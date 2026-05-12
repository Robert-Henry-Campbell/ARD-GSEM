run_munge <- function(config, sex) {
  log_info("munge", sprintf("=== Munge stage: %s ===", sex))

  variants_manifest <- load_variants_manifest(config$paths$variants_manifest)
  traits <- discover_traits(config, sex)
  munged_files <- character(0)
  warnings_list <- character(0)

  for (trait in traits) {
    log_info("munge", sprintf("Processing %s %s...", trait, sex))

    tryCatch({
      sumstats_file <- file.path(config$paths$sumstats_dir, sex,
                                 paste0(trait, ".gwas.imputed_v3.", sex, ".tsv.bgz"))
      dt <- fread(cmd = paste("zcat", shQuote(sumstats_file)))

      n_raw <- nrow(dt)
      # MAF filtering is applied here, then again inside GenomicSEM::munge() (HM3 join),
      # and again inside GenomicSEM::sumstats() (1000G alignment). Redundant but cheap;
      # each downstream step expects its own MAF guard.
      dt <- filter_variants(dt,
                            maf_threshold = config$munge$maf_threshold,
                            info_threshold = NULL)
      log_debug("munge", sprintf("After filtering: %s/%s variants",
                                 format(nrow(dt), big.mark = ","),
                                 format(n_raw, big.mark = ",")))

      dt <- parse_variant_column(dt)
      dt <- map_rsids(dt, variants_manifest)
      log_debug("munge", sprintf("rsid mapped: %s variants", format(nrow(dt), big.mark = ",")))

      dt[, A1 := alt]
      dt[, A2 := ref]

      cc <- get_case_control(config, sex, trait)
      neff <- compute_neff(cc$n_cases, cc$n_controls)

      if (is.na(neff)) {
        log_warn("munge", sprintf("%s %s: Neff is NA (cases=%s controls=%s); skipping",
                                  trait, sex,
                                  as.character(cc$n_cases), as.character(cc$n_controls)))
        warnings_list <- c(warnings_list, sprintf("%s: Neff NA, skipped", trait))
        next
      }
      if (neff < 1000) {
        log_warn("munge", sprintf("%s %s: Neff=%.0f (very low; expect underpowered h2)", trait, sex, neff))
        warnings_list <- c(warnings_list, sprintf("%s: Neff=%.0f", trait, neff))
      }

      K <- cc$n_cases / (cc$n_cases + cc$n_controls)
      log_info("munge", sprintf(
        "%s %s: Neff=%.0f (cases=%d, controls=%d), K=%.4f -> converting BETA/SE to log(OR) via beta/(K(1-K))",
        trait, sex, neff, cc$n_cases, cc$n_controls, K))

      K_warn <- config$munge$K_warn_threshold %||% 0.01
      K_reject <- config$munge$K_reject_threshold %||% 1e-9
      conv <- linear_to_logor(beta = dt$beta, se = dt$se, K = K,
                              reject_threshold = K_reject,
                              warn_threshold = K_warn)
      if (isTRUE(conv$warn)) {
        log_warn("munge", sprintf(
          "%s %s: K=%.4f outside [%.3f, %.3f]; linear->log(OR) approximation has growing higher-order bias for SNPs with non-trivial effect",
          trait, sex, K, K_warn, 1 - K_warn))
        warnings_list <- c(warnings_list, sprintf("%s: K=%.4f extreme", trait, K))
      }
      dt[, `:=`(beta = conv$beta, se = conv$se)]

      sentinel_cfg <- config$munge$polarity_sentinel
      polarity <- check_polarity_sentinel(dt,
                                          sentinel_rsid = sentinel_cfg$rsid,
                                          expected_a1 = sentinel_cfg$a1)
      if (!is.na(polarity$polarity_correct)) {
        if (polarity$polarity_correct) {
          log_info("munge", sprintf("Polarity check PASSED (%s)", polarity$message))
        } else {
          log_warn("munge", sprintf("Polarity check FAILED (%s)", polarity$message))
        }
      } else {
        log_debug("munge", sprintf("Polarity: %s", polarity$message))
      }

      out_dir <- file.path(config$paths$output_dir, sex, "munge")
      dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
      prefix <- paste0(sex, "_")
      munge_input <- file.path(out_dir, paste0(prefix, trait, "_pre_munge.tsv"))
      munge_dt <- dt[, .(SNP, CHR = chr, BP = pos, A1, A2, MAF = minor_AF,
                         effect = beta, SE = se, P = pval, N = neff)]
      fwrite(munge_dt, munge_input, sep = "\t")

      # GenomicSEM::munge() writes paste0(trait.names[i], ".sumstats(.gz)") and the log file
      # as paths relative to CWD (verified by reading munge_main.R: there is even a literal
      # log line "saved ... in the current working directory"). Without forcing CWD to
      # out_dir during the call, the munged file lands wherever R happened to be, i.e. at the
      # project root, exactly the symptom seen on the previous server run.
      hm3_abs <- normalizePath(config$paths$hm3_snplist, mustWork = TRUE)
      munge_input_abs <- normalizePath(munge_input, mustWork = TRUE)
      munge_trait <- paste0(prefix, trait)
      stale <- file.path(out_dir, paste0(munge_trait, ".sumstats.gz"))
      if (file.exists(stale)) file.remove(stale)
      saved_wd <- getwd()
      tryCatch({
        setwd(out_dir)
        GenomicSEM::munge(
          files = munge_input_abs,
          hm3 = hm3_abs,
          trait.names = munge_trait,
          N = neff,
          maf.filter = config$munge$maf_threshold,
          log.name = paste0(munge_trait, "_munge")
        )
      }, finally = setwd(saved_wd))

      munged_path <- file.path(out_dir, paste0(munge_trait, ".sumstats.gz"))
      if (!file.exists(munged_path)) {
        legacy <- paste0(munge_input, ".sumstats.gz")
        if (file.exists(legacy)) munged_path <- legacy
      }

      if (file.exists(munged_path)) {
        munged_files <- c(munged_files, munged_path)
        n_snps <- count_gz_lines(munged_path) - 1L
        log_info("munge", sprintf("%s %s: %s HM3 SNPs retained",
                                  trait, sex, format(n_snps, big.mark = ",")))
      } else {
        log_warn("munge", sprintf("%s %s: munged file not found at expected path", trait, sex))
      }

      # Keep munge_input (pre-munge wide TSV): the sumstats stage reads it for the
      # raw per-trait BETA/SE on the log(OR) scale.

    }, error = function(e) {
      log_error("munge", sprintf("%s %s FAILED: %s", trait, sex, e$message))
      warnings_list <<- c(warnings_list, sprintf("%s: %s", trait, e$message))
    })
  }

  log_info("munge", sprintf("Munge complete: %d/%d %s traits processed",
                            length(munged_files), length(traits), sex))

  if (length(munged_files) == length(traits)) {
    write_stage_manifest("munge", sex, config, munged_files, warnings_list)
  } else {
    munged_traits <- sub(paste0("^", sex, "_"), "",
                         sub("\\.sumstats\\.gz$", "", basename(munged_files)))
    munged_traits <- sub("_pre_munge\\.tsv$", "", munged_traits)
    failed <- setdiff(traits, munged_traits)
    log_warn("munge", sprintf(
      "Only %d/%d traits munged successfully; NOT writing stage manifest so --resume will retry. Failed: %s",
      length(munged_files), length(traits),
      paste(failed, collapse = ", ")))
  }

  list(munged_files = munged_files, traits = traits, n_munged = length(munged_files))
}
