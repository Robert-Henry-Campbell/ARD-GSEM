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

      if (neff < 1000) {
        log_warn("munge", sprintf("%s %s: Neff=%.0f (very low; expect underpowered h2)", trait, sex, neff))
        warnings_list <- c(warnings_list, sprintf("%s: Neff=%.0f", trait, neff))
      }

      log_info("munge", sprintf("%s %s: Neff=%.0f (cases=%d, controls=%d)",
                                trait, sex, neff, cc$n_cases, cc$n_controls))

      polarity <- check_polarity_sentinel(dt)
      if (!is.na(polarity$polarity_correct)) {
        if (polarity$polarity_correct) {
          log_info("munge", sprintf("Polarity check PASSED (%s)", polarity$message))
        } else {
          log_warn("munge", sprintf("Polarity check FAILED (%s)", polarity$message))
        }
      }

      out_dir <- file.path(config$paths$output_dir, sex, "munge")
      munge_input <- file.path(out_dir, paste0(trait, "_pre_munge.tsv"))
      munge_dt <- dt[, .(SNP, A1, A2, effect = beta, SE = se, P = pval, N = neff)]
      fwrite(munge_dt, munge_input, sep = "\t")

      GenomicSEM::munge(
        files = munge_input,
        hm3 = config$paths$hm3_snplist,
        trait.names = trait,
        N = neff,
        maf.filter = config$munge$maf_threshold,
        log.name = file.path(out_dir, paste0(trait, "_munge"))
      )

      munged_path <- paste0(munge_input, ".sumstats.gz")
      if (!file.exists(munged_path)) {
        alt_path <- file.path(out_dir, paste0(trait, ".sumstats.gz"))
        munged_candidates <- list.files(out_dir, pattern = paste0(trait, ".*\\.sumstats\\.gz$"),
                                        full.names = TRUE)
        if (length(munged_candidates) > 0) munged_path <- munged_candidates[1]
      }

      if (file.exists(munged_path)) {
        munged_files <- c(munged_files, munged_path)
        n_snps <- length(readLines(gzcon(file(munged_path, "rb")), n = -1)) - 1
        log_info("munge", sprintf("%s %s: %s HM3 SNPs retained",
                                  trait, sex, format(n_snps, big.mark = ",")))
      } else {
        log_warn("munge", sprintf("%s %s: munged file not found at expected path", trait, sex))
      }

      file.remove(munge_input)

    }, error = function(e) {
      log_error("munge", sprintf("%s %s FAILED: %s", trait, sex, e$message))
      warnings_list <<- c(warnings_list, sprintf("%s: %s", trait, e$message))
    })
  }

  log_info("munge", sprintf("Munge complete: %d/%d %s traits processed",
                            length(munged_files), length(traits), sex))

  write_stage_manifest("munge", sex, config, munged_files, warnings_list)

  list(munged_files = munged_files, traits = traits, n_munged = length(munged_files))
}
