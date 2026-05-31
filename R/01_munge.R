run_munge <- function(config, sex) {
  if (identical(sex, "bothsex")) return(run_munge_panukb(config))
  if (identical(sex, "bothsex_meta")) return(run_munge_bothsex_meta(config))
  run_munge_neale(config, sex)
}

run_munge_neale <- function(config, sex) {
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

      verify_genome_build(dt, config$genome_build_check,
                          build_label = config$gwas$genome_build %||% "GRCh37")

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
        "%s %s: Neff=%.0f (cases=%d, controls=%d), K_sample=%.4f -> raw linear BETA/SE retained; linprob=TRUE handles conversion in sumstats()",
        trait, sex, neff, cc$n_cases, cc$n_controls, K))

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

run_munge_panukb <- function(config) {
  sex <- "bothsex"
  log_info("munge", "=== Munge stage: bothsex (Pan-UKB) ===")

  vm <- load_panukb_variants_manifest(
    config$paths$panukb_variants_manifest,
    cols = config$panukb$variants_columns)

  traits <- discover_traits(config, sex)
  if (length(traits) == 0L) {
    log_fatal("munge", sprintf("No Pan-UKB bothsex sumstats found under %s",
                                file.path(config$paths$sumstats_dir, sex)))
  }
  log_info("munge", sprintf("Pan-UKB traits discovered: %s",
                            paste(traits, collapse = ", ")))

  col_defaults <- list(chr = "chr", pos = "pos", ref = "ref", alt = "alt",
                       beta = "beta_EUR", se = "se_EUR",
                       neglog10_pval = "neglog10_pval_EUR",
                       low_confidence = "low_confidence_EUR")
  cols <- resolve_column_map(config$panukb$columns, col_defaults,
                             label = "config$panukb$columns")
  required_keys <- names(col_defaults)

  munged_files <- character(0)
  warnings_list <- character(0)

  for (trait in traits) {
    log_info("munge", sprintf("Processing %s bothsex...", trait))
    tryCatch({
      sumstats_file <- file.path(config$paths$sumstats_dir, sex,
                                 sprintf("icd10-%s-both_sexes.tsv.bgz", trait))
      if (!file.exists(sumstats_file)) {
        log_fatal("munge", sprintf("Pan-UKB sumstat missing: %s", sumstats_file))
      }

      need <- vapply(required_keys, function(k) cols[[k]], character(1))
      header <- names(data.table::fread(cmd = paste("zcat", shQuote(sumstats_file)),
                                          nrows = 0L))
      missing_cols <- setdiff(need, header)
      if (length(missing_cols) > 0L) {
        log_fatal("munge", sprintf(
          "Pan-UKB sumstat %s missing column(s): %s. Actual header: %s",
          basename(sumstats_file), paste(missing_cols, collapse = ", "),
          paste(header, collapse = ", ")))
      }

      classes <- list()
      classes[[cols$chr]] <- "character"
      classes[[cols$low_confidence]] <- "character"
      dt <- data.table::fread(cmd = paste("zcat", shQuote(sumstats_file)),
                              select = unname(need), colClasses = classes,
                              na.strings = c("", "NA"))
      n_raw <- nrow(dt)

      data.table::setnames(dt, old = unname(need), new = required_keys)

      dt[, low_confidence := tolower(low_confidence) == "true"]
      dt <- dt[low_confidence == FALSE & !is.na(beta) & !is.na(se) &
               !is.na(neglog10_pval)]
      log_debug("munge", sprintf("After low_confidence/NA filter: %s/%s variants",
                                  format(nrow(dt), big.mark = ","),
                                  format(n_raw, big.mark = ",")))

      dt[, pval := 10^(-neglog10_pval)]

      dt <- merge_panukb_variants(dt, vm)

      info_thr <- config$munge$info_threshold %||% 0.8
      dt <- dt[info >= info_thr]
      log_debug("munge", sprintf("After INFO>=%.2f: %s variants",
                                  info_thr, format(nrow(dt), big.mark = ",")))

      dt[, MAF := pmin(af_EUR, 1 - af_EUR)]
      dt <- dt[MAF >= config$munge$maf_threshold]
      log_debug("munge", sprintf("After MAF>=%.2f: %s variants",
                                  config$munge$maf_threshold,
                                  format(nrow(dt), big.mark = ",")))

      dt[, `:=`(A1 = alt, A2 = ref, SNP = rsid)]

      verify_genome_build(dt, config$genome_build_check,
                          build_label = config$gwas$genome_build %||% "GRCh37")

      cc <- get_case_control(config, sex, trait)
      neff <- compute_neff(cc$n_cases, cc$n_controls)
      if (is.na(neff)) {
        log_warn("munge", sprintf("%s bothsex: Neff NA (cases=%s controls=%s); skipping",
                                  trait, as.character(cc$n_cases),
                                  as.character(cc$n_controls)))
        warnings_list <- c(warnings_list, sprintf("%s: Neff NA, skipped", trait))
        next
      }
      if (neff < 1000) {
        log_warn("munge", sprintf("%s bothsex: Neff=%.0f (very low)", trait, neff))
        warnings_list <- c(warnings_list, sprintf("%s: Neff=%.0f", trait, neff))
      }

      K <- cc$n_cases / (cc$n_cases + cc$n_controls)
      log_info("munge", sprintf(
        "%s bothsex: Neff=%.0f (cases=%d, controls=%d), K_sample=%.4f; SAIGE log(OR) input -> se.logit=TRUE in sumstats()",
        trait, neff, cc$n_cases, cc$n_controls, K))

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
      munge_dt <- dt[, .(SNP, CHR = chr, BP = pos, A1, A2, MAF,
                          effect = beta, SE = se, P = pval, N = neff)]
      data.table::fwrite(munge_dt, munge_input, sep = "\t")

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
      if (file.exists(munged_path)) {
        munged_files <- c(munged_files, munged_path)
        n_snps <- count_gz_lines(munged_path) - 1L
        log_info("munge", sprintf("%s bothsex: %s HM3 SNPs retained",
                                    trait, format(n_snps, big.mark = ",")))
      } else {
        log_warn("munge", sprintf("%s bothsex: munged file not found at expected path", trait))
      }
    }, error = function(e) {
      log_error("munge", sprintf("%s bothsex FAILED: %s", trait, e$message))
      warnings_list <<- c(warnings_list, sprintf("%s: %s", trait, e$message))
    })
  }

  log_info("munge", sprintf("Munge complete: %d/%d bothsex traits processed",
                            length(munged_files), length(traits)))

  if (length(munged_files) == length(traits)) {
    write_stage_manifest("munge", sex, config, munged_files, warnings_list)
  } else {
    log_warn("munge", sprintf(
      "Only %d/%d bothsex traits munged; NOT writing stage manifest so --resume will retry.",
      length(munged_files), length(traits)))
  }

  list(munged_files = munged_files, traits = traits, n_munged = length(munged_files))
}

run_munge_bothsex_meta <- function(config) {
  sex <- "bothsex_meta"
  log_info("munge", "=== Munge stage: bothsex_meta (FinnGen + UKBB EUR meta + FG-only) ===")

  meta_cfg <- config$bothsex_meta
  if (is.null(meta_cfg)) {
    log_fatal("munge", "config$bothsex_meta block missing from pipeline.yaml")
  }

  het_p_min <- meta_cfg$het_p_min %||% 1.0e-3
  chain_path <- meta_cfg$chain_file
  if (!is_absolute_path(chain_path)) {
    chain_path <- file.path(config$project$root %||% getwd(), chain_path)
  }

  # R9: intersect discovered sumstats files with the manifest's phenotype list.
  # The manifest generator drops traits without sufficient cohort coverage;
  # files for those traits may still sit in sumstats_dir. Without this filter,
  # the munge loop would error on each missing-from-manifest trait, refuse to
  # write the stage manifest because length(munged_files) != length(traits),
  # and --resume would retry forever.
  rda_path <- meta_cfg$manifest_rda
  if (!is_absolute_path(rda_path)) {
    rda_path <- file.path(config$project$root %||% getwd(), rda_path)
  }
  env <- new.env()
  load(rda_path, envir = env)
  manifest <- get(ls(env)[1L], envir = env)
  manifest_traits <- as.character(manifest$phenotype)

  # Discover trait files + their type (meta vs fg). discover_traits_with_type
  # returns data.table(trait, file_type, file_path).
  traits_dt <- discover_traits_with_type(config)
  if (nrow(traits_dt) == 0L) {
    log_fatal("munge", sprintf("No bothsex_meta sumstats found under %s",
                                file.path(config$paths$sumstats_dir, sex)))
  }
  dropped <- setdiff(traits_dt$trait, manifest_traits)
  if (length(dropped) > 0L) {
    log_warn("munge", sprintf(
      "bothsex_meta: dropping %d trait file(s) not in manifest (likely insufficient cohort coverage): %s",
      length(dropped), paste(dropped, collapse = ", ")))
    traits_dt <- traits_dt[trait %in% manifest_traits]
  }
  if (nrow(traits_dt) == 0L) {
    log_fatal("munge", "No bothsex_meta traits remain after manifest intersection")
  }
  type_counts <- table(traits_dt$file_type)
  log_info("munge", sprintf("bothsex_meta traits to process: %d (%s)",
                            nrow(traits_dt),
                            paste(sprintf("%d %s", as.integer(type_counts), names(type_counts)),
                                  collapse = ", ")))

  munged_files <- character(0)
  warnings_list <- character(0)

  for (i in seq_len(nrow(traits_dt))) {
    trait <- traits_dt$trait[i]
    file_type <- traits_dt$file_type[i]
    file_path <- traits_dt$file_path[i]
    log_info("munge", sprintf("Processing %s bothsex_meta [%s]...", trait, file_type))

    result <- tryCatch({
      cc <- get_case_control(config, sex, trait)
      if (identical(file_type, "meta")) {
        .munge_one_meta_trait(trait, file_path, config, meta_cfg, cc,
                              het_p_min = het_p_min, chain_path = chain_path)
      } else if (identical(file_type, "fg")) {
        .munge_one_fg_trait(trait, file_path, config, meta_cfg, cc,
                            chain_path = chain_path)
      } else {
        log_warn("munge", sprintf("%s: unknown file_type '%s'; skipping", trait,
                                  as.character(file_type)))
        list(success = FALSE, warning = sprintf("%s: unknown file_type", trait))
      }
    }, error = function(e) {
      log_error("munge", sprintf("%s bothsex_meta [%s] FAILED: %s",
                                  trait, file_type, e$message))
      list(success = FALSE, warning = sprintf("%s: %s", trait, e$message))
    })

    if (isTRUE(result$success)) {
      munged_files <- c(munged_files, result$munged_path)
    }
    if (!is.null(result$warning)) {
      warnings_list <- c(warnings_list, result$warning)
    }
  }

  log_info("munge", sprintf("Munge complete: %d/%d bothsex_meta traits processed",
                            length(munged_files), nrow(traits_dt)))

  if (length(munged_files) == nrow(traits_dt)) {
    write_stage_manifest("munge", sex, config, munged_files, warnings_list)
  } else {
    log_warn("munge", sprintf(
      "Only %d/%d bothsex_meta traits munged; NOT writing stage manifest so --resume will retry.",
      length(munged_files), nrow(traits_dt)))
  }

  list(munged_files = munged_files, traits = traits_dt$trait,
       n_munged = length(munged_files))
}

# Helper: shared post-processing tail. Takes a cleaned dt (with chr, pos in
# GRCh37 already; A1/A2/SNP set; beta/se/pval canonical) + scalar Neff +
# pooled cohort counts. Performs polarity check, writes pre_munge.tsv, calls
# GenomicSEM::munge, returns list(success, munged_path, warning).
.write_munge_output <- function(trait, dt, neff_scalar, n_cases_pooled,
                                n_controls_pooled, config, sex = "bothsex_meta",
                                file_type_label = "meta") {
  K <- if (!is.na(n_cases_pooled) && !is.na(n_controls_pooled) &&
           (n_cases_pooled + n_controls_pooled) > 0) {
    n_cases_pooled / (n_cases_pooled + n_controls_pooled)
  } else NA_real_
  log_info("munge", sprintf(
    "%s bothsex_meta [%s]: pooled n_cases=%s, n_controls=%s, K_sample=%s, scalar N=%.0f; betas already log(OR) -> se.logit=TRUE in sumstats()",
    trait, file_type_label,
    format(as.numeric(n_cases_pooled), big.mark = ","),
    format(as.numeric(n_controls_pooled), big.mark = ","),
    if (is.na(K)) "NA" else sprintf("%.4f", K),
    neff_scalar))

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
  } else if (isTRUE(polarity$sentinel_missing)) {
    log_warn("munge", sprintf(
      "Polarity sentinel %s not found post-liftover (%s); polarity unverified for %s",
      sentinel_cfg$rsid, polarity$message, trait))
  } else {
    log_debug("munge", sprintf("Polarity: %s", polarity$message))
  }

  out_dir <- file.path(config$paths$output_dir, sex, "munge")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  prefix <- paste0(sex, "_")
  munge_input <- file.path(out_dir, paste0(prefix, trait, "_pre_munge.tsv"))
  munge_dt <- dt[, .(SNP, CHR = chr, BP = pos, A1, A2, MAF,
                      effect = beta, SE = se, P = pval, N = neff_scalar)]
  data.table::fwrite(munge_dt, munge_input, sep = "\t")

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
      N = neff_scalar,
      maf.filter = config$munge$maf_threshold,
      log.name = paste0(munge_trait, "_munge")
    )
  }, finally = setwd(saved_wd))

  munged_path <- file.path(out_dir, paste0(munge_trait, ".sumstats.gz"))
  if (file.exists(munged_path)) {
    n_snps <- count_gz_lines(munged_path) - 1L
    log_info("munge", sprintf("%s bothsex_meta [%s]: %s HM3 SNPs retained",
                                trait, file_type_label, format(n_snps, big.mark = ",")))
    return(list(success = TRUE, munged_path = munged_path, warning = NULL))
  }
  log_warn("munge", sprintf("%s bothsex_meta [%s]: munged file not found at expected path",
                              trait, file_type_label))
  list(success = FALSE,
       warning = sprintf("%s: munged file not produced", trait))
}

# Helper: process one _meta_out.tsv.gz file (2-way FG+UKBB inverse-variance meta).
.munge_one_meta_trait <- function(trait, sumstats_file, config, meta_cfg, cc,
                                  het_p_min, chain_path) {
  col_defaults <- list(chr = "#CHR", pos = "POS", ref = "REF", alt = "ALT",
                       snp = "SNP", rsid = "rsid",
                       beta = "all_inv_var_meta_beta",
                       se = "all_inv_var_meta_sebeta",
                       pval = "all_inv_var_meta_p",
                       mlogp = "all_inv_var_meta_mlogp",
                       het_p = "all_inv_var_het_p",
                       all_meta_N = "all_meta_N")
  cols <- resolve_column_map(meta_cfg$columns, col_defaults,
                             label = "config$bothsex_meta$columns")

  af_cohorts <- meta_cfg$af_cohorts
  if (length(af_cohorts) == 0L) {
    log_fatal("munge", "config$bothsex_meta$af_cohorts is empty; must list at least one cohort")
  }
  af_cols <- vapply(af_cohorts, function(co) co$af_col, character(1))
  beta_cols <- vapply(af_cohorts, function(co) co$beta_col, character(1))

  need <- unique(c(unname(unlist(cols)), unname(af_cols), unname(beta_cols)))
  header <- names(data.table::fread(cmd = paste("zcat", shQuote(sumstats_file)),
                                      nrows = 0L))
  missing_cols <- setdiff(need, header)
  if (length(missing_cols) > 0L) {
    log_fatal("munge", sprintf(
      "bothsex_meta meta sumstat %s missing column(s): %s. Actual header: %s",
      basename(sumstats_file), paste(missing_cols, collapse = ", "),
      paste(header, collapse = ", ")))
  }

  classes <- list()
  classes[[cols$chr]]  <- "character"
  classes[[cols$pos]]  <- "integer"   # R4
  classes[[cols$rsid]] <- "character"
  classes[[cols$snp]]  <- "character"
  dt <- data.table::fread(cmd = paste("zcat", shQuote(sumstats_file)),
                          select = need, colClasses = classes,
                          na.strings = c("", "NA", "."))
  n_raw <- nrow(dt)

  core_renames <- c(chr = cols$chr, pos = cols$pos, ref = cols$ref,
                    alt = cols$alt, snp = cols$snp, rsid = cols$rsid,
                    beta = cols$beta, se = cols$se, pval = cols$pval,
                    mlogp = cols$mlogp, het_p = cols$het_p,
                    all_meta_N = cols$all_meta_N)
  for (new_nm in names(core_renames)) {
    old_nm <- core_renames[[new_nm]]
    if (!identical(old_nm, new_nm) && old_nm %in% names(dt)) {
      data.table::setnames(dt, old_nm, new_nm)
    }
  }
  dt[, het_p := as.numeric(het_p)]                # R5
  dt[, all_meta_N := as.integer(all_meta_N)]

  # R6: cohort-count sanity check.
  max_cohorts_in_data <- suppressWarnings(max(dt$all_meta_N, na.rm = TRUE))
  if (is.finite(max_cohorts_in_data) && max_cohorts_in_data > length(af_cohorts)) {
    log_warn("munge", sprintf(
      "%s bothsex_meta [meta]: max(all_meta_N)=%d but only %d cohort(s) configured in af_cohorts (%s). The all_inv_var_meta_* betas reflect the wider meta; AF / per-SNP Neff will reflect only the configured cohorts. Verify the sumstats source is EUR-only.",
      trait, max_cohorts_in_data, length(af_cohorts),
      paste(vapply(af_cohorts, function(co) co$name, character(1)), collapse = ", ")))
  }

  dt <- dt[!is.na(beta) & !is.na(se) & !is.na(pval)]
  log_debug("munge", sprintf("[meta] After NA filter: %s/%s variants",
                              format(nrow(dt), big.mark = ","),
                              format(n_raw, big.mark = ",")))
  dt <- dt[!is.na(rsid) & nzchar(rsid)]

  manifest_row <- cc$manifest_row
  n_weights <- vapply(af_cohorts, function(co) {
    nca <- manifest_row[[co$manifest_n_cases]]
    nco <- manifest_row[[co$manifest_n_controls]]
    if (is.null(nca) || is.null(nco) || length(nca) == 0L || length(nco) == 0L ||
        is.na(nca) || is.na(nco)) return(0)
    as.numeric(nca) + as.numeric(nco)
  }, numeric(1))
  if (sum(n_weights) <= 0) {
    log_warn("munge", sprintf("%s bothsex_meta [meta]: all cohort weights zero/NA; skipping", trait))
    return(list(success = FALSE,
                warning = sprintf("%s: no cohort weights", trait)))
  }
  af_mat <- as.matrix(dt[, af_cols, with = FALSE])
  mode(af_mat) <- "numeric"
  meta_af_vec <- nweighted_meta_af(af_mat, n_weights)
  dt[, meta_af := meta_af_vec]
  dt[, MAF := pmin(meta_af, 1 - meta_af)]
  n_before_maf <- nrow(dt)
  dt <- dt[!is.na(MAF) & MAF >= config$munge$maf_threshold]
  log_debug("munge", sprintf("[meta] After MAF>=%.2f: %s/%s variants",
                              config$munge$maf_threshold,
                              format(nrow(dt), big.mark = ","),
                              format(n_before_maf, big.mark = ",")))

  # Heterogeneity QC.
  n_before_het <- nrow(dt)
  dt <- dt[is.na(het_p) | het_p >= het_p_min]
  log_info("munge", sprintf(
    "[meta] het_filter: %s -> %s SNPs (dropped %s with het_p < %.1e)",
    format(n_before_het, big.mark = ","), format(nrow(dt), big.mark = ","),
    format(n_before_het - nrow(dt), big.mark = ","), het_p_min))

  dt <- liftover_grch38_to_grch37(dt, chain_path)
  dt[, `:=`(A1 = alt, A2 = ref, SNP = rsid)]
  verify_genome_build(dt, config$genome_build_check,
                      build_label = config$gwas$genome_build %||% "GRCh37")

  # Scalar N from per-SNP Neff median (S4).
  neff_vec <- per_snp_neff(dt, af_cohorts, manifest_row)
  neff_pos <- neff_vec[neff_vec > 0]
  if (length(neff_pos) == 0L) {
    log_warn("munge", sprintf("%s bothsex_meta [meta]: no SNP has any cohort contribution; skipping", trait))
    return(list(success = FALSE,
                warning = sprintf("%s: per-SNP Neff all zero", trait)))
  }
  neff_scalar <- median(neff_pos, na.rm = TRUE)
  if (is.na(neff_scalar) || neff_scalar <= 0) {
    log_warn("munge", sprintf("%s bothsex_meta [meta]: median per-SNP Neff is NA/0; skipping", trait))
    return(list(success = FALSE, warning = sprintf("%s: Neff NA", trait)))
  }
  if (neff_scalar < 1000) {
    log_warn("munge", sprintf("%s bothsex_meta [meta]: median per-SNP Neff=%.0f (very low)",
                              trait, neff_scalar))
  }
  log_info("munge", sprintf(
    "[meta] per-SNP Neff: min=%.0f, median=%.0f (passed to GenomicSEM as scalar N), max=%.0f",
    min(neff_pos), neff_scalar, max(neff_vec, na.rm = TRUE)))

  .write_munge_output(trait, dt, neff_scalar,
                       n_cases_pooled = cc$n_cases,
                       n_controls_pooled = cc$n_controls,
                       config = config, file_type_label = "meta")
}

# Helper: process one _fg.tsv.gz file (raw single-cohort FinnGen GWAS).
# Schema is entirely different from the meta files -- no all_inv_var_meta_*,
# no het_p, lowercase column names, rsids (plural, may be comma-separated).
.munge_one_fg_trait <- function(trait, sumstats_file, config, meta_cfg, cc,
                                chain_path) {
  fg_defaults <- list(chr = "#chrom", pos = "pos", ref = "ref", alt = "alt",
                      rsids = "rsids", beta = "beta", se = "sebeta",
                      pval = "pval", mlogp = "mlogp", af_alt = "af_alt")
  fg_cols <- resolve_column_map(meta_cfg$fg_only_columns, fg_defaults,
                                label = "config$bothsex_meta$fg_only_columns")
  need <- unique(unname(unlist(fg_cols)))
  header <- names(data.table::fread(cmd = paste("zcat", shQuote(sumstats_file)),
                                      nrows = 0L))
  missing_cols <- setdiff(need, header)
  if (length(missing_cols) > 0L) {
    log_fatal("munge", sprintf(
      "bothsex_meta fg sumstat %s missing column(s): %s. Actual header: %s",
      basename(sumstats_file), paste(missing_cols, collapse = ", "),
      paste(header, collapse = ", ")))
  }

  classes <- list()
  classes[[fg_cols$chr]]   <- "character"
  classes[[fg_cols$pos]]   <- "integer"    # R4
  classes[[fg_cols$rsids]] <- "character"
  dt <- data.table::fread(cmd = paste("zcat", shQuote(sumstats_file)),
                          select = need, colClasses = classes,
                          na.strings = c("", "NA", "."))
  n_raw <- nrow(dt)

  core_renames <- c(chr = fg_cols$chr, pos = fg_cols$pos, ref = fg_cols$ref,
                    alt = fg_cols$alt, rsids = fg_cols$rsids,
                    beta = fg_cols$beta, se = fg_cols$se, pval = fg_cols$pval,
                    mlogp = fg_cols$mlogp, af_alt = fg_cols$af_alt)
  for (new_nm in names(core_renames)) {
    old_nm <- core_renames[[new_nm]]
    if (!identical(old_nm, new_nm) && old_nm %in% names(dt)) {
      data.table::setnames(dt, old_nm, new_nm)
    }
  }

  # rsids -> first rsid only (some rows have comma-separated multi-rsids).
  dt[, rsid := sub(",.*", "", rsids)]
  dt <- dt[!is.na(rsid) & nzchar(rsid)]

  dt <- dt[!is.na(beta) & !is.na(se) & !is.na(pval)]
  log_debug("munge", sprintf("[fg] After NA + rsid filter: %s/%s variants",
                              format(nrow(dt), big.mark = ","),
                              format(n_raw, big.mark = ",")))

  # MAF: single-cohort AF, no weighted mean needed.
  dt[, af_alt := as.numeric(af_alt)]
  dt[, MAF := pmin(af_alt, 1 - af_alt)]
  n_before_maf <- nrow(dt)
  dt <- dt[!is.na(MAF) & MAF >= config$munge$maf_threshold]
  log_debug("munge", sprintf("[fg] After MAF>=%.2f: %s/%s variants",
                              config$munge$maf_threshold,
                              format(nrow(dt), big.mark = ","),
                              format(n_before_maf, big.mark = ",")))

  # No heterogeneity filter (single cohort -> no het_p).
  # No R6 cohort-count check (no all_meta_N column).

  dt <- liftover_grch38_to_grch37(dt, chain_path)
  dt[, `:=`(A1 = alt, A2 = ref, SNP = rsid)]
  verify_genome_build(dt, config$genome_build_check,
                      build_label = config$gwas$genome_build %||% "GRCh37")

  # Scalar N from manifest's FG-only cohort counts. No per-SNP variation since
  # only one cohort contributes. fg cohort is the first entry in af_cohorts by
  # convention; pull its manifest_n_cases / _n_controls keys.
  manifest_row <- cc$manifest_row
  af_cohorts <- meta_cfg$af_cohorts
  fg_co <- NULL
  for (co in af_cohorts) {
    if (identical(co$name, "fg")) { fg_co <- co; break }
  }
  if (is.null(fg_co)) {
    log_fatal("munge", "config$bothsex_meta$af_cohorts must contain a cohort named 'fg' for fg_only file handling")
  }
  n_cases_fg    <- as.numeric(manifest_row[[fg_co$manifest_n_cases]])
  n_controls_fg <- as.numeric(manifest_row[[fg_co$manifest_n_controls]])
  neff_scalar <- compute_neff(n_cases_fg, n_controls_fg)
  if (is.na(neff_scalar) || neff_scalar <= 0) {
    log_warn("munge", sprintf("%s bothsex_meta [fg]: Neff NA/0 (n_cases=%s, n_controls=%s); skipping",
                              trait, as.character(n_cases_fg), as.character(n_controls_fg)))
    return(list(success = FALSE, warning = sprintf("%s: Neff NA", trait)))
  }
  if (neff_scalar < 1000) {
    log_warn("munge", sprintf("%s bothsex_meta [fg]: Neff=%.0f (very low)",
                              trait, neff_scalar))
  }
  log_info("munge", sprintf("[fg] single-cohort N: n_cases_fg=%.0f, n_controls_fg=%.0f, Neff=%.0f",
                            n_cases_fg, n_controls_fg, neff_scalar))

  .write_munge_output(trait, dt, neff_scalar,
                       n_cases_pooled = n_cases_fg,
                       n_controls_pooled = n_controls_fg,
                       config = config, file_type_label = "fg")
}
