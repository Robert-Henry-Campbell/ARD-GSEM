run_gwas_smoke <- function(config, sex) {
  log_info("smoke", sprintf("=== userGWAS smoke test: %s ===", sex))

  ldsc_path <- file.path(config$paths$output_dir, sex, "ldsc", "ldsc_full.rds")
  ldsc_output <- readRDS(ldsc_path)

  snp_path <- file.path(config$paths$output_dir, sex, "sumstats",
                        paste0(sex, "_snp_sumstats.rds"))
  if (!file.exists(snp_path)) {
    legacy <- file.path(config$paths$output_dir, sex, "sumstats", "snp_sumstats.rds")
    if (file.exists(legacy)) snp_path <- legacy
    else log_fatal("smoke", sprintf("snp_sumstats not found at %s", snp_path))
  }
  snp_sumstats <- readRDS(snp_path)

  n_test <- as.integer(config$smoke$n_snps %||% 100L)
  if (nrow(snp_sumstats) <= n_test) {
    snp_test <- snp_sumstats
  } else {
    set.seed(20260516L)
    snp_test <- snp_sumstats[sample.int(nrow(snp_sumstats), n_test), ]
  }
  log_info("smoke", sprintf("Sampling %d SNPs from %d for smoke test",
                             nrow(snp_test), nrow(snp_sumstats)))

  # Prefer the split-sample-locked bifactor spec (04_cfa) so the smoke gate tests the
  # model the GWAS will actually run; fall back to the legacy chapter model otherwise.
  spec_path <- file.path(config$paths$output_dir, sex, "cfa", "model_spec.rds")
  cov_for_gwas <- ldsc_output
  if (file.exists(spec_path)) {
    bspec <- readRDS(spec_path)
    log_info("smoke", sprintf("Bifactor smoke: %d factor(s) [%s] over %d indicators",
                              length(bspec$factors), paste(bspec$factors, collapse = ", "),
                              length(bspec$general_traits)))
    cov_for_gwas <- subset_covstruc(ldsc_output, bspec$general_traits)
    snp_test <- subset_snps(snp_test, colnames(cov_for_gwas$S))
    spec <- list(model = bspec$gwas_model, sub = bspec$sub, factors = bspec$factors)
    std_lv <- TRUE
  } else {
    h2 <- fread(file.path(config$paths$output_dir, sex, "ldsc", "h2_qc.csv"))
    retained <- h2[pass == TRUE]$trait
    categories <- load_trait_categories(config, sex)
    loadings_path <- file.path(config$paths$output_dir, sex, "cfa", "factor_loadings.csv")
    ordering_table <- if (file.exists(loadings_path)) fread(loadings_path) else NULL
    apriori_model <- build_apriori_model(
      retained, categories,
      min_indicators = config$cfa$min_indicators_per_factor,
      ordering_table = ordering_table,
      code_format = if (identical(sex, "bothsex_meta")) "free" else "icd3")
    spec <- build_gwas_model(apriori_model)
    std_lv <- isTRUE(config$gwas$std_lv %||% FALSE)
  }

  fix_measurement <- isTRUE(config$gwas$fix_measurement %||% TRUE)
  gc_mode <- config$gwas$genomic_control %||% "standard"
  q_snp_flag <- isTRUE(config$gwas$q_snp %||% TRUE)

  log_info("smoke", "Calling userGWAS on smoke subset...")
  t0 <- Sys.time()
  out <- GenomicSEM::userGWAS(
    covstruc = cov_for_gwas,
    SNPs = snp_test,
    estimation = config$cfa$estimator,
    model = spec$model,
    sub = spec$sub,
    cores = min(4L, as.integer(config$parallel$n_workers %||% 1L)),
    parallel = TRUE,
    smooth_check = TRUE,
    std.lv = std_lv,
    fix_measurement = fix_measurement,
    GC = gc_mode,
    Q_SNP = q_snp_flag,
    # printwarn MUST be TRUE: with Q_SNP + smooth_check, GenomicSEM's per-SNP
    # column names always include the error/warning columns, but those columns
    # are only appended when printwarn=TRUE. printwarn=FALSE yields a column
    # count 2 short of the name vector -> "'names' attribute [26] must be the
    # same length as the vector [24]" and the whole run aborts.
    printwarn = TRUE
  )
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  log_info("smoke", sprintf("userGWAS smoke complete in %.1fs", elapsed))

  if (length(spec$sub) == 1L && is.data.frame(out)) out <- list(out)
  first <- out[[1]]

  q_cols <- intersect(c("Q_SNP", "Q"), names(first))
  if (length(q_cols) == 0L) {
    log_fatal("smoke",
              sprintf("Q_SNP column missing from userGWAS output. Columns present: %s. Q_SNP=TRUE is set; check GenomicSEM version (need post-2024).",
                      paste(names(first), collapse = ", ")))
  }
  q_vals <- first[[q_cols[1]]]
  nonna_frac <- mean(!is.na(q_vals))
  min_frac <- as.numeric(config$smoke$min_qsnp_nonna_frac %||% 0.10)
  log_info("smoke", sprintf("Q_SNP column '%s': %.1f%% non-NA (threshold %.0f%%)",
                             q_cols[1], 100 * nonna_frac, 100 * min_frac))
  if (nonna_frac < min_frac) {
    log_fatal("smoke",
              sprintf("Q_SNP populated for only %.1f%% of test SNPs (need >=%.0f%%). Check model identifiability.",
                      100 * nonna_frac, 100 * min_frac))
  }

  log_info("smoke", "Smoke test PASSED")
  invisible(list(n_snps = nrow(snp_test),
                 q_nonna_frac = nonna_frac,
                 elapsed_sec = elapsed))
}
