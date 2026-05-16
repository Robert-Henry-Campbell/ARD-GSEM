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

  h2 <- fread(file.path(config$paths$output_dir, sex, "ldsc", "h2_qc.csv"))
  retained <- h2[pass == TRUE]$trait
  categories <- fread(file.path(config$paths$meta_dir, "icd10_categories.csv"))
  loadings_path <- file.path(config$paths$output_dir, sex, "cfa", "factor_loadings.csv")
  ordering_table <- if (file.exists(loadings_path)) fread(loadings_path) else NULL
  apriori_model <- build_apriori_model(
    retained, categories,
    min_indicators = config$cfa$min_indicators_per_factor,
    ordering_table = ordering_table)
  spec <- build_gwas_model(apriori_model)

  fix_measurement <- isTRUE(config$gwas$fix_measurement %||% TRUE)
  gc_mode <- config$gwas$genomic_control %||% "standard"
  std_lv <- isTRUE(config$gwas$std_lv %||% FALSE)
  q_snp_flag <- isTRUE(config$gwas$q_snp %||% TRUE)

  log_info("smoke", "Calling userGWAS on smoke subset...")
  t0 <- Sys.time()
  out <- GenomicSEM::userGWAS(
    covstruc = ldsc_output,
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
    printwarn = FALSE
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
