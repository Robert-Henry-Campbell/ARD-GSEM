parse_factor_names <- function(model_str) {
  if (nchar(model_str) == 0L) return(character(0))
  lines <- strsplit(model_str, "\n", fixed = TRUE)[[1]]
  lines <- lines[grepl("=~", lines, fixed = TRUE)]
  trimws(sub("\\s*=~.*$", "", lines))
}

build_gwas_model <- function(apriori_model) {
  factor_names <- parse_factor_names(apriori_model)
  if (length(factor_names) == 0L) {
    stop("build_gwas_model: no factors found in a priori model")
  }
  snp_regressions <- paste0(factor_names, " ~ SNP")
  gwas_model <- paste(c(apriori_model, snp_regressions), collapse = "\n")
  list(model = gwas_model, sub = snp_regressions, factors = factor_names)
}

run_gwas <- function(config, sex) {
  log_info("gwas", sprintf("=== GWAS stage: %s ===", sex))

  ldsc_path <- file.path(config$paths$output_dir, sex, "ldsc", "ldsc_full.rds")
  ldsc_output <- readRDS(ldsc_path)

  snp_path <- file.path(config$paths$output_dir, sex, "sumstats",
                        paste0(sex, "_snp_sumstats.rds"))
  if (!file.exists(snp_path)) {
    legacy <- file.path(config$paths$output_dir, sex, "sumstats", "snp_sumstats.rds")
    if (file.exists(legacy)) {
      log_warn("gwas", sprintf(
        "Using legacy (unprefixed) snp_sumstats at %s (mtime=%s) -- expected %s",
        legacy, format(file.mtime(legacy)), snp_path))
      snp_path <- legacy
    } else {
      log_fatal("gwas", sprintf(
        "snp_sumstats not found at %s (run --stage sumstats first)", snp_path))
    }
  }
  snp_sumstats <- readRDS(snp_path)

  h2_path <- file.path(config$paths$output_dir, sex, "ldsc", "h2_qc.csv")
  retained <- fread(h2_path)[pass == TRUE]$trait

  categories <- fread(file.path(config$paths$meta_dir, "icd10_categories.csv"))
  loadings_path <- file.path(config$paths$output_dir, sex, "cfa", "factor_loadings.csv")
  ordering_table <- if (file.exists(loadings_path)) {
    log_info("gwas", sprintf("Indicator ordering: descending |std_loading| from %s", loadings_path))
    fread(loadings_path)
  } else {
    log_warn("gwas", sprintf("factor_loadings.csv not found at %s; falling back to alphabetical ordering within each chapter",
                              loadings_path))
    NULL
  }
  apriori_model <- build_apriori_model(
    retained, categories,
    min_indicators = config$cfa$min_indicators_per_factor,
    ordering_table = ordering_table,
    add_factor_covariances = TRUE,
    add_heywood_constraints = TRUE)
  if (nchar(apriori_model) == 0L) {
    log_fatal("gwas", "A priori model is empty; cannot run factor GWAS")
  }

  spec <- build_gwas_model(apriori_model)
  log_info("gwas", sprintf("Factor GWAS over %d chapter factor(s): %s",
                            length(spec$factors), paste(spec$factors, collapse = ", ")))
  log_info("gwas", sprintf("Augmented model:\n%s", spec$model))

  cores <- as.integer(config$parallel$n_workers %||% 1L)
  smooth_check <- isTRUE(config$gwas$smooth_check)

  log_info("gwas", sprintf("Calling GenomicSEM::userGWAS (cores=%d, smooth_check=%s)...",
                            cores, smooth_check))
  t0 <- Sys.time()

  std_lv <- isTRUE(config$gwas$std_lv %||% FALSE)
  log_info("gwas", sprintf("userGWAS std.lv=%s (%s)",
                            std_lv,
                            if (std_lv) "factor variance fixed to 1; SNP betas are per-SD-of-factor"
                            else "marker-variable identification; first indicator loading fixed to 1; SNP betas in units of the first indicator"))

  fix_measurement <- isTRUE(config$gwas$fix_measurement %||% TRUE)
  gc_mode <- config$gwas$genomic_control %||% "standard"
  q_snp_flag <- isTRUE(config$gwas$q_snp %||% TRUE)

  log_info("gwas", sprintf(
    "userGWAS args: Q_SNP=%s, fix_measurement=%s, GC=%s, smooth_check=%s, std.lv=%s",
    q_snp_flag, fix_measurement, gc_mode, smooth_check, std_lv))

  gwas_result <- GenomicSEM::userGWAS(
    covstruc = ldsc_output,
    SNPs = snp_sumstats,
    estimation = config$cfa$estimator,
    model = spec$model,
    sub = spec$sub,
    cores = cores,
    parallel = TRUE,
    smooth_check = smooth_check,
    std.lv = std_lv,
    fix_measurement = fix_measurement,
    GC = gc_mode,
    Q_SNP = q_snp_flag,
    printwarn = FALSE
  )

  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  log_info("gwas", sprintf("userGWAS complete in %.1f min", elapsed))

  # userGWAS returns a bare data.frame (not a length-1 list) when sub has 1 entry.
  if (length(spec$sub) == 1L && is.data.frame(gwas_result)) {
    gwas_result <- list(gwas_result)
  }
  if (!is.list(gwas_result) || length(gwas_result) != length(spec$sub)) {
    log_fatal("gwas", sprintf(
      "userGWAS returned %d elements; expected %d (one per sub regression)",
      if (is.list(gwas_result)) length(gwas_result) else 1L, length(spec$sub)))
  }
  names(gwas_result) <- spec$factors

  out_dir <- file.path(config$paths$output_dir, sex, "gwas")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(gwas_result, file.path(out_dir, paste0(sex, "_userGWAS_raw.rds")))
  saveRDS(spec, file.path(out_dir, paste0(sex, "_model_spec.rds")))

  for (f in names(gwas_result)) {
    g <- gwas_result[[f]]
    n_rows <- if (is.data.frame(g)) nrow(g) else NA_integer_
    log_info("gwas", sprintf("  factor %s: %s SNPs", f, format(n_rows, big.mark = ",")))
  }

  write_stage_manifest("gwas", sex, config, list.files(out_dir, full.names = TRUE))

  list(gwas_result = gwas_result, factors = spec$factors, model = spec$model)
}
