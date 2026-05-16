run_cfa <- function(config, sex, ldsc_results = NULL, efa_results = NULL) {
  log_info("cfa", sprintf("=== CFA stage: %s ===", sex))

  if (is.null(ldsc_results)) {
    ldsc_path <- file.path(config$paths$output_dir, sex, "ldsc", "ldsc_full.rds")
    ldsc_output <- readRDS(ldsc_path)
    h2_table <- fread(file.path(config$paths$output_dir, sex, "ldsc", "h2_qc.csv"))
    retained <- h2_table[pass == TRUE]$trait
    trait_names <- h2_table$trait
  } else {
    ldsc_output <- ldsc_results$ldsc_output
    retained <- ldsc_results$retained_traits
    trait_names <- ldsc_results$trait_names
  }

  if (is.null(efa_results)) {
    loadings_path <- file.path(config$paths$output_dir, sex, "efa", "loadings.rds")
    if (file.exists(loadings_path)) {
      efa_loadings <- readRDS(loadings_path)
    } else {
      efa_loadings <- NULL
    }
  } else {
    efa_loadings <- efa_results$loadings
  }

  results <- list(converged = FALSE)
  out_dir <- file.path(config$paths$output_dir, sex, "cfa")

  # EFA-derived model
  if (!is.null(efa_loadings)) {
    efa_model <- efa_to_lavaan(efa_loadings, threshold = config$efa$loading_threshold)
    if (nchar(efa_model) == 0) {
      log_warn("cfa", sprintf(
        "EFA-derived model is empty: no indicator passed loading_threshold=%.2f; skipping EFA CFA fit",
        config$efa$loading_threshold))
    }
    if (nchar(efa_model) > 0) {
      log_info("cfa", sprintf("Fitting EFA-derived model (%d factor(s), %d indicators)...",
                              ncol(efa_loadings), sum(apply(abs(efa_loadings) >= config$efa$loading_threshold, 2, sum))))
      tryCatch({
        efa_fit <- GenomicSEM::usermodel(
          covstruc = ldsc_output,
          model = efa_model,
          estimation = config$cfa$estimator
        )
        results$efa_fit <- efa_fit
        results$converged <- TRUE

        fit_indices <- extract_fit(efa_fit)
        log_info("cfa", sprintf("EFA model: converged=TRUE | CFI=%.3f | RMSEA=%.3f | df=%d",
                                fit_indices$cfi, fit_indices$rmsea, fit_indices$df))
        if (fit_indices$df == 0) {
          log_warn("cfa", "Model is just-identified (0 df); fit indices are trivially perfect")
        }

        saveRDS(efa_fit, file.path(out_dir, "model_efa.rds"))
      }, error = function(e) {
        log_error("cfa", sprintf("EFA model failed: %s", e$message))
      })
    }
  }

  # A priori model
  if (config$cfa$model_source %in% c("both", "apriori")) {
    log_info("cfa", "Building a priori model from ICD-10 chapters...")
    categories <- fread(file.path(config$paths$meta_dir, "icd10_categories.csv"))
    apriori_model <- build_apriori_model(retained, categories,
                                         min_indicators = config$cfa$min_indicators_per_factor)

    if (nchar(apriori_model) > 0) {
      log_info("cfa", sprintf("A priori model:\n%s", apriori_model))
      tryCatch({
        apriori_fit <- GenomicSEM::usermodel(
          covstruc = ldsc_output,
          model = apriori_model,
          estimation = config$cfa$estimator
        )
        results$apriori_fit <- apriori_fit
        if (!results$converged) results$converged <- TRUE

        fit_indices_ap <- extract_fit(apriori_fit)
        log_info("cfa", sprintf("A priori model: converged=TRUE | CFI=%.3f | RMSEA=%.3f | df=%d",
                                fit_indices_ap$cfi, fit_indices_ap$rmsea, fit_indices_ap$df))

        saveRDS(apriori_fit, file.path(out_dir, "model_apriori.rds"))
      }, error = function(e) {
        log_error("cfa", sprintf("A priori model failed: %s", e$message))
      })
    } else {
      log_warn("cfa", "A priori model: no factors with >=2 indicators; skipping")
    }
  }

  # Save fit comparison
  fit_comparison <- data.table(model = character(), cfi = numeric(),
                               rmsea = numeric(), srmr = numeric(),
                               chisq = numeric(), df = integer())
  if (!is.null(results$efa_fit)) {
    f <- extract_fit(results$efa_fit)
    fit_comparison <- rbind(fit_comparison, data.table(
      model = "efa", cfi = f$cfi, rmsea = f$rmsea, srmr = f$srmr,
      chisq = f$chisq, df = f$df))
  }
  if (!is.null(results$apriori_fit)) {
    f <- extract_fit(results$apriori_fit)
    fit_comparison <- rbind(fit_comparison, data.table(
      model = "apriori", cfi = f$cfi, rmsea = f$rmsea, srmr = f$srmr,
      chisq = f$chisq, df = f$df))
  }
  fwrite(fit_comparison, file.path(out_dir, "fit_comparison.csv"))

  # Save standardized loadings
  if (!is.null(results$efa_fit)) {
    sol <- results$efa_fit$results
    if (!is.null(sol)) {
      fwrite(as.data.table(sol), file.path(out_dir, "loadings.csv"))
      log_info("cfa", "Standardized EFA loadings saved")
    }
  }
  if (!is.null(results$apriori_fit)) {
    sol_ap <- results$apriori_fit$results
    if (!is.null(sol_ap)) {
      fwrite(as.data.table(sol_ap), file.path(out_dir, "loadings_apriori.csv"))
      log_info("cfa", "Standardized a priori (chapter) loadings saved")
      factor_loadings <- extract_apriori_loadings(results$apriori_fit, categories)
      if (!is.null(factor_loadings) && nrow(factor_loadings) > 0L) {
        fwrite(factor_loadings, file.path(out_dir, "factor_loadings.csv"))
        log_info("cfa", sprintf(
          "factor_loadings.csv written (%d indicator rows); will be used by R/08_gwas.R for descending-loading indicator ordering",
          nrow(factor_loadings)))
      } else {
        log_warn("cfa", "extract_apriori_loadings returned empty; userGWAS will fall back to alphabetical ordering")
      }
    }
  }

  results$fit <- if (!is.null(results$efa_fit)) extract_fit(results$efa_fit) else
                 if (!is.null(results$apriori_fit)) extract_fit(results$apriori_fit) else
                 list(cfi = NA, rmsea = NA, srmr = NA, chisq = NA, df = NA)

  write_stage_manifest("cfa", sex, config, list.files(out_dir, full.names = TRUE))
  results
}

extract_fit <- function(model_result) {
  if (is.null(model_result$modelfit)) {
    return(list(cfi = NA, rmsea = NA, srmr = NA, chisq = NA, df = NA))
  }
  mf <- model_result$modelfit
  pick <- function(...) {
    for (nm in c(...)) {
      if (nm %in% colnames(mf)) return(mf[1, nm])
    }
    NA
  }
  list(
    cfi = as.numeric(pick("CFI", "cfi")),
    rmsea = as.numeric(pick("RMSEA", "rmsea")),
    srmr = as.numeric(pick("SRMR", "srmr")),
    chisq = as.numeric(pick("chisq", "Chisq", "CHISQ")),
    df = as.integer(pick("df", "DF"))
  )
}
