run_cfa <- function(config, sex, ldsc_results = NULL, efa_results = NULL) {
  log_info("cfa", sprintf("=== CFA stage: %s ===", sex))

  # Defensible split-sample path: if the EFA discovery structure + even-half covstruc +
  # pinned domains are present, run the held-out bake-off and lock the bifactor spec.
  # Otherwise fall through to the legacy EFA/a-priori CFA (used for sexes without a split).
  bakeoff <- run_cfa_bakeoff(config, sex)
  if (!is.null(bakeoff)) return(bakeoff)

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
    categories <- load_trait_categories(config, sex)
    code_format <- if (identical(sex, "bothsex_meta")) "free" else "icd3"
    apriori_model <- build_apriori_model(retained, categories,
                                         min_indicators = config$cfa$min_indicators_per_factor,
                                         code_format = code_format)

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

# --- Split-sample confirmatory bake-off + model lock (defensible ARD GSEM) ---
# Structure discovered on ODD (03_efa) is CONFIRMED on the held-out EVEN half, then the
# winning bifactor is re-fit on FULL and its GWAS spec saved. Returns NULL if the
# split-sample inputs are absent (caller falls back to the legacy CFA path).
# Skill basis: 03_interpretation (same-S comparison, CFI/SRMR/AIC + Δχ², anomaly checks).

extract_fit_row <- function(fit, key) {
  if (is.null(fit)) return(NULL)
  mf <- fit$modelfit
  r <- as.data.table(fit$results); setnames(r, tolower(names(r)))
  ec <- grep("std_genotype$|std_all$", names(r), value = TRUE)[1]
  sl <- if (!is.na(ec)) r[op == "=~" & grepl("^S[0-9]", lhs), as.numeric(get(ec))] else numeric(0)
  gg <- if (!is.na(ec)) r[op == "=~" & lhs == "ARD", as.numeric(get(ec))] else numeric(0)
  data.table(model = key, CFI = mf$CFI, SRMR = mf$SRMR, chisq = mf$chisq, df = mf$df, AIC = mf$AIC,
             spec_neg = sum(sl < 0, na.rm = TRUE), spec_tiny = sum(abs(sl) < 0.1, na.rm = TRUE),
             gen_load_mean = if (length(gg)) mean(gg, na.rm = TRUE) else NA_real_)
}

run_cfa_bakeoff <- function(config, sex) {
  ldsc_dir <- file.path(config$paths$output_dir, sex, "ldsc")
  out_dir  <- file.path(config$paths$output_dir, sex, "cfa")
  cand_path <- file.path(config$paths$output_dir, sex, "efa", "candidate_structures.rds")
  even_path <- file.path(ldsc_dir, "ldsc_even.rds")
  if (!file.exists(cand_path) || !file.exists(even_path) || is.null(config$cfa$domains)) return(NULL)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  cand <- readRDS(cand_path); broad <- cand$broad
  ld_even <- readRDS(even_path); ld_full <- readRDS(file.path(ldsc_dir, "ldsc_full.rds"))
  est <- config$cfa$estimator %||% "DWLS"
  fit_even <- function(model, label) tryCatch(
    suppressWarnings(GenomicSEM::usermodel(covstruc = ld_even, model = model, estimation = est, std.lv = TRUE)),
    error = function(e) { log_warn("cfa", sprintf("%s failed: %s", label, conditionMessage(e))); NULL })

  domains <- lapply(config$cfa$domains, function(v) intersect(as.character(v), broad))
  domains <- domains[vapply(domains, length, integer(1)) >= 2L]
  core <- unlist(domains, use.names = FALSE)
  log_info("cfa", sprintf("Split-sample confirm on EVEN half: general set=%d traits; %d domain(s): %s",
                          length(broad), length(domains),
                          paste(sprintf("%s(%d)", names(domains), lengths(domains)), collapse = ", ")))

  rvv <- function(v) c(sprintf("%s ~~ rvv_%s*%s", v, v, v), sprintf("rvv_%s > 0.001", v))
  umb_b <- paste(c(paste0("ARD =~ ", paste(broad, collapse = " + ")), rvv(broad)), collapse = "\n")
  bif   <- build_bifactor_model(broad, domains, "ARD", add_snp = FALSE, heywood = TRUE)
  umb_c <- paste(c(paste0("ARD =~ ", paste(core, collapse = " + ")), rvv(core)), collapse = "\n")
  corr  <- paste(c(vapply(seq_along(domains), function(i)
                   paste0("F", i, " =~ ", paste(domains[[i]], collapse = " + ")), character(1)), rvv(core)),
                 collapse = "\n")
  hier  <- build_hierarchical_model(setNames(domains, paste0("F", seq_along(domains))), "ARD", heywood = TRUE)

  tab <- rbindlist(list(
    extract_fit_row(fit_even(umb_b, "umbrella (broad)"),        "umbrella_broad"),
    extract_fit_row(fit_even(bif,   "bifactor (broad)"),        "bifactor_broad"),
    extract_fit_row(fit_even(umb_c, "umbrella (core)"),         "umbrella_core"),
    extract_fit_row(fit_even(corr,  "correlated (core)"),       "correlated_core"),
    extract_fit_row(fit_even(hier,  "hierarchical (core)"),     "hierarchical_core")), fill = TRUE)
  log_info("cfa", "Held-out (EVEN) model comparison:")
  for (i in seq_len(nrow(tab)))
    log_info("cfa", sprintf("  %-18s CFI=%.3f SRMR=%.3f AIC=%.0f | spec neg/tiny=%d/%d",
                            tab$model[i], tab$CFI[i], tab$SRMR[i], tab$AIC[i], tab$spec_neg[i], tab$spec_tiny[i]))
  ub <- tab[model == "umbrella_broad"]; bf <- tab[model == "bifactor_broad"]
  if (nrow(ub) && nrow(bf)) {
    dchi <- ub$chisq - bf$chisq; ddf <- ub$df - bf$df; pv <- pchisq(dchi, ddf, lower.tail = FALSE)
    log_info("cfa", sprintf("  Delta-chisq umbrella vs bifactor: dchi=%.1f ddf=%d p=%.2e -> structure %s the umbrella",
                            dchi, ddf, pv, if (pv < 0.05) "significantly beats" else "does NOT beat"))
  }
  fwrite(tab, file.path(out_dir, "fit_comparison.csv"))

  bif_full <- tryCatch(suppressWarnings(GenomicSEM::usermodel(covstruc = ld_full, model = bif, estimation = est, std.lv = TRUE)),
                       error = function(e) { log_error("cfa", sprintf("full bifactor fit failed: %s", conditionMessage(e))); NULL })
  if (!is.null(bif_full)) {
    fr <- extract_fit_row(bif_full, "bifactor_full")
    log_info("cfa", sprintf("Locked bifactor on FULL: CFI=%.3f SRMR=%.3f spec neg/tiny=%d/%d gen-load mean=%.2f",
                            fr$CFI, fr$SRMR, fr$spec_neg, fr$spec_tiny, fr$gen_load_mean))
    saveRDS(bif_full, file.path(out_dir, "model_bifactor_full.rds"))
    fwrite(fr, file.path(out_dir, "bifactor_full_fit.csv"))
  }
  factors <- c("ARD", names(domains))
  spec <- list(general_name = "ARD", general_traits = broad, domains = domains, factors = factors,
               measurement_model = bif,
               gwas_model = build_bifactor_model(broad, domains, "ARD", add_snp = TRUE, heywood = TRUE),
               sub = paste0(factors, " ~ SNP"))
  saveRDS(spec, file.path(out_dir, "model_spec.rds"))
  log_info("cfa", sprintf("Saved model_spec.rds -- bifactor GWAS over %d factor(s): %s",
                          length(factors), paste(factors, collapse = ", ")))
  write_stage_manifest("cfa", sex, config, list.files(out_dir, full.names = TRUE))
  invisible(list(fit_comparison = tab, spec = spec))
}
