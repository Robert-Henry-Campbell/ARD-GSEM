# EFA stage — split-sample structure DISCOVERY.
#
# Defensible-analysis design (GSEM skill non-negotiable #5: "EFA on an independent
# split before committing to a CFA"): when config$ldsc$chr_split is set, the EFA runs
# on the ODD-chromosome half (ldsc_odd.rds); 04_cfa then confirms the structure on the
# EVEN half. Factor number is chosen by CONVERGENT criteria (scree, parallel analysis,
# Velicer MAP, VSS) rather than a single rule, and indicator assignment is pre-specified
# (dominant loading >= specific_loading_min, non-cross-loading), never forced. A
# per-factor genetic-correlation coherence metric flags grab-bag factors. Outputs the
# candidate structures (per k) + the broad general-factor indicator set for the bake-off.

run_efa <- function(config, sex, ldsc_results = NULL) {
  log_info("efa", sprintf("=== EFA stage: %s ===", sex))
  ldsc_dir <- file.path(config$paths$output_dir, sex, "ldsc")
  out_dir  <- file.path(config$paths$output_dir, sex, "efa")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  # Discovery half: odd chromosomes when chr_split, else the full-sample S.
  odd_path <- file.path(ldsc_dir, "ldsc_odd.rds")
  use_split <- isTRUE(config$ldsc$chr_split) && file.exists(odd_path)
  ldsc_output <- if (use_split) readRDS(odd_path) else readRDS(file.path(ldsc_dir, "ldsc_full.rds"))
  log_info("efa", sprintf("Discovery covstruc: %s", if (use_split) "ODD-chromosome half (split-sample)" else "full sample"))

  h2 <- fread(file.path(ldsc_dir, "h2_qc.csv"))
  retained <- h2[pass == TRUE]$trait
  trait_names <- colnames(ldsc_output$S)
  dimnames(ldsc_output$S) <- list(trait_names, trait_names)   # ldsc() sets colnames only
  retained <- intersect(retained, trait_names)
  S <- ldsc_output$S[retained, retained, drop = FALSE]
  rg <- cov2cor(S); dimnames(rg) <- list(retained, retained)
  rg_pd <- as.matrix(Matrix::nearPD(rg, corr = TRUE)$mat); dimnames(rg_pd) <- list(retained, retained)
  gmean <- mean(rg[upper.tri(rg)])

  neffs <- vapply(retained, function(tr) {
    cc <- get_case_control(config, sex, tr); compute_neff(cc$n_cases, cc$n_controls)
  }, numeric(1))
  n_obs <- harmonic_neff(neffs); if (is.na(n_obs)) n_obs <- 1000
  log_info("efa", sprintf("n.obs = harmonic Neff = %.0f | global mean rg = %.3f", n_obs, gmean))

  # --- Convergent factor-number criteria ---
  ev <- eigen(rg_pd, symmetric = TRUE, only.values = TRUE)$values
  pa <- tryCatch(psych::fa.parallel(rg_pd, n.obs = n_obs, fa = "fa", plot = FALSE, fm = "minres")$nfact,
                 error = function(e) NA_integer_)
  nf <- tryCatch({ invisible(capture.output(
           r <- psych::nfactors(rg_pd, n = min(10L, nrow(rg_pd) - 1L), n.obs = n_obs,
                                fm = "minres", rotate = "promax"))); r },
                 error = function(e) NULL)
  map_k <- if (!is.null(nf)) which.min(nf$map) else NA_integer_
  vss1_k <- if (!is.null(nf)) which.max(nf$vss.stats$cfit.1) else NA_integer_
  bic_k <- if (!is.null(nf)) which.min(nf$vss.stats$BIC) else NA_integer_
  log_info("efa", sprintf("k-criteria: scree(top5)=%s Kaiser>1=%d | parallel=%s MAP=%s VSS1=%s BIC=%s",
                          paste(round(ev[1:5], 2), collapse = ","), sum(ev > 1),
                          pa, map_k, vss1_k, bic_k))
  fwrite(data.table(criterion = c("parallel", "MAP", "VSS1", "BIC", "Kaiser"),
                    k = c(pa, map_k, vss1_k, bic_k, sum(ev > 1))),
         file.path(out_dir, "n_factors_criteria.csv"))
  fwrite(data.table(component = seq_along(ev), eigenvalue = ev), file.path(out_dir, "scree.csv"))

  # --- Broad general-factor indicator set (1-factor EFA loading >= general_loading_min) ---
  gmin <- config$efa$general_loading_min %||% 0.30
  gload <- setNames(as.numeric(psych::fa(rg_pd, nfactors = 1, rotate = "none", fm = "ml",
                                         n.obs = n_obs)$loadings[, 1]), retained)
  broad <- names(gload)[abs(gload) >= gmin]
  log_info("efa", sprintf("Broad general-factor set (|loading|>=%.2f): %d/%d traits; excluded: %s",
                          gmin, length(broad), length(retained), paste(setdiff(retained, broad), collapse = ",")))

  # --- Candidate domain structures across k, with indicator discipline + coherence ---
  smin <- config$efa$specific_loading_min %||% 0.40
  cross <- config$efa$cross_load_ratio %||% 1.30
  min_ind <- config$cfa$min_indicators_per_factor %||% 3L
  k_range <- config$efa$candidate_k %||% 3:6
  assign_struct <- function(k) {
    fa <- tryCatch(psych::fa(rg_pd, nfactors = k, rotate = "promax", fm = "ml", n.obs = n_obs),
                   error = function(e) psych::fa(rg_pd, nfactors = k, rotate = "promax", fm = "minres", n.obs = n_obs))
    L <- as.matrix(fa$loadings[, , drop = FALSE]); a <- rep(NA_integer_, nrow(L))
    for (i in seq_len(nrow(L))) {
      o <- order(abs(L[i, ]), decreasing = TRUE)
      if (abs(L[i, o[1]]) >= smin && abs(L[i, o[1]]) >= cross * abs(L[i, o[2]])) a[i] <- o[1]
    }
    g <- split(retained, a); g[vapply(g, length, integer(1)) >= min_ind]
  }
  structs <- list()
  for (k in k_range) {
    g <- assign_struct(k)
    if (length(g)) names(g) <- paste0("F", seq_along(g))
    log_info("efa", sprintf("k=%d -> %d coherent factor(s), %d/%d traits assigned",
                            k, length(g), length(unlist(g)), length(retained)))
    for (nm in names(g)) {
      co <- factor_coherence(rg, g[[nm]])
      log_info("efa", sprintf("   %s (n=%d within-rg mean=%.2f min=%.2f%s): %s", nm, length(g[[nm]]),
                              co$mean, co$min, if (co$mean < gmean + 0.15) " [INCOHERENT]" else "",
                              paste(g[[nm]], collapse = ",")))
    }
    structs[[as.character(k)]] <- g
  }

  saveRDS(list(structs = structs, broad = broad, retained = retained, n_obs = n_obs,
               gmean = gmean, gload = gload, half = if (use_split) "odd" else "full"),
          file.path(out_dir, "candidate_structures.rds"))
  write_stage_manifest("efa", sex, config, list.files(out_dir, full.names = TRUE))
  log_info("efa", "EFA discovery complete; candidate_structures.rds written")
  invisible(list(structs = structs, broad = broad, retained = retained, n_obs = n_obs))
}
