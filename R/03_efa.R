run_efa <- function(config, sex, ldsc_results = NULL) {
  log_info("efa", sprintf("=== EFA stage: %s ===", sex))

  if (is.null(ldsc_results)) {
    ldsc_path <- file.path(config$paths$output_dir, sex, "ldsc", "ldsc_full.rds")
    h2_path <- file.path(config$paths$output_dir, sex, "ldsc", "h2_qc.csv")
    ldsc_output <- readRDS(ldsc_path)
    h2_table <- fread(h2_path)
    retained <- h2_table[pass == TRUE]$trait
    S <- ldsc_output$S
    trait_names <- h2_table$trait
  } else {
    S <- ldsc_results$S
    retained <- ldsc_results$retained_traits
    trait_names <- ldsc_results$trait_names
  }

  keep_idx <- which(trait_names %in% retained)
  S_sub <- S[keep_idx, keep_idx, drop = FALSE]
  colnames(S_sub) <- retained
  rownames(S_sub) <- retained

  diag_sqrt <- sqrt(diag(S_sub))
  rg <- S_sub / outer(diag_sqrt, diag_sqrt)
  diag(rg) <- 1

  log_info("efa", sprintf("Parallel analysis on %dx%d genetic correlation matrix (%s)",
                          nrow(rg), ncol(rg), sex))

  neffs <- vapply(retained, function(trait) {
    cc <- get_case_control(config, sex, trait)
    compute_neff(cc$n_cases, cc$n_controls)
  }, numeric(1))
  n_obs <- harmonic_neff(neffs)
  if (is.na(n_obs)) {
    log_warn("efa", "harmonic Neff is NA; falling back to n.obs=1000 placeholder")
    n_obs <- 1000
  }
  log_info("efa", sprintf("EFA n.obs = harmonic mean of per-trait Neff = %.0f", n_obs))

  n_factors <- 1
  if (nrow(rg) >= 4 && config$efa$use_parallel_analysis) {
    if (!is.null(ldsc_results) && !is.null(ldsc_results$ldsc_output)) {
      pa_covstruc <- ldsc_results$ldsc_output
    } else {
      pa_covstruc <- readRDS(file.path(config$paths$output_dir, sex, "ldsc", "ldsc_full.rds"))
    }
    n_factors <- tryCatch({
      # Signature unknown in advance; try defaults first. Confirm on doraemon6
      # via ?GenomicSEM::paLDSC and adjust args if the function takes named iters.
      pa <- GenomicSEM::paLDSC(covstruc = pa_covstruc)
      pa_n <- if (is.list(pa) && !is.null(pa$nfact)) pa$nfact
              else if (is.list(pa) && !is.null(pa$n_factors)) pa$n_factors
              else if (is.numeric(pa) && length(pa) == 1L) as.integer(pa)
              else NA_integer_
      if (is.na(pa_n) || pa_n < 1L) {
        stop(sprintf("paLDSC returned no usable factor count (class=%s)", class(pa)[1]))
      }
      log_info("efa", sprintf("paLDSC suggested %d factor(s)", as.integer(pa_n)))
      max(1L, as.integer(pa_n))
    }, error = function(e) {
      log_warn("efa", sprintf(
        "GenomicSEM::paLDSC unavailable/failed (%s); falling back to psych::fa.parallel on rg",
        conditionMessage(e)))
      pa <- psych::fa.parallel(rg, n.obs = n_obs, fa = "fa", fm = "ml",
                               plot = FALSE, n.iter = 20)
      max(1L, as.integer(pa$nfact))
    })
  } else if (nrow(rg) == 3) {
    n_factors <- 1
  }

  n_factors <- min(n_factors, floor(nrow(rg) / 2))
  log_info("efa", sprintf("Suggested factors: %d", n_factors))

  rotation <- config$efa$rotation
  if (n_factors == 1) rotation <- "none"

  fa_result <- psych::fa(rg, nfactors = n_factors, rotate = rotation,
                         fm = "ml", n.obs = n_obs)

  loadings_matrix <- matrix(fa_result$loadings[], nrow = nrow(rg), ncol = n_factors)
  rownames(loadings_matrix) <- retained
  colnames(loadings_matrix) <- paste0("F", seq_len(n_factors))

  for (f in seq_len(n_factors)) {
    high_loading <- retained[abs(loadings_matrix[, f]) >= config$efa$loading_threshold]
    loading_vals <- round(loadings_matrix[, f], 3)
    names(loading_vals) <- retained
    log_info("efa", sprintf("Factor %d loadings: %s",
                            f, paste(sprintf("%s=%.3f", retained, loading_vals), collapse = ", ")))
  }

  out_dir <- file.path(config$paths$output_dir, sex, "efa")
  saveRDS(fa_result, file.path(out_dir, "fa_result.rds"))
  saveRDS(loadings_matrix, file.path(out_dir, "loadings.rds"))
  writeLines(as.character(n_factors), file.path(out_dir, "n_factors.txt"))

  loadings_df <- as.data.frame(loadings_matrix)
  loadings_df$trait <- rownames(loadings_matrix)
  fwrite(loadings_df, file.path(out_dir, "loadings.csv"))

  write_stage_manifest("efa", sex, config, list.files(out_dir, full.names = TRUE))

  list(
    fa_result = fa_result,
    loadings = loadings_matrix,
    n_factors = n_factors,
    rg = rg,
    retained = retained
  )
}
