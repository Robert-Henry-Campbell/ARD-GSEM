run_ldsc <- function(config, sex) {
  log_info("ldsc", sprintf("=== LDSC stage: %s ===", sex))

  munge_dir <- file.path(config$paths$output_dir, sex, "munge")
  munged_files <- list.files(munge_dir, pattern = "\\.sumstats\\.gz$", full.names = TRUE)

  if (length(munged_files) < 2) {
    log_fatal("ldsc", sprintf("Need >=2 munged traits for LDSC, found %d", length(munged_files)))
  }

  trait_names <- sub("_pre_munge\\.tsv\\.sumstats\\.gz$", "", basename(munged_files))
  trait_names <- sub("\\.sumstats\\.gz$", "", trait_names)
  log_info("ldsc", sprintf("Running LDSC for %s (%d traits: %s)",
                           sex, length(trait_names), paste(trait_names, collapse = ", ")))

  sample_prevs <- vapply(trait_names, function(trait) {
    cc <- get_case_control(config, sex, trait)
    cc$n_cases / (cc$n_cases + cc$n_controls)
  }, numeric(1))

  pop_prevs <- rep(NA, length(trait_names))

  ld_path <- config$paths$ld_scores
  wld_path <- config$paths$ld_scores

  ldsc_output <- GenomicSEM::ldsc(
    traits = munged_files,
    sample.prev = sample_prevs,
    population.prev = pop_prevs,
    ld = ld_path,
    wld = wld_path,
    trait.names = trait_names,
    stand = TRUE
  )

  S <- ldsc_output$S
  V <- ldsc_output$V
  I <- ldsc_output$I

  k <- length(trait_names)
  V_dim <- nrow(V)
  if (V_dim == k * k) {
    se_idx <- seq(1L, length.out = k, by = k + 1L)
    log_debug("ldsc", sprintf("V is vec-ordered (%dx%d); h2 SE at stride k+1", V_dim, V_dim))
  } else if (V_dim == k * (k + 1L) / 2L) {
    se_idx <- vech_diag_index(k)
    log_debug("ldsc", sprintf("V is vech-ordered (%dx%d); h2 SE at cumulative positions", V_dim, V_dim))
  } else {
    log_fatal("ldsc", sprintf("V has unexpected dimension %d for k=%d traits (expected %d=vec or %d=vech)",
                              V_dim, k, k * k, k * (k + 1L) / 2L))
  }
  stopifnot(length(se_idx) == k, all(se_idx <= V_dim))

  log_info("ldsc", "h2 estimates:")
  h2_table <- data.table(
    trait = trait_names,
    h2 = diag(S),
    se = sqrt(diag(V)[se_idx])
  )
  h2_table[, z := h2 / se]
  h2_table[, pass := z >= config$ldsc$h2_z_threshold]


  for (i in seq_len(nrow(h2_table))) {
    row <- h2_table[i]
    status <- ifelse(row$pass, "pass", "FAIL")
    log_info("ldsc", sprintf("  %s: h2=%.4f, SE=%.4f, Z=%.2f [%s]",
                             row$trait, row$h2, row$se, row$z, status))
    if (!row$pass) {
      log_warn("ldsc", sprintf("%s h2_z=%.2f < threshold %.1f; dropping from model",
                               row$trait, row$z, config$ldsc$h2_z_threshold))
    }
  }

  retained <- h2_table[pass == TRUE]$trait
  log_info("ldsc", sprintf("Traits retained: %d/%d (%s)",
                           length(retained), length(trait_names), paste(retained, collapse = ", ")))

  if (length(retained) < 2) {
    log_fatal("ldsc", "Fewer than 2 traits passed h2 filter; cannot proceed")
  }

  keep_idx <- which(trait_names %in% retained)
  S_filtered <- S[keep_idx, keep_idx, drop = FALSE]
  S_raw <- S_filtered

  psd_ok <- check_psd(S_filtered)
  smoothed <- FALSE
  if (psd_ok) {
    log_info("ldsc", sprintf("S matrix (%dx%d): all eigenvalues positive (PSD confirmed)",
                             nrow(S_filtered), ncol(S_filtered)))
  } else {
    log_warn("ldsc", "S matrix is NOT positive semi-definite; applying Matrix::nearPD smoothing")
    smoothed_result <- Matrix::nearPD(S_filtered, corr = FALSE, keepDiag = TRUE,
                                       ensureSymmetry = TRUE)
    S_filtered <- as.matrix(smoothed_result$mat)
    frob_delta <- sqrt(sum((S_filtered - S_raw)^2))
    log_info("ldsc", sprintf("nearPD smoothing: Frobenius-norm delta = %.4g; %d iterations",
                             frob_delta, smoothed_result$iterations))
    smoothed <- TRUE
  }

  if (!is.null(I)) {
    max_offdiag <- max(abs(I[upper.tri(I)]))
    log_info("ldsc", sprintf("V intercept matrix -- max |off-diagonal|: %.3f", max_offdiag))
    if (max_offdiag > 0.1) {
      log_warn("ldsc", "Large off-diagonal intercepts detected; check sample overlap assumptions")
    }
    intercept_path <- file.path(config$paths$output_dir, sex, "ldsc", "intercepts.csv")
    int_df <- as.data.frame(I)
    colnames(int_df) <- trait_names
    rownames(int_df) <- trait_names
    fwrite(int_df, intercept_path)
    log_debug("ldsc", sprintf("Full intercept matrix saved to %s", intercept_path))
  }

  if (smoothed) {
    ldsc_output$S[keep_idx, keep_idx] <- S_filtered
  }

  out_dir <- file.path(config$paths$output_dir, sex, "ldsc")
  saveRDS(ldsc_output, file.path(out_dir, "ldsc_full.rds"))
  saveRDS(S, file.path(out_dir, "S.rds"))
  saveRDS(V, file.path(out_dir, "V.rds"))
  if (smoothed) {
    saveRDS(S_raw, file.path(out_dir, "S_filtered_raw.rds"))
  }
  fwrite(h2_table, file.path(out_dir, "h2_qc.csv"))

  write_stage_manifest("ldsc", sex, config, list.files(out_dir, full.names = TRUE))

  list(
    S = S,
    V = V,
    I = I,
    ldsc_output = ldsc_output,
    h2_table = h2_table,
    retained_traits = retained,
    n_valid = length(retained),
    trait_names = trait_names
  )
}
