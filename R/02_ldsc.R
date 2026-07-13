run_ldsc <- function(config, sex) {
  log_info("ldsc", sprintf("=== LDSC stage: %s ===", sex))

  munge_dir <- file.path(config$paths$output_dir, sex, "munge")
  munged_files <- list.files(munge_dir, pattern = "\\.sumstats\\.gz$", full.names = TRUE)

  if (length(munged_files) < 2) {
    log_fatal("ldsc", sprintf("Need >=2 munged traits for LDSC, found %d", length(munged_files)))
  }

  trait_names <- basename(munged_files)
  trait_names <- sub("_pre_munge\\.tsv\\.sumstats\\.gz$", "", trait_names)
  trait_names <- sub("\\.sumstats\\.gz$", "", trait_names)
  trait_names <- sub(paste0("^", sex, "_"), "", trait_names)
  log_info("ldsc", sprintf("Running LDSC for %s (%d traits: %s)",
                           sex, length(trait_names), paste(trait_names, collapse = ", ")))

  sample_prevs <- rep(0.5, length(trait_names))

  pop_prev_source <- config$ldsc$pop_prev_source %||% "ukb_cohort"
  if (isTRUE(config$ldsc$liability_scale) && identical(pop_prev_source, "ukb_cohort")) {
    pop_prevs <- vapply(trait_names, function(trait) {
      cc <- get_case_control(config, sex, trait)
      cc$n_cases / (cc$n_cases + cc$n_controls)
    }, numeric(1))
    log_info("ldsc", sprintf(
      "population.prev sourced from UKB cohort case fraction (proxy for K_pop; see methods caveat)"))
  } else {
    pop_prev_csv <- file.path(config$paths$meta_dir,
                              sprintf("ukb_pop_prev_%s.csv", sex))
    if (file.exists(pop_prev_csv)) {
      kp <- fread(pop_prev_csv)
      pop_prevs <- kp$K_pop[match(trait_names, kp$code)]
      log_info("ldsc", sprintf("population.prev loaded from %s", basename(pop_prev_csv)))
    } else {
      pop_prevs <- rep(NA_real_, length(trait_names))
      log_warn("ldsc", "population.prev=NA for all traits (observed-scale h2)")
    }
  }

  ld_path <- config$paths$ld_scores
  wld_path <- config$paths$ld_scores

  # GenomicSEM::ldsc() also writes a "*_ldsc.log" to CWD (relative path, same trap as
  # munge/sumstats). Force CWD to the per-sex ldsc dir for the duration of the call.
  out_dir <- file.path(config$paths$output_dir, sex, "ldsc")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  munged_files_abs <- normalizePath(munged_files, mustWork = TRUE)

  # --- Published-h2-Z pre-filter (before ldsc; keeps the block-jackknife valid) ---
  # GenomicSEM inflates n.blocks to ~k(k+1)/2 when n.traits>18, so >~43 traits -> tiny
  # blocks -> corrupted SEs. Select on the PUBLISHED h2-Z (Pan-UKB EUR / Neale both_sexes),
  # floored at published_h2_z_min and capped to the block budget. The existing post-ldsc
  # own-z gate + Frobenius loop below then do final QC on this valid set.
  h2_source  <- config$ldsc$h2_z_source %||% "published"
  max_blocks <- config$ldsc$max_jackknife_blocks %||% 1000L
  k_max      <- max_traits_for_blocks(max_blocks)
  if (identical(h2_source, "published")) {
    z_min <- config$ldsc$published_h2_z_min %||% 4.0
    pub   <- published_h2z(config, sex, trait_names)
    pubz  <- pub$pub_h2_z[match(trait_names, pub$trait)]
    if (all(is.na(pubz))) {
      log_warn("ldsc", sprintf(
        "No published h2-Z source for sex=%s; skipping published pre-filter (%d traits)", sex, length(trait_names)))
      if (length(trait_names) > k_max) {
        log_warn("ldsc", sprintf(
          "%d traits > k_max=%d -> GenomicSEM will set n.blocks=%d (>%d); SEs may be biased",
          length(trait_names), k_max, n_jackknife_blocks(length(trait_names)), max_blocks))
      }
    } else {
      keep   <- which(!is.na(pubz) & pubz >= z_min)
      capped <- length(keep) > k_max
      if (capped) keep <- sort(keep[order(pubz[keep], decreasing = TRUE)][seq_len(k_max)])
      fwrite(data.table(trait = trait_names, pub_h2_z = pubz,
                        pub_h2_confidence = pub$pub_h2_confidence[match(trait_names, pub$trait)],
                        selected = seq_along(trait_names) %in% keep),
             file.path(out_dir, "h2_screen.csv"))
      src_lbl <- if (identical(sex, "bothsex")) "Pan-UKB EUR" else "Neale both_sexes"
      if (length(keep) < 2L) {
        log_fatal("ldsc", sprintf(
          "Published-h2-Z filter (%s, z>=%.2f) kept %d trait(s); need >=2. See %s.",
          src_lbl, z_min, length(keep), file.path(out_dir, "h2_screen.csv")))
      }
      log_info("ldsc", sprintf(
        "Published-h2-Z pre-filter (%s): %d -> %d traits (z>=%.2f, effective z=%.2f%s; n.blocks=%d, budget=%d)",
        src_lbl, length(trait_names), length(keep), z_min, min(pubz[keep]),
        if (capped) sprintf(", capped to k_max=%d", k_max) else "",
        n_jackknife_blocks(length(keep)), max_blocks))
      munged_files_abs <- munged_files_abs[keep]
      trait_names      <- trait_names[keep]
      sample_prevs     <- sample_prevs[keep]
      pop_prevs        <- pop_prevs[keep]
    }
  } else {
    log_info("ldsc", sprintf("h2_z_source=computed: no published-z pre-filter (%d traits)", length(trait_names)))
    if (length(trait_names) > k_max) {
      log_warn("ldsc", sprintf(
        "computed mode: %d traits > k_max=%d -> GenomicSEM n.blocks=%d (>%d); SEs may be biased",
        length(trait_names), k_max, n_jackknife_blocks(length(trait_names)), max_blocks))
    }
  }

  saved_wd <- getwd()
  # Attempt stand=TRUE for standardized loadings (S_Stand / V_Stand). On the V_Stand
  # crash (occurs when any h2 estimate is near-zero / negative), fall back to FALSE.
  call_ldsc <- function(stand_flag, select = FALSE) {
    GenomicSEM::ldsc(
      traits = munged_files_abs,
      sample.prev = sample_prevs,
      population.prev = pop_prevs,
      ld = ld_path,
      wld = wld_path,
      trait.names = trait_names,
      stand = stand_flag,
      select = select
    )
  }
  ldsc_output <- tryCatch({
    setwd(out_dir)
    log_info("ldsc", "Attempting GenomicSEM::ldsc(stand=TRUE) ...")
    out <- tryCatch(call_ldsc(TRUE), error = function(e) {
      log_warn("ldsc", sprintf(
        "stand=TRUE failed (%s); falling back to stand=FALSE",
        conditionMessage(e)))
      call_ldsc(FALSE)
    })
    if (!is.null(out$S_Stand) && !is.null(out$V_Stand)) {
      log_info("ldsc", "stand=TRUE succeeded; S_Stand + V_Stand available")
    } else {
      log_warn("ldsc", "stand path completed without S_Stand/V_Stand; loadings reported on raw scale")
    }
    out
  }, finally = setwd(saved_wd))

  # --- Split-sample halves for defensible EFA/CFA (structure from odd, confirmed on
  # even; see docs plan). GenomicSEM::ldsc(select=) restricts LD scores + weights + M
  # to odd/even chromosomes, so each half is an independent-genome covstruc over the
  # SAME trait set. Saved beside ldsc_full.rds for 03_efa (odd) and 04_cfa (even). ---
  if (isTRUE(config$ldsc$chr_split)) {
    for (half in c("ODD", "EVEN")) {
      log_info("ldsc", sprintf("chr_split: %s-chromosome LDSC (%d traits) for split-sample EFA/CFA...",
                               half, length(trait_names)))
      half_out <- tryCatch({
        setwd(out_dir)
        tryCatch(call_ldsc(TRUE, half), error = function(e) {
          log_warn("ldsc", sprintf("chr_split %s stand=TRUE failed (%s); stand=FALSE",
                                    half, conditionMessage(e)))
          call_ldsc(FALSE, half)
        })
      }, finally = setwd(saved_wd))
      saveRDS(half_out, file.path(out_dir, sprintf("ldsc_%s.rds", tolower(half))))
      log_info("ldsc", sprintf("chr_split: saved ldsc_%s.rds", tolower(half)))
    }
  }

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

  halt_thresh <- config$munge$smoothing_frob_halt %||% 0.025
  z_thresh    <- config$ldsc$h2_z_threshold
  z_step      <- config$ldsc$h2_z_step %||% 0.5
  z_max       <- config$ldsc$h2_z_max  %||% 10.0

  smoothed   <- FALSE
  keep_idx   <- NULL
  S_filtered <- NULL
  S_raw      <- NULL

  repeat {
    retained   <- h2_table[z >= z_thresh]$trait
    keep_idx   <- which(trait_names %in% retained)
    S_filtered <- S[keep_idx, keep_idx, drop = FALSE]
    S_raw      <- S_filtered

    if (length(retained) < 2) {
      log_fatal("ldsc", sprintf(
        "Fewer than 2 traits pass h2_z_threshold=%.2f; cannot proceed", z_thresh))
    }

    if (check_psd(S_filtered)) {
      log_info("ldsc", sprintf("S matrix (%dx%d): all eigenvalues positive (PSD confirmed) at z_thresh=%.2f",
                               nrow(S_filtered), ncol(S_filtered), z_thresh))
      break
    }

    log_warn("ldsc", "S matrix is NOT positive semi-definite; applying Matrix::nearPD smoothing")
    smoothed_result <- Matrix::nearPD(S_filtered, corr = FALSE, keepDiag = TRUE,
                                       ensureSymmetry = TRUE)
    S_filtered  <- as.matrix(smoothed_result$mat)
    frob_delta  <- sqrt(sum((S_filtered - S_raw)^2))
    log_info("ldsc", sprintf("nearPD smoothing: Frobenius-norm delta = %.4g; %d iterations",
                             frob_delta, smoothed_result$iterations))
    log_warn("ldsc",
             "PSD smoothing applied to S only; V (sampling-covariance matrix of vech(S)) is unchanged. Downstream sandwich SEs in userGWAS may be slightly biased -- see GenomicSEM wiki on Matrix::nearPD use.")

    if (frob_delta <= halt_thresh) {
      smoothed <- TRUE
      break
    }

    new_thresh <- z_thresh + z_step
    if (new_thresh > z_max) {
      log_fatal("ldsc", sprintf(
        "Frobenius-norm smoothing delta %.4g > halt threshold %.4g even after raising h2_z_threshold to %.2f (z_max=%.2f). Consider lowering smoothing_frob_halt or checking input data.",
        frob_delta, halt_thresh, z_thresh, z_max))
    }
    log_warn("ldsc", sprintf(
      "Frobenius-norm delta %.4g > halt threshold %.4g; raising h2_z_threshold %.2f -> %.2f and retrying",
      frob_delta, halt_thresh, z_thresh, new_thresh))
    z_thresh <- new_thresh
  }

  h2_table[, pass := z >= z_thresh]
  log_info("ldsc", sprintf("Traits retained: %d/%d (%s) [final h2_z_threshold=%.2f]",
                           length(retained), length(trait_names),
                           paste(retained, collapse = ", "), z_thresh))

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

  categories <- load_trait_categories(config, sex)
  survival <- check_chapter_survival(
    retained, categories,
    min_indicators = config$cfa$min_indicators_per_factor %||% 3L,
    min_factors_required = config$ldsc$min_factors_required %||% 3L)
  fwrite(data.table(chapter = names(survival$counts),
                    n_retained = as.integer(survival$counts)),
         file.path(out_dir, "chapter_survival.csv"))

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
