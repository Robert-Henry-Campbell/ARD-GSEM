run_comparison <- function(config, male_cfa = NULL, female_cfa = NULL) {
  log_info("comparison", "=== Sex comparison stage ===")

  out_dir <- file.path(config$paths$output_dir, "comparison")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  male_h2 <- fread(file.path(config$paths$output_dir, "male", "ldsc", "h2_qc.csv"))
  female_h2 <- fread(file.path(config$paths$output_dir, "female", "ldsc", "h2_qc.csv"))

  male_retained <- male_h2[pass == TRUE]$trait
  female_retained <- female_h2[pass == TRUE]$trait
  shared_valid <- intersect(male_retained, female_retained)

  log_info("comparison", sprintf("Shared valid traits: %s (intersection of h2-passing sets)",
                                 paste(shared_valid, collapse = ", ")))

  if (length(shared_valid) < 2) {
    log_warn("comparison", "Fewer than 2 shared valid traits; comparison limited")
  }

  # Load CFA results
  male_model <- load_cfa_result(config, "male")
  female_model <- load_cfa_result(config, "female")

  if (is.null(male_model) || is.null(female_model)) {
    log_warn("comparison", "Cannot compare: missing CFA results for one or both sexes")
    return(list(shared_valid = shared_valid, comparison_table = NULL))
  }

  male_loadings <- extract_loadings_table(male_model)
  female_loadings <- extract_loadings_table(female_model)

  if (is.null(male_loadings) || is.null(female_loadings)) {
    log_warn("comparison", "Cannot extract loadings from CFA results")
    return(list(shared_valid = shared_valid, comparison_table = NULL))
  }

  alignment_warnings <- check_factor_alignment(male_loadings, female_loadings)
  if (length(alignment_warnings) > 0) {
    for (w in alignment_warnings) {
      log_warn("comparison", paste("Factor alignment:", w))
    }
    log_warn("comparison", "Factor labels may represent different constructs across sexes; interpret Z-tests with caution")
  }

  shared_traits_in_model <- intersect(male_loadings$trait, female_loadings$trait)
  shared_traits_in_model <- intersect(shared_traits_in_model, shared_valid)

  if (length(shared_traits_in_model) == 0) {
    log_warn("comparison", "No shared traits in both CFA models")
    return(list(shared_valid = shared_valid, comparison_table = NULL))
  }

  log_info("comparison", "Loading differences:")
  comparison <- data.table(
    trait = character(), factor = character(),
    male_loading = numeric(), male_se = numeric(),
    female_loading = numeric(), female_se = numeric(),
    diff = numeric(), z_diff = numeric(), p_diff = numeric()
  )

  for (trait in shared_traits_in_model) {
    m_row <- male_loadings[male_loadings$trait == trait, ]
    f_row <- female_loadings[female_loadings$trait == trait, ]

    if (nrow(m_row) > 0 && nrow(f_row) > 0) {
      result <- compute_loading_diff(
        male_loading = m_row$loading[1],
        male_se = m_row$se[1],
        female_loading = f_row$loading[1],
        female_se = f_row$se[1]
      )
      comparison <- rbind(comparison, data.table(
        trait = trait,
        factor = m_row$factor[1],
        male_loading = m_row$loading[1],
        male_se = m_row$se[1],
        female_loading = f_row$loading[1],
        female_se = f_row$se[1],
        diff = result$diff,
        z_diff = result$z_diff,
        p_diff = result$p_diff
      ))
      log_info("comparison", sprintf("  %s: male=%.3f, female=%.3f, diff=%.3f, Z=%.2f, p=%.3f",
                                     trait, m_row$loading[1], f_row$loading[1],
                                     result$diff, result$z_diff, result$p_diff))
    }
  }

  if (nrow(comparison) > 0) {
    comparison[, p_bonferroni := pmin(p_diff * .N, 1)]
    comparison[, p_fdr := p.adjust(p_diff, method = "BH")]

    sig_bonf <- comparison[p_bonferroni < 0.05]
    if (nrow(sig_bonf) > 0) {
      log_info("comparison", sprintf("%d significant sex differences at p<0.05 (Bonferroni)", nrow(sig_bonf)))
    } else {
      log_info("comparison", "No significant sex differences at p<0.05 (Bonferroni-adjusted)")
    }
  }

  fwrite(comparison, file.path(out_dir, "loading_diff.csv"))

  write_stage_manifest("comparison", "both", config, list.files(out_dir, full.names = TRUE))

  list(shared_valid = shared_valid, comparison_table = comparison)
}

load_cfa_result <- function(config, sex) {
  efa_path <- file.path(config$paths$output_dir, sex, "cfa", "model_efa.rds")
  apriori_path <- file.path(config$paths$output_dir, sex, "cfa", "model_apriori.rds")
  if (file.exists(efa_path)) return(readRDS(efa_path))
  if (file.exists(apriori_path)) return(readRDS(apriori_path))
  NULL
}

extract_loadings_table <- function(model_result) {
  if (is.null(model_result$results)) return(NULL)
  res <- as.data.table(model_result$results)
  if (!"op" %in% names(res)) return(NULL)
  loadings <- res[op == "=~"]
  if (nrow(loadings) == 0) return(NULL)
  data.table(
    trait = loadings$rhs,
    factor = loadings$lhs,
    loading = loadings$Unstand_Est,
    se = loadings$Unstand_SE
  )
}

check_factor_alignment <- function(male_loadings, female_loadings) {
  male_factors <- unique(male_loadings$factor)
  female_factors <- unique(female_loadings$factor)
  warnings <- character(0)

  if (length(male_factors) != length(female_factors)) {
    warnings <- c(warnings, sprintf(
      "Different number of factors: male=%d, female=%d",
      length(male_factors), length(female_factors)))
  }

  shared_factors <- intersect(male_factors, female_factors)
  for (f in shared_factors) {
    male_traits <- sort(male_loadings[factor == f]$trait)
    female_traits <- sort(female_loadings[factor == f]$trait)
    if (!identical(male_traits, female_traits)) {
      warnings <- c(warnings, sprintf(
        "Factor %s has different indicators: male={%s}, female={%s}",
        f, paste(male_traits, collapse = ","), paste(female_traits, collapse = ",")))
    }
  }

  warnings
}

