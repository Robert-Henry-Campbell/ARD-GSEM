#!/usr/bin/env Rscript
# Smoke test: full pipeline on shared traits (D50, D64, E10, E11) both sexes
# Requires: reference data downloaded, renv restored, GenomicSEM installed
# Expected runtime: <15 min on 24 cores
# Exit 0 = pass, exit 1 = fail

setwd("/mnt/sdg/robert/ardmr/GSEM")
source("R/00_setup.R")
source("R/01_munge.R")
source("R/02_ldsc.R")
source("R/03_efa.R")
source("R/04_cfa.R")
source("R/05_comparison.R")
source("R/06_report.R")

config <- read_config("config/pipeline.yaml")
config <- setup_pipeline(config, sex = "both", mode = "smoke", threads = 4)

errors <- character(0)
ldsc_results_cache <- list()

for (sex in c("male", "female")) {
  log_info("smoke", paste("=== Starting", sex, "==="))

  # Stage 1: Munge
  tryCatch({
    munge_results <- run_munge(config, sex = sex)
    n_munged <- munge_results$n_munged
    if (n_munged < 2) stop(paste("Only", n_munged, "traits munged; need >=2"))
    log_info("smoke", paste(n_munged, "traits munged successfully"))
  }, error = function(e) {
    errors <<- c(errors, paste(sex, "munge:", e$message))
    log_error("smoke", paste(sex, "munge FAILED:", e$message))
  })

  # Stage 2: LDSC
  tryCatch({
    ldsc_results <- run_ldsc(config, sex = sex)
    ldsc_results_cache[[sex]] <- ldsc_results
    if (!check_psd(ldsc_results$S)) log_warn("smoke", "S matrix not PSD")
    if (any(is.na(ldsc_results$S))) stop("NAs in S matrix")
    log_info("smoke", paste("LDSC complete;", ldsc_results$n_valid, "traits pass h2 filter"))
  }, error = function(e) {
    errors <<- c(errors, paste(sex, "ldsc:", e$message))
    log_error("smoke", paste(sex, "ldsc FAILED:", e$message))
  })

  # Stage 3: EFA (only if >=3 traits pass)
  tryCatch({
    if (!is.null(ldsc_results_cache[[sex]]) && ldsc_results_cache[[sex]]$n_valid >= 3) {
      efa_results <- run_efa(config, sex = sex, ldsc_results = ldsc_results_cache[[sex]])
      log_info("smoke", paste("EFA:", efa_results$n_factors, "factor(s) extracted"))
    } else {
      log_warn("smoke", paste("Skipping EFA for", sex, ": <3 valid traits"))
    }
  }, error = function(e) {
    errors <<- c(errors, paste(sex, "efa:", e$message))
    log_error("smoke", paste(sex, "efa FAILED:", e$message))
  })

  # Stage 4: CFA
  tryCatch({
    cfa_results <- run_cfa(config, sex = sex, ldsc_results = ldsc_results_cache[[sex]])
    if (!cfa_results$converged) log_warn("smoke", paste(sex, "CFA did not converge"))
    else log_info("smoke", paste(sex, "CFA converged"))
  }, error = function(e) {
    errors <<- c(errors, paste(sex, "cfa:", e$message))
    log_error("smoke", paste(sex, "cfa FAILED:", e$message))
  })
}

# Stage 5: Comparison
tryCatch({
  comp_results <- run_comparison(config)
  log_info("smoke", "Sex comparison complete")
}, error = function(e) {
  errors <<- c(errors, paste("comparison:", e$message))
  log_error("smoke", paste("comparison FAILED:", e$message))
})

# Stage 6: Report
tryCatch({
  suppressPackageStartupMessages(library(rmarkdown))
  run_report(config)
  log_info("smoke", "Report rendered")
}, error = function(e) {
  log_warn("smoke", paste("Report rendering failed (non-fatal):", e$message))
})

# Final verdict
cat("\n")
if (length(errors) == 0) {
  log_info("smoke", "=== ALL STAGES PASSED ===")
  close_logging()
  cat("SMOKE TEST PASSED\n")
  quit(status = 0)
} else {
  log_error("smoke", paste("FAILURES:", length(errors)))
  for (err in errors) log_error("smoke", err)
  close_logging()
  cat("SMOKE TEST FAILED\n")
  cat(paste(" -", errors, collapse = "\n"), "\n")
  quit(status = 1)
}
