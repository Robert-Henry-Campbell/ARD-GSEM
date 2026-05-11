#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

parse_args <- function(args) {
  opts <- list(sex = "both", mode = "smoke", stage = "all", threads = 24, resume = FALSE)
  i <- 1
  while (i <= length(args)) {
    if (args[i] == "--sex") { opts$sex <- args[i + 1]; i <- i + 2 }
    else if (args[i] == "--mode") { opts$mode <- args[i + 1]; i <- i + 2 }
    else if (args[i] == "--stage") { opts$stage <- args[i + 1]; i <- i + 2 }
    else if (args[i] == "--threads") { opts$threads <- as.integer(args[i + 1]); i <- i + 2 }
    else if (args[i] == "--resume") { opts$resume <- TRUE; i <- i + 1 }
    else if (args[i] == "--help" || args[i] == "-h") {
      cat("Usage: Rscript scripts/run_pipeline.R [options]\n\n")
      cat("Options:\n")
      cat("  --sex       male|female|both  (default: both)\n")
      cat("  --mode      smoke|full        (default: smoke)\n")
      cat("  --stage     all|munge|ldsc|efa|cfa|comparison|report  (default: all)\n")
      cat("  --threads   N                 (default: 24)\n")
      cat("  --resume    Skip stages with valid output manifests\n")
      cat("  --help      Show this message\n")
      quit(status = 0)
    }
    else { cat("Unknown argument:", args[i], "\n"); quit(status = 1) }
  }
  opts
}

opts <- parse_args(args)

setwd("/mnt/sdg/robert/ardmr/GSEM")
source("R/00_setup.R")
source("R/01_munge.R")
source("R/02_ldsc.R")
source("R/03_efa.R")
source("R/04_cfa.R")
source("R/05_comparison.R")
source("R/06_report.R")

config <- read_config("config/pipeline.yaml")
config <- setup_pipeline(config, sex = opts$sex, mode = opts$mode, threads = opts$threads)

sexes <- if (opts$sex == "both") c("male", "female") else opts$sex
stages <- if (opts$stage == "all") c("munge", "ldsc", "efa", "cfa", "comparison", "report") else opts$stage

should_run <- function(stage_name) {
  stage_name %in% stages
}

should_skip <- function(stage_name, sex_val) {
  if (!opts$resume) return(FALSE)
  manifest <- read_stage_manifest(stage_name, sex_val, config)
  config_hash <- digest::digest(config[[stage_name]] %||% config, algo = "sha256")
  fresh <- stage_is_fresh(manifest, config_hash)
  if (fresh) {
    log_info("setup", sprintf("--resume: %s/%s manifest fresh; skipping", stage_name, sex_val))
  }
  fresh
}

tryCatch({

  if (should_run("munge")) {
    for (sex in sexes) {
      if (!should_skip("munge", sex)) run_munge(config, sex)
    }
  }

  if (should_run("ldsc")) {
    for (sex in sexes) {
      if (!should_skip("ldsc", sex)) run_ldsc(config, sex)
    }
  }

  if (should_run("efa")) {
    for (sex in sexes) {
      if (!should_skip("efa", sex)) {
        h2 <- fread(file.path(config$paths$output_dir, sex, "ldsc", "h2_qc.csv"))
        n_valid <- sum(h2$pass)
        if (n_valid >= 3) {
          run_efa(config, sex)
        } else {
          log_warn("efa", sprintf("Skipping EFA for %s: only %d valid traits (need >=3)", sex, n_valid))
        }
      }
    }
  }

  if (should_run("cfa")) {
    for (sex in sexes) {
      if (!should_skip("cfa", sex)) run_cfa(config, sex)
    }
  }

  if (should_run("comparison") && opts$sex == "both") {
    if (!should_skip("comparison", "both")) run_comparison(config)
  }

  if (should_run("report")) {
    suppressPackageStartupMessages(library(rmarkdown))
    run_report(config)
  }

}, error = function(e) {
  log_error("pipeline", sprintf("Pipeline failed: %s", e$message))
  close_logging()
  quit(status = 1)
})

close_logging()
cat("\nPipeline completed successfully.\n")
