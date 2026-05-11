suppressPackageStartupMessages({
  library(data.table)
  library(yaml)
  library(jsonlite)
})

project_root <- "/mnt/sdg/robert/ardmr/GSEM"
setwd(project_root)

source("R/utils.R")

setup_pipeline <- function(config, sex = "both", mode = "smoke", threads = NULL) {
  if (!is.null(threads)) config$parallel$n_workers <- threads

  init_logging(config)
  log_info("setup", sprintf("Config: %s", config$project$name))
  log_info("setup", sprintf("Threads: %d | Mode: %s | Sex: %s",
                            config$parallel$n_workers, mode, sex))

  validate_reference(config)

  sexes <- if (sex == "both") c("male", "female") else sex
  for (s in sexes) {
    traits <- discover_traits(config, s)
    log_info("setup", sprintf("Traits detected (%s): %s", s, paste(traits, collapse = ", ")))
  }

  if (sex == "both") {
    shared <- get_shared_traits(config)
    log_info("setup", sprintf("Shared traits: %s", paste(shared, collapse = ", ")))
  }

  for (s in sexes) {
    for (stage in c("munge", "ldsc", "efa", "cfa")) {
      dir.create(file.path(config$paths$output_dir, s, stage),
                 recursive = TRUE, showWarnings = FALSE)
    }
  }
  dir.create(file.path(config$paths$output_dir, "comparison"),
             recursive = TRUE, showWarnings = FALSE)

  config
}

validate_reference <- function(config) {
  ld_dir <- config$paths$ld_scores
  hm3 <- config$paths$hm3_snplist

  if (!dir.exists(ld_dir)) {
    log_fatal("setup", sprintf("LD scores directory not found: %s", ld_dir))
  }
  ld_files <- list.files(ld_dir, pattern = "\\.l2\\.ldscore\\.gz$")
  if (length(ld_files) < 22) {
    log_fatal("setup", sprintf("Expected 22 LD score files, found %d in %s", length(ld_files), ld_dir))
  }

  if (!file.exists(hm3)) {
    log_fatal("setup", sprintf("HM3 SNP list not found: %s", hm3))
  }

  log_info("setup", sprintf("Reference validated: eur_w_ld_chr (%d chr), w_hm3.snplist", length(ld_files)))
}
