library(data.table)

find_project_root <- function() {
  start_d <- normalizePath(".", winslash = "/", mustWork = FALSE)
  d <- start_d
  visited <- character(0)
  while (d != dirname(d)) {
    visited <- c(visited, d)
    if (file.exists(file.path(d, "config", "pipeline.yaml"))) return(d)
    d <- dirname(d)
  }
  stop(sprintf(
    "project root not found.\n  CWD: %s\n  walked: %s\n  hint: run tests from the repo or any subdirectory containing 'config/pipeline.yaml'.",
    getwd(), paste(visited, collapse = " -> ")))
}

project_root <- find_project_root()
source(file.path(project_root, "R/utils.R"))
source(file.path(project_root, "R/05_comparison.R"))

fixtures_dir <- file.path(project_root, "tests/fixtures")

skip_if_no_reference <- function() {
  ref_dir <- file.path(project_root, "reference")
  if (!dir.exists(file.path(ref_dir, "eur_w_ld_chr"))) {
    skip("Reference data not downloaded")
  }
}

skip_if_no_sumstats <- function() {
  if (!dir.exists(file.path(project_root, "sumstats/male"))) {
    skip("Sumstats not available")
  }
}
