library(data.table)

find_project_root <- function() {
  d <- normalizePath(".", winslash = "/", mustWork = FALSE)
  while (d != dirname(d)) {
    if (file.exists(file.path(d, "config", "pipeline.yaml"))) return(d)
    d <- dirname(d)
  }
  stop("project root not found (no config/pipeline.yaml ancestor of ", getwd(), ")")
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
