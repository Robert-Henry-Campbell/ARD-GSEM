library(data.table)

project_root <- "/mnt/sdg/robert/ardmr/GSEM"
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
