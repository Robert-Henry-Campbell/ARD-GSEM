#!/usr/bin/env Rscript

setwd("/mnt/sdg/robert/ardmr/GSEM")

if (!requireNamespace("renv", quietly = TRUE)) {
  install.packages("renv")
}

renv::init(bare = TRUE, restart = FALSE)

cat("Installing CRAN packages...\n")
renv::install(c(
  "data.table",
  "lavaan",
  "psych",
  "GPArotation",
  "furrr",
  "future",
  "yaml",
  "Matrix",
  "ggplot2",
  "dplyr",
  "R.utils",
  "rmarkdown",
  "knitr",
  "pheatmap",
  "remotes",
  "jsonlite",
  "digest",
  "testthat"
))

cat("Installing GenomicSEM from GitHub...\n")
renv::install("GenomicSEM/GenomicSEM")

cat("Snapshotting renv...\n")
renv::snapshot(prompt = FALSE)

cat("\n=== Installation complete ===\n")
cat("Verify: Rscript -e 'library(GenomicSEM); cat(\"GenomicSEM OK\\n\")'\n")
