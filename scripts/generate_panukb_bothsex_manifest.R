#!/usr/bin/env Rscript
# Derives manifest/panukb_bothsex_manifest.rda from the Pan-UKB phenotype
# manifest CSV, keeping only icd10/both_sexes rows with non-NA EUR counts.

find_project_root <- function() {
  full_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- full_args[grep("^--file=", full_args)]
  start_d <- if (length(file_arg) > 0L) {
    dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = FALSE))
  } else normalizePath(getwd(), mustWork = FALSE)
  d <- start_d
  while (d != dirname(d)) {
    if (file.exists(file.path(d, "config", "pipeline.yaml"))) return(d)
    d <- dirname(d)
  }
  stop("project root not found")
}

project_root <- find_project_root()
setwd(project_root)
source("R/00_setup.R")
config <- read_config("config/pipeline.yaml", root = project_root)

manifest_csv <- config$panukb$phenotype_manifest
if (is.null(manifest_csv) || !nzchar(manifest_csv) || !file.exists(manifest_csv)) {
  stop("Pan-UKB phenotype manifest not found: ", as.character(manifest_csv))
}

m <- data.table::fread(manifest_csv)
m <- m[trait_type == "icd10" & pheno_sex == "both_sexes"]

n_cases <- suppressWarnings(as.numeric(m$n_cases_EUR))
n_controls <- suppressWarnings(as.numeric(m$n_controls_EUR))
ok <- !is.na(n_cases) & !is.na(n_controls) & n_cases > 0

panukb_bothsex_manifest <- data.frame(
  phenotype  = as.character(m$phenocode[ok]),
  n_cases    = n_cases[ok],
  n_controls = n_controls[ok],
  stringsAsFactors = FALSE
)

out_path <- file.path(config$paths$manifest_dir, "panukb_bothsex_manifest.rda")
dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
save(panukb_bothsex_manifest, file = out_path)
message(sprintf("wrote %s (%d rows)", out_path, nrow(panukb_bothsex_manifest)))
