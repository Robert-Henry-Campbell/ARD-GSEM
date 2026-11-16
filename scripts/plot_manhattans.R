#!/usr/bin/env Rscript
# Renders one Manhattan plot per latent factor of each finished factor GWAS into
# output/<sex>/plots/<sex>_<factor>_manhattan.png, without running the rest of the
# plots stage (QQ, rg heatmap, loadings bar, path diagram).
#
# Usage: Rscript scripts/plot_manhattans.R [--sex bothsex[,male,...]] [--config path]
# Default --sex: every sex with an output/<sex>/gwas/<sex>_userGWAS_raw.rds.

suppressPackageStartupMessages({
  library(data.table)
})

parse_args <- function(args) {
  opts <- list(sex = NULL, config = "config/pipeline.yaml")
  i <- 1
  while (i <= length(args)) {
    if (args[i] == "--sex") { opts$sex <- args[i + 1]; i <- i + 2 }
    else if (args[i] == "--config") { opts$config <- args[i + 1]; i <- i + 2 }
    else if (args[i] == "--help" || args[i] == "-h") {
      cat("Usage: Rscript scripts/plot_manhattans.R [--sex <list>] [--config <path>]\n\n")
      cat("  --sex     Comma-separated subset of {male, female, bothsex, bothsex_meta}.\n")
      cat("            Default: every sex with a finished factor GWAS.\n")
      cat("  --config  path to YAML config (default: config/pipeline.yaml)\n")
      quit(status = 0)
    }
    else { cat("Unknown argument:", args[i], "\n"); quit(status = 1) }
  }
  opts
}

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

opts <- parse_args(commandArgs(trailingOnly = TRUE))

project_root <- find_project_root()
setwd(project_root)
source("R/00_setup.R")
source("R/09_write_vcf.R")
source("R/10_plots.R")
config <- read_config(opts$config, root = project_root)

finished <- function(sex) {
  file.exists(file.path(config$paths$output_dir, sex, "gwas",
                        paste0(sex, "_userGWAS_raw.rds")))
}
sexes <- if (is.null(opts$sex)) {
  Filter(finished, c("male", "female", "bothsex", "bothsex_meta"))
} else {
  trimws(strsplit(opts$sex, ",", fixed = TRUE)[[1L]])
}
if (length(sexes) == 0L) {
  stop(sprintf("No finished factor GWAS under %s", config$paths$output_dir), call. = FALSE)
}

written <- unlist(lapply(sexes, function(sex) run_manhattans(config, sex)))
cat(sprintf("\nWrote %d Manhattan plot(s):\n%s\n", length(written),
            paste0("  ", written, collapse = "\n")))
