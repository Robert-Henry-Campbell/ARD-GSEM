#!/usr/bin/env Rscript
# Derives meta/ukb_pop_prev_{sex}.csv (columns: code, K_pop) from the Neale
# manifests. K_pop here is UKB cohort case fraction -- a documented proxy
# for population prevalence, not literature K_pop. See methods caveat.

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

for (sex in c("male", "female")) {
  manifest_file <- file.path(config$paths$manifest_dir,
                             sprintf("neale_%s_manifest.rda", sex))
  if (!file.exists(manifest_file)) {
    message(sprintf("manifest missing for %s; skipping", sex)); next
  }
  env <- new.env()
  load(manifest_file, envir = env)
  objs <- ls(env)
  stopifnot(length(objs) == 1L)
  m <- get(objs[1], envir = env)
  m <- as.data.frame(m)
  out <- data.frame(
    code = m$phenotype,
    K_pop = m$n_cases / (m$n_cases + m$n_controls),
    stringsAsFactors = FALSE
  )
  out_path <- file.path(config$paths$meta_dir, sprintf("ukb_pop_prev_%s.csv", sex))
  data.table::fwrite(out, out_path)
  message(sprintf("wrote %s (%d rows)", out_path, nrow(out)))
}
