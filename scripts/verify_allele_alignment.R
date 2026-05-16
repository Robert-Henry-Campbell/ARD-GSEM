#!/usr/bin/env Rscript
# Iterates meta/anchor_snps_by_chapter.csv; for the first chapter whose
# anchor_trait_code is in the retained-traits set AND the factor's VCF exists,
# checks ES sign at the anchor SNP matches expected_sign. Halts (status 1) on
# sign mismatch. Robust to z>=4 trait-dropout.

suppressPackageStartupMessages({
  library(data.table)
})

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

args <- commandArgs(trailingOnly = TRUE)
sex <- if (length(args) >= 1L) args[1] else "male"

project_root <- find_project_root()
setwd(project_root)
source("R/00_setup.R")
config <- read_config("config/pipeline.yaml", root = project_root)

anchors <- fread(file.path(config$paths$meta_dir, "anchor_snps_by_chapter.csv"))
h2 <- fread(file.path(config$paths$output_dir, sex, "ldsc", "h2_qc.csv"))
retained <- h2[pass == TRUE]$trait

vcf <- file.path(config$paths$output_dir, sex, "gwas",
                 sprintf("%s_chapter_gsem.vcf.gz", sex))
if (!file.exists(vcf)) stop(sprintf("VCF not found: %s", vcf))

bcftools <- Sys.which("bcftools")
if (!nzchar(bcftools)) stop("bcftools not in $PATH")

samples <- system2(bcftools, c("query", "-l", vcf), stdout = TRUE)

picked <- NULL
for (i in seq_len(nrow(anchors))) {
  a <- anchors[i]
  if (!(a$anchor_trait_code %in% retained)) next
  if (!(a$chapter %in% samples)) next
  picked <- a
  break
}
if (is.null(picked)) {
  stop("No anchor SNP applicable: no chapter has both an entry in anchor_snps_by_chapter.csv AND a surviving anchor_trait_code AND a SAMPLE column in the VCF")
}

cat(sprintf("Using anchor: chapter=%s, rsid=%s, expected_a1=%s, expected_sign=%s\n",
            picked$chapter, picked$rsid, picked$expected_a1, picked$expected_sign))

q <- system2(bcftools,
             c("view", "-s", picked$chapter, "-i", sprintf('ID=="%s"', picked$rsid),
               "-Ov", vcf),
             stdout = TRUE)
data_lines <- q[!startsWith(q, "#") & nzchar(q)]
if (length(data_lines) == 0L) stop(sprintf("SNP %s not in VCF", picked$rsid))

fields <- strsplit(data_lines[1], "\t", fixed = TRUE)[[1]]
ref <- fields[4]; alt <- fields[5]; fmt <- fields[9]; sample_val <- fields[10]
fmt_keys <- strsplit(fmt, ":", fixed = TRUE)[[1]]
sample_vals <- strsplit(sample_val, ":", fixed = TRUE)[[1]]
names(sample_vals) <- fmt_keys
es <- suppressWarnings(as.numeric(sample_vals["ES"]))
if (is.na(es)) stop(sprintf("ES missing for %s/%s", picked$rsid, picked$chapter))

cat(sprintf("Variant: REF=%s ALT=%s ES=%s\n", ref, alt, format(es)))

# Convention: VCF stores effect on ALT. If expected_a1 == ALT, expected_sign applies directly;
# if expected_a1 == REF, flip the sign expectation.
expected_sign_for_alt <- picked$expected_sign
if (picked$expected_a1 == ref) {
  expected_sign_for_alt <- ifelse(picked$expected_sign == "+", "-", "+")
}
observed_sign <- ifelse(es >= 0, "+", "-")
cat(sprintf("Expected sign on ALT=%s: %s; observed: %s\n",
            alt, expected_sign_for_alt, observed_sign))

if (observed_sign != expected_sign_for_alt) {
  cat("FAIL: allele-alignment direction mismatch\n", file = stderr())
  quit(status = 1L)
}
cat("PASS: allele-alignment direction matches expectation\n")
