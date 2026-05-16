run_split_vcf <- function(config, sex) {
  log_info("split_vcf", sprintf("=== Split VCF stage: %s ===", sex))

  in_vcf <- file.path(config$paths$output_dir, sex, "gwas",
                      sprintf("%s_chapter_gsem.vcf.gz", sex))
  if (!file.exists(in_vcf)) {
    log_fatal("split_vcf", sprintf("input VCF not found: %s", in_vcf))
  }

  out_dir <- file.path(config$paths$output_dir, sex, "gwas_split")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  bcftools <- Sys.which("bcftools")
  if (!nzchar(bcftools)) {
    log_fatal("split_vcf",
              "bcftools not found in $PATH; required for VCF splitting. Install via apt/conda.")
  }

  samples_raw <- system2(bcftools, c("query", "-l", in_vcf), stdout = TRUE)
  samples <- samples_raw[nzchar(samples_raw)]
  if (length(samples) == 0L) {
    log_fatal("split_vcf", "no SAMPLE columns found in input VCF")
  }
  log_info("split_vcf", sprintf("Splitting %d sample column(s): %s",
                                 length(samples), paste(samples, collapse = ", ")))

  written <- character(0)
  for (s in samples) {
    out_vcf <- file.path(out_dir, sprintf("%s_%s_gsem.vcf.gz", sex, s))
    args <- c("view", "-Oz", "-s", s, "-o", out_vcf, in_vcf)
    status <- system2(bcftools, args, stdout = NULL, stderr = NULL)
    if (status != 0L || !file.exists(out_vcf)) {
      log_warn("split_vcf", sprintf("bcftools view failed for sample %s (status %d)", s, status))
      next
    }
    idx_status <- system2(bcftools, c("index", "-t", out_vcf),
                          stdout = NULL, stderr = NULL)
    if (idx_status != 0L) {
      log_warn("split_vcf", sprintf("bcftools index failed for %s (status %d)", out_vcf, idx_status))
    }
    written <- c(written, out_vcf)
    log_info("split_vcf", sprintf("  wrote %s", out_vcf))
  }

  write_stage_manifest("split_vcf", sex, config, written)
  invisible(written)
}
