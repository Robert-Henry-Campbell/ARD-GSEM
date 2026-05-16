run_factor_ldsc <- function(config, sex) {
  log_info("factor_ldsc", sprintf("=== factor LDSC stage: %s ===", sex))

  vcf_path <- file.path(config$paths$output_dir, sex, "gwas",
                        sprintf("%s_chapter_gsem.vcf.gz", sex))
  if (!file.exists(vcf_path)) {
    log_fatal("factor_ldsc", sprintf("VCF not found: %s", vcf_path))
  }

  out_dir <- file.path(config$paths$output_dir, sex, "factor_ldsc")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  gwas_raw_path <- file.path(config$paths$output_dir, sex, "gwas",
                             paste0(sex, "_userGWAS_raw.rds"))
  gwas_raw <- readRDS(gwas_raw_path)
  factor_names <- names(gwas_raw)

  ld_path <- config$paths$ld_scores
  wld_path <- config$paths$ld_scores

  summary_rows <- list()

  for (f in factor_names) {
    g <- gwas_raw[[f]]
    g <- as.data.table(g)
    snp_col <- intersect(c("SNP", "rsid"), names(g))[1]
    a1_col  <- intersect(c("A1", "EA"), names(g))[1]
    a2_col  <- intersect(c("A2", "OA"), names(g))[1]
    est_col <- intersect(c("est", "Estimate", "Unstand_Est", "Effect"), names(g))[1]
    se_col  <- intersect(c("SE", "Std_Error", "Unstand_SE"), names(g))[1]
    p_col   <- intersect(c("Pval_Estimate", "P", "pval", "p_value"), names(g))[1]
    maf_col <- intersect(c("MAF", "EAF"), names(g))[1]
    if (any(is.na(c(snp_col, a1_col, a2_col, est_col, se_col, p_col, maf_col)))) {
      log_warn("factor_ldsc", sprintf("factor %s missing required columns; skipping", f))
      next
    }

    nhat <- compute_factor_nhat(g[[maf_col]], g[[se_col]],
                                maf_threshold = config$gwas$nhat_maf_threshold %||% 0.10)
    if (is.na(nhat) || nhat <= 0) {
      log_warn("factor_ldsc", sprintf("factor %s: N_hat is NA/<=0; skipping", f))
      next
    }

    dt <- data.table(
      SNP = g[[snp_col]],
      A1  = g[[a1_col]],
      A2  = g[[a2_col]],
      Z   = g[[est_col]] / g[[se_col]],
      P   = g[[p_col]],
      N   = round(nhat)
    )
    dt <- dt[!is.na(SNP) & !is.na(Z) & !is.na(P) & is.finite(Z) & is.finite(P)]
    dt <- dt[A1 %in% c("A","C","G","T") & A2 %in% c("A","C","G","T")]
    sumstats_path <- file.path(out_dir, sprintf("%s.sumstats.gz", f))
    gz <- gzfile(sumstats_path, "wt")
    fwrite(dt, gz, sep = "\t")
    close(gz)
    log_info("factor_ldsc", sprintf("factor %s: %s SNPs written; N_hat=%.0f",
                                     f, format(nrow(dt), big.mark = ","), nhat))

    saved_wd <- getwd()
    res <- tryCatch({
      setwd(out_dir)
      GenomicSEM::ldsc(
        traits = normalizePath(sumstats_path, mustWork = TRUE),
        sample.prev = NA,
        population.prev = NA,
        ld = ld_path,
        wld = wld_path,
        trait.names = f,
        stand = FALSE
      )
    }, error = function(e) {
      log_warn("factor_ldsc", sprintf("ldsc() failed for factor %s: %s", f, conditionMessage(e)))
      NULL
    }, finally = setwd(saved_wd))

    if (is.null(res)) next

    h2 <- as.numeric(diag(res$S))[1]
    h2_se <- sqrt(as.numeric(diag(res$V)))[1]
    intercept <- if (!is.null(res$I)) as.numeric(diag(res$I))[1] else NA_real_
    summary_rows[[f]] <- data.table(
      factor = f, n_snps = nrow(dt), N_hat = nhat,
      h2 = h2, h2_se = h2_se, intercept = intercept)

    writeLines(c(
      sprintf("factor: %s", f),
      sprintf("n_snps: %d", nrow(dt)),
      sprintf("N_hat: %.0f", nhat),
      sprintf("h2: %.4f (SE %.4f)", h2, h2_se),
      sprintf("intercept: %.4f", intercept)
    ), file.path(out_dir, sprintf("%s_intercept.txt", f)))

    lo <- config$factor_ldsc$intercept_warn_low %||% 0.9
    hi <- config$factor_ldsc$intercept_warn_high %||% 1.1
    if (!is.na(intercept) && (intercept < lo || intercept > hi)) {
      log_warn("factor_ldsc", sprintf(
        "factor %s: intercept=%.4f outside [%.2f, %.2f] -- possible confounding",
        f, intercept, lo, hi))
    }
  }

  if (length(summary_rows) > 0L) {
    summary_dt <- rbindlist(summary_rows)
    fwrite(summary_dt, file.path(out_dir, "summary.csv"))
    log_info("factor_ldsc", sprintf("Wrote summary.csv with %d factor row(s)", nrow(summary_dt)))
  }

  write_stage_manifest("factor_ldsc", sex, config, list.files(out_dir, full.names = TRUE))
  invisible(summary_rows)
}
