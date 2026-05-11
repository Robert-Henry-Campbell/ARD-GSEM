factor_summary_row <- function(g, factor_name, sig_threshold = 5e-8) {
  n_snps <- nrow(g)
  z <- ifelse(!is.na(g$ES) & !is.na(g$SE) & g$SE > 0, g$ES / g$SE, NA_real_)
  mean_chisq <- mean(z^2, na.rm = TRUE)
  lambda <- compute_lambda_gc(g$P)
  n_gws <- sum(!is.na(g$P) & g$P < sig_threshold)
  n_qsnp <- sum(!is.na(g$Q_SNP_pval) & g$Q_SNP_pval < sig_threshold)
  data.table(
    factor = factor_name,
    n_snps = n_snps,
    mean_chisq = mean_chisq,
    lambda_gc = lambda,
    n_genome_wide_sig = n_gws,
    n_qsnp_het_sig = n_qsnp
  )
}

top_snps_table <- function(g, factor_name, n_top = 20L) {
  d <- as.data.table(g)
  d <- d[!is.na(P)]
  if (nrow(d) == 0L) return(data.table())
  setorder(d, P)
  d <- head(d, n_top)
  d[, factor := factor_name]
  cols <- intersect(c("factor", "SNP", "CHR", "BP", "A1", "A2", "AF",
                      "ES", "SE", "P", "Q_SNP", "Q_SNP_pval"),
                    names(d))
  d[, ..cols]
}

run_factor_summary <- function(config, sex, sig_threshold = 5e-8, n_top = 20L) {
  log_info("factor_summary", sprintf("=== Factor summary stage: %s ===", sex))

  gwas_path <- file.path(config$paths$output_dir, sex, "gwas",
                         paste0(sex, "_userGWAS_raw.rds"))
  if (!file.exists(gwas_path)) {
    log_warn("factor_summary",
             sprintf("userGWAS_raw not found at %s; skipping", gwas_path))
    return(invisible(NULL))
  }
  gwas_result <- readRDS(gwas_path)

  snp_path <- file.path(config$paths$output_dir, sex, "sumstats",
                        paste0(sex, "_snp_sumstats.rds"))
  snp_sumstats <- if (file.exists(snp_path)) readRDS(snp_path) else NULL

  out_dir <- file.path(config$paths$output_dir, sex, "gwas")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  summary_rows <- list()
  top_rows <- list()
  for (f in names(gwas_result)) {
    g <- normalize_gwas_columns(gwas_result[[f]], snp_sumstats)
    summary_rows[[f]] <- factor_summary_row(g, f, sig_threshold = sig_threshold)
    top_rows[[f]] <- top_snps_table(g, f, n_top = n_top)
    log_info("factor_summary",
             sprintf("  %s: n_snps=%s, lambda_GC=%.3f, n_sig=%d, n_qsnp_het=%d",
                     f, format(summary_rows[[f]]$n_snps, big.mark = ","),
                     summary_rows[[f]]$lambda_gc,
                     summary_rows[[f]]$n_genome_wide_sig,
                     summary_rows[[f]]$n_qsnp_het_sig))
  }

  summary_dt <- rbindlist(summary_rows)
  summary_dt[, sex := sex]
  setcolorder(summary_dt, c("sex", setdiff(names(summary_dt), "sex")))
  summary_path <- file.path(out_dir, paste0(sex, "_factor_summary.csv"))
  fwrite(summary_dt, summary_path)

  top_dt <- rbindlist(top_rows, fill = TRUE)
  if (nrow(top_dt) > 0L) top_dt[, sex := sex]
  top_path <- file.path(out_dir, paste0(sex, "_top_snps.csv"))
  fwrite(top_dt, top_path)

  log_info("factor_summary", sprintf("Wrote %s and %s",
                                      basename(summary_path), basename(top_path)))

  write_stage_manifest("factor_summary", sex, config, c(summary_path, top_path))

  list(summary = summary_dt, top_snps = top_dt)
}
