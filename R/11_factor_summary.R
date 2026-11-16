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

# Post-GWAS (bifactor) summary: one row per factor from the qc_filters, loci, factor_ldsc
# and mr_export outputs, plus top SNPs of the full and qfilter tiers with their flags.
run_factor_summary_postgwas <- function(config, sex, n_top = 20L) {
  stage <- "factor_summary"
  out_dir <- stage_dir(config, sex, stage)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  read_opt <- function(st, file) {
    p <- file.path(stage_dir(config, sex, st), file)
    if (file.exists(p)) fread(p) else NULL
  }
  s <- read_opt("qc_filters", "qc_summary.csv")
  if (is.null(s)) log_fatal(stage, "qc_filters/qc_summary.csv missing (run qc_filters)")
  calib <- read_opt("qc_filters", "qsnp_calibration.csv")
  if (!is.null(calib)) s <- merge(s, calib[, .(factor, df_eff, lambda_Q_gsem, lambda_Q_eff)], by = "factor", all.x = TRUE)
  loci <- read_opt("loci", "loci_summary.csv")
  if (!is.null(loci)) s <- merge(s, loci[, !"n_gws"], by = "factor", all.x = TRUE)
  ld <- read_opt("factor_ldsc", "summary.csv")
  if (!is.null(ld)) {
    w <- dcast(ld[tier %in% c("qfilter", "qsnp")], factor ~ tier,
               value.var = c("intercept", "intercept_se", "h2", "h2_se", "mean_chisq"))
    s <- merge(s, w, by = "factor", all.x = TRUE)
  }
  mr <- read_opt("mr_export", "mr_export_summary.csv")
  if (!is.null(mr)) s <- merge(s, mr[, .(factor, n_instruments)], by = "factor", all.x = TRUE)
  excl <- config$postgwas$excluded_factors %||% list()
  if (length(excl)) {
    s <- rbind(s, data.table(factor = names(excl), note = unlist(excl)), fill = TRUE)
  }
  s[, sex := sex]
  fwrite(s, file.path(out_dir, paste0(sex, "_factor_summary.csv")))

  top <- rbindlist(lapply(postgwas_factors(config), function(f) {
    d <- load_postgwas(config, sex, f, "full",
                       cols = c("SNP", "CHR", "BP", "A1", "A2", "EAF_A1", "beta", "SE", "P", "LOG10P",
                                "Q_SNP_pval", "qsnp_hit", "qsnp_ld", "pass_qfilter", "is_mhc"))
    rbind(d[order(-LOG10P)][seq_len(min(.N, n_top))][, tier := "full"],
          d[pass_qfilter == TRUE][order(-LOG10P)][seq_len(min(.N, n_top))][, tier := "qfilter"])[, factor := f]
  }))
  setcolorder(top, c("factor", "tier"))
  fwrite(top, file.path(out_dir, paste0(sex, "_top_snps.csv")))
  clear_postgwas_cache()
  log_info(stage, sprintf("Wrote %s_factor_summary.csv and %s_top_snps.csv", sex, sex))
  write_stage_manifest(stage, sex, config, list.files(out_dir, full.names = TRUE))
  list(summary = s, top_snps = top)
}

run_factor_summary <- function(config, sex, sig_threshold = 5e-8, n_top = 20L) {
  log_info("factor_summary", sprintf("=== Factor summary stage: %s ===", sex))
  if (exists("factor_registry_path", mode = "function") &&
      file.exists(factor_registry_path(config, sex))) {
    return(run_factor_summary_postgwas(config, sex, n_top))
  }

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
