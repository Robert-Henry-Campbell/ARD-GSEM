find_col <- function(dt, candidates) {
  hit <- intersect(candidates, names(dt))
  if (length(hit) == 0L) return(NA_character_)
  hit[1]
}

clump_qsnp <- function(qsnp_tsv, bfile, out_prefix, plink = "plink",
                        p1 = 5e-8, p2 = 1.0, r2 = 0.1, kb = 250) {
  cmd <- c(
    "--bfile", bfile,
    "--clump", qsnp_tsv,
    "--clump-p1", as.character(p1),
    "--clump-p2", as.character(p2),
    "--clump-r2", as.character(r2),
    "--clump-kb", as.character(kb),
    "--clump-field", "Q_SNP_pval",
    "--clump-snp-field", "SNP",
    "--out", out_prefix
  )
  status <- system2(plink, cmd, stdout = NULL, stderr = NULL)
  clumped_file <- paste0(out_prefix, ".clumped")
  if (!file.exists(clumped_file)) {
    if (status != 0L) {
      log_warn("qc_filters", sprintf("plink --clump returned status %d and no .clumped file; treating as no hits", status))
    } else {
      log_info("qc_filters", "plink --clump produced no .clumped file (no hits)")
    }
    return(character(0))
  }
  clumped <- tryCatch(fread(clumped_file), error = function(e) NULL)
  if (is.null(clumped) || nrow(clumped) == 0L) return(character(0))
  index_snps <- clumped$SNP
  buddies <- unique(unlist(strsplit(paste(clumped$SP2, collapse = ","), ",", fixed = TRUE)))
  buddies <- sub("\\(.*\\)$", "", buddies)
  buddies <- buddies[nzchar(buddies) & buddies != "NONE"]
  unique(c(index_snps, buddies))
}

run_qc_filters <- function(config, sex) {
  log_info("qc_filters", sprintf("=== QC filters stage: %s ===", sex))

  gwas_raw_path <- file.path(config$paths$output_dir, sex, "gwas",
                             paste0(sex, "_userGWAS_raw.rds"))
  if (!file.exists(gwas_raw_path)) {
    log_fatal("qc_filters", sprintf("userGWAS raw not found at %s", gwas_raw_path))
  }
  gwas_raw <- readRDS(gwas_raw_path)
  factor_names <- names(gwas_raw)

  out_dir <- file.path(config$paths$output_dir, sex, "qc_filters")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  qcfg <- config$qc_filters %||% list()
  qsnp_thresh <- as.numeric(qcfg$qsnp_pval_threshold %||% 5e-8)
  zsm_thresh <- as.numeric(qcfg$zsmooth_threshold %||% 1.96)
  r2 <- as.numeric(qcfg$clump_r2 %||% 0.1)
  kb <- as.integer(qcfg$clump_kb %||% 250L)
  clump_p2 <- as.numeric(qcfg$clump_p2 %||% 1.0)
  plink_bin <- qcfg$plink_binary %||% "plink"
  bfile <- config$paths$thousand_g_plink_bfile %||% ""
  mhc_chr <- as.character(config$gwas$mhc_chr %||% "6")
  mhc_start <- as.integer(config$gwas$mhc_start %||% 25000000L)
  mhc_end <- as.integer(config$gwas$mhc_end %||% 35000000L)

  have_plink <- nzchar(bfile) && file.exists(paste0(bfile, ".bim"))
  if (!have_plink) {
    log_warn("qc_filters", sprintf(
      "PLINK bfile not configured/found (%s.bim); LD-buddy clumping disabled. Only exact-Q_SNP hits will be dropped.", bfile))
  }

  summary_rows <- list()
  for (f in factor_names) {
    g <- as.data.table(gwas_raw[[f]])
    snp_col <- find_col(g, c("SNP", "rsid"))
    chr_col <- find_col(g, c("CHR", "chromosome", "chr"))
    bp_col  <- find_col(g, c("BP", "position", "pos"))
    p_col   <- find_col(g, c("Pval_Estimate", "P", "pval", "p_value"))
    q_col   <- find_col(g, c("Q_SNP", "Q"))
    qp_col  <- find_col(g, c("Q_SNP_pval", "Q_pval"))
    sm_col  <- find_col(g, c("smooth_check", "Z_smooth_max"))
    if (any(is.na(c(snp_col, chr_col, bp_col, p_col)))) {
      log_warn("qc_filters", sprintf("factor %s missing core columns; skipping", f))
      next
    }

    base <- data.table(
      SNP = g[[snp_col]],
      CHR = as.character(g[[chr_col]]),
      BP  = as.integer(g[[bp_col]]),
      P   = as.numeric(g[[p_col]]),
      Q_SNP_pval = if (!is.na(qp_col)) as.numeric(g[[qp_col]]) else NA_real_,
      Z_smooth   = if (!is.na(sm_col)) as.numeric(g[[sm_col]]) else NA_real_
    )
    n0 <- nrow(base)
    base <- base[!is.na(SNP) & !is.na(P)]

    qsnp_hits <- base[!is.na(Q_SNP_pval) & Q_SNP_pval < qsnp_thresh, SNP]
    qsnp_drop <- qsnp_hits

    if (have_plink && length(qsnp_hits) > 0L) {
      tsv_path <- file.path(out_dir, sprintf("%s_qsnp_input.tsv", f))
      fwrite(base[!is.na(Q_SNP_pval), .(SNP, Q_SNP_pval)],
             tsv_path, sep = "\t")
      log_info("qc_filters", sprintf(
        "factor %s: %d Q_SNP hits; running plink --clump (r2=%.2f, kb=%d)",
        f, length(qsnp_hits), r2, kb))
      qsnp_drop <- clump_qsnp(tsv_path, bfile,
                              file.path(out_dir, sprintf("%s_qsnp_clump", f)),
                              plink = plink_bin, p1 = qsnp_thresh, p2 = clump_p2,
                              r2 = r2, kb = kb)
      log_info("qc_filters", sprintf(
        "factor %s: clump expanded to %d SNP(s) (index + LD-buddies)",
        f, length(qsnp_drop)))
    } else if (length(qsnp_hits) > 0L) {
      log_info("qc_filters", sprintf("factor %s: %d Q_SNP hits, no clumping (PLINK unavailable)",
                                       f, length(qsnp_hits)))
    }

    zsm_drop <- if (!is.na(sm_col)) base[!is.na(Z_smooth) & abs(Z_smooth) > zsm_thresh, SNP] else character(0)
    mhc_drop <- base[CHR == mhc_chr & BP >= mhc_start & BP <= mhc_end, SNP]

    dropped_all <- unique(c(qsnp_drop, zsm_drop, mhc_drop))
    pass <- base[!SNP %in% dropped_all]

    qc_pass_path <- file.path(out_dir, sprintf("%s_qc_pass.tsv.gz", f))
    gz <- gzfile(qc_pass_path, "wt")
    fwrite(pass, gz, sep = "\t")
    close(gz)

    summary_rows[[f]] <- data.table(
      factor = f,
      n_in = n0,
      n_dropped_qsnp = length(qsnp_drop),
      n_dropped_smooth = length(zsm_drop),
      n_dropped_mhc = length(mhc_drop),
      n_dropped_total = length(dropped_all),
      n_qc_pass = nrow(pass)
    )
    log_info("qc_filters", sprintf(
      "factor %s: in=%s, Q_SNP+LD drop=%d, smooth drop=%d, MHC drop=%d, total drop=%d, pass=%s",
      f, format(n0, big.mark = ","),
      length(qsnp_drop), length(zsm_drop), length(mhc_drop),
      length(dropped_all), format(nrow(pass), big.mark = ",")))
  }

  if (length(summary_rows) > 0L) {
    fwrite(rbindlist(summary_rows), file.path(out_dir, "summary.csv"))
  }
  write_stage_manifest("qc_filters", sex, config, list.files(out_dir, full.names = TRUE))
  invisible(summary_rows)
}
