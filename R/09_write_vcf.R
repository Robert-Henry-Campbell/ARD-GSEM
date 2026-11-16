normalize_gwas_columns <- function(g, snp_sumstats) {
  est_col <- intersect(c("est", "Estimate", "Unstand_Est", "Effect"), names(g))[1]
  se_col  <- intersect(c("SE", "Std_Error", "Unstand_SE"), names(g))[1]
  p_col   <- intersect(c("Pval_Estimate", "P", "pval", "p_value"), names(g))[1]
  if (is.na(est_col) || is.na(se_col) || is.na(p_col)) {
    stop(sprintf(
      "userGWAS output is missing estimate/SE/P column. Have: %s",
      paste(names(g), collapse = ", ")))
  }

  q_col    <- intersect(c("Q_SNP", "Q"), names(g))[1]
  qdf_col  <- intersect(c("Q_SNP_df", "Q_df"), names(g))[1]
  qp_col   <- intersect(c("Q_SNP_pval", "Q_pval"), names(g))[1]
  fail_col <- intersect(c("fail", "failure", "Failed", "Failure"), names(g))[1]

  snp_col <- intersect(c("SNP", "rsid"), names(g))[1]
  chr_col <- intersect(c("CHR", "chromosome", "chr"), names(g))[1]
  bp_col  <- intersect(c("BP", "position", "pos"), names(g))[1]
  a1_col  <- intersect(c("A1", "EA"), names(g))[1]
  a2_col  <- intersect(c("A2", "OA"), names(g))[1]
  maf_col <- intersect(c("MAF", "EAF"), names(g))[1]

  result <- data.table(
    SNP = if (!is.na(snp_col)) g[[snp_col]] else NA_character_,
    CHR = if (!is.na(chr_col)) g[[chr_col]] else NA_character_,
    BP  = if (!is.na(bp_col))  g[[bp_col]]  else NA_integer_,
    A1  = if (!is.na(a1_col))  g[[a1_col]]  else NA_character_,
    A2  = if (!is.na(a2_col))  g[[a2_col]]  else NA_character_,
    AF  = if (!is.na(maf_col)) g[[maf_col]] else NA_real_,
    ES  = g[[est_col]],
    SE  = g[[se_col]],
    P   = g[[p_col]],
    Q_SNP      = if (!is.na(q_col))   g[[q_col]]   else NA_real_,
    Q_SNP_df   = if (!is.na(qdf_col)) g[[qdf_col]] else NA_real_,
    Q_SNP_pval = if (!is.na(qp_col))  g[[qp_col]]  else NA_real_,
    fail       = if (!is.na(fail_col)) as.logical(g[[fail_col]]) else NA
  )

  if (any(is.na(result$CHR) | is.na(result$BP) | is.na(result$A1) | is.na(result$A2)) &&
      !is.null(snp_sumstats)) {
    keys <- intersect(c("SNP", "rsid"), names(snp_sumstats))
    if (length(keys) > 0L) {
      key <- keys[1]
      ss <- as.data.table(snp_sumstats)
      ss_cols <- intersect(c("SNP", "CHR", "BP", "A1", "A2", "MAF"), names(ss))
      ss_sub <- unique(ss[, ..ss_cols])
      result <- merge(result[, !c("CHR","BP","A1","A2","AF"), with = FALSE],
                      ss_sub, by.x = "SNP", by.y = key, all.x = TRUE, sort = FALSE)
      if ("MAF" %in% names(result)) setnames(result, "MAF", "AF")
    }
  }

  result[]
}

format_vcf_field <- function(x, digits = 6) {
  ifelse(is.na(x), ".", formatC(x, format = "g", digits = digits))
}

compute_lambda_gc <- function(p) {
  p <- p[!is.na(p) & p > 0 & p < 1]
  if (length(p) < 10L) return(NA_real_)
  chisq <- qchisq(p, df = 1L, lower.tail = FALSE)
  median(chisq) / qchisq(0.5, df = 1L)
}

bgzip_and_index <- function(plain_vcf, bgz_path) {
  if (!requireNamespace("Rsamtools", quietly = TRUE)) {
    log_warn("vcf", "Rsamtools not installed; falling back to plain gzip (no .tbi index will be created).")
    return(FALSE)
  }
  tryCatch({
    Rsamtools::bgzip(file = plain_vcf, dest = bgz_path, overwrite = TRUE)
    Rsamtools::indexTabix(bgz_path, format = "vcf")
    TRUE
  }, error = function(e) {
    log_warn("vcf", sprintf("bgzip/tabix failed (%s); falling back to plain gzip.",
                             conditionMessage(e)))
    FALSE
  })
}

write_gwas_vcf <- function(gwas_list, snp_sumstats, output_path, metadata) {
  stopifnot(is.list(gwas_list), length(gwas_list) >= 1L)
  factor_names <- names(gwas_list)
  stopifnot(!is.null(factor_names), all(nzchar(factor_names)))

  per_factor <- lapply(gwas_list, normalize_gwas_columns, snp_sumstats = snp_sumstats)
  base <- per_factor[[1]][, .(SNP, CHR, BP, A1, A2, AF)]
  for (i in seq_along(per_factor)) {
    g <- per_factor[[i]]
    if (!identical(g$SNP, base$SNP)) {
      stop(sprintf("factor %s has a different SNP set than the first factor; cannot interleave",
                   factor_names[i]))
    }
  }

  base <- base[!is.na(CHR) & !is.na(BP)]
  ord <- order(base$CHR, base$BP)
  base <- base[ord]
  snp_keep <- base$SNP
  per_factor <- lapply(per_factor, function(g) g[match(snp_keep, g$SNP)])

  lambda_by_factor <- vapply(per_factor, function(g) compute_lambda_gc(g$P),
                              numeric(1))
  names(lambda_by_factor) <- factor_names

  n_fail_by_factor <- vapply(per_factor,
                              function(g) sum(g$fail, na.rm = TRUE),
                              integer(1))
  names(n_fail_by_factor) <- factor_names

  nhat_maf <- metadata$nhat_maf_threshold %||% 0.10
  nhat_by_factor <- vapply(per_factor, function(g) {
    compute_factor_nhat(g$AF, g$SE, maf_threshold = nhat_maf)
  }, numeric(1))
  names(nhat_by_factor) <- factor_names
  if (is.null(metadata$neff_by_factor)) metadata$neff_by_factor <- as.list(nhat_by_factor)

  hdr <- c(
    "##fileformat=VCFv4.2",
    sprintf("##fileDate=%s", format(Sys.Date(), "%Y%m%d")),
    sprintf("##source=ARD-GSEM_pipeline_v1"),
    sprintf("##gsem_scale=log_odds_ratio"),
    sprintf("##gsem_source=Neale_UKBB_round2_linear_regression"),
    sprintf("##gsem_conversion=GenomicSEM_sumstats_linprob_TRUE"),
    sprintf("##gsem_identification=%s",
            if (isTRUE(metadata$std_lv)) "std.lv (factor variance fixed to 1; SNP betas per-SD-of-factor)"
            else "marker-variable (first indicator loading fixed to 1; SNP betas in units of the first indicator)"),
    sprintf("##gsem_model=one_factor_per_ICD10_chapter"),
    sprintf("##gsem_sex=%s", metadata$sex),
    sprintf("##gsem_genome_build=%s", metadata$genome_build %||% "GRCh37"),
    sprintf("##gsem_study_id_prefix=%s", metadata$study_id_prefix %||% "ARD-GSEM"),
    sprintf("##gsem_n_factors=%d", length(factor_names)),
    sprintf("##gsem_nhat_formula=mean(1/((2*MAF*(1-MAF))*SE^2))_at_MAF_geq_%.2f", nhat_maf)
  )
  for (f in factor_names) {
    hdr <- c(hdr,
             sprintf("##gsem_lambda_gc_%s=%s", f,
                     if (is.na(lambda_by_factor[[f]])) "NA" else
                       formatC(lambda_by_factor[[f]], format = "f", digits = 4)),
             sprintf("##gsem_nhat_%s=%s", f,
                     if (is.na(nhat_by_factor[[f]])) "NA" else
                       formatC(nhat_by_factor[[f]], format = "f", digits = 1)))
    if (!is.na(n_fail_by_factor[[f]]) && n_fail_by_factor[[f]] > 0L) {
      hdr <- c(hdr, sprintf("##gsem_n_sem_fail_%s=%d", f, n_fail_by_factor[[f]]))
    }
  }
  hdr <- c(hdr,
    "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Reference panel allele frequency (alt allele)\">",
    "##FORMAT=<ID=ES,Number=A,Type=Float,Description=\"Effect size estimate (log OR per allele copy of ALT)\">",
    "##FORMAT=<ID=SE,Number=A,Type=Float,Description=\"Standard error of effect size estimate\">",
    "##FORMAT=<ID=LP,Number=A,Type=Float,Description=\"-log10 p-value for effect estimate\">",
    "##FORMAT=<ID=AF,Number=A,Type=Float,Description=\"Alternative allele frequency\">",
    "##FORMAT=<ID=SS,Number=A,Type=Float,Description=\"Effective sample size used in factor GWAS\">",
    "##FORMAT=<ID=HET,Number=A,Type=Float,Description=\"Q_SNP heterogeneity chi-square statistic (per-factor)\">",
    "##FORMAT=<ID=HP,Number=A,Type=Float,Description=\"Q_SNP heterogeneity p-value (per-factor); downstream filter on this rather than a variant-level FILTER\">",
    "##FORMAT=<ID=SF,Number=A,Type=Integer,Description=\"Per-factor SEM convergence failure (1 = userGWAS did not converge for this SNP/factor, 0 = converged)\">"
  )
  for (f in factor_names) {
    hdr <- c(hdr, sprintf("##SAMPLE=<ID=%s,Description=\"GSEM chapter factor %s in %s\">",
                          f, f, metadata$sex))
  }
  hdr <- c(hdr,
           paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT",
                   factor_names), collapse = "\t"))

  format_id <- "ES:SE:LP:AF:SS:HET:HP:SF"
  body <- character(nrow(base))
  for (i in seq_len(nrow(base))) {
    chrom <- as.character(base$CHR[i])
    pos   <- as.integer(base$BP[i])
    id    <- as.character(base$SNP[i])
    ref   <- as.character(base$A2[i])
    alt   <- as.character(base$A1[i])
    af    <- format_vcf_field(base$AF[i])
    info  <- sprintf("AF=%s", af)
    sample_vals <- vapply(factor_names, function(f) {
      g <- per_factor[[f]]
      es <- format_vcf_field(g$ES[i])
      se <- format_vcf_field(g$SE[i])
      p  <- g$P[i]
      lp <- if (is.na(p) || p <= 0) "." else format_vcf_field(-log10(p))
      ss <- if (is.null(metadata$neff_by_factor)) "." else format_vcf_field(metadata$neff_by_factor[[f]])
      het <- format_vcf_field(g$Q_SNP[i])
      hp  <- format_vcf_field(g$Q_SNP_pval[i])
      sf <- if (is.na(g$fail[i])) "." else if (isTRUE(g$fail[i])) "1" else "0"
      paste(es, se, lp, af, ss, het, hp, sf, sep = ":")
    }, character(1))
    body[i] <- paste(c(chrom, pos, id, ref, alt, ".", "PASS", info, format_id,
                       sample_vals), collapse = "\t")
  }

  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)

  # Write plain VCF first, then bgzip + tabix index. On failure of either,
  # fall back to plain gzip so the pipeline still ships *a* VCF.
  plain_path <- sub("\\.gz$", "", output_path)
  if (!grepl("\\.vcf$", plain_path)) plain_path <- paste0(plain_path, ".vcf")
  con <- file(plain_path, open = "wt")
  writeLines(c(hdr, body), con)
  close(con)

  bgz_path <- paste0(plain_path, ".gz")
  indexed <- bgzip_and_index(plain_path, bgz_path)
  if (indexed) {
    if (file.exists(plain_path)) file.remove(plain_path)
    final_path <- bgz_path
    tbi_path <- paste0(bgz_path, ".tbi")
  } else {
    # Fallback: plain gzip (no .tbi).
    fallback_con <- gzfile(bgz_path, open = "wt")
    writeLines(c(hdr, body), fallback_con)
    close(fallback_con)
    if (file.exists(plain_path)) file.remove(plain_path)
    final_path <- bgz_path
    tbi_path <- NA_character_
  }

  list(
    path = final_path,
    tbi_path = tbi_path,
    n_variants = nrow(base),
    n_factors = length(factor_names),
    n_sem_fail = n_fail_by_factor,
    lambda_gc = lambda_by_factor,
    nhat_by_factor = nhat_by_factor,
    indexed = indexed
  )
}

# --- GWAS-VCF from the post-GWAS QC artefacts (bifactor path) ------------------------------
# REF/ALT follow GRCh37: the factor GWAS effect allele A1 is usually the reference allele, so
# ES is the ALT-allele effect (-beta when A1 is REF) and AF is the ALT frequency (Pan-UKB EUR).
# One multi-sample VCF (a sample per factor) plus single-sample VCFs per factor (ardmr reads
# only the first sample column).

.vcf_num <- function(x, digits = 6) {
  out <- as.character(signif(x, digits))
  out[is.na(x)] <- "."
  out
}

postgwas_vcf_header <- function(config, sex, factors, reg, fai, sample_desc) {
  pg <- config$postgwas %||% list()
  contigs <- fread(fai, header = FALSE, select = 1:2, col.names = c("id", "len"))
  contigs[, id := sub("^chr", "", id)]
  contigs <- contigs[id %in% c(as.character(1:22), "X")]
  contigs <- contigs[order(match(id, c(as.character(1:22), "X")))]
  hdr <- c(
    "##fileformat=VCFv4.2",
    sprintf("##fileDate=%s", format(Sys.Date(), "%Y%m%d")),
    "##source=ARD-GSEM_post-GWAS_pipeline",
    "##reference=GRCh37 (hg19)",
    sprintf("##contig=<ID=%s,length=%s,assembly=GRCh37>", contigs$id, contigs$len),
    "##gsem_model=bifactor: general factor ARD + orthogonal domain factors (GenomicSEM userGWAS, DWLS, fix_measurement)",
    "##gsem_inputs=Pan-UKB EUR both-sexes ICD-10 GWAS (SAIGE logistic; liability-scaled by GenomicSEM sumstats)",
    "##gsem_scale=std.lv: SNP effects per SD of the latent genetic factor",
    sprintf("##gsem_tiers=full: converged and |Z_smooth|<=%s; qfilter: full minus Q_SNP p<%s and their LD partners (r2>%s within %s kb, 1000G EUR)",
            config$qc_filters$zsmooth_threshold %||% 1.96, pg$qsnp_p %||% 5e-8,
            config$qc_filters$qfilter_r2 %||% 0.1, config$qc_filters$qfilter_kb %||% 1000)
  )
  for (f in names(pg$excluded_factors %||% list())) {
    hdr <- c(hdr, sprintf("##gsem_excluded_factor=%s,reason=\"%s\"", f,
                          gsub("\\s+", " ", pg$excluded_factors[[f]])))
  }
  for (f in factors) {
    r <- reg[[f]]
    hdr <- c(hdr, sprintf(paste0("##gsem_factor=<ID=%s,role=%s,n_indicators=%d,N_hat_std=%.0f,marker=%s,",
                                 "marker_lambda=%.6g,N_hat_marker=%.0f,qsnp_df_gsem=%s,qsnp_df_used=%s,lambda_gc=%.4f>"),
                          f, r$role, length(r$indicators), r$N_hat_std, r$marker, r$marker_lambda,
                          r$N_hat_marker, paste(r$qsnp_df_gsem, collapse = "/"), r$qsnp_df_eff, r$lambda_gc_full))
  }
  c(hdr,
    "##INFO=<ID=AF,Number=A,Type=Float,Description=\"ALT allele frequency, Pan-UKB EUR\">",
    "##INFO=<ID=SI,Number=1,Type=Float,Description=\"Imputation INFO, Pan-UKB\">",
    "##FORMAT=<ID=ES,Number=A,Type=Float,Description=\"Effect of the ALT allele per SD of the latent factor\">",
    "##FORMAT=<ID=SE,Number=A,Type=Float,Description=\"Standard error of ES\">",
    "##FORMAT=<ID=LP,Number=A,Type=Float,Description=\"-log10 p-value of ES\">",
    "##FORMAT=<ID=AF,Number=A,Type=Float,Description=\"ALT allele frequency, Pan-UKB EUR\">",
    "##FORMAT=<ID=SS,Number=A,Type=Float,Description=\"Effective sample size N_hat (GenomicSEM wiki formula, std.lv scale)\">",
    "##FORMAT=<ID=HET,Number=A,Type=Float,Description=\"Q_SNP heterogeneity chi-square\">",
    "##FORMAT=<ID=HP,Number=A,Type=Float,Description=\"Q_SNP p-value at the df used (see gsem_factor qsnp_df_used)\">",
    "##FORMAT=<ID=QF,Number=A,Type=Integer,Description=\"1 = fails the qfilter tier (Q_SNP hit or LD partner of one)\">",
    "##FORMAT=<ID=SF,Number=A,Type=Integer,Description=\"1 = userGWAS fit unusable (no estimate / error / smoothing)\">",
    sprintf("##SAMPLE=<ID=%s,Description=\"%s\">", factors, sample_desc[factors]))
}

write_postgwas_vcf <- function(config, sex, out_dir) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  reg <- read_factor_registry(config, sex)
  factors <- postgwas_factors(config)
  first <- load_postgwas(config, sex, factors[1], "all",
                         cols = c("CHR", "BP", "SNP", "A1", "A2", "a1_is_ref", "EAF_A1", "INFO"))
  keep <- !is.na(first$a1_is_ref)
  base <- first[keep]
  base[, `:=`(CHROM = fifelse(CHR == 23L, "X", as.character(CHR)),
              REF = fifelse(a1_is_ref, A1, A2), ALT = fifelse(a1_is_ref, A2, A1),
              AF = fifelse(a1_is_ref, 1 - EAF_A1, EAF_A1))]
  sgn <- fifelse(base$a1_is_ref, -1, 1)
  info_col <- paste0("AF=", .vcf_num(base$AF), ";SI=", .vcf_num(base$INFO, 4))
  sample_cols <- list()
  for (f in factors) {
    d <- load_postgwas(config, sex, f, "all",
                       cols = c("SNP", "beta", "SE", "LOG10P", "Q_SNP", "Q_SNP_pval", "pass_full", "pass_qfilter"))[keep]
    stopifnot(identical(d$SNP, base$SNP))
    sample_cols[[f]] <- paste(.vcf_num(sgn * d$beta), .vcf_num(d$SE), .vcf_num(d$LOG10P, 5),
                              .vcf_num(base$AF), round(reg[[f]]$N_hat_std),
                              .vcf_num(d$Q_SNP, 5), .vcf_num(d$Q_SNP_pval, 4),
                              as.integer(!d$pass_qfilter), as.integer(!d$pass_full), sep = ":")
  }
  fixed <- data.table(base$CHROM, base$BP, base$SNP, base$REF, base$ALT, ".", "PASS", info_col,
                      "ES:SE:LP:AF:SS:HET:HP:QF:SF")
  desc <- setNames(sprintf("GenomicSEM bifactor factor %s (%s)", factors,
                           vapply(factors, function(f) reg[[f]]$role, character(1))), factors)
  fai <- paste0(config$refs$hg19_fasta, ".fai")
  write_one <- function(path, fs) {
    plain <- sub("\\.gz$", "", path)
    hdr <- postgwas_vcf_header(config, sex, fs, reg, fai, desc)
    writeLines(c(hdr, paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO",
                              "FORMAT", fs), collapse = "\t")), plain)
    fwrite(cbind(fixed, as.data.table(sample_cols[fs])), plain, sep = "\t", quote = FALSE,
           col.names = FALSE, append = TRUE)
    Rsamtools::bgzip(plain, dest = path, overwrite = TRUE)
    unlink(plain)
    Rsamtools::indexTabix(path, format = "vcf")
    path
  }
  basename_vcf <- config$vcf$basename %||% paste0(sex, "_bifactor_gsem")
  paths <- write_one(file.path(out_dir, paste0(basename_vcf, ".vcf.gz")), factors)
  for (f in factors) paths <- c(paths, write_one(file.path(out_dir, sprintf("%s_%s.vcf.gz", basename_vcf, f)), f))
  list(paths = paths, n_variants = nrow(base), n_dropped_unoriented = sum(!keep))
}

run_write_vcf <- function(config, sex) {
  log_info("vcf", sprintf("=== VCF write stage: %s ===", sex))

  if (exists("factor_registry_path", mode = "function") &&
      file.exists(factor_registry_path(config, sex))) {
    out_dir <- stage_dir(config, sex, "vcf")
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    res <- write_postgwas_vcf(config, sex, out_dir)
    log_info("vcf", sprintf("Wrote %d VCF(s), %s variants (%d dropped: REF not resolvable)",
                            length(res$paths), format(res$n_variants, big.mark = ","),
                            res$n_dropped_unoriented))
    write_stage_manifest("vcf", sex, config, c(res$paths, paste0(res$paths, ".tbi")))
    return(invisible(res))
  }

  gwas_dir <- file.path(config$paths$output_dir, sex, "gwas")
  gwas_path <- file.path(gwas_dir, paste0(sex, "_userGWAS_raw.rds"))
  if (!file.exists(gwas_path)) {
    legacy <- file.path(gwas_dir, "userGWAS_raw.rds")
    if (file.exists(legacy)) {
      log_warn("vcf", sprintf("Using legacy (unprefixed) userGWAS_raw at %s (mtime=%s) -- expected %s",
                              legacy, format(file.mtime(legacy)), gwas_path))
      gwas_path <- legacy
    } else {
      log_fatal("vcf", sprintf("userGWAS output not found at %s (run --stage gwas first)", gwas_path))
    }
  }
  gwas_result <- readRDS(gwas_path)

  ss_dir <- file.path(config$paths$output_dir, sex, "sumstats")
  snp_path <- file.path(ss_dir, paste0(sex, "_snp_sumstats.rds"))
  if (!file.exists(snp_path)) {
    legacy <- file.path(ss_dir, "snp_sumstats.rds")
    if (file.exists(legacy)) {
      log_warn("vcf", sprintf("Using legacy (unprefixed) snp_sumstats at %s (mtime=%s)",
                              legacy, format(file.mtime(legacy))))
      snp_path <- legacy
    } else snp_path <- NA_character_
  }
  snp_sumstats <- if (!is.na(snp_path)) readRDS(snp_path) else NULL

  output_path <- file.path(gwas_dir, sprintf("%s_chapter_gsem.vcf.gz", sex))

  metadata <- list(
    sex = sex,
    genome_build = config$gwas$genome_build %||% "GRCh37",
    study_id_prefix = config$gwas$study_id_prefix %||% "ARD-GSEM",
    std_lv = isTRUE(config$gwas$std_lv %||% FALSE),
    nhat_maf_threshold = config$gwas$nhat_maf_threshold %||% 0.10,
    neff_by_factor = NULL
  )

  res <- write_gwas_vcf(gwas_result, snp_sumstats, output_path, metadata)
  log_info("vcf", sprintf("Wrote %s: %s variants x %d factor(s); bgzip+tabix=%s",
                          res$path, format(res$n_variants, big.mark = ","),
                          res$n_factors, res$indexed))
  if (!is.na(res$tbi_path)) {
    log_info("vcf", sprintf("Tabix index: %s", res$tbi_path))
  }
  for (f in names(res$lambda_gc)) {
    nf <- res$n_sem_fail[[f]]
    nh <- res$nhat_by_factor[[f]]
    log_info("vcf", sprintf("  factor %s: lambda_GC = %.4f, N_hat = %s, SEM_fail = %d",
                             f, res$lambda_gc[[f]],
                             if (is.na(nh)) "NA" else format(round(nh), big.mark = ","),
                             if (is.na(nf)) 0L else nf))
  }

  write_stage_manifest("vcf", sex, config,
                       c(res$path, if (!is.na(res$tbi_path)) res$tbi_path))

  res
}
