normalize_gwas_columns <- function(g, snp_sumstats) {
  est_col <- intersect(c("est", "Estimate", "Unstand_Est", "Effect"), names(g))[1]
  se_col  <- intersect(c("SE", "Std_Error", "Unstand_SE"), names(g))[1]
  p_col   <- intersect(c("Pval_Estimate", "P", "pval", "p_value"), names(g))[1]
  if (is.na(est_col) || is.na(se_col) || is.na(p_col)) {
    stop(sprintf(
      "userGWAS output is missing estimate/SE/P column. Have: %s",
      paste(names(g), collapse = ", ")))
  }

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
    P   = g[[p_col]]
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
  out <- ifelse(is.na(x), ".", formatC(x, format = "g", digits = digits))
  out
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

  hdr <- c(
    "##fileformat=VCFv4.2",
    sprintf("##fileDate=%s", format(Sys.Date(), "%Y%m%d")),
    sprintf("##source=ARD-GSEM_pipeline_v1"),
    sprintf("##gsem_scale=log_odds_ratio"),
    sprintf("##gsem_source=Neale_UKBB_round2_linear_regression_converted"),
    sprintf("##gsem_conversion=beta_logOR_eq_beta_linear_div_K_times_1_minus_K"),
    sprintf("##gsem_model=one_factor_per_ICD10_chapter"),
    sprintf("##gsem_sex=%s", metadata$sex),
    sprintf("##gsem_genome_build=%s", metadata$genome_build %||% "GRCh37"),
    sprintf("##gsem_study_id_prefix=%s", metadata$study_id_prefix %||% "ARD-GSEM"),
    sprintf("##gsem_n_factors=%d", length(factor_names)),
    "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Reference panel allele frequency (alt allele)\">",
    "##FORMAT=<ID=ES,Number=A,Type=Float,Description=\"Effect size estimate (log OR per allele copy of ALT)\">",
    "##FORMAT=<ID=SE,Number=A,Type=Float,Description=\"Standard error of effect size estimate\">",
    "##FORMAT=<ID=LP,Number=A,Type=Float,Description=\"-log10 p-value for effect estimate\">",
    "##FORMAT=<ID=AF,Number=A,Type=Float,Description=\"Alternative allele frequency\">",
    "##FORMAT=<ID=SS,Number=A,Type=Float,Description=\"Effective sample size used in factor GWAS\">"
  )
  for (f in factor_names) {
    hdr <- c(hdr, sprintf("##SAMPLE=<ID=%s,Description=\"GSEM chapter factor %s in %s\">",
                          f, f, metadata$sex))
  }
  hdr <- c(hdr,
           paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT",
                   factor_names), collapse = "\t"))

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
      paste(es, se, lp, af, ss, sep = ":")
    }, character(1))
    body[i] <- paste(c(chrom, pos, id, ref, alt, ".", "PASS", info, "ES:SE:LP:AF:SS",
                       sample_vals), collapse = "\t")
  }

  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  con <- gzfile(output_path, open = "wt")
  on.exit(close(con), add = TRUE)
  writeLines(c(hdr, body), con)

  list(path = output_path, n_variants = nrow(base), n_factors = length(factor_names))
}

run_write_vcf <- function(config, sex) {
  log_info("vcf", sprintf("=== VCF write stage: %s ===", sex))

  gwas_path <- file.path(config$paths$output_dir, sex, "gwas", "userGWAS_raw.rds")
  if (!file.exists(gwas_path)) {
    log_fatal("vcf", sprintf("userGWAS_raw.rds not found at %s (run --stage gwas first)", gwas_path))
  }
  gwas_result <- readRDS(gwas_path)

  snp_path <- file.path(config$paths$output_dir, sex, "sumstats", "snp_sumstats.rds")
  snp_sumstats <- if (file.exists(snp_path)) readRDS(snp_path) else NULL

  scale_meta_path <- file.path(config$paths$output_dir, sex, "sumstats", "scale_metadata.json")
  scale_meta <- if (file.exists(scale_meta_path)) jsonlite::read_json(scale_meta_path) else list()

  output_path <- file.path(config$paths$output_dir, sex, "gwas",
                           sprintf("%s_chapter_gsem.vcf.gz", sex))

  metadata <- list(
    sex = sex,
    genome_build = config$gwas$genome_build %||% "GRCh37",
    study_id_prefix = config$gwas$study_id_prefix %||% "ARD-GSEM",
    neff_by_factor = NULL
  )

  res <- write_gwas_vcf(gwas_result, snp_sumstats, output_path, metadata)
  log_info("vcf", sprintf("Wrote %s: %s variants x %d factor(s)",
                          res$path, format(res$n_variants, big.mark = ","), res$n_factors))

  write_stage_manifest("vcf", sex, config, output_path)

  res
}
