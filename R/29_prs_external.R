# External (comparator) GWAS summary statistics for the prs stage: each row of the manifest
# extdata/comparator_sumstats.tsv (config prs$external_manifest) becomes a source with the same
# input table as a factor (SNP CHR BP A1 A2 EAF_A1 beta SE P LOG10P; per-allele effects of A1,
# oriented so that a higher score means worse ageing / more disease), so PRS_METHODS apply
# unchanged. Sample size per the manifest `n` column: a number, "col" (median of the file's N
# column) or "nhat" (GenomicSEM factors: compute_factor_nhat on MAF and SE, as for our factors).

read_external_manifest <- function(config) {
  p <- config$prs$external_manifest %||% ""
  if (!nzchar(p)) return(NULL)
  if (!is_absolute_path(p)) p <- file.path(config$project$root, p)
  if (!file.exists(p)) log_fatal("prs", sprintf("external manifest not found: %s", p))
  m <- fread(p, sep = "\t", colClasses = "character", na.strings = c("", "NA"))
  m[, flip_sign := toupper(flip_sign) %in% c("TRUE", "T", "1", "YES")]
  abs <- vapply(m$file, is_absolute_path, logical(1))
  m[, file := ifelse(abs, file, file.path(dirname(p), file))]
  m[]
}

# The pipeline's variant universe with Pan-UKB EUR frequencies of A1 (rsID assignment for
# chr:pos identifiers and frequencies for sources without an EAF column).
external_reference <- function(config, sex) {
  f <- postgwas_factors(config)[1]
  load_postgwas(config, sex, f, "all", cols = c("SNP", "CHR", "BP", "A1", "A2", "EAF_A1"))
}

# One manifest row -> list(d = standard input table, n = sample size, summary = one-row table).
read_external_sumstats <- function(spec, ref) {
  stage <- "prs"
  sep <- switch(spec$sep %||% "auto", tab = "\t", comma = ",", space = " ", "auto")
  col <- function(nm) { v <- spec[[paste0("col_", nm)]]; if (is.null(v) || is.na(v) || !nzchar(v)) NULL else v }
  want <- unique(unlist(lapply(c("snp", "chr", "bp", "a1", "a2", "eaf", "beta", "se", "p", "mlogp", "n"), col)))
  d <- fread(spec$file, sep = sep, select = want, colClasses = "character", na.strings = c("", "NA", "."))
  out <- data.table(SNP = d[[col("snp")]], CHR = d[[col("chr")]], BP = d[[col("bp")]],
                    A1 = toupper(d[[col("a1")]]), A2 = toupper(d[[col("a2")]]),
                    beta = as.numeric(d[[col("beta")]]), SE = as.numeric(d[[col("se")]]))
  out[, CHR := suppressWarnings(as.integer(sub("^chr", "", CHR, ignore.case = TRUE)))]
  out[, BP := suppressWarnings(as.integer(BP))]
  out[, EAF_A1 := if (!is.null(col("eaf"))) as.numeric(d[[col("eaf")]]) else NA_real_]
  if (!is.null(col("p"))) {
    out[, P := as.numeric(d[[col("p")]])]
    out[, LOG10P := -log10(P)]
  } else if (!is.null(col("mlogp"))) {
    out[, LOG10P := as.numeric(d[[col("mlogp")]])]
    out[, P := 10^-LOG10P]
  } else {
    out[, LOG10P := log10p_from_z(beta / SE)]
    out[, P := 10^-LOG10P]
  }
  out[, N_snp := if (identical(spec$n, "col") && !is.null(col("n"))) as.numeric(d[[col("n")]]) else NA_real_]
  rm(d)
  n0 <- nrow(out)
  out <- out[CHR %in% 1:22 & !is.na(BP) & A1 %in% c("A", "C", "G", "T") & A2 %in% c("A", "C", "G", "T") &
             A1 != A2 & is.finite(beta) & is.finite(SE) & SE > 0 & is.finite(P)]
  # chr:pos identifiers (or missing rsIDs) -> rsID from the reference by position + alleles
  no_rs <- !grepl("^rs[0-9]+$", out$SNP)
  if (any(no_rs)) {
    r <- ref[, .(CHR, BP, rsid = SNP, rA1 = A1, rA2 = A2)]
    m <- merge(out[no_rs, .(i = which(no_rs), CHR, BP, A1, A2)], r, by = c("CHR", "BP"), allow.cartesian = TRUE)
    m <- m[(A1 == rA1 & A2 == rA2) | (A1 == rA2 & A2 == rA1)][!duplicated(i)]
    out[m$i, SNP := m$rsid]
    out <- out[grepl("^rs[0-9]+$", SNP)]
  }
  out <- out[SNP %in% ref$SNP][!duplicated(SNP)]
  if (isTRUE(spec$flip_sign)) out[, beta := -beta]
  # frequencies: the file's EAF, else the reference frequency aligned to this A1
  eaf_source <- "file"
  if (all(is.na(out$EAF_A1))) {
    f <- ref[match(out$SNP, SNP)]
    same <- out$A1 == f$A1 & out$A2 == f$A2
    swap <- out$A1 == f$A2 & out$A2 == f$A1
    out[, EAF_A1 := fifelse(same, f$EAF_A1, fifelse(swap, 1 - f$EAF_A1, NA_real_))]
    eaf_source <- "panukb_eur"
  }
  out <- out[!is.na(EAF_A1)]
  n_rule <- spec$n %||% "nhat"
  n <- switch(n_rule,
              col = stats::median(out$N_snp, na.rm = TRUE),
              nhat = compute_factor_nhat(pmin(out$EAF_A1, 1 - out$EAF_A1), out$SE, maf_threshold = 0.10),
              suppressWarnings(as.numeric(n_rule)))
  if (!is.finite(n)) log_fatal(stage, sprintf("%s: sample size could not be determined (n = %s)", spec$id, n_rule))
  log_info(stage, sprintf("%s: %s rows -> %s usable SNPs in the reference; N = %s%s; sign %s",
                          spec$id, format(n0, big.mark = ","), format(nrow(out), big.mark = ","),
                          format(round(n), big.mark = ","),
                          if (identical(n_rule, "nhat")) " (N_hat from MAF and SE)" else "",
                          if (isTRUE(spec$flip_sign)) "flipped so higher = worse" else "kept (higher = worse)"))
  list(d = out[, .(SNP, CHR, BP, A1, A2, EAF_A1, beta, SE, P, LOG10P)], n = n,
       summary = data.table(id = spec$id, n_rows = n0, n_usable = nrow(out), n_used = round(n), n_rule = n_rule,
                            flip_sign = isTRUE(spec$flip_sign), eaf_source = eaf_source))
}
