# QC filters stage: the GenomicSEM authors' post-GWAS QC, applied once to the factor GWAS.
# Writes the canonical per-factor artefacts every later stage reads (output/<sex>/qc_filters/
# <factor>.qc.rds + bgzipped/tabixed <factor>.qc.tsv.gz) and a factor registry.
#
# Tiers (both kept as flags on every SNP):
#   pass_full    converged (estimate present, no lavaan error, not "solution has NOT been
#                found") and |Z_smooth| <= zsmooth_threshold (GenomicSEM wiki: drop SNPs that
#                needed heavy smoothing, e.g. Z > 1.96)
#   pass_qfilter pass_full minus Q_SNP p < postgwas$qsnp_p and their LD partners (wiki: "remove
#                SNPs at some Q_SNP_pval threshold ... along with any SNPs that are in LD with
#                Q_SNP significant results"). Partners: r2 > qfilter_r2 within qfilter_kb in
#                1000G EUR; hits without LD information (chrX, absent from the panel) lose a
#                +-nold_window_kb window instead.
# The excluded factor (Osteoarthritis, under-identified) is not processed; the registry
# records why.

# Per-SNP lavaan messages -> class. "fail" SNPs have no usable estimate; "gradient" (NLMINB
# gradient not ~0, seen when I25's residual sits on its 0.001 bound) and "ok" are kept.
classify_gsem_warnings <- function(warning, fail_regex, benign_regex) {
  w <- as.character(warning)
  cls <- rep("other", length(w))
  cls[is.na(w) | w == "0"] <- "ok"
  for (rx in benign_regex) cls[cls == "other" & grepl(rx, w)] <- "ok"
  cls[grepl("gradient are \\(near\\) zero", w)] <- "gradient"
  cls[grepl(fail_regex, w)] <- "fail"
  cls
}

# Reads the userGWAS output for `keep_factors`: from the per-chunk checkpoints in parallel
# when they exist (fast), else from the combined RDS. Returns list(info, warning, error,
# factors = named list of data.tables).
read_factor_gwas <- function(config, sex, all_factors, keep_factors, workers = 8L) {
  gwas_dir <- file.path(config$paths$output_dir, sex, "gwas")
  cols <- c("est", "SE", "Z_Estimate", "Pval_Estimate", "Q_SNP", "Q_SNP_df", "Q_SNP_pval",
            "Z_smooth")
  idx <- match(keep_factors, all_factors)
  pick <- function(x) {
    if (length(x) != length(all_factors)) stop("userGWAS result has ", length(x),
                                               " elements; model spec has ", length(all_factors))
    first <- x[[1]]
    need <- c("SNP", "CHR", "BP", "MAF", "A1", "A2", cols, "error", "warning")
    miss <- setdiff(need, names(first))
    if (length(miss)) stop("userGWAS columns missing: ", paste(miss, collapse = ", "))
    for (j in idx) if (!identical(x[[j]]$SNP, first$SNP)) stop("SNP order differs across factors")
    list(info = as.data.table(first[, c("SNP", "CHR", "BP", "MAF", "A1", "A2")]),
         warning = as.character(first$warning),
         error = as.character(first$error),
         factors = lapply(x[idx], function(d) as.data.table(d[, cols])))
  }
  chunk_files <- sort(list.files(file.path(gwas_dir, "chunks"),
                                 pattern = sprintf("^%s_chunk_[0-9]+\\.rds$", sex), full.names = TRUE))
  if (length(chunk_files)) {
    log_info("qc_filters", sprintf("Reading %d userGWAS chunks (%d workers)", length(chunk_files), workers))
    parts <- parallel::mclapply(chunk_files, function(f) pick(readRDS(f)), mc.cores = workers)
    bad <- vapply(parts, inherits, logical(1), "try-error")
    if (any(bad)) log_fatal("qc_filters", sprintf("failed to read chunk(s): %s",
                                                  paste(basename(chunk_files[bad]), collapse = ", ")))
    out <- list(info = rbindlist(lapply(parts, `[[`, "info")),
                warning = unlist(lapply(parts, `[[`, "warning"), use.names = FALSE),
                error = unlist(lapply(parts, `[[`, "error"), use.names = FALSE),
                factors = setNames(lapply(seq_along(idx), function(j)
                  rbindlist(lapply(parts, function(p) p$factors[[j]]))), keep_factors))
  } else {
    raw <- file.path(gwas_dir, paste0(sex, "_userGWAS_raw.rds"))
    if (!file.exists(raw)) log_fatal("qc_filters", sprintf("no userGWAS output in %s", gwas_dir))
    log_info("qc_filters", sprintf("Reading %s", raw))
    out <- pick(readRDS(raw))
    names(out$factors) <- keep_factors
  }
  out
}

# Effect-allele (A1) orientation and frequency from the Pan-UKB variant manifest (GRCh37,
# EUR allele frequencies of the GWAS sample itself). GenomicSEM's MAF column is the 1000G
# minor-allele frequency, not A1's frequency.
panukb_allele_annotation <- function(info, path) {
  log_info("qc_filters", sprintf("Joining Pan-UKB variant manifest %s", basename(path)))
  pan <- fread(cmd = paste("zcat", shQuote(path)),
               select = c("chrom", "pos", "ref", "alt", "rsid", "nearest_genes", "info", "af_EUR"),
               colClasses = list(character = c("chrom", "ref", "alt", "rsid", "nearest_genes")),
               showProgress = FALSE)
  pan[, CHR := suppressWarnings(as.integer(fifelse(chrom == "X", "23", chrom)))]
  pan <- pan[!is.na(CHR)]
  setnames(pan, c("pos", "info"), c("BP", "INFO"))
  key <- info[, .(row = .I, CHR = as.integer(CHR), BP = as.integer(BP), A1, A2)]
  m <- pan[key, on = .(CHR, BP), nomatch = NULL, allow.cartesian = TRUE]
  m <- m[(A1 == ref & A2 == alt) | (A1 == alt & A2 == ref)]
  m <- m[!duplicated(row)]
  out <- data.table(a1_is_ref = rep(NA, nrow(info)), EAF_A1 = NA_real_, INFO = NA_real_,
                    nearest_genes = NA_character_, panukb_rsid = NA_character_)
  out[m$row, `:=`(a1_is_ref = m$A1 == m$ref,
                  EAF_A1 = fifelse(m$A1 == m$alt, m$af_EUR, 1 - m$af_EUR),
                  INFO = m$INFO, nearest_genes = m$nearest_genes, panukb_rsid = m$rsid)]
  out
}

# REF base check against the hg19 FASTA for SNPs Pan-UKB did not resolve, plus a seeded
# audit sample. Returns list(a1_is_ref, audit) where audit is the fraction with A1 == REF.
fasta_ref_check <- function(info, a1_is_ref, fasta, n_audit = 20000L, seed = 20261108L) {
  fa <- Rsamtools::FaFile(fasta)
  autos <- which(as.integer(info$CHR) <= 22L)
  todo <- intersect(which(is.na(a1_is_ref)), autos)
  set.seed(seed)
  audit <- sort(sample(intersect(which(!is.na(a1_is_ref)), autos), min(n_audit, length(autos))))
  base_at <- function(rows) {
    if (!length(rows)) return(character(0))
    gr <- GenomicRanges::GRanges(paste0("chr", info$CHR[rows]),
                                 IRanges::IRanges(as.integer(info$BP[rows]), width = 1L))
    toupper(as.character(Rsamtools::getSeq(fa, gr)))
  }
  b <- base_at(todo)
  a1_is_ref[todo] <- fifelse(b == info$A1[todo], TRUE, fifelse(b == info$A2[todo], FALSE, NA))
  ba <- base_at(audit)
  list(a1_is_ref = a1_is_ref, n_fasta_resolved = sum(!is.na(a1_is_ref[todo])),
       audit_frac_a1_ref = mean(ba == info$A1[audit]), audit_concordance = mean((ba == info$A1[audit]) == a1_is_ref[audit]))
}

# Loadings of the full-sample bifactor (cfa/model_bifactor_full.rds) -> marker per factor.
marker_loadings <- function(config, sex) {
  fit <- readRDS(file.path(config$paths$output_dir, sex, "cfa", "model_bifactor_full.rds"))
  r <- as.data.table(fit$results)[op == "=~"]
  r[, .(factor = lhs, indicator = rhs, lambda = as.numeric(Unstand_Est),
        std = as.numeric(STD_All))]
}

choose_marker <- function(load_tab, fac, rule = "max_std") {
  lt <- load_tab[factor == fac]
  if (!nrow(lt)) stop("no loadings for factor ", fac)
  if (identical(rule, "max_std")) lt[which.max(std)]
  else if (rule %in% lt$indicator) lt[indicator == rule]
  else stop("marker_indicator must be 'max_std' or an indicator of ", fac)
}

# Flag Q_SNP LD partners (exact r2 in 1000G EUR) and window-flag hits that have no LD info.
flag_qsnp_partners <- function(config, dt, hits, bfile, r2, kb, nold_kb, prefix) {
  in_panel <- dt$in_1000g & dt$CHR <= 22L
  hit_rows <- which(dt$SNP %in% hits)
  ld_hits <- dt$SNP[hit_rows[in_panel[hit_rows]]]
  nold_rows <- hit_rows[!in_panel[hit_rows]]
  partner <- rep(FALSE, nrow(dt)); method <- rep(NA_character_, nrow(dt))
  if (length(ld_hits)) {
    ld <- plink_ld_partners(config, bfile, ld_hits, r2, kb, prefix, "qc_filters",
                            threads = as.integer(config$postgwas$workers %||% 8L))
    p_rows <- which(dt$SNP %in% ld$SNP_B)
    partner[p_rows] <- TRUE; method[p_rows] <- "r2"
  }
  if (length(nold_rows)) {
    w <- data.table(CHR = dt$CHR[nold_rows], lo = dt$BP[nold_rows] - nold_kb * 1000,
                    hi = dt$BP[nold_rows] + nold_kb * 1000)
    pos <- dt[, .(row = .I, CHR, BP)]
    hit_w <- pos[w, on = .(CHR, BP >= lo, BP <= hi), nomatch = NULL, .(row = x.row)]$row
    hit_w <- setdiff(unique(hit_w), which(partner))
    partner[hit_w] <- TRUE; method[hit_w] <- "window"
  }
  is_hit <- seq_len(nrow(dt)) %in% hit_rows
  partner[is_hit] <- FALSE; method[is_hit] <- NA_character_
  list(qsnp_ld = partner, qsnp_ld_method = method, n_ld_hits = length(ld_hits),
       n_window_hits = length(nold_rows))
}

run_qc_filters <- function(config, sex) {
  stage <- "qc_filters"
  log_info(stage, sprintf("=== Post-GWAS QC filters: %s ===", sex))
  pg <- config$postgwas %||% list()
  qc <- config$qc_filters %||% list()
  out_dir <- stage_dir(config, sex, stage)
  work <- file.path(out_dir, "work")
  dir.create(work, recursive = TRUE, showWarnings = FALSE)
  workers <- as.integer(pg$workers %||% 8L)
  zthr <- as.numeric(qc$zsmooth_threshold %||% 1.96)
  qsnp_p <- as.numeric(pg$qsnp_p %||% 5e-8)
  maf_nhat <- as.numeric(pg$nhat_maf_threshold %||% 0.10)
  bfile <- config$refs$ld_1000g_bfile
  if (!file.exists(paste0(bfile, ".bim"))) log_fatal(stage, sprintf("1000G bfile missing: %s.bim", bfile))

  spec <- readRDS(file.path(config$paths$output_dir, sex, "gwas", paste0(sex, "_model_spec.rds")))
  factors <- postgwas_factors(config)
  miss <- setdiff(factors, spec$factors)
  if (length(miss)) log_fatal(stage, sprintf("factor(s) not in the GWAS model: %s", paste(miss, collapse = ", ")))
  excluded <- pg$excluded_factors %||% list()
  for (f in names(excluded)) log_info(stage, sprintf("Excluded factor %s: %s", f, excluded[[f]]))

  gw <- read_factor_gwas(config, sex, spec$factors, factors, workers)
  info <- gw$info
  info[, `:=`(CHR = as.integer(CHR), BP = as.integer(BP))]
  n <- nrow(info)
  log_info(stage, sprintf("%s SNPs x %d factor(s): %s", format(n, big.mark = ","), length(factors),
                          paste(factors, collapse = ", ")))

  # --- SNP-level annotation shared by all factors ---
  cls <- classify_gsem_warnings(gw$warning, qc$fail_warning_regex %||% "solution has NOT been found",
                                unlist(qc$benign_warning_regex))
  warn_tab <- data.table(class = cls)[, .N, by = class]
  log_info(stage, sprintf("lavaan message classes: %s",
                          paste(sprintf("%s=%s", warn_tab$class, format(warn_tab$N, big.mark = ",")), collapse = ", ")))
  if (any(cls == "other")) {
    ex <- unique(gw$warning[cls == "other"])[1:min(3, sum(cls == "other"))]
    log_warn(stage, sprintf("%d SNP(s) with an unrecognised lavaan message, e.g.: %s",
                            sum(cls == "other"), paste(substr(ex, 1, 160), collapse = " | ")))
  }
  err <- !(is.na(gw$error) | gw$error == "0")

  alle <- panukb_allele_annotation(info, config$refs$panukb_variant_qc)
  log_info(stage, sprintf("Pan-UKB allele match: %.2f%% of SNPs", 100 * mean(!is.na(alle$a1_is_ref))))
  fr <- fasta_ref_check(info, alle$a1_is_ref, config$refs$hg19_fasta)
  alle$a1_is_ref <- fr$a1_is_ref
  log_info(stage, sprintf("hg19 audit: A1 is REF for %.2f%% of a 20k sample; Pan-UKB vs FASTA concordance %.4f; FASTA resolved %d unmatched SNPs",
                          100 * fr$audit_frac_a1_ref, fr$audit_concordance, fr$n_fasta_resolved))

  bim <- fread(paste0(bfile, ".bim"), header = FALSE, select = 2L, col.names = "SNP")
  in_1000g <- info$SNP %in% bim$SNP
  rm(bim); gc(verbose = FALSE)
  is_mhc <- in_mhc(info$CHR, info$BP, config)

  annot <- cbind(info[, .(SNP, CHR, BP, A1, A2, MAF)], alle, in_1000g = in_1000g)
  saveRDS(annot, file.path(out_dir, "variant_annot.rds"), compress = FALSE)

  load_tab <- marker_loadings(config, sex)
  registry <- list(); summary_rows <- list(); calib_rows <- list(); sens_rows <- list()

  for (f in factors) {
    log_info(stage, sprintf("--- factor %s ---", f))
    g <- gw$factors[[f]]
    dt <- data.table(SNP = info$SNP, CHR = info$CHR, BP = info$BP, A1 = info$A1, A2 = info$A2,
                     a1_is_ref = alle$a1_is_ref, MAF = info$MAF, EAF_A1 = alle$EAF_A1, INFO = alle$INFO,
                     beta = g$est, SE = suppressWarnings(as.numeric(g$SE)), Z = g$Z_Estimate,
                     P = g$Pval_Estimate, Q_SNP = g$Q_SNP, Q_SNP_df = g$Q_SNP_df,
                     Q_SNP_pval_gsem = g$Q_SNP_pval, Z_smooth = g$Z_smooth, warn_class = cls,
                     flag_error = err, is_mhc = is_mhc, is_chrX = info$CHR == 23L, in_1000g = in_1000g)
    dt[, LOG10P := log10p_from_z(Z)]
    dt[, fail_fit := is.na(beta) | is.na(SE) | warn_class == "fail"]
    dt[, flag_smooth := !is.na(Z_smooth) & abs(Z_smooth) > zthr]
    dt[, pass_full := !fail_fit & !flag_error & !flag_smooth & SE > 0]
    dt[is.na(pass_full), pass_full := FALSE]

    df_gsem <- unique(na.omit(dt$Q_SNP_df))
    df_eff <- as.numeric(pg$qsnp_df_override[[f]] %||% df_gsem)
    dt[, Q_SNP_df_eff := df_eff]
    dt[, Q_SNP_pval := pchisq(Q_SNP, Q_SNP_df_eff, lower.tail = FALSE)]
    dt[, Q_LOG10P := log10p_chisq(Q_SNP, Q_SNP_df_eff)]

    full <- dt$pass_full
    n_hat_std <- compute_factor_nhat(dt$MAF[full], dt$SE[full], maf_threshold = maf_nhat)
    mk <- choose_marker(load_tab, f, pg$marker_indicator %||% "max_std")
    dt[, `:=`(beta_marker = mk$lambda * beta, SE_marker = mk$lambda * SE)]
    n_hat_marker <- n_hat_std / mk$lambda^2
    sens_rows[[f]] <- load_tab[factor == f, .(factor, indicator, lambda, std,
                                             N_hat_marker = n_hat_std / lambda^2)]
    log_info(stage, sprintf("N_hat_std = %.0f; marker %s (lambda %.4f, std %.3f) -> N_hat_marker = %.0f",
                            n_hat_std, mk$indicator, mk$lambda, mk$std, n_hat_marker))

    hits <- dt[full & !is.na(Q_SNP_pval) & Q_SNP_pval < qsnp_p, SNP]
    dt[, qsnp_hit := SNP %in% hits]
    fl <- flag_qsnp_partners(config, dt, hits, bfile, as.numeric(qc$qfilter_r2 %||% 0.1),
                             as.integer(qc$qfilter_kb %||% 1000L), as.integer(qc$nold_window_kb %||% 500L),
                             file.path(work, paste0(f, "_qsnp")))
    dt[, `:=`(qsnp_ld = fl$qsnp_ld, qsnp_ld_method = fl$qsnp_ld_method)]
    dt[, pass_qfilter := pass_full & !qsnp_hit & !qsnp_ld]
    log_info(stage, sprintf("Q_SNP (df %s; GenomicSEM df %s): %s hits, %s LD partners (%d hits by r2, %d by window); qfilter keeps %s of %s",
                            df_eff, paste(df_gsem, collapse = "/"), format(length(hits), big.mark = ","),
                            format(sum(fl$qsnp_ld), big.mark = ","), fl$n_ld_hits, fl$n_window_hits,
                            format(sum(dt$pass_qfilter), big.mark = ","), format(sum(full), big.mark = ",")))

    setcolorder(dt, c("CHR", "BP", "SNP", "A1", "A2", "a1_is_ref", "MAF", "EAF_A1", "INFO",
                      "beta", "SE", "Z", "P", "LOG10P", "beta_marker", "SE_marker",
                      "Q_SNP", "Q_SNP_df", "Q_SNP_df_eff", "Q_SNP_pval_gsem", "Q_SNP_pval", "Q_LOG10P",
                      "Z_smooth", "warn_class", "fail_fit", "flag_error", "flag_smooth",
                      "qsnp_hit", "qsnp_ld", "qsnp_ld_method", "is_mhc", "is_chrX", "in_1000g",
                      "pass_full", "pass_qfilter"))
    setorder(dt, CHR, BP)
    saveRDS(dt, qc_rds_path(config, sex, f), compress = FALSE)
    write_bgz_tsv(dt, file.path(out_dir, paste0(f, ".qc.tsv.gz")))
    fwrite(dt[qsnp_hit == TRUE, .(CHR, BP, SNP, A1, A2, beta, SE, P, Q_SNP, Q_SNP_df_eff, Q_SNP_pval)],
           file.path(out_dir, paste0(f, "_qsnp_hits.tsv")), sep = "\t")

    q_full <- dt$Q_SNP[full]
    lam <- function(df) median(q_full, na.rm = TRUE) / qchisq(0.5, df)
    calib_rows[[f]] <- data.table(factor = f, df_gsem = paste(df_gsem, collapse = "/"), df_eff = df_eff,
                                  median_Q = median(q_full, na.rm = TRUE),
                                  lambda_Q_gsem = lam(df_gsem[1]), lambda_Q_eff = lam(df_eff),
                                  n_hits_gsem = sum(full & dt$Q_SNP_pval_gsem < qsnp_p, na.rm = TRUE),
                                  n_hits_eff = length(hits))
    gws <- full & dt$P < as.numeric(pg$p_gws %||% 5e-8)
    summary_rows[[f]] <- data.table(
      factor = f, n_snps = n, n_fail_fit = sum(dt$fail_fit), n_error = sum(dt$flag_error),
      n_gradient_warning = sum(dt$warn_class == "gradient"), n_smooth = sum(dt$flag_smooth),
      n_full = sum(full), n_qsnp_hit = length(hits), n_qsnp_ld = sum(dt$qsnp_ld),
      n_qfilter = sum(dt$pass_qfilter), lambda_gc_full = compute_lambda_gc(dt$P[full]),
      mean_chisq_full = mean(dt$Z[full]^2), n_gws_full = sum(gws),
      n_gws_qfilter = sum(gws & dt$pass_qfilter), N_hat_std = n_hat_std, marker = mk$indicator,
      marker_lambda = mk$lambda, N_hat_marker = n_hat_marker)
    registry[[f]] <- list(role = if (f == spec$general_name %||% "ARD") "general" else "specific",
                          identified = TRUE, excluded = FALSE,
                          indicators = if (f == (spec$general_name %||% "ARD")) spec$general_traits else spec$domains[[f]],
                          std_lv = isTRUE(pg$std_lv %||% TRUE), marker = mk$indicator,
                          marker_lambda = mk$lambda, marker_std = mk$std,
                          N_hat_std = n_hat_std, N_hat_marker = n_hat_marker,
                          qsnp_df_gsem = df_gsem, qsnp_df_eff = df_eff,
                          lambda_gc_full = summary_rows[[f]]$lambda_gc_full,
                          n_full = sum(full), n_qfilter = sum(dt$pass_qfilter))
    rm(dt, g); gc(verbose = FALSE)
  }
  for (f in names(excluded)) {
    registry[[f]] <- list(role = "specific", identified = FALSE, excluded = TRUE,
                          indicators = spec$domains[[f]] %||% character(0), reason = excluded[[f]])
  }

  fwrite(rbindlist(summary_rows), file.path(out_dir, "qc_summary.csv"))
  fwrite(warn_tab, file.path(out_dir, "warning_classes.csv"))
  fwrite(rbindlist(calib_rows), file.path(out_dir, "qsnp_calibration.csv"))
  fwrite(rbindlist(sens_rows), file.path(out_dir, "nhat_marker_sensitivity.csv"))
  jsonlite::write_json(list(audit_frac_a1_ref = fr$audit_frac_a1_ref, audit_concordance = fr$audit_concordance,
                            panukb_match = mean(!is.na(alle$a1_is_ref))),
                       file.path(out_dir, "allele_audit.json"), auto_unbox = TRUE, pretty = TRUE)
  jsonlite::write_json(registry, factor_registry_path(config, sex), auto_unbox = TRUE, pretty = TRUE, digits = NA)
  clear_postgwas_cache()

  write_stage_manifest(stage, sex, config, list.files(out_dir, full.names = TRUE))
  invisible(registry)
}
