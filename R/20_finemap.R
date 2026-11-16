# Fine-mapping stage: SuSiE-RSS (susieR) on every fine-mapping region, with UK Biobank
# in-sample LD from the ld_panel stage (GWASTutorial "Fine-mapping by susieR").
#  - SNPs: full-QC tier (no holes in LD blocks) present in the LD panel; SNPs whose panel A1
#    frequency differs from Pan-UKB's by more than max_af_diff are dropped (allele mismatch)
#  - susie_rss(z, R, n = N_hat_std, L, coverage 0.95, min_abs_corr 0.5, refine) plus an L = 1
#    sensitivity fit; estimate_s_rss() and kriging_rss() diagnose LD mismatch, and SNPs flagged
#    by kriging (logLR > 2 and |z| > 2) are removed and the region refitted (at most twice)
#  - each 95% credible set gets a PIP lead (highest PIP); a set is "qflag" when any member fails
#    the Q_SNP filter, and carries its locus's heterogeneous label

finemap_one <- function(z, R, n, cfg) {
  L <- as.integer(cfg$L %||% 10L)
  fit <- susieR::susie_rss(z = z, R = R, n = n, L = L, coverage = as.numeric(cfg$coverage %||% 0.95),
                           min_abs_corr = as.numeric(cfg$min_abs_corr %||% 0.5),
                           refine = isTRUE(cfg$refine %||% TRUE), max_iter = as.integer(cfg$max_iter %||% 1000L),
                           estimate_residual_variance = FALSE)
  fit
}

FINEMAP_COLS <- c("SNP", "CHR", "BP", "A1", "A2", "EAF_A1", "beta", "SE", "Z", "LOG10P", "pass_qfilter", "qsnp_hit")

# g_all: the factor's full-QC tier (FINEMAP_COLS), loaded once in the parent before forking
run_finemap_region <- function(config, sex, reg_row, cfg, reg_info, out_dir, g_all) {
  f <- reg_row$factor
  rdir <- file.path(out_dir, f, reg_row$ld_id)
  done <- file.path(rdir, "done")
  if (file.exists(done)) return(readRDS(file.path(rdir, "summary.rds")))
  dir.create(rdir, recursive = TRUE, showWarnings = FALSE)
  ld <- read_region_ld(config, sex, reg_row$ld_id)
  g <- g_all[CHR == reg_row$chr & BP >= reg_row$start & BP <= reg_row$end]
  n_region <- nrow(g)
  top_snp <- g$SNP[which.max(g$LOG10P)]
  i <- match(g$SNP, ld$snps)
  g <- g[!is.na(i)]; i <- i[!is.na(i)]
  af_ok <- abs(ld$eaf_panel_a1[i] - g$EAF_A1) <= as.numeric(cfg$max_af_diff %||% 0.1)
  af_ok[is.na(af_ok)] <- FALSE
  n_af_drop <- sum(!af_ok)
  g <- g[af_ok]; i <- i[af_ok]
  R <- ld$R[i, i, drop = FALSE]
  bad_diag <- !is.finite(diag(R)) | rowSums(!is.finite(R)) > 0
  g <- g[!bad_diag]; R <- R[!bad_diag, !bad_diag, drop = FALSE]
  z <- setNames(g$Z, g$SNP)
  n <- reg_info$N_hat_std

  kriged <- character(0); s_est <- NA_real_
  for (round in 0:2) {
    # one eigendecomposition per round, shared by estimate_s_rss and kriging_rss
    # (subsetting R below drops the attribute, so it is recomputed after removals)
    attr(R, "eigen") <- eigen(R, symmetric = TRUE)
    s_est <- susieR::estimate_s_rss(z, R, n = n)
    kr <- susieR::kriging_rss(z, R, n = n, s = s_est)$conditional_dist
    out <- which(kr$logLR > 2 & abs(kr$z) > 2)
    if (!length(out) || round == 2) break
    kriged <- c(kriged, names(z)[out])
    keep <- -out
    z <- z[keep]; R <- R[keep, keep, drop = FALSE]; g <- g[keep]
  }
  attr(R, "eigen") <- NULL
  fit <- tryCatch(finemap_one(z, R, n, cfg), error = function(e) e)
  if (inherits(fit, "error")) {
    summ <- data.table(factor = f, fm_id = reg_row$fm_id, n_region = n_region, n_used = length(z),
                       error = conditionMessage(fit))
    saveRDS(summ, file.path(rdir, "summary.rds")); file.create(done)
    return(summ)
  }
  fit1 <- tryCatch(finemap_one(z, R, n, modifyList(cfg, list(L = 1L, refine = FALSE))), error = function(e) NULL)
  cs <- susieR::susie_get_cs(fit, coverage = as.numeric(cfg$coverage %||% 0.95),
                             min_abs_corr = as.numeric(cfg$min_abs_corr %||% 0.5), Xcorr = R)
  pip <- susieR::susie_get_pip(fit)
  pt <- data.table(SNP = names(z), CHR = g$CHR, BP = g$BP, A1 = g$A1, A2 = g$A2, z = unname(z),
                   LOG10P = g$LOG10P, pip = pip, pip_L1 = if (is.null(fit1)) NA_real_ else susieR::susie_get_pip(fit1),
                   cs = NA_integer_, pass_qfilter = g$pass_qfilter)
  cs_rows <- list()
  if (length(cs$cs)) {
    for (k in seq_along(cs$cs)) {
      mem <- cs$cs[[k]]
      pt[mem, cs := k]
      lead <- mem[which.max(pip[mem])]
      cs_rows[[k]] <- data.table(factor = f, fm_id = reg_row$fm_id, locus_id = reg_row$locus_id,
                                 cs = k, cs_id = sprintf("%s_cs%d", reg_row$fm_id, k),
                                 susie_l = as.integer(sub("^L", "", names(cs$cs)[k])), size = length(mem),
                                 coverage = cs$coverage[k], purity_min = cs$purity$min.abs.corr[k],
                                 purity_mean = cs$purity$mean.abs.corr[k],
                                 log10bf = fit$lbf[as.integer(sub("^L", "", names(cs$cs)[k]))] / log(10),
                                 lead_snp = names(z)[lead], lead_chr = g$CHR[lead], lead_bp = g$BP[lead],
                                 lead_pip = pip[lead], lead_z = unname(z[lead]), lead_log10p = g$LOG10P[lead],
                                 qflag = any(!g$pass_qfilter[mem]),
                                 members = paste(names(z)[mem], collapse = ","))
    }
  }
  csd <- rbindlist(cs_rows)
  fit$lbf_variable_snps <- names(z)
  saveRDS(fit, file.path(rdir, "fit.rds"))
  fwrite(pt, file.path(rdir, "pip.tsv.gz"), sep = "\t")
  fwrite(csd, file.path(rdir, "cs.tsv"), sep = "\t")
  summ <- data.table(factor = f, fm_id = reg_row$fm_id, locus_id = reg_row$locus_id, ld_id = reg_row$ld_id,
                     n_region = n_region, n_used = length(z), snp_coverage = length(z) / n_region,
                     n_af_drop = n_af_drop, n_kriging_removed = length(kriged), s = s_est,
                     top_snp = top_snp, top_snp_used = top_snp %in% names(z),
                     converged = isTRUE(fit$converged), niter = fit$niter, n_cs = nrow(csd),
                     L_saturated = nrow(csd) >= as.integer(cfg$L %||% 10L),
                     max_pip = max(pip), error = NA_character_)
  saveRDS(summ, file.path(rdir, "summary.rds"))
  writeLines(kriged, file.path(rdir, "kriging_removed.txt"))
  file.create(done)
  summ
}

run_finemap <- function(config, sex) {
  stage <- "finemap"
  log_info(stage, sprintf("=== Fine-mapping (SuSiE-RSS, UKB in-sample LD): %s ===", sex))
  if (!requireNamespace("susieR", quietly = TRUE)) log_fatal(stage, "susieR not installed (scripts/install_postgwas_tools.sh)")
  cfg <- config$finemap %||% list()
  out_dir <- stage_dir(config, sex, stage)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  reg <- read_factor_registry(config, sex)
  regs <- fread(file.path(stage_dir(config, sex, "ld_panel"), "finemap_regions.tsv"))
  loci <- fread(file.path(stage_dir(config, sex, "loci"), "regions.tsv"))
  regs <- merge(regs, loci[, .(locus_id = region_id, lead_snp, heterogeneous)], by = "locus_id")
  workers <- as.integer(cfg$workers %||% 12L)
  log_info(stage, sprintf("%d regions, %d workers", nrow(regs), workers))
  # each factor's full tier is read once here; the forked workers share it copy-on-write
  g_all <- lapply(setNames(nm = unique(regs$factor)), function(f) load_postgwas(config, sex, f, "full", cols = FINEMAP_COLS))
  clear_postgwas_cache()
  res <- parallel::mclapply(seq_len(nrow(regs)), function(i) tryCatch(
    run_finemap_region(config, sex, regs[i], cfg, reg[[regs$factor[i]]], out_dir, g_all[[regs$factor[i]]]),
    error = function(e) data.table(factor = regs$factor[i], fm_id = regs$fm_id[i], error = conditionMessage(e))),
    mc.cores = workers, mc.preschedule = FALSE)
  diag <- rbindlist(res, fill = TRUE)
  for (col in c("s", "snp_coverage")) if (!col %in% names(diag)) set(diag, j = col, value = NA_real_)
  diag <- merge(diag, regs[, .(fm_id, heterogeneous)], by = "fm_id", all.x = TRUE)
  diag[, ld_reliability := fifelse(is.na(s), NA_character_, fifelse(s > 0.3, "unreliable", fifelse(s > 0.1, "caution", "ok")))]
  fwrite(diag, file.path(out_dir, "region_diagnostics.tsv"), sep = "\t")
  n_err <- sum(!is.na(diag$error))
  if (n_err) log_warn(stage, sprintf("%d region(s) failed; see region_diagnostics.tsv", n_err))

  cs <- rbindlist(lapply(seq_len(nrow(regs)), function(i) {
    p <- file.path(out_dir, regs$factor[i], regs$ld_id[i], "cs.tsv")
    if (file.exists(p) && file.size(p) > 0) fread(p) else NULL
  }), fill = TRUE)
  if (nrow(cs)) {
    cs <- merge(cs, regs[, .(fm_id, heterogeneous)], by = "fm_id", all.x = TRUE)
    cs <- merge(cs, diag[, .(fm_id, ld_reliability, snp_coverage)], by = "fm_id", all.x = TRUE)
    cs[, factor_level := !heterogeneous & !qflag]
  }
  fwrite(cs, file.path(out_dir, "all_cs.tsv"), sep = "\t")
  if (nrow(cs)) fwrite(cs[, .(factor, cs_id, locus_id, lead_snp, lead_chr, lead_bp, lead_pip, lead_z, lead_log10p,
                              size, purity_min, heterogeneous, qflag, factor_level, ld_reliability)],
                       file.path(out_dir, "pip_leads.tsv"), sep = "\t")
  log_info(stage, sprintf("%d credible sets in %d regions (%d factor-level); median size %s",
                          nrow(cs), uniqueN(cs$fm_id), sum(cs$factor_level %in% TRUE),
                          if (nrow(cs)) median(cs$size) else NA))
  clear_postgwas_cache()
  write_stage_manifest(stage, sex, config, file.path(out_dir, c("all_cs.tsv", "pip_leads.tsv", "region_diagnostics.tsv")))
  invisible(cs)
}
