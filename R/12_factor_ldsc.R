# Factor LDSC stage (GenomicSEM paper checks on the factor GWAS):
# - h2 / intercept / ratio per factor on the qfilter tier (primary) and the full tier;
#   an intercept near 1 means the factor's lambda_GC is polygenic signal, not inflation
# - LDSC of the Q_SNP statistics (intercept = calibration, slope = genetic signal in
#   heterogeneity), as in Grotzinger et al. 2019
# - rg between factors (the bifactor makes them orthogonal: rg should be ~0) and between
#   each factor and its indicators (should track the standardised loadings)
# N for every factor file is N_hat_std (std.lv scale); intercept, ratio and rg do not
# depend on N, and h2 on the marker scale is h2_std * N_hat_std / N_hat_marker.

# Munged (LDSC) format restricted to HapMap3, alleles aligned to w_hm3.snplist.
write_munged_sumstats <- function(dt, hm3, n, path) {
  m <- merge(dt[, .(SNP, A1, A2, Z)], hm3, by = "SNP", suffixes = c("", ".h"))
  same <- m$A1 == m$A1.h & m$A2 == m$A2.h
  flip <- m$A1 == m$A2.h & m$A2 == m$A1.h
  m <- m[same | flip]
  flip <- flip[same | flip]
  m[flip, Z := -Z]
  out <- m[, .(SNP, N = n, Z, A1 = A1.h, A2 = A2.h)][is.finite(Z)]
  fwrite(out, path, sep = "\t")
  nrow(out)
}

# GenomicSEM::ldsc in out_dir (it writes its log to the working directory).
run_gsem_ldsc <- function(config, files, names, out_dir, log_name, stand = FALSE) {
  saved <- getwd()
  on.exit(setwd(saved), add = TRUE)
  setwd(out_dir)
  call <- function(st) GenomicSEM::ldsc(traits = normalizePath(files), sample.prev = rep(NA, length(files)),
                                        population.prev = rep(NA, length(files)),
                                        ld = config$paths$ld_scores, wld = config$paths$ld_scores,
                                        trait.names = names, stand = st, ldsc.log = log_name)
  res <- tryCatch(call(stand), error = function(e) {
    if (!stand) stop(e)
    log_warn("factor_ldsc", sprintf("ldsc(stand=TRUE) failed (%s); rerunning with stand=FALSE", conditionMessage(e)))
    call(FALSE)
  })
  list(res = res, log = file.path(out_dir, paste0(log_name, "_ldsc.log")))
}

# Per-trait lines of a GenomicSEM ldsc log -> one row per trait, in log order.
parse_gsem_ldsc_log <- function(log_file) {
  x <- readLines(log_file, warn = FALSE)
  starts <- grep("^Heritability Results for trait:", x)
  num <- function(line, pat) {
    m <- regmatches(line, regexec(pat, line))[[1]]
    if (length(m) < 2) NA_real_ else as.numeric(m[2])
  }
  rbindlist(lapply(starts, function(s) {
    blk <- x[s:min(length(x), s + 12)]
    get <- function(key) blk[grep(key, blk)[1]]
    data.table(trait_file = basename(sub("^Heritability Results for trait: ", "", x[s])),
               mean_chisq = num(get("^Mean Chi"), ":\\s*([-0-9.eE]+)"),
               lambda_gc = num(get("^Lambda GC"), ":\\s*([-0-9.eE]+)"),
               intercept = num(get("^Intercept"), ":\\s*([-0-9.eE]+)"),
               intercept_se = num(get("^Intercept"), "\\(([-0-9.eE]+)\\)"),
               ratio = num(get("^Ratio"), ":\\s*([-0-9.eE]+)"),
               ratio_se = num(get("^Ratio"), "\\(([-0-9.eE]+)\\)"),
               h2 = num(get("Scale h2"), ":\\s*([-0-9.eE]+)"),
               h2_se = num(get("Scale h2"), "\\(([-0-9.eE]+)\\)"),
               h2_z = num(get("^h2 Z"), ":\\s*([-0-9.eE]+)"))
  }))
}

# rg and SE between `rows` and every trait, from a stand=TRUE GenomicSEM ldsc result.
rg_table <- function(res, rows) {
  if (is.null(res$S_Stand) || is.null(res$V_Stand)) return(data.table())
  tr <- colnames(res$S)
  k <- length(tr)
  low <- which(lower.tri(matrix(0, k, k), diag = TRUE), arr.ind = TRUE)
  se <- sqrt(diag(res$V_Stand))
  rbindlist(lapply(rows, function(a) {
    i <- match(a, tr)
    rbindlist(lapply(setdiff(seq_len(k), i), function(j) {
      v <- which((low[, 1] == max(i, j)) & (low[, 2] == min(i, j)))
      data.table(trait_a = a, trait_b = tr[j], rg = res$S_Stand[i, j], rg_se = se[v])
    }))
  }))
}

run_factor_ldsc <- function(config, sex) {
  stage <- "factor_ldsc"
  log_info(stage, sprintf("=== Factor LDSC: %s ===", sex))
  cfg <- config$factor_ldsc %||% list()
  out_dir <- stage_dir(config, sex, stage)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  reg <- read_factor_registry(config, sex)
  factors <- postgwas_factors(config)
  hm3 <- fread(config$paths$hm3_snplist, select = c("SNP", "A1", "A2"))

  files <- list()
  for (tier in c("qfilter", "full")) {
    for (f in factors) {
      dt <- load_postgwas(config, sex, f, tier, cols = c("SNP", "A1", "A2", "Z", "Q_SNP_pval"))
      p <- file.path(out_dir, sprintf("%s_%s.sumstats.gz", f, tier))
      n_w <- write_munged_sumstats(dt, hm3, reg[[f]]$N_hat_std, p)
      files[[paste(f, tier)]] <- data.table(factor = f, tier = tier, file = p, n_hm3 = n_w)
      if (tier == "full" && isTRUE(cfg$qsnp_ldsc %||% TRUE)) {
        q <- dt[!is.na(Q_SNP_pval) & Q_SNP_pval > 0, .(SNP, A1, A2, Z = qnorm(Q_SNP_pval / 2, lower.tail = FALSE))]
        pq <- file.path(out_dir, sprintf("%s_qsnp.sumstats.gz", f))
        files[[paste(f, "qsnp")]] <- data.table(factor = f, tier = "qsnp", file = pq,
                                                n_hm3 = write_munged_sumstats(q, hm3, reg[[f]]$N_hat_std, pq))
      }
    }
  }
  files <- rbindlist(files)
  clear_postgwas_cache()

  res_rows <- list()
  for (tr in c("qfilter", "full", "qsnp")) {
    tier <- tr
    ft <- files[files$tier == tr]
    if (!nrow(ft)) next
    r <- run_gsem_ldsc(config, ft$file, paste0(ft$factor, "_", tier), out_dir, paste0("factors_", tier),
                       stand = tier == "qfilter")
    pr <- parse_gsem_ldsc_log(r$log)
    pr[, `:=`(factor = ft$factor[match(trait_file, basename(ft$file))], tier = tier)]
    res_rows[[tier]] <- pr
    if (tier == "qfilter") {
      fwrite(rg_table(r$res, colnames(r$res$S)), file.path(out_dir, "rg_factors.csv"))
    }
  }
  s <- rbindlist(res_rows)
  s[, N_hat_std := vapply(factor, function(f) reg[[f]]$N_hat_std, numeric(1))]
  s[, h2_marker_scale := h2 * N_hat_std / vapply(factor, function(f) reg[[f]]$N_hat_marker, numeric(1))]
  lo <- as.numeric(cfg$intercept_warn_low %||% 0.9); hi <- as.numeric(cfg$intercept_warn_high %||% 1.1)
  for (i in which(s$tier != "qsnp" & (s$intercept < lo | s$intercept > hi))) {
    log_warn(stage, sprintf("%s (%s): LDSC intercept %.3f outside [%.2f, %.2f]", s$factor[i], s$tier[i],
                            s$intercept[i], lo, hi))
  }
  for (i in seq_len(nrow(s))) {
    log_info(stage, sprintf("%-16s %-7s intercept %.3f (%.3f) ratio %s h2(std) %.3f (%.3f) mean chi2 %.3f",
                            s$factor[i], s$tier[i], s$intercept[i], s$intercept_se[i],
                            format(round(s$ratio[i], 3)), s$h2[i], s$h2_se[i], s$mean_chisq[i]))
  }
  fwrite(s, file.path(out_dir, "summary.csv"))

  if (isTRUE(cfg$rg_with_indicators %||% TRUE)) {
    inds <- unique(unlist(lapply(factors, function(f) reg[[f]]$indicators)))
    ind_files <- file.path(config$paths$output_dir, sex, "munge", sprintf("%s_%s.sumstats.gz", sex, inds))
    have <- file.exists(ind_files)
    if (!all(have)) log_warn(stage, sprintf("munged indicator files missing: %s", paste(inds[!have], collapse = ", ")))
    fq <- files[tier == "qfilter"]
    log_info(stage, sprintf("rg: %d factors x %d indicators (multivariable LDSC; slow)", nrow(fq), sum(have)))
    r <- run_gsem_ldsc(config, c(fq$file, ind_files[have]), c(fq$factor, inds[have]), out_dir,
                       "factors_vs_indicators", stand = TRUE)
    rg <- rg_table(r$res, fq$factor)
    fwrite(rg, file.path(out_dir, "rg_indicators.csv"))
    fit <- readRDS(file.path(config$paths$output_dir, sex, "cfa", "model_bifactor_full.rds"))
    ld <- as.data.table(fit$results)[op == "=~", .(trait_a = lhs, trait_b = rhs, std_loading = as.numeric(STD_All))]
    cmp <- merge(rg, ld, by = c("trait_a", "trait_b"))
    fwrite(cmp, file.path(out_dir, "rg_vs_loadings.csv"))
    for (f in unique(cmp$trait_a)) {
      d <- cmp[trait_a == f]
      log_info(stage, sprintf("%s: cor(rg with indicator, standardised loading) = %.3f over %d indicators",
                              f, cor(d$rg, d$std_loading), nrow(d)))
    }
  }
  write_stage_manifest(stage, sex, config, list.files(out_dir, full.names = TRUE))
  invisible(s)
}
