# FinnGen summary-level PRS evaluation (European check). For a score S = sum_j w_j G_j and a
# FinnGen endpoint Y with per-allele log-odds beta_j (alt allele) and alt frequency p_j:
#   cov(S, Y) ~ K(1-K) * sum_j w_j beta_j 2 p_j (1-p_j)        (linear-probability approximation)
#   var(S)    = sample variance of S in an LD reference (plink --score on 1000G EUR, or a UKB
#               EUR subset when plink2 is available)
#   R2_obs    = cov^2 / (var(S) K(1-K)),  then Lee et al. 2012 to the liability scale with K = P
#               (FinnGen is a cohort, so the case fraction is the prevalence)
# CIs by leave-one-chromosome-out jackknife of the numerator. Per endpoint and per weight set;
# the mean liability R2 over a factor's mapped endpoints is the summary number. This is the
# lassosum / PUMAS pseudo-validation idea; it needs no individual-level data.

FINNGEN_URL <- "https://storage.googleapis.com/finngen-public-data-r12/summary_stats/release/finngen_R12_%s.gz"

read_finngen_endpoint_map <- function(path) {
  m <- fread(path, sep = "\t", colClasses = "character")
  m[!is.na(endpoint) & nzchar(endpoint)]
}

finngen_file <- function(config, ep) {
  dir <- config$refs$finngen_dir; dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  local <- file.path(config$refs$finngen_local_dir %||% "", paste0(ep, "_fg.tsv.gz"))
  if (nzchar(config$refs$finngen_local_dir %||% "") && file.exists(local)) return(local)
  f <- file.path(dir, sprintf("finngen_R12_%s.gz", ep))
  if (!file.exists(f)) {
    run_tool("curl", c("-sS", "-L", "--fail", "--retry", "5", "-C", "-", "-o", paste0(f, ".part"), sprintf(FINNGEN_URL, ep)),
             "prs_eval", paste0(f, ".curl.log"), expect = paste0(f, ".part"))
    file.rename(paste0(f, ".part"), f)
  }
  f
}

# One endpoint file -> (rsids, ref, alt, beta, se, af), multi-rsID cells split, restricted to
# `keep_snps` when given. Read once per endpoint; aligned to each weight set by finngen_align().
finngen_read <- function(f, keep_snps = NULL) {
  x <- fread(f, select = c("#chrom", "rsids", "ref", "alt", "beta", "sebeta", "af_alt"), colClasses = list(character = c("#chrom", "rsids", "ref", "alt")))
  setnames(x, c("chr", "rsids", "ref", "alt", "beta", "se", "af"))
  x <- x[rsids != "" & !is.na(beta)]
  multi <- grepl(",", x$rsids, fixed = TRUE)
  if (any(multi)) {
    sp <- strsplit(x$rsids[multi], ",", fixed = TRUE)
    x <- rbind(x[!multi], x[rep(which(multi), lengths(sp))][, rsids := unlist(sp)])
  }
  if (!is.null(keep_snps)) x <- x[rsids %in% keep_snps]
  x
}

# Endpoint statistics (finngen_read output) restricted to the weight SNPs, effect oriented to A1 of `w`.
finngen_align <- function(x, w) {
  if (is.character(x)) x <- finngen_read(x, w$SNP)
  m <- merge(w, x[rsids %in% w$SNP], by.x = "SNP", by.y = "rsids")
  m <- m[(A1 == alt & A2 == ref) | (A1 == ref & A2 == alt)]
  m[, beta_a1 := fifelse(A1 == alt, beta, -beta)]
  m[, p_a1 := fifelse(A1 == alt, af, 1 - af)]
  m[, CHR := suppressWarnings(as.integer(sub("^chr", "", chr)))]   # chromosome from FinnGen (weight tables carry none)
  m[!duplicated(SNP) & !is.na(CHR), .(SNP, CHR, weight, beta_a1, se, p_a1)]
}

# var(S) in the LD reference: plink 1.9 --score (sum of per-allele weights) on 1000G EUR.
finngen_score_variance <- function(config, w, prs_id, work, weights_mtime = NULL) {
  pre <- file.path(work, paste0("var_", prs_id))
  prof <- paste0(pre, ".profile")
  stale <- file.exists(prof) && is.finite(weights_mtime %||% NA_real_) && as.numeric(file.mtime(prof)) < weights_mtime
  if (!file.exists(prof) || stale) {
    wf <- paste0(pre, ".weights.txt")
    fwrite(w[, .(SNP, A1, weight)], wf, sep = "\t", col.names = FALSE)
    run_tool(tool_path(config, "plink"), c("--bfile", config$refs$ld_1000g_bfile, "--score", wf, "1", "2", "3", "sum", "--out", pre),
             "prs_eval", paste0(pre, ".stdout.log"), expect = paste0(pre, ".profile"))
  }
  pr <- fread(paste0(pre, ".profile"))
  list(var = stats::var(pr$SCORESUM), n_ref = nrow(pr), n_snps_scored = pr$CNT[1] / 2)
}

finngen_r2 <- function(m, var_s, K) {
  num_chr <- m[, .(c = sum(weight * beta_a1 * 2 * p_a1 * (1 - p_a1))), by = CHR]
  cov_all <- sum(num_chr$c)
  r2 <- function(cv) K * (1 - K) * cv^2 / var_s
  r2_obs <- r2(cov_all)
  jack <- vapply(seq_len(nrow(num_chr)), function(i) r2(cov_all - num_chr$c[i]), numeric(1))
  nj <- length(jack); se <- sqrt((nj - 1) / nj * sum((jack - mean(jack))^2))
  t <- stats::qnorm(1 - K); z <- stats::dnorm(t)
  liab <- function(r) r * K * (1 - K) / z^2
  data.table(cov = cov_all, direction = sign(cov_all), r2_obs = r2_obs, r2_liability = liab(r2_obs),
             ci_low = liab(max(0, r2_obs - 1.96 * se)), ci_high = liab(r2_obs + 1.96 * se), jackknife_se_obs = se, n_snps = nrow(m))
}

run_prs_eval_finngen <- function(config, sex, sets, out_dir) {
  stage <- "prs_eval"
  cfg <- config$prs_eval$finngen %||% list()
  fd <- file.path(out_dir, "finngen"); work <- file.path(fd, "work"); dir.create(work, recursive = TRUE, showWarnings = FALSE)
  map_path <- cfg$endpoints %||% file.path(config$project$root, "manifest", "indicator_finngen_endpoints.tsv")
  if (!is_absolute_path(map_path)) map_path <- file.path(config$project$root, map_path)
  map <- read_finngen_endpoint_map(map_path)
  man <- fread(config$refs$finngen_manifest, sep = "\t")
  man[, K := num_cases / (num_cases + num_controls)]
  reg <- read_factor_registry(config, sex)
  eps <- unique(map$endpoint)
  # Results are cached per weight set (results.tsv + done_sets.tsv): a rerun recomputes only sets
  # that are new or whose weight file is newer than their cached result, because every endpoint
  # file (~750 MB) has to be re-read for any set that is missing.
  res_path <- file.path(fd, "results.tsv"); done_path <- file.path(fd, "done_sets.tsv")
  cur <- data.table(prs_id = vapply(sets, `[[`, "", "prs_id"),
                    weights_mtime = vapply(sets, function(s) as.numeric(s$mtime %||% Inf), numeric(1)))
  have_cache <- file.exists(res_path) && file.exists(done_path)
  done <- if (have_cache) fread(done_path, colClasses = list(character = "prs_id")) else data.table(prs_id = character(), weights_mtime = numeric(), done_at = numeric())
  is_done <- cur$prs_id %in% done$prs_id & cur$weights_mtime <= done$done_at[match(cur$prs_id, done$prs_id)]
  cached <- if (have_cache && any(is_done)) fread(res_path, colClasses = list(character = c("prs_id", "threshold")))[prs_id %in% cur$prs_id[is_done]] else data.table()
  todo <- sets[!is_done]
  log_info(stage, sprintf("FinnGen: %d endpoints for %d indicators, %d weight sets (%d cached, %d to compute)",
                          length(eps), nrow(map), length(sets), sum(is_done), length(todo)))
  rows <- list()
  vars <- list()
  for (s in todo) {
    v <- tryCatch(finngen_score_variance(config, s$dt, s$prs_id, work, s$mtime), error = function(e) NULL)
    if (is.null(v)) { log_warn(stage, sprintf("FinnGen: variance of %s failed", s$prs_id)); next }
    vars[[s$prs_id]] <- v
  }
  all_snps <- unique(unlist(lapply(todo, function(s) s$dt$SNP)))
  for (ep in if (length(todo)) eps else character(0)) {
    f <- tryCatch(finngen_file(config, ep), error = function(e) { log_warn(stage, sprintf("FinnGen %s: %s", ep, conditionMessage(e))); NULL })
    if (is.null(f)) next
    K <- man[phenocode == ep, K][1]
    if (!is.finite(K)) { log_warn(stage, sprintf("FinnGen %s: not in the manifest", ep)); next }
    codes <- map[endpoint == ep, icd10]
    x <- finngen_read(f, all_snps)
    log_info(stage, sprintf("FinnGen %s: %s weight SNPs matched by rsID", ep, format(nrow(x), big.mark = ",")))
    for (s in todo) {
      if (is.null(vars[[s$prs_id]])) next
      fac_codes <- reg[[s$factor]]$indicators
      if (!any(codes %in% fac_codes) && s$kind == "factor") next
      m <- finngen_align(x, s$dt)
      if (nrow(m) < 100) next
      r <- finngen_r2(m, vars[[s$prs_id]]$var, K)
      rows[[length(rows) + 1]] <- cbind(data.table(endpoint = ep, icd10 = paste(codes, collapse = ";"), prs_id = s$prs_id, kind = s$kind,
                                                   name = s$name, factor = s$factor, method = s$method, threshold = s$threshold,
                                                   n_cases = man[phenocode == ep, num_cases][1], K = K, n_ref = vars[[s$prs_id]]$n_ref), r)
    }
  }
  res <- rbindlist(c(list(cached), rows), fill = TRUE)
  fwrite(res, res_path, sep = "\t")
  fwrite(rbind(done[prs_id %in% cur$prs_id[is_done]],
               data.table(prs_id = cur$prs_id[!is_done], weights_mtime = cur$weights_mtime[!is_done], done_at = as.numeric(Sys.time()))),
         done_path, sep = "\t")
  if (nrow(res)) {
    summ <- res[, .(n_endpoints = .N, mean_r2_liability = mean(r2_liability), median_r2_liability = stats::median(r2_liability),
                    n_positive_direction = sum(direction > 0)), by = .(prs_id, kind, name, factor, method, threshold)][order(factor, -mean_r2_liability)]
    fwrite(summ, file.path(fd, "summary.tsv"), sep = "\t")
    for (f in unique(summ$factor)) {
      b <- summ[factor == f & kind == "factor"][1]
      if (nrow(b)) log_info(stage, sprintf("FinnGen %s: best set %s, mean liability R2 %.4f over %d endpoints", f, b$prs_id, b$mean_r2_liability, b$n_endpoints))
    }
  }
  invisible(res)
}
