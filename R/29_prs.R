# PRS stage (GWASTutorial "Polygenic Risk Score"): per-factor PRS weights from the qfilter tier
# (autosomal), by a registry of methods. Each method writes weights/<method>/<factor>.weights.tsv
# with columns SNP A1 A2 weight (+ method extras). Scoring and evaluation in a target cohort run
# only when prs$target is configured (not yet decided).
#   ct            clumping + thresholding (plink 1.9, 1000G EUR, r2 0.1 / 250 kb as in the guide),
#                 one weight set per p-value threshold
#   prscs         PRS-CS-auto (UKB EUR LD blocks), n_gwas = N_hat_std
#   ldpred2_auto  LDpred2-auto (bigsnpr; UKB HapMap3 LD of Prive et al. 2022), no tuning sample
#   sbayesrc      SBayesRC (UKB imputed LD eigen-decomposition + baselineLD annotations), via its
#                 container

prs_input <- function(config, sex, f) {
  d <- load_postgwas(config, sex, f, "qfilter",
                     cols = c("SNP", "CHR", "BP", "A1", "A2", "EAF_A1", "beta", "SE", "P", "LOG10P"))
  d <- d[CHR <= 22L & !is.na(EAF_A1)]
  d[, P := pmax(P, 10^-pmin(LOG10P, 300))]
  d[]
}

prs_ct <- function(config, d, f, dir, cfg, n) {
  thr <- as.numeric(unlist(cfg$ct_thresholds %||% c(5e-8, 1e-6, 1e-4, 1e-3, 0.01, 0.05, 0.1, 0.2, 0.5, 1)))
  cl <- plink_clump(config, config$refs$ld_1000g_bfile, d[SNP %in% fread(paste0(config$refs$ld_1000g_bfile, ".bim"),
                                                                       header = FALSE, select = 2L)$V2],
                    p1 = 1, p2 = 1, r2 = as.numeric(cfg$ct_r2 %||% 0.1), kb = as.integer(cfg$ct_kb %||% 250L),
                    out_prefix = file.path(dir, paste0(f, "_clump")), stage = "prs")
  w <- d[SNP %in% cl$SNP, .(SNP, A1, A2, weight = beta, P)]
  fwrite(data.table(name = paste0("pT_", format(thr, scientific = TRUE)), lower = 0, upper = thr),
         file.path(dir, paste0(f, ".ranges")), sep = " ", col.names = FALSE)
  w[, n_thresholds_passed := vapply(P, function(p) sum(p <= thr), integer(1))]
  w
}

prs_prscs <- function(config, d, f, dir, cfg, n) {
  py <- tool_path(config, "py")
  prscs <- tool_path(config, "prscs")
  sst <- file.path(dir, paste0(f, ".prscs.sst"))
  fwrite(d[, .(SNP, A1, A2, BETA = beta, SE)], sst, sep = "\t")
  outp <- file.path(dir, paste0(f, "_prscs"))
  env <- c("MKL_NUM_THREADS=1", "NUMEXPR_NUM_THREADS=1", "OMP_NUM_THREADS=1")
  r <- parallel::mclapply(1:22, function(c) tryCatch(run_tool(py, c(prscs, paste0("--ref_dir=", config$refs$prscs_ref),
      paste0("--bim_prefix=", config$refs$ld_1000g_bfile), paste0("--sst_file=", sst), paste0("--n_gwas=", round(n)),
      paste0("--out_dir=", outp), paste0("--chrom=", c), paste0("--seed=", cfg$seed %||% 20261108)),
      "prs", paste0(outp, "_chr", c, ".log"), env = env), error = function(e) e),
    mc.cores = as.integer(cfg$workers %||% 22L))
  bad <- vapply(r, inherits, logical(1), "error")
  if (any(bad)) stop("PRS-CS failed on ", sum(bad), " chromosome(s): ", conditionMessage(r[bad][[1]]))
  files <- list.files(dirname(outp), pattern = paste0("^", basename(outp), "_pst_eff_.*_chr[0-9]+\\.txt$"), full.names = TRUE)
  w <- rbindlist(lapply(files, fread, header = FALSE, col.names = c("CHR", "SNP", "BP", "A1", "A2", "weight")))
  w[, .(SNP, A1, A2, weight)]
}

prs_ldpred2_auto <- function(config, d, f, dir, cfg, n) {
  if (!requireNamespace("bigsnpr", quietly = TRUE)) stop("bigsnpr not installed")
  ld_dir <- config$refs$ldpred2_ld_dir
  map <- as.data.table(readRDS(file.path(ld_dir, "map.rds")))
  ss <- d[, .(chr = CHR, pos = BP, a0 = A2, a1 = A1, beta, beta_se = SE, n_eff = n, rsid = SNP)]
  info <- as.data.table(bigsnpr::snp_match(as.data.frame(ss), as.data.frame(map[, .(chr, pos, a0, a1)]),
                                           return_flip_and_rev = TRUE))
  tmp <- file.path(dir, paste0(f, "_ldpred2_corr"))
  corr <- NULL
  for (c in 1:22) {
    ind <- info[chr == c, `_NUM_ID_`]
    blk <- readRDS(file.path(ld_dir, sprintf("LD_with_blocks_chr%d.rds", c)))
    first <- which(map$chr == c)[1] - 1L
    sub <- blk[ind - first, ind - first]
    if (is.null(corr)) corr <- bigsparser::as_SFBM(sub, tmp, compact = TRUE) else corr$add_columns(sub, nrow(corr))
  }
  ld <- map$ld[info$`_NUM_ID_`]
  chi2 <- (info$beta / info$beta_se)^2
  h2_ldsc <- bigsnpr::snp_ldsc(ld, ld_size = nrow(map), chi2 = chi2, sample_size = info$n_eff, blocks = NULL)[["h2"]]
  set.seed(as.integer(cfg$seed %||% 20261108L))
  auto <- bigsnpr::snp_ldpred2_auto(corr, as.data.frame(info[, .(beta, beta_se, n_eff)]), h2_init = h2_ldsc,
                                    vec_p_init = bigsnpr::seq_log(1e-4, 0.2, length.out = as.integer(cfg$ldpred2_chains %||% 30L)),
                                    allow_jump_sign = FALSE, shrink_corr = 0.95,
                                    ncores = as.integer(cfg$workers %||% 22L))
  rng <- vapply(auto, function(a) diff(range(a$corr_est)), numeric(1))
  keep <- which(rng > 0.95 * stats::quantile(rng, 0.95, na.rm = TRUE))
  w <- rowMeans(sapply(auto[keep], function(a) a$beta_est))
  data.table(SNP = info$rsid, A1 = info$a1, A2 = info$a0, weight = w,
             h2_est = mean(vapply(auto[keep], `[[`, numeric(1), "h2_est")),
             p_est = mean(vapply(auto[keep], `[[`, numeric(1), "p_est")), n_chains_kept = length(keep))
}

prs_sbayesrc <- function(config, d, f, dir, cfg, n) {
  sif <- config$tools$sbayesrc_sif %||% ""
  if (!file.exists(sif)) stop("SBayesRC container missing (tools$sbayesrc_sif)")
  ma <- file.path(dir, paste0(f, ".ma"))
  fwrite(d[, .(SNP, A1, A2, freq = EAF_A1, b = beta, se = SE, p = P, N = round(n))], ma, sep = "\t")
  pre <- file.path(dir, paste0(f, "_sbrc"))
  thr <- as.integer(cfg$sbayesrc_threads %||% 16L)
  rcode <- sprintf(paste0("SBayesRC::tidy(mafile='%s', LDdir='%s', output='%s_tidy.ma', log2file=TRUE);",
                          "SBayesRC::impute(mafile='%s_tidy.ma', LDdir='%s', output='%s_imp.ma', log2file=TRUE);",
                          "SBayesRC::sbayesrc(mafile='%s_imp.ma', LDdir='%s', outPrefix='%s', annot='%s', log2file=TRUE)"),
                   ma, config$refs$sbayesrc_ld, pre, pre, config$refs$sbayesrc_ld, pre, pre, config$refs$sbayesrc_ld,
                   pre, config$refs$sbayesrc_annot)
  run_tool("apptainer", c("exec", "-B", "/mnt/sdg", sif, "Rscript", "--vanilla", "-e", shQuote(rcode)), "prs",
           paste0(pre, ".stdout.log"), expect = paste0(pre, ".txt"), env = paste0("OMP_NUM_THREADS=", thr))
  # tidy() aligns alleles to the LD reference and impute() adds SNPs absent from the GWAS, so A2
  # comes from the reference's snp.info (A1 is the effect allele of BETA)
  w <- fread(paste0(pre, ".txt"))
  info <- fread(file.path(config$refs$sbayesrc_ld, "snp.info"), select = c("ID", "A1", "A2"))
  i <- match(w$SNP, info$ID)
  w[, .(SNP, A1, A2 = fifelse(A1 == info$A1[i], info$A2[i], info$A1[i]), weight = BETA)]
}

PRS_METHODS <- list(ct = prs_ct, prscs = prs_prscs, ldpred2_auto = prs_ldpred2_auto, sbayesrc = prs_sbayesrc)

run_prs <- function(config, sex) {
  stage <- "prs"
  log_info(stage, sprintf("=== PRS weights: %s ===", sex))
  cfg <- config$prs %||% list()
  out_dir <- stage_dir(config, sex, stage); dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  reg <- read_factor_registry(config, sex)
  methods <- unlist(cfg$methods %||% names(PRS_METHODS))
  unknown <- setdiff(methods, names(PRS_METHODS))
  if (length(unknown)) log_fatal(stage, sprintf("unknown PRS method(s): %s", paste(unknown, collapse = ", ")))
  summ <- list()
  for (f in postgwas_factors(config)) {
    d <- prs_input(config, sex, f)
    for (m in methods) {
      dir <- file.path(out_dir, "weights", m); dir.create(dir, recursive = TRUE, showWarnings = FALSE)
      out <- file.path(dir, paste0(f, ".weights.tsv"))
      if (file.exists(out)) { log_info(stage, sprintf("%s/%s: weights exist; skipping", f, m)); next }
      t0 <- Sys.time()
      w <- tryCatch(PRS_METHODS[[m]](config, d, f, dir, cfg, reg[[f]]$N_hat_std), error = function(e) e)
      if (inherits(w, "error")) {
        log_warn(stage, sprintf("%s/%s failed: %s", f, m, conditionMessage(w)))
        summ[[paste(f, m)]] <- data.table(factor = f, method = m, n_snps = NA_integer_, error = conditionMessage(w))
        next
      }
      fwrite(w, out, sep = "\t")
      summ[[paste(f, m)]] <- data.table(factor = f, kind = "factor", method = m, n_snps = nrow(w), error = NA_character_,
                                        minutes = round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1))
      log_info(stage, sprintf("%s/%s: %s SNPs", f, m, format(nrow(w), big.mark = ",")))
    }
  }
  # comparator scores from external summary statistics (extdata/comparator_sumstats.tsv)
  ext <- read_external_manifest(config)
  if (!is.null(ext) && nrow(ext)) {
    ref <- external_reference(config, sex)
    ext_methods <- unlist(cfg$external_methods %||% methods)
    ext_summ <- list()
    for (i in seq_len(nrow(ext))) {
      spec <- as.list(ext[i])
      todo <- ext_methods[!file.exists(file.path(out_dir, "weights", ext_methods, paste0(spec$id, ".weights.tsv")))]
      if (!length(todo)) { log_info(stage, sprintf("%s: all weights exist; skipping", spec$id)); next }
      x <- tryCatch(read_external_sumstats(spec, ref), error = function(e) e)
      if (inherits(x, "error")) {
        log_warn(stage, sprintf("%s: could not read summary statistics: %s", spec$id, conditionMessage(x)))
        summ[[spec$id]] <- data.table(factor = spec$id, kind = "external", method = NA_character_, n_snps = NA_integer_,
                                      error = conditionMessage(x))
        next
      }
      ext_summ[[spec$id]] <- x$summary
      for (m in todo) {
        dir <- file.path(out_dir, "weights", m); dir.create(dir, recursive = TRUE, showWarnings = FALSE)
        t0 <- Sys.time()
        w <- tryCatch(PRS_METHODS[[m]](config, x$d, spec$id, dir, cfg, x$n), error = function(e) e)
        if (inherits(w, "error")) {
          log_warn(stage, sprintf("%s/%s failed: %s", spec$id, m, conditionMessage(w)))
          summ[[paste(spec$id, m)]] <- data.table(factor = spec$id, kind = "external", method = m, n_snps = NA_integer_,
                                                  error = conditionMessage(w))
          next
        }
        fwrite(w, file.path(dir, paste0(spec$id, ".weights.tsv")), sep = "\t")
        summ[[paste(spec$id, m)]] <- data.table(factor = spec$id, kind = "external", method = m, n_snps = nrow(w),
                                                error = NA_character_,
                                                minutes = round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1))
        log_info(stage, sprintf("%s/%s: %s SNPs", spec$id, m, format(nrow(w), big.mark = ",")))
      }
    }
    if (length(ext_summ)) fwrite(rbindlist(ext_summ), file.path(out_dir, "external_sumstats_summary.csv"))
    rm(ref)
  }
  if (is.null(cfg$target)) log_info(stage, "prs$target not set: weights only (scoring/evaluation in prs_eval)")
  fwrite(rbindlist(summ, fill = TRUE), file.path(out_dir, "prs_summary.csv"))
  clear_postgwas_cache()
  write_stage_manifest(stage, sex, config, file.path(out_dir, "prs_summary.csv"))
  invisible(rbindlist(summ, fill = TRUE))
}
