# PRS evaluation stage: scores every weight set (factor x method, plus comparator scores built
# from external summary statistics) in the UK Biobank South Asian and Black participants
# (/mnt/sdh/upamanyu/UKB/Genomics_PGEN/Asian_Black, read-only; zero overlap with the Pan-UKB
# EUR discovery) and fits, per group separately:
#  - primary: negative binomial of the ARD indicator-diagnosis count (ICD-10 three-character
#    unit) on the score -> IRR per SD; domain scores on their own indicator counts
#  - sensitivity: the same on the GBD-cause roll-up count; all-cause mortality (Cox, age scale)
#  - robustness (ARD scores): GCTA unrelated lists instead of KING; global UKB PCs; injury count
#  - FinnGen summary-level check (R/31_prs_eval_finngen.R)
# Sample QC and within-group PCs use plink2 on the array genotypes; scoring uses plink2 --score
# on the imputed dosages with alleles pre-aligned in R. Every plink2 --out is asserted to lie
# under output/ so nothing is ever written next to the genotypes.

LONG_RANGE_LD_HG19 <- data.table(
  CHR = c(1, 2, 2, 2, 3, 3, 3, 5, 5, 5, 5, 6, 6, 6, 7, 8, 8, 8, 10, 11, 11, 12, 12, 20),
  START = c(48e6, 86e6, 134.5e6, 183e6, 47.5e6, 83.5e6, 89e6, 44.5e6, 98e6, 129e6, 135.5e6, 25.5e6, 57e6, 140e6, 55e6, 8e6, 43e6, 112e6, 37e6, 46e6, 87.5e6, 33e6, 109.5e6, 32e6),
  END = c(52e6, 100.5e6, 138e6, 190e6, 50e6, 87e6, 97.5e6, 50.5e6, 100.5e6, 132e6, 138.5e6, 33.5e6, 64e6, 142.5e6, 66e6, 12e6, 50e6, 115e6, 43e6, 57e6, 90.5e6, 40e6, 112e6, 34.5e6))

pe_plink2 <- function(config, args, stage, log_file, expect = character(0)) {
  i <- which(args == "--out")
  if (!length(i)) stop("pe_plink2: --out missing")
  out_dir <- normalizePath(dirname(args[i + 1]), mustWork = FALSE)
  if (!startsWith(out_dir, normalizePath(config$paths$output_dir, mustWork = FALSE)))
    stop(sprintf("pe_plink2: refusing to write outside output/: %s", args[i + 1]))
  cfg <- config$prs_eval %||% list()
  run_tool(tool_path(config, "plink2"), c(args, "--threads", as.character(cfg$plink_threads %||% 8),
                                           "--memory", as.character(cfg$plink_memory_mb %||% 32000)),
           stage, log_file, expect = expect)
}

.pe_write_keep <- function(ids, path) {
  fwrite(data.table(FID = ids, IID = ids), path, sep = "\t", col.names = FALSE)
  Sys.chmod(path, mode = "0600")
  path
}

# Sample QC per group: ethnicity, genetic-QC flags, relatedness (KING primary, GCTA lists
# sensitivity), within-group PCs. Returns list(group -> list(keep, keep_gcta, pcs, steps)).
pe_sample_qc <- function(config, ph, out_dir) {
  stage <- "prs_eval"
  cfg <- config$prs_eval; refs <- config$refs
  qd <- file.path(out_dir, "sample_qc"); dir.create(qd, recursive = TRUE, showWarnings = FALSE)
  pgen <- refs$ukb_asian_black_pgen_dir
  psam <- fread(file.path(pgen, "UKB_chr22.psam"), colClasses = "character")
  in_pgen <- as.character(psam[[2]])
  ph <- copy(ph)[eid %in% in_pgen]
  ph[, group := pe_assign_group(eth, cfg$groups)]
  sexnum <- function(x) suppressWarnings(as.integer(round(as.numeric(x))))
  res <- list()
  lrld <- file.path(qd, "long_range_ld_hg19.bed")
  fwrite(LONG_RANGE_LD_HG19[, .(CHR, START = as.integer(START), END = as.integer(END), LABEL = paste0("lrld", .I))], lrld, sep = "\t", col.names = FALSE)
  for (g in names(cfg$groups)) {
    d <- ph[group == g]
    steps <- list(in_pgen = length(in_pgen), in_group = nrow(d))
    ok <- rep(TRUE, nrow(d))
    sq <- cfg$sample_qc %||% list()
    if (isTRUE(sq$require_genetic_qc %||% TRUE)) { ok <- ok & !is.na(d$genetic_sex) & nzchar(d$genetic_sex); steps$has_genetic_qc <- sum(ok) }
    if (isTRUE(sq$require_sex_match %||% TRUE)) { ok <- ok & sexnum(d$sex) == sexnum(d$genetic_sex); steps$sex_match <- sum(ok) }
    if (isTRUE(sq$exclude_sex_aneuploidy %||% TRUE)) { ok <- ok & (is.na(d$aneuploidy) | !nzchar(d$aneuploidy)); steps$no_aneuploidy <- sum(ok) }
    if (isTRUE(sq$exclude_het_missing_outliers %||% TRUE)) { ok <- ok & (is.na(d$het_miss) | !nzchar(d$het_miss)); steps$no_het_missing_outlier <- sum(ok) }
    ok[is.na(ok)] <- FALSE
    basic <- d$eid[ok]
    steps$basic_keep <- length(basic)
    basic_file <- .pe_write_keep(basic, file.path(qd, paste0(g, ".basic_keep.txt")))
    work <- file.path(qd, "work"); dir.create(work, showWarnings = FALSE)
    # LD pruning on the array genotypes, long-range LD regions excluded
    pre <- file.path(work, paste0(g, "_prune"))
    if (!file.exists(paste0(pre, ".prune.in"))) {
      pe_plink2(config, c("--pfile", file.path(pgen, "UKB_genotyped"), "--keep", basic_file, "--maf", "0.01", "--geno", "0.02",
                          "--hwe", "1e-6", "--exclude", "bed1", lrld, "--indep-pairwise", "1000kb", "1", "0.1", "--out", pre),
                stage, paste0(pre, ".stdout.log"), expect = paste0(pre, ".prune.in"))
    }
    # relatedness: KING 3rd degree (primary)
    kp <- file.path(work, paste0(g, "_king"))
    if (!file.exists(paste0(kp, ".king.cutoff.in.id"))) {
      pe_plink2(config, c("--pfile", file.path(pgen, "UKB_genotyped"), "--keep", basic_file, "--extract", paste0(pre, ".prune.in"),
                          "--king-cutoff", as.character(cfg$relatedness$king_cutoff %||% 0.0442), "--out", kp),
                stage, paste0(kp, ".stdout.log"), expect = paste0(kp, ".king.cutoff.in.id"))
    }
    king_in <- fread(paste0(kp, ".king.cutoff.in.id"), colClasses = "character")
    keep_king <- intersect(basic, as.character(king_in[[ncol(king_in)]]))
    steps$unrelated_king <- length(keep_king)
    # GCTA singleton lists (sensitivity)
    gf <- refs[[paste0("ukb_unrelated_", tolower(gsub("([a-z])([A-Z])", "\\1_\\2", g)))]]
    keep_gcta <- if (!is.null(gf) && file.exists(gf)) intersect(basic, as.character(fread(gf, header = FALSE, colClasses = "character")[[2]])) else character(0)
    steps$unrelated_gcta_list <- length(keep_gcta)
    primary <- if (identical(cfg$relatedness$primary %||% "king", "gcta_list")) keep_gcta else keep_king
    keep_file <- .pe_write_keep(primary, file.path(qd, paste0(g, ".keep.txt")))
    if (length(keep_gcta)) .pe_write_keep(keep_gcta, file.path(qd, paste0(g, ".gcta_keep.txt")))
    # within-group PCs on the unrelated pruned set
    pp <- file.path(work, paste0(g, "_pca"))
    if (!file.exists(paste0(pp, ".eigenvec"))) {
      pe_plink2(config, c("--pfile", file.path(pgen, "UKB_genotyped"), "--keep", keep_file, "--extract", paste0(pre, ".prune.in"),
                          "--pca", as.character(cfg$pcs$n %||% 10), "--out", pp), stage, paste0(pp, ".stdout.log"), expect = paste0(pp, ".eigenvec"))
    }
    ev <- fread(paste0(pp, ".eigenvec"), colClasses = list(character = 1:2))
    pcs <- data.table(eid = as.character(ev[[2]]), ev[, -(1:2)]); setnames(pcs, -1, paste0("PC", seq_len(ncol(pcs) - 1)))
    steps$n_final <- length(primary)
    log_info(stage, sprintf("%s: %s", g, paste(sprintf("%s=%s", names(steps), unlist(steps)), collapse = ", ")))
    res[[g]] <- list(keep = primary, keep_gcta = keep_gcta, pcs = pcs, steps = steps)
  }
  jsonlite::write_json(lapply(res, function(r) r$steps), file.path(qd, "sample_qc.json"), auto_unbox = TRUE, pretty = TRUE)
  res
}

# Target variants (SNP CHR BP REF ALT) from the 22 .pvar files, cached.
pe_target_variants <- function(config, out_dir) {
  p <- file.path(out_dir, "weights", "target_variants.rds")
  if (file.exists(p)) return(readRDS(p))
  dir.create(dirname(p), recursive = TRUE, showWarnings = FALSE)
  pgen <- config$refs$ukb_asian_black_pgen_dir
  tv <- rbindlist(lapply(1:22, function(c) fread(file.path(pgen, sprintf("UKB_chr%d.pvar", c)), skip = "#CHROM",
                                                 select = 1:5, col.names = c("CHR", "BP", "SNP", "REF", "ALT"))))
  tv[, CHR := as.integer(CHR)]
  saveRDS(tv, paste0(p, ".tmp")); file.rename(paste0(p, ".tmp"), p)   # atomic: a killed run leaves no partial cache
  tv
}

# All weight sets: factor x method (C+T expanded per threshold) and comparators.
pe_load_weight_sets <- function(config, sex, cfg) {
  wdir <- file.path(stage_dir(config, sex, "prs"), "weights")
  methods <- unlist(cfg$weights$methods %||% c("ct", "prscs", "ldpred2_auto", "sbayesrc"))
  ext <- read_external_manifest(config)
  sources <- c(lapply(postgwas_factors(config), function(f) list(name = f, kind = "factor", factor = f)),
               if (!is.null(ext)) lapply(seq_len(nrow(ext)), function(i) list(name = ext$id[i], kind = "comparator", factor = "ARD", label = ext$label[i])))
  sets <- list()
  for (s in sources) for (m in methods) {
    f <- file.path(wdir, m, paste0(s$name, ".weights.tsv"))
    if (!file.exists(f)) next
    w <- fread(f)
    if (m == "ct") {
      rg <- fread(file.path(wdir, m, paste0(s$name, ".ranges")), header = FALSE, col.names = c("name", "lower", "upper"))
      for (k in seq_len(nrow(rg))) {
        sets[[length(sets) + 1]] <- c(s, list(prs_id = sprintf("%s__ct__%s", s$name, rg$name[k]), method = "ct", threshold = rg$name[k],
                                              file = f, mtime = as.numeric(file.mtime(f)), dt = w[P <= rg$upper[k], .(SNP, A1, A2, weight)]))
      }
    } else {
      sets[[length(sets) + 1]] <- c(s, list(prs_id = sprintf("%s__%s", s$name, m), method = m, threshold = NA_character_,
                                            file = f, mtime = as.numeric(file.mtime(f)), dt = w[, .(SNP, A1, A2, weight)]))
    }
  }
  sets
}

# Per-chromosome score matrices (SNP A1=ALT <prs_id...>), alleles aligned to the target .pvar.
pe_build_score_matrix <- function(sets, target, out_dir) {
  md <- file.path(out_dir, "weights", "score_matrix"); dir.create(md, recursive = TRUE, showWarnings = FALSE)
  cov <- list(); aligned <- list()
  tgt <- target[, .(SNP, CHR, BP, REF, ALT)]; setkey(tgt, SNP)   # keyed once; joins are then ~1 s per set
  for (s in sets) {
    w <- tgt[s$dt, on = "SNP", nomatch = NULL]
    ok_same <- w$A1 == w$ALT & w$A2 == w$REF
    ok_flip <- w$A1 == w$REF & w$A2 == w$ALT
    w <- w[ok_same | ok_flip]
    w[, weight := fifelse(A1 == REF, -weight, weight)]
    cov[[s$prs_id]] <- data.table(prs_id = s$prs_id, kind = s$kind, name = s$name, method = s$method, threshold = s$threshold,
                                  n_weight = nrow(s$dt), n_aligned = nrow(w), pct_found = round(100 * nrow(w) / max(1, nrow(s$dt)), 1),
                                  n_sign_flipped = sum(ok_flip[ok_same | ok_flip]))
    aligned[[s$prs_id]] <- w[weight != 0, .(SNP, CHR, weight)]
  }
  ids <- names(aligned)
  for (c in 1:22) {
    base <- target[CHR == c, .(SNP, A1 = ALT)]
    use <- rbindlist(lapply(ids, function(id) aligned[[id]][CHR == c, .(SNP, prs_id = id, weight)]))
    if (!nrow(use)) next
    wide <- dcast(use, SNP ~ prs_id, value.var = "weight", fill = 0)
    for (id in setdiff(ids, names(wide))) wide[, (id) := 0]
    setcolorder(wide, c("SNP", ids))
    wide <- merge(base, wide, by = "SNP")
    fwrite(wide, file.path(md, sprintf("chr%d.tsv", c)), sep = "\t")
    writeLines(wide$SNP, file.path(md, sprintf("chr%d.snps", c)))
  }
  covt <- rbindlist(cov)
  covt[, column_index := match(prs_id, ids) + 2L]
  fwrite(covt, file.path(out_dir, "weights", "coverage.tsv"), sep = "\t")
  fwrite(covt[, .(prs_id, kind, name, method, threshold, column_index)], file.path(out_dir, "weights", "score_sets.tsv"), sep = "\t")
  covt
}

# A scored chromosome is reused only while its score matrix is byte-identical to the one that was
# scored (md5 kept in the .done file), so weight sets added by a later prs run are never skipped.
.pe_checkpoint_ok <- function(done, mat) {
  if (!file.exists(done) || !file.exists(mat)) return(FALSE)
  h <- readLines(done, n = 1, warn = FALSE)
  length(h) == 1 && identical(unname(tools::md5sum(mat)), h)
}

# plink2 --score per chromosome (dosages), summed to one row per participant.
pe_score <- function(config, out_dir, keep_ids, n_sets) {
  stage <- "prs_eval"
  cfg <- config$prs_eval
  sd <- file.path(out_dir, "scores", "work"); dir.create(sd, recursive = TRUE, showWarnings = FALSE)
  keep_file <- .pe_write_keep(keep_ids, file.path(out_dir, "scores", "union_keep.txt"))
  md <- file.path(out_dir, "weights", "score_matrix")
  pgen <- config$refs$ukb_asian_black_pgen_dir
  one <- function(c) {
    pre <- file.path(sd, sprintf("chr%d", c))
    mat <- file.path(md, sprintf("chr%d.tsv", c))
    if (!file.exists(mat)) { unlink(paste0(pre, c(".done", ".sscore"))); return(FALSE) }
    if (.pe_checkpoint_ok(paste0(pre, ".done"), mat)) return(TRUE)
    unlink(paste0(pre, c(".done", ".sscore")))
    pe_plink2(config, c("--pfile", file.path(pgen, sprintf("UKB_chr%d", c)), "--keep", keep_file,
                        "--extract", file.path(md, sprintf("chr%d.snps", c)),
                        "--score", mat, "1", "2", "header-read", "cols=+scoresums,-scoreavgs", "list-variants", "ignore-dup-ids",
                        "--score-col-nums", sprintf("3-%d", n_sets + 2L), "--out", pre),
              stage, paste0(pre, ".stdout.log"), expect = paste0(pre, ".sscore"))
    writeLines(unname(tools::md5sum(mat)), paste0(pre, ".done")); TRUE
  }
  r <- parallel::mclapply(1:22, function(c) tryCatch(one(c), error = function(e) e), mc.cores = as.integer(cfg$score_parallel_chr %||% 4L))
  bad <- vapply(r, inherits, logical(1), "error")
  if (any(bad)) log_fatal(stage, sprintf("scoring failed: %s", conditionMessage(r[bad][[1]])))
  tot <- NULL
  for (c in 1:22) {
    f <- file.path(sd, sprintf("chr%d.sscore", c))
    if (!file.exists(f)) next
    x <- fread(f, colClasses = list(character = 1:2))
    x <- x[, c(2, grep("_SUM$", names(x))), with = FALSE]; setnames(x, 1, "eid")
    tot <- if (is.null(tot)) x else { m <- merge(tot, x, by = "eid", suffixes = c("", ".y")); for (nm in setdiff(names(x), "eid")) { set(m, j = nm, value = m[[nm]] + m[[paste0(nm, ".y")]]); set(m, j = paste0(nm, ".y"), value = NULL) }; m }
  }
  setnames(tot, sub("_SUM$", "", names(tot)))
  fwrite(tot, file.path(out_dir, "scores", "raw_sums.tsv"), sep = "\t")
  tot
}

# Models for one group. Returns list(results, deciles).
pe_evaluate_group <- function(config, g, pheno, raw, qc, sets, cfg, gbd = TRUE) {
  stage <- "prs_eval"
  gi <- qc[[g]]
  ids <- gi$keep
  res <- list(); dec <- list()
  n_gpc <- as.integer(cfg$pcs$sensitivity_global_pcs %||% 4L)
  one_sample <- function(sample_ids, sample_set, covariate_set) {
    dd <- pheno[eid %in% sample_ids]
    if (covariate_set != "global_pcs") dd <- merge(dd, gi$pcs, by = "eid")
    dd <- merge(dd, raw[eid %in% sample_ids], by = "eid")
    cov_cols <- if (covariate_set == "global_pcs") paste0("gpc", seq_len(n_gpc)) else grep("^PC[0-9]+$", names(dd), value = TRUE)
    if (covariate_set == "contact") {   # + log(1 + distinct hospital diagnosis dates): hospital-contact adjustment
      dd[, log_hes_dates := log1p(n_hes_dates)]
      cov_cols <- c(cov_cols, "log_hes_dates")
    }
    X <- as.data.frame(dd[, c("age", "age2", "sex", cov_cols), with = FALSE])
    for (s in sets) {
      if (!s$prs_id %in% names(dd)) next
      z <- as.numeric(scale(dd[[s$prs_id]]))
      base <- data.table(group = g, prs_id = s$prs_id, kind = s$kind, name = s$name, method = s$method, threshold = s$threshold,
                         factor = s$factor, sample_set = sample_set, covariate_set = covariate_set)
      f <- s$factor
      outs <- list(list(outcome = paste0("count_", f), unit = "icd3", type = "count"))
      if (gbd && paste0("count_", f, "_gbd") %in% names(dd)) outs <- c(outs, list(list(outcome = paste0("count_", f, "_gbd"), unit = "gbd_cause", type = "count")))
      outs <- c(outs, list(list(outcome = "mortality", unit = "death", type = "cox")))
      if (f == "ARD") outs <- c(outs, list(list(outcome = "injury_count", unit = "icd3", type = "count")))
      for (o in outs) {
        r <- if (o$type == "count") fit_count(dd[[o$outcome]], z, X, min_events = as.integer(cfg$min_events %||% 20L))
             else fit_cox_age(dd$age_entry, dd$age_exit, dd$death, z, X[, setdiff(names(X), c("age", "age2")), drop = FALSE],
                              min_events = as.integer(cfg$min_events %||% 20L))
        res[[length(res) + 1]] <<- cbind(base, data.table(outcome = o$outcome, unit = o$unit), r)
      }
      if (sample_set == "king" && covariate_set == "primary" && f == "ARD")
        dec[[length(dec) + 1]] <<- cbind(data.table(group = g, prs_id = s$prs_id), decile_table(z, dd[[paste0("count_", f)]]))
    }
  }
  one_sample(ids, "king", "primary")
  if (length(gi$keep_gcta)) one_sample(gi$keep_gcta, "gcta_list", "primary")
  if (any(grepl("^gpc", names(pheno)))) one_sample(ids, "king", "global_pcs")
  if ("n_hes_dates" %in% names(pheno)) one_sample(ids, "king", "contact")
  list(results = rbindlist(res, fill = TRUE), deciles = rbindlist(dec))
}

run_prs_eval <- function(config, sex) {
  stage <- "prs_eval"
  log_info(stage, sprintf("=== PRS evaluation (UKB South Asian / Black): %s ===", sex))
  if (sex != "bothsex") { log_warn(stage, "PRS evaluation runs on the both-sex factors only"); return(invisible(NULL)) }
  cfg <- config$prs_eval %||% list(); refs <- config$refs
  out_dir <- stage_dir(config, sex, stage); dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  for (k in c("ukb_asian_black_pgen_dir", "ukb_pheno_extract")) if (!file.exists(refs[[k]] %||% "")) log_fatal(stage, sprintf("refs$%s missing", k))
  reg <- read_factor_registry(config, sex)
  indicators <- lapply(setNames(nm = postgwas_factors(config)), function(f) reg[[f]]$indicators)
  # phenotypes (cached)
  pd <- file.path(out_dir, "pheno"); dir.create(pd, showWarnings = FALSE)
  ids <- as.character(fread(file.path(refs$ukb_asian_black_pgen_dir, "UKB_chr22.psam"), colClasses = "character")[[2]])
  ph <- pe_read_pheno_extract(refs$ukb_pheno_extract, ids)
  gbd_path <- cfg$gbd_map %||% file.path(config$project$root, "manifest", "indicator_gbd_causes.tsv")
  gbd <- if (file.exists(gbd_path)) pe_read_gbd_map(gbd_path) else { log_warn(stage, "GBD roll-up manifest missing; sensitivity unit skipped"); NULL }
  pheno <- pe_build_pheno(ph, indicators, gbd, censor = cfg$censor %||% list(hes = "2022-10-31", death = "2022-12-19"),
                          injury_chapters = unlist(cfg$injury_chapters %||% "S"))
  pheno <- merge(pheno, ph[, c("eid", grep("^gpc", names(ph), value = TRUE)), with = FALSE], by = "eid")
  for (col in grep("^gpc", names(pheno), value = TRUE)) set(pheno, j = col, value = suppressWarnings(as.numeric(pheno[[col]])))
  saveRDS(pheno, file.path(pd, "pheno_all.rds"))
  log_info(stage, sprintf("phenotypes for %d participants (%d with HES records, %d deaths)", nrow(pheno), sum(pheno$any_hes), sum(pheno$death)))
  # weight sets and score matrices (no plink2 needed)
  sets <- pe_load_weight_sets(config, sex, cfg)
  if (!length(sets)) log_fatal(stage, "no weight files under prs/weights (run the prs stage)")
  log_info(stage, sprintf("%d weight sets (%d factor, %d comparator)", length(sets), sum(vapply(sets, function(s) s$kind == "factor", logical(1))),
                          sum(vapply(sets, function(s) s$kind == "comparator", logical(1)))))
  target <- pe_target_variants(config, out_dir)
  covt <- pe_build_score_matrix(sets, target, out_dir)
  for (i in seq_len(nrow(covt))) if (covt$pct_found[i] < 80) log_warn(stage, sprintf("%s: only %.1f%% of weight SNPs found in the target", covt$prs_id[i], covt$pct_found[i]))
  # FinnGen summary-level check needs only plink 1.9 and downloads, so it runs before the plink2 gate
  if (isTRUE((cfg$finngen %||% list())$run %||% FALSE)) {
    fr <- tryCatch(run_prs_eval_finngen(config, sex, sets, out_dir), error = function(e) { log_warn(stage, sprintf("FinnGen module failed: %s", conditionMessage(e))); NULL })
  }
  # plink2-dependent part
  if (!nzchar(tool_path(config, "plink2", required = FALSE)))
    log_fatal(stage, "plink2 not installed (scripts/install_postgwas_tools.sh); phenotypes, score matrices and the FinnGen results are cached, rerun resumes here")
  qc <- pe_sample_qc(config, ph, out_dir)
  raw <- pe_score(config, out_dir, unique(unlist(lapply(qc, function(q) union(q$keep, q$keep_gcta)))), nrow(covt))
  # models
  all_res <- list(); all_dec <- list()
  for (g in names(qc)) {
    e <- pe_evaluate_group(config, g, pheno, raw, qc, sets, cfg, gbd = !is.null(gbd))
    all_res[[g]] <- e$results; all_dec[[g]] <- e$deciles
    fwrite(e$results, file.path(out_dir, sprintf("results_%s.tsv", g)), sep = "\t")
    fwrite(e$deciles, file.path(out_dir, sprintf("decile_%s.tsv", g)), sep = "\t")
    z <- merge(pheno[eid %in% qc[[g]]$keep, .(eid)], raw, by = "eid")
    for (col in setdiff(names(z), "eid")) set(z, j = col, value = as.numeric(scale(z[[col]])))
    fwrite(z, file.path(out_dir, "scores", sprintf("%s.tsv", g)), sep = "\t")
  }
  res <- rbindlist(all_res, fill = TRUE)
  prim <- res[sample_set == "king" & covariate_set == "primary" & unit == "icd3" & outcome == paste0("count_", factor)]
  fwrite(ct_select_threshold(prim), file.path(out_dir, "ct_selected_thresholds.tsv"), sep = "\t")
  fwrite(prim[, .(group, factor, kind, name, method, threshold, prs_id, n, irr = effect, ci_low, ci_high, p, lr_chisq, pseudo_r2)][order(factor, group, -lr_chisq)],
         file.path(out_dir, "method_comparison.tsv"), sep = "\t")
  fwrite(res[kind == "comparator" & sample_set == "king" & covariate_set == "primary"], file.path(out_dir, "comparators.tsv"), sep = "\t")
  gs <- names(qc)
  if (length(gs) == 2) {
    a <- res[group == gs[1]]; b <- res[group == gs[2]]
    m <- merge(a[, .(prs_id, outcome, unit, sample_set, covariate_set, e1 = estimate, s1 = se)],
               b[, .(prs_id, outcome, unit, sample_set, covariate_set, e2 = estimate, s2 = se)],
               by = c("prs_id", "outcome", "unit", "sample_set", "covariate_set"))
    m <- m[is.finite(e1) & is.finite(e2)]
    meta <- cbind(m[, .(prs_id, outcome, unit, sample_set, covariate_set)], ivw_meta(m$e1, m$s1, m$e2, m$s2))
    meta[, note := "supplementary inverse-variance meta of the two groups; groups are reported separately"]
    fwrite(meta, file.path(out_dir, "results_meta.tsv"), sep = "\t")
  }
  for (g in gs) {
    h <- prim[group == g & kind == "factor" & factor == "ARD"][order(-lr_chisq)][1]
    if (nrow(h) && !is.na(h$effect)) log_info(stage, sprintf("%s: ARD count, best ARD score %s: IRR %.3f (%.3f-%.3f) per SD, p = %.2g, n = %d",
                                                               g, h$prs_id, h$effect, h$ci_low, h$ci_high, h$p, h$n))
  }
  clear_postgwas_cache()
  write_stage_manifest(stage, sex, config, file.path(out_dir, c(sprintf("results_%s.tsv", gs), "method_comparison.tsv", "comparators.tsv")))
  invisible(res)
}
