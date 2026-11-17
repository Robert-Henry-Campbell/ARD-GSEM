# Single-cell stage (GWASTutorial "Single-cell genomics integration"), on MAGMA gene results:
#  - MAGMA gene-property analysis with Tabula Sapiens cell-type mean expression (log1p CPM),
#    conditioned on the average, one-sided; Bonferroni over cell types
#  - scDRS: top 1,000 MAGMA genes (z-score weights) scored per cell on a seeded Tabula Sapiens
#    subsample (<= max_per_type cells per cell_ontology_class), then cell-type association
#    (Monte-Carlo p, FDR over cell types)
# (LDSC-SEG runs in the sldsc stage.)

# compute-score is the slow step (~2 h per factor on one core), so a rerun (e.g. after a downstream
# failure) reuses its output when nothing feeding it has changed: the log shows the same options (the
# CLI echoes each one under "Call:") and "Trait=<factor>, ..." for every factor (printed only after both
# score files are written), the gene sets just munged from MAGMA are unchanged, and the scores are newer
# than the Tabula Sapiens inputs.
scdrs_scores_current <- function(args, factors, log_file, inputs, gs_prev, gs_now) {
  out <- args[match("--out-folder", args) + 1L]
  sf <- file.path(out, paste0(rep(factors, each = 2L), c(".score.gz", ".full_score.gz")))
  if (is.null(gs_prev) || !identical(gs_prev, gs_now) || !all(file.exists(c(sf, log_file, inputs)))) return(FALSE)
  lg <- sub("\\s*\\\\$", "", readLines(log_file, warn = FALSE))
  opts <- args[-1L]
  same_call <- all(paste(opts[c(TRUE, FALSE)], opts[c(FALSE, TRUE)]) %in% lg)
  scored <- all(vapply(factors, function(f) any(grepl(paste0("Trait=", f, ", n_gene="), lg, fixed = TRUE)), logical(1)))
  same_call && scored && min(file.mtime(sf)) > max(file.mtime(inputs))
}

run_scdrs <- function(config, sex) {
  stage <- "scdrs"
  log_info(stage, sprintf("=== Single-cell integration: %s ===", sex))
  cfg <- config$scdrs %||% list()
  out_dir <- stage_dir(config, sex, stage); dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  py_sc <- tool_path(config, "scdrs_python")
  scdrs <- tool_path(config, "scdrs")
  magma <- tool_path(config, "magma")
  ts_dir <- file.path(out_dir, "tabula_sapiens")
  ts_out <- file.path(ts_dir, c("ts_subsample.h5ad", "ts_cov.tsv", "ts_celltype_mean_log1pcpm.tsv"))
  if (!all(file.exists(ts_out))) {   # all three outputs, not just the h5ad (a crash once left a partial one)
    run_tool(py_sc, c(file.path(config$project$root, "scripts", "py", "ts_prepare.py"), "--h5ad", config$refs$tabula_sapiens_h5ad,
                      "--out-dir", ts_dir, "--max-per-type", as.character(cfg$max_per_type %||% 1000),
                      "--seed", as.character(cfg$seed %||% 20261108)),
             stage, file.path(out_dir, "ts_prepare.log"), expect = ts_out)
  }
  # MAGMA cell-type covariates (symbols -> Entrez via the MAGMA gene locations)
  loc <- fread(config$refs$magma_gene_loc, header = FALSE, col.names = c("entrez", "chr", "start", "end", "strand", "symbol"))
  ct <- fread(file.path(ts_dir, "ts_celltype_mean_log1pcpm.tsv"))
  setnames(ct, 1, "symbol")
  ct <- ct[symbol %in% loc$symbol][!duplicated(symbol)]
  m <- as.matrix(ct[, -1])
  cv <- data.table(GENE = loc$entrez[match(ct$symbol, loc$symbol)], m)
  setnames(cv, make.names(names(cv), unique = TRUE))
  cv[, Average := rowMeans(m)]
  covar <- file.path(out_dir, "ts_celltype_covar.txt")
  fwrite(cv, covar, sep = " ")
  res_ct <- list(); res_sc <- list()
  zfile <- file.path(out_dir, "magma_z.tsv")
  z <- NULL
  for (f in postgwas_factors(config)) {
    raw <- file.path(stage_dir(config, sex, "magma"), paste0(f, ".genes.raw"))
    gt <- file.path(stage_dir(config, sex, "magma"), paste0(f, "_genes.tsv"))
    if (!file.exists(raw) || !file.exists(gt)) log_fatal(stage, sprintf("MAGMA results for %s missing (run magma)", f))
    pre <- file.path(out_dir, paste0(f, "_ts"))
    run_tool(magma, c("--gene-results", raw, "--gene-covar", covar, "--model", "condition-hide=Average",
                      "direction=greater", "--out", pre), stage, paste0(pre, ".stdout.log"), expect = paste0(pre, ".gsa.out"))
    r <- magma_read_out(paste0(pre, ".gsa.out"))
    r[, `:=`(factor = f, sig_bonf = P < 0.05 / .N, q_bh = p.adjust(P, "BH"))]
    res_ct[[f]] <- r
    g <- fread(gt)[!is.na(SYMBOL), .(GENE = SYMBOL, z = ZSTAT)]
    setnames(g, "z", f)
    z <- if (is.null(z)) g else merge(z, g, by = "GENE", all = TRUE)
  }
  fwrite(z, zfile, sep = "\t", na = "NA")
  gs <- file.path(out_dir, "magma_top1000.gs")
  gs_prev <- if (file.exists(gs)) readLines(gs, warn = FALSE)
  run_tool(scdrs, c("munge-gs", "--out-file", gs, "--zscore-file", zfile, "--weight", "zscore",
                    "--n-max", as.character(cfg$n_genes %||% 1000)), stage, file.path(out_dir, "munge_gs.log"), expect = gs)
  sc_out <- file.path(out_dir, "scdrs"); dir.create(sc_out, showWarnings = FALSE)
  h5 <- file.path(ts_dir, "ts_subsample.h5ad")
  cov_file <- file.path(ts_dir, "ts_cov.tsv")
  cs_args <- c("compute-score", "--h5ad-file", h5, "--h5ad-species", "human", "--gs-file", gs,
               "--gs-species", "human", "--cov-file", cov_file, "--flag-filter-data", "True",
               "--flag-raw-count", "True", "--n-ctrl", as.character(cfg$n_ctrl %||% 1000),
               "--flag-return-ctrl-raw-score", "False", "--flag-return-ctrl-norm-score", "True",
               "--out-folder", sc_out)
  cs_log <- file.path(out_dir, "compute_score.log")
  if (scdrs_scores_current(cs_args, postgwas_factors(config), cs_log, c(h5, cov_file), gs_prev, readLines(gs, warn = FALSE))) {
    log_info(stage, "compute-score: reusing the last run's scores (same call and gene sets, every factor scored)")
  } else {
    run_tool(scdrs, cs_args, stage, cs_log)
  }
  for (f in postgwas_factors(config)) {
    sf <- file.path(sc_out, paste0(f, ".full_score.gz"))
    run_tool(scdrs, c("perform-downstream", "--h5ad-file", h5, "--score-file", sf, "--out-folder", sc_out,
                      "--group-analysis", "cell_ontology_class", "--flag-filter-data", "True", "--flag-raw-count", "True"),
             stage, file.path(out_dir, paste0(f, "_downstream.log")))
    grp <- file.path(sc_out, paste0(f, ".scdrs_group.cell_ontology_class"))
    if (file.exists(grp)) {
      x <- fread(grp)
      setnames(x, 1, "cell_type")
      x[, `:=`(factor = f, assoc_q = p.adjust(assoc_mcp, "BH"))]
      res_sc[[f]] <- x
      log_info(stage, sprintf("%s: %d cell types FDR<5%% (scDRS); %d Bonferroni (MAGMA cell type, top %s)", f,
                              sum(x$assoc_q < 0.05), sum(res_ct[[f]]$sig_bonf), res_ct[[f]][which.min(P), VARIABLE]))
    }
  }
  fwrite(rbindlist(res_ct), file.path(out_dir, "magma_celltype_results.tsv"), sep = "\t")
  fwrite(rbindlist(res_sc, fill = TRUE), file.path(out_dir, "scdrs_celltype_results.tsv"), sep = "\t")
  write_stage_manifest(stage, sex, config, file.path(out_dir, c("magma_celltype_results.tsv", "scdrs_celltype_results.tsv")))
  invisible(rbindlist(res_sc, fill = TRUE))
}
