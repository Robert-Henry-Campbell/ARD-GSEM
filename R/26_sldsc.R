# Stratified LDSC stage (GWASTutorial "LD score regression by LDSC": partitioned heritability
# and LDSC-SEG), on the factor_ldsc qfilter sumstats (N = N_hat_std; enrichments do not
# depend on N):
#  - partitioned h2 with the baselineLD v2.2 model (enrichment and tau; FDR over annotations)
#  - LDSC-SEG (--h2-cts) on GTEx v8 tissue-specific expression annotations built here (the
#    official Multi_tissue_gene_expr scores are only in a requester-pays bucket): for each tissue
#    the top 10% of genes by specificity (TPM in tissue / sum over tissues), GENCODE v19 gene
#    bodies +-100 kb; control = every gene in the matrix. LD scores are computed once per
#    chromosome for all tissues together (1000G Phase 3 EUR, 1 cM window, HapMap3 SNPs printed)
#    and split into one file set per tissue. Conditioned on baseline v1.2; FDR over tissues.

ldsc_py <- function(config, args, stage, log_file, expect = character(0)) {
  py <- tool_path(config, "ldsc_python")
  run_tool(py, c(file.path(config$tools$ldsc_dir, "ldsc.py"), args), stage, log_file, expect = expect)
}

# Enrichment is defined only for binary annotations; continuous ones (GERP, MaxCPP, LLD, ...)
# are read from their tau. Classified from the chr22 .annot values.
baseline_binary <- function(prefix) {
  a <- fread(paste0(prefix, "22.annot.gz"))
  a <- a[, !intersect(c("CHR", "BP", "SNP", "CM"), names(a)), with = FALSE]
  vapply(a, function(v) all(v %in% c(0, 1)), logical(1))
}

# re-use an ldsc result only when it is newer than its input sumstats
stale <- function(out, input) !file.exists(out) || file.mtime(out) < file.mtime(input)

build_gtex_ldsc_seg <- function(config, ann_dir, top_frac = 0.1, flank = 100000L, workers = 22L) {
  stage <- "sldsc"
  done <- file.path(ann_dir, "gtex.ldcts")
  if (file.exists(done)) return(done)
  dir.create(ann_dir, recursive = TRUE, showWarnings = FALSE)
  tpm <- fread(cmd = paste("zcat", shQuote(config$refs$gtex_median_tpm)), skip = 2)
  genes <- read_gencode_genes(config$refs$gencode_gtf)
  tpm <- tpm[Description %in% genes$symbol][!duplicated(Description)]
  m <- as.matrix(tpm[, !c("Name", "Description")]); rownames(m) <- tpm$Description
  m <- m[rowSums(m) > 0, , drop = FALSE]
  spec <- m / rowSums(m)
  tissues <- sub("_+$", "", gsub("[^A-Za-z0-9]+", "_", colnames(m)))
  sets <- lapply(seq_len(ncol(m)), function(j) rownames(m)[spec[, j] >= quantile(spec[, j], 1 - top_frac)])
  names(sets) <- tissues
  sets$control <- rownames(m)
  g <- genes[symbol %in% rownames(m), .(symbol, CHR, lo = pmax(1L, start - flank), hi = end + flank)]
  plink_pre <- config$refs$ldsc_plink
  one_chr <- function(c) {
    out <- file.path(ann_dir, sprintf("all.%d", c))
    if (file.exists(paste0(out, ".l2.M_5_50"))) return(c)   # written last by ldsc --l2
    bim <- fread(sprintf("%s%d.bim", plink_pre, c), header = FALSE, col.names = c("CHR", "SNP", "CM", "BP", "A1", "A2"))
    ann <- as.data.table(lapply(sets, function(s) {
      w <- g[CHR == c & symbol %in% s]
      hit <- if (nrow(w)) unique(bim[, .(i = .I, BP)][w, on = .(BP >= lo, BP <= hi), nomatch = NULL, x.i]) else integer(0)
      v <- integer(nrow(bim)); v[hit] <- 1L; v
    }))
    fwrite(ann, paste0(out, ".annot.gz"), sep = "\t")
    ldsc_py(config, c("--l2", "--bfile", sprintf("%s%d", plink_pre, c), "--ld-wind-cm", "1",
                      "--annot", paste0(out, ".annot.gz"), "--thin-annot", "--out", out,
                      "--print-snps", config$refs$ldsc_hm3_list),
            stage, paste0(out, ".l2.log"), expect = paste0(out, ".l2.M_5_50"))
    c
  }
  r <- parallel::mclapply(1:22, function(c) tryCatch(one_chr(c), error = function(e) e), mc.cores = workers)
  bad <- vapply(r, inherits, logical(1), "error")
  if (any(bad)) log_fatal(stage, sprintf("LD score computation failed: %s", conditionMessage(r[bad][[1]])))
  # split the multi-annotation LD scores into one prefix per annotation
  for (c in 1:22) {
    src <- file.path(ann_dir, sprintf("all.%d", c))
    ld <- fread(paste0(src, ".l2.ldscore.gz"))
    m_all <- scan(paste0(src, ".l2.M"), quiet = TRUE)
    m550 <- scan(paste0(src, ".l2.M_5_50"), quiet = TRUE)
    for (k in seq_along(sets)) {
      nm <- names(sets)[k]
      dir.create(file.path(ann_dir, nm), showWarnings = FALSE)
      pre <- file.path(ann_dir, nm, sprintf("%s.%d", nm, c))
      col <- paste0(nm, "L2")
      if (!col %in% names(ld)) log_fatal(stage, sprintf("%s missing from %s.l2.ldscore.gz", col, src))
      fwrite(ld[, c("CHR", "SNP", "BP", col), with = FALSE][, setnames(.SD, col, "L2")], paste0(pre, ".l2.ldscore.gz"), sep = "\t")
      writeLines(as.character(m_all[k]), paste0(pre, ".l2.M"))
      writeLines(as.character(m550[k]), paste0(pre, ".l2.M_5_50"))
    }
  }
  ctl <- file.path(ann_dir, "control", "control.")
  writeLines(sprintf("%s\t%s,%s", tissues, file.path(ann_dir, tissues, paste0(tissues, ".")), ctl), done)
  done
}

run_sldsc <- function(config, sex) {
  stage <- "sldsc"
  log_info(stage, sprintf("=== Stratified LDSC: %s ===", sex))
  cfg <- config$sldsc %||% list()
  refs <- config$refs
  out_dir <- stage_dir(config, sex, stage); dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  ldcts <- build_gtex_ldsc_seg(config, file.path(out_dir, "gtex_annotations"),
                               as.numeric(cfg$top_frac %||% 0.1), as.integer(cfg$flank_kb %||% 100L) * 1000L,
                               as.integer(cfg$workers %||% 22L))
  is_binary <- baseline_binary(refs$ldsc_baselineld)
  # ldsc is single-threaded: the factors run in parallel
  one_factor <- function(f) {
    ss <- file.path(stage_dir(config, sex, "factor_ldsc"), sprintf("%s_qfilter.sumstats.gz", f))
    if (!file.exists(ss)) stop(sprintf("%s missing (run factor_ldsc)", ss))
    pre <- file.path(out_dir, paste0(f, "_baselineLD"))
    if (stale(paste0(pre, ".results"), ss)) {
      ldsc_py(config, c("--h2", ss, "--ref-ld-chr", refs$ldsc_baselineld, "--w-ld-chr", refs$ldsc_weights,
                        "--overlap-annot", "--frqfile-chr", refs$ldsc_frq, "--print-coefficients", "--out", pre),
              stage, paste0(pre, ".stdout.log"), expect = paste0(pre, ".results"))
    }
    r <- fread(paste0(pre, ".results"))
    r[, `:=`(factor = f, annotation = sub("L2_0$", "", Category))]
    r[, binary := unname(is_binary[annotation])]
    r[, tau_p := 2 * pnorm(-abs(`Coefficient_z-score`))]
    r[annotation != "base", tau_q := p.adjust(tau_p, "BH")]
    r[binary == TRUE & annotation != "base", enrichment_q := p.adjust(Enrichment_p, "BH")]
    r[binary == FALSE, `:=`(Enrichment = NA_real_, Enrichment_std_error = NA_real_, Enrichment_p = NA_real_)]
    pc <- file.path(out_dir, paste0(f, "_gtex"))
    if (stale(paste0(pc, ".cell_type_results.txt"), ss)) {
      ldsc_py(config, c("--h2-cts", ss, "--ref-ld-chr", refs$ldsc_baseline, "--out", pc,
                        "--ref-ld-chr-cts", ldcts, "--w-ld-chr", refs$ldsc_weights),
              stage, paste0(pc, ".stdout.log"), expect = paste0(pc, ".cell_type_results.txt"))
    }
    ct <- fread(paste0(pc, ".cell_type_results.txt"))
    ct[, `:=`(factor = f, q_bh = p.adjust(Coefficient_P_value, "BH"))]
    log_info(stage, sprintf("%s: %d binary baselineLD annotations enriched and %d with non-zero tau (FDR<5%%); top GTEx tissue %s (p = %.2g; %d FDR<5%%)",
                            f, sum(r$enrichment_q < 0.05, na.rm = TRUE), sum(r$tau_q < 0.05, na.rm = TRUE),
                            ct$Name[which.min(ct$Coefficient_P_value)], min(ct$Coefficient_P_value), sum(ct$q_bh < 0.05)))
    list(res = r, cts = ct)
  }
  factors <- postgwas_factors(config)
  out <- parallel::mclapply(factors, function(f) tryCatch(one_factor(f), error = function(e) e),
                            mc.cores = length(factors))
  bad <- vapply(out, inherits, logical(1), "error")
  if (any(bad)) log_fatal(stage, sprintf("%s: %s", factors[bad][1], conditionMessage(out[bad][[1]])))
  res <- rbindlist(lapply(out, `[[`, "res")); cts <- rbindlist(lapply(out, `[[`, "cts"))
  fwrite(res, file.path(out_dir, "baselineLD_results.tsv"), sep = "\t")
  fwrite(cts, file.path(out_dir, "gtex_ldsc_seg_results.tsv"), sep = "\t")
  write_stage_manifest(stage, sex, config, file.path(out_dir, c("baselineLD_results.tsv", "gtex_ldsc_seg_results.tsv")))
  invisible(cts)
}
