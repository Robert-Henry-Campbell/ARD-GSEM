# TWAS stage (GWASTutorial "TWAS (using FUSION)"): FUSION.assoc_test.R with GTEx v8 expression
# weights and the 1000G EUR LDREF, on the qfilter tier (SNP A1 A2 Z, matched by rsID), every
# tissue x chromosome in parallel; Bonferroni over all gene-tissue models tested per factor;
# FUSION.post_process.R (joint / conditional) on tissue-chromosomes with significant genes.
# GTEx v8 FUSION weight positions are checked against GENCODE v19 and the .pos gene windows
# lifted to hg19 when they are GRCh38 (post_process groups loci on LDREF hg19 positions).

fusion_pos_files <- function(config, out_dir) {
  wdir <- config$refs$fusion_weights_dir
  pos <- list.files(wdir, pattern = "\\.pos$", recursive = TRUE, full.names = TRUE)
  # the standard panels only: .nofilter.pos repeats every heritable model and adds
  # non-heritable ones (double-counted tests); .hg19.pos are this function's own output
  pos <- pos[!grepl("\\.(hg19|nofilter)\\.pos$", pos)]
  if (!length(pos)) log_fatal("twas", sprintf("no FUSION .pos files under %s", wdir))
  # build check on the first panel: GTEx v8 positions vs GENCODE v19 gene starts
  p1 <- fread(pos[1])
  genes <- read_gencode_genes(config$refs$gencode_gtf)
  sym <- if (all(startsWith(p1$ID, "ENSG"))) genes$gene_id else genes$symbol
  key <- sub("\\..*$", "", p1$ID)
  d <- p1[, .(P0, g = genes$start[match(key, sym)])][!is.na(g)]
  hg38 <- nrow(d) > 20 && mean(abs(d$P0 - d$g) < 5000) < 0.5
  # WGT paths inside a .pos file are relative to the directory holding it
  if (!hg38) return(data.table(pos = pos, wdir = dirname(pos)))
  log_info("twas", "FUSION weight positions are GRCh38: lifting .pos gene windows to hg19")
  chain <- file.path(config$paths$reference_dir, "chains", "hg38ToHg19.over.chain")
  out <- vapply(pos, function(p) {
    q <- sub("\\.pos$", ".hg19.pos", p)
    if (!file.exists(q)) {
      x <- fread(p)
      x[, `:=`(P0 = lift_positions(CHR, P0, chain, out_dir), P1 = lift_positions(CHR, P1, chain, out_dir))]
      fwrite(x[!is.na(P0) & !is.na(P1)], q, sep = "\t")
    }
    q
  }, character(1))
  data.table(pos = unname(out), wdir = dirname(pos))
}

run_twas <- function(config, sex) {
  stage <- "twas"
  log_info(stage, sprintf("=== TWAS (FUSION, GTEx v8): %s ===", sex))
  cfg <- config$twas %||% list()
  fdir <- config$tools$fusion_dir
  rlib <- config$tools$fusion_rlib
  for (p in c(file.path(fdir, "FUSION.assoc_test.R"), rlib)) if (!file.exists(p)) log_fatal(stage, paste(p, "missing (install script)"))
  out_dir <- stage_dir(config, sex, stage); dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  pos <- fusion_pos_files(config, out_dir)
  ldref <- config$refs$fusion_ldref
  ldsnps <- unique(unlist(lapply(1:22, function(c) fread(sprintf("%s%d.bim", ldref, c), header = FALSE, select = 2L)$V2)))
  env <- c(paste0("R_LIBS=", rlib), "OMP_NUM_THREADS=1", "OPENBLAS_NUM_THREADS=1")
  all <- list()
  for (f in postgwas_factors(config)) {
    fd <- file.path(out_dir, f); dir.create(fd, showWarnings = FALSE)
    ss <- file.path(fd, paste0(f, ".fusion.sumstats"))
    d <- load_postgwas(config, sex, f, "qfilter", cols = c("SNP", "A1", "A2", "Z"))
    fwrite(d[SNP %in% ldsnps], ss, sep = "\t")
    jobs <- CJ(i = seq_len(nrow(pos)), chr = 1:22)
    res <- parallel::mclapply(seq_len(nrow(jobs)), function(k) {
      p <- pos[jobs$i[k]]
      tissue <- sub("(\\.hg19)?\\.pos$", "", basename(p$pos))
      out <- file.path(fd, sprintf("%s.%d.dat", tissue, jobs$chr[k]))
      if (!file.exists(out)) {
        tryCatch(run_tool("Rscript", c("--vanilla", file.path(fdir, "FUSION.assoc_test.R"), "--sumstats", ss,
                                       "--weights", p$pos, "--weights_dir", p$wdir, "--ref_ld_chr", ldref,
                                       "--chr", jobs$chr[k], "--out", out),
                          stage, paste0(out, ".log"), env = env, wd = fdir), error = function(e) NULL)
      }
      if (file.exists(out)) fread(out)[, tissue := tissue] else NULL
    }, mc.cores = as.integer(cfg$workers %||% 48L), mc.preschedule = FALSE)
    x <- rbindlist(res, fill = TRUE)
    # every FUSION job failing is a broken installation, not "no TWAS signal": stop the stage
    if (!nrow(x)) {
      lg <- list.files(fd, "\\.dat\\.log$", full.names = TRUE)
      log_fatal(stage, sprintf("%s: no TWAS output from %d FUSION jobs; first job log:\n%s", f, nrow(jobs),
                               if (length(lg)) paste(utils::head(readLines(lg[1], warn = FALSE), 10), collapse = "\n") else "<none>"))
    }
    x <- x[!is.na(TWAS.P)]
    x[, `:=`(factor = f, n_tests = .N, sig_bonf = TWAS.P < 0.05 / .N)]
    if (all(startsWith(x$ID, "ENSG"))) {
      sym <- gencode_gene_symbols(config$refs$gencode_gtf)
      x[, ID := fifelse(is.na(sym$symbol[match(sub("\\..*$", "", ID), sym$gene_id)]), ID,
                        sym$symbol[match(sub("\\..*$", "", ID), sym$gene_id)])]
    }
    # joint / conditional post-processing where genes are significant
    sig <- unique(x[sig_bonf == TRUE, .(tissue, CHR)])
    for (k in seq_len(nrow(sig))) {
      top <- file.path(fd, sprintf("%s.%d.top", sig$tissue[k], sig$CHR[k]))
      raw <- fread(file.path(fd, sprintf("%s.%d.dat", sig$tissue[k], sig$CHR[k])))
      fwrite(raw[TWAS.P < 0.05 / x$n_tests[1]], top, sep = "\t")
      tryCatch(run_tool("Rscript", c("--vanilla", file.path(fdir, "FUSION.post_process.R"), "--sumstats", ss,
                                     "--input", top, "--out", sub("\\.top$", ".analysis", top),
                                     "--ref_ld_chr", ldref, "--chr", sig$CHR[k], "--save_loci", "--plot",
                                     "--locus_win", "100000"),
                        stage, paste0(top, ".post.log"), env = env, wd = fdir),
               error = function(e) log_warn(stage, sprintf("post_process %s: %s", basename(top), conditionMessage(e))))
    }
    all[[f]] <- x
    log_info(stage, sprintf("%s: %s gene-tissue models, Bonferroni %.2e -> %d significant (%d genes)",
                            f, format(nrow(x), big.mark = ","), 0.05 / nrow(x), sum(x$sig_bonf), uniqueN(x[sig_bonf == TRUE, ID])))
  }
  tw <- rbindlist(all, fill = TRUE)
  fwrite(tw, file.path(out_dir, "twas_all.tsv"), sep = "\t")
  joint <- list.files(out_dir, pattern = "joint_included\\.dat$", recursive = TRUE, full.names = TRUE)
  fwrite(rbindlist(lapply(joint, function(j) fread(j)[, source := basename(j)]), fill = TRUE),
         file.path(out_dir, "twas_joint.tsv"), sep = "\t")
  clear_postgwas_cache()
  write_stage_manifest(stage, sex, config, file.path(out_dir, c("twas_all.tsv", "twas_joint.tsv")))
  invisible(tw)
}
