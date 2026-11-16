# Conditional analysis stage (GWASTutorial "Conditional analysis"): GCTA-COJO stepwise selection
# of independent signals in each fine-mapping region, with a 20,000-sample UK Biobank hard-call
# panel as the LD reference (the guide asks for > 4,000 samples), as a cross-check on the number
# of SuSiE credible sets. .ma: SNP A1 A2 freq b se p N with freq = A1 frequency (Pan-UKB EUR),
# b/se on the std.lv scale and N = N_hat_std (mutually consistent, so Vp ~ 1).

run_cojo <- function(config, sex) {
  stage <- "cojo"
  log_info(stage, sprintf("=== GCTA-COJO: %s ===", sex))
  cfg <- config$cojo %||% list()
  gcta <- tool_path(config, "gcta")
  out_dir <- stage_dir(config, sex, stage)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  reg <- read_factor_registry(config, sex)
  ldp <- stage_dir(config, sex, "ld_panel")
  regs <- fread(file.path(ldp, "finemap_regions.tsv"))
  p_gws <- as.numeric(config$postgwas$p_gws %||% 5e-8)
  # each factor's full tier is read once here; the forked workers share it copy-on-write
  g_all <- lapply(setNames(nm = unique(regs$factor)), function(f) load_postgwas(config, sex, f, "full",
    cols = c("SNP", "CHR", "BP", "A1", "A2", "EAF_A1", "beta", "SE", "P", "LOG10P")))
  clear_postgwas_cache()

  one <- function(i) {
    x <- regs[i]
    pre <- file.path(out_dir, x$factor, x$ld_id)
    dir.create(dirname(pre), recursive = TRUE, showWarnings = FALSE)
    bed <- file.path(ldp, "cojo_bed", x$ld_id)
    g <- g_all[[x$factor]][CHR == x$chr & BP >= x$start & BP <= x$end & !is.na(EAF_A1)]
    bim <- fread(paste0(bed, ".bim"), header = FALSE, select = 2L, col.names = "SNP")
    g <- g[SNP %in% bim$SNP]
    if (!any(g$LOG10P > -log10(p_gws))) return(data.table(fm_id = x$fm_id, factor = x$factor, n_signals = 0L))
    ma <- g[, .(SNP, A1, A2, freq = EAF_A1, b = beta, se = SE, p = pmax(P, 1e-300),
                N = round(reg[[x$factor]]$N_hat_std))]
    fwrite(ma, paste0(pre, ".ma"), sep = " ")
    run_tool(gcta, c("--bfile", bed, "--cojo-file", paste0(pre, ".ma"), "--cojo-slct",
                     "--cojo-p", format(p_gws), "--cojo-wind", as.character(cfg$wind_kb %||% 10000),
                     "--cojo-collinear", as.character(cfg$collinear %||% 0.9),
                     "--diff-freq", as.character(cfg$diff_freq %||% 0.2), "--maf", "0.001",
                     "--thread-num", as.character(cfg$threads %||% 4), "--out", pre),
             stage, paste0(pre, ".gcta.log"))
    jma <- paste0(pre, ".jma.cojo")
    if (!file.exists(jma)) return(data.table(fm_id = x$fm_id, factor = x$factor, n_signals = 0L))
    j <- fread(jma)
    j[, `:=`(fm_id = x$fm_id, factor = x$factor)]
    fwrite(j, paste0(pre, ".signals.tsv"), sep = "\t")
    data.table(fm_id = x$fm_id, factor = x$factor, n_signals = nrow(j), signals = paste(j$SNP, collapse = ","))
  }
  res <- rbindlist(parallel::mclapply(seq_len(nrow(regs)), function(i) tryCatch(one(i), error = function(e)
    data.table(fm_id = regs$fm_id[i], factor = regs$factor[i], error = conditionMessage(e))),
    mc.cores = as.integer(cfg$workers %||% 8L), mc.preschedule = FALSE), fill = TRUE)
  diag_f <- file.path(stage_dir(config, sex, "finemap"), "region_diagnostics.tsv")
  diag <- if (file.exists(diag_f)) fread(diag_f) else NULL
  if ("n_cs" %in% names(diag)) res <- merge(res, diag[, .(fm_id, n_cs)], by = "fm_id", all.x = TRUE)
  fwrite(res, file.path(out_dir, "cojo_signals.tsv"), sep = "\t")
  if (all(c("n_cs", "n_signals") %in% names(res))) {
    ok <- res[!is.na(n_cs) & !is.na(n_signals)]
    log_info(stage, sprintf("COJO signals vs SuSiE credible sets: equal in %d/%d regions, within 1 in %d",
                            sum(ok$n_signals == ok$n_cs), nrow(ok), sum(abs(ok$n_signals - ok$n_cs) <= 1)))
  }
  clear_postgwas_cache()
  write_stage_manifest(stage, sex, config, file.path(out_dir, "cojo_signals.tsv"))
  invisible(res)
}
