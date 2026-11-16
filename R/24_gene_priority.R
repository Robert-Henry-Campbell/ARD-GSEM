# Gene prioritisation stage: candidate genes per credible set with transparent evidence columns
# and a weighted evidence score (config gene_priority$weights). Candidates: protein-coding genes
# (GENCODE v19) with TSS within window_kb of the credible set's PIP lead, plus genes whose body
# overlaps any credible-set member. Evidence sources are optional; a missing stage output simply
# contributes no evidence (logged).

run_gene_priority <- function(config, sex) {
  stage <- "gene_priority"
  log_info(stage, sprintf("=== Gene prioritisation: %s ===", sex))
  cfg <- config$gene_priority %||% list()
  wts <- modifyList(list(coding_pip = 2, qtl_coloc = 2, twas = 1, magma = 1, nearest_tss = 1, cs_overlap = 1,
                         qtl_lookup = 0.5), cfg$weights %||% list())
  out_dir <- stage_dir(config, sex, stage); dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  rd <- function(st, file) { p <- file.path(stage_dir(config, sex, st), file); if (file.exists(p)) fread(p) else NULL }
  cs <- rd("finemap", "all_cs.tsv")
  if (is.null(cs) || !nrow(cs)) log_fatal(stage, "finemap/all_cs.tsv missing or empty")
  genes <- read_gencode_genes(config$refs$gencode_gtf)
  win <- as.integer(cfg$window_kb %||% 500L) * 1000L

  cand <- genes[cs[, .(cs_id, factor, CHR = lead_chr, lo = lead_bp - win, hi = lead_bp + win, lead_bp)],
                on = .(CHR, TSS >= lo, TSS <= hi), nomatch = NULL,
                .(cs_id, factor, symbol, gene_id, CHR, TSS = x.TSS, dist_tss = abs(x.TSS - lead_bp))]
  cand[, nearest_tss := dist_tss == min(dist_tss), by = cs_id]

  members <- cs[, .(SNP = unlist(strsplit(members, ","))), by = .(cs_id, factor)]
  ann <- rd("annotate", "variants_annotated.tsv.gz")
  if (!is.null(ann)) {
    members <- merge(members, ann, by = "SNP", all.x = TRUE)
    ov <- genes[members[!is.na(BP), .(cs_id, CHR, BP)], on = .(CHR, start <= BP, end >= BP), nomatch = NULL,
                .(cs_id, symbol)][, .(cs_overlap = TRUE), by = .(cs_id, symbol)]
    cand <- merge(cand, unique(ov), by = c("cs_id", "symbol"), all.x = TRUE)
  } else log_warn(stage, "annotate outputs missing: no overlap/coding evidence")
  cand[is.na(cs_overlap), cs_overlap := FALSE]

  # coding / splice members (PIP >= 0.1): VEP gene when available, else Neale consequence + nearest gene
  pipf <- rbindlist(lapply(seq_len(nrow(cs)), function(i) {
    p <- file.path(stage_dir(config, sex, "finemap"), cs$factor[i], sub(paste0("^", cs$factor[i], "_"), "", cs$fm_id[i]), "pip.tsv.gz")
    k <- cs$cs[i]; id <- cs$cs_id[i]   # pip.tsv.gz has its own `cs` column
    if (file.exists(p)) fread(p)[cs == k, .(cs_id = id, SNP, pip)] else NULL
  }))
  cod <- data.table()
  if (!is.null(ann) && nrow(pipf)) {
    x <- merge(pipf[pip >= 0.1], ann, by = "SNP")
    coding_rx <- "missense|stop_gained|stop_lost|start_lost|frameshift|splice_acceptor|splice_donor|inframe"
    if ("vep_Consequence" %in% names(x)) {
      cod <- x[grepl(coding_rx, vep_Consequence), .(coding_pip = max(pip)), by = .(cs_id, symbol = vep_SYMBOL)]
    } else {
      cod <- x[grepl(coding_rx, neale_consequence), .(coding_pip = max(pip)), by = .(cs_id, symbol = nearest_tss_gene)]
    }
  }
  if (nrow(cod)) cand <- merge(cand, cod, by = c("cs_id", "symbol"), all.x = TRUE) else cand[, coding_pip := NA_real_]

  mg <- rbindlist(lapply(postgwas_factors(config), function(f) {
    p <- file.path(stage_dir(config, sex, "magma"), paste0(f, "_genes.tsv"))
    if (file.exists(p)) fread(p)[, .(factor = f, symbol = SYMBOL, magma_p = P, magma_sig = sig_bonf)] else NULL
  }))
  if (nrow(mg)) cand <- merge(cand, mg, by = c("factor", "symbol"), all.x = TRUE) else cand[, `:=`(magma_p = NA_real_, magma_sig = NA)]
  tw <- rd("twas", "twas_all.tsv")
  if (!is.null(tw)) {
    t2 <- tw[, .(twas_best_p = min(TWAS.P, na.rm = TRUE), twas_sig = any(sig_bonf), twas_tissues = paste(unique(tissue[sig_bonf]), collapse = ",")),
             by = .(factor, symbol = ID)]
    cand <- merge(cand, t2, by = c("factor", "symbol"), all.x = TRUE)
  } else cand[, `:=`(twas_best_p = NA_real_, twas_sig = NA, twas_tissues = NA_character_)]
  lk <- rd("qtl_lookup", "lookup_all.tsv")
  if (!is.null(lk)) {
    l2 <- lk[sig_bonf == TRUE, .(qtl_lookup_sources = paste(unique(paste0(family, ":", dataset_id)), collapse = ";")),
             by = .(cs_id, symbol = gene_symbol)]
    cand <- merge(cand, l2, by = c("cs_id", "symbol"), all.x = TRUE)
  } else cand[, qtl_lookup_sources := NA_character_]
  cl <- rd("coloc", "coloc_results.tsv")
  if (!is.null(cl) && nrow(cl)) {
    c2 <- cl[, .(eqtl_coloc_h4 = suppressWarnings(max(PP.H4.abf[family == "eQTL"], na.rm = TRUE)),
                 pqtl_coloc_h4 = suppressWarnings(max(PP.H4.abf[family == "pQTL"], na.rm = TRUE))),
             by = .(cs_id, symbol = gene_symbol)]
    for (col in c("eqtl_coloc_h4", "pqtl_coloc_h4")) set(c2, which(!is.finite(c2[[col]])), col, NA_real_)
    cand <- merge(cand, c2, by = c("cs_id", "symbol"), all.x = TRUE)
  } else cand[, `:=`(eqtl_coloc_h4 = NA_real_, pqtl_coloc_h4 = NA_real_)]

  t <- function(x) !is.na(x) & x
  cand[, score := wts$nearest_tss * nearest_tss + wts$cs_overlap * cs_overlap +
         wts$coding_pip * fifelse(is.na(coding_pip), 0, coding_pip) +
         wts$magma * t(magma_sig) + wts$twas * t(twas_sig) + wts$qtl_lookup * (!is.na(qtl_lookup_sources)) +
         wts$qtl_coloc * (t(eqtl_coloc_h4 >= 0.8) | t(pqtl_coloc_h4 >= 0.8))]
  cand[, evidence := paste(c(if (nearest_tss) "nearest-TSS", if (cs_overlap) "CS-overlap",
                             if (!is.na(coding_pip)) sprintf("coding(PIP %.2f)", coding_pip),
                             if (t(magma_sig)) "MAGMA", if (t(twas_sig)) "TWAS",
                             if (t(eqtl_coloc_h4 >= 0.8)) "eQTL-coloc", if (t(pqtl_coloc_h4 >= 0.8)) "pQTL-coloc",
                             if (!is.na(qtl_lookup_sources)) "QTL-lookup"), collapse = ";"), by = seq_len(nrow(cand))]
  cand <- merge(cand, cs[, .(cs_id, lead_snp, lead_pip, heterogeneous, factor_level)], by = "cs_id")
  setorder(cand, factor, cs_id, -score, dist_tss)
  cand[, rank := seq_len(.N), by = cs_id]
  fwrite(cand, file.path(out_dir, "gene_priority.tsv"), sep = "\t")
  top <- cand[rank == 1L]
  fwrite(top, file.path(out_dir, "cs_top_genes.tsv"), sep = "\t")
  log_info(stage, sprintf("%d candidate genes for %d credible sets; top gene has >1 line of evidence in %d sets",
                          nrow(cand), uniqueN(cand$cs_id), sum(lengths(strsplit(top$evidence, ";")) > 1)))
  write_stage_manifest(stage, sex, config, file.path(out_dir, c("gene_priority.tsv", "cs_top_genes.tsv")))
  invisible(cand)
}
