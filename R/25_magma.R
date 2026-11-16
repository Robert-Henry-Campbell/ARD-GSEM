# MAGMA stage (GWASTutorial "Gene/Gene-set tests by MAGMA"), on the qfilter tier:
#  - SNP -> gene annotation (NCBI37.3, window 35 kb upstream / 10 kb downstream, as in the guide)
#  - gene analysis on SNP p-values (g1000_eur LD, SNP-wise mean model, N = N_hat_std per SNP)
#  - gene-set analysis (MSigDB 2024.1 Entrez GMT)
#  - gene-property (tissue) analysis: GTEx v8 median TPM -> log2(TPM + 1) per tissue, conditioned
#    on the average across tissues (FUMA convention), one-sided (higher expression)
# Bonferroni across genes and tissues; BH and Bonferroni across gene sets. Cell-type gene
# properties (Tabula Sapiens) run in the scdrs stage on this stage's gene results.

magma_read_out <- function(f) {
  x <- readLines(f)
  x <- x[!startsWith(x, "#")]
  fread(text = x, header = TRUE)
}

gtex_gene_covar <- function(config, out_file) {
  tpm <- fread(cmd = paste("zcat", shQuote(config$refs$gtex_median_tpm)), skip = 2)
  loc <- fread(config$refs$magma_gene_loc, header = FALSE, col.names = c("entrez", "chr", "start", "end", "strand", "symbol"))
  tpm <- tpm[Description %in% loc$symbol][!duplicated(Description)]
  m <- as.matrix(tpm[, !c("Name", "Description")])
  m <- log2(m + 1)
  cv <- data.table(GENE = loc$entrez[match(tpm$Description, loc$symbol)], m)
  setnames(cv, make.names(names(cv)))
  cv[, Average := rowMeans(m)]
  fwrite(cv, out_file, sep = " ")
  out_file
}

run_magma <- function(config, sex) {
  stage <- "magma"
  log_info(stage, sprintf("=== MAGMA: %s ===", sex))
  cfg <- config$magma %||% list()
  magma <- tool_path(config, "magma")
  out_dir <- stage_dir(config, sex, stage); dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  reg <- read_factor_registry(config, sex)
  win <- cfg$window %||% "35,10"
  covar <- gtex_gene_covar(config, file.path(out_dir, "gtex_v8_covar.txt"))
  summ <- list()
  for (f in postgwas_factors(config)) {
    pre <- file.path(out_dir, f)
    d <- load_postgwas(config, sex, f, "qfilter", cols = c("SNP", "CHR", "BP", "P", "LOG10P"))
    fwrite(d[, .(SNP, CHR, BP)], paste0(pre, ".snploc"), sep = "\t", col.names = FALSE)
    fwrite(d[, .(SNP, P = pmax(P, 10^-pmin(LOG10P, 300)), N = round(reg[[f]]$N_hat_std))], paste0(pre, ".pval"), sep = "\t")
    run_tool(magma, c("--annotate", paste0("window=", win), "--snp-loc", paste0(pre, ".snploc"),
                      "--gene-loc", config$refs$magma_gene_loc, "--out", pre),
             stage, paste0(pre, ".annot.log"), expect = paste0(pre, ".genes.annot"))
    run_tool(magma, c("--bfile", config$refs$magma_bfile, "--pval", paste0(pre, ".pval"), "ncol=N",
                      "--gene-annot", paste0(pre, ".genes.annot"), "--gene-model", "snp-wise=mean", "--out", pre),
             stage, paste0(pre, ".gene.log"), expect = paste0(pre, ".genes.out"))
    run_tool(magma, c("--gene-results", paste0(pre, ".genes.raw"), "--set-annot", config$refs$msigdb_gmt,
                      "--out", paste0(pre, "_msigdb")),
             stage, paste0(pre, "_msigdb.stdout.log"), expect = paste0(pre, "_msigdb.gsa.out"))
    run_tool(magma, c("--gene-results", paste0(pre, ".genes.raw"), "--gene-covar", covar,
                      "--model", "condition-hide=Average", "direction=greater", "--out", paste0(pre, "_gtex")),
             stage, paste0(pre, "_gtex.stdout.log"), expect = paste0(pre, "_gtex.gsa.out"))

    loc <- fread(config$refs$magma_gene_loc, header = FALSE, col.names = c("GENE", "chr", "start", "end", "strand", "SYMBOL"))
    g <- merge(magma_read_out(paste0(pre, ".genes.out")), loc[, .(GENE, SYMBOL)], by = "GENE", all.x = TRUE)
    g[, `:=`(factor = f, p_bonf = pmin(1, P * .N), sig_bonf = P < 0.05 / .N)]
    s <- magma_read_out(paste0(pre, "_msigdb.gsa.out"))
    s[, `:=`(factor = f, q_bh = p.adjust(P, "BH"), sig_bonf = P < 0.05 / .N)]
    t <- magma_read_out(paste0(pre, "_gtex.gsa.out"))
    t[, `:=`(factor = f, sig_bonf = P < 0.05 / .N)]
    fwrite(g[order(P)], paste0(pre, "_genes.tsv"), sep = "\t")
    fwrite(s[order(P)], paste0(pre, "_genesets.tsv"), sep = "\t")
    fwrite(t[order(P)], paste0(pre, "_tissues.tsv"), sep = "\t")
    summ[[f]] <- data.table(factor = f, n_genes = nrow(g), n_genes_bonf = sum(g$sig_bonf),
                            n_sets_fdr = sum(s$q_bh < 0.05), n_sets_bonf = sum(s$sig_bonf),
                            n_tissues_bonf = sum(t$sig_bonf), top_tissue = t[which.min(P), VARIABLE])
    log_info(stage, sprintf("%s: %d Bonferroni genes; %d gene sets FDR<5%%; %d tissues Bonferroni (top: %s)",
                            f, sum(g$sig_bonf), sum(s$q_bh < 0.05), sum(t$sig_bonf), summ[[f]]$top_tissue))
  }
  fwrite(rbindlist(summ), file.path(out_dir, "magma_summary.csv"))
  clear_postgwas_cache()
  write_stage_manifest(stage, sex, config, list.files(out_dir, pattern = "\\.(tsv|csv)$", full.names = TRUE))
  invisible(rbindlist(summ))
}
