for (f in c("R/14_split_vcf.R", "R/15_postgwas_common.R", "R/18_annotate.R", "R/19_ld_panel.R",
            "R/20_finemap.R", "R/22_qtl_lookup.R", "R/23_coloc.R", "R/24_gene_priority.R", "R/25_magma.R",
            "R/26_sldsc.R", "R/28_scdrs.R")) {
  source(file.path(project_root, f), local = TRUE)
}

test_that("baselineLD annotations are classed binary or continuous from their values", {
  pre <- file.path(tempfile(), "baselineLD."); dir.create(dirname(pre))
  fwrite(data.table(CHR = 22L, BP = 1:4, SNP = paste0("rs", 1:4), CM = 0, base = 1, Coding = c(0, 1, 0, 0),
                    GERP.NS = c(-1.2, 0, 3.4, 0.5), MaxCPP = c(0, 0.2, 1, 0)),
         paste0(pre, "22.annot.gz"), sep = "\t")
  expect_identical(baseline_binary(pre), c(base = TRUE, Coding = TRUE, GERP.NS = FALSE, MaxCPP = FALSE))
})

test_that("split_vcf records the per-factor VCFs written by the vcf stage", {
  tmp <- tempfile()
  cfg <- list(paths = list(output_dir = tmp), postgwas = list(factors = c("F1", "F2")), vcf = list(basename = "t"))
  qd <- file.path(tmp, "bothsex", "qc_filters"); dir.create(qd, recursive = TRUE)
  writeLines("{}", file.path(qd, "factor_registry.json"))
  vd <- file.path(tmp, "bothsex", "vcf"); dir.create(vd)
  file.create(file.path(vd, c("t_F1.vcf.gz", "t_F1.vcf.gz.tbi", "t_F2.vcf.gz")))
  expect_error(run_split_vcf(cfg, "bothsex"), "t_F2.vcf.gz")
  file.create(file.path(vd, "t_F2.vcf.gz.tbi"))
  expect_identical(basename(run_split_vcf(cfg, "bothsex")), c("t_F1.vcf.gz", "t_F2.vcf.gz"))
})

test_that("UKB LD samples: QC filters applied, COJO subset nested, keep files private", {
  tmp <- tempfile(); dir.create(tmp)
  pg <- file.path(tmp, "pgen"); dir.create(pg)
  fwrite(data.table(`#FID` = 1001:1030, IID = 1001:1030, SEX = 1L), file.path(pg, "UKB_chr22.psam"), sep = "\t")
  qc <- data.table(eid = 1001:1032, `22006-0.0` = 1L, `22021-0.0` = 0L, `22027-0.0` = NA_integer_,
                   `22019-0.0` = NA_integer_)
  qc[eid == 1001, `22006-0.0` := NA_integer_]        # not White British
  qc[eid == 1002, `22021-0.0` := 1L]                 # has relatives
  qc[eid == 1003, `22027-0.0` := 1L]                 # het/missingness outlier
  qc[eid == 1004, `22019-0.0` := 1L]                 # sex-chromosome aneuploidy
  fwrite(qc, file.path(tmp, "qc.tsv"), sep = "\t")
  cfg <- list(refs = list(ukb_pgen_dir = pg, ukb_sample_qc = file.path(tmp, "qc.tsv"), ukb_withdrawals = ""))
  out <- file.path(tmp, "out"); dir.create(out)
  steps <- select_ld_samples(cfg, out, n_ld = 10L, n_cojo = 4L, seed = 7L)
  expect_identical(c(steps$in_qc_table, steps$no_sex_aneuploidy), c(30L, 26L))
  ld <- fread(file.path(out, "keep_ld.txt"), header = FALSE)
  cojo <- fread(file.path(out, "keep_cojo.txt"), header = FALSE)
  expect_identical(c(nrow(ld), nrow(cojo)), c(10L, 4L))
  expect_false(any(ld$V2 %in% 1001:1004))
  expect_true(all(cojo$V2 %in% ld$V2))
  expect_identical(as.character(file.mode(file.path(out, "keep_ld.txt"))), "600")
  out2 <- file.path(tmp, "out2"); dir.create(out2)
  select_ld_samples(cfg, out2, n_ld = 10L, n_cojo = 4L, seed = 7L)
  expect_identical(readLines(file.path(out2, "keep_ld.txt")), readLines(file.path(out, "keep_ld.txt")))
})

test_that("fine-mapping regions wider than max_mb are split at midpoints between hits", {
  tmp <- tempfile()
  cfg <- list(paths = list(output_dir = tmp))
  ld <- file.path(tmp, "bothsex", "loci"); dir.create(ld, recursive = TRUE)
  fwrite(data.table(region_id = c("F_chr1_a", "F_chr2_b", "F_chr6_mhc"), factor = "F", chr = c(1L, 2L, 6L),
                    start = c(0L, 1000000L, 29000000L), end = c(5000000L, 2000000L, 30000000L),
                    finemap_eligible = c(TRUE, TRUE, FALSE)), file.path(ld, "regions.tsv"), sep = "\t")
  fwrite(data.table(CHR = c(1L, 1L, 1L, 2L), BP = c(1000000L, 2000000L, 4000000L, 1500000L)),
         file.path(ld, "F_independent_hits.tsv"), sep = "\t")
  r <- finemap_regions(cfg, "bothsex", max_mb = 3)
  expect_identical(r$locus_id, c("F_chr1_a", "F_chr1_a", "F_chr2_b"))
  expect_identical(r$start, c(0L, 3000000L, 1000000L))
  expect_identical(r$end, c(3000000L, 5000000L, 2000000L))
  expect_true(all(r$end - r$start <= 3e6))
  expect_identical(r$fm_id[1], "F_chr1_0_3000000")
})

test_that("QTL alleles align to A1; the lead is preferred; Bonferroni counts tests per family", {
  expect_identical(allele_sign(c("A", "A", "A"), c("G", "G", "G"), c("A", "G", "C"), c("G", "A", "T")),
                   c(1, -1, NA_real_))
  x <- data.table(cs_id = "cs1", dataset_id = "d", molecular_trait_id = "g", SNP = c("rs1", "rs2", "rs3"),
                  is_lead = c(FALSE, FALSE, TRUE), pip = c(0.2, 0.5, 0.3))
  expect_identical(pick_tested_snp(x)$SNP, "rs3")
  expect_identical(pick_tested_snp(x[is_lead == FALSE])$SNP, "rs2")
  y <- data.table(family = c(rep("eQTL", 4), "pQTL"), p = c(0.01, 0.0125, 0.02, 1e-6, 0.04),
                  z_qtl = c(2, -2, 1, 5, 2), z_factor = c(1, 1, 1, 1, -1))
  y <- qtl_multiple_testing(y)
  expect_identical(y$n_tests_family, c(4L, 4L, 4L, 4L, 1L))
  expect_identical(y$sig_bonf, c(TRUE, TRUE, FALSE, TRUE, TRUE))
  expect_identical(y$sig_fdr, c(TRUE, TRUE, TRUE, TRUE, TRUE))
  expect_identical(y$concordant, c(TRUE, FALSE, TRUE, TRUE, FALSE))
})

test_that("eQTL Catalogue reads: tabix CLI records, de-duplication, and failures that raise", {
  tabix <- "/mnt/sdg/robert/miniforge3/envs/gsem-cli/bin/tabix"; bgzip <- sub("tabix$", "bgzip", tabix)
  skip_if_not(file.exists(tabix) && file.exists(bgzip))
  td <- tempfile(); dir.create(td)
  mk <- function(trait, chr, pos) data.table(molecular_trait_id = trait, chromosome = chr, position = pos, ref = "A", alt = "G",
                                             variant = sprintf("chr%s_%d_A_G", chr, pos), ma_samples = 10L, maf = 0.2, pvalue = 1e-3,
                                             beta = 0.1, se = 0.02, type = "SNP", ac = 40L, an = 400L, r2 = 1, molecular_trait_object_id = trait,
                                             gene_id = trait, median_tpm = 5, rsid = sprintf("rs%d", pos))
  d <- rbind(mk("G1", "1", c(100L, 200L, 300L)), mk("G2", "1", c(100L, 250L)), mk("G1", "2", 500L))[order(chromosome, position)]
  raw <- file.path(td, "T1.all.tsv"); fwrite(d, raw, sep = "\t")
  system2(bgzip, c("-f", raw)); system2(tabix, c("-f", "-s2", "-b3", "-e3", "-S1", paste0(raw, ".gz")))
  src <- paste0(raw, ".gz"); cache <- file.path(td, "tbi"); dir.create(cache)
  file.copy(paste0(src, ".tbi"), file.path(cache, "T1.all.tsv.gz.tbi"))
  idx <- file.path(cache, "T1.all.tsv.gz.tbi")
  # exact positions; a duplicated / overlapping region does not duplicate records
  l <- tabix_read(tabix, src, idx, data.table(chr = c("1", "1", "1", "2"), start = c(100L, 100L, 250L, 500L), end = c(100L, 100L, 300L, 500L)))
  x <- parse_eqtlcat_lines(l)
  expect_equal(nrow(x), 5L)                                    # G1/G2 at 100, G2 at 250, G1 at 300, G1 at chr2:500
  expect_identical(sort(unique(x$position)), c(100L, 250L, 300L, 500L))
  expect_identical(tabix_seqnames(tabix, src, idx), c("1", "2"))
  # dataset query: GRCh38 positions, chromosome naming, query count
  ds <- data.table(dataset_id = "T1", ftp_path = src)
  q <- query_eqtlcat_dataset(ds, data.table(CHR = c(1L, 1L, 2L, 3L), BP38 = c(200L, 999L, 500L, 7L)), cache, tabix)
  expect_identical(sort(q$position), c(200L, 500L)); expect_identical(attr(q, "n_query"), 3L)   # chr3 absent from the file
  # a damaged file must raise, not return an empty result
  bad <- file.path(td, "T2.all.tsv.gz"); writeBin(readBin(src, "raw", file.size(src))[seq_len(file.size(src) %/% 2)], bad)
  expect_error(tabix_read(tabix, bad, idx, data.table(chr = "1", start = 100L, end = 300L), retries = 2L, wait = 0), "failed after 2 attempts")
})

test_that("eQTL Catalogue lbf extracts: streamed once per dataset, one matrix per needed trait", {
  skip_if_not(nzchar(Sys.which("curl")))
  td <- tempfile(); dir.create(td)
  lbf <- data.table(molecular_trait_id = rep(c("G1", "G2", "G3"), each = 3), region = "r", variant = rep(c("v1", "v2", "v3"), 3),
                    chromosome = "1", position = rep(1:3, 3), lbf_variable1 = c(5, 0, 0, 0, 7, 0, 1, 1, 1), lbf_variable2 = 0)
  f <- file.path(td, "T1.lbf_variable.txt"); fwrite(lbf, f, sep = "\t"); system2("gzip", c("-f", f))
  ds <- data.table(dataset_id = "T1", ftp_lbf_path = paste0("file://", f, ".gz"))
  m <- eqtlcat_lbf(ds, "G2", file.path(td, "cache"), traits = c("G2", "G1"))
  expect_equal(dim(m), c(2L, 3L)); expect_identical(colnames(m), c("v1", "v2", "v3")); expect_equal(m[1, "v2"], 7)
  expect_null(eqtlcat_lbf(ds, "G9", file.path(td, "cache"), traits = c("G1", "G2")))         # not extracted
  expect_false(file.exists(file.path(td, "cache", "T1", paste0(digest::digest("G3"), ".rds"))))
  ds_bad <- data.table(dataset_id = "T9", ftp_lbf_path = paste0("file://", td, "/missing.gz"))
  expect_error(eqtlcat_lbf_subset(ds_bad, "G1", file.path(td, "cache")), "lbf stream for T9 failed")
})

test_that("SuSiE-coloc aligns log Bayes factors by SNP name, not column position", {
  skip_if_not_installed("coloc")
  snps <- paste0("rs", 1:60)
  lbf1 <- matrix(c(rep(0, 9), 30, rep(0, 50)), 1, dimnames = list(NULL, snps))         # signal at rs10
  lbf2 <- matrix(c(rep(0, 9), 200, rep(0, 50)), 1, dimnames = list(NULL, snps))
  perm <- c(60:1)                                                                       # same SNPs, other order
  pri <- list(p1 = 1e-4, p2 = 1e-4, p12 = 5e-6)
  r <- bfbf_pair(lbf1, lbf2[, perm, drop = FALSE], pri)
  expect_gt(r$PP.H4.abf, 0.99); expect_identical(r$top_snp, "rs10")
  # the trap itself: coloc pairs equal-sized matrices by position
  raw <- coloc::coloc.bf_bf(lbf1, lbf2[, perm, drop = FALSE], p1 = 1e-4, p2 = 1e-4, p12 = 5e-6)$summary
  expect_lt(raw$PP.H4.abf, 0.01)
  # extra, unshared SNPs in either matrix
  lbf3 <- cbind(lbf2[, perm, drop = FALSE], matrix(0, 1, 5, dimnames = list(NULL, paste0("x", 1:5))))
  expect_gt(bfbf_pair(lbf1, lbf3, pri)$PP.H4.abf, 0.99)
  expect_gt(bfbf_pair(lbf3, lbf1, pri)$PP.H4.abf, 0.99)
  # different signals stay distinct
  lbf4 <- matrix(c(rep(0, 39), 200, rep(0, 20)), 1, dimnames = list(NULL, snps))
  expect_gt(bfbf_pair(lbf1, lbf4[, perm, drop = FALSE], pri)$PP.H3.abf, 0.99)
})

test_that("eQTLGen Z converts to per-allele beta/SE on the A1 scale (SMR formula)", {
  tmp <- tempfile(); dir.create(tmp)
  f <- file.path(tmp, "eqtlgen.txt.gz")
  e <- data.table(Pvalue = 1e-5, SNP = c("rs1", "rs2", "rs3", "rs4"), SNPChr = 1L, SNPPos = c(100L, 200L, 300L, 250L),
                  AssessedAllele = c("A", "G", "C", "A"), OtherAllele = c("G", "A", "T", "G"),
                  Zscore = c(4, 4, 3, 9), Gene = c("ENSG1", "ENSG1", "ENSG1", "ENSG2"), GeneSymbol = "X",
                  GeneChr = 1L, GenePos = 1L, NrCohorts = 1L, NrSamples = 30000L, FDR = 0, BonferroniP = 1)
  fwrite(e, f, sep = "\t")
  gw <- data.table(SNP = c("rs1", "rs2", "rs3", "rs4"), CHR = 1L, BP = c(100L, 200L, 300L, 250L),
                   A1 = c("A", "A", "A", "A"), A2 = c("G", "G", "G", "G"))
  af <- data.table(SNP = c("rs1", "rs2", "rs3", "rs4"), AlleleB = c("G", "G", "T", "G"), AlleleB_all = c(0.7, 0.7, 0.5, 0.5))
  r <- eqtlgen_region(list(refs = list(eqtlgen_cis = f)), "ENSG1", gw, af)
  expect_identical(r$SNP, c("rs1", "rs2"))                        # rs3 alleles mismatch, rs4 other gene
  den <- sqrt(2 * 0.3 * 0.7 * (30000 + 16))
  expect_equal(r$beta, c(4 / den, -4 / den))
  # the one-pass extract (used by run_coloc) gives the same rows; its cache is keyed to the gene list
  eg <- eqtlgen_subset(list(refs = list(eqtlgen_cis = f)), c("ENSG2", "ENSG1"), tmp)
  expect_identical(sort(unique(eg$Gene)), c("ENSG1", "ENSG2"))
  r2 <- eqtlgen_region(list(refs = list(eqtlgen_cis = f)), "ENSG1", gw, af, subset = eg)
  expect_equal(as.data.frame(r2), as.data.frame(r))
  expect_null(eqtlgen_region(list(refs = list(eqtlgen_cis = f)), "ENSG9", gw, af, subset = eg))
  eg2 <- eqtlgen_subset(list(refs = list(eqtlgen_cis = f)), "ENSG2", tmp)   # new gene list -> rebuilt
  expect_identical(unique(eg2$Gene), "ENSG2")
  expect_equal(r$se, c(1 / den, 1 / den))
  expect_equal(r$maf, c(0.3, 0.3))
})

test_that("MAGMA output tables are read past the comment header", {
  f <- tempfile(fileext = ".gsa.out")
  writeLines(c("# TOTAL_GENES = 2", "# CONDITIONED_INTERNAL = Average",
               "VARIABLE  TYPE  NGENES  BETA  BETA_STD  SE  P  FULL_NAME",
               "Liver  COVAR  18000  0.05  0.02  0.01  1e-05  Liver",
               "Lung  COVAR  18000  0.01  0.004  0.01  0.3  Lung"), f)
  x <- magma_read_out(f)
  expect_identical(x$VARIABLE, c("Liver", "Lung"))
  expect_equal(x$P, c(1e-5, 0.3))
})

test_that("scDRS scores are reused only for the same call and gene sets, with every factor scored", {
  tmp <- tempfile(); dir.create(tmp)
  inputs <- file.path(tmp, c("ts.h5ad", "ts_cov.tsv")); file.create(inputs)
  Sys.setFileTime(inputs, Sys.time() - 3600)
  out <- file.path(tmp, "scdrs"); dir.create(out)
  file.create(file.path(out, paste0(rep(c("F1", "F2", "F3"), each = 2), c(".score.gz", ".full_score.gz"))))
  log <- file.path(tmp, "compute_score.log")
  args <- c("compute-score", "--h5ad-file", inputs[1], "--n-ctrl", "1000", "--out-folder", out)
  # the CLI echoes every option (defaults too) under "Call:"; stdout can land after a progress bar's \r
  writeLines(c("Call: scdrs compute-score \\", paste("--h5ad-file", inputs[1], "\\"), "--weight-opt vs \\",
               "--n-ctrl 1000 \\", paste("--out-folder", out), "",
               "Trait=F1, n_gene=900: 5/100 FDR<0.1 cells, 9/100 FDR<0.2 cells (sys_time=10.0s)",
               "\rComputing control scores: 100%|#| 1000/1000Trait=F2, n_gene=880: 0/100 FDR<0.1 cells"), log)
  gs <- c("TRAIT\tGENESET", "F1\tG1:2.5,G2:1.9", "F2\tG3:3.1")
  expect_true(scdrs_scores_current(args, c("F1", "F2"), log, inputs, gs, gs))
  expect_false(scdrs_scores_current(args, c("F1", "F2"), log, inputs, NULL, gs))                  # first run
  expect_false(scdrs_scores_current(args, c("F1", "F2"), log, inputs, gs, c(gs, "F3\tG4:2.2")))   # gene sets changed
  expect_false(scdrs_scores_current(replace(args, 5, "500"), c("F1", "F2"), log, inputs, gs, gs)) # other options
  expect_false(scdrs_scores_current(args, c("F1", "F2", "F3"), log, inputs, gs, gs))              # F3 never finished
  expect_false(scdrs_scores_current(args, c("F1", "F4"), log, inputs, gs, gs))                    # no F4 score files
  Sys.setFileTime(inputs, Sys.time() + 60)
  expect_false(scdrs_scores_current(args, c("F1", "F2"), log, inputs, gs, gs))                    # inputs rebuilt since
})

test_that("gene priority scores and labels each line of evidence", {
  tmp <- tempfile()
  cfg <- list(paths = list(output_dir = tmp), postgwas = list(factors = "F1"), refs = list(gencode_gtf = "unused"))
  sd <- function(st) { d <- file.path(tmp, "bothsex", st); dir.create(d, recursive = TRUE, showWarnings = FALSE); d }
  fwrite(data.table(cs_id = "F1_chr1_1000_5000_cs1", factor = "F1", cs = 1L, fm_id = "F1_chr1_1000_5000",
                    lead_chr = 1L, lead_bp = 2000L, lead_snp = "rs1", lead_pip = 0.9, members = "rs1,rs2",
                    heterogeneous = FALSE, factor_level = TRUE), file.path(sd("finemap"), "all_cs.tsv"), sep = "\t")
  dir.create(file.path(tmp, "bothsex", "finemap", "F1", "chr1_1000_5000"), recursive = TRUE)
  fwrite(data.table(SNP = c("rs1", "rs2"), pip = c(0.9, 0.08), cs = 1L),
         file.path(tmp, "bothsex", "finemap", "F1", "chr1_1000_5000", "pip.tsv.gz"), sep = "\t")
  fwrite(data.table(SNP = c("rs1", "rs2"), CHR = 1L, BP = c(2000L, 3000L),
                    neale_consequence = c("missense_variant", "intron_variant"), nearest_tss_gene = c("GENEC", "GENEB")),
         file.path(sd("annotate"), "variants_annotated.tsv.gz"), sep = "\t")
  fwrite(data.table(SYMBOL = c("GENEA", "GENEB"), P = c(1e-8, 0.5), sig_bonf = c(TRUE, FALSE)),
         file.path(sd("magma"), "F1_genes.tsv"), sep = "\t")
  fwrite(data.table(cs_id = "F1_chr1_1000_5000_cs1", gene_symbol = "GENEB", family = "eQTL", dataset_id = "QTD1",
                    sig_bonf = TRUE), file.path(sd("qtl_lookup"), "lookup_all.tsv"), sep = "\t")
  fwrite(data.table(cs_id = "F1_chr1_1000_5000_cs1", gene_symbol = "GENEB", family = "eQTL", PP.H4.abf = 0.92),
         file.path(sd("coloc"), "coloc_results.tsv"), sep = "\t")
  read_gencode_genes <- function(gtf) data.table(gene_id = c("E1", "E2", "E3"), symbol = c("GENEA", "GENEB", "GENEC"),
                                                 CHR = 1L, start = c(1500L, 2900L, 1900L), end = c(2500L, 4000L, 2100L),
                                                 strand = c("+", "-", "+"), TSS = c(1500L, 4000L, 1900L))
  environment(run_gene_priority) <- environment()
  g <- run_gene_priority(cfg, "bothsex")
  expect_identical(g$symbol, c("GENEC", "GENEB", "GENEA"))
  expect_equal(g$score, c(1 + 1 + 2 * 0.9, 1 + 0.5 + 2, 1 + 1))
  expect_identical(g$evidence, c("nearest-TSS;CS-overlap;coding(PIP 0.90)", "CS-overlap;eQTL-coloc;QTL-lookup",
                                 "CS-overlap;MAGMA"))
  expect_true(file.exists(file.path(tmp, "bothsex", "gene_priority", "cs_top_genes.tsv")))
})

test_that("SuSiE-RSS wrapper recovers a simulated causal variant", {
  skip_if_not_installed("susieR")
  p <- 60
  R <- 0.7^abs(outer(seq_len(p), seq_len(p), "-"))
  z <- setNames(as.numeric(8 * R[, 25]), paste0("rs", seq_len(p)))
  fit <- finemap_one(z, R, n = 20000, list(L = 5L))
  cs <- susieR::susie_get_cs(fit, Xcorr = R)$cs
  expect_length(cs, 1L)
  expect_true(25L %in% cs[[1]])
  expect_identical(unname(which.max(susieR::susie_get_pip(fit))), 25L)
})

test_that("coloc.abf supports a shared causal variant", {
  skip_if_not_installed("coloc")
  p <- 60
  R <- 0.7^abs(outer(seq_len(p), seq_len(p), "-"))
  snp <- paste0("rs", seq_len(p))
  gw <- data.table(SNP = snp, beta = 0.04 * R[, 30], SE = 0.005, EAF_A1 = 0.3)
  qtl <- data.table(SNP = snp, beta = 0.3 * R[, 30], se = 0.02, n = 30000, maf = 0.3)
  r <- abf_pair(gw, qtl, n_gwas = 20000, priors = list(p1 = 1e-4, p2 = 1e-4, p12 = 5e-6))
  expect_gt(r$PP.H4.abf, 0.8)
  expect_identical(r$top_snp, "rs30")
})
