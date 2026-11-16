for (f in c("R/09_write_vcf.R", "R/11_factor_summary.R", "R/12_factor_ldsc.R", "R/13_qc_filters.R",
            "R/15_postgwas_common.R", "R/16_loci.R", "R/17_mr_export.R", "R/18_annotate.R")) {
  source(file.path(project_root, f), local = TRUE)
}

test_that("upstream stage hashes ignore the postgwas/refs/tools blocks; post-GWAS stages include them", {
  cfg <- list(paths = list(output_dir = "/x"), gwas = list(chunk_size = 1),
              postgwas = list(factors = "ARD"), refs = list(a = "/r"), tools = list(plink = "/p"))
  up <- stage_hash_payload("gwas", "bothsex", cfg)
  expect_identical(names(up), c("stage", "paths", "source"))
  cfg2 <- cfg; cfg2$postgwas$factors <- c("ARD", "CM"); cfg2$tools$plink <- "/q"
  expect_identical(digest::digest(up), digest::digest(stage_hash_payload("gwas", "bothsex", cfg2)))
  pg <- stage_hash_payload("loci", "bothsex", cfg)
  expect_true(all(c("postgwas", "refs", "tools") %in% names(pg)))
  expect_false(identical(digest::digest(pg), digest::digest(stage_hash_payload("loci", "bothsex", cfg2))))
})

test_that("lavaan messages are classed as ok / gradient / fail / other", {
  w <- c("0", NA,
         "lavaan->lav_lavaan_step02_options():  \n   the following argument(s) override(s) the options in slotOptions: se",
         "lavaan->lav_lavaan_step11_estoptim():  \n   the optimizer (NLMINB) claimed the model converged, but not all elements \n   of the gradient are (near) zero; ...",
         "lavaan->lav_lavaan_step11_estoptim():  \n   the optimizer warns that a solution has NOT been found!",
         "something new")
  cls <- classify_gsem_warnings(w, "solution has NOT been found",
                                c("override\\(s\\) the options in slotOptions", "gradient are \\(near\\) zero"))
  expect_identical(cls, c("ok", "ok", "ok", "gradient", "fail", "other"))
})

test_that("log-space p-values match the direct ones and survive underflow", {
  z <- c(0.5, 3, 8)
  expect_equal(log10p_from_z(z), -log10(2 * pnorm(-abs(z))), tolerance = 1e-10)
  expect_true(is.finite(log10p_from_z(60)) && log10p_from_z(60) > 700)
  expect_equal(log10p_chisq(50, 30), -log10(pchisq(50, 30, lower.tail = FALSE)), tolerance = 1e-10)
})

test_that("merge_intervals merges overlaps per chromosome only", {
  m <- merge_intervals(c(1, 1, 1, 2), c(100, 150, 500, 120), c(200, 300, 600, 130))
  expect_equal(nrow(m), 3)
  expect_equal(m[chr == 1, start], c(100, 500))
  expect_equal(sort(unlist(m[chr == 1][1, rows])), c(1L, 2L))
})

test_that("choose_marker picks the highest standardised loading or a named indicator", {
  lt <- data.table(factor = c("A", "A", "B"), indicator = c("x", "y", "z"),
                   lambda = c(0.1, 0.3, 0.2), std = c(0.9, 0.5, 0.7))
  expect_identical(choose_marker(lt, "A")$indicator, "x")
  expect_identical(choose_marker(lt, "A", "y")$lambda, 0.3)
  expect_error(choose_marker(lt, "A", "z"))
})

test_that("Q_SNP partners: exact LD for panel SNPs, a window for SNPs without LD information", {
  dt <- data.table(SNP = paste0("rs", 1:8), CHR = c(rep(1L, 6), 23L, 23L),
                   BP = c(1e6, 1.1e6, 1.2e6, 5e6, 5.3e6, 9e6, 2e6, 2.4e6),
                   in_1000g = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE))
  plink_ld_partners <- function(config, bfile, snps, r2_min, kb, out_prefix, stage, threads = 8L) {
    expect_identical(sort(snps), "rs1")
    data.table(SNP_A = "rs1", SNP_B = c("rs2", "rs1"), R2 = c(0.5, 1))
  }
  environment(flag_qsnp_partners) <- environment()
  fl <- flag_qsnp_partners(list(postgwas = list(workers = 1)), dt, c("rs1", "rs7"), "bfile",
                           0.1, 1000, 500, tempfile())
  expect_identical(which(fl$qsnp_ld), c(2L, 8L))            # rs2 by r2; rs8 within 500 kb of chrX hit rs7
  expect_identical(fl$qsnp_ld_method[c(2, 8)], c("r2", "window"))
  expect_false(any(fl$qsnp_ld[c(1, 7)]))                    # hits themselves are not partners
})

test_that("munged LDSC files align alleles to HapMap3 and drop mismatches", {
  hm3 <- data.table(SNP = c("rs1", "rs2", "rs3"), A1 = c("A", "C", "G"), A2 = c("G", "T", "T"))
  dt <- data.table(SNP = c("rs1", "rs2", "rs3", "rs4"), A1 = c("A", "T", "G", "A"),
                   A2 = c("G", "C", "C", "C"), Z = c(2, 3, 4, 5))
  p <- tempfile(fileext = ".sumstats.gz")
  expect_equal(write_munged_sumstats(dt, hm3, 1000, p), 2)
  out <- fread(p)
  expect_equal(out[SNP == "rs1", Z], 2)
  expect_equal(out[SNP == "rs2", c(Z, A1)], c("-3", "C"))
})

test_that("ardmr export helpers: palindromes and the column check", {
  expect_identical(palindrome(c("A", "C", "A"), c("T", "G", "G")), c(TRUE, TRUE, FALSE))
  x <- as.data.table(setNames(as.list(rep(1, length(ARDMR_EXPOSURE_COLS))), ARDMR_EXPOSURE_COLS))
  expect_identical(validate_ardmr_exposure(x, ""), "built-in column check")
  expect_error(validate_ardmr_exposure(x[, !"pval.exposure"], ""))
})

test_that("distance clumping keeps the strongest SNP per window", {
  d <- data.table(SNP = c("a", "b", "c", "d"), CHR = c(1L, 1L, 1L, 2L),
                  BP = c(1e6, 1.3e6, 2e6, 1e6), LOG10P = c(9, 12, 8, 10))
  expect_identical(sort(distance_clump_simple(d, 500)$SNP), c("b", "c", "d"))
})

test_that("nearest TSS and overlapping genes from a toy GENCODE table", {
  genes <- data.table(gene_id = c("G1", "G2"), symbol = c("AAA", "BBB"), CHR = c(1L, 1L),
                      start = c(1000L, 5000L), end = c(2000L, 9000L), strand = c("+", "-"),
                      TSS = c(1000L, 9000L))
  v <- data.table(CHR = c(1L, 1L), BP = c(1500L, 8000L))
  ng <- nearest_genes_gencode(v, genes)
  expect_identical(ng$nearest_tss_gene, c("AAA", "BBB"))
  expect_identical(ng$nearest_tss_dist, c(500L, 1000L))
  expect_identical(ng$overlapping_genes, c("AAA", "BBB"))
})

test_that("GenomicSEM ldsc logs are parsed per trait", {
  log <- tempfile()
  writeLines(c("Heritability Results for trait: /x/ARD_qfilter.sumstats.gz",
               "Mean Chi^2 across remaining SNPs: 1.5461", "Lambda GC: 1.428",
               "Intercept: 0.9904 (0.0091)", "Ratio: -0.0176 (0.0167)",
               "Total Observed Scale h2: 1.6922 (0.0598)", "h2 Z: 28.3", ""), log)
  p <- parse_gsem_ldsc_log(log)
  expect_identical(p$trait_file, "ARD_qfilter.sumstats.gz")
  expect_equal(c(p$intercept, p$intercept_se, p$ratio, p$h2, p$h2_z), c(0.9904, 0.0091, -0.0176, 1.6922, 28.3))
})

test_that("post-GWAS VCF orients REF/ALT to GRCh37 and flips ES/AF when A1 is REF", {
  skip_if_not_installed("Rsamtools")
  tmp <- tempfile(); dir.create(tmp)
  fasta <- file.path(tmp, "fake.fa")
  writeLines(sprintf("chr%s\t1000000\t0\t60\t61", c(1:22, "X")), paste0(fasta, ".fai"))
  cfg <- list(paths = list(output_dir = file.path(tmp, "out")), refs = list(hg19_fasta = fasta),
              postgwas = list(factors = "F1", excluded_factors = list(F2 = "two indicators"), qsnp_p = 5e-8),
              qc_filters = list(), vcf = list(basename = "t"))
  qd <- file.path(tmp, "out", "bothsex", "qc_filters"); dir.create(qd, recursive = TRUE)
  dt <- data.table(CHR = c(1L, 2L), BP = c(100L, 200L), SNP = c("rs1", "rs2"), A1 = c("A", "C"), A2 = c("G", "T"),
                   a1_is_ref = c(TRUE, FALSE), EAF_A1 = c(0.3, 0.2), INFO = c(0.99, 0.95),
                   beta = c(0.1, 0.2), SE = c(0.01, 0.02), LOG10P = c(20, 5), Q_SNP = c(1, 2),
                   Q_SNP_pval = c(0.5, 0.4), pass_full = c(TRUE, TRUE), pass_qfilter = c(TRUE, FALSE))
  saveRDS(dt, file.path(qd, "F1.qc.rds"))
  jsonlite::write_json(list(F1 = list(role = "general", indicators = c("a", "b"), N_hat_std = 1000,
                                      marker = "a", marker_lambda = 0.2, N_hat_marker = 25000,
                                      qsnp_df_gsem = 1, qsnp_df_eff = 1, lambda_gc_full = 1.01)),
                       file.path(qd, "factor_registry.json"), auto_unbox = TRUE)
  res <- write_postgwas_vcf(cfg, "bothsex", file.path(tmp, "vcf"))
  x <- readLines(gzfile(res$paths[1]))
  expect_true(any(grepl("^##gsem_excluded_factor=F2", x)))
  body <- strsplit(grep("^[^#]", x, value = TRUE), "\t")
  expect_identical(body[[1]][4:5], c("A", "G"))                         # A1 is REF
  expect_identical(strsplit(body[[1]][10], ":")[[1]][c(1, 4)], c("-0.1", "0.7"))
  expect_identical(body[[2]][4:5], c("T", "C"))                         # A1 is ALT
  expect_identical(strsplit(body[[2]][10], ":")[[1]][c(1, 4, 8)], c("0.2", "0.2", "1"))
  expect_true(file.exists(paste0(res$paths[1], ".tbi")))
  expect_length(res$paths, 2)
})
