source(file.path(project_root, "R/09_write_vcf.R"))
source(file.path(project_root, "R/10_plots.R"))
source(file.path(project_root, "R/11_factor_summary.R"))

test_that("compute_lambda_gc returns ~1 for uniform p-values", {
  set.seed(1)
  p <- runif(50000)
  expect_equal(compute_lambda_gc(p), 1, tolerance = 0.05)
})

test_that("compute_lambda_gc is > 1 for inflated z-scores", {
  set.seed(2)
  z <- rnorm(5000, mean = 1.5)
  p <- 2 * pnorm(-abs(z))
  expect_gt(compute_lambda_gc(p), 1.5)
})

test_that("compute_lambda_gc is NA for too-small input", {
  expect_true(is.na(compute_lambda_gc(c(0.1, 0.5))))
})

test_that("neg_log10_p recovers P values that underflowed to 0 from Z", {
  expect_equal(neg_log10_p(c(0.1, 1e-10)), c(1, 10))
  expect_equal(neg_log10_p(0, 40), -(pnorm(-40, log.p = TRUE) + log(2)) / log(10))
  expect_gt(neg_log10_p(0, 40), 300)
})

test_that("chrom_labels maps PLINK codes 23-26 and strips chr prefixes", {
  expect_equal(chrom_labels(c(1, 22, 23, 24, 25, 26)), c("1", "22", "X", "Y", "XY", "MT"))
  expect_equal(chrom_labels(c("chr6", "chrX")), c("6", "X"))
})

test_that("p_label formats thresholds for plotmath and plain text", {
  expect_equal(p_label(5e-8), "5 %*% 10^-8")
  expect_equal(p_label(1e-5), "10^-5")
  expect_equal(p_label(5e-8, plotmath = FALSE), "5×10⁻⁸")
})

test_that("parse_indicators reads each factor's indicators from lavaan syntax", {
  ind <- parse_indicators("ARD =~ D12 + E03\nCM =~ NA*E11 + E66\nARD ~~ 0*CM\nD12 ~~ rvv_D12*D12")
  expect_equal(ind, list(ARD = c("D12", "E03"), CM = c("E11", "E66")))
})

test_that("distance_clump merges SNPs within the window and counts the MHC once", {
  d <- data.table(CHR = c("1", "1", "1", "6", "6", "6"),
                  BP = c(1e6, 1.3e6, 3e6, 26e6, 30e6, 34e6),
                  lp = c(10, 9, 8, 12, 11, 9))
  leads <- distance_clump(d, window = 5e5, mhc = list(chr = 6, start = 25e6, end = 35e6))
  expect_equal(sort(leads$BP), c(1e6, 3e6, 26e6))
  expect_equal(nrow(distance_clump(d, window = 5e5)), 5L)
})

test_that("compress_above_transform is linear below the break and invertible above it", {
  tr <- compress_above_transform(30, 0.1)
  expect_equal(tr$transform(c(10, 30, 130)), c(10, 30, 40))
  expect_equal(tr$inverse(tr$transform(c(5, 50, 97))), c(5, 50, 97))
})

test_that("plot_manhattan draws every SNP below keep_threshold and thins the rest", {
  set.seed(5)
  n <- 60000
  g <- data.table(CHR = sample(1:23, n, replace = TRUE),
                  BP = sample.int(5e7, n, replace = TRUE),
                  ES = rnorm(n), SE = 1, Q_SNP_pval = runif(n))
  g[, P := 2 * pnorm(-abs(ES / SE))]
  # An underflowed hit (P = 0, Z = 40) and a Q_SNP-heterogeneous genome-wide hit.
  g[1:2, `:=`(CHR = 2L, BP = c(1e6, 2e6), ES = c(40, 9), P = c(0, 2 * pnorm(-9)),
              Q_SNP_pval = c(0.5, 1e-9))]
  p <- plot_manhattan(g, factor_name = "F1", sex = "bothsex", width = 3, height = 1, dpi = 50)
  lp <- neg_log10_p(g$P, g$ES / g$SE)
  expect_equal(sum(p$data$lp >= 3), sum(lp >= 3))
  expect_lt(nrow(p$data), n)
  expect_gt(max(p$data$lp), 300)
  expect_true("X" %in% p$scales$get_scales("x")$labels)
  het <- Filter(function(l) inherits(l$geom, "GeomPoint") && is.data.frame(l$data), p$layers)
  expect_equal(nrow(het[[1L]]$data), 1L)
  expect_no_error(ggplot2::ggplot_build(p))
})

test_that("factor_summary_row produces expected columns", {
  set.seed(3)
  n <- 100
  g <- data.table(
    SNP = paste0("rs", 1:n),
    CHR = "1", BP = 1:n,
    A1 = "A", A2 = "G",
    AF = runif(n, 0.05, 0.5),
    ES = rnorm(n, 0, 0.02),
    SE = runif(n, 0.005, 0.015),
    P = runif(n),
    Q_SNP = abs(rnorm(n, 0, 2)),
    Q_SNP_df = 1L,
    Q_SNP_pval = runif(n)
  )
  row <- factor_summary_row(g, "TestFactor")
  expect_equal(row$factor, "TestFactor")
  expect_equal(row$n_snps, n)
  expect_true(is.finite(row$mean_chisq))
  expect_true(is.finite(row$lambda_gc))
  expect_true(row$n_genome_wide_sig >= 0)
  expect_true(row$n_qsnp_het_sig >= 0)
})

test_that("top_snps_table returns the smallest-P rows", {
  g <- data.table(
    SNP = paste0("rs", 1:50),
    CHR = "1", BP = 1:50,
    A1 = "A", A2 = "G",
    AF = 0.1, ES = 0, SE = 0.01,
    P = (50:1) / 100,
    Q_SNP = 0, Q_SNP_pval = 0.5
  )
  top <- top_snps_table(g, "F1", n_top = 5)
  expect_equal(nrow(top), 5L)
  expect_equal(top$SNP[1], "rs50")
  expect_lt(top$P[5], top$P[1] + 1)
})

test_that("dot_path_diagram emits an edge per loading row", {
  model <- "Blood =~ D50 + D64\nEndocrine =~ E10 + E11"
  ld <- data.table(
    trait = c("D50", "D64", "E10", "E11"),
    factor = c("Blood", "Blood", "Endocrine", "Endocrine"),
    loading = c(0.7, 0.6, 0.8, 0.75)
  )
  dot <- dot_path_diagram(model, ld, "male")
  expect_match(dot, "^digraph")
  expect_match(dot, "Blood")
  expect_match(dot, "Endocrine")
  expect_true(grepl("\"Blood\" -> \"D50\"", dot, fixed = TRUE))
  expect_true(grepl("\"Endocrine\" -> \"E11\"", dot, fixed = TRUE))
  expect_true(grepl("label=\"0.70\"", dot, fixed = TRUE))
})

test_that("VCF FILTER is always PASS (no variant-level Q_SNP_HET); HP is exposed per-factor", {
  set.seed(4)
  snps <- paste0("rs", 1:30)
  mk <- function(qp) data.table(
    SNP = snps, CHR = "1", BP = seq_along(snps),
    A1 = "A", A2 = "G", MAF = 0.2,
    est = rnorm(30, 0, 0.01), SE = runif(30, 0.005, 0.01),
    Pval_Estimate = runif(30),
    Q_SNP = abs(rnorm(30)), Q_SNP_df = 1L, Q_SNP_pval = qp
  )
  qp <- runif(30)
  qp[c(5, 10, 15)] <- 1e-10
  g1 <- mk(qp); g2 <- mk(runif(30))
  tmp <- tempfile(fileext = ".vcf.gz")
  on.exit(unlink(c(tmp, paste0(tmp, ".tbi"),
                   sub("\\.gz$", "", tmp))), add = TRUE)
  res <- write_gwas_vcf(list(A = g1, B = g2), NULL, tmp,
                        metadata = list(sex = "male", std_lv = TRUE))

  lines <- readLines(gzfile(res$path))
  hdr <- lines[startsWith(lines, "#")]
  data_rows <- lines[!startsWith(lines, "#")]

  # No variant-level FILTER for Q_SNP_HET anywhere.
  expect_false(any(grepl("^##FILTER=<ID=Q_SNP_HET", hdr)))
  filters <- vapply(strsplit(data_rows, "\t", fixed = TRUE), `[[`, character(1), 7L)
  expect_true(all(filters == "PASS"))

  # The companion _qc.vcf.gz must NOT be emitted.
  qc_companion <- sub("\\.vcf\\.gz$", "_qc.vcf.gz", res$path)
  expect_false(file.exists(qc_companion))

  # HP (per-factor Q_SNP_pval) is still a FORMAT field carrying the per-(SNP,factor) info.
  expect_true(any(grepl("^##FORMAT=<ID=HP,", hdr)))
})

test_that("VCF FORMAT field includes HET, HP, and SF", {
  snps <- paste0("rs", 1:5)
  g <- data.table(
    SNP = snps, CHR = "1", BP = 1:5, A1 = "A", A2 = "G", MAF = 0.1,
    est = rnorm(5), SE = runif(5, 0.01, 0.02), Pval_Estimate = runif(5),
    Q_SNP = c(1, 2, 3, 4, 5), Q_SNP_df = 1L,
    Q_SNP_pval = c(0.5, 0.3, 0.1, 0.01, 1e-9),
    fail = c(FALSE, FALSE, TRUE, FALSE, FALSE)
  )
  tmp <- tempfile(fileext = ".vcf.gz")
  on.exit(unlink(c(tmp, paste0(tmp, ".tbi"),
                   sub("\\.gz$", "", tmp))), add = TRUE)
  res <- write_gwas_vcf(list(A = g), NULL, tmp,
                        metadata = list(sex = "male", std_lv = TRUE))
  lines <- readLines(gzfile(res$path))
  hdr <- lines[startsWith(lines, "#")]
  expect_true(any(grepl("^##FORMAT=<ID=HET,", hdr)))
  expect_true(any(grepl("^##FORMAT=<ID=HP,", hdr)))
  expect_true(any(grepl("^##FORMAT=<ID=SF,", hdr)))
  expect_true(any(grepl("^##gsem_lambda_gc_A=", hdr)))
})
