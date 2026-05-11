source(file.path(project_root, "R/09_write_vcf.R"))
source(file.path(project_root, "R/10_plots.R"))
source(file.path(project_root, "R/11_factor_summary.R"))

test_that("compute_lambda_gc returns ~1 for uniform p-values", {
  set.seed(1)
  p <- runif(5000)
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

test_that("VCF Q_SNP_HET FILTER fires when Q_SNP_pval below threshold", {
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
  on.exit(unlink(c(tmp, sub("\\.vcf\\.gz$", "_qc.vcf.gz", tmp))), add = TRUE)
  res <- write_gwas_vcf(list(A = g1, B = g2), NULL, tmp,
                        metadata = list(sex = "male"),
                        qsnp_threshold = 5e-8,
                        write_qc_filtered = TRUE)
  expect_equal(res$n_qsnp_het, 3L)
  expect_true(!is.null(res$qc_path) && file.exists(res$qc_path))

  lines <- readLines(gzfile(tmp))
  data_rows <- lines[!startsWith(lines, "#")]
  filters <- vapply(strsplit(data_rows, "\t", fixed = TRUE), `[[`, character(1), 7L)
  expect_equal(sum(filters == "Q_SNP_HET"), 3L)

  qc_lines <- readLines(gzfile(res$qc_path))
  qc_data <- qc_lines[!startsWith(qc_lines, "#")]
  expect_equal(length(qc_data), nrow(g1) - 3L)
})

test_that("VCF FORMAT field includes HET and HP", {
  snps <- paste0("rs", 1:5)
  g <- data.table(
    SNP = snps, CHR = "1", BP = 1:5, A1 = "A", A2 = "G", MAF = 0.1,
    est = rnorm(5), SE = runif(5, 0.01, 0.02), Pval_Estimate = runif(5),
    Q_SNP = c(1, 2, 3, 4, 5), Q_SNP_df = 1L,
    Q_SNP_pval = c(0.5, 0.3, 0.1, 0.01, 1e-9)
  )
  tmp <- tempfile(fileext = ".vcf.gz")
  on.exit(unlink(c(tmp, sub("\\.vcf\\.gz$", "_qc.vcf.gz", tmp))), add = TRUE)
  write_gwas_vcf(list(A = g), NULL, tmp,
                 metadata = list(sex = "male"),
                 qsnp_threshold = 5e-8,
                 write_qc_filtered = TRUE)
  lines <- readLines(gzfile(tmp))
  hdr <- lines[startsWith(lines, "#")]
  expect_true(any(grepl("^##FORMAT=<ID=HET,", hdr)))
  expect_true(any(grepl("^##FORMAT=<ID=HP,", hdr)))
  expect_true(any(grepl("^##FILTER=<ID=Q_SNP_HET,", hdr)))
  expect_true(any(grepl("^##gsem_lambda_gc_A=", hdr)))
})
