source(file.path(project_root, "R/09_write_vcf.R"))

make_factor_df <- function(snps, with_fail = FALSE) {
  set.seed(42)
  d <- data.table(
    SNP = snps,
    CHR = sample(c("1","2","3"), length(snps), replace = TRUE),
    BP  = sample.int(1e8, length(snps)),
    A1  = sample(c("A","C","G","T"), length(snps), replace = TRUE),
    A2  = sample(c("A","C","G","T"), length(snps), replace = TRUE),
    MAF = runif(length(snps), 0.05, 0.45),
    est = rnorm(length(snps), 0, 0.05),
    SE  = runif(length(snps), 0.005, 0.02),
    Pval_Estimate = runif(length(snps))
  )
  if (with_fail) {
    d[, fail := sample(c(TRUE, FALSE), length(snps), replace = TRUE, prob = c(0.1, 0.9))]
  }
  d
}

# read a possibly-bgzipped or plain-gzipped VCF (round-trip both forms)
read_vcf_lines <- function(path) {
  con <- gzfile(path, open = "rt")
  on.exit(close(con), add = TRUE)
  readLines(con)
}

test_that("write_gwas_vcf produces a parseable GWAS-VCF round-trip", {
  snps <- paste0("rs", 1:20)
  gwas_list <- list(
    Blood = make_factor_df(snps),
    Endocrine = make_factor_df(snps)
  )

  tmp <- tempfile(fileext = ".vcf.gz")
  on.exit(unlink(c(tmp, paste0(tmp, ".tbi"),
                   sub("\\.gz$", "", tmp))), add = TRUE)

  res <- write_gwas_vcf(gwas_list,
                        snp_sumstats = NULL,
                        output_path = tmp,
                        metadata = list(sex = "male",
                                        genome_build = "GRCh37",
                                        study_id_prefix = "ARD-GSEM",
                                        std_lv = TRUE))
  expect_true(file.exists(res$path))
  expect_equal(res$n_factors, 2L)
  expect_equal(res$n_variants, 20L)

  lines <- read_vcf_lines(res$path)
  hdr <- lines[startsWith(lines, "#")]
  body <- lines[!startsWith(lines, "#")]

  expect_true(any(grepl("^##fileformat=VCFv4\\.2$", hdr)))
  expect_true(any(grepl("^##gsem_scale=log_odds_ratio$", hdr)))
  expect_true(any(grepl("^##gsem_model=one_factor_per_ICD10_chapter$", hdr)))
  expect_true(any(grepl("^##gsem_sex=male$", hdr)))
  expect_true(any(grepl("^##gsem_identification=std\\.lv", hdr)))

  # Variant-level FILTER for Q_SNP_HET must be GONE; HP remains as a per-factor FORMAT field.
  expect_false(any(grepl("^##FILTER=<ID=Q_SNP_HET", hdr)))

  format_ids <- c("ES", "SE", "LP", "AF", "SS", "HET", "HP", "SF")
  for (id in format_ids) {
    expect_true(any(grepl(sprintf("^##FORMAT=<ID=%s,", id), hdr)),
                info = sprintf("missing FORMAT for %s", id))
  }

  sample_lines <- hdr[startsWith(hdr, "##SAMPLE=")]
  expect_equal(length(sample_lines), 2L)
  expect_true(any(grepl("ID=Blood", sample_lines)))
  expect_true(any(grepl("ID=Endocrine", sample_lines)))

  col_line <- hdr[startsWith(hdr, "#CHROM")]
  expect_length(col_line, 1L)
  cols <- strsplit(col_line, "\t", fixed = TRUE)[[1]]
  expect_equal(tail(cols, 2), c("Blood", "Endocrine"))

  expect_equal(length(body), 20L)
  parts <- strsplit(body[1], "\t", fixed = TRUE)[[1]]
  expect_equal(length(parts), 9L + 2L)
  # FILTER is always PASS (no variant-level Q_SNP_HET).
  expect_equal(parts[7], "PASS")
  expect_equal(parts[9], "ES:SE:LP:AF:SS:HET:HP:SF")
  sample_field <- strsplit(parts[10], ":", fixed = TRUE)[[1]]
  expect_equal(length(sample_field), 8L)
})

test_that("write_gwas_vcf rejects mismatched SNP sets across factors", {
  g1 <- make_factor_df(paste0("rs", 1:10))
  g2 <- make_factor_df(paste0("rs", 11:20))
  gwas_list <- list(A = g1, B = g2)
  tmp <- tempfile(fileext = ".vcf.gz")
  on.exit(unlink(c(tmp, paste0(tmp, ".tbi"))), add = TRUE)
  expect_error(write_gwas_vcf(gwas_list, NULL, tmp,
                               list(sex = "male")), "different SNP set")
})

test_that("normalize_gwas_columns handles missing position fields by joining snp_sumstats", {
  snps <- paste0("rs", 1:5)
  g <- data.table(
    SNP = snps,
    est = rnorm(5),
    SE  = runif(5, 0.01, 0.05),
    Pval_Estimate = runif(5)
  )
  ss <- data.table(
    SNP = snps,
    CHR = rep("1", 5),
    BP  = 1:5 * 1000L,
    A1  = rep("A", 5),
    A2  = rep("G", 5),
    MAF = rep(0.1, 5)
  )
  out <- normalize_gwas_columns(g, ss)
  expect_equal(out$CHR, rep("1", 5))
  expect_equal(out$BP, 1:5 * 1000L)
})

test_that("normalize_gwas_columns extracts the per-factor `fail` column", {
  snps <- paste0("rs", 1:10)
  g <- make_factor_df(snps, with_fail = TRUE)
  out <- normalize_gwas_columns(g, NULL)
  expect_true("fail" %in% names(out))
  expect_equal(sum(out$fail, na.rm = TRUE), sum(g$fail, na.rm = TRUE))
})

test_that("SF FORMAT field is per-factor: failed rows write 1, converged rows 0", {
  snps <- paste0("rs", 1:10)
  g <- make_factor_df(snps, with_fail = TRUE)
  # Make sure CHR/BP are deterministic so we can match SNP -> output row order
  g[, CHR := "1"]
  g[, BP := 1:10]
  gwas_list <- list(A = g)
  tmp <- tempfile(fileext = ".vcf.gz")
  on.exit(unlink(c(tmp, paste0(tmp, ".tbi"))), add = TRUE)
  res <- write_gwas_vcf(gwas_list, NULL, tmp, list(sex = "male", std_lv = TRUE))

  lines <- read_vcf_lines(res$path)
  body <- lines[!startsWith(lines, "#")]
  # Parse SF from each row (last colon-separated field of the single sample column).
  sf <- vapply(body, function(L) {
    parts <- strsplit(L, "\t", fixed = TRUE)[[1]]
    sample_fields <- strsplit(parts[10], ":", fixed = TRUE)[[1]]
    sample_fields[length(sample_fields)]
  }, character(1))
  # Output is sorted by CHR then BP; we set BP=1..10 so order is preserved.
  expect_equal(sum(sf == "1"), sum(g$fail))
  expect_equal(sum(sf == "0"), sum(!g$fail))
})
