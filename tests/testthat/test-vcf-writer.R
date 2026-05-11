source(file.path(project_root, "R/09_write_vcf.R"))

make_factor_df <- function(snps) {
  set.seed(42)
  data.table(
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
}

test_that("write_gwas_vcf produces a parseable GWAS-VCF round-trip", {
  snps <- paste0("rs", 1:20)
  gwas_list <- list(
    Blood = make_factor_df(snps),
    Endocrine = make_factor_df(snps)
  )

  tmp <- tempfile(fileext = ".vcf.gz")
  on.exit(unlink(tmp), add = TRUE)

  res <- write_gwas_vcf(gwas_list,
                        snp_sumstats = NULL,
                        output_path = tmp,
                        metadata = list(sex = "male",
                                        genome_build = "GRCh37",
                                        study_id_prefix = "ARD-GSEM"))
  expect_true(file.exists(tmp))
  expect_equal(res$n_factors, 2L)
  expect_equal(res$n_variants, 20L)

  lines <- readLines(gzfile(tmp))
  hdr <- lines[startsWith(lines, "#")]
  body <- lines[!startsWith(lines, "#")]

  expect_true(any(grepl("^##fileformat=VCFv4\\.2$", hdr)))
  expect_true(any(grepl("^##gsem_scale=log_odds_ratio$", hdr)))
  expect_true(any(grepl("^##gsem_model=one_factor_per_ICD10_chapter$", hdr)))
  expect_true(any(grepl("^##gsem_sex=male$", hdr)))

  format_ids <- c("ES", "SE", "LP", "AF", "SS", "HET", "HP")
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
  expect_equal(parts[9], "ES:SE:LP:AF:SS:HET:HP")
  sample_field <- strsplit(parts[10], ":", fixed = TRUE)[[1]]
  expect_equal(length(sample_field), 7L)
})

test_that("write_gwas_vcf rejects mismatched SNP sets across factors", {
  g1 <- make_factor_df(paste0("rs", 1:10))
  g2 <- make_factor_df(paste0("rs", 11:20))
  gwas_list <- list(A = g1, B = g2)
  tmp <- tempfile(fileext = ".vcf.gz")
  on.exit(unlink(tmp), add = TRUE)
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
