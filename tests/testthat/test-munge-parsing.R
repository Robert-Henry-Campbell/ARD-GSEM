test_that("variant column splits correctly into chr:pos:ref:alt", {
  dt <- data.table(variant = c("1:69487:G:A", "22:50000:C:T"))
  dt <- parse_variant_column(dt)
  expect_equal(dt$chr, c("1", "22"))
  expect_equal(dt$pos, c(69487L, 50000L))
  expect_equal(dt$ref, c("G", "C"))
  expect_equal(dt$alt, c("A", "T"))
})

test_that("multi-character alleles parse correctly", {
  dt <- data.table(variant = "1:100:AT:G")
  dt <- parse_variant_column(dt)
  expect_equal(dt$ref, "AT")
  expect_equal(dt$alt, "G")
})

test_that("A1 = alt (effect allele), A2 = ref after parsing", {
  dt <- data.table(variant = "1:69487:G:A", beta = 0.05, se = 0.01)
  dt <- parse_variant_column(dt)
  dt[, A1 := alt]
  dt[, A2 := ref]
  expect_equal(dt$A1, "A")
  expect_equal(dt$A2, "G")
})

test_that("rsid mapping works and drops unmapped variants", {
  sumstats <- data.table(variant = c("1:69487:G:A", "1:99999:X:Y"))
  manifest <- data.table(variant = c("1:69487:G:A"), rsid = c("rs568226429"))
  setkey(manifest, "variant")
  result <- map_rsids(sumstats, manifest)
  expect_equal(nrow(result), 1)
  expect_equal(result$SNP, "rs568226429")
})

test_that("rsid mapping preserves all columns", {
  sumstats <- data.table(variant = "1:69487:G:A", beta = 0.05, se = 0.01, pval = 0.001)
  manifest <- data.table(variant = "1:69487:G:A", rsid = "rs568226429")
  setkey(manifest, "variant")
  result <- map_rsids(sumstats, manifest)
  expect_true("beta" %in% names(result))
  expect_equal(result$beta, 0.05)
})

test_that("sentinel SNP polarity check passes for known variant", {
  dt <- data.table(SNP = "rs568226429", A1 = "A", A2 = "G", beta = 0.05)
  result <- check_polarity_sentinel(dt, sentinel_rsid = "rs568226429", expected_a1 = "A")
  expect_true(result$polarity_correct)
})

test_that("sentinel SNP polarity check fails when A1 is wrong", {
  dt <- data.table(SNP = "rs568226429", A1 = "G", A2 = "A", beta = -0.05)
  result <- check_polarity_sentinel(dt, sentinel_rsid = "rs568226429", expected_a1 = "A")
  expect_false(result$polarity_correct)
})

test_that("sentinel SNP polarity returns NA when not found", {
  dt <- data.table(SNP = "rs999999", A1 = "A", A2 = "G", beta = 0.05)
  result <- check_polarity_sentinel(dt, sentinel_rsid = "rs568226429", expected_a1 = "A")
  expect_true(is.na(result$polarity_correct))
})

test_that("parse_variant_column errors on variants without 4 colon-separated fields", {
  dt <- data.table(variant = c("1:69487:G", "22:50000:C:T"))
  expect_error(parse_variant_column(dt), "4 fields")
})

test_that("parse_variant_column errors on empty fields", {
  dt <- data.table(variant = c("1:69487::A"))
  expect_error(parse_variant_column(dt), "NA/empty")
})
