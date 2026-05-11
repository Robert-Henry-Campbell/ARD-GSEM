test_that("polarity check is skipped when sentinel_rsid is NULL", {
  dt <- data.table(SNP = c("rs1", "rs2"), A1 = c("A", "G"))
  out <- check_polarity_sentinel(dt, sentinel_rsid = NULL, expected_a1 = NULL)
  expect_true(is.na(out$polarity_correct))
  expect_match(out$message, "skipped")
})

test_that("polarity check reports correct when A1 matches", {
  dt <- data.table(SNP = c("rs1", "rs2"), A1 = c("A", "G"))
  out <- check_polarity_sentinel(dt, sentinel_rsid = "rs1", expected_a1 = "A")
  expect_true(out$polarity_correct)
})

test_that("polarity check reports incorrect when A1 mismatches", {
  dt <- data.table(SNP = c("rs1", "rs2"), A1 = c("T", "G"))
  out <- check_polarity_sentinel(dt, sentinel_rsid = "rs1", expected_a1 = "A")
  expect_false(out$polarity_correct)
})

test_that("polarity check returns NA when sentinel SNP missing", {
  dt <- data.table(SNP = c("rs2"), A1 = c("G"))
  out <- check_polarity_sentinel(dt, sentinel_rsid = "rs1", expected_a1 = "A")
  expect_true(is.na(out$polarity_correct))
  expect_match(out$message, "not found")
})
