source(file.path(project_root, "R/04_cfa.R"))

test_that("extract_fit reads uppercase CFI/RMSEA/SRMR columns", {
  mr <- list(modelfit = data.frame(
    chisq = 12.3, df = 4L, CFI = 0.95, SRMR = 0.04, RMSEA = 0.06,
    stringsAsFactors = FALSE
  ))
  out <- extract_fit(mr)
  expect_equal(out$cfi, 0.95)
  expect_equal(out$rmsea, 0.06)
  expect_equal(out$srmr, 0.04)
  expect_equal(out$chisq, 12.3)
  expect_equal(out$df, 4L)
})

test_that("extract_fit falls back to lowercase if uppercase missing", {
  mr <- list(modelfit = data.frame(
    chisq = 12.3, df = 4L, cfi = 0.91, srmr = 0.05, rmsea = 0.07,
    stringsAsFactors = FALSE
  ))
  out <- extract_fit(mr)
  expect_equal(out$cfi, 0.91)
  expect_equal(out$rmsea, 0.07)
  expect_equal(out$srmr, 0.05)
})

test_that("extract_fit returns NA list when modelfit is NULL", {
  out <- extract_fit(list())
  expect_true(is.na(out$cfi))
  expect_true(is.na(out$rmsea))
  expect_true(is.na(out$srmr))
})

test_that("extract_fit returns NA per-field when a column is missing", {
  mr <- list(modelfit = data.frame(chisq = 12.3, df = 4L, stringsAsFactors = FALSE))
  out <- extract_fit(mr)
  expect_equal(out$chisq, 12.3)
  expect_equal(out$df, 4L)
  expect_true(is.na(out$cfi))
})
