test_that("low_confidence_variant = TRUE rows are removed", {
  dt <- data.table(
    low_confidence_variant = c(TRUE, FALSE, FALSE),
    minor_AF = c(0.1, 0.1, 0.1),
    beta = c(0.1, 0.2, 0.3),
    se = c(0.01, 0.01, 0.01)
  )
  result <- filter_variants(dt)
  expect_equal(nrow(result), 2)
  expect_equal(result$beta, c(0.2, 0.3))
})

test_that("MAF < 0.01 is filtered", {
  dt <- data.table(
    minor_AF = c(0.005, 0.01, 0.3),
    low_confidence_variant = c(FALSE, FALSE, FALSE),
    beta = 1:3,
    se = c(0.01, 0.01, 0.01)
  )
  result <- filter_variants(dt, maf_threshold = 0.01)
  expect_equal(nrow(result), 2)
  expect_equal(result$minor_AF, c(0.01, 0.3))
})

test_that("MAF exactly at threshold is retained", {
  dt <- data.table(
    minor_AF = 0.01,
    low_confidence_variant = FALSE,
    beta = 0.5,
    se = 0.01
  )
  result <- filter_variants(dt, maf_threshold = 0.01)
  expect_equal(nrow(result), 1)
})

test_that("variants with missing beta are dropped", {
  dt <- data.table(
    beta = c(0.1, NA, 0.2),
    se = c(0.01, 0.01, 0.01),
    minor_AF = c(0.1, 0.1, 0.1),
    low_confidence_variant = c(FALSE, FALSE, FALSE)
  )
  result <- filter_variants(dt)
  expect_equal(nrow(result), 2)
})

test_that("variants with missing se are dropped", {
  dt <- data.table(
    beta = c(0.1, 0.2, 0.3),
    se = c(0.01, NA, 0.01),
    minor_AF = c(0.1, 0.1, 0.1),
    low_confidence_variant = c(FALSE, FALSE, FALSE)
  )
  result <- filter_variants(dt)
  expect_equal(nrow(result), 2)
})

test_that("info filter works when column exists", {
  dt <- data.table(
    info = c(0.8, 0.95, 0.99),
    minor_AF = c(0.1, 0.1, 0.1),
    low_confidence_variant = c(FALSE, FALSE, FALSE),
    beta = c(0.1, 0.2, 0.3),
    se = c(0.01, 0.01, 0.01)
  )
  result <- filter_variants(dt, info_threshold = 0.9)
  expect_equal(nrow(result), 2)
})

test_that("all filters applied together", {
  dt <- data.table(
    variant = c("1:1:A:T", "1:2:A:T", "1:3:A:T", "1:4:A:T", "1:5:A:T"),
    minor_AF = c(0.001, 0.1, 0.1, 0.1, 0.1),
    low_confidence_variant = c(FALSE, TRUE, FALSE, FALSE, FALSE),
    beta = c(0.1, 0.1, NA, 0.1, 0.1),
    se = c(0.01, 0.01, 0.01, 0.01, 0.01)
  )
  result <- filter_variants(dt, maf_threshold = 0.01)
  expect_equal(nrow(result), 2)
})
