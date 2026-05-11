test_that("loading difference Z-score is computed correctly", {
  result <- compute_loading_diff(
    male_loading = 0.5, male_se = 0.05,
    female_loading = 0.7, female_se = 0.06
  )
  expected_z <- (0.7 - 0.5) / sqrt(0.05^2 + 0.06^2)
  expect_equal(result$z_diff, expected_z, tolerance = 1e-6)
})

test_that("loading difference is female - male", {
  result <- compute_loading_diff(
    male_loading = 0.8, male_se = 0.05,
    female_loading = 0.6, female_se = 0.05
  )
  expect_true(result$diff < 0)
  expect_true(result$z_diff < 0)
})

test_that("equal loadings produce z_diff = 0", {
  result <- compute_loading_diff(
    male_loading = 0.5, male_se = 0.05,
    female_loading = 0.5, female_se = 0.06
  )
  expect_equal(result$z_diff, 0)
  expect_equal(result$p_diff, 1)
})

test_that("p_diff is two-tailed", {
  result <- compute_loading_diff(
    male_loading = 0.3, male_se = 0.05,
    female_loading = 0.8, female_se = 0.05
  )
  expect_true(result$p_diff < 0.05)
  expect_true(result$p_diff > 0)

  result_neg <- compute_loading_diff(
    male_loading = 0.8, male_se = 0.05,
    female_loading = 0.3, female_se = 0.05
  )
  expect_equal(result$p_diff, result_neg$p_diff, tolerance = 1e-10)
})

test_that("se_diff combines in quadrature", {
  result <- compute_loading_diff(
    male_loading = 0.5, male_se = 0.03,
    female_loading = 0.6, female_se = 0.04
  )
  expect_equal(result$se_diff, sqrt(0.03^2 + 0.04^2), tolerance = 1e-10)
})

test_that("extract_loadings_table uses Unstand_Est column", {
  mock_results <- data.table::data.table(
    lhs = c("F1", "F1"),
    op = c("=~", "=~"),
    rhs = c("D50", "D64"),
    Unstand_Est = c(0.45, 0.62),
    Unstand_SE = c(0.05, 0.04),
    STD_Genotype = c(0.50, 0.68),
    STD_Genotype_SE = c(0.06, 0.05),
    STD_All = c(0.55, 0.72),
    p_value = c(1e-10, 1e-15)
  )
  result <- extract_loadings_table(list(results = mock_results))
  expect_equal(result$loading, c(0.45, 0.62))
  expect_equal(result$se, c(0.05, 0.04))
})

test_that("extract_loadings_table returns NULL for missing results", {
  expect_null(extract_loadings_table(list(results = NULL)))
})

test_that("extract_loadings_table returns NULL when no =~ rows", {
  mock_results <- data.table::data.table(
    lhs = "F1", op = "~~", rhs = "F1",
    Unstand_Est = 1, Unstand_SE = 0.1,
    STD_Genotype = 1, STD_Genotype_SE = 0.1,
    STD_All = 1, p_value = 1e-10
  )
  expect_null(extract_loadings_table(list(results = mock_results)))
})

test_that("check_factor_alignment detects different indicators", {
  male <- data.table::data.table(
    trait = c("D50", "D64"), factor = c("F1", "F1"),
    loading = c(0.5, 0.6), se = c(0.05, 0.04))
  female <- data.table::data.table(
    trait = c("E10", "E11"), factor = c("F1", "F1"),
    loading = c(0.7, 0.8), se = c(0.03, 0.04))
  warnings <- check_factor_alignment(male, female)
  expect_true(length(warnings) > 0)
  expect_true(any(grepl("different indicators", warnings)))
})

test_that("check_factor_alignment detects different factor counts", {
  male <- data.table::data.table(
    trait = c("D50", "D64", "E10"), factor = c("F1", "F1", "F2"),
    loading = c(0.5, 0.6, 0.7), se = c(0.05, 0.04, 0.03))
  female <- data.table::data.table(
    trait = c("D50", "D64"), factor = c("F1", "F1"),
    loading = c(0.5, 0.6), se = c(0.05, 0.04))
  warnings <- check_factor_alignment(male, female)
  expect_true(any(grepl("Different number of factors", warnings)))
})

test_that("check_factor_alignment is silent for matching structures", {
  male <- data.table::data.table(
    trait = c("D50", "D64"), factor = c("F1", "F1"),
    loading = c(0.5, 0.6), se = c(0.05, 0.04))
  female <- data.table::data.table(
    trait = c("D50", "D64"), factor = c("F1", "F1"),
    loading = c(0.7, 0.8), se = c(0.03, 0.04))
  warnings <- check_factor_alignment(male, female)
  expect_equal(length(warnings), 0)
})
