test_that("h2 Z-score filter works at threshold = 2", {
  h2_table <- data.table(
    trait = c("D50", "E10", "E14"),
    h2 = c(0.05, 0.03, 0.01),
    se = c(0.01, 0.015, 0.02)
  )
  h2_table[, z := h2 / se]
  result <- filter_h2(h2_table, z_threshold = 2)
  expect_equal(result$trait, c("D50", "E10"))
})

test_that("h2 filter retains exact threshold", {
  h2_table <- data.table(trait = "A", h2 = 0.04, se = 0.02, z = 2.0)
  result <- filter_h2(h2_table, z_threshold = 2.0)
  expect_equal(nrow(result), 1)
})

test_that("h2 filter computes z if not present", {
  h2_table <- data.table(trait = c("A", "B"), h2 = c(0.1, 0.01), se = c(0.02, 0.02))
  result <- filter_h2(h2_table, z_threshold = 2)
  expect_equal(result$trait, "A")
})

test_that("S matrix PSD check passes for valid matrix", {
  good_S <- matrix(c(1, 0.5, 0.5, 1), 2, 2)
  expect_true(check_psd(good_S))
})

test_that("S matrix PSD check fails for invalid matrix", {
  bad_S <- matrix(c(1, 1.5, 1.5, 1), 2, 2)
  expect_false(check_psd(bad_S))
})

test_that("S matrix PSD check works for 3x3", {
  S <- matrix(c(
    0.05, 0.02, 0.01,
    0.02, 0.04, 0.015,
    0.01, 0.015, 0.03
  ), 3, 3)
  expect_true(check_psd(S))
})

test_that("PSD check handles near-zero eigenvalues", {
  S <- matrix(c(1, 0.9999, 0.9999, 1), 2, 2)
  expect_true(check_psd(S))
})
