test_that("at K=0.5 the scale factor is 0.25 so logOR = 4 * beta_linear", {
  out <- linear_to_logor(beta = 0.01, se = 0.002, K = 0.5)
  expect_equal(out$beta, 0.04)
  expect_equal(out$se, 0.008)
  expect_equal(out$scale, 0.25)
})

test_that("Z = beta/SE is invariant under the conversion", {
  for (K in c(0.05, 0.2, 0.5, 0.8)) {
    out <- linear_to_logor(beta = 0.005, se = 0.001, K = K)
    expect_equal(out$beta / out$se, 5)
  }
})

test_that("vectorised beta/se with scalar K", {
  out <- linear_to_logor(beta = c(0.01, 0.02, -0.005),
                         se = c(0.002, 0.003, 0.001),
                         K = 0.1)
  expect_length(out$beta, 3)
  expect_equal(out$beta[1], 0.01 / (0.1 * 0.9))
  expect_equal(out$se[2], 0.003 / (0.1 * 0.9))
})

test_that("invalid K is rejected", {
  expect_error(linear_to_logor(0.01, 0.002, K = 0), "invalid K")
  expect_error(linear_to_logor(0.01, 0.002, K = 1), "invalid K")
  expect_error(linear_to_logor(0.01, 0.002, K = -0.1), "invalid K")
  expect_error(linear_to_logor(0.01, 0.002, K = NA), "invalid K")
})

test_that("extreme K with K(1-K) below threshold is rejected", {
  expect_error(linear_to_logor(0.01, 0.002, K = 1e-7), "too imbalanced")
  expect_error(linear_to_logor(0.01, 0.002, K = 1 - 1e-7), "too imbalanced")
})

test_that("vector K is rejected", {
  expect_error(linear_to_logor(c(0.01, 0.02), c(0.002, 0.003), K = c(0.1, 0.2)),
               "scalar")
})
