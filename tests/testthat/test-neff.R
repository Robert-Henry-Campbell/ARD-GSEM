test_that("Neff formula is correct for known values", {
  expect_equal(compute_neff(1000, 99000), 4 * 1000 * 99000 / (1000 + 99000), tolerance = 0.01)
})

test_that("Neff for equal case/control", {
  expect_equal(compute_neff(5000, 5000), 10000)
})

test_that("Neff for very imbalanced (E14 male: 145 cases, 166875 controls)", {
  expected <- 4 * 145 * 166875 / (145 + 166875)
  expect_equal(compute_neff(145, 166875), expected, tolerance = 0.01)
  expect_true(compute_neff(145, 166875) < 600)
})

test_that("Neff errors on zero cases", {
  expect_error(compute_neff(0, 10000))
})

test_that("Neff errors on negative controls", {

  expect_error(compute_neff(100, -5))
})

test_that("Neff is vectorized", {
  result <- compute_neff(c(100, 200), c(9900, 9800))
  expect_length(result, 2)
  expect_equal(result[1], 4 * 100 * 9900 / 10000)
})
