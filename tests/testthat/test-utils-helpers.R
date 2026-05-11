test_that("harmonic_neff matches the harmonic mean", {
  x <- c(1000, 2000, 4000)
  expected <- length(x) / sum(1 / x)
  expect_equal(harmonic_neff(x), expected)
})

test_that("harmonic_neff drops NA and non-positive entries", {
  expect_equal(harmonic_neff(c(1000, NA, 0, -5, 2000)),
               2 / (1/1000 + 1/2000))
})

test_that("harmonic_neff returns NA on all-bad input", {
  expect_true(is.na(harmonic_neff(c(NA, NA))))
  expect_true(is.na(harmonic_neff(c(0, -1))))
  expect_true(is.na(harmonic_neff(numeric(0))))
})

test_that("vech_diag_index gives the leading-diagonal positions of vech(S)", {
  expect_equal(vech_diag_index(1), 1L)
  expect_equal(vech_diag_index(2), c(1L, 3L))
  expect_equal(vech_diag_index(3), c(1L, 4L, 6L))
  expect_equal(vech_diag_index(4), c(1L, 5L, 8L, 10L))
})

test_that("vech_diag_index length matches k", {
  for (k in 1:8) expect_length(vech_diag_index(k), k)
})

test_that("vech_diag_index rejects k < 1", {
  expect_error(vech_diag_index(0))
})
