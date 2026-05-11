test_that("efa_to_lavaan produces valid syntax for 2-factor model", {
  loadings <- matrix(c(0.8, 0.7, 0.1, 0.05, 0.1, 0.05, 0.9, 0.85), 4, 2)
  rownames(loadings) <- c("D50", "D64", "E10", "E11")
  model <- efa_to_lavaan(loadings, threshold = 0.3)
  expect_match(model, "F1 =~ D50 \\+ D64")
  expect_match(model, "F2 =~ E10 \\+ E11")
})

test_that("efa_to_lavaan drops traits below threshold", {
  loadings <- matrix(c(0.8, 0.2, 0.1), 3, 1)
  rownames(loadings) <- c("D50", "D64", "E10")
  model <- efa_to_lavaan(loadings, threshold = 0.3)
  expect_match(model, "D50")
  expect_no_match(model, "D64")
  expect_no_match(model, "E10")
})

test_that("efa_to_lavaan handles single factor", {
  loadings <- matrix(c(0.8, 0.7, 0.6), 3, 1)
  rownames(loadings) <- c("D50", "D64", "E10")
  model <- efa_to_lavaan(loadings, threshold = 0.3)
  expect_match(model, "F1 =~ D50 \\+ D64 \\+ E10")
  expect_no_match(model, "F2")
})

test_that("efa_to_lavaan handles cross-loading traits", {
  loadings <- matrix(c(0.6, 0.5, 0.4, 0.3, 0.4, 0.5, 0.6, 0.7), 4, 2)
  rownames(loadings) <- c("A", "B", "C", "D")
  model <- efa_to_lavaan(loadings, threshold = 0.3)
  expect_match(model, "F1 =~")
  expect_match(model, "F2 =~")
})

test_that("efa_to_lavaan returns empty string when no loadings pass", {
  loadings <- matrix(c(0.1, 0.1, 0.1), 3, 1)
  rownames(loadings) <- c("A", "B", "C")
  model <- efa_to_lavaan(loadings, threshold = 0.3)
  expect_equal(model, "")
})
