test_that("a priori model groups traits by ICD-10 chapter (plain core)", {
  traits <- c("D50", "D64", "E10", "E11")
  categories <- data.table(
    code = c("D50", "D64", "E10", "E11"),
    chapter = c("Blood diseases", "Blood diseases", "Endocrine", "Endocrine")
  )
  model <- build_apriori_model(traits, categories,
                                add_factor_covariances = FALSE,
                                add_heywood_constraints = FALSE)
  expect_match(model, "D50")
  expect_match(model, "D64")
  expect_match(model, "E10")
  expect_match(model, "E11")
  lines <- strsplit(model, "\n")[[1]]
  expect_equal(length(lines), 2)
})

test_that("a priori model with defaults adds factor covariances and Heywood constraints", {
  traits <- c("D50", "D64", "E10", "E11")
  categories <- data.table(
    code = c("D50", "D64", "E10", "E11"),
    chapter = c("Blood_diseases", "Blood_diseases", "Endocrine", "Endocrine")
  )
  model <- build_apriori_model(traits, categories)
  expect_match(model, "Blood_diseases =~")
  expect_match(model, "Endocrine =~")
  expect_match(model, "Blood_diseases ~~ Endocrine")
  for (code in traits) {
    expect_match(model, sprintf("%s ~~ rv_%s\\*%s", code, code, code))
    expect_match(model, sprintf("rv_%s > 0\\.001", code))
  }
})

test_that("factors with <2 indicators are dropped", {
  traits <- c("D50", "D64", "E10")
  categories <- data.table(
    code = c("D50", "D64", "E10"),
    chapter = c("Blood diseases", "Blood diseases", "Endocrine")
  )
  model <- build_apriori_model(traits, categories, min_indicators = 2)
  expect_no_match(model, "E10")
  expect_match(model, "D50")
  expect_match(model, "D64")
})

test_that("all single-indicator factors produce empty model", {
  traits <- c("D50", "E10", "G35")
  categories <- data.table(
    code = c("D50", "E10", "G35"),
    chapter = c("Blood", "Endocrine", "Nervous")
  )
  model <- build_apriori_model(traits, categories, min_indicators = 2)
  expect_equal(model, "")
})

test_that("factor names are sanitized for lavaan", {
  traits <- c("D50", "D64")
  categories <- data.table(
    code = c("D50", "D64"),
    chapter = c("Diseases of the blood", "Diseases of the blood")
  )
  model <- build_apriori_model(traits, categories)
  factor_name <- sub(" =~.*", "", model)
  expect_no_match(factor_name, " ")
  expect_match(model, "=~ D50 \\+ D64")
})

test_that("min_indicators = 3 raises bar correctly", {
  traits <- c("D50", "D64", "E10", "E11", "E14")
  categories <- data.table(
    code = c("D50", "D64", "E10", "E11", "E14"),
    chapter = c("Blood", "Blood", "Endocrine", "Endocrine", "Endocrine")
  )
  model <- build_apriori_model(traits, categories, min_indicators = 3)
  expect_match(model, "E10")
  expect_no_match(model, "D50")
})
