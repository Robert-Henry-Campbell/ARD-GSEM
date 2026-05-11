test_that("build_apriori_model asserts 3-character trait codes", {
  categories <- data.table(code = c("D50", "D64"), chapter = c("Blood", "Blood"))
  expect_error(build_apriori_model(c("D50", "D6"), categories), "nchar")
  expect_error(build_apriori_model(c("D50", "D6400"), categories), "nchar")
})

test_that("no first-letter fallback: unmapped D50 does NOT fall into a 'D' chapter", {
  categories <- data.table(
    code = c("C00", "C01", "C02", "E10", "E11"),
    chapter = c("Neoplasms", "Neoplasms", "Neoplasms", "Endocrine", "Endocrine")
  )
  model <- build_apriori_model(c("D50", "D64", "E10", "E11"), categories)
  expect_no_match(model, "Neoplasms")
  expect_no_match(model, "D50")
  expect_no_match(model, "D64")
  expect_match(model, "E10")
  expect_match(model, "E11")
})

test_that("all-unmapped traits produce empty model", {
  categories <- data.table(code = c("C00"), chapter = c("Neoplasms"))
  expect_equal(build_apriori_model(c("D50", "D64"), categories), "")
})

test_that("duplicated factor names after 20-char truncation are rejected", {
  categories <- data.table(
    code = c("D50", "D64", "E10", "E11"),
    chapter = c(
      "AAAAAAAAAAAAAAAAAAAA_long_one",
      "AAAAAAAAAAAAAAAAAAAA_long_one",
      "AAAAAAAAAAAAAAAAAAAA_long_two",
      "AAAAAAAAAAAAAAAAAAAA_long_two"
    )
  )
  expect_error(build_apriori_model(c("D50", "D64", "E10", "E11"), categories),
               "duplicated")
})
