source(file.path(project_root, "R/08_gwas.R"))

test_that("parse_factor_names extracts LHS of =~ lines", {
  m <- "Blood =~ D50 + D64\nEndocrine =~ E10 + E11"
  expect_equal(parse_factor_names(m), c("Blood", "Endocrine"))
})

test_that("parse_factor_names trims whitespace and ignores non-loading lines", {
  m <- "  Blood   =~ D50 + D64  \nEndocrine =~ E10\nBlood ~~ Endocrine"
  expect_equal(parse_factor_names(m), c("Blood", "Endocrine"))
})

test_that("parse_factor_names returns empty for empty model string", {
  expect_equal(parse_factor_names(""), character(0))
})

test_that("build_gwas_model appends Factor ~ SNP lines and emits sub vector", {
  m <- "Blood =~ D50 + D64\nEndocrine =~ E10 + E11"
  spec <- build_gwas_model(m)
  expect_equal(spec$factors, c("Blood", "Endocrine"))
  expect_equal(spec$sub, c("Blood ~ SNP", "Endocrine ~ SNP"))
  expect_match(spec$model, "Blood ~ SNP")
  expect_match(spec$model, "Endocrine ~ SNP")
  expect_true(grepl("=~ D50", spec$model))
})

test_that("build_gwas_model errors on empty model", {
  expect_error(build_gwas_model(""), "no factors")
})
