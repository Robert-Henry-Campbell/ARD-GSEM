library(data.table)

test_that("build_apriori_model ordering_table sorts indicators by descending |sort_metric|", {
  traits <- c("E10", "E11", "E14")
  categories <- data.table(
    code = c("E10", "E11", "E14"),
    chapter = c("Endocrine", "Endocrine", "Endocrine")
  )
  ordering <- data.table(
    code = c("E10", "E11", "E14"),
    chapter = c("Endocrine", "Endocrine", "Endocrine"),
    sort_metric = c(0.4, 0.9, -0.7)
  )
  model <- build_apriori_model(traits, categories,
                                ordering_table = ordering,
                                add_factor_covariances = FALSE,
                                add_heywood_constraints = FALSE)
  loading_line <- grep("=~", strsplit(model, "\n")[[1]], value = TRUE)[1]
  expect_match(loading_line, "Endocrine =~ E11 \\+ E14 \\+ E10$")
})

test_that("build_apriori_model with NULL ordering_table falls back to alphabetical", {
  traits <- c("E14", "E10", "E11")
  categories <- data.table(
    code = traits,
    chapter = rep("Endocrine", 3)
  )
  model <- build_apriori_model(traits, categories,
                                add_factor_covariances = FALSE,
                                add_heywood_constraints = FALSE)
  loading_line <- grep("=~", strsplit(model, "\n")[[1]], value = TRUE)[1]
  expect_match(loading_line, "Endocrine =~ E10 \\+ E11 \\+ E14$")
})

test_that("check_chapter_survival reports counts and surviving chapters", {
  retained <- c("D50", "D64", "E10", "E11", "E14", "G35")
  categories <- data.table(
    code = c("D50", "D64", "E10", "E11", "E14", "G35"),
    chapter = c("Blood", "Blood", "Endocrine", "Endocrine", "Endocrine", "Nervous")
  )
  res <- check_chapter_survival(retained, categories,
                                 min_indicators = 3L,
                                 min_factors_required = 1L)
  expect_setequal(res$surviving_chapters, "Endocrine")
  expect_setequal(res$dropped_chapters, c("Blood", "Nervous"))
})

test_that("check_chapter_survival halts when too few factors survive", {
  retained <- c("D50", "D64")
  categories <- data.table(code = c("D50", "D64"), chapter = c("Blood", "Blood"))
  expect_error(check_chapter_survival(retained, categories,
                                       min_indicators = 3L,
                                       min_factors_required = 1L),
               "check_chapter_survival")
})

test_that("verify_genome_build passes when SNPs match expected coordinates", {
  dt <- data.table(
    SNP = c("rs1421085", "rs429358", "rs7903146"),
    chr = c("16", "19", "10"),
    pos = c(53800954L, 45411941L, 114758349L)
  )
  expected <- list(
    list(rsid = "rs1421085", chr = "16", pos = 53800954L),
    list(rsid = "rs429358",  chr = "19", pos = 45411941L)
  )
  expect_silent(verify_genome_build(dt, expected))
})

test_that("verify_genome_build fatals on coordinate mismatch", {
  dt <- data.table(SNP = "rs1421085", chr = "16", pos = 99999999L)
  expected <- list(list(rsid = "rs1421085", chr = "16", pos = 53800954L))
  expect_error(verify_genome_build(dt, expected), "verify_genome_build")
})

test_that("verify_genome_build is a no-op when no SNPs are configured", {
  dt <- data.table(SNP = "rs1", chr = "1", pos = 1L)
  expect_silent(verify_genome_build(dt, list()))
})

test_that("compute_factor_nhat matches Wiki 4 formula on MAF >= threshold subset", {
  maf <- c(0.05, 0.20, 0.50)
  se  <- c(0.01, 0.02, 0.03)
  # default threshold = 0.10 drops the first
  expected <- mean(1 / ((2 * maf[2:3] * (1 - maf[2:3])) * se[2:3]^2))
  expect_equal(compute_factor_nhat(maf, se), expected)
})

test_that("compute_factor_nhat returns NA when no SNPs pass MAF filter", {
  expect_true(is.na(compute_factor_nhat(c(0.01, 0.05), c(0.01, 0.01))))
})

test_that("extract_apriori_loadings shapes a usermodel-style results frame", {
  apriori_fit <- list(results = data.frame(
    lhs = c("Blood", "Blood", "Endocrine"),
    op  = c("=~", "=~", "=~"),
    rhs = c("D50", "D64", "E11"),
    STD_Genotype = c(0.6, 0.4, 0.8),
    stringsAsFactors = FALSE
  ))
  categories <- data.table(
    code = c("D50", "D64", "E11"),
    chapter = c("Blood", "Blood", "Endocrine")
  )
  out <- extract_apriori_loadings(apriori_fit, categories)
  expect_equal(nrow(out), 3L)
  expect_setequal(names(out), c("code", "chapter", "factor", "std_loading", "sort_metric"))
  expect_equal(out[code == "D50"]$std_loading, 0.6)
})

test_that("extract_apriori_loadings returns NULL on empty input", {
  expect_null(extract_apriori_loadings(list(results = NULL), data.table()))
})
