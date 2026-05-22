test_that("discover_traits bothsex extracts 3-char codes from icd10-XXX-both_sexes.tsv.bgz", {
  tmp <- tempfile("disc_"); dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  bothsex_dir <- file.path(tmp, "bothsex"); dir.create(bothsex_dir)
  for (code in c("E11", "A02", "Z00")) {
    file.create(file.path(bothsex_dir, sprintf("icd10-%s-both_sexes.tsv.bgz", code)))
    file.create(file.path(bothsex_dir, sprintf("icd10-%s-both_sexes.tsv.bgz.tbi", code)))
  }
  file.create(file.path(bothsex_dir, "readme.txt"))
  cfg <- list(paths = list(sumstats_dir = tmp),
              discover = list(bothsex_filename_pattern = "^icd10-.*-both_sexes\\.tsv\\.bgz$"))
  out <- discover_traits(cfg, "bothsex")
  expect_setequal(out, c("E11", "A02", "Z00"))
})

test_that("discover_traits bothsex ignores tbi index files and unrelated files", {
  tmp <- tempfile("disc_"); dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  bothsex_dir <- file.path(tmp, "bothsex"); dir.create(bothsex_dir)
  file.create(file.path(bothsex_dir, "icd10-E11-both_sexes.tsv.bgz"))
  file.create(file.path(bothsex_dir, "icd10-E11-both_sexes.tsv.bgz.tbi"))
  file.create(file.path(bothsex_dir, "biomarkers-30600-both_sexes-irnt.tsv.bgz"))
  cfg <- list(paths = list(sumstats_dir = tmp), discover = list())
  out <- discover_traits(cfg, "bothsex")
  expect_equal(out, "E11")
})

test_that("discover_traits Neale path unchanged for male/female", {
  tmp <- tempfile("disc_"); dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  male_dir <- file.path(tmp, "male"); dir.create(male_dir)
  file.create(file.path(male_dir, "D50.gwas.imputed_v3.male.tsv.bgz"))
  file.create(file.path(male_dir, "E11.gwas.imputed_v3.male.tsv.bgz"))
  cfg <- list(paths = list(sumstats_dir = tmp), discover = list())
  out <- discover_traits(cfg, "male")
  expect_setequal(out, c("D50", "E11"))
})
