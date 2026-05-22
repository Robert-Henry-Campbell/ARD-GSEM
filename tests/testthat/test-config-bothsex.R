test_that("config/pipeline.yaml loads new Pan-UKB keys", {
  cfg <- read_config(file.path(project_root, "config/pipeline.yaml"),
                     root = project_root)
  expect_true(!is.null(cfg$panukb))
  expect_equal(cfg$panukb$ancestry, "EUR")
  # cfg$panukb is not path-resolved by read_config (only cfg$paths is); the
  # value is a relative path resolved against CWD by callers.
  expect_equal(cfg$panukb$phenotype_manifest,
               "manifest/panukb_phenotype_manifest.csv")
  expect_equal(cfg$paths$panukb_variants_manifest,
               file.path(project_root, "manifest/full_variant_qc_metrics.txt.bgz"))
  expect_equal(cfg$panukb$columns$beta, "beta_EUR")
  expect_equal(cfg$panukb$columns$neglog10_pval, "neglog10_pval_EUR")
  expect_equal(cfg$panukb$columns$low_confidence, "low_confidence_EUR")
  expect_equal(cfg$panukb$variants_columns$chr, "chrom")
  expect_equal(cfg$panukb$variants_columns$af, "af_EUR")
})

test_that("config has info_threshold (munge) and info_filter (sumstats)", {
  cfg <- read_config(file.path(project_root, "config/pipeline.yaml"),
                     root = project_root)
  expect_equal(cfg$munge$info_threshold, 0.8)
  expect_equal(cfg$sumstats$info_filter, 0.6)
})

test_that("config has bothsex_filename_pattern under discover", {
  cfg <- read_config(file.path(project_root, "config/pipeline.yaml"),
                     root = project_root)
  expect_equal(cfg$discover$bothsex_filename_pattern,
               "^icd10-.*-both_sexes\\.tsv\\.bgz$")
})
