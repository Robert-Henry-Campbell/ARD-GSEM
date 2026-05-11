test_that("manifest hash for stages without a config block is stable across unrelated config changes", {
  cfg_a <- list(munge = list(maf_threshold = 0.01),
                logging = list(level = "INFO"),
                ldsc = list(h2_z_threshold = 2.0))
  cfg_b <- cfg_a
  cfg_b$logging$level <- "DEBUG"
  cfg_b$plots <- list(top_snps_n = 50)

  for (stage in c("sumstats", "vcf", "factor_summary", "report")) {
    hash_a <- digest::digest(cfg_a[[stage]] %||% list(), algo = "sha256")
    hash_b <- digest::digest(cfg_b[[stage]] %||% list(), algo = "sha256")
    expect_equal(hash_a, hash_b,
                 info = sprintf("stage=%s hash should not depend on unrelated config keys", stage))
  }
})

test_that("manifest hash for stages WITH a config block changes when that block changes", {
  cfg_a <- list(munge = list(maf_threshold = 0.01))
  cfg_b <- list(munge = list(maf_threshold = 0.05))
  hash_a <- digest::digest(cfg_a[["munge"]] %||% list(), algo = "sha256")
  hash_b <- digest::digest(cfg_b[["munge"]] %||% list(), algo = "sha256")
  expect_false(hash_a == hash_b)
})
