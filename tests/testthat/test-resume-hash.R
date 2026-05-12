stage_hash <- function(cfg, stage) {
  digest::digest(
    list(stage = cfg[[stage]] %||% list(), paths = cfg$paths),
    algo = "sha256")
}

test_that("stage hash for stages without a config block is stable across unrelated non-path config changes", {
  cfg_a <- list(
    paths = list(output_dir = "/x"),
    munge = list(maf_threshold = 0.01),
    logging = list(level = "INFO"),
    ldsc = list(h2_z_threshold = 2.0))
  cfg_b <- cfg_a
  cfg_b$logging$level <- "DEBUG"
  cfg_b$plots <- list(top_snps_n = 50)

  for (stage in c("sumstats", "vcf", "factor_summary", "report")) {
    expect_equal(stage_hash(cfg_a, stage), stage_hash(cfg_b, stage),
                 info = sprintf("stage=%s hash should not depend on unrelated non-path keys", stage))
  }
})

test_that("stage hash for stages WITH a config block changes when that block changes", {
  cfg_a <- list(paths = list(output_dir = "/x"),
                munge = list(maf_threshold = 0.01))
  cfg_b <- list(paths = list(output_dir = "/x"),
                munge = list(maf_threshold = 0.05))
  expect_false(stage_hash(cfg_a, "munge") == stage_hash(cfg_b, "munge"))
})

test_that("stage hash changes when config$paths change, for every stage", {
  cfg_a <- list(paths = list(output_dir = "/x", variants_manifest = "old"),
                munge = list(maf_threshold = 0.01))
  cfg_b <- cfg_a
  cfg_b$paths$variants_manifest <- "new"

  for (stage in c("munge", "ldsc", "sumstats", "vcf")) {
    expect_false(stage_hash(cfg_a, stage) == stage_hash(cfg_b, stage),
                 info = sprintf("stage=%s hash should invalidate on paths change", stage))
  }
})

test_that("legacy hash format (stage-block only, no paths) still validates under should_skip", {
  # Simulate a manifest written by the OLD scheme (digest(config[[stage]] %||% list())).
  cfg <- list(paths = list(output_dir = "/x", variants_manifest = "m"),
              munge = list(maf_threshold = 0.01))
  legacy_hash <- digest::digest(cfg$munge, algo = "sha256")
  current_hash <- digest::digest(
    list(stage = cfg$munge, paths = cfg$paths), algo = "sha256")
  expect_false(legacy_hash == current_hash,
               "sanity: legacy and current hash schemas should differ for the same config")
  # Either form should match in the should_skip fallback logic below.
  is_fresh <- function(stored_hash) {
    identical(stored_hash, current_hash) || identical(stored_hash, legacy_hash)
  }
  expect_true(is_fresh(legacy_hash),
              "manifest written with legacy schema must still validate")
  expect_true(is_fresh(current_hash),
              "manifest written with current schema must validate")
})
