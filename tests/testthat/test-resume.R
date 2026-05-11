test_that("stage is skipped when config hash matches", {
  manifest <- list(stage = "munge", config_hash = "abc123", output_hash = "def456")
  expect_true(stage_is_fresh(manifest, "abc123"))
})

test_that("stage re-runs when config changes", {
  manifest <- list(stage = "munge", config_hash = "abc123", output_hash = "def456")
  expect_false(stage_is_fresh(manifest, "CHANGED"))
})

test_that("NULL manifest means stage needs to run", {
  expect_false(stage_is_fresh(NULL, "abc123"))
})

test_that("empty config hash doesn't match", {
  manifest <- list(stage = "munge", config_hash = "abc123")
  expect_false(stage_is_fresh(manifest, ""))
})
