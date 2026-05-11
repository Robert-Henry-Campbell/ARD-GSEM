test_that("config resolves relative paths against project root", {
  cfg <- read_config(file.path(project_root, "config/pipeline.yaml"))
  expect_true(startsWith(cfg$paths$sumstats_dir, "/"))
  expect_true(startsWith(cfg$paths$output_dir, "/"))
  expect_true(startsWith(cfg$paths$reference_dir, "/"))
})

test_that("config has expected default values", {
  cfg <- read_config(file.path(project_root, "config/pipeline.yaml"))
  expect_equal(cfg$parallel$n_workers, 24)
  expect_equal(cfg$ldsc$h2_z_threshold, 2.0)
  expect_equal(cfg$munge$maf_threshold, 0.01)
  expect_equal(cfg$efa$rotation, "geomin")
  expect_equal(cfg$cfa$estimator, "DWLS")
})

test_that("config errors on missing file", {
  expect_error(read_config("/nonexistent/path.yaml"))
})

test_that("config paths point to existing directories (where expected)", {
  cfg <- read_config(file.path(project_root, "config/pipeline.yaml"))
  expect_true(dir.exists(cfg$paths$sumstats_dir))
  expect_true(dir.exists(cfg$paths$manifest_dir))
  expect_true(dir.exists(cfg$paths$meta_dir))
})
