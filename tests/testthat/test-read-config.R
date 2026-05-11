make_tmp_yaml <- function(root_value = NULL) {
  yaml_path <- tempfile(fileext = ".yaml")
  on.exit(NULL)
  root_line <- if (is.null(root_value)) "  root: null" else sprintf("  root: %s", root_value)
  writeLines(c(
    "project:",
    "  name: test",
    root_line,
    "paths:",
    "  sumstats_dir: sumstats",
    "  reference_dir: reference",
    "  ld_scores: reference/eur_w_ld_chr"
  ), yaml_path)
  yaml_path
}

test_that("read_config absolutises relative paths against the supplied root", {
  root <- tempfile("proj_")
  dir.create(file.path(root, "config"), recursive = TRUE)
  yaml_path <- file.path(root, "config", "pipeline.yaml")
  file.copy(make_tmp_yaml("/nonexistent/path"), yaml_path, overwrite = TRUE)

  cfg <- read_config(yaml_path, root = root)
  expect_equal(cfg$project$root, normalizePath(root, winslash = "/", mustWork = TRUE))
  expect_match(cfg$paths$sumstats_dir, "sumstats$")
  expect_true(startsWith(cfg$paths$sumstats_dir, normalizePath(root, winslash = "/")))
})

test_that("read_config falls back to the YAML directory when project.root is unusable", {
  root <- tempfile("proj_")
  dir.create(file.path(root, "config"), recursive = TRUE)
  yaml_path <- file.path(root, "config", "pipeline.yaml")
  file.copy(make_tmp_yaml("/nonexistent/path"), yaml_path, overwrite = TRUE)
  cfg <- read_config(yaml_path)
  expect_equal(cfg$project$root, normalizePath(root, winslash = "/", mustWork = TRUE))
})

test_that("read_config leaves absolute paths unchanged", {
  root <- tempfile("proj_")
  dir.create(file.path(root, "config"), recursive = TRUE)
  yaml_path <- file.path(root, "config", "pipeline.yaml")
  writeLines(c("project:", "  root: null",
               "paths:", "  ld_scores: /abs/ld",
               "  hm3_snplist: C:/abs/hm3"), yaml_path)
  cfg <- read_config(yaml_path, root = root)
  expect_equal(cfg$paths$ld_scores, "/abs/ld")
  expect_equal(cfg$paths$hm3_snplist, "C:/abs/hm3")
})

test_that("is_absolute_path recognises POSIX, Windows, and UNC paths", {
  expect_true(is_absolute_path("/foo/bar"))
  expect_true(is_absolute_path("C:/Users"))
  expect_true(is_absolute_path("D:\\Users"))
  expect_true(is_absolute_path("\\\\server\\share"))
  expect_false(is_absolute_path("foo/bar"))
  expect_false(is_absolute_path("./foo"))
})
