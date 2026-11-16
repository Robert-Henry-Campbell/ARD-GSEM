for (f in c("R/15_postgwas_common.R", "R/32_figures.R")) source(file.path(project_root, f), local = TRUE)

test_that("theme_nature is a 7-pt theme with bold panel tags", {
  th <- theme_nature()
  expect_s3_class(th, "theme")
  expect_equal(th$text$size, 7)
  expect_identical(th$plot.tag$face, "bold")
  expect_equal(theme_nature(grid = "y")$panel.grid.major.y$linewidth, 0.25)
})

test_that("save_nature_fig writes PDF, PNG and source data, and warns on non-standard widths", {
  skip_if_not(capabilities("cairo"))
  tmp <- tempfile(); dir.create(tmp)
  p <- with_source_data(ggplot2::ggplot(data.frame(x = 1:3, y = 1:3), ggplot2::aes(x, y)) + ggplot2::geom_point() + theme_nature(),
                        data.frame(x = 1:3, y = 1:3))
  files <- save_nature_fig(p, "t", 89, 60, tmp, dpi = 72)
  expect_true(all(file.exists(files)))
  expect_setequal(basename(files), c("t.pdf", "t.png", "t__source_data.csv"))
  expect_equal(nrow(fread(files[3])), 3)
  invisible(capture.output(f2 <- save_nature_fig(p, "u", 100, 300, tmp, dpi = 72, source_data = list(a = data.frame(z = 1), b = data.frame(z = 2))), type = "message"))
  expect_length(f2, 4)
  expect_true(file.exists(file.path(tmp, "u__source_data_a.csv")) && file.exists(file.path(tmp, "u__source_data_b.csv")))
})

test_that("assemble_figure skips NULL panels and reports errors without stopping", {
  fig_ok <- function(config, sex) with_source_data(ggplot2::ggplot(data.frame(x = 1, y = 1), ggplot2::aes(x, y)) + ggplot2::geom_point(), data.frame(x = 1))
  fig_none <- function(config, sex) NULL
  fig_bad <- function(config, sex) stop("boom")
  environment(assemble_figure) <- environment()
  a <- assemble_figure(list(name = "f", panels = c("fig_ok", "fig_none", "fig_bad", "fig_missing"), ncol = 2), list(), "bothsex")
  expect_identical(a$status$fig_ok, "ok")
  expect_identical(a$status$fig_none, "skipped: inputs missing")
  expect_match(a$status$fig_bad, "^error: boom")
  expect_match(a$status$fig_missing, "no such panel")
  expect_true(inherits(a$plot, "ggplot"))
  expect_identical(names(a$source_data), "a_ok")
})

test_that("run_figures writes an index and a manifest even when every panel is skipped", {
  tmp <- tempfile()
  cfg <- list(paths = list(output_dir = tmp), figures = list(main = list(list(name = "fig1", panels = c("fig_nothing"), width_mm = 89, height_mm = 60))))
  fig_nothing <- function(config, sex) NULL
  environment(run_figures) <- environment(); environment(assemble_figure) <- environment()
  idx <- run_figures(cfg, "bothsex")
  expect_identical(idx$status, "skipped: inputs missing")
  expect_true(file.exists(file.path(tmp, "bothsex", "figures", "figure_index.tsv")))
  expect_true(file.exists(file.path(tmp, "bothsex", "figures", ".manifest.json")))
})
