# scripts/run_pipeline.R auto-runs parse_args() at top level so we can't source it.
# Mirror the parse_sexes helper here verbatim; the orchestrator script is the
# single source of truth and a divergence is caught by smoke runs.
parse_sexes_local <- function(sex_arg) {
  valid <- c("male", "female", "bothsex", "bothsex_meta")
  sexes <- trimws(strsplit(sex_arg, ",", fixed = TRUE)[[1L]])
  sexes <- sexes[nzchar(sexes)]
  if (length(sexes) == 0L) stop("--sex parsed to an empty list", call. = FALSE)
  unknown <- setdiff(sexes, valid)
  if (length(unknown) > 0L) {
    stop(sprintf("Unknown --sex value(s): %s (valid: %s)",
                 paste(unknown, collapse = ","),
                 paste(valid, collapse = ",")), call. = FALSE)
  }
  if (anyDuplicated(sexes)) stop(sprintf("--sex contains duplicates: %s",
                                          paste(sexes, collapse = ",")),
                                  call. = FALSE)
  sexes
}

test_that("parse_sexes accepts comma-separated subset of {male,female,bothsex,bothsex_meta}", {
  expect_equal(parse_sexes_local("male"), "male")
  expect_equal(parse_sexes_local("female"), "female")
  expect_equal(parse_sexes_local("bothsex"), "bothsex")
  expect_equal(parse_sexes_local("bothsex_meta"), "bothsex_meta")
  expect_equal(parse_sexes_local("male,female"), c("male", "female"))
  expect_equal(parse_sexes_local("male,female,bothsex"), c("male", "female", "bothsex"))
  expect_equal(parse_sexes_local("bothsex,bothsex_meta"), c("bothsex", "bothsex_meta"))
})

test_that("parse_sexes trims whitespace inside list", {
  expect_equal(parse_sexes_local(" male , bothsex "), c("male", "bothsex"))
})

test_that("parse_sexes rejects the legacy 'both' alias", {
  expect_error(parse_sexes_local("both"), "Unknown --sex value")
})

test_that("parse_sexes rejects unknown values", {
  expect_error(parse_sexes_local("foobar"), "Unknown --sex value")
  expect_error(parse_sexes_local("male,foobar"), "foobar")
})

test_that("parse_sexes rejects duplicates", {
  expect_error(parse_sexes_local("male,male"), "duplicates")
  expect_error(parse_sexes_local("male,bothsex,bothsex"), "duplicates")
})
