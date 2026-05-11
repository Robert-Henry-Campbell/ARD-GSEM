test_that("log messages include timestamp, level, and stage", {
  msg <- format_log_message("INFO", "munge", "Processing D50 male")
  expect_match(msg, "^\\d{4}-\\d{2}-\\d{2} \\d{2}:\\d{2}:\\d{2}")
  expect_match(msg, "\\[INFO\\]")
  expect_match(msg, "\\[munge\\]")
  expect_match(msg, "Processing D50 male")
})

test_that("log level filtering suppresses DEBUG in INFO mode", {
  expect_true(should_log("WARN", min_level = "INFO"))
  expect_true(should_log("ERROR", min_level = "INFO"))
  expect_true(should_log("INFO", min_level = "INFO"))
  expect_false(should_log("DEBUG", min_level = "INFO"))
})

test_that("log level filtering passes everything in DEBUG mode", {
  expect_true(should_log("DEBUG", min_level = "DEBUG"))
  expect_true(should_log("INFO", min_level = "DEBUG"))
  expect_true(should_log("WARN", min_level = "DEBUG"))
})

test_that("log level filtering only passes FATAL in FATAL mode", {
  expect_false(should_log("ERROR", min_level = "FATAL"))
  expect_true(should_log("FATAL", min_level = "FATAL"))
})

test_that("format_log_message output is single-line", {
  msg <- format_log_message("WARN", "ldsc", "h2_z below threshold")
  expect_false(grepl("\n", msg))
})

test_that("LEVELS ordering is correct", {
  expect_true(LEVELS["DEBUG"] < LEVELS["INFO"])
  expect_true(LEVELS["INFO"] < LEVELS["WARN"])
  expect_true(LEVELS["WARN"] < LEVELS["ERROR"])
  expect_true(LEVELS["ERROR"] < LEVELS["FATAL"])
})
