test_that("get_case_control dispatches to panukb_bothsex_manifest.rda for sex=bothsex", {
  tmp <- tempfile("cc_"); dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  manifest_dir <- file.path(tmp, "manifest"); dir.create(manifest_dir)
  panukb_bothsex_manifest <- data.frame(
    phenotype = c("E11", "A02"),
    n_cases = c(20000, 112),
    n_controls = c(400000, 420419),
    stringsAsFactors = FALSE
  )
  save(panukb_bothsex_manifest,
       file = file.path(manifest_dir, "panukb_bothsex_manifest.rda"))
  cfg <- list(paths = list(manifest_dir = manifest_dir))
  cc <- get_case_control(cfg, "bothsex", "E11")
  expect_equal(cc$n_cases, 20000)
  expect_equal(cc$n_controls, 400000)
})

test_that("get_case_control still reads Neale manifest for sex=male", {
  tmp <- tempfile("cc_"); dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  manifest_dir <- file.path(tmp, "manifest"); dir.create(manifest_dir)
  neale_male_manifest <- data.frame(
    phenotype = c("D50", "E11"),
    n_cases = c(5000, 30000),
    n_controls = c(45000, 370000),
    stringsAsFactors = FALSE
  )
  save(neale_male_manifest,
       file = file.path(manifest_dir, "neale_male_manifest.rda"))
  cfg <- list(paths = list(manifest_dir = manifest_dir))
  cc <- get_case_control(cfg, "male", "D50")
  expect_equal(cc$n_cases, 5000)
})

test_that("resolve_column_map: NULL user_map returns defaults unchanged", {
  defaults <- list(chr = "chr", beta = "beta_EUR", se = "se_EUR")
  out <- resolve_column_map(NULL, defaults, label = "test")
  expect_equal(out, defaults)
})

test_that("resolve_column_map: partial override fills missing keys from defaults", {
  defaults <- list(chr = "chr", pos = "pos", beta = "beta_EUR",
                   se = "se_EUR", neglog10_pval = "neglog10_pval_EUR")
  user <- list(beta = "beta_meta")
  out <- resolve_column_map(user, defaults, label = "test")
  expect_equal(out$beta, "beta_meta")
  expect_equal(out$chr, "chr")
  expect_equal(out$pos, "pos")
  expect_equal(out$se, "se_EUR")
  expect_equal(out$neglog10_pval, "neglog10_pval_EUR")
})

test_that("resolve_column_map: NULL/empty values in user map ignored, defaults retained", {
  defaults <- list(chr = "chr", beta = "beta_EUR")
  user <- list(chr = NULL, beta = "")
  out <- resolve_column_map(user, defaults, label = "test")
  expect_equal(out$chr, "chr")
  expect_equal(out$beta, "beta_EUR")
})

test_that("resolve_column_map: reordered user keys do not affect result mapping", {
  defaults <- list(chr = "chrom", pos = "pos", af = "af_EUR")
  user <- list(af = "af_AFR", chr = "chromosome")
  out <- resolve_column_map(user, defaults, label = "test")
  expect_equal(out$chr, "chromosome")
  expect_equal(out$pos, "pos")
  expect_equal(out$af, "af_AFR")
})

test_that("merge_panukb_variants performs inner join on (chr,pos,ref,alt) and drops rsid-less rows", {
  vm <- data.table::data.table(
    chr = c("1", "1", "1"),
    pos = c(100L, 200L, 300L),
    ref = c("A", "C", "G"),
    alt = c("T", "G", "A"),
    rsid = c("rs1", NA_character_, "rs3"),
    info = c(0.95, 0.85, 0.91),
    af_EUR = c(0.20, 0.40, 0.55)
  )
  data.table::setkey(vm, chr, pos, ref, alt)
  dt <- data.table::data.table(
    chr = c(1L, 1L, 1L, 1L),
    pos = c(100L, 200L, 300L, 999L),
    ref = c("A", "C", "G", "T"),
    alt = c("T", "G", "A", "C"),
    beta = c(0.01, 0.02, 0.03, 0.04)
  )
  res <- merge_panukb_variants(dt, vm)
  expect_equal(nrow(res), 2)
  expect_setequal(res$rsid, c("rs1", "rs3"))
  expect_true("af_EUR" %in% names(res))
  expect_true("info" %in% names(res))
  expect_true("beta" %in% names(res))
})
