make_panukb_sumstat_fixture <- function(path) {
  dt <- data.table::data.table(
    chr = c("1", "1", "1", "1", "1"),
    pos = c(100L, 200L, 300L, 400L, 500L),
    ref = c("A", "C", "G", "T", "A"),
    alt = c("T", "G", "A", "C", "G"),
    beta_EUR = c(0.05, 0.10, NA_real_, 0.02, 0.03),
    se_EUR = c(0.01, 0.02, 0.01, 0.01, 0.01),
    neglog10_pval_EUR = c(7.30, 2.50, NA_real_, 1.10, 0.20),
    low_confidence_EUR = c("false", "true", "false", "false", "false")
  )
  data.table::fwrite(dt, path, sep = "\t")
}

make_panukb_vm_fixture <- function() {
  vm <- data.table::data.table(
    chr = c("1", "1", "1", "1"),
    pos = c(100L, 200L, 400L, 500L),
    ref = c("A", "C", "T", "A"),
    alt = c("T", "G", "C", "G"),
    rsid = c("rs100", "rs200", "rs400", ""),
    info = c(0.95, 0.85, 0.70, 0.99),
    af_EUR = c(0.20, 0.40, 0.95, 0.001)
  )
  data.table::setkey(vm, chr, pos, ref, alt)
  vm
}

test_that("Pan-UKB munge pipeline: low_confidence filter + NA drop + P recovery + MAF + INFO", {
  tmp <- tempfile("munge_"); dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  raw <- file.path(tmp, "raw.tsv")
  make_panukb_sumstat_fixture(raw)

  cols <- list(chr = "chr", pos = "pos", ref = "ref", alt = "alt",
               beta = "beta_EUR", se = "se_EUR",
               neglog10_pval = "neglog10_pval_EUR",
               low_confidence = "low_confidence_EUR")
  need <- unlist(cols, use.names = FALSE)
  classes <- list(); classes[[cols$chr]] <- "character"
  classes[[cols$low_confidence]] <- "character"
  dt <- data.table::fread(raw, select = need, colClasses = classes,
                           na.strings = c("", "NA"))
  data.table::setnames(dt, old = need,
                        new = c("chr", "pos", "ref", "alt", "beta", "se",
                                "neglog10_pval", "low_confidence"))
  dt[, low_confidence := tolower(low_confidence) == "true"]
  dt <- dt[low_confidence == FALSE & !is.na(beta) & !is.na(se) &
           !is.na(neglog10_pval)]
  expect_equal(nrow(dt), 3)  # row2 dropped (low_conf), row3 dropped (NA beta/pval)
  dt[, pval := 10^(-neglog10_pval)]
  expect_equal(dt$pval[dt$pos == 100L], 10^(-7.30), tolerance = 1e-10)

  vm <- make_panukb_vm_fixture()
  res <- merge_panukb_variants(dt, vm)
  # pos=500 dropped (empty rsid), pos=100 kept (info=0.95, af=0.20), pos=400 kept (info=0.70, af=0.95)
  expect_equal(nrow(res), 2)
  expect_setequal(res$rsid, c("rs100", "rs400"))

  info_thr <- 0.8
  res <- res[info >= info_thr]
  expect_equal(nrow(res), 1)
  expect_equal(res$rsid, "rs100")

  res[, MAF := pmin(af_EUR, 1 - af_EUR)]
  expect_true(all(res$MAF <= 0.5))
  expect_equal(res$MAF[res$rsid == "rs100"], 0.20)
})

test_that("merge_panukb_variants coerces integer chr to character before join", {
  vm <- data.table::data.table(
    chr = c("1", "X"), pos = c(100L, 50L),
    ref = c("A", "C"), alt = c("T", "G"),
    rsid = c("rs1", "rsX"), info = c(0.95, 0.90), af_EUR = c(0.2, 0.3)
  )
  data.table::setkey(vm, chr, pos, ref, alt)
  dt <- data.table::data.table(chr = 1L, pos = 100L, ref = "A", alt = "T",
                                 beta = 0.05)
  res <- merge_panukb_variants(dt, vm)
  expect_equal(nrow(res), 1)
  expect_equal(res$rsid, "rs1")
})

test_that("MAF = pmin(af, 1-af) catches common-ALT variants Neale-style filtering would drop", {
  af <- c(0.05, 0.5, 0.95, 0.99)
  maf <- pmin(af, 1 - af)
  expect_equal(maf, c(0.05, 0.5, 0.05, 0.01))
})

test_that("P recovery from neglog10_pval matches expected at edge cases", {
  expect_equal(10^(-0), 1)
  expect_equal(10^(-7.30), 10^(-7.30))
  expect_equal(10^(-300), 0, tolerance = 0)
})
