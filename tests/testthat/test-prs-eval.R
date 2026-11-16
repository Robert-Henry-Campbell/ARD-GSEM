for (f in c("R/15_postgwas_common.R", "R/29_prs_external.R", "R/31_prs_eval_pheno.R", "R/31_prs_eval_stats.R",
            "R/31_prs_eval_finngen.R", "R/31_prs_eval.R")) {
  source(file.path(project_root, f), local = TRUE)
}

test_that("score matrices align weights to the target ALT allele and record coverage", {
  tmp <- tempfile(); dir.create(tmp)
  target <- data.table(SNP = c("rs1", "rs2", "rs3", "rs4"), CHR = c(1L, 1L, 2L, 2L), BP = 1:4, REF = c("A", "C", "G", "T"), ALT = c("G", "T", "A", "C"))
  sets <- list(list(prs_id = "F__ct__pT_1", kind = "factor", name = "F", factor = "F", method = "ct", threshold = "pT_1",
                    dt = data.table(SNP = c("rs1", "rs2", "rs3", "rs9"), A1 = c("G", "C", "C", "A"), A2 = c("A", "T", "A", "G"), weight = c(0.5, 0.2, 0.1, 0.3))),
               list(prs_id = "C__prscs", kind = "comparator", name = "C", factor = "ARD", method = "prscs", threshold = NA_character_,
                    dt = data.table(SNP = c("rs2", "rs4"), A1 = c("T", "T"), A2 = c("C", "C"), weight = c(-0.4, 0.6))))
  cov <- pe_build_score_matrix(sets, target, tmp)
  expect_identical(cov$n_aligned, c(2L, 2L))                        # rs3 allele mismatch, rs9 absent
  expect_identical(cov$n_sign_flipped, c(1L, 1L))
  m1 <- fread(file.path(tmp, "weights", "score_matrix", "chr1.tsv"))
  expect_identical(names(m1), c("SNP", "A1", "F__ct__pT_1", "C__prscs"))
  expect_equal(m1[SNP == "rs1", F__ct__pT_1], 0.5)                  # A1 was ALT: unchanged
  expect_equal(m1[SNP == "rs2", F__ct__pT_1], -0.2)                 # A1 was REF: flipped, now per ALT
  expect_equal(m1[SNP == "rs2", C__prscs], -0.4)
  expect_identical(m1$A1, c("G", "T"))
  m2 <- fread(file.path(tmp, "weights", "score_matrix", "chr2.tsv"))
  expect_equal(m2[SNP == "rs4", C__prscs], -0.6)
  expect_equal(cov$column_index, c(3L, 4L))
})

test_that("weight sets are loaded per method, with C+T expanded per threshold and comparators tagged", {
  tmp <- tempfile()
  cfg <- list(paths = list(output_dir = tmp), postgwas = list(factors = c("ARD", "UpperGI")),
              prs = list(external_manifest = ""), prs_eval = list(weights = list(methods = c("ct", "prscs"))))
  wd <- file.path(tmp, "bothsex", "prs", "weights"); dir.create(file.path(wd, "ct"), recursive = TRUE); dir.create(file.path(wd, "prscs"))
  fwrite(data.table(SNP = c("rs1", "rs2", "rs3"), A1 = "A", A2 = "G", weight = c(0.1, 0.2, 0.3), P = c(1e-9, 1e-5, 0.5), n_thresholds_passed = 1L),
         file.path(wd, "ct", "ARD.weights.tsv"), sep = "\t")
  writeLines(c("pT_5e-08 0 5e-08", "pT_1e-04 0 1e-04", "pT_1e+00 0 1"), file.path(wd, "ct", "ARD.ranges"))
  fwrite(data.table(SNP = c("rs1", "rs2"), A1 = "A", A2 = "G", weight = c(0.05, 0.02)), file.path(wd, "prscs", "UpperGI.weights.tsv"), sep = "\t")
  s <- pe_load_weight_sets(cfg, "bothsex", cfg$prs_eval)
  expect_length(s, 4)
  ids <- vapply(s, `[[`, character(1), "prs_id")
  expect_setequal(ids, c("ARD__ct__pT_5e-08", "ARD__ct__pT_1e-04", "ARD__ct__pT_1e+00", "UpperGI__prscs"))
  expect_identical(nrow(s[[which(ids == "ARD__ct__pT_1e-04")]]$dt), 2L)
  expect_identical(s[[which(ids == "UpperGI__prscs")]]$factor, "UpperGI")
  expect_true(all(vapply(s, `[[`, character(1), "kind") == "factor"))
})

test_that("sample QC: group assignment, QC chain, KING and GCTA lists, within-group PCs (plink2 mocked)", {
  tmp <- tempfile(); dir.create(tmp)
  pgen <- file.path(tmp, "pgen"); dir.create(pgen)
  ids <- sprintf("%d", 1001:1012)
  fwrite(data.table(`#FID` = ids, IID = ids, SEX = "NA"), file.path(pgen, "UKB_chr22.psam"), sep = "\t")
  writeLines(c("1001 1001", "1002 1002", "1003 1003", "1004 1004"), file.path(tmp, "sa_unrelated.txt"))
  cfg <- list(paths = list(output_dir = tmp),
              refs = list(ukb_asian_black_pgen_dir = pgen, ukb_unrelated_south_asian = file.path(tmp, "sa_unrelated.txt")),
              prs_eval = list(groups = list(SouthAsian = list(codes = c(3, 3001:3004)), Black = list(codes = c(4, 4001:4003))),
                              relatedness = list(primary = "king", king_cutoff = 0.0442), pcs = list(n = 10)))
  ph <- data.table(eid = ids, eth = c(rep("3001.0", 8), rep("4001.0", 4)), sex = c("0", "1", "0", "1", "1", "0", "1", "0", "1", "0", "1", "0"),
                   genetic_sex = c("0.0", "1.0", "1.0", "1.0", "1.0", "0.0", "", "0.0", "1.0", "0.0", "1.0", "0.0"),
                   aneuploidy = c("", "", "", "", "1", "", "", "", "", "", "", ""), het_miss = "")
  # 1003: sex mismatch; 1005: aneuploidy; 1007: no genetic QC -> 5 South Asian pass basic QC; all 4 Black pass
  calls <- list()
  pe_plink2 <- function(config, args, stage, log_file, expect = character(0)) {
    calls[[length(calls) + 1]] <<- args
    out <- args[which(args == "--out") + 1]
    keep <- fread(args[which(args == "--keep") + 1], header = FALSE, colClasses = "character")[[2]]
    if ("--indep-pairwise" %in% args) writeLines(paste0("rs", 1:50), paste0(out, ".prune.in"))
    if ("--king-cutoff" %in% args) fwrite(data.table(`#FID` = keep[-1], IID = keep[-1]), paste0(out, ".king.cutoff.in.id"), sep = "\t")
    if ("--pca" %in% args) fwrite(cbind(data.table(`#FID` = keep, IID = keep), as.data.table(matrix(rnorm(length(keep) * 10), ncol = 10, dimnames = list(NULL, paste0("PC", 1:10))))),
                                  paste0(out, ".eigenvec"), sep = "\t")
    invisible(NULL)
  }
  environment(pe_sample_qc) <- environment()
  qc <- pe_sample_qc(cfg, ph, tmp)
  expect_identical(names(qc), c("SouthAsian", "Black"))
  expect_identical(qc$SouthAsian$steps$in_group, 8L); expect_identical(qc$SouthAsian$steps$basic_keep, 5L)
  expect_identical(qc$SouthAsian$steps$unrelated_king, 4L)                 # mock drops the first of the 5
  expect_identical(qc$SouthAsian$steps$unrelated_gcta_list, 3L)            # 1001, 1002, 1004 of the list pass basic QC
  expect_identical(qc$Black$steps$basic_keep, 4L); expect_identical(qc$Black$steps$unrelated_gcta_list, 0L)
  expect_identical(ncol(qc$SouthAsian$pcs), 11L); expect_identical(nrow(qc$SouthAsian$pcs), 4L)
  kf <- file.path(tmp, "sample_qc", "SouthAsian.keep.txt")
  expect_true(file.exists(kf)); expect_identical(as.character(file.mode(kf)), "600")
  expect_true(file.exists(file.path(tmp, "sample_qc", "sample_qc.json")))
  expect_true(all(vapply(calls, function(a) startsWith(normalizePath(dirname(a[which(a == "--out") + 1]), mustWork = FALSE), normalizePath(tmp)), logical(1))))
  expect_error(pe_plink2_real <- NULL, NA)
})

test_that("pe_plink2 refuses an --out outside the output directory", {
  cfg <- list(paths = list(output_dir = tempfile()), prs_eval = list(), tools = list(plink2 = ""))
  expect_error(pe_plink2(cfg, c("--pfile", "x", "--out", "/mnt/sdh/upamanyu/UKB/forbidden"), "prs_eval", tempfile()), "outside output")
})

test_that("FinnGen summary-level R2 recovers a simulated score-trait covariance", {
  set.seed(3)
  n_snp <- 200
  w <- data.table(SNP = paste0("rs", 1:n_snp), CHR = rep(1:4, each = 50), A1 = "A", A2 = "G", weight = rnorm(n_snp, 0, 0.05))
  p <- runif(n_snp, 0.1, 0.9)
  # independent SNPs: var(S) = sum w^2 2p(1-p); a liability with R2 = 0.02 from the score gives
  # marginal betas beta_j = w_j * sqrt(R2 / var(S)) on the (unit-variance) liability scale
  var_s <- sum(w$weight^2 * 2 * p * (1 - p))
  K <- 0.1
  beta_liab <- w$weight * sqrt(0.02 / var_s)
  # logistic beta ~ liability beta * (K(1-K))^-1 * dnorm(qnorm(1-K)) is the inverse of the Lee factor; feed the implied log-odds
  z <- stats::dnorm(stats::qnorm(1 - K))
  m <- data.table(SNP = w$SNP, CHR = w$CHR, weight = w$weight, beta_a1 = beta_liab * z / (K * (1 - K)), se = 0.01, p_a1 = p)
  r <- finngen_r2(m, var_s, K)
  expect_lt(abs(r$r2_liability - 0.02), 0.002)
  expect_true(r$ci_low <= r$r2_liability && r$r2_liability <= r$ci_high)
  expect_equal(r$n_snps, n_snp)
  # alignment: beta oriented to A1, multi-rsid cells split, allele mismatches dropped
  tmp <- tempfile(fileext = ".gz")
  fwrite(data.table(`#chrom` = "1", rsids = c("rs1", "rs2,rs3", "rs4"), ref = c("G", "A", "C"), alt = c("A", "G", "T"),
                    beta = c(0.1, -0.2, 0.3), sebeta = 0.01, af_alt = c(0.3, 0.6, 0.5)), tmp, sep = "\t")
  a <- finngen_align(finngen_read(tmp), w[1:4, .(SNP, A1, A2, weight)])   # pipeline weight tables have no CHR
  expect_identical(finngen_align(tmp, w[1:4, .(SNP, A1, A2, weight)]), a) # a path is read on the fly
  expect_identical(a$CHR, c(1L, 1L, 1L))
  expect_identical(a$SNP, c("rs1", "rs2", "rs3"))
  expect_equal(a$beta_a1, c(0.1, 0.2, 0.2))                           # rs2/rs3: A1 = ref -> sign flipped
  expect_equal(a$p_a1, c(0.3, 0.4, 0.4))
})

test_that("ICD-10 arrays parse position by position, with sentinel and unequal-length handling", {
  x <- pe_parse_icd(eid = c(1, 2, 3, 4),
                    codes = c("D125|E880|I10", "", "K219|M1999", "E119|I25"),
                    dates = c("2017-11-24|2021-06-25|2017-11-24", "", "1900-01-01|2008-08-25", "2010-01-01"))
  expect_equal(nrow(x), 7)
  expect_identical(x[eid == "1", code3], c("D12", "E88", "I10"))
  expect_identical(as.character(x[eid == "1", date]), c("2017-11-24", "2021-06-25", "2017-11-24"))
  expect_true(is.na(x[eid == "3" & code == "K219", date]))            # sentinel date -> NA
  expect_identical(as.character(x[eid == "3" & code == "M1999", date]), "2008-08-25")
  expect_true(all(is.na(x[eid == "4", date])))                        # unequal lengths -> NA dates
  expect_false("2" %in% x$eid)
})

test_that("per-person counts: ICD three-character unit, GBD roll-up, injuries, mortality follow-up", {
  ph <- data.table(eid = c("a", "b", "c"), sex = c("0", "1", "1"), age = c("55.0", "60", "48"),
                   assess = c("2010-01-01", "2008-06-15", "2009-03-03"),
                   death_date = c(NA, "2016-06-30", NA),
                   icd = c("I201|I251|K210|Z950|S720", "", "E119|E112|J459|T814"),
                   icd_dates = c("2005-01-01|2012-01-01|2015-01-01|2015-01-01|2011-05-05", "", "2000-01-01|2001-01-01|2020-01-01|2020-01-01"))
  ind <- list(ARD = c("I20", "I25", "K21", "Z95", "E11", "J45"), CM = c("I20", "I25", "Z95", "E11"), UGI = "K21")
  gbd <- data.table(icd10 = c("I20", "I25", "K21", "Z95", "E11", "J45"),
                    gbd_cause = c("IHD", "IHD", "Upper digestive", "Adverse effects", "Diabetes", "Asthma"))
  p <- pe_build_pheno(ph, ind, gbd)
  expect_identical(p[eid == "a", count_ARD], 4L); expect_identical(p[eid == "a", count_ARD_gbd], 3L)
  expect_identical(p[eid == "a", count_CM], 3L);  expect_identical(p[eid == "a", count_CM_gbd], 2L)
  expect_identical(p[eid == "a", count_UGI], 1L); expect_identical(p[eid == "a", injury_count], 1L)
  expect_identical(p[eid == "b", c(count_ARD, count_ARD_gbd, injury_count)], c(0L, 0L, 0L))
  expect_identical(p[eid == "c", count_ARD], 2L)                      # E119 and E112 are one code
  expect_identical(p[eid == "c", injury_count], 0L)                   # chapter T (T814, complication of care) is not an injury
  expect_identical(pe_build_pheno(ph, ind, gbd, injury_chapters = c("S", "T"))[eid == "c", injury_count], 1L)
  expect_identical(p$n_hes_dates, c(4L, 0L, 3L))                      # distinct first-diagnosis dates (hospital-contact proxy)
  expect_identical(p$death, c(0L, 1L, 0L))
  expect_equal(p[eid == "b", t_death], as.numeric(as.Date("2016-06-30") - as.Date("2008-06-15")) / 365.25)
  expect_equal(p[eid == "a", t_death], as.numeric(as.Date("2022-12-19") - as.Date("2010-01-01")) / 365.25)
  expect_equal(p[eid == "b", age_exit], 60 + p[eid == "b", t_death])
  expect_identical(p$any_hes, c(TRUE, FALSE, TRUE))
})

test_that("ethnicity codes map to groups", {
  g <- list(SouthAsian = list(codes = c(3, 3001, 3002, 3003, 3004)), Black = list(codes = c(4, 4001, 4002, 4003)))
  expect_identical(pe_assign_group(c("3001.0", "4", "1001.0", "", "4002.0"), g), c("SouthAsian", "Black", NA, NA, "Black"))
})

test_that("negative binomial recovers a simulated IRR per SD", {
  skip_if_not_installed("MASS")
  set.seed(1)
  n <- 4000; prs <- rnorm(n); age <- rnorm(n, 55, 8); sex <- rbinom(n, 1, 0.5)
  mu <- exp(0.5 + 0.2 * prs + 0.01 * (age - 55))
  y <- MASS::rnegbin(n, mu = mu, theta = 2)
  r <- fit_count(y, prs, data.frame(age = age, sex = sex))
  expect_identical(r$model, "nb")
  expect_lt(abs(r$estimate - 0.2), 0.05)
  expect_lt(r$p, 1e-10)
  expect_true(r$lr_chisq > 50 && r$pseudo_r2 > 0 && r$ci_low < r$effect && r$effect < r$ci_high)
  r0 <- fit_count(rep(0L, 100), rnorm(100), data.frame(age = rnorm(100)))
  expect_true(is.na(r0$estimate) && grepl("too few", r0$note))
})

test_that("Cox with age as timescale recovers a simulated hazard ratio", {
  skip_if_not_installed("survival")
  set.seed(2)
  n <- 5000; prs <- rnorm(n); sex <- rbinom(n, 1, 0.5); entry <- runif(n, 40, 70)
  t <- rexp(n, rate = 0.02 * exp(0.3 * prs))
  exit <- entry + pmin(t, 12); event <- as.integer(t <= 12)
  r <- fit_cox_age(entry, exit, event, prs, data.frame(sex = sex))
  expect_identical(r$model, "cox")
  expect_lt(abs(r$estimate - 0.3), 0.06)
  expect_gt(r$n_events, 20)
})

test_that("cross-group C+T threshold selection, deciles and IVW meta", {
  res <- data.table(group = rep(c("A", "B"), each = 3), name = "ARD", kind = "factor", factor = "ARD", method = "ct",
                    threshold = rep(c("pT_1e-04", "pT_0.05", "pT_1"), 2), lr_chisq = c(10, 30, 20, 25, 15, 5))
  cmp <- data.table(group = c("A", "B"), name = "mvARD", kind = "comparator", factor = "ARD", method = "ct",
                    threshold = "pT_1", lr_chisq = c(500, 500))        # a comparator on the same factor must not be mixed in
  s <- ct_select_threshold(rbind(res, cmp))
  expect_identical(s[group == "A" & name == "ARD", threshold], "pT_1e-04")   # best in B
  expect_identical(s[group == "B" & name == "ARD", threshold], "pT_0.05")    # best in A
  expect_identical(s[name == "mvARD", threshold], c("pT_1", "pT_1"))
  d <- decile_table(seq(1, 1000), seq(1, 1000) / 100)
  expect_equal(nrow(d), 10); expect_true(all(diff(d$mean) > 0))
  # tied scores (a C+T set with a handful of SNPs): quantile breaks coincide, deciles must still form
  d2 <- decile_table(rep(c(0, 1, 2), c(500, 400, 100)), seq_len(1000))
  expect_equal(nrow(d2), 10); expect_equal(sum(d2$n), 1000L); expect_true(all(d2$n == 100L))
  expect_equal(nrow(decile_table(c(NA, 1), c(1, NA))), 0)
  m <- ivw_meta(0.2, 0.05, 0.1, 0.05)
  expect_equal(m$estimate, 0.15); expect_equal(m$se, sqrt(1 / (2 / 0.05^2)))
  expect_gt(m$q_p, 0.05)
})

test_that("a scored chromosome is reused only while its score matrix is unchanged", {
  td <- tempfile(); dir.create(td)
  mat <- file.path(td, "chr1.tsv"); done <- file.path(td, "chr1.done")
  writeLines(c("SNP\tA1\tx", "rs1\tA\t0.1"), mat)
  expect_false(.pe_checkpoint_ok(done, mat))
  file.create(done)                                   # the old empty marker is not trusted
  expect_false(.pe_checkpoint_ok(done, mat))
  writeLines(unname(tools::md5sum(mat)), done)
  expect_true(.pe_checkpoint_ok(done, mat))
  writeLines(c("SNP\tA1\tx\ty", "rs1\tA\t0.1\t0.2"), mat)   # a new weight set adds a column
  expect_false(.pe_checkpoint_ok(done, mat))
})

test_that("external summary statistics are harmonised: chr:pos IDs, lower-case alleles, sign flip, N rules", {
  tmp <- tempfile(fileext = ".tsv"); dir.create(dirname(tmp), showWarnings = FALSE)
  fwrite(data.table(SNP = c("1:100", "1:200", "rs3", "1:400", "rs5"), Chr = 1L, Position = c(100L, 200L, 300L, 400L, 500L),
                    EA = c("a", "g", "C", "t", "A"), NEA = c("g", "a", "T", "c", "G"), EAF = c(0.3, 0.6, 0.2, 0.5, 0.1),
                    Beta = c(0.1, -0.2, 0.3, 0.4, 0.5), SE = 0.05, `P-value` = c(0.05, 1e-4, 1e-8, 0.5, 0.01),
                    Effective_N = c(1000, 1000, 2000, 3000, 1000)), tmp, sep = "\t")
  ref <- data.table(SNP = c("rs1", "rs2", "rs3", "rs5"), CHR = 1L, BP = c(100L, 200L, 300L, 500L),
                    A1 = c("A", "A", "C", "A"), A2 = c("G", "G", "T", "G"), EAF_A1 = c(0.3, 0.4, 0.2, 0.1))
  spec <- list(id = "t", file = tmp, sep = "tab", col_snp = "SNP", col_chr = "Chr", col_bp = "Position", col_a1 = "EA",
               col_a2 = "NEA", col_eaf = "EAF", col_beta = "Beta", col_se = "SE", col_p = "P-value", col_n = "Effective_N",
               n = "col", flip_sign = TRUE)
  x <- read_external_sumstats(spec, ref)
  expect_identical(x$d$SNP, c("rs1", "rs2", "rs3", "rs5"))             # 1:400 is not in the reference
  expect_identical(x$d$A1, c("A", "G", "C", "A"))                       # upper-cased, orientation untouched
  expect_equal(x$d$beta, -c(0.1, -0.2, 0.3, 0.5))                       # flipped
  expect_equal(x$n, 1000)                                               # median of the N column
  spec2 <- modifyList(spec, list(n = "nhat", col_eaf = NA, flip_sign = FALSE))
  x2 <- read_external_sumstats(spec2, ref)
  expect_equal(x2$d[SNP == "rs2", EAF_A1], 0.6)                         # reference frequency aligned to the file's A1 (G)
  expect_equal(x2$d$beta[1], 0.1)
  expect_true(is.finite(x2$n) && x2$n > 0)
})
