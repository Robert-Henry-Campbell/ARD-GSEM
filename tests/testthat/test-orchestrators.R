source(file.path(project_root, "R/07_sumstats.R"))
source(file.path(project_root, "R/08_gwas.R"))
source(file.path(project_root, "R/09_write_vcf.R"))
source(file.path(project_root, "R/10_plots.R"))
source(file.path(project_root, "R/11_factor_summary.R"))

make_fake_config <- function(tmpdir, sex = "male") {
  list(
    project = list(name = "test", root = tmpdir),
    paths = list(
      output_dir = file.path(tmpdir, "output"),
      meta_dir = file.path(tmpdir, "meta"),
      manifest_dir = file.path(tmpdir, "manifest"),
      thousand_g_reference = file.path(tmpdir, "reference.1000G.maf.0.005.txt")
    ),
    munge = list(maf_threshold = 0.01),
    cfa = list(estimator = "DWLS", min_indicators_per_factor = 2),
    gwas = list(smooth_check = TRUE, genome_build = "GRCh37", std_lv = TRUE),
    parallel = list(n_workers = 1)
  )
}

write_fake_ldsc <- function(config, sex, traits) {
  dir.create(file.path(config$paths$output_dir, sex, "ldsc"),
             recursive = TRUE, showWarnings = FALSE)
  k <- length(traits)
  S <- diag(0.1, k); dimnames(S) <- list(traits, traits)
  V <- diag(1e-4, k * (k + 1L) / 2L)
  ldsc_output <- list(S = S, V = V, I = diag(1, k))
  saveRDS(ldsc_output,
          file.path(config$paths$output_dir, sex, "ldsc", "ldsc_full.rds"))
  fwrite(data.table(trait = traits, h2 = diag(S), se = 0.01,
                    z = diag(S) / 0.01,
                    pass = TRUE),
         file.path(config$paths$output_dir, sex, "ldsc", "h2_qc.csv"))
}

write_fake_munge_inputs <- function(config, sex, traits, n_snps = 50L) {
  dir.create(file.path(config$paths$output_dir, sex, "munge"),
             recursive = TRUE, showWarnings = FALSE)
  for (trait in traits) {
    dt <- data.table(
      SNP = paste0("rs", seq_len(n_snps)),
      CHR = "1", BP = seq_len(n_snps),
      A1 = "A", A2 = "G", MAF = 0.2,
      effect = rnorm(n_snps, 0, 0.01),
      SE = runif(n_snps, 0.005, 0.01),
      P = runif(n_snps), N = 50000
    )
    fwrite(dt, file.path(config$paths$output_dir, sex, "munge",
                         sprintf("%s_%s_pre_munge.tsv", sex, trait)),
           sep = "\t")
  }
}

write_fake_manifest <- function(config, sex, traits) {
  dir.create(config$paths$manifest_dir, recursive = TRUE, showWarnings = FALSE)
  manifest <- data.frame(phenotype = traits,
                         n_cases = 5000,
                         n_controls = 45000)
  manifest_name <- sprintf("neale_%s_manifest", sex)
  assign(manifest_name, manifest)
  save(list = manifest_name,
       file = file.path(config$paths$manifest_dir,
                        sprintf("neale_%s_manifest.rda", sex)))
}

write_fake_plink_ref <- function(config) {
  # GenomicSEM::sumstats() reads `ref` as a single text file via fread().
  # Create a stub file so resolve_1000g_reference() can find it.
  ref_path <- config$paths$thousand_g_reference
  dir.create(dirname(ref_path), recursive = TRUE, showWarnings = FALSE)
  writeLines(c("SNP\tA1\tA2\tMAF", "rs1\tA\tG\t0.2"), ref_path)
}

skip_if_no_mockery <- function() {
  if (!requireNamespace("mockery", quietly = TRUE)) skip("mockery not installed")
}

test_that("run_sumstats wires inputs and writes the prefixed output", {
  skip_if_no_mockery()
  tmp <- tempfile("orch_"); dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  cfg <- make_fake_config(tmp, sex = "male")
  traits <- c("D50", "D64")
  write_fake_manifest(cfg, "male", traits)
  write_fake_ldsc(cfg, "male", traits)
  write_fake_munge_inputs(cfg, "male", traits)
  write_fake_plink_ref(cfg)

  fake_snp <- data.table(SNP = paste0("rs", 1:50), A1 = "A", A2 = "G",
                          CHR = "1", BP = 1:50, MAF = 0.2)
  captured <- list()
  mockery::stub(run_sumstats, "GenomicSEM::sumstats", function(...) {
    captured$args <<- list(...)
    fake_snp
  })

  res <- run_sumstats(cfg, sex = "male")

  expect_true(file.exists(file.path(cfg$paths$output_dir, "male", "sumstats",
                                     "male_snp_sumstats.rds")))
  expect_true(file.exists(file.path(cfg$paths$output_dir, "male", "sumstats",
                                     "male_scale_metadata.json")))
  args <- captured$args
  expect_true(all(args$se.logit), "se.logit must be TRUE for log-odds input")
  expect_true(all(!args$OLS), "OLS must be FALSE")
  expect_true(all(!args$linprob), "linprob must be FALSE (already log-odds)")
  expect_equal(args$trait.names, traits)
})

test_that("run_gwas builds the augmented model and saves the prefixed outputs", {
  skip_if_no_mockery()
  tmp <- tempfile("orch_"); dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  cfg <- make_fake_config(tmp, sex = "male")
  traits <- c("D50", "D64")
  write_fake_manifest(cfg, "male", traits)
  write_fake_ldsc(cfg, "male", traits)

  dir.create(file.path(cfg$paths$output_dir, "male", "sumstats"),
             recursive = TRUE, showWarnings = FALSE)
  fake_snp <- data.table(SNP = paste0("rs", 1:30), A1 = "A", A2 = "G",
                          CHR = "1", BP = 1:30, MAF = 0.2)
  saveRDS(fake_snp, file.path(cfg$paths$output_dir, "male", "sumstats",
                               "male_snp_sumstats.rds"))

  dir.create(cfg$paths$meta_dir, recursive = TRUE, showWarnings = FALSE)
  fwrite(data.table(code = c("D50", "D64"),
                    chapter = c("Diseases of the blood", "Diseases of the blood")),
         file.path(cfg$paths$meta_dir, "icd10_categories.csv"))

  captured <- list()
  fake_factor_df <- data.table(SNP = fake_snp$SNP, CHR = "1", BP = 1:30,
                                 A1 = "A", A2 = "G", MAF = 0.2,
                                 est = rnorm(30), SE = runif(30, 0.01, 0.02),
                                 Pval_Estimate = runif(30),
                                 Q_SNP = abs(rnorm(30)), Q_SNP_df = 1L,
                                 Q_SNP_pval = runif(30))
  mockery::stub(run_gwas, "GenomicSEM::userGWAS", function(model, sub, ..., std.lv) {
    captured$model <<- model
    captured$sub <<- sub
    captured$std.lv <<- std.lv
    rep(list(fake_factor_df), length(sub))
  })

  res <- run_gwas(cfg, sex = "male")

  expect_true(file.exists(file.path(cfg$paths$output_dir, "male", "gwas",
                                     "male_userGWAS_raw.rds")))
  expect_true(file.exists(file.path(cfg$paths$output_dir, "male", "gwas",
                                     "male_model_spec.rds")))
  expect_match(captured$model, "=~ D50 \\+ D64")
  expect_true(any(grepl("~ SNP$", captured$sub)))
  expect_true(isTRUE(captured$std.lv),
              "std.lv = TRUE must be passed to GenomicSEM::userGWAS by default")
})

test_that("run_gwas wraps a single-factor userGWAS return value (length(sub)==1) into a list", {
  skip_if_no_mockery()
  tmp <- tempfile("orch_"); dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  cfg <- make_fake_config(tmp, sex = "male")
  traits <- c("D50", "D64")   # one chapter -> one factor
  write_fake_manifest(cfg, "male", traits)
  write_fake_ldsc(cfg, "male", traits)
  dir.create(file.path(cfg$paths$output_dir, "male", "sumstats"),
             recursive = TRUE, showWarnings = FALSE)
  fake_snp <- data.table(SNP = paste0("rs", 1:30), A1 = "A", A2 = "G",
                          CHR = "1", BP = 1:30, MAF = 0.2)
  saveRDS(fake_snp, file.path(cfg$paths$output_dir, "male", "sumstats",
                               "male_snp_sumstats.rds"))

  dir.create(cfg$paths$meta_dir, recursive = TRUE, showWarnings = FALSE)
  fwrite(data.table(code = c("D50", "D64"),
                    chapter = c("Diseases of the blood", "Diseases of the blood")),
         file.path(cfg$paths$meta_dir, "icd10_categories.csv"))

  fake_factor_df <- data.table(SNP = fake_snp$SNP, CHR = "1", BP = 1:30,
                                 A1 = "A", A2 = "G", MAF = 0.2,
                                 est = rnorm(30), SE = runif(30, 0.01, 0.02),
                                 Pval_Estimate = runif(30),
                                 Q_SNP = abs(rnorm(30)), Q_SNP_df = 1L,
                                 Q_SNP_pval = runif(30))
  # Stub userGWAS to return a bare data.frame (not a length-1 list), as
  # the real userGWAS does when sub has length 1.
  mockery::stub(run_gwas, "GenomicSEM::userGWAS", function(...) {
    fake_factor_df
  })

  res <- run_gwas(cfg, sex = "male")

  out <- readRDS(file.path(cfg$paths$output_dir, "male", "gwas",
                            "male_userGWAS_raw.rds"))
  expect_true(is.list(out))
  expect_equal(length(out), 1L)
  expect_true(is.data.frame(out[[1]]))
})

test_that("run_write_vcf reads sex-prefixed inputs and produces sex-prefixed VCF", {
  skip_if_no_mockery()
  tmp <- tempfile("orch_"); dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  cfg <- make_fake_config(tmp, sex = "male")

  set.seed(5)
  fake_factor_df <- data.table(SNP = paste0("rs", 1:30), CHR = "1", BP = 1:30,
                                 A1 = "A", A2 = "G", MAF = 0.2,
                                 est = rnorm(30, 0, 0.01), SE = runif(30, 0.005, 0.01),
                                 Pval_Estimate = runif(30),
                                 Q_SNP = abs(rnorm(30)), Q_SNP_df = 1L,
                                 Q_SNP_pval = runif(30))
  dir.create(file.path(cfg$paths$output_dir, "male", "gwas"),
             recursive = TRUE, showWarnings = FALSE)
  saveRDS(list(Blood = fake_factor_df),
          file.path(cfg$paths$output_dir, "male", "gwas",
                    "male_userGWAS_raw.rds"))

  res <- run_write_vcf(cfg, sex = "male")
  expect_true(file.exists(file.path(cfg$paths$output_dir, "male", "gwas",
                                     "male_chapter_gsem.vcf.gz")))
  expect_equal(res$n_factors, 1L)
  expect_equal(res$n_variants, 30L)
})
