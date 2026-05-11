resolve_1000g_prefix <- function(dir_path) {
  if (!dir.exists(dir_path)) {
    log_fatal("sumstats", sprintf("1000G PLINK reference dir not found: %s", dir_path))
  }
  bed_files <- list.files(dir_path, pattern = "\\.bed$", full.names = TRUE)
  if (length(bed_files) == 0L) {
    log_fatal("sumstats", sprintf("No .bed file found in 1000G PLINK dir: %s", dir_path))
  }
  prefixes <- unique(sub("\\.bed$", "", bed_files))
  if (length(prefixes) != 1L) {
    log_fatal("sumstats", sprintf(
      "Expected exactly one .bed/.bim/.fam trio in %s, found %d candidate prefix(es)",
      dir_path, length(prefixes)))
  }
  prefixes[1L]
}

run_sumstats <- function(config, sex) {
  log_info("sumstats", sprintf("=== Sumstats stage: %s ===", sex))

  h2_path <- file.path(config$paths$output_dir, sex, "ldsc", "h2_qc.csv")
  if (!file.exists(h2_path)) {
    log_fatal("sumstats", sprintf("h2_qc.csv not found at %s (run LDSC first)", h2_path))
  }
  h2_table <- fread(h2_path)
  retained <- h2_table[pass == TRUE]$trait

  if (length(retained) < 2L) {
    log_fatal("sumstats", sprintf("Need >=2 retained traits, found %d", length(retained)))
  }

  munge_dir <- file.path(config$paths$output_dir, sex, "munge")
  files <- file.path(munge_dir, paste0(retained, "_pre_munge.tsv"))
  missing <- files[!file.exists(files)]
  if (length(missing) > 0L) {
    log_fatal("sumstats", sprintf(
      "Pre-munge TSV(s) missing for retained traits: %s (re-run --stage munge)",
      paste(basename(missing), collapse = ", ")))
  }

  cc_list <- lapply(retained, function(trait) get_case_control(config, sex, trait))
  sample_prev <- vapply(cc_list, function(cc) cc$n_cases / (cc$n_cases + cc$n_controls),
                        numeric(1))
  neffs <- vapply(cc_list, function(cc) compute_neff(cc$n_cases, cc$n_controls),
                  numeric(1))

  log_info("sumstats", sprintf("Calling GenomicSEM::sumstats() for %d traits (%s)",
                                length(retained), paste(retained, collapse = ", ")))
  log_info("sumstats", sprintf(
    "Input scale: log(OR) (converted at munge time). se.logit=TRUE, OLS=FALSE, linprob=FALSE."))

  ref_prefix <- resolve_1000g_prefix(config$paths$thousand_g_plink)
  log_info("sumstats", sprintf("1000G reference prefix: %s", ref_prefix))

  snp_sumstats <- GenomicSEM::sumstats(
    files = files,
    ref = ref_prefix,
    trait.names = retained,
    se.logit = rep(TRUE, length(retained)),
    OLS = rep(FALSE, length(retained)),
    linprob = rep(FALSE, length(retained)),
    N = neffs,
    betas = "effect",
    ses = "SE",
    info.filter = 0.6,
    maf.filter = config$munge$maf_threshold,
    keep.indel = FALSE,
    parallel = TRUE,
    cores = config$parallel$n_workers
  )

  out_dir <- file.path(config$paths$output_dir, sex, "sumstats")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(snp_sumstats, file.path(out_dir, "snp_sumstats.rds"))

  scale_meta <- list(
    scale = "log_odds_ratio",
    source = "Neale UKBB round 2 linear regression on 0/1 phenotype",
    conversion = "beta_logOR = beta_linear / (K * (1 - K))",
    sample_prev = setNames(as.list(sample_prev), retained),
    neff = setNames(as.list(neffs), retained),
    sumstats_call = list(
      se.logit = TRUE, OLS = FALSE, linprob = FALSE,
      info.filter = 0.6, maf.filter = config$munge$maf_threshold
    )
  )
  jsonlite::write_json(scale_meta, file.path(out_dir, "scale_metadata.json"),
                       auto_unbox = TRUE, pretty = TRUE)

  n_snps <- if (is.data.frame(snp_sumstats)) nrow(snp_sumstats) else NA_integer_
  log_info("sumstats", sprintf("Aligned SNP-level sumstats saved: %s SNPs x %d traits",
                                format(n_snps, big.mark = ","), length(retained)))

  write_stage_manifest("sumstats", sex, config, list.files(out_dir, full.names = TRUE))

  list(
    snp_sumstats = snp_sumstats,
    retained = retained,
    sample_prev = sample_prev,
    neffs = neffs
  )
}
