resolve_1000g_reference <- function(path) {
  # GenomicSEM::sumstats() reads `ref` via fread() as a single tab-delimited file
  # (typically reference.1000G.maf.0.005.txt -- SNP, A1, A2, MAF columns).
  # It is NOT a PLINK .bed/.bim/.fam prefix. fread() transparently decompresses .gz.
  candidates <- unique(c(path, paste0(path, ".gz"), sub("\\.gz$", "", path)))
  hit <- candidates[file.exists(candidates)][1]
  if (is.na(hit)) {
    log_fatal("sumstats",
              sprintf("1000G reference file not found (tried: %s)",
                      paste(candidates, collapse = ", ")))
  }
  normalizePath(hit, winslash = "/", mustWork = TRUE)
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
  prefix <- paste0(sex, "_")
  files <- file.path(munge_dir, paste0(prefix, retained, "_pre_munge.tsv"))
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

  ref_path <- resolve_1000g_reference(
    config$paths$thousand_g_reference %||% config$paths$thousand_g_plink)
  log_info("sumstats", sprintf("1000G reference file: %s", ref_path))

  out_dir <- file.path(config$paths$output_dir, sex, "sumstats")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  # GenomicSEM::sumstats() writes paste0(<traits>, "_sumstats.log") relative to CWD
  # (same pattern as munge). Force CWD = out_dir for the call so the log lands in the
  # sex-specific sumstats directory.
  files_abs <- normalizePath(files, mustWork = TRUE)
  saved_wd <- getwd()
  snp_sumstats <- tryCatch({
    setwd(out_dir)
    GenomicSEM::sumstats(
      files = files_abs,
      ref = ref_path,
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
  }, finally = setwd(saved_wd))

  saveRDS(snp_sumstats, file.path(out_dir, paste0(sex, "_snp_sumstats.rds")))

  # Informational only: a JSON provenance record of the scale conversion + sumstats() call.
  # Not consumed by any downstream stage in this pipeline; intended for human/audit use and
  # for any external consumer that needs to know the units of the VCF effect estimates.
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
  jsonlite::write_json(scale_meta,
                       file.path(out_dir, paste0(sex, "_scale_metadata.json")),
                       auto_unbox = TRUE, pretty = TRUE)

  n_snps <- if (is.data.frame(snp_sumstats)) nrow(snp_sumstats) else NA_integer_
  log_info("sumstats", sprintf("Aligned SNP-level sumstats saved: %s SNPs x %d traits",
                                format(n_snps, big.mark = ","), length(retained)))

  if (isTRUE(config$sumstats$cleanup_pre_munge)) {
    removed <- 0L
    for (f in files) {
      if (file.exists(f) && file.remove(f)) removed <- removed + 1L
    }
    log_info("sumstats", sprintf(
      "cleanup_pre_munge=true: removed %d/%d pre-munge TSV(s) from %s",
      removed, length(files), file.path(config$paths$output_dir, sex, "munge")))
  }

  write_stage_manifest("sumstats", sex, config, list.files(out_dir, full.names = TRUE))

  list(
    snp_sumstats = snp_sumstats,
    retained = retained,
    sample_prev = sample_prev,
    neffs = neffs
  )
}
