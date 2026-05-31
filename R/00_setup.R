suppressPackageStartupMessages({
  library(data.table)
  library(yaml)
  library(jsonlite)
})

source("R/utils.R")

setup_pipeline <- function(config, sex, threads = NULL) {
  if (!is.null(threads)) config$parallel$n_workers <- threads

  sexes <- as.character(sex)
  valid <- c("male", "female", "bothsex", "bothsex_meta")
  bad <- setdiff(sexes, valid)
  if (length(bad) > 0L) stop(sprintf("setup_pipeline: invalid sex value(s): %s",
                                      paste(bad, collapse = ",")))

  init_logging(config)
  log_info("setup", sprintf("Config: %s", config$project$name))
  log_info("setup", sprintf("Threads: %d | Sex: %s",
                            config$parallel$n_workers,
                            paste(sexes, collapse = ",")))
  report_tempdir()

  validate_reference(config, sexes)

  for (s in sexes) {
    traits <- discover_traits(config, s)
    log_info("setup", sprintf("Traits detected (%s): %s", s, paste(traits, collapse = ", ")))
  }

  both_strat <- all(c("male", "female") %in% sexes)
  if (both_strat) {
    shared <- get_shared_traits(config)
    log_info("setup", sprintf("Shared traits (male ∩ female): %s",
                              paste(shared, collapse = ", ")))
  }

  for (s in sexes) {
    for (stage in c("munge", "ldsc", "efa", "cfa", "sumstats", "gwas")) {
      dir.create(file.path(config$paths$output_dir, s, stage),
                 recursive = TRUE, showWarnings = FALSE)
    }
  }
  dir.create(file.path(config$paths$output_dir, "comparison"),
             recursive = TRUE, showWarnings = FALSE)

  config
}

report_tempdir <- function(warn_threshold_gb = 200) {
  td <- tempdir()
  td_env <- Sys.getenv("TMPDIR", unset = NA)
  avail_gb <- tryCatch({
    out <- suppressWarnings(system2(
      "df", c("-B1", "--output=avail", shQuote(td)),
      stdout = TRUE, stderr = FALSE))
    if (length(out) >= 2L) as.numeric(out[2]) / 1024^3 else NA_real_
  }, error = function(e) NA_real_)
  env_str <- if (is.na(td_env) || !nzchar(td_env)) "(TMPDIR env unset)"
             else sprintf("TMPDIR=%s", td_env)
  if (is.na(avail_gb)) {
    log_info("setup", sprintf("R tempdir: %s  [%s; free space: unknown]", td, env_str))
  } else if (avail_gb < warn_threshold_gb) {
    log_warn("setup", sprintf(
      "R tempdir: %s  [%s; only %.1f GB free]. fread of Neale bgz files needs ~1-2 GB per trait. If munge fails with 'External command failed' / 'disk is full in the temporary directory', launch R with TMPDIR=/path/on/big/mount Rscript ...",
      td, env_str, avail_gb))
  } else {
    log_info("setup", sprintf("R tempdir: %s  [%s; %.1f GB free]",
                              td, env_str, avail_gb))
  }
}

validate_reference <- function(config, sexes = c("male", "female")) {
  ld_dir <- config$paths$ld_scores
  hm3 <- config$paths$hm3_snplist
  vm <- config$paths$variants_manifest
  meta <- config$paths$meta_dir
  manifest_dir <- config$paths$manifest_dir
  sumstats_dir <- config$paths$sumstats_dir

  if (!dir.exists(ld_dir)) {
    log_fatal("setup", sprintf("LD scores directory not found: %s", ld_dir))
  }
  ld_files <- list.files(ld_dir, pattern = "\\.l2\\.ldscore\\.gz$")
  if (length(ld_files) < 22) {
    log_fatal("setup", sprintf("Expected 22 LD score files, found %d in %s", length(ld_files), ld_dir))
  }

  if (!file.exists(hm3)) {
    log_fatal("setup", sprintf("HM3 SNP list not found: %s", hm3))
  }
  needs_neale <- any(c("male", "female") %in% sexes)
  if (needs_neale && !file.exists(vm)) {
    log_fatal("setup", sprintf("Variants manifest not found: %s", vm))
  }
  # bothsex_meta sources chapter info from the FinnGen manifest, not icd10_categories.csv
  needs_icd_csv <- any(sexes %in% c("male", "female", "bothsex"))
  icd_path <- file.path(meta, "icd10_categories.csv")
  if (needs_icd_csv && !file.exists(icd_path)) {
    log_fatal("setup", sprintf("ICD-10 categories file not found: %s", icd_path))
  }

  for (s in sexes) {
    if (identical(s, "bothsex")) {
      rda <- file.path(manifest_dir, "panukb_bothsex_manifest.rda")
      if (!file.exists(rda)) {
        log_fatal("setup", sprintf("Pan-UKB bothsex manifest not found: %s", rda))
      }
      pvm <- config$paths$panukb_variants_manifest
      if (is.null(pvm) || !nzchar(pvm) || !file.exists(pvm)) {
        log_fatal("setup", sprintf("Pan-UKB variants manifest not found: %s",
                                    as.character(pvm)))
      }
    } else if (identical(s, "bothsex_meta")) {
      rda <- file.path(manifest_dir, "bothsex_meta_manifest.rda")
      if (!file.exists(rda)) {
        log_fatal("setup", sprintf("bothsex_meta manifest not found: %s. Run scripts/generate_bothsex_meta_manifest.R.", rda))
      }
      meta_tsv <- config$bothsex_meta$manifest
      if (!is.null(meta_tsv) && !is_absolute_path(meta_tsv)) {
        meta_tsv <- file.path(config$project$root %||% getwd(), meta_tsv)
      }
      if (is.null(meta_tsv) || !nzchar(meta_tsv) || !file.exists(meta_tsv)) {
        log_fatal("setup", sprintf("bothsex_meta phenotype manifest TSV not found: %s",
                                    as.character(meta_tsv)))
      }
      chain <- config$bothsex_meta$chain_file
      if (!is.null(chain) && !is_absolute_path(chain)) {
        chain <- file.path(config$project$root %||% getwd(), chain)
      }
      if (is.null(chain) || !nzchar(chain) || !file.exists(chain)) {
        log_fatal("setup", sprintf("bothsex_meta liftover chain file not found: %s. Run bash scripts/download_chain.sh.",
                                    as.character(chain)))
      }
      if (!requireNamespace("rtracklayer", quietly = TRUE)) {
        log_fatal("setup", "rtracklayer is required for bothsex_meta liftover but is not installed. Install with: Rscript -e 'BiocManager::install(\"rtracklayer\")'")
      }
      # R12: load the .rda once and verify every configured af_cohort's
      # manifest_n_cases / manifest_n_controls column exists. A typo in the
      # YAML would otherwise silently zero that cohort's weight at runtime.
      env_rda <- new.env()
      load(rda, envir = env_rda)
      objs <- ls(env_rda)
      if (length(objs) != 1L) {
        log_fatal("setup", sprintf("bothsex_meta manifest .rda must contain exactly 1 object, found %d (%s)",
                                    length(objs), paste(objs, collapse = ", ")))
      }
      mf <- get(objs[1L], envir = env_rda)
      mf_cols <- names(mf)
      af_cohorts_cfg <- config$bothsex_meta$af_cohorts
      if (length(af_cohorts_cfg) == 0L) {
        log_fatal("setup", "config$bothsex_meta$af_cohorts is empty; need at least one cohort")
      }
      for (co in af_cohorts_cfg) {
        for (key in c("manifest_n_cases", "manifest_n_controls")) {
          col_name <- co[[key]]
          if (is.null(col_name) || !nzchar(col_name) || !(col_name %in% mf_cols)) {
            log_fatal("setup", sprintf(
              "bothsex_meta af_cohort '%s': %s='%s' not found in %s (available: %s). Check config$bothsex_meta$af_cohorts and re-run scripts/generate_bothsex_meta_manifest.R if needed.",
              co$name %||% "<unnamed>", key,
              as.character(col_name), basename(rda),
              paste(mf_cols, collapse = ", ")))
          }
        }
      }
    } else {
      rda <- file.path(manifest_dir, sprintf("neale_%s_manifest.rda", s))
      if (!file.exists(rda)) {
        log_fatal("setup", sprintf("Neale %s manifest not found: %s", s, rda))
      }
    }
    sd <- file.path(sumstats_dir, s)
    if (!dir.exists(sd) || length(list.files(sd)) == 0L) {
      log_fatal("setup", sprintf("Sumstats dir empty or missing: %s", sd))
    }
  }

  tg <- config$paths$thousand_g_reference %||% config$paths$thousand_g_plink
  if (!is.null(tg) && nzchar(tg)) {
    if (!file.exists(tg)) {
      log_warn("setup", sprintf(
        "1000G reference file not found: %s (sumstats/gwas/vcf stages will fail). Download from https://utexas.box.com/s/vkd36n197m8klbaio3yzoxsee6sxo11v",
        tg))
    } else {
      log_info("setup", sprintf("1000G reference file present: %s", tg))
    }
  }

  log_info("setup", sprintf("Reference validated: eur_w_ld_chr (%d chr), w_hm3.snplist, variants_manifest, icd10_categories, neale manifests for %s",
                            length(ld_files), paste(sexes, collapse = "+")))
}
