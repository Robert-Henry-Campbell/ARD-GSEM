suppressPackageStartupMessages({
  library(data.table)
  library(yaml)
  library(jsonlite)
})

source("R/utils.R")

setup_pipeline <- function(config, sex = "both", mode = "smoke", threads = NULL) {
  if (!is.null(threads)) config$parallel$n_workers <- threads

  init_logging(config)
  log_info("setup", sprintf("Config: %s", config$project$name))
  log_info("setup", sprintf("Threads: %d | Mode: %s | Sex: %s",
                            config$parallel$n_workers, mode, sex))

  sexes <- if (sex == "both") c("male", "female") else sex
  validate_reference(config, sexes)

  for (s in sexes) {
    traits <- discover_traits(config, s)
    log_info("setup", sprintf("Traits detected (%s): %s", s, paste(traits, collapse = ", ")))
  }

  if (sex == "both") {
    shared <- get_shared_traits(config)
    log_info("setup", sprintf("Shared traits: %s", paste(shared, collapse = ", ")))
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
  if (!file.exists(vm)) {
    log_fatal("setup", sprintf("Variants manifest not found: %s", vm))
  }
  icd_path <- file.path(meta, "icd10_categories.csv")
  if (!file.exists(icd_path)) {
    log_fatal("setup", sprintf("ICD-10 categories file not found: %s", icd_path))
  }

  for (s in sexes) {
    rda <- file.path(manifest_dir, sprintf("neale_%s_manifest.rda", s))
    if (!file.exists(rda)) {
      log_fatal("setup", sprintf("Neale %s manifest not found: %s", s, rda))
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
