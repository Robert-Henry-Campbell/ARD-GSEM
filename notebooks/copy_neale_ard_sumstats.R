#!/usr/bin/env Rscript
# Copy Neale sex-stratified ARD sumstats from the cache into GSEM sumstats dirs.
# Run from the GSEM project root:
#   Rscript notebooks/copy_neale_ard_sumstats.R

# ---- 0. Paths ---------------------------------------------------------------
cache_dir  <- "/mnt/sdg/robert/ardmr/ardmr_cache/neale_sumstats"
dest_male  <- "sumstats/male"
dest_fem   <- "sumstats/female"

dir.create(dest_male, recursive = TRUE, showWarnings = FALSE)
dir.create(dest_fem,  recursive = TRUE, showWarnings = FALSE)


# ---- 1. Load manifests ------------------------------------------------------
e <- new.env()

load("manifest/male_ARD.rda",          envir = e)
load("manifest/female_ARD.rda",        envir = e)
load("manifest/neale_file_manifest.rda", envir = e)

male_ard   <- e$male_ARD
female_ard <- e$female_ARD
file_mf    <- e$neale_file_manifest


# ---- 2. Extract unique 3-char ICD10 codes from each ARD ---------------------
ard_codes <- function(ard) {
  codes <- ard$ICD10_explo
  sort(unique(trimws(codes[!is.na(codes) & nzchar(trimws(codes))])))
}

male_codes   <- ard_codes(male_ard)
female_codes <- ard_codes(female_ard)

cat(sprintf("ARD ICD10 codes  -- male: %d  female: %d\n",
            length(male_codes), length(female_codes)))


# ---- 3. Match codes to files via neale_file_manifest ------------------------
# neale_file_manifest columns: "Phenotype Code", "Sex", "File", ...
# Sex values: "male", "female", "both_sexes"

get_files <- function(codes, sex_label) {
  sub <- file_mf[
    file_mf[["Phenotype Code"]] %in% codes &
    file_mf[["Sex"]] == sex_label,
  ]
  sub[, c("Phenotype Code", "File")]
}

male_files   <- get_files(male_codes,   "male")
female_files <- get_files(female_codes, "female")

# Codes in ARD but absent from file manifest
report_missing_codes <- function(codes, files_df, sex_label) {
  found  <- unique(files_df[["Phenotype Code"]])
  absent <- setdiff(codes, found)
  if (length(absent) > 0) {
    cat(sprintf("[WARN] %d %s ARD code(s) have no entry in neale_file_manifest: %s\n",
                length(absent), sex_label, paste(absent, collapse = ", ")))
  }
}
report_missing_codes(male_codes,   male_files,   "male")
report_missing_codes(female_codes, female_files, "female")

cat(sprintf("Files to copy    -- male: %d  female: %d\n",
            nrow(male_files), nrow(female_files)))


# ---- 4. Copy loop ------------------------------------------------------------
copy_files <- function(files_df, dest) {
  n_copied  <- 0L
  n_skipped <- 0L
  n_missing <- 0L

  for (i in seq_len(nrow(files_df))) {
    fname  <- files_df[["File"]][i]
    code   <- files_df[["Phenotype Code"]][i]
    src    <- file.path(cache_dir, fname)
    dst    <- file.path(dest, fname)

    if (file.exists(dst)) {
      cat(sprintf("  [skip]   %s\n", fname))
      n_skipped <- n_skipped + 1L
      next
    }

    if (!file.exists(src)) {
      cat(sprintf("  [MISS]   %s  (not in cache)\n", fname))
      n_missing <- n_missing + 1L
      next
    }

    ok <- file.copy(src, dst, overwrite = FALSE)
    if (ok) {
      cat(sprintf("  [copy]   %s\n", fname))
      n_copied <- n_copied + 1L
    } else {
      cat(sprintf("  [ERROR]  %s  (file.copy returned FALSE)\n", fname))
    }
  }

  cat(sprintf("  -> copied: %d  skipped: %d  missing from cache: %d\n\n",
              n_copied, n_skipped, n_missing))
  invisible(list(copied = n_copied, skipped = n_skipped, missing = n_missing))
}

cat("\n=== male -> sumstats/male ===\n")
r_m <- copy_files(male_files, dest_male)

cat("=== female -> sumstats/female ===\n")
r_f <- copy_files(female_files, dest_fem)


# ---- 5. Summary -------------------------------------------------------------
cat("=== Summary ===\n")
cat(sprintf("male   : %d copied, %d already present, %d not in cache\n",
            r_m$copied, r_m$skipped, r_m$missing))
cat(sprintf("female : %d copied, %d already present, %d not in cache\n",
            r_f$copied, r_f$skipped, r_f$missing))
cat(sprintf("sumstats/male   now has %d files\n",   length(list.files(dest_male))))
cat(sprintf("sumstats/female now has %d files\n\n", length(list.files(dest_fem))))

if (r_m$missing + r_f$missing > 0) {
  cat("[WARN] Some files were not in the cache. Check the [MISS] lines above.\n")
} else {
  cat("All files accounted for.\n")
}
