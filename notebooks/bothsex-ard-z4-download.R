# bothsex ARD -> Pan-UKB ICD10 sumstats download
#
# Step through this top-to-bottom (Ctrl/Cmd-Enter line-by-line).
# Working dir assumed = GSEM/notebooks/.

# ---- 0. Setup ---------------------------------------------------------------
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
})

setwd("notebooks")
manifest_path <- "../manifest/panukb_phenotype_manifest.csv"
rda_path      <- "../manifest/bothsex_ARD.rda"
dest_dir      <- "../sumstats/bothsex"
csv_out       <- "bothsex_ARD_panukb_icd10_filtered.csv"

dir.create(dest_dir, showWarnings = FALSE, recursive = TRUE)
cat("Destination for sumstats:", normalizePath(dest_dir), "\n")


# ---- 1. Load Pan-UKB phenotype manifest -------------------------------------
manifest <- readr::read_csv(
  manifest_path,
  show_col_types = FALSE,
  guess_max      = 50000
)

cat("manifest dim:", dim(manifest), "\n")
print(sort(table(manifest$trait_type), decreasing = TRUE))

# Peek at the columns we'll use
manifest |>
  dplyr::filter(trait_type == "icd10") |>
  dplyr::select(trait_type, phenocode, description,
                sldsc_25bin_h2_liability_EUR,
                sldsc_25bin_h2_z_EUR,
                filename, wget) |>
  head(5) |>
  print()


# ---- 2. Load bothsex_ARD.rda ------------------------------------------------
loaded <- load(rda_path)
cat("Objects loaded:", paste(loaded, collapse = ", "), "\n")

ard <- get(loaded[1])
cat("class:", paste(class(ard), collapse = " / "), "\n")
cat("dim:",   dim(ard), "\n")
print(colnames(ard))
print(head(ard, 5))


# ---- 3. Pull ICD10 codes from bothsex_ARD -----------------------------------
# ICD10_explo holds the 3-char codes (B17, C92, ...) -- the form Pan-UKB
# indexes ICD10 GWAS by (trait_type == "icd10", phenocode == 3-char code).
ard_icd10 <- ard$ICD10_explo |>
  stringr::str_trim() |>
  (\(x) x[!is.na(x) & nzchar(x)])() |>
  unique() |>
  sort()

cat("Unique ICD10 (3-char) codes in bothsex_ARD:", length(ard_icd10), "\n")
print(ard_icd10)


# ---- 4. Filter the manifest -------------------------------------------------
# Keep only Pan-UKB ICD10 GWAS that (a) match an ARD ICD10 code and
# (b) have EUR h2 Z-score > 4 (sldsc_25bin_h2_z_EUR) -- the Pan-UKB-recommended
# heritability-significance threshold.
Z_MIN <- 0.01

matched <- manifest |>
  dplyr::filter(trait_type == "icd10", phenocode %in% ard_icd10) |>
  dplyr::select(
    trait_type, phenocode, pheno_sex, description,
    h2_liability_EUR = sldsc_25bin_h2_liability_EUR,
    h2_z_EUR         = sldsc_25bin_h2_z_EUR,
    filename, aws_link, wget,
    md5_hex, size_in_bytes
  ) |>
  dplyr::arrange(phenocode)

filtered <- matched |>
  dplyr::filter(!is.na(h2_z_EUR), h2_z_EUR > Z_MIN)

cat("ARD codes:                          ", length(ard_icd10), "\n")
cat("Pan-UKB ICD10 GWAS matching codes:  ", nrow(matched), "\n")
cat("...with EUR h2 Z >", Z_MIN, ":             ", nrow(filtered), "\n\n")

missing_codes  <- setdiff(ard_icd10, matched$phenocode)
dropped_for_Z  <- matched |>
  dplyr::filter(is.na(h2_z_EUR) | h2_z_EUR <= Z_MIN) |>
  dplyr::select(phenocode, description, h2_z_EUR)

cat("ARD codes with no Pan-UKB ICD10 GWAS (", length(missing_codes), "): ",
    paste(missing_codes, collapse = ", "), "\n", sep = "")
cat("Dropped by Z <=", Z_MIN, "or NA (", nrow(dropped_for_Z), "):\n", sep = "")
print(dropped_for_Z, n = Inf)


# ---- 5. Print filename + h2 + Z, and save CSV -------------------------------
filtered |>
  dplyr::select(phenocode, description, filename, h2_liability_EUR, h2_z_EUR) |>
  print(n = Inf)

readr::write_csv(filtered, csv_out)
cat("wrote:", normalizePath(csv_out), "  (", nrow(filtered), "rows)\n")


# ---- 6. Download the sumstats -----------------------------------------------
download_one <- function(wget_cmd, fname, dest) {
  out_path <- file.path(dest, fname)
  if (file.exists(out_path)) {
    cat("[skip] ", fname, " (already present)\n", sep = "")
    return(invisible(0L))
  }
  # manifest's `wget` column looks like: 'wget https://...'
  cmd <- sub("^wget\\s+",
             paste0("wget --no-clobber --directory-prefix=",
                    shQuote(dest), " "),
             wget_cmd)
  cat("[get ] ", fname, "\n", sep = "")
  status <- system(cmd)
  if (status != 0) warning("wget failed for ", fname, " (exit ", status, ")")
  invisible(status)
}

# Dry-run preview
print(head(filtered$wget, 3))

# Actual download loop (idempotent -- skips files already in dest_dir)
for (i in seq_len(nrow(filtered))) {
  download_one(filtered$wget[i], filtered$filename[i], dest_dir)
}

cat("\nFiles now in", dest_dir, ":\n")
print(list.files(dest_dir))


# ---- 7. (optional) verify sizes against the manifest ------------------------
check <- filtered |>
  dplyr::mutate(
    local_path   = file.path(dest_dir, filename),
    exists       = file.exists(local_path),
    local_bytes  = ifelse(exists, file.info(local_path)$size, NA_real_),
    size_matches = !is.na(local_bytes) & local_bytes == size_in_bytes
  ) |>
  dplyr::select(phenocode, filename, exists, local_bytes, size_in_bytes, size_matches)

print(check, n = Inf)
cat("\nall downloaded & size-matched:", all(check$size_matches, na.rm = TRUE), "\n")












