#!/usr/bin/env Rscript
# Derives manifest/bothsex_meta_manifest.rda from the FinnGen+UKBB phenotype
# manifest TSV. The sumstats dir mixes two file types (from the same FinnGen
# release): <TRAIT>_meta_out.tsv.gz (2-way FG+UKBB inverse-variance meta) and
# <TRAIT>_fg.tsv.gz (raw single-cohort FinnGen GWAS). file_type is detected
# from the suffix on disk and recorded in the .rda so run_munge_bothsex_meta
# can dispatch per file.
#
# Inclusion filter:
#   - <TRAIT>_meta_out.tsv.gz: require BOTH fg_n_cases AND ukbb_n_cases non-NA
#   - <TRAIT>_fg.tsv.gz:       require ONLY fg_n_cases non-NA
#
# The .rda holds:
#   phenotype       : FinnGen endpoint id (= fg_phenotype)
#   name            : display name
#   category        : raw FinnGen category text
#   chapter         : canonical ICD chapter (Roman-numeral-stripped) or "Unclassified"
#   chapter_source  : "category-table" | "prefix-fallback" | "unclassified"
#   file_type       : "meta" or "fg"
#   fg_n_cases, fg_n_controls          : always populated
#   ukbb_n_cases, ukbb_n_controls      : NA for fg-only traits
#   n_cases, n_controls                : pooled (fg + ukbb, NA->0) -- logging only

find_project_root <- function() {
  full_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- full_args[grep("^--file=", full_args)]
  start_d <- if (length(file_arg) > 0L) {
    dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = FALSE))
  } else normalizePath(getwd(), mustWork = FALSE)
  d <- start_d
  while (d != dirname(d)) {
    if (file.exists(file.path(d, "config", "pipeline.yaml"))) return(d)
    d <- dirname(d)
  }
  stop("project root not found")
}

project_root <- find_project_root()
setwd(project_root)
source("R/00_setup.R")
config <- read_config("config/pipeline.yaml", root = project_root)

meta_cfg <- config$bothsex_meta
if (is.null(meta_cfg)) stop("config$bothsex_meta block missing")

tsv_path <- meta_cfg$manifest
if (!is.null(tsv_path) && !is_absolute_path(tsv_path)) {
  tsv_path <- file.path(project_root, tsv_path)
}
if (!file.exists(tsv_path)) {
  stop("bothsex_meta manifest TSV not found: ", tsv_path)
}

sumstats_dir <- file.path(config$paths$sumstats_dir, "bothsex_meta")
if (!dir.exists(sumstats_dir)) {
  stop("bothsex_meta sumstats dir not found: ", sumstats_dir)
}
all_files <- list.files(sumstats_dir,
                         pattern = "_(meta_out|fg)\\.tsv\\.gz$")
file_type_by_trait <- ifelse(grepl("_meta_out\\.tsv\\.gz$", all_files), "meta", "fg")
trait_by_file <- sub("_(meta_out|fg)\\.tsv\\.gz$", "", all_files)
present_dt <- data.table::data.table(phenotype = trait_by_file,
                                      file_type = file_type_by_trait)
# Defensive: if a trait somehow appears as both meta and fg, prefer meta.
present_dt <- present_dt[order(phenotype, factor(file_type, levels = c("meta", "fg")))]
present_dt <- unique(present_dt, by = "phenotype")
message(sprintf("Found %d sumstats files under %s (%d meta, %d fg)",
                length(all_files), sumstats_dir,
                sum(present_dt$file_type == "meta"),
                sum(present_dt$file_type == "fg")))

m <- data.table::fread(tsv_path)
required_cols <- c("fg_phenotype", "name", "category",
                   "fg_n_cases", "fg_n_controls",
                   "ukbb_n_cases", "ukbb_n_controls")
missing_cols <- setdiff(required_cols, names(m))
if (length(missing_cols) > 0L) {
  stop("Manifest missing column(s): ", paste(missing_cols, collapse = ", "))
}

# Restrict to phenotypes whose sumstats file is present, and attach file_type.
m <- m[fg_phenotype %in% present_dt$phenotype]
m <- merge(m, present_dt, by.x = "fg_phenotype", by.y = "phenotype",
            all.x = TRUE, sort = FALSE)
message(sprintf("After sumstats-presence filter: %d rows", nrow(m)))

m[, fg_n_cases      := suppressWarnings(as.numeric(fg_n_cases))]
m[, fg_n_controls   := suppressWarnings(as.numeric(fg_n_controls))]
m[, ukbb_n_cases    := suppressWarnings(as.numeric(ukbb_n_cases))]
m[, ukbb_n_controls := suppressWarnings(as.numeric(ukbb_n_controls))]

# Inclusion filter is file-type-aware:
#   meta -> require both FG and UKBB
#   fg   -> require only FG
fg_ok   <- !is.na(m$fg_n_cases)   & !is.na(m$fg_n_controls)   & m$fg_n_cases   > 0
ukbb_ok <- !is.na(m$ukbb_n_cases) & !is.na(m$ukbb_n_controls) & m$ukbb_n_cases > 0
keep <- (m$file_type == "meta" & fg_ok & ukbb_ok) |
        (m$file_type == "fg"   & fg_ok)
skipped <- m[!keep, .(fg_phenotype, file_type,
                      fg_ok = fg_ok[!keep], ukbb_ok = ukbb_ok[!keep])]
if (nrow(skipped) > 0L) {
  message(sprintf("Skipping %d trait(s) failing per-file-type inclusion:",
                  nrow(skipped)))
  for (i in seq_len(nrow(skipped))) {
    message(sprintf("  %s [%s] fg_ok=%s, ukbb_ok=%s",
                    skipped$fg_phenotype[i], skipped$file_type[i],
                    skipped$fg_ok[i], skipped$ukbb_ok[i]))
  }
}
m <- m[keep]

# Canonical chapter assignment using the extractor + fallback in utils.R.
src_utils <- file.path(project_root, "R", "utils.R")
source(src_utils)
chap <- extract_chapter_with_fallback(m$category, m$fg_phenotype)

# Pooled n_cases / n_controls: reflects what's ACTUALLY USED by the munge step.
# For meta traits: FG + UKBB (both contribute to the meta).
# For fg traits: FG only (the UKBB columns in the manifest TSV may be non-NA but
# UKBB sumstats are not part of the input file -- summing would be misleading).
# Used for logging / K computation; NOT for the scalar N passed to
# GenomicSEM::munge (which comes from per-SNP Neff median for meta files, or
# single-cohort Neff for fg files).
fg_cases   <- ifelse(is.na(m$fg_n_cases),    0, m$fg_n_cases)
fg_ctrls   <- ifelse(is.na(m$fg_n_controls), 0, m$fg_n_controls)
ukbb_cases <- ifelse(m$file_type == "fg", 0,
                     ifelse(is.na(m$ukbb_n_cases),    0, m$ukbb_n_cases))
ukbb_ctrls <- ifelse(m$file_type == "fg", 0,
                     ifelse(is.na(m$ukbb_n_controls), 0, m$ukbb_n_controls))

bothsex_meta_manifest <- data.frame(
  phenotype        = as.character(m$fg_phenotype),
  name             = as.character(m$name),
  category         = as.character(m$category),
  chapter          = chap$chapter,
  chapter_source   = chap$chapter_source,
  file_type        = as.character(m$file_type),
  # Per-cohort N columns keep the source-TSV naming (prefix style) so the
  # config$bothsex_meta$af_cohorts$manifest_n_cases / _n_controls keys read
  # them directly. DO NOT rename to suffix style -- it would silently break
  # per_snp_neff (lookup returns NULL -> all cohorts skipped -> trait dropped).
  fg_n_cases       = m$fg_n_cases,
  fg_n_controls    = m$fg_n_controls,
  ukbb_n_cases     = m$ukbb_n_cases,     # NA for fg-only traits
  ukbb_n_controls  = m$ukbb_n_controls,  # NA for fg-only traits
  n_cases          = fg_cases + ukbb_cases,
  n_controls       = fg_ctrls + ukbb_ctrls,
  stringsAsFactors = FALSE
)

# Per-trait log.
type_summary <- table(bothsex_meta_manifest$file_type)
message("File-type breakdown:")
for (t in names(type_summary)) {
  message(sprintf("  %-6s : %d", t, type_summary[[t]]))
}
src_summary <- table(bothsex_meta_manifest$chapter_source)
message("Chapter-source breakdown:")
for (s in names(src_summary)) {
  message(sprintf("  %-18s: %d", s, src_summary[[s]]))
}
unclassified <- bothsex_meta_manifest[
  bothsex_meta_manifest$chapter_source == "unclassified", "phenotype"]
if (length(unclassified) > 0L) {
  message(sprintf("Unclassified traits (will be dropped from CFA model by min_indicators): %s",
                  paste(unclassified, collapse = ", ")))
}

out_path <- file.path(config$paths$manifest_dir, "bothsex_meta_manifest.rda")
dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
save(bothsex_meta_manifest, file = out_path)
message(sprintf("wrote %s (%d rows)", out_path, nrow(bothsex_meta_manifest)))
