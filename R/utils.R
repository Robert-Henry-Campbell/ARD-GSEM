# --- Logging ---

.log_env <- new.env(parent = emptyenv())
.log_env$file <- NULL
.log_env$json_file <- NULL
.log_env$min_level <- "INFO"
.log_env$run_id <- NULL
.log_env$start_time <- NULL

LEVELS <- c(DEBUG = 0L, INFO = 1L, WARN = 2L, ERROR = 3L, FATAL = 4L)

init_logging <- function(config) {
  op <- options(digits.secs = 3)
  on.exit(options(op), add = TRUE)
  run_id <- format(Sys.time(), "%Y%m%d_%H%M%OS3")
  run_id <- gsub("\\.", "_", run_id)
  log_dir <- file.path(config$paths$output_dir, "logs")
  dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

  log_path <- file.path(log_dir, paste0("gsem_", run_id, ".log"))
  .log_env$file <- file(log_path, open = "wt")
  .log_env$min_level <- config$logging$level %||% "INFO"
  .log_env$run_id <- run_id
  .log_env$start_time <- Sys.time()

  if (isTRUE(config$logging$json)) {
    json_path <- file.path(log_dir, paste0("gsem_", run_id, ".json"))
    .log_env$json_file <- file(json_path, open = "wt")
  }

  if (.Platform$OS.type == "unix") {
    latest <- file.path(log_dir, "latest.log")
    if (file.exists(latest) || is.symlink(latest)) file.remove(latest)
    file.symlink(log_path, latest)
  }

  invisible(run_id)
}

is.symlink <- function(path) {
  info <- file.info(path)
  !is.na(info$isdir) && nchar(Sys.readlink(path)) > 0
}

should_log <- function(level, min_level = .log_env$min_level) {
  LEVELS[level] >= LEVELS[min_level]
}

format_log_message <- function(level, stage, msg) {
  ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  sprintf("%s [%s] [%s] %s", ts, level, stage, msg)
}

write_log <- function(level, stage, msg, data = NULL) {
  if (!should_log(level)) return(invisible())
  formatted <- format_log_message(level, stage, msg)

  if (isTRUE(.log_env$min_level == "DEBUG") || level %in% c("WARN", "ERROR", "FATAL")) {
    message(formatted)
  } else if (should_log(level)) {
    cat(formatted, "\n")
  }

  if (!is.null(.log_env$file)) {
    writeLines(formatted, .log_env$file)
    flush(.log_env$file)
  }

  if (!is.null(.log_env$json_file)) {
    entry <- list(
      ts = format(Sys.time(), "%Y-%m-%dT%H:%M:%S"),
      level = level,
      stage = stage,
      msg = msg
    )
    if (!is.null(data)) entry$data <- data
    writeLines(jsonlite::toJSON(entry, auto_unbox = TRUE), .log_env$json_file)
    flush(.log_env$json_file)
  }
}

log_debug <- function(stage, msg, data = NULL) write_log("DEBUG", stage, msg, data)
log_info  <- function(stage, msg, data = NULL) write_log("INFO",  stage, msg, data)
log_warn  <- function(stage, msg, data = NULL) write_log("WARN",  stage, msg, data)
log_error <- function(stage, msg, data = NULL) write_log("ERROR", stage, msg, data)
log_fatal <- function(stage, msg, data = NULL) {
  write_log("FATAL", stage, msg, data)
  stop(msg, call. = FALSE)
}

close_logging <- function() {
  if (!is.null(.log_env$start_time)) {
    elapsed <- as.numeric(difftime(Sys.time(), .log_env$start_time, units = "secs"))
    mins <- floor(elapsed / 60)
    secs <- round(elapsed %% 60)
    log_info("setup", sprintf("Pipeline finished. Total time: %dm %ds", mins, secs))
  }
  if (!is.null(.log_env$file)) { close(.log_env$file); .log_env$file <- NULL }
  if (!is.null(.log_env$json_file)) { close(.log_env$json_file); .log_env$json_file <- NULL }
}

# --- Config ---

is_absolute_path <- function(p) {
  startsWith(p, "/") || grepl("^[A-Za-z]:[/\\\\]", p) || startsWith(p, "\\\\")
}

read_config <- function(path = "config/pipeline.yaml", root = NULL) {
  if (!file.exists(path)) stop("Config not found: ", path)
  cfg <- yaml::read_yaml(path)
  effective_root <- root %||% cfg$project$root
  if (is.null(effective_root) || !nzchar(effective_root) ||
      !dir.exists(effective_root)) {
    effective_root <- normalizePath(dirname(dirname(path)),
                                    winslash = "/", mustWork = TRUE)
  }
  for (nm in names(cfg$paths)) {
    p <- cfg$paths[[nm]]
    if (!is_absolute_path(p)) cfg$paths[[nm]] <- file.path(effective_root, p)
  }
  cfg$project$root <- effective_root
  cfg
}

# --- Neff ---

compute_neff <- function(n_cases, n_controls) {
  bad <- n_cases <= 0 | n_controls <= 0 | is.na(n_cases) | is.na(n_controls)
  out <- 4 * n_cases * n_controls / (n_cases + n_controls)
  if (any(bad)) {
    log_warn("munge", sprintf("compute_neff: %d row(s) with non-positive/NA counts -> NA",
                              sum(bad)))
    out[bad] <- NA_real_
  }
  out
}

linear_to_logor <- function(beta, se, K,
                            reject_threshold = 1e-9,
                            warn_threshold = 0.01) {
  if (length(K) != 1L) stop("linear_to_logor: K must be a scalar (per-trait sample prevalence)")
  if (is.na(K) || K <= 0 || K >= 1) {
    stop(sprintf("linear_to_logor: invalid K=%s", as.character(K)))
  }
  scale <- K * (1 - K)
  if (scale < reject_threshold) {
    stop(sprintf("linear_to_logor: K(1-K)=%g below reject threshold %g; trait too imbalanced for stable conversion",
                 scale, reject_threshold))
  }
  warn <- (K < warn_threshold || K > 1 - warn_threshold)
  list(beta = beta / scale, se = se / scale, scale = scale, warn = warn)
}

harmonic_neff <- function(neffs) {
  neffs <- neffs[!is.na(neffs) & neffs > 0]
  if (length(neffs) == 0L) return(NA_real_)
  length(neffs) / sum(1 / neffs)
}

vech_diag_index <- function(k) {
  if (k < 1L) stop("vech_diag_index: k must be >= 1")
  if (k == 1L) return(1L)
  cumsum(c(1L, seq(as.integer(k), 2L, by = -1L)))
}

count_gz_lines <- function(path, chunk = 10000L) {
  con <- gzfile(path, "r")
  on.exit(close(con), add = TRUE)
  n <- 0L
  repeat {
    lines <- readLines(con, n = chunk, warn = FALSE)
    if (length(lines) == 0L) break
    n <- n + length(lines)
  }
  n
}

# --- Variant Parsing ---

parse_variant_column <- function(dt) {
  # Count colons per row first: tstrsplit pads short rows with NA up to the max field
  # count across the table, so a length-based check on its output would let a short
  # row slip through and fail later with a misleading "NA/empty" message.
  n_colons <- nchar(dt$variant) - nchar(gsub(":", "", dt$variant, fixed = TRUE))
  wrong_count <- n_colons != 3L
  if (any(wrong_count)) {
    stop(sprintf(
      "parse_variant_column: expected chr:pos:ref:alt (4 fields, 3 colons); %d row(s) had a different colon count. sample bad row: %s",
      sum(wrong_count), dt$variant[which(wrong_count)[1]]))
  }
  parts <- data.table::tstrsplit(dt$variant, ":", fixed = TRUE)
  empty_mask <- Reduce(`|`, lapply(parts, function(p) is.na(p) | !nzchar(p)))
  if (any(empty_mask)) {
    stop(sprintf(
      "parse_variant_column: %d row(s) have NA/empty fields; sample row: %s",
      sum(empty_mask), dt$variant[which(empty_mask)[1]]))
  }
  dt[, `:=`(chr = parts[[1]], pos = as.integer(parts[[2]]),
            ref = parts[[3]], alt = parts[[4]])]
  dt
}

# --- GRCh38 -> GRCh37 liftover (bothsex_meta) ---

# rtracklayer::liftOver returns a GRangesList: each input range maps to 0, 1,
# or N output ranges. A naive as.data.table() expansion silently misaligns rows
# with their associated effect estimates -- which is data corruption, not a
# fixable warning. The helper below tags each input row with .row_id, lifts,
# and keeps only rows whose .row_id maps exactly once. Multi-mapped and
# zero-mapped rows are dropped (counted separately for the log).
liftover_grch38_to_grch37 <- function(dt, chain_path) {
  if (!requireNamespace("rtracklayer", quietly = TRUE)) {
    log_fatal("munge", sprintf(
      "rtracklayer is required for liftOver but is not installed. Install with: Rscript -e 'BiocManager::install(\"rtracklayer\")'"))
  }
  if (!requireNamespace("GenomicRanges", quietly = TRUE)) {
    log_fatal("munge", "GenomicRanges (rtracklayer dependency) is missing")
  }
  if (is.null(chain_path) || !nzchar(chain_path) || !file.exists(chain_path)) {
    log_fatal("munge", sprintf("liftover chain file not found: %s", as.character(chain_path)))
  }
  n_in <- nrow(dt)
  if (n_in == 0L) return(dt)

  dt[, .row_id := .I]
  chr_in <- as.character(dt$chr)
  pos_in <- as.integer(dt$pos)
  gr <- GenomicRanges::GRanges(
    seqnames = paste0("chr", chr_in),
    ranges   = IRanges::IRanges(start = pos_in, end = pos_in),
    .row_id  = dt$.row_id
  )
  chain_obj <- rtracklayer::import.chain(chain_path)
  lifted_grl <- rtracklayer::liftOver(gr, chain_obj)
  lifted <- unlist(lifted_grl)

  row_ids <- S4Vectors::mcols(lifted)$.row_id
  tab <- table(row_ids)
  unique_ids <- as.integer(names(tab)[tab == 1L])
  multi_ids  <- as.integer(names(tab)[tab >  1L])
  n_multi <- length(multi_ids)
  n_unmap <- n_in - length(unique_ids) - n_multi

  keep_mask <- row_ids %in% unique_ids
  lifted_keep <- lifted[keep_mask]
  lookup <- data.table::data.table(
    .row_id = S4Vectors::mcols(lifted_keep)$.row_id,
    chr_new = sub("^chr", "", as.character(GenomicRanges::seqnames(lifted_keep))),
    pos_new = as.integer(GenomicRanges::start(lifted_keep))
  )
  data.table::setkey(lookup, .row_id)

  out <- merge(dt, lookup, by = ".row_id", all = FALSE, sort = FALSE)
  out[, chr := chr_new]
  out[, pos := pos_new]
  out[, c(".row_id", "chr_new", "pos_new") := NULL]
  dt[, .row_id := NULL]

  n_out <- nrow(out)
  pct <- if (n_in == 0L) 0 else 100 * n_out / n_in
  log_info("munge", sprintf(
    "liftover_grch38_to_grch37: %s -> %s variants (%.1f%% retained, %s multi-mapped, %s unmapped)",
    format(n_in, big.mark = ","), format(n_out, big.mark = ","), pct,
    format(n_multi, big.mark = ","), format(max(0L, n_unmap), big.mark = ",")))
  out
}

# --- Meta-analysis helpers (bothsex_meta) ---

# N-weighted mean of per-cohort effect-allele frequencies.
# af_mat: numeric matrix, rows = SNPs, cols = cohorts (in the order of
# config$bothsex_meta$af_cohorts).
# n_weights: numeric vector, length = ncol(af_mat); per-cohort trait-level
# (n_cases + n_controls) sum from the manifest.
# Returns: numeric vector length nrow(af_mat). For each row, NA AFs zero out
# their cohort's weight. Returns NA only if every cohort is NA on that row.
nweighted_meta_af <- function(af_mat, n_weights) {
  if (!is.matrix(af_mat)) af_mat <- as.matrix(af_mat)
  if (length(n_weights) != ncol(af_mat)) {
    stop(sprintf("nweighted_meta_af: n_weights length %d != ncol(af_mat) %d",
                 length(n_weights), ncol(af_mat)))
  }
  na_mask <- is.na(af_mat)
  af_filled <- af_mat
  af_filled[na_mask] <- 0
  w <- matrix(n_weights, nrow = nrow(af_mat), ncol = ncol(af_mat), byrow = TRUE)
  w[na_mask] <- 0
  num <- rowSums(af_filled * w)
  den <- rowSums(w)
  out <- ifelse(den > 0, num / den, NA_real_)
  n_all_na <- sum(den == 0)
  if (n_all_na > 0L) {
    log_info("munge", sprintf(
      "nweighted_meta_af: %s row(s) had no cohort AF (will be dropped by MAF filter)",
      format(n_all_na, big.mark = ",")))
  }
  out
}

# Per-SNP effective N for a multi-cohort meta-analysis.
# Each SNP's Neff = sum over cohorts of (cohort_Neff * I[cohort contributed]).
# Cohort contribution per row is inferred from non-NA beta_col.
# Returns: numeric vector length nrow(dt). Caller takes median() as the
# scalar trait-level N passed to GenomicSEM::munge().
per_snp_neff <- function(dt, af_cohorts, manifest_row) {
  if (length(af_cohorts) == 0L) {
    stop("per_snp_neff: af_cohorts is empty; check config$bothsex_meta$af_cohorts")
  }
  neff_per_snp <- numeric(nrow(dt))
  for (co in af_cohorts) {
    nca <- manifest_row[[co$manifest_n_cases]]
    nco <- manifest_row[[co$manifest_n_controls]]
    if (length(nca) == 0L || length(nco) == 0L ||
        is.na(nca) || is.na(nco) || nca <= 0 || nco <= 0) {
      next
    }
    cohort_neff <- 4 * nca * nco / (nca + nco)
    beta_col <- co$beta_col
    if (!beta_col %in% names(dt)) {
      log_warn("munge", sprintf(
        "per_snp_neff: cohort %s beta_col '%s' missing; skipping cohort",
        co$name, beta_col))
      next
    }
    present <- !is.na(dt[[beta_col]])
    neff_per_snp <- neff_per_snp + cohort_neff * as.numeric(present)
  }
  neff_per_snp
}

# --- Chapter assignment from FinnGen category (bothsex_meta) ---

# Canonical ICD chapter strings (from meta/icd10_categories.csv `chapter`
# column). Kept here as a constant so we don't re-read the CSV every time the
# extractor is called.
.canonical_icd_chapters <- c(
  "Certain infectious and parasitic diseases",
  "Neoplasms",
  "Diseases of the blood and blood-forming organs and certain disorders involving the immune mechanism",
  "Endocrine, nutritional and metabolic diseases",
  "Mental, Behavioral and Neurodevelopmental disorders",
  "Diseases of the nervous system",
  "Diseases of the eye and adnexa",
  "Diseases of the ear and mastoid process",
  "Diseases of the circulatory system",
  "Diseases of the respiratory system",
  "Diseases of the digestive system",
  "Diseases of the skin and subcutaneous tissue",
  "Diseases of the musculoskeletal system and connective tissue",
  "Diseases of the genitourinary system",
  "Pregnancy, childbirth and the puerperium",
  "Certain conditions originating in the perinatal period",
  "Congenital malformations, deformations and chromosomal abnormalities",
  "Symptoms, signs and abnormal clinical and laboratory findings, not elsewhere classified",
  "Injury, poisoning and certain other consequences of external causes",
  "External causes of morbidity",
  "Factors influencing health status and contact with health services"
)

# FinnGen endpoint-prefix -> canonical chapter. Used as fallback when the
# manifest `category` column is non-ICD (e.g. "Comorbidities of COPD",
# "Interstitial lung disease endpoints"). Keys match `^[A-Z]+[0-9]+_` parsed
# from the trait ID.
.finngen_prefix_chapter <- list(
  AB1  = "Certain infectious and parasitic diseases",
  CD2  = "Neoplasms",
  C3   = "Neoplasms",
  D3   = "Diseases of the blood and blood-forming organs and certain disorders involving the immune mechanism",
  E4   = "Endocrine, nutritional and metabolic diseases",
  F5   = "Mental, Behavioral and Neurodevelopmental disorders",
  G6   = "Diseases of the nervous system",
  H7   = "Diseases of the eye and adnexa",
  H8   = "Diseases of the ear and mastoid process",
  I9   = "Diseases of the circulatory system",
  J10  = "Diseases of the respiratory system",
  K11  = "Diseases of the digestive system",
  L12  = "Diseases of the skin and subcutaneous tissue",
  M13  = "Diseases of the musculoskeletal system and connective tissue",
  N14  = "Diseases of the genitourinary system",
  Q17  = "Congenital malformations, deformations and chromosomal abnormalities",
  R18  = "Symptoms, signs and abnormal clinical and laboratory findings, not elsewhere classified",
  ST19 = "Injury, poisoning and certain other consequences of external causes",
  Z21  = "Factors influencing health status and contact with health services"
)

# Strip Roman-numeral chapter prefix + trailing parenthetical, trim.
# Vectorised.
.clean_finngen_category <- function(s) {
  s <- as.character(s)
  s <- sub("^(I{1,3}|IV|VI{0,3}|IX|X{1,3}|XI{1,3}|XIV|XV|XVI{0,3}|XIX|XX|XXI)\\s+",
           "", s, perl = TRUE)
  s <- sub("\\s*\\([^)]*\\)\\s*$", "", s, perl = TRUE)
  trimws(s)
}

# For each (category_text, trait_id) pair, return a list with:
#   chapter:        canonical chapter string (or "Unclassified")
#   chapter_source: "category-table" | "prefix-fallback" | "unclassified"
# Vectorised. trait_ids must be same length as category_text.
extract_chapter_with_fallback <- function(category_text, trait_ids) {
  if (length(category_text) != length(trait_ids)) {
    stop("extract_chapter_with_fallback: category_text and trait_ids must be same length")
  }
  cleaned <- .clean_finngen_category(category_text)
  canon_lc <- tolower(.canonical_icd_chapters)
  cleaned_lc <- tolower(cleaned)
  chapter <- character(length(trait_ids))
  source  <- character(length(trait_ids))
  for (i in seq_along(trait_ids)) {
    hit <- match(cleaned_lc[i], canon_lc)
    if (!is.na(hit)) {
      chapter[i] <- .canonical_icd_chapters[hit]
      source[i]  <- "category-table"
      next
    }
    pref <- regmatches(trait_ids[i],
                       regexpr("^[A-Z]+[0-9]+", trait_ids[i], perl = TRUE))
    if (length(pref) == 1L && nzchar(pref) && !is.null(.finngen_prefix_chapter[[pref]])) {
      chapter[i] <- .finngen_prefix_chapter[[pref]]
      source[i]  <- "prefix-fallback"
      next
    }
    chapter[i] <- "Unclassified"
    source[i]  <- "unclassified"
  }
  list(chapter = chapter, chapter_source = source)
}

# Source-aware trait-categories loader. Replaces all direct reads of
# meta/icd10_categories.csv so the CFA/LDSC/GWAS/plots stages can stay
# source-agnostic. Returned data.table schema: code, title, chapter_range,
# chapter (matching the current icd10_categories.csv schema).
load_trait_categories <- function(config, sex) {
  if (identical(sex, "bothsex_meta")) {
    rda <- config$bothsex_meta$manifest_rda %||%
           file.path(config$paths$manifest_dir, "bothsex_meta_manifest.rda")
    if (!is_absolute_path(rda)) {
      rda <- file.path(config$project$root %||% getwd(), rda)
    }
    if (!file.exists(rda)) {
      log_fatal("setup", sprintf("bothsex_meta manifest .rda not found: %s. Run scripts/generate_bothsex_meta_manifest.R.", rda))
    }
    env <- new.env()
    load(rda, envir = env)
    objs <- ls(env)
    if (length(objs) != 1L) {
      stop(sprintf("load_trait_categories: expected exactly 1 object in %s, got %d",
                   basename(rda), length(objs)))
    }
    m <- as.data.table(get(objs[1], envir = env))
    prefix <- regmatches(m$phenotype,
                         regexpr("^[A-Z]+[0-9]+_?", m$phenotype, perl = TRUE))
    if (length(prefix) != nrow(m)) {
      # regmatches with vector input returns only matched entries; pad NAs for
      # rows without a parseable prefix
      prefix_full <- character(nrow(m))
      raw <- regexpr("^[A-Z]+[0-9]+_?", m$phenotype, perl = TRUE)
      for (i in seq_len(nrow(m))) {
        if (raw[i] > 0L) {
          prefix_full[i] <- substr(m$phenotype[i], 1L,
                                   raw[i] + attr(raw, "match.length")[i] - 1L)
        } else {
          prefix_full[i] <- NA_character_
        }
      }
      prefix <- prefix_full
    }
    return(data.table(
      code          = m$phenotype,
      title         = if ("name" %in% names(m)) m$name else m$phenotype,
      chapter_range = prefix,
      chapter       = m$chapter
    ))
  }
  fread(file.path(config$paths$meta_dir, "icd10_categories.csv"))
}

# --- rsID Mapping ---

load_variants_manifest <- function(path) {
  log_info("munge", "Loading variants manifest...")
  vt <- data.table::fread(
    cmd = paste("zcat", shQuote(path)),
    select = c("variant", "rsid"),
    key = "variant"
  )
  log_debug("munge", sprintf("Variants manifest: %s rows loaded", format(nrow(vt), big.mark = ",")))
  vt
}

map_rsids <- function(sumstats, manifest) {
  n_before <- nrow(sumstats)
  result <- manifest[sumstats, on = "variant", nomatch = NULL]
  result[, SNP := rsid]
  result[, rsid := NULL]
  pct <- if (n_before == 0L) 0 else 100 * nrow(result) / n_before
  log_info("munge",
           sprintf("rsid map: %s -> %s variants (%.1f%% mapped)",
                   format(n_before, big.mark = ","),
                   format(nrow(result), big.mark = ","),
                   pct))
  if (n_before > 0L && pct < 50) {
    log_warn("munge",
             sprintf("rsid map: only %.1f%% of variants mapped; check variants_manifest", pct))
  }
  result
}

# --- Pan-UKB variants manifest ---

# Merge a user-supplied column map (from YAML) with defaults, validate that
# every required key is present after the merge, and return a list keyed by
# internal-name where values are the source-column names to read from disk.
# Robust to YAML key reorderings and partial overrides (the original
# unlist+positional-rename approach silently misaligned columns).
resolve_column_map <- function(user_map, defaults, label = "column map") {
  out <- defaults
  if (!is.null(user_map)) {
    for (k in names(user_map)) {
      v <- user_map[[k]]
      if (!is.null(v) && nzchar(v)) out[[k]] <- v
    }
  }
  missing <- setdiff(names(defaults), names(out))
  if (length(missing) > 0L) {
    stop(sprintf("%s missing key(s): %s", label, paste(missing, collapse = ",")),
         call. = FALSE)
  }
  out
}

load_panukb_variants_manifest <- function(path, cols = NULL) {
  log_info("munge", "Loading Pan-UKB variants manifest...")
  defaults <- list(chr = "chrom", pos = "pos", ref = "ref", alt = "alt",
                   rsid = "rsid", info = "info", af = "af_EUR")
  internal <- list(chr = "chr", pos = "pos", ref = "ref", alt = "alt",
                   rsid = "rsid", info = "info", af = "af_EUR")
  cm <- resolve_column_map(cols, defaults, label = "variants_columns")
  source_names <- vapply(names(internal), function(k) cm[[k]], character(1))
  target_names <- vapply(names(internal), function(k) internal[[k]], character(1))
  classes <- list()
  classes[[cm$chr]]  <- "character"
  classes[[cm$rsid]] <- "character"
  classes[[cm$info]] <- "numeric"
  classes[[cm$af]]   <- "numeric"
  vt <- data.table::fread(
    cmd = paste("zcat", shQuote(path)),
    select = unname(source_names), colClasses = classes,
    na.strings = c("", "NA")
  )
  data.table::setnames(vt, old = unname(source_names), new = unname(target_names))
  data.table::setkey(vt, chr, pos, ref, alt)
  log_debug("munge", sprintf("Pan-UKB variants manifest: %s rows loaded",
                              format(nrow(vt), big.mark = ",")))
  vt
}

merge_panukb_variants <- function(dt, vm) {
  n_before <- nrow(dt)
  dt[, chr := as.character(chr)]
  result <- vm[dt, on = c("chr", "pos", "ref", "alt"), nomatch = NULL]
  result <- result[!is.na(rsid) & nzchar(rsid)]
  pct <- if (n_before == 0L) 0 else 100 * nrow(result) / n_before
  log_info("munge",
           sprintf("Pan-UKB variant join: %s -> %s variants (%.1f%% matched + rsid'd)",
                   format(n_before, big.mark = ","),
                   format(nrow(result), big.mark = ","), pct))
  result
}

# --- Filtering ---

filter_variants <- function(dt, maf_threshold = 0.01, info_threshold = NULL) {
  dt <- dt[low_confidence_variant == FALSE]
  dt <- dt[minor_AF >= maf_threshold]
  # info_threshold preserved for non-Neale inputs that carry an 'info' column.
  # Neale UKBB sumstats do not have one; the branch is a no-op there.
  if (!is.null(info_threshold) && "info" %in% names(dt)) {
    dt <- dt[info >= info_threshold]
  }
  dt <- dt[!is.na(beta) & !is.na(se)]
  dt
}

# --- Polarity Check ---

check_polarity_sentinel <- function(dt, sentinel_rsid = NULL, expected_a1 = NULL) {
  # Returns list with polarity_correct (logical/NA), message (chr), and
  # sentinel_missing (logical). The sentinel_missing flag lets callers escalate
  # the "not found" case to WARN (e.g. post-liftover for bothsex_meta) instead
  # of accepting the silent DEBUG default.
  if (is.null(sentinel_rsid) || is.null(expected_a1)) {
    return(list(polarity_correct = NA,
                message = "polarity check skipped (no sentinel configured)",
                sentinel_missing = FALSE))
  }
  row <- dt[SNP == sentinel_rsid]
  if (nrow(row) == 0) {
    return(list(polarity_correct = NA,
                message = sprintf("Sentinel SNP %s not found", sentinel_rsid),
                sentinel_missing = TRUE))
  }
  correct <- row$A1[1] == expected_a1
  list(polarity_correct = correct,
       message = sprintf("Sentinel %s: A1=%s (expected %s)",
                         sentinel_rsid, row$A1[1], expected_a1),
       sentinel_missing = FALSE)
}

# --- A Priori Model Building ---

build_apriori_model <- function(traits, categories, min_indicators = 2,
                                ordering_table = NULL,
                                add_factor_covariances = TRUE,
                                add_heywood_constraints = TRUE,
                                code_format = c("icd3", "free")) {
  code_format <- match.arg(code_format)
  if (identical(code_format, "icd3")) {
    bad_nchar <- nchar(traits) != 3L
    if (any(bad_nchar)) {
      log_warn("cfa", sprintf(
        "build_apriori_model: dropping %d non-3-char trait(s): %s",
        sum(bad_nchar), paste(traits[bad_nchar], collapse = ", ")))
      traits <- traits[!bad_nchar]
    }
  }
  if (length(traits) == 0L) return("")
  mapping <- categories[code %in% traits, .(code, chapter)]
  unmapped <- setdiff(traits, mapping$code)
  if (length(unmapped) > 0) {
    log_warn("cfa", sprintf("build_apriori_model: %d trait(s) not in icd10_categories (no fallback): %s",
                             length(unmapped), paste(unmapped, collapse = ", ")))
  }
  if (nrow(mapping) == 0) return("")

  order_codes <- function(codes) {
    if (!is.null(ordering_table)) {
      ord <- as.data.table(ordering_table)
      hit <- ord[code %in% codes]
      if (nrow(hit) > 0L) {
        hit <- hit[order(-abs(sort_metric))]
        ordered <- hit$code
        # Append any codes missing from ordering_table at the tail in original order.
        missing <- setdiff(codes, ordered)
        return(c(ordered, missing))
      }
    }
    sort(codes)
  }

  by_factor <- split(mapping$code, mapping$chapter)
  by_factor <- by_factor[vapply(by_factor, length, integer(1)) >= min_indicators]

  if (length(by_factor) == 0) return("")

  by_factor <- lapply(by_factor, order_codes)

  factor_names <- vapply(names(by_factor), function(ch) {
    fn <- gsub("[^A-Za-z0-9]", "_", ch)
    substr(fn, 1, 20)
  }, character(1))
  stopifnot(!any(duplicated(factor_names)))

  loading_lines <- mapply(function(fn, codes) {
    paste0(fn, " =~ ", paste(codes, collapse = " + "))
  }, factor_names, by_factor, USE.NAMES = FALSE)

  parts <- loading_lines

  if (isTRUE(add_factor_covariances) && length(factor_names) >= 2L) {
    pairs <- combn(factor_names, 2L)
    cov_lines <- apply(pairs, 2, function(p) sprintf("%s ~~ %s", p[1], p[2]))
    parts <- c(parts, cov_lines)
  }

  if (isTRUE(add_heywood_constraints)) {
    all_codes <- unique(unlist(by_factor, use.names = FALSE))
    rv_lines <- sprintf("%s ~~ rv_%s*%s", all_codes, all_codes, all_codes)
    cn_lines <- sprintf("rv_%s > 0.001", all_codes)
    parts <- c(parts, rv_lines, cn_lines)
  }

  paste(parts, collapse = "\n")
}

# --- Chapter-survival pre-check ---

check_chapter_survival <- function(retained, categories,
                                   min_indicators = 3L,
                                   min_factors_required = 3L) {
  cats <- as.data.table(categories)
  mapping <- cats[code %in% retained, .(code, chapter)]
  by_chapter <- split(mapping$code, mapping$chapter)
  counts <- vapply(by_chapter, length, integer(1))
  surviving <- names(counts)[counts >= min_indicators]
  dropped <- names(counts)[counts < min_indicators]

  log_info("ldsc", sprintf(
    "check_chapter_survival: %d/%d chapters have >=%d retained indicators",
    length(surviving), length(counts), min_indicators))

  for (ch in dropped) {
    log_warn("ldsc", sprintf(
      "  chapter dropped: '%s' (only %d retained: %s)",
      ch, counts[[ch]], paste(by_chapter[[ch]], collapse = ", ")))
  }

  if (length(surviving) < min_factors_required) {
    log_fatal("ldsc", sprintf(
      "check_chapter_survival: only %d chapters survive at min_indicators=%d (need >=%d). Lower the h2_z_threshold or reduce min_indicators_per_factor.",
      length(surviving), min_indicators, min_factors_required))
  }

  list(n_surviving = length(surviving),
       n_dropped = length(dropped),
       surviving_chapters = surviving,
       dropped_chapters = dropped,
       counts = counts)
}

# --- Genome-build coordinate sanity check ---

verify_genome_build <- function(dt, expected, build_label = "GRCh37") {
  if (length(expected) == 0L) {
    log_debug("munge", "verify_genome_build: no expected SNPs configured; skipping")
    return(invisible(TRUE))
  }
  rsids <- vapply(expected, function(e) e$rsid, character(1))
  hits <- dt[SNP %in% rsids, .(SNP, chr, pos)]
  mismatches <- character(0)
  for (e in expected) {
    row <- hits[SNP == e$rsid]
    if (nrow(row) == 0L) next
    if (as.character(row$chr[1]) != as.character(e$chr) ||
        as.integer(row$pos[1]) != as.integer(e$pos)) {
      mismatches <- c(mismatches, sprintf(
        "%s: expected %s:%d, got %s:%d", e$rsid, e$chr, e$pos,
        row$chr[1], row$pos[1]))
    }
  }
  if (length(mismatches) > 0L) {
    log_fatal("munge", sprintf(
      "verify_genome_build (%s): %d coordinate mismatch(es) -- %s",
      build_label, length(mismatches), paste(mismatches, collapse = "; ")))
  }
  log_debug("munge", sprintf("verify_genome_build (%s): all %d seed SNP(s) at expected coordinates",
                             build_label, nrow(hits)))
  invisible(TRUE)
}

# --- Standardized loading extraction from usermodel() ---

extract_apriori_loadings <- function(apriori_fit, categories) {
  sol <- apriori_fit$results
  if (is.null(sol) || nrow(sol) == 0L) return(NULL)
  sol <- as.data.table(sol)
  std_col <- intersect(c("STD_All", "STD_Genotype", "Standardized_Est"), names(sol))[1]
  if (is.na(std_col)) return(NULL)
  loadings <- sol[op == "=~", .(factor = lhs, code = rhs,
                                std_loading = as.numeric(get(std_col)))]
  if (nrow(loadings) == 0L) return(NULL)
  cats <- as.data.table(categories)[, .(code, chapter)]
  out <- merge(loadings, cats, by = "code", all.x = TRUE, sort = FALSE)
  out[, sort_metric := std_loading]
  out[, .(code, chapter, factor, std_loading, sort_metric)]
}

# --- Factor effective sample size (Wiki 4 N_hat) ---

compute_factor_nhat <- function(maf, se, maf_threshold = 0.10) {
  ok <- !is.na(maf) & !is.na(se) & maf >= maf_threshold &
        maf <= (1 - maf_threshold) & se > 0
  if (!any(ok)) return(NA_real_)
  mean(1 / ((2 * maf[ok] * (1 - maf[ok])) * se[ok]^2))
}

# --- EFA to lavaan ---

efa_to_lavaan <- function(loadings_matrix, threshold = 0.30) {
  n_factors <- ncol(loadings_matrix)
  trait_names <- rownames(loadings_matrix)
  lines <- character(n_factors)
  for (f in seq_len(n_factors)) {
    indicators <- trait_names[abs(loadings_matrix[, f]) >= threshold]
    if (length(indicators) > 0) {
      lines[f] <- paste0("F", f, " =~ ", paste(indicators, collapse = " + "))
    }
  }
  paste(lines[nchar(lines) > 0], collapse = "\n")
}

# --- S Matrix Checks ---

check_psd <- function(S) {
  eigenvalues <- eigen(S, symmetric = TRUE, only.values = TRUE)$values
  all(eigenvalues >= -1e-10)
}

# --- h2 Filtering ---

filter_h2 <- function(h2_table, z_threshold = 2.0) {
  if (!"z" %in% names(h2_table)) {
    h2_table[, z := h2 / se]
  }
  h2_table[z >= z_threshold]
}

# --- Loading Differences ---

compute_loading_diff <- function(male_loading, male_se, female_loading, female_se) {
  diff <- female_loading - male_loading
  se_diff <- sqrt(male_se^2 + female_se^2)
  z_diff <- diff / se_diff
  p_diff <- 2 * pnorm(-abs(z_diff))
  list(diff = diff, se_diff = se_diff, z_diff = z_diff, p_diff = p_diff)
}

# --- Stage Manifests ---

write_stage_manifest <- function(stage, sex, config, output_files, warnings = character(0)) {
  # Hash includes source-specific config (panukb / bothsex_meta) so editing the
  # source block correctly invalidates --resume on the next run. should_skip()
  # in run_pipeline.R falls back to compat / legacy hashes for older manifests.
  source_key <- switch(sex,
                       bothsex      = "panukb",
                       bothsex_meta = "bothsex_meta",
                       NULL)
  source_cfg <- if (is.null(source_key)) NULL else config[[source_key]]
  manifest <- list(
    stage = stage,
    sex = sex,
    completed_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S"),
    run_id = .log_env$run_id %||% "unknown",
    config_hash = digest::digest(
      list(stage = config[[stage]] %||% list(),
           paths = config$paths,
           source = source_cfg),
      algo = "sha256"),
    output_files = basename(output_files),
    traits_processed = length(output_files),
    warnings = warnings
  )
  out_dir <- file.path(config$paths$output_dir, sex, stage)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  jsonlite::write_json(manifest, file.path(out_dir, ".manifest.json"),
                       auto_unbox = TRUE, pretty = TRUE)
}

stage_is_fresh <- function(manifest, current_config_hash) {
  if (is.null(manifest)) return(FALSE)
  identical(manifest$config_hash, current_config_hash)
}

read_stage_manifest <- function(stage, sex, config) {
  path <- file.path(config$paths$output_dir, sex, stage, ".manifest.json")
  if (!file.exists(path)) return(NULL)
  jsonlite::read_json(path)
}

# --- Trait Discovery ---

discover_traits <- function(config, sex) {
  dir_path <- file.path(config$paths$sumstats_dir, sex)
  if (identical(sex, "bothsex")) {
    pat <- config$discover$bothsex_filename_pattern %||%
           "^icd10-.*-both_sexes\\.tsv\\.bgz$"
    files <- list.files(dir_path, pattern = pat, full.names = FALSE)
    return(sub("^icd10-(.*)-both_sexes\\.tsv\\.bgz$", "\\1", files))
  }
  if (identical(sex, "bothsex_meta")) {
    # Bothsex_meta dir mixes two file types from the same FinnGen release:
    # <TRAIT>_meta_out.tsv.gz (2-way FG+UKBB meta) and <TRAIT>_fg.tsv.gz
    # (FG-only single-cohort GWAS). Strip whichever suffix to recover trait ID.
    # See discover_traits_with_type() for the (trait, file_type) shape used by
    # run_munge_bothsex_meta to dispatch per-file.
    pat <- config$discover$bothsex_meta_filename_pattern %||%
           "^[A-Z][A-Z0-9_]+_(meta_out|fg)\\.tsv\\.gz$"
    files <- list.files(dir_path, pattern = pat, full.names = FALSE)
    return(sub("_(meta_out|fg)\\.tsv\\.gz$", "", files))
  }
  # Filename pattern is config-driven so non-Neale inputs can be plugged in.
  # Stored under config$discover (not config$munge) so the munge stage hash is
  # unaffected by adding/changing this key -- prevents stale-manifest invalidation.
  pat <- config$discover$sumstats_filename_pattern %||% "\\.gwas\\.imputed_v3\\."
  files <- list.files(dir_path, pattern = pat, full.names = FALSE)
  sub(paste0(pat, ".*"), "", files)
}

# bothsex_meta-only helper. Returns data.table(trait, file_type, file_path)
# where file_type is "meta" (file ends in _meta_out.tsv.gz, i.e. 2-way
# FG+UKBB inverse-variance meta) or "fg" (file ends in _fg.tsv.gz, i.e. raw
# single-cohort FinnGen GWAS). The two file types have different schemas and
# need different munge code paths (see run_munge_bothsex_meta).
discover_traits_with_type <- function(config) {
  dir_path <- file.path(config$paths$sumstats_dir, "bothsex_meta")
  pat <- config$discover$bothsex_meta_filename_pattern %||%
         "^[A-Z][A-Z0-9_]+_(meta_out|fg)\\.tsv\\.gz$"
  files <- list.files(dir_path, pattern = pat, full.names = FALSE)
  if (length(files) == 0L) {
    return(data.table::data.table(trait = character(0),
                                   file_type = character(0),
                                   file_path = character(0)))
  }
  is_meta <- grepl("_meta_out\\.tsv\\.gz$", files)
  is_fg   <- grepl("_fg\\.tsv\\.gz$", files)
  file_type <- ifelse(is_meta, "meta", ifelse(is_fg, "fg", NA_character_))
  trait <- sub("_(meta_out|fg)\\.tsv\\.gz$", "", files)
  data.table::data.table(trait     = trait,
                          file_type = file_type,
                          file_path = file.path(dir_path, files))
}

get_shared_traits <- function(config) {
  male <- discover_traits(config, "male")
  female <- discover_traits(config, "female")
  intersect(male, female)
}

# --- Case/Control Lookup ---

get_case_control <- function(config, sex, trait) {
  manifest_file <- if (identical(sex, "bothsex")) {
    file.path(config$paths$manifest_dir, "panukb_bothsex_manifest.rda")
  } else if (identical(sex, "bothsex_meta")) {
    file.path(config$paths$manifest_dir, "bothsex_meta_manifest.rda")
  } else {
    file.path(config$paths$manifest_dir, paste0("neale_", sex, "_manifest.rda"))
  }
  env <- new.env()
  load(manifest_file, envir = env)
  objs <- ls(env)
  if (length(objs) != 1L) {
    stop(sprintf(
      "get_case_control: expected exactly 1 object in %s, found %d (%s). Refusing to guess.",
      basename(manifest_file), length(objs), paste(objs, collapse = ", ")))
  }
  manifest <- get(objs[1], envir = env)
  row <- manifest[manifest$phenotype == trait, ]
  if (nrow(row) == 0) stop(sprintf("Trait %s not found in %s manifest", trait, sex))
  # For bothsex_meta, also return the full row so callers can compute per-SNP
  # Neff using per-cohort counts. row[1,] is a single-row data.frame.
  out <- list(n_cases = row$n_cases[1], n_controls = row$n_controls[1])
  if (identical(sex, "bothsex_meta")) {
    out$manifest_row <- as.list(row[1, , drop = FALSE])
  }
  out
}
