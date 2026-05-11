# --- Logging ---

.log_env <- new.env(parent = emptyenv())
.log_env$file <- NULL
.log_env$json_file <- NULL
.log_env$min_level <- "INFO"
.log_env$run_id <- NULL
.log_env$start_time <- NULL

LEVELS <- c(DEBUG = 0L, INFO = 1L, WARN = 2L, ERROR = 3L, FATAL = 4L)

init_logging <- function(config) {
  run_id <- format(Sys.time(), "%Y%m%d_%H%M%S")
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

read_config <- function(path = "config/pipeline.yaml") {
  if (!file.exists(path)) stop("Config not found: ", path)
  cfg <- yaml::read_yaml(path)
  root <- cfg$project$root
  for (nm in names(cfg$paths)) {
    p <- cfg$paths[[nm]]
    if (!startsWith(p, "/")) cfg$paths[[nm]] <- file.path(root, p)
  }
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

linear_to_logor <- function(beta, se, K) {
  if (length(K) != 1L) stop("linear_to_logor: K must be a scalar (per-trait sample prevalence)")
  if (is.na(K) || K <= 0 || K >= 1) {
    stop(sprintf("linear_to_logor: invalid K=%s", as.character(K)))
  }
  scale <- K * (1 - K)
  if (scale < 1e-6) {
    stop(sprintf("linear_to_logor: K(1-K)=%g below 1e-6; trait too imbalanced for stable conversion", scale))
  }
  list(beta = beta / scale, se = se / scale, scale = scale)
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
  parts <- data.table::tstrsplit(dt$variant, ":", fixed = TRUE)
  dt[, `:=`(chr = parts[[1]], pos = as.integer(parts[[2]]),
            ref = parts[[3]], alt = parts[[4]])]
  dt
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
  result <- manifest[sumstats, on = "variant", nomatch = NULL]
  result[, SNP := rsid]
  result[, rsid := NULL]
  result
}

# --- Filtering ---

filter_variants <- function(dt, maf_threshold = 0.01, info_threshold = NULL) {
  dt <- dt[low_confidence_variant == FALSE]
  dt <- dt[minor_AF >= maf_threshold]
  if (!is.null(info_threshold) && "info" %in% names(dt)) {
    dt <- dt[info >= info_threshold]
  }
  dt <- dt[!is.na(beta) & !is.na(se)]
  dt
}

# --- Polarity Check ---

check_polarity_sentinel <- function(dt, sentinel_rsid = NULL, expected_a1 = NULL) {
  if (is.null(sentinel_rsid) || is.null(expected_a1)) {
    return(list(polarity_correct = NA,
                message = "polarity check skipped (no sentinel configured)"))
  }
  row <- dt[SNP == sentinel_rsid]
  if (nrow(row) == 0) {
    return(list(polarity_correct = NA,
                message = sprintf("Sentinel SNP %s not found", sentinel_rsid)))
  }
  correct <- row$A1[1] == expected_a1
  list(polarity_correct = correct,
       message = sprintf("Sentinel %s: A1=%s (expected %s)",
                         sentinel_rsid, row$A1[1], expected_a1))
}

# --- A Priori Model Building ---

build_apriori_model <- function(traits, categories, min_indicators = 2) {
  stopifnot(all(nchar(traits) == 3L))
  mapping <- categories[code %in% traits, .(code, chapter)]
  unmapped <- setdiff(traits, mapping$code)
  if (length(unmapped) > 0) {
    log_warn("cfa", sprintf("build_apriori_model: %d trait(s) not in icd10_categories (no fallback): %s",
                             length(unmapped), paste(unmapped, collapse = ", ")))
  }
  if (nrow(mapping) == 0) return("")

  by_factor <- split(mapping$code, mapping$chapter)
  by_factor <- by_factor[vapply(by_factor, length, integer(1)) >= min_indicators]

  if (length(by_factor) == 0) return("")

  factor_names <- vapply(names(by_factor), function(ch) {
    fn <- gsub("[^A-Za-z0-9]", "_", ch)
    substr(fn, 1, 20)
  }, character(1))
  stopifnot(!any(duplicated(factor_names)))

  lines <- mapply(function(fn, codes) {
    paste0(fn, " =~ ", paste(codes, collapse = " + "))
  }, factor_names, by_factor, USE.NAMES = FALSE)

  paste(lines, collapse = "\n")
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
  manifest <- list(
    stage = stage,
    sex = sex,
    completed_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S"),
    run_id = .log_env$run_id %||% "unknown",
    config_hash = digest::digest(config[[stage]] %||% config, algo = "sha256"),
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
  files <- list.files(dir_path, pattern = "\\.gwas\\.imputed_v3\\.", full.names = FALSE)
  sub("\\.gwas\\.imputed_v3\\..*", "", files)
}

get_shared_traits <- function(config) {
  male <- discover_traits(config, "male")
  female <- discover_traits(config, "female")
  intersect(male, female)
}

# --- Case/Control Lookup ---

get_case_control <- function(config, sex, trait) {
  manifest_file <- file.path(config$paths$manifest_dir,
                             paste0("neale_", sex, "_manifest.rda"))
  env <- new.env()
  load(manifest_file, envir = env)
  manifest_name <- ls(env)[1]
  manifest <- get(manifest_name, envir = env)
  row <- manifest[manifest$phenotype == trait, ]
  if (nrow(row) == 0) stop(sprintf("Trait %s not found in %s manifest", trait, sex))
  list(n_cases = row$n_cases[1], n_controls = row$n_controls[1])
}
