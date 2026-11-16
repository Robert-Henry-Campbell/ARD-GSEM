# Shared helpers for the post-GWAS stages (qc_filters onwards). Every post-GWAS stage
# reads the per-factor artefacts written once by qc_filters (output/<sex>/qc_filters/
# <factor>.qc.rds) instead of the 2.7 GB raw userGWAS RDS.

stage_dir <- function(config, sex, stage) {
  file.path(config$paths$output_dir, sex, stage)
}

# Factors analysed after the GWAS: config$postgwas$factors minus the excluded ones, in
# the order given in config.
postgwas_factors <- function(config) {
  f <- unlist(config$postgwas$factors %||% character(0))
  setdiff(f, names(config$postgwas$excluded_factors %||% list()))
}

in_mhc <- function(chr, bp, config) {
  m <- config$postgwas$mhc %||% list(chr = 6, start = 25e6, end = 35e6)
  as.integer(chr) == as.integer(m$chr) & bp >= m$start & bp <= m$end
}

# -log10(P) for a two-sided normal test, computed in log space so it survives P < 1e-308.
log10p_from_z <- function(z) {
  -(pnorm(-abs(z), log.p = TRUE) + log(2)) / log(10)
}

# -log10 of a chi-square upper-tail p-value, in log space.
log10p_chisq <- function(q, df) {
  -pchisq(q, df, lower.tail = FALSE, log.p = TRUE) / log(10)
}

# --- canonical per-factor artefacts --------------------------------------------------------

qc_rds_path <- function(config, sex, factor) {
  file.path(stage_dir(config, sex, "qc_filters"), paste0(factor, ".qc.rds"))
}

.postgwas_cache <- new.env(parent = emptyenv())

# Canonical QC'd sumstats for one factor. tier: "all" (every SNP, with flags), "full"
# (converged, smoothing OK) or "qfilter" (full minus Q_SNP hits and their LD partners).
load_postgwas <- function(config, sex, factor, tier = c("all", "full", "qfilter"), cols = NULL) {
  tier <- match.arg(tier)
  path <- qc_rds_path(config, sex, factor)
  key <- paste(path, file.mtime(path))
  dt <- .postgwas_cache[[key]]
  if (is.null(dt)) {
    if (!file.exists(path)) {
      log_fatal("postgwas", sprintf("%s not found (run --stage qc_filters first)", path))
    }
    dt <- readRDS(path)
    .postgwas_cache[[key]] <- dt
  }
  keep <- switch(tier, all = rep(TRUE, nrow(dt)), full = dt$pass_full, qfilter = dt$pass_qfilter)
  if (is.null(cols)) dt[keep] else dt[keep, ..cols]
}

clear_postgwas_cache <- function() rm(list = ls(.postgwas_cache), envir = .postgwas_cache)

factor_registry_path <- function(config, sex) {
  file.path(stage_dir(config, sex, "qc_filters"), "factor_registry.json")
}

read_factor_registry <- function(config, sex) {
  path <- factor_registry_path(config, sex)
  if (!file.exists(path)) log_fatal("postgwas", sprintf("%s not found (run qc_filters)", path))
  jsonlite::read_json(path, simplifyVector = TRUE)
}

# --- external tools ------------------------------------------------------------------------

# Absolute path of a configured tool (config$tools$<name>), else the same name on PATH.
tool_path <- function(config, name, required = TRUE) {
  p <- config$tools[[name]] %||% ""
  if (!nzchar(p) || !file.exists(p)) p <- unname(Sys.which(name))
  # a plain error: callers with optional methods (prs) catch it; otherwise the pipeline reports it
  if (!nzchar(p) && required) {
    stop(sprintf("tool '%s' not found (config$tools$%s = %s); run scripts/install_postgwas_tools.sh",
                 name, name, config$tools[[name]] %||% "<unset>"), call. = FALSE)
  }
  p
}

# Run an external program, sending stdout+stderr to log_file. Fails the stage, quoting the
# log's last lines, if the exit status is non-zero or an expected output is missing.
run_tool <- function(bin, args, stage, log_file, expect = character(0), env = character(0), wd = NULL) {
  if (!is.null(wd)) {   # paths in args must then be absolute; forked workers change only their own wd
    log_file <- normalizePath(log_file, mustWork = FALSE)
    owd <- setwd(wd); on.exit(setwd(owd), add = TRUE)
  }
  status <- system2(bin, args, stdout = log_file, stderr = log_file, env = env)
  missing <- expect[!file.exists(expect)]
  if (!identical(as.integer(status), 0L) || length(missing)) {
    tail_txt <- if (file.exists(log_file)) paste(utils::tail(readLines(log_file, warn = FALSE), 30),
                                                 collapse = "\n") else "<no log>"
    log_fatal(stage, sprintf("%s exited with status %s%s\n--- %s (tail) ---\n%s",
                             basename(bin), status,
                             if (length(missing)) paste0("; missing: ", paste(basename(missing), collapse = ", ")) else "",
                             basename(log_file), tail_txt))
  }
  invisible(status)
}

# --- plink 1.9 wrappers (LD reference = 1000G EUR unless stated) ---------------------------

# All SNPs with r2 > r2_min to any SNP in `snps`, within kb (exact; unlike clumping with
# p2 = 1, a hit absorbed into another hit's clump still contributes its own partners).
# Returns data.table(SNP_A, SNP_B, R2), self-pairs excluded.
plink_ld_partners <- function(config, bfile, snps, r2_min, kb, out_prefix, stage, threads = 8L) {
  plink <- tool_path(config, "plink")
  snp_file <- paste0(out_prefix, ".snps")
  writeLines(unique(snps), snp_file)
  run_tool(plink, c("--bfile", bfile, "--r2", "--ld-snp-list", snp_file,
                    "--ld-window-kb", as.character(kb), "--ld-window", "999999",
                    "--ld-window-r2", as.character(r2_min), "--threads", as.character(threads),
                    "--out", out_prefix),
           stage, paste0(out_prefix, ".plink.log"), expect = paste0(out_prefix, ".ld"))
  ld <- fread(paste0(out_prefix, ".ld"), select = c("SNP_A", "SNP_B", "R2"))
  ld[SNP_A != SNP_B]
}

# PLINK 1.9 clumping. `assoc` has columns SNP and P. Returns one row per index SNP with
# its clumped members (SP2, parsed). Returns an empty table if nothing passes p1.
plink_clump <- function(config, bfile, assoc, p1, p2, r2, kb, out_prefix, stage,
                        chr = NULL, threads = 8L) {
  plink <- tool_path(config, "plink")
  in_file <- paste0(out_prefix, ".assoc")
  fwrite(assoc[, .(SNP, P)], in_file, sep = "\t")
  args <- c("--bfile", bfile, "--clump", in_file, "--clump-p1", format(p1), "--clump-p2", format(p2),
            "--clump-r2", format(r2), "--clump-kb", format(kb), "--threads", as.character(threads),
            "--out", out_prefix)
  if (!is.null(chr)) args <- c(args, "--chr", paste(chr, collapse = ","))
  run_tool(plink, args, stage, paste0(out_prefix, ".plink.log"))
  f <- paste0(out_prefix, ".clumped")
  if (!file.exists(f)) return(data.table(SNP = character(), CHR = integer(), BP = integer(),
                                          P = numeric(), members = list()))
  cl <- fread(f, select = c("CHR", "BP", "SNP", "P", "SP2"))
  cl[, members := lapply(strsplit(SP2, ",", fixed = TRUE), function(x) {
    x <- sub("\\(1\\)$", "", x); x[x != "NONE"]
  })]
  cl[, SP2 := NULL]
  cl[]
}

# --- intervals -----------------------------------------------------------------------------

# Merge overlapping [start, end] intervals per chromosome; returns chr, start, end and a
# list column of the row indices merged into each interval.
merge_intervals <- function(chr, start, end) {
  d <- data.table(chr = chr, start = start, end = end, i = seq_along(chr))
  setorder(d, chr, start)
  d[, grp := cumsum(c(TRUE, start[-1] > cummax(end)[-.N])), by = chr]
  d[, .(start = min(start), end = max(end), rows = list(i)), by = .(chr, grp)][, grp := NULL][]
}

# --- liftover ------------------------------------------------------------------------------

# Positions lifted with a UCSC chain (e.g. hg19ToHg38): returns the new position, NA when a
# position is unmapped, maps to several places or to another chromosome. chr: 1..22, 23 = X.
lift_positions <- function(chr, pos, chain_path, cache_dir = tempdir()) {
  if (grepl("\\.gz$", chain_path)) {
    plain <- file.path(cache_dir, sub("\\.gz$", "", basename(chain_path)))
    if (!file.exists(plain)) R.utils::gunzip(chain_path, plain, remove = FALSE)
    chain_path <- plain
  }
  chain <- rtracklayer::import.chain(chain_path)
  seq <- paste0("chr", ifelse(as.integer(chr) == 23L, "X", as.character(chr)))
  gr <- GenomicRanges::GRanges(seq, IRanges::IRanges(as.integer(pos), width = 1L), id = seq_along(pos))
  l <- unlist(rtracklayer::liftOver(gr, chain))
  ids <- S4Vectors::mcols(l)$id
  ok <- !(ids %in% ids[duplicated(ids)]) & as.character(GenomicRanges::seqnames(l)) == seq[ids]
  out <- rep(NA_integer_, length(pos))
  out[ids[ok]] <- as.integer(GenomicRanges::start(l)[ok])
  out
}

# --- writers -------------------------------------------------------------------------------

# Write a CHR/BP-sorted table as bgzipped TSV with a tabix index. The header line starts
# with '#' so tabix skips it.
write_bgz_tsv <- function(dt, path, seq_col = "CHR", pos_col = "BP") {
  stopifnot(grepl("\\.gz$", path))
  plain <- sub("\\.gz$", "", path)
  out <- copy(dt)
  setnames(out, seq_col, paste0("#", seq_col))
  fwrite(out, plain, sep = "\t", na = "NA", quote = FALSE)
  Rsamtools::bgzip(plain, dest = path, overwrite = TRUE)
  unlink(plain)
  cols <- names(out)
  Rsamtools::indexTabix(path, seq = match(paste0("#", seq_col), cols),
                        start = match(pos_col, cols), end = match(pos_col, cols), comment = "#")
  invisible(path)
}
