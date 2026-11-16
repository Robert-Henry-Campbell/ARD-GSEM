# LD panel stage: in-sample (UK Biobank) LD for fine-mapping and COJO, computed from the
# local UKB imputed genotypes (plink2 PGEN, read-only, another user's directory).
#  - samples: unrelated White British (22006 == 1, 22021 == 0, not a 22027 het/missing
#    outlier, no 22019 sex-chromosome aneuploidy), minus withdrawals when a list is supplied;
#    a seeded random 50k (dosage LD for SuSiE) and a nested 20k (hard calls for GCTA-COJO)
#  - regions: fine-mapping-eligible loci (autosomal, not MHC). A merged locus wider than
#    max_region_mb is split at the midpoints between its hits.
#  - per region: plink2 --r-unphased square ref-based (dosage r between REF-allele counts),
#    --freq, and a 20k-sample hard-call .bed. Nothing is written outside output/.

select_ld_samples <- function(config, out_dir, n_ld, n_cojo, seed) {
  stage <- "ld_panel"
  refs <- config$refs
  psam <- fread(file.path(refs$ukb_pgen_dir, "UKB_chr22.psam"), select = 1:2, col.names = c("FID", "IID"))
  qc <- fread(refs$ukb_sample_qc, colClasses = list(character = "eid"))
  col <- function(field) qc[[grep(paste0("^", field, "-0\\.0$"), names(qc))[1]]]
  steps <- list(in_pgen = nrow(psam))
  ok <- qc$eid %in% as.character(psam$IID); steps$in_qc_table <- sum(ok)
  ok <- ok & col("22006") %in% 1; steps$white_british <- sum(ok)
  ok <- ok & col("22021") %in% 0; steps$no_relatives <- sum(ok)
  ok <- ok & is.na(col("22027")); steps$not_het_missing_outlier <- sum(ok)
  ok <- ok & is.na(col("22019")); steps$no_sex_aneuploidy <- sum(ok)
  if (nzchar(refs$ukb_withdrawals %||% "") && file.exists(refs$ukb_withdrawals)) {
    wd <- fread(refs$ukb_withdrawals, header = FALSE, colClasses = "character")$V1
    ok <- ok & !(qc$eid %in% wd); steps$after_withdrawals <- sum(ok)
  } else {
    log_warn(stage, "refs$ukb_withdrawals not set: withdrawn participants are not removed (supply the current list for application 15181)")
  }
  pool <- sort(qc$eid[ok])
  if (length(pool) < n_ld) log_fatal(stage, sprintf("only %d eligible samples (< %d)", length(pool), n_ld))
  set.seed(seed)
  ld <- sort(sample(pool, n_ld))
  cojo <- sort(sample(ld, n_cojo))
  psam[, IID := as.character(IID)]
  fwrite(psam[IID %in% ld], file.path(out_dir, "keep_ld.txt"), sep = "\t", col.names = FALSE)
  fwrite(psam[IID %in% cojo], file.path(out_dir, "keep_cojo.txt"), sep = "\t", col.names = FALSE)
  Sys.chmod(file.path(out_dir, c("keep_ld.txt", "keep_cojo.txt")), mode = "0600")
  steps$ld_sample <- n_ld; steps$cojo_sample <- n_cojo; steps$seed <- seed
  jsonlite::write_json(steps, file.path(out_dir, "sample_selection.json"), auto_unbox = TRUE, pretty = TRUE)
  log_info(stage, sprintf("samples: %s", paste(sprintf("%s=%s", names(steps), unlist(steps)), collapse = ", ")))
  invisible(steps)
}

# Fine-mapping regions from loci/regions.tsv: eligible loci, split when wider than max_mb.
finemap_regions <- function(config, sex, max_mb) {
  r <- fread(file.path(stage_dir(config, sex, "loci"), "regions.tsv"))[finemap_eligible == TRUE]
  hits <- rbindlist(lapply(unique(r$factor), function(f)
    fread(file.path(stage_dir(config, sex, "loci"), paste0(f, "_independent_hits.tsv")))[, factor := f]))
  out <- list()
  for (i in seq_len(nrow(r))) {
    x <- r[i]
    if (x$end - x$start <= max_mb * 1e6) {
      out[[length(out) + 1]] <- x[, .(factor, locus_id = region_id, chr, start, end)]
      next
    }
    h <- sort(hits[factor == x$factor & CHR == x$chr & BP >= x$start & BP <= x$end, BP])
    cuts <- c(x$start, (h[-1] + h[-length(h)]) %/% 2L, x$end)
    # merge consecutive pieces until each is as wide as possible without exceeding max_mb
    pieces <- list(); s <- cuts[1]
    for (k in 2:length(cuts)) {
      if (cuts[k] - s > max_mb * 1e6 && cuts[k - 1] > s) { pieces[[length(pieces) + 1]] <- c(s, cuts[k - 1]); s <- cuts[k - 1] }
    }
    pieces[[length(pieces) + 1]] <- c(s, x$end)
    for (p in pieces) out[[length(out) + 1]] <- data.table(factor = x$factor, locus_id = x$region_id,
                                                            chr = x$chr, start = p[1], end = p[2])
  }
  out <- rbindlist(out)
  out[, ld_id := sprintf("chr%d_%d_%d", chr, start, end)]
  out[, fm_id := sprintf("%s_%s", factor, ld_id)]
  out[]
}

run_ld_panel <- function(config, sex) {
  stage <- "ld_panel"
  log_info(stage, sprintf("=== UKB LD panel: %s ===", sex))
  cfg <- config$ld_panel %||% list()
  out_dir <- stage_dir(config, sex, stage)
  for (d in c("ld", "cojo_bed", "snps", "logs")) dir.create(file.path(out_dir, d), recursive = TRUE, showWarnings = FALSE)
  Sys.chmod(out_dir, mode = "0700")
  plink2 <- tool_path(config, "plink2")
  pgen_dir <- config$refs$ukb_pgen_dir
  workers <- as.integer(cfg$workers %||% 4L)
  threads <- as.integer(cfg$plink_threads %||% 8L)

  select_ld_samples(config, out_dir, as.integer(cfg$n_ld %||% 50000L), as.integer(cfg$n_cojo %||% 20000L),
                    as.integer(cfg$seed %||% 20261108L))
  regs <- finemap_regions(config, sex, as.numeric(cfg$max_region_mb %||% 3))
  fwrite(regs, file.path(out_dir, "finemap_regions.tsv"), sep = "\t")
  ld_regs <- unique(regs[, .(ld_id, chr, start, end)])
  log_info(stage, sprintf("%d fine-mapping regions (%d distinct LD windows)", nrow(regs), nrow(ld_regs)))

  # Full-QC SNPs are the same for every factor (convergence is per SNP): take the first factor.
  f1 <- postgwas_factors(config)[1]
  snps <- load_postgwas(config, sex, f1, "full", cols = c("SNP", "CHR", "BP", "A1", "A2"))
  one <- function(i) {
    x <- ld_regs[i]
    pre <- file.path(out_dir, "ld", x$ld_id)
    if (file.exists(paste0(pre, ".done"))) return(x$ld_id)
    s <- snps[CHR == x$chr & BP >= x$start & BP <= x$end]
    pvar <- fread(file.path(pgen_dir, sprintf("UKB_chr%d.pvar", x$chr)), skip = "#CHROM",
                  col.names = c("CHROM", "POS", "ID", "REF", "ALT"))[POS >= x$start & POS <= x$end]
    m <- merge(s, pvar[, .(SNP = ID, BP = POS, REF, ALT)], by = c("SNP", "BP"))
    m <- m[(A1 == REF & A2 == ALT) | (A1 == ALT & A2 == REF)]
    snp_file <- file.path(out_dir, "snps", paste0(x$ld_id, ".snps"))
    writeLines(m$SNP, snp_file)
    fwrite(m[, .(SNP, BP, A1, A2, REF, ALT, a1_is_panel_ref = A1 == REF)],
           file.path(out_dir, "snps", paste0(x$ld_id, ".map.tsv")), sep = "\t")
    pf <- file.path(pgen_dir, sprintf("UKB_chr%d", x$chr))
    run_tool(plink2, c("--pfile", pf, "--keep", file.path(out_dir, "keep_ld.txt"), "--extract", snp_file,
                       "--r-unphased", "square", "ref-based", "bin4", "--freq", "--threads", threads,
                       "--memory", as.character(cfg$plink_memory_mb %||% 32000), "--out", pre),
             stage, file.path(out_dir, "logs", paste0(x$ld_id, ".ld.log")),
             expect = paste0(pre, c(".unphased.vcor1.bin", ".unphased.vcor1.bin.vars", ".afreq")))
    run_tool(plink2, c("--pfile", pf, "--keep", file.path(out_dir, "keep_cojo.txt"), "--extract", snp_file,
                       "--make-bed", "--threads", threads, "--out", file.path(out_dir, "cojo_bed", x$ld_id)),
             stage, file.path(out_dir, "logs", paste0(x$ld_id, ".bed.log")),
             expect = file.path(out_dir, "cojo_bed", paste0(x$ld_id, ".bed")))
    file.create(paste0(pre, ".done"))
    x$ld_id
  }
  res <- parallel::mclapply(seq_len(nrow(ld_regs)), function(i) tryCatch(one(i), error = function(e) e),
                            mc.cores = workers)
  bad <- vapply(res, inherits, logical(1), "error")
  if (any(bad)) log_fatal(stage, sprintf("%d LD window(s) failed, e.g. %s", sum(bad), conditionMessage(res[bad][[1]])))
  cov <- rbindlist(lapply(ld_regs$ld_id, function(id) {
    n_map <- nrow(fread(file.path(out_dir, "snps", paste0(id, ".map.tsv"))))
    x <- ld_regs[ld_id == id]
    data.table(ld_id = id, n_full_qc = nrow(snps[CHR == x$chr & BP >= x$start & BP <= x$end]), n_in_panel = n_map)
  }))
  cov[, coverage := n_in_panel / n_full_qc]
  fwrite(cov, file.path(out_dir, "panel_coverage.tsv"), sep = "\t")
  log_info(stage, sprintf("panel coverage of full-QC SNPs: median %.3f (min %.3f)", median(cov$coverage), min(cov$coverage)))
  clear_postgwas_cache()
  write_stage_manifest(stage, sex, config, c(file.path(out_dir, c("finemap_regions.tsv", "panel_coverage.tsv",
                                                                   "sample_selection.json"))))
  invisible(regs)
}

# Signed LD matrix for one window (SNP order of the .vars file), oriented to each SNP's A1.
read_region_ld <- function(config, sex, ld_id) {
  d <- file.path(stage_dir(config, sex, "ld_panel"))
  vars <- readLines(file.path(d, "ld", paste0(ld_id, ".unphased.vcor1.bin.vars")))
  n <- length(vars)
  con <- file(file.path(d, "ld", paste0(ld_id, ".unphased.vcor1.bin")), "rb")
  on.exit(close(con))
  R <- matrix(readBin(con, "numeric", n = n * n, size = 4), n, n)
  map <- fread(file.path(d, "snps", paste0(ld_id, ".map.tsv")))[match(vars, SNP)]
  sgn <- ifelse(map$a1_is_panel_ref, 1, -1)       # ref-based r -> r between A1 dosages
  R <- R * outer(sgn, sgn)
  dimnames(R) <- list(vars, vars)
  af <- fread(file.path(d, "ld", paste0(ld_id, ".afreq")))
  setnames(af, sub("^#", "", names(af)))
  alt_f <- af$ALT_FREQS[match(vars, af$ID)]
  list(R = R, snps = vars, eaf_panel_a1 = ifelse(map$a1_is_panel_ref, 1 - alt_f, alt_f))
}
