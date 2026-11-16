# QTL lookup stage: eQTL and pQTL associations of the fine-mapped PIP leads.
#  - query SNPs: each 95% credible set's lead plus its top-PIP members (max_cs_members); for a
#    dataset where the lead was not tested, the highest-PIP member that was is used instead
#  - eQTL family: eQTLGen phase 1 cis (blood, GRCh37 local file) and eQTL Catalogue ge and
#    microarray datasets (GRCh38, remote tabix); pQTL family: UKB-PPP (Olink, Synapse) and the
#    eQTL Catalogue aptamer dataset (Sun_2018)
#  - QTL effects are aligned to the factor's effect allele A1; "concordant" = same sign as the
#    factor z
#  - multiple testing: Bonferroni within each family over every (credible set, molecular trait,
#    dataset) test actually performed (primary), BH-FDR within family (secondary). A significant
#    lookup is only "supported" once coloc finds a shared variant (coloc stage).

# column order of the current (v7) *.all.tsv.gz files; scanTabix returns no header line
EQTLCAT_COLS <- c("molecular_trait_id", "chromosome", "position", "ref", "alt", "variant", "ma_samples",
                  "maf", "pvalue", "beta", "se", "type", "ac", "an", "r2", "molecular_trait_object_id",
                  "gene_id", "median_tpm", "rsid")

# Ensembl gene ID (unversioned) -> symbol for every GENCODE gene (any biotype).
gencode_gene_symbols <- function(gtf) {
  g <- fread(cmd = paste("zcat", shQuote(gtf), "| awk -F'\\t' '$3==\"gene\"'"), header = FALSE, sep = "\t",
             select = 9L, col.names = "attr")
  g[, .(gene_id = sub('.*gene_id "([^".]+)[^"]*".*', "\\1", attr), symbol = sub('.*gene_name "([^"]+)".*', "\\1", attr))]
}

qtl_query_snps <- function(config, sex, max_members) {
  cs <- fread(file.path(stage_dir(config, sex, "finemap"), "all_cs.tsv"))
  rbindlist(lapply(seq_len(nrow(cs)), function(i) {
    x <- cs[i]
    pip <- fread(file.path(stage_dir(config, sex, "finemap"), x$factor,
                           sub(paste0("^", x$factor, "_"), "", x$fm_id), "pip.tsv.gz"))[cs == x$cs]
    pip <- pip[order(-pip)][seq_len(min(.N, max_members))]
    pip[, .(cs_id = x$cs_id, factor = x$factor, SNP, CHR, BP, A1, A2, z, pip,
            is_lead = SNP == x$lead_snp, heterogeneous = x$heterogeneous, qflag = x$qflag,
            factor_level = x$factor_level)]
  }))
}

# Align a QTL row's effect to A1: +1, -1, or NA when the alleles do not match.
allele_sign <- function(a1, a2, eff, oth) {
  fifelse(eff == a1 & oth == a2, 1, fifelse(eff == a2 & oth == a1, -1, NA_real_))
}

lookup_eqtlgen <- function(config, q, work) {
  f <- config$refs$eqtlgen_cis
  if (is.null(f) || !file.exists(f)) { log_warn("qtl_lookup", "eQTLGen file missing; skipped"); return(NULL) }
  snp_file <- file.path(work, "eqtlgen_query.snps")
  writeLines(unique(q$SNP), snp_file)
  cmd <- sprintf("zcat %s | awk -F'\\t' -v f=%s 'BEGIN{while((getline l < f)>0) a[l]=1} NR==1 || ($2 in a)'",
                 shQuote(f), shQuote(snp_file))
  e <- fread(cmd = cmd, sep = "\t")
  if (!nrow(e)) return(NULL)
  m <- merge(q, e[, .(SNP, eff = AssessedAllele, oth = OtherAllele, Zscore, Pvalue, gene_id = Gene,
                      gene_symbol = GeneSymbol, n = NrSamples)], by = "SNP", allow.cartesian = TRUE)
  m[, s := allele_sign(A1, A2, eff, oth)]
  m <- m[!is.na(s)]
  m[, .(source = "eQTLGen", family = "eQTL", dataset_id = "eQTLGen_phase1", study = "eQTLGen",
        tissue = "blood", quant_method = "ge", cs_id, factor, SNP, is_lead, pip,
        molecular_trait_id = gene_id, gene_id, gene_symbol, beta = NA_real_, se = NA_real_,
        z_qtl = s * Zscore, p = Pvalue, n_qtl = n, z_factor = z, heterogeneous, factor_level)]
}

eqtlcat_datasets <- function(config, quant) {
  p <- fread(config$refs$eqtlcat_paths)[quant_method %in% quant]
  for (col in c("ftp_path", "ftp_cs_path", "ftp_lbf_path")) set(p, j = col, value = sub("^ftp://", "https://", p[[col]]))
  p[, family := fifelse(quant_method == "aptamer", "pQTL", "eQTL")]
  p[]
}

# Remote tabix reads use the htslib command-line tool. Rsamtools::scanTabix printed htslib read
# errors on remote files ("[E::hts_itr_next] Failed to seek ...: Illegal seek") without raising
# them, so failed regions came back empty and the 2026-11-10 lookup kept ~2% of the records. A
# read counts only if tabix exits 0 AND prints no "[E::" line; otherwise it is retried and then
# raised. The local .tbi is passed with htslib's "##idx##" syntax, so nothing is downloaded into
# the working directory. `regions`: data.table(chr, start, end), 1-based inclusive.
tabix_read <- function(tabix, src, index, regions, retries = 4L, timeout = 3600, wait = 20) {
  rf <- tempfile(fileext = ".regions.txt"); on.exit(unlink(rf), add = TRUE)
  fwrite(unique(regions)[order(chr, start)], rf, sep = "\t", col.names = FALSE)
  target <- if (nzchar(index %||% "")) paste0(src, "##idx##", index) else src
  last <- ""
  for (k in seq_len(retries)) {
    r <- tryCatch(processx::run(tabix, c(target, "-R", rf), error_on_status = FALSE, timeout = timeout),
                  error = function(e) list(status = -1L, stdout = "", stderr = conditionMessage(e)))
    if (identical(as.integer(r$status), 0L) && !grepl("[E::", r$stderr, fixed = TRUE)) {
      out <- strsplit(r$stdout, "\n", fixed = TRUE)[[1]]
      return(unique(out[nzchar(out)]))   # unique: overlapping regions print a record twice
    }
    last <- sprintf("status %s: %s", r$status, substr(gsub("\\s+", " ", r$stderr), 1, 300))
    if (k < retries) Sys.sleep(wait * k)
  }
  stop(sprintf("tabix read of %s failed after %d attempts (%s)", basename(src), retries, last), call. = FALSE)
}

# Sequence names in a tabix index (reads only the local index).
tabix_seqnames <- function(tabix, src, index) {
  r <- processx::run(tabix, c("-l", paste0(src, "##idx##", index)), error_on_status = FALSE, timeout = 300)
  if (!identical(as.integer(r$status), 0L) || grepl("[E::", r$stderr, fixed = TRUE))
    stop(sprintf("tabix -l %s failed: %s", basename(index), substr(r$stderr, 1, 300)), call. = FALSE)
  s <- strsplit(r$stdout, "\n", fixed = TRUE)[[1]]
  s[nzchar(s)]
}

# Local copy of a dataset's .tbi (downloaded once).
eqtlcat_index <- function(ds, cache_dir) {
  tbi <- file.path(cache_dir, paste0(ds$dataset_id, ".all.tsv.gz.tbi"))
  if (!file.exists(tbi) || file.size(tbi) < 100) {
    tmp <- paste0(tbi, ".tmp")
    utils::download.file(paste0(ds$ftp_path, ".tbi"), tmp, mode = "wb", quiet = TRUE)
    file.rename(tmp, tbi)
  }
  tbi
}

# eQTL Catalogue chromosome names for GWAS chromosomes (23 = X; "chr" prefix when the file uses it).
eqtlcat_chr <- function(chr_num, seqs) {
  chr <- ifelse(chr_num == 23L, "X", as.character(chr_num))
  if (length(seqs) && all(startsWith(seqs, "chr"))) chr <- paste0("chr", chr)
  chr
}

parse_eqtlcat_lines <- function(lines) {
  if (!length(lines)) return(data.table())
  fread(text = lines, header = FALSE, sep = "\t", col.names = EQTLCAT_COLS,
        colClasses = list(character = c("chromosome", "rsid", "molecular_trait_id", "gene_id")))
}

# Records of one eQTL Catalogue dataset at the GRCh38 positions in q38 (one remote read). Exact
# positions by default: htslib coalesces nearby point reads (GTEx liver, 1,570 positions: 13 s),
# whereas merging positions into windows (merge_gap > 0) pulls every variant x gene in between and
# was ~10x slower on these dense files. Records are cut back to the exact query positions.
query_eqtlcat_dataset <- function(ds, q38, cache_dir, tabix, retries = 4L, merge_gap = 0) {
  tbi <- eqtlcat_index(ds, cache_dir)
  seqs <- tabix_seqnames(tabix, ds$ftp_path, tbi)
  p <- unique(data.table(chr = eqtlcat_chr(q38$CHR, seqs), pos = as.integer(q38$BP38)))[chr %in% seqs][order(chr, pos)]
  if (!nrow(p)) { x <- data.table(); attr(x, "n_query") <- 0L; return(x) }
  p[, cl := cumsum(c(TRUE, diff(pos) > merge_gap | chr[-1] != chr[-.N]))]
  regions <- p[, .(start = min(pos), end = max(pos)), by = .(cl, chr)][, .(chr, start, end)]
  x <- parse_eqtlcat_lines(tabix_read(tabix, ds$ftp_path, tbi, regions, retries))
  if (nrow(x)) x <- x[paste(chromosome, position) %in% paste(p$chr, p$pos)]
  attr(x, "n_query") <- nrow(p)
  x
}

lookup_eqtlcat <- function(config, q, work, quant, workers, min_frac = 0.2) {
  ds <- eqtlcat_datasets(config, quant)
  cache <- file.path(work, "eqtlcat_tbi"); dir.create(cache, showWarnings = FALSE)
  q38 <- unique(q[, .(SNP, CHR, BP, A1, A2)])
  q38[, BP38 := lift_positions(CHR, BP, config$refs$chain_hg19_to_hg38, work)]
  log_info("qtl_lookup", sprintf("eQTL Catalogue: %d datasets x %d query SNPs (%d lifted to GRCh38)",
                                 nrow(ds), nrow(q38), sum(!is.na(q38$BP38))))
  q38 <- q38[!is.na(BP38)]
  tabix <- tool_path(config, "tabix")
  ckpt <- file.path(work, "eqtlcat_ckpt"); dir.create(ckpt, showWarnings = FALSE)
  one <- function(i) {
    d <- ds[i]
    f <- file.path(ckpt, paste0(d$dataset_id, ".rds"))
    if (file.exists(f)) return(readRDS(f))
    x <- query_eqtlcat_dataset(d, q38, cache, tabix)
    cov <- data.table(dataset_id = d$dataset_id, study = d$study_label, tissue = d$tissue_label,
                      quant_method = d$quant_method, n_query = attr(x, "n_query"), n_records = nrow(x),
                      n_query_found = if (nrow(x)) uniqueN(x$position) else 0L)
    if (nrow(x)) {
      x <- merge(q38, x, by.x = c("BP38"), by.y = c("position"), allow.cartesian = TRUE)
      x <- x[as.character(fifelse(CHR == 23L, "X", as.character(CHR))) == sub("^chr", "", chromosome)]
      x[, s := allele_sign(A1, A2, alt, ref)]
      x <- x[!is.na(s)][!duplicated(paste(SNP, molecular_trait_id))]
      x <- x[, .(dataset_id = d$dataset_id, study = d$study_label, tissue = d$tissue_label,
                 quant_method = d$quant_method, family = d$family, SNP, molecular_trait_id, gene_id,
                 beta = s * beta, se, p = pvalue, n_qtl = an / 2, maf, rsid_qtl = rsid)]
    }
    cov[, n_snps_matched := if (nrow(x)) uniqueN(x$SNP) else 0L]
    # backstop for reads that fail without an htslib error: a dataset covering < min_frac of the
    # query positions is not checkpointed (raised -> logged -> retried on the next run)
    if (cov$n_query > 0 && cov$n_query_found / cov$n_query < min_frac)
      stop(sprintf("%s: only %d of %d query positions returned", d$dataset_id, cov$n_query_found, cov$n_query))
    out <- list(x = x, coverage = cov)
    saveRDS(out, paste0(f, ".tmp")); file.rename(paste0(f, ".tmp"), f)
    out
  }
  is_bad <- function(r) inherits(r, "try-error") || !is.list(r) || is.null(r$coverage)
  res <- parallel::mclapply(seq_len(nrow(ds)), one, mc.cores = workers, mc.preschedule = FALSE)
  bad <- vapply(res, is_bad, logical(1))
  if (any(bad)) {   # second, gentler pass before giving up
    log_warn("qtl_lookup", sprintf("%d eQTL Catalogue dataset(s) failed; retrying with 2 connections", sum(bad)))
    res[bad] <- parallel::mclapply(which(bad), one, mc.cores = 2L, mc.preschedule = FALSE)
    bad <- vapply(res, is_bad, logical(1))
  }
  # incomplete eQTL Catalogue results would silently shrink the Bonferroni denominators again, so a
  # dataset that still fails stops the stage (finished datasets are checkpointed; a rerun resumes)
  if (any(bad)) log_fatal("qtl_lookup", sprintf("%d eQTL Catalogue dataset(s) still failing: %s", sum(bad),
                                                paste(utils::head(vapply(res[bad], function(r) substr(as.character(r), 1, 160), ""), 3), collapse = " | ")))
  cov <- rbindlist(lapply(res[!bad], `[[`, "coverage"))
  if (nrow(cov)) {
    cov[, frac_found := round(n_query_found / pmax(n_query, 1L), 3)]
    fwrite(cov[order(frac_found)], file.path(work, "eqtlcat_coverage.tsv"), sep = "\t")
    log_info("qtl_lookup", sprintf("eQTL Catalogue coverage: %d datasets read; query positions found per dataset median %.2f (min %.2f, %d below 0.5)",
                                   nrow(cov), stats::median(cov$frac_found), min(cov$frac_found), sum(cov$frac_found < 0.5)))
  }
  hits <- rbindlist(lapply(res[!bad], `[[`, "x"), fill = TRUE)
  if (!nrow(hits)) return(NULL)
  m <- merge(q[, .(cs_id, factor, SNP, is_lead, pip, z, heterogeneous, factor_level)], hits, by = "SNP",
             allow.cartesian = TRUE)
  m[, .(source = "eQTL_Catalogue", family, dataset_id, study, tissue, quant_method, cs_id, factor, SNP, is_lead,
        pip, molecular_trait_id, gene_id, gene_symbol = NA_character_, beta, se, z_qtl = beta / se, p,
        n_qtl, z_factor = z, heterogeneous, factor_level)]
}

# UKB-PPP (Synapse) via scripts/py/ukbppp_fetch.py: proteins whose gene TSS is within cis_kb of
# a credible-set lead, cis windows only. Needs SYNAPSE_AUTH_TOKEN.
lookup_ukbppp <- function(config, q, work, cis_kb) {
  # without a token, proceed only on cached regions (the fetch logs in only for uncached proteins)
  if (!nzchar(Sys.getenv("SYNAPSE_AUTH_TOKEN")) && !length(list.files(file.path(work, "ukbppp", "regions")))) {
    log_warn("qtl_lookup", "SYNAPSE_AUTH_TOKEN not set and no cached UKB-PPP regions: pQTL lookup skipped")
    return(NULL)
  }
  py <- tool_path(config, "py")
  genes <- read_gencode_genes(config$refs$gencode_gtf)
  leads <- unique(q[is_lead == TRUE, .(cs_id, CHR, BP)])
  near <- genes[leads[, .(cs_id, CHR, lo = BP - cis_kb * 1000, hi = BP + cis_kb * 1000, BP)],
                on = .(CHR, TSS >= lo, TSS <= hi), nomatch = NULL,
                .(symbol, chr = CHR, start = pmax(1, BP - cis_kb * 1000), end = BP + cis_kb * 1000)]
  proteins <- lapply(split(near, near$symbol), function(d)
    list(gene = d$symbol[1], chr = d$chr[1], windows = unname(as.matrix(unique(d[, .(start, end)])))))
  spec <- list(out_dir = file.path(work, "ukbppp"), cache_dir = config$refs$ukbppp_cache %||% file.path(work, "ukbppp_cache"),
               proteins = unname(proteins), snps = unique(q$SNP),
               keep_tars = isTRUE(config$qtl_lookup$ukbppp_keep_tars %||% FALSE))
  log_info("qtl_lookup", sprintf("UKB-PPP: %d candidate genes in cis of %d leads", length(proteins), nrow(leads)))
  spec_file <- file.path(work, "ukbppp_spec.json")
  jsonlite::write_json(spec, spec_file, auto_unbox = TRUE, digits = NA)
  run_tool(py, c(file.path(config$project$root, "scripts", "py", "ukbppp_fetch.py"), "--spec", spec_file),
           "qtl_lookup", file.path(work, "ukbppp_fetch.log"), expect = file.path(spec$out_dir, "lookup.tsv"))
  u <- fread(file.path(spec$out_dir, "lookup.tsv"))
  if (!nrow(u)) return(NULL)
  m <- merge(q, u, by.x = "SNP", by.y = "rsid", allow.cartesian = TRUE)
  m[, s := allele_sign(A1, A2, ALLELE1, ALLELE0)]
  m <- m[!is.na(s)]
  m[, .(source = "UKB-PPP", family = "pQTL", dataset_id = "UKB-PPP_EUR_discovery", study = "UKB-PPP",
        tissue = "plasma", quant_method = "olink", cs_id, factor, SNP, is_lead, pip,
        molecular_trait_id = protein, gene_id = NA_character_, gene_symbol = gene, beta = s * BETA, se = SE,
        z_qtl = s * BETA / SE, p = 10^-LOG10P, n_qtl = N, z_factor = z, heterogeneous, factor_level)]
}

# Per (credible set, molecular trait, dataset): the lead if tested, else the highest-PIP member.
pick_tested_snp <- function(x) {
  x <- x[order(cs_id, dataset_id, molecular_trait_id, -is_lead, -pip)]
  x[!duplicated(x[, .(cs_id, dataset_id, molecular_trait_id)])]
}

# Bonferroni within each QTL family over the tests actually performed (primary call) and BH-FDR
# per family (secondary); direction concordance between the QTL and factor z.
qtl_multiple_testing <- function(x) {
  x[, n_tests_family := .N, by = family]
  x[, `:=`(alpha_bonf = 0.05 / n_tests_family, sig_bonf = p <= 0.05 / n_tests_family,
           q_bh = p.adjust(p, "BH")), by = family]
  x[, sig_fdr := q_bh <= 0.05]
  x[, concordant := sign(z_qtl) == sign(z_factor)]
  x[]
}

run_qtl_lookup <- function(config, sex) {
  stage <- "qtl_lookup"
  log_info(stage, sprintf("=== eQTL / pQTL lookup of fine-mapped leads: %s ===", sex))
  cfg <- config$qtl_lookup %||% list()
  out_dir <- stage_dir(config, sex, stage); work <- file.path(out_dir, "work")
  dir.create(work, recursive = TRUE, showWarnings = FALSE)
  q <- qtl_query_snps(config, sex, as.integer(cfg$max_cs_members %||% 20L))
  log_info(stage, sprintf("%d credible sets, %d query SNPs", uniqueN(q$cs_id), uniqueN(q$SNP)))
  # the per-dataset / per-protein checkpoints are valid only for this query set (and cis window)
  qhash <- digest::digest(list(q[order(cs_id, SNP), .(cs_id, SNP, CHR, BP, A1, A2, is_lead)], cfg$cis_kb),
                          algo = "sha256")
  hf <- file.path(work, "query_hash.txt")
  if (!file.exists(hf) || !identical(readLines(hf), qhash)) {
    unlink(file.path(work, c("eqtlcat_ckpt", "ukbppp")), recursive = TRUE)
    writeLines(qhash, hf)
  }
  parts <- list(
    eqtlgen = lookup_eqtlgen(config, q, work),
    eqtlcat = lookup_eqtlcat(config, q, work, unlist(cfg$eqtlcat_quant %||% c("ge", "microarray", "aptamer")),
                             as.integer(cfg$workers %||% 8L), as.numeric(cfg$eqtlcat_min_found %||% 0.2)),
    ukbppp = lookup_ukbppp(config, q, work, as.integer(cfg$cis_kb %||% 1000L)))
  x <- rbindlist(parts, fill = TRUE)
  if (!nrow(x)) log_fatal(stage, "no QTL results from any source (see the warnings above)")
  x <- qtl_multiple_testing(pick_tested_snp(x))
  sym <- gencode_gene_symbols(config$refs$gencode_gtf)
  x[is.na(gene_symbol), gene_symbol := sym$symbol[match(sub("\\..*$", "", gene_id), sym$gene_id)]]
  setorder(x, family, p)
  fwrite(x, file.path(out_dir, "lookup_all.tsv"), sep = "\t")
  fwrite(x[sig_bonf == TRUE], file.path(out_dir, "lookup_sig.tsv"), sep = "\t")
  tests <- x[, .(n_tests = .N, n_tests_family = n_tests_family[1], alpha_bonf = alpha_bonf[1],
                 n_sig_bonf = sum(sig_bonf), n_sig_fdr = sum(sig_fdr),
                 n_cs_with_sig = uniqueN(cs_id[sig_bonf])), by = .(family, source)]
  fwrite(tests, file.path(out_dir, "tests.csv"))
  for (i in seq_len(nrow(tests))) {
    log_info(stage, sprintf("%s / %s: %s tests, Bonferroni alpha %.2e -> %d significant (%d credible sets); FDR<5%%: %d",
                            tests$family[i], tests$source[i], format(tests$n_tests[i], big.mark = ","),
                            tests$alpha_bonf[i], tests$n_sig_bonf[i], tests$n_cs_with_sig[i], tests$n_sig_fdr[i]))
  }
  write_stage_manifest(stage, sex, config, file.path(out_dir, c("lookup_all.tsv", "lookup_sig.tsv", "tests.csv")))
  invisible(x)
}
