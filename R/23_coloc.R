# Colocalisation stage (GWASTutorial "Colocalization"): does a factor credible set share its
# causal variant with the QTL it was significant for in qtl_lookup?
#  - coloc.abf (one causal variant per trait) for every pair: factor beta/SE (sdY = 1, std.lv)
#    against the QTL's region statistics
#  - SuSiE-coloc (coloc.bf_bf) where both sides have SuSiE log Bayes factors: the factor fit from
#    the finemap stage against eQTL Catalogue's published lbf_variable files, or against a
#    susie_rss fit of the UKB-PPP pQTL using the same UKB in-sample LD (same population)
#  - eQTLGen (Z only): beta/SE from Z, allele frequency and N (SMR convention)
#  - priors p1 = p2 = 1e-4, p12 = 5e-6; colocalised = PP.H4 >= 0.8, suggestive 0.5-0.8,
#    distinct variants = PP.H3 >= 0.8

coloc_region_gwas <- function(config, sex, cs_row, fm_dir) {
  pip <- fread(file.path(fm_dir, cs_row$factor, cs_row$ld_id, "pip.tsv.gz"))
  g <- load_postgwas(config, sex, cs_row$factor, "full", cols = c("SNP", "beta", "SE", "EAF_A1"))
  pip <- merge(pip, g, by = "SNP")
  fit <- readRDS(file.path(fm_dir, cs_row$factor, cs_row$ld_id, "fit.rds"))
  lbf <- fit$lbf_variable
  colnames(lbf) <- fit$lbf_variable_snps
  # only this credible set's SuSiE component, so a neighbouring set's coloc is not reported here
  list(snps = pip, lbf = lbf[cs_row$susie_l, , drop = FALSE])
}

abf_pair <- function(gw, qtl, n_gwas, priors) {
  s <- intersect(gw$SNP, qtl$SNP)
  if (length(s) < 20) return(NULL)
  gw <- gw[match(s, SNP)]; qtl <- qtl[match(s, SNP)]
  d1 <- list(snp = s, beta = gw$beta, varbeta = gw$SE^2, type = "quant", sdY = 1, N = n_gwas,
             MAF = pmin(gw$EAF_A1, 1 - gw$EAF_A1))
  d2 <- list(snp = s, beta = qtl$beta, varbeta = qtl$se^2, type = "quant", sdY = 1, N = max(qtl$n, na.rm = TRUE),
             MAF = qtl$maf)
  r <- suppressWarnings(coloc::coloc.abf(d1, d2, p1 = priors$p1, p2 = priors$p2, p12 = priors$p12))
  top <- r$results$snp[which.max(r$results$SNP.PP.H4)]
  data.table(method = "abf", nsnps = length(s), t(r$summary[c("PP.H0.abf", "PP.H1.abf", "PP.H2.abf",
                                                              "PP.H3.abf", "PP.H4.abf")]), top_snp = top)
}

bfbf_pair <- function(lbf1, lbf2, priors) {
  s <- intersect(colnames(lbf1), colnames(lbf2))
  if (length(s) < 20 || !nrow(lbf2)) return(NULL)
  # coloc.bf_bf (5.2.3) re-indexes a matrix by SNP name only when it has columns beyond the shared
  # set; a matrix holding exactly the shared SNPs is used as is, so two such matrices in different
  # orders are paired by POSITION. Put both in one order (shared SNPs first, identically ordered,
  # then each matrix's own extra SNPs), so positions and names agree whichever branch coloc takes.
  # Full matrices are still passed: coloc checks how much posterior mass the shared SNPs carry.
  canon <- function(m) m[, c(s, setdiff(colnames(m), s)), drop = FALSE]
  r <- suppressWarnings(coloc::coloc.bf_bf(canon(lbf1), canon(lbf2), p1 = priors$p1, p2 = priors$p2, p12 = priors$p12))
  if (is.null(r$summary) || !nrow(r$summary)) return(NULL)
  sm <- as.data.table(r$summary)
  best <- sm[which.max(PP.H4.abf)]
  data.table(method = "susie_bf", nsnps = length(s), PP.H0.abf = best$PP.H0.abf, PP.H1.abf = best$PP.H1.abf,
             PP.H2.abf = best$PP.H2.abf, PP.H3.abf = best$PP.H3.abf, PP.H4.abf = best$PP.H4.abf,
             top_snp = best$hit1)
}

# All eQTL Catalogue records of one dataset in one GRCh38 region, cached per (dataset, region):
# every trait and credible set in that region reuses the read. Returns the cache file path.
eqtlcat_region_file <- function(ds, chr_num, start, end, tbi_dir, region_dir, tabix) {
  f <- file.path(region_dir, sprintf("%s_%s_%d_%d.tsv.gz", ds$dataset_id, chr_num, start, end))
  if (!file.exists(f)) {
    tbi <- eqtlcat_index(ds, tbi_dir)
    seqs <- tabix_seqnames(tabix, ds$ftp_path, tbi)
    chr <- eqtlcat_chr(chr_num, seqs)
    x <- if (chr %in% seqs) parse_eqtlcat_lines(tabix_read(tabix, ds$ftp_path, tbi, data.table(chr = chr, start = start, end = end))) else data.table()
    # a region of a genome-wide eQTL file on a chromosome it contains is never empty: treat that as
    # a failed read (not cached, retried) rather than as "no eQTL"
    if (!nrow(x) && chr %in% seqs) stop(sprintf("%s: empty read of %s:%d-%d", ds$dataset_id, chr, start, end), call. = FALSE)
    if (!nrow(x)) x <- as.data.table(setNames(replicate(length(EQTLCAT_COLS), character(0), simplify = FALSE), EQTLCAT_COLS))
    fwrite(x, paste0(f, ".tmp.gz"), sep = "\t"); file.rename(paste0(f, ".tmp.gz"), f)
  }
  f
}

# eQTL Catalogue region statistics for one molecular trait (GRCh38), aligned to A1.
eqtlcat_region <- function(ds, trait, gw, tbi_dir, region_dir, tabix, span = NULL) {
  gw38 <- gw[!is.na(BP38)]
  if (is.null(span)) span <- range(gw38$BP38)   # run_coloc passes the fine-mapping region's span (= prefetch key)
  f <- eqtlcat_region_file(ds, gw38$CHR[1], as.integer(span[1]), as.integer(span[2]), tbi_dir, region_dir, tabix)
  x <- fread(f, colClasses = list(character = c("chromosome", "rsid", "molecular_trait_id", "gene_id")))
  if (!nrow(x)) return(NULL)
  x <- x[molecular_trait_id == trait]
  x <- merge(gw38[, .(SNP, BP38, A1, A2)], x, by.x = "BP38", by.y = "position")
  x[, s := allele_sign(A1, A2, alt, ref)]
  x <- x[!is.na(s)][!duplicated(SNP)]
  x[, .(SNP, beta = s * beta, se, n = an / 2, maf, variant)]
}

# SuSiE lbf matrices (10 components x variants) of the needed traits of one dataset. The lbf files
# are 0.6-1.8 GB and not indexed, so each is streamed once (curl | zcat | awk) keeping only the
# needed traits; every stage of the pipe is checked (pipefail), and one .rds per trait is written.
# The trait list the cache was built from is stored; a different list rebuilds it.
eqtlcat_lbf_subset <- function(ds, traits, cache_dir) {
  traits <- sort(unique(traits))
  d <- file.path(cache_dir, ds$dataset_id); dir.create(d, recursive = TRUE, showWarnings = FALSE)
  gl <- file.path(d, "traits.txt")
  if (file.exists(gl) && identical(readLines(gl), traits)) return(d)
  tl <- tempfile(fileext = ".traits"); writeLines(traits, tl)
  tmp <- tempfile(fileext = ".tsv"); on.exit(unlink(c(tl, tmp)), add = TRUE)
  cmd <- sprintf("curl -sS --fail --retry 5 --retry-delay 20 --connect-timeout 60 %s | zcat | awk -F'\\t' 'NR==FNR{t[$1]=1; next} FNR==1 || ($1 in t)' %s - > %s",
                 shQuote(ds$ftp_lbf_path), shQuote(tl), shQuote(tmp))
  st <- system2("bash", c("-o", "pipefail", "-c", shQuote(cmd)))
  if (!identical(as.integer(st), 0L)) stop(sprintf("lbf stream for %s failed (exit %s)", ds$dataset_id, st), call. = FALSE)
  x <- fread(tmp, colClasses = list(character = c("molecular_trait_id", "variant")))
  unlink(list.files(d, "\\.rds$", full.names = TRUE))
  for (tr in unique(x$molecular_trait_id)) {
    y <- x[molecular_trait_id == tr]
    m <- t(as.matrix(y[, grep("^lbf_variable", names(y)), with = FALSE])); colnames(m) <- y$variant
    saveRDS(m, file.path(d, paste0(digest::digest(tr), ".rds")))
  }
  writeLines(traits, gl)
  d
}

# lbf matrix of one trait (NULL when the dataset has no SuSiE result for it).
eqtlcat_lbf <- function(ds, trait, cache_dir, traits = trait) {
  d <- eqtlcat_lbf_subset(ds, traits, cache_dir)
  f <- file.path(d, paste0(digest::digest(trait), ".rds"))
  if (file.exists(f)) readRDS(f) else NULL
}

# All eQTLGen rows for a set of genes in ONE pass over the 3.9 GB cis file (a per-pair zcat|awk scan
# takes ~2 min; 300+ pairs would take half a day). Cached in work/ with the gene list it was built from.
eqtlgen_subset <- function(config, genes, work) {
  genes <- sort(unique(genes))
  out <- file.path(work, "eqtlgen_subset.tsv"); gl <- file.path(work, "eqtlgen_subset.genes")
  if (!(file.exists(out) && file.exists(gl) && identical(readLines(gl), genes))) {
    writeLines(genes, paste0(gl, ".tmp"))
    cmd <- sprintf("zcat %s | awk -F'\\t' 'NR==FNR{g[$1]=1; next} FNR==1 || ($8 in g)' %s - > %s",
                   shQuote(config$refs$eqtlgen_cis), shQuote(paste0(gl, ".tmp")), shQuote(paste0(out, ".tmp")))
    st <- system(cmd)
    if (st != 0) stop("eQTLGen extraction failed (exit status ", st, ")")
    file.rename(paste0(out, ".tmp"), out); file.rename(paste0(gl, ".tmp"), gl)
  }
  x <- fread(out)
  setkey(x, Gene)
  x
}

# eQTLGen rows for one gene in the factor region: from the one-pass `subset` (eqtlgen_subset())
# when given, otherwise by scanning the full cis file.
eqtlgen_region <- function(config, gene, gw, af, subset = NULL) {
  rng <- range(gw$BP)
  if (!is.null(subset)) {
    e <- subset[.(gene), nomatch = NULL][SNPChr == gw$CHR[1] & SNPPos >= rng[1] & SNPPos <= rng[2]]
  } else {
    cmd <- sprintf("zcat %s | awk -F'\\t' 'NR==1 || ($8==\"%s\" && $3==%d && $4>=%d && $4<=%d)'",
                   shQuote(config$refs$eqtlgen_cis), gene, gw$CHR[1], rng[1], rng[2])
    e <- fread(cmd = cmd)
  }
  if (!nrow(e)) return(NULL)
  e <- merge(e, af, by = "SNP")
  e <- merge(gw[, .(SNP, A1, A2)], e, by = "SNP")
  e[, s := allele_sign(A1, A2, AssessedAllele, OtherAllele)]
  e <- e[!is.na(s)]
  p <- fifelse(e$AssessedAllele == e$AlleleB, e$AlleleB_all, 1 - e$AlleleB_all)
  den <- sqrt(2 * p * (1 - p) * (e$NrSamples + e$Zscore^2))
  e[, .(SNP, beta = s * Zscore / den, se = 1 / den, n = NrSamples, maf = pmin(p, 1 - p))]
}

ukbppp_region <- function(work_lookup, protein, gw) {
  f <- file.path(work_lookup, "ukbppp", "regions", paste0(protein, ".tsv.gz"))
  if (!file.exists(f)) return(NULL)
  u <- merge(gw[, .(SNP, A1, A2)], fread(f), by.x = "SNP", by.y = "rsid")
  u[, s := allele_sign(A1, A2, ALLELE1, ALLELE0)]
  u <- u[!is.na(s)][!duplicated(SNP)]
  u[, .(SNP, beta = s * BETA, se = SE, n = N, maf = pmin(A1FREQ, 1 - A1FREQ))]
}

run_coloc <- function(config, sex) {
  stage <- "coloc"
  log_info(stage, sprintf("=== Colocalisation: %s ===", sex))
  for (p in c("coloc", "susieR")) if (!requireNamespace(p, quietly = TRUE)) log_fatal(stage, paste(p, "not installed"))
  cfg <- config$coloc %||% list()
  priors <- list(p1 = as.numeric(cfg$p1 %||% 1e-4), p2 = as.numeric(cfg$p2 %||% 1e-4), p12 = as.numeric(cfg$p12 %||% 5e-6))
  out_dir <- stage_dir(config, sex, stage); work <- file.path(out_dir, "work")
  dir.create(work, recursive = TRUE, showWarnings = FALSE)
  fm_dir <- stage_dir(config, sex, "finemap")
  lk_dir <- stage_dir(config, sex, "qtl_lookup")
  reg <- read_factor_registry(config, sex)
  look <- fread(file.path(lk_dir, "lookup_all.tsv"))
  pairs <- if (identical(cfg$pairs %||% "bonferroni", "fdr")) look[sig_fdr == TRUE] else look[sig_bonf == TRUE]
  cs <- fread(file.path(fm_dir, "all_cs.tsv"))
  cs[, ld_id := sub(paste0("^", factor, "_"), "", fm_id), by = factor]
  ds_all <- eqtlcat_datasets(config, c("ge", "microarray", "aptamer"))
  af <- NULL
  log_info(stage, sprintf("%d lookup-significant pairs", nrow(pairs)))
  eg_genes <- unique(pairs[source == "eQTLGen", molecular_trait_id])
  eg <- if (length(eg_genes)) eqtlgen_subset(config, eg_genes, work) else NULL
  if (!is.null(eg)) log_info(stage, sprintf("eQTLGen: %s rows for %d genes extracted in one pass", format(nrow(eg), big.mark = ","), length(eg_genes)))
  pd <- file.path(work, "pairs"); dir.create(pd, showWarnings = FALSE)   # one checkpoint per pair
  # eQTL Catalogue keys carry a version: pairs computed before the 2026-11-12 tabix fix read
  # truncated regions and are never reused
  ck_path <- function(pr) file.path(pd, paste0(digest::digest(paste(pr$cs_id, pr$source, pr$dataset_id, pr$molecular_trait_id,
                                                                     priors$p1, priors$p2, priors$p12,
                                                                     if (pr$source == "eQTL_Catalogue") "tabix-cli-v2" else NULL, sep = "|")), ".rds"))
  # eQTL Catalogue inputs are fetched up front, in parallel, into per-(dataset, region) region files
  # and per-dataset lbf extracts; the pair loop then reads only local files.
  tabix <- tool_path(config, "tabix")
  tbi_dir <- file.path(lk_dir, "work", "eqtlcat_tbi"); region_dir <- file.path(work, "eqtlcat_regions")
  lbf_dir <- file.path(work, "eqtlcat_lbf"); dir.create(region_dir, showWarnings = FALSE)
  ec <- pairs[source == "eQTL_Catalogue"]
  ec <- ec[!file.exists(vapply(seq_len(nrow(ec)), function(i) ck_path(ec[i]), ""))]
  ec <- merge(ec, cs[, .(cs_id, ld_id)], by = "cs_id")
  lbf_traits <- split(pairs[source == "eQTL_Catalogue", molecular_trait_id], pairs[source == "eQTL_Catalogue", dataset_id])
  bp38 <- list()   # SNP -> GRCh38 position, per fine-mapping region
  region_bp38 <- function(f, ld_id) {
    if (is.null(bp38[[ld_id]])) {
      pip <- fread(file.path(fm_dir, f, ld_id, "pip.tsv.gz"), select = c("SNP", "CHR", "BP"))
      v <- setNames(lift_positions(pip$CHR, pip$BP, config$refs$chain_hg19_to_hg38, work), pip$SNP)
      attr(v, "chr") <- as.integer(pip$CHR[1])
      bp38[[ld_id]] <<- v
    }
    bp38[[ld_id]]
  }
  workers <- as.integer(cfg$workers %||% 6L)
  if (nrow(ec)) {
    rg <- unique(ec[, .(factor, ld_id)])
    for (k in seq_len(nrow(rg))) region_bp38(rg$factor[k], rg$ld_id[k])
    jobs <- unique(ec[, .(dataset_id, factor, ld_id)])
    jobs[, `:=`(chr = vapply(ld_id, function(l) attr(bp38[[l]], "chr"), 1L),
                start = vapply(ld_id, function(l) min(bp38[[l]], na.rm = TRUE), 1L),
                end = vapply(ld_id, function(l) max(bp38[[l]], na.rm = TRUE), 1L))]
    log_info(stage, sprintf("eQTL Catalogue: %d pairs to test; fetching %d (dataset, region) reads and %d lbf extracts (%d workers)",
                            nrow(ec), nrow(jobs), uniqueN(ec$dataset_id), workers))
    r1 <- parallel::mclapply(seq_len(nrow(jobs)), function(k) tryCatch({
      eqtlcat_region_file(ds_all[dataset_id == jobs$dataset_id[k]], jobs$chr[k], jobs$start[k], jobs$end[k], tbi_dir, region_dir, tabix); NULL
    }, error = function(e) conditionMessage(e)), mc.cores = workers, mc.preschedule = FALSE)
    e1 <- unlist(r1)
    if (length(e1)) log_warn(stage, sprintf("%d region reads failed (retried in the pair loop): %s", length(e1), e1[1]))
    dsl <- unique(ec$dataset_id)
    r2 <- parallel::mclapply(dsl, function(id) tryCatch({
      ds <- ds_all[dataset_id == id]
      if (!is.na(ds$ftp_lbf_path) && nzchar(ds$ftp_lbf_path)) eqtlcat_lbf_subset(ds, lbf_traits[[id]], lbf_dir); NULL
    }, error = function(e) conditionMessage(e)), mc.cores = min(workers, 4L), mc.preschedule = FALSE)
    e2 <- unlist(r2)
    if (length(e2)) log_warn(stage, sprintf("%d lbf extracts failed (retried in the pair loop): %s", length(e2), e2[1]))
  }
  rows <- list()
  for (i in seq_len(nrow(pairs))) {
    pr <- pairs[i]
    ck <- ck_path(pr)
    if (file.exists(ck)) { r0 <- readRDS(ck); if (nrow(r0)) rows[[length(rows) + 1]] <- r0; next }
    cr <- cs[cs_id == pr$cs_id]
    gwr <- coloc_region_gwas(config, sex, cr, fm_dir)
    gw <- gwr$snps
    n_gwas <- reg[[cr$factor]]$N_hat_std
    res <- list()
    ok <- tryCatch({   # one failing pair (download, LD, SuSiE) is logged and retried on rerun, not fatal
    if (pr$source == "eQTL_Catalogue") {
      ds <- ds_all[dataset_id == pr$dataset_id]
      rb <- region_bp38(cr$factor, cr$ld_id)
      gw[, BP38 := unname(rb[SNP])]
      qtl <- eqtlcat_region(ds, pr$molecular_trait_id, gw, tbi_dir, region_dir, tabix, span = range(rb, na.rm = TRUE))
      if (!is.null(qtl)) res$abf <- abf_pair(gw, qtl, n_gwas, priors)
      # an lbf failure now fails the pair (logged, not checkpointed, retried on rerun) instead of
      # silently dropping the SuSiE comparison
      lb <- if (is.na(ds$ftp_lbf_path) || !nzchar(ds$ftp_lbf_path)) NULL else
        eqtlcat_lbf(ds, pr$molecular_trait_id, lbf_dir, lbf_traits[[pr$dataset_id]])
      if (!is.null(lb) && !is.null(qtl)) {
        # rename variants matched to factor SNPs to rsIDs; unmatched ones keep their IDs (unshared)
        m <- match(colnames(lb), qtl$variant)
        colnames(lb) <- ifelse(is.na(m), colnames(lb), qtl$SNP[m])
        res$bf <- bfbf_pair(gwr$lbf, lb, priors)
      }
    } else if (pr$source == "eQTLGen") {
      if (is.null(af)) af <- fread(config$refs$eqtlgen_af, select = c("SNP", "AlleleB", "AlleleB_all"))
      qtl <- eqtlgen_region(config, pr$molecular_trait_id, gw, af, subset = eg)
      if (!is.null(qtl)) res$abf <- abf_pair(gw, qtl, n_gwas, priors)
    } else if (pr$source == "UKB-PPP") {
      qtl <- ukbppp_region(file.path(lk_dir, "work"), pr$molecular_trait_id, gw)
      if (!is.null(qtl)) {
        res$abf <- abf_pair(gw, qtl, n_gwas, priors)
        ld <- read_region_ld(config, sex, cr$ld_id)
        s <- intersect(qtl$SNP, ld$snps)
        if (length(s) >= 20) {
          q2 <- qtl[match(s, SNP)]
          sd2 <- file.path(work, "pqtl_susie"); dir.create(sd2, showWarnings = FALSE)
          fc <- file.path(sd2, paste0(digest::digest(list(pr$molecular_trait_id, cr$ld_id, s, q2$beta, q2$se)), ".rds"))
          lb2 <- if (file.exists(fc)) readRDS(fc) else {
            fit2 <- tryCatch(susieR::susie_rss(z = setNames(q2$beta / q2$se, s), R = ld$R[s, s], n = max(q2$n),
                                               L = 10, estimate_residual_variance = FALSE, max_iter = 1000),
                             error = function(e) NULL)
            l <- if (!is.null(fit2) && length(fit2$sets$cs)) {
              m <- fit2$lbf_variable[fit2$sets$cs_index, , drop = FALSE]; colnames(m) <- s; m
            } else matrix(numeric(0), 0, length(s), dimnames = list(NULL, s))
            if (!is.null(fit2)) { saveRDS(l, paste0(fc, ".tmp")); file.rename(paste0(fc, ".tmp"), fc) }
            l
          }
          if (nrow(lb2)) res$bf <- bfbf_pair(gwr$lbf, lb2, priors)
        }
      }
    }
    TRUE }, error = function(e) {
      log_warn(stage, sprintf("pair %d (%s, %s, %s) failed: %s", i, pr$source, pr$molecular_trait_id, pr$cs_id, conditionMessage(e))); FALSE })
    pr_rows <- rbindlist(lapply(res[!vapply(res, is.null, logical(1))], function(r)
      cbind(pr[, .(factor, cs_id, SNP, family, source, dataset_id, study, tissue,
                   molecular_trait_id, gene_symbol, p_lookup = p, heterogeneous, factor_level)], r)), fill = TRUE)
    if (ok) { saveRDS(pr_rows, paste0(ck, ".tmp")); file.rename(paste0(ck, ".tmp"), ck) }
    if (nrow(pr_rows)) rows[[length(rows) + 1]] <- pr_rows
    if (i %% 25 == 0 || i == nrow(pairs)) log_info(stage, sprintf("pair %d of %d", i, nrow(pairs)))
  }
  out <- rbindlist(rows, fill = TRUE)
  if (nrow(out)) {
    out[, call := fifelse(PP.H4.abf >= 0.8, "colocalised",
                          fifelse(PP.H4.abf >= 0.5, "suggestive", fifelse(PP.H3.abf >= 0.8, "distinct", "inconclusive")))]
  }
  fwrite(out, file.path(out_dir, "coloc_results.tsv"), sep = "\t")
  n_col <- if (nrow(out)) out[call == "colocalised", .(.N, uniqueN(cs_id))] else data.table(0L, 0L)
  log_info(stage, sprintf("%d coloc tests; %d colocalised (PP.H4 >= 0.8) across %d credible sets",
                          nrow(out), n_col[[1]], n_col[[2]]))
  clear_postgwas_cache()
  write_stage_manifest(stage, sex, config, file.path(out_dir, "coloc_results.tsv"))
  invisible(out)
}
