# Loci stage: independent factor hits, loci and fine-mapping regions, "unique" loci and
# Q_SNP (heterogeneity) loci.
#
# - Independent hits: plink clumping of the full-QC tier (1000G EUR), r2 < 0.1 within 500 kb,
#   1 Mb on chr6 and chr8 (Grotzinger et al. 2019). GWS SNPs absent from the panel, and chrX,
#   are added by greedy distance clumping (+-500 kb).
# - Loci: [hit +- locus_flank_kb] intervals merged per factor; anything touching the MHC is one
#   locus. A locus is "heterogeneous" when its lead fails the Q_SNP filter (a Q_SNP hit or an
#   LD partner of one): fine-mapped and looked up, but kept out of factor-level conclusions.
# - Unique: no lead of the locus is, or is in LD (r2 >= unique_r2) with, a genome-wide
#   significant SNP of any of the factor's indicator GWASs (the paper's definition).
# - Q_SNP loci: clumps of Q_SNP p < qsnp_p, with the indicator that drives each.

# Greedy distance clumping: most significant first, absorb everything within window_kb.
distance_clump_simple <- function(d, window_kb = 500) {
  d <- d[order(-LOG10P)]
  taken <- logical(nrow(d)); lead <- logical(nrow(d))
  for (i in seq_len(nrow(d))) {
    if (taken[i]) next
    lead[i] <- TRUE
    taken <- taken | (d$CHR == d$CHR[i] & abs(d$BP - d$BP[i]) <= window_kb * 1000)
  }
  d[lead]
}

# Independent hits for one factor's GWS SNPs (data.table with SNP CHR BP P LOG10P in_1000g).
independent_hits <- function(config, gws, bfile, cfg, prefix, stage) {
  r2 <- as.numeric(cfg$clump_r2 %||% 0.1)
  kb <- as.integer(cfg$clump_kb %||% 500L)
  kb68 <- as.integer(cfg$clump_kb_chr6_8 %||% 1000L)
  auto <- gws[CHR <= 22L & in_1000g == TRUE]
  # plink needs P > 0; recover underflowed P from LOG10P
  auto[, P := pmax(P, 10^-pmin(LOG10P, 300))]
  hits <- list()
  for (part in list(list(chr = setdiff(1:22, c(6L, 8L)), kb = kb, tag = "main"),
                    list(chr = c(6L, 8L), kb = kb68, tag = "chr6_8"))) {
    a <- auto[CHR %in% part$chr]
    if (!nrow(a)) next
    cl <- plink_clump(config, bfile, a, p1 = 5e-8, p2 = 1, r2 = r2, kb = part$kb,
                      out_prefix = paste0(prefix, "_", part$tag), stage = stage, chr = part$chr)
    if (nrow(cl)) hits[[part$tag]] <- cl[, .(SNP, method = "plink", n_members = lengths(members))]
  }
  hits <- rbindlist(hits)
  # GWS SNPs plink could not use (chrX, absent from the panel): distance-clump the ones not
  # already within a clump window of a plink lead.
  rest <- gws[!(CHR <= 22L & in_1000g == TRUE)]
  if (nrow(rest)) {
    leads <- gws[SNP %in% hits$SNP]
    near <- vapply(seq_len(nrow(rest)), function(i)
      any(leads$CHR == rest$CHR[i] & abs(leads$BP - rest$BP[i]) <= kb * 1000), logical(1))
    rest <- rest[!near]
    if (nrow(rest)) {
      dc <- distance_clump_simple(rest, kb)
      hits <- rbind(hits, dc[, .(SNP, method = "distance", n_members = NA_integer_)])
    }
  }
  merge(hits, gws, by = "SNP")[order(CHR, BP)]
}

# SNP pairs (lead -> partner) with r2 >= r2_min within kb, from 1000G EUR.
lead_ld_partners <- function(config, leads, bfile, r2_min, kb, prefix) {
  in_panel <- leads[in_1000g == TRUE & CHR <= 22L, SNP]
  if (!length(in_panel)) return(data.table(SNP_A = character(), SNP_B = character(), R2 = numeric()))
  plink_ld_partners(config, bfile, in_panel, r2_min, kb, prefix, "loci",
                    threads = as.integer(config$postgwas$workers %||% 8L))
}

# Genome-wide significant SNPs of each indicator, from the sumstats() output (Z = beta/se,
# oriented to the same A1 as the factor GWAS).
indicator_gws <- function(config, sex, indicators, p_gws) {
  path <- file.path(config$paths$output_dir, sex, "sumstats", paste0(sex, "_snp_sumstats.rds"))
  log_info("loci", sprintf("Reading indicator sumstats %s", basename(path)))
  ss <- as.data.table(readRDS(path))
  zthr <- qnorm(p_gws / 2, lower.tail = FALSE)
  zmat <- vapply(indicators, function(t) ss[[paste0("beta.", t)]] / ss[[paste0("se.", t)]],
                 numeric(nrow(ss)))
  out <- list(info = ss[, .(SNP, CHR = as.integer(CHR), BP = as.integer(BP))], z = zmat)
  rm(ss); gc(verbose = FALSE)
  hit_any <- rowSums(abs(zmat) > zthr, na.rm = TRUE) > 0
  out$gws_snps <- out$info$SNP[hit_any]
  out
}

run_loci <- function(config, sex) {
  stage <- "loci"
  log_info(stage, sprintf("=== Loci: %s ===", sex))
  pg <- config$postgwas %||% list(); cfg <- config$loci %||% list()
  out_dir <- stage_dir(config, sex, stage); work <- file.path(out_dir, "work")
  dir.create(work, recursive = TRUE, showWarnings = FALSE)
  bfile <- config$refs$ld_1000g_bfile
  p_gws <- as.numeric(pg$p_gws %||% 5e-8)
  qsnp_p <- as.numeric(pg$qsnp_p %||% 5e-8)
  flank <- as.integer(cfg$locus_flank_kb %||% 500L) * 1000L
  reg <- read_factor_registry(config, sex)
  factors <- postgwas_factors(config)

  all_ind <- unique(unlist(lapply(factors, function(f) reg[[f]]$indicators)))
  ind <- indicator_gws(config, sex, all_ind, p_gws)
  zthr <- qnorm(p_gws / 2, lower.tail = FALSE)

  regions <- list(); summ <- list()
  for (f in factors) {
    log_info(stage, sprintf("--- factor %s ---", f))
    dt <- load_postgwas(config, sex, f, "full",
                        cols = c("SNP", "CHR", "BP", "A1", "A2", "P", "LOG10P", "beta", "SE",
                                 "Q_SNP_pval", "Q_LOG10P", "qsnp_hit", "qsnp_ld", "pass_qfilter",
                                 "is_mhc", "in_1000g"))
    gws <- dt[P < p_gws | LOG10P > -log10(p_gws)]
    if (!nrow(gws)) {
      log_warn(stage, sprintf("%s: no genome-wide significant SNPs", f))
      next
    }
    hits <- independent_hits(config, gws, bfile, cfg, file.path(work, paste0(f, "_clump")), stage)
    log_info(stage, sprintf("%s: %s GWS SNPs -> %d independent hits (%d by distance)", f,
                            format(nrow(gws), big.mark = ","), nrow(hits), sum(hits$method == "distance")))

    # loci: hit +- flank, merged; the MHC is one locus
    iv <- merge_intervals(hits$CHR, pmax(1L, hits$BP - flank), hits$BP + flank)
    m <- config$postgwas$mhc %||% list(chr = 6, start = 25e6, end = 35e6)
    touch_mhc <- iv$chr == m$chr & iv$end >= m$start & iv$start <= m$end
    if (sum(touch_mhc) > 1L) {
      mhc_iv <- iv[touch_mhc][, .(chr = chr[1], start = min(start), end = max(end),
                                  rows = list(unlist(rows)))]
      iv <- rbind(iv[!touch_mhc], mhc_iv)[order(chr, start)]
    }
    loci <- rbindlist(lapply(seq_len(nrow(iv)), function(i) {
      h <- hits[iv$rows[[i]]]
      in_locus <- gws[CHR == iv$chr[i] & BP >= iv$start[i] & BP <= iv$end[i]]
      lead <- h[which.max(LOG10P)]
      data.table(factor = f, chr = iv$chr[i], start = iv$start[i], end = iv$end[i],
                 lead_snp = lead$SNP, lead_bp = lead$BP, lead_log10p = lead$LOG10P,
                 lead_beta = lead$beta, lead_se = lead$SE, lead_q_log10p = lead$Q_LOG10P,
                 n_hits = nrow(h), hits = paste(h$SNP, collapse = ","),
                 n_gws = nrow(in_locus), n_gws_qfail = sum(!in_locus$pass_qfilter),
                 lead_qsnp_hit = lead$qsnp_hit, lead_qsnp_ld = lead$qsnp_ld,
                 mhc = iv$chr[i] == m$chr & iv$end[i] >= m$start & iv$start[i] <= m$end,
                 chrX = iv$chr[i] == 23L)
    }))
    loci[, heterogeneous := lead_qsnp_hit | lead_qsnp_ld]
    loci[, finemap_eligible := !mhc & !chrX]

    # unique vs the factor's indicators
    f_ind <- reg[[f]]$indicators
    ind_hit <- rowSums(abs(ind$z[, f_ind, drop = FALSE]) > zthr, na.rm = TRUE) > 0
    ind_gws <- ind$info$SNP[ind_hit]
    ld <- lead_ld_partners(config, hits, bfile, as.numeric(cfg$unique_r2 %||% 0.1), 1000L,
                           file.path(work, paste0(f, "_unique")))
    hit_tagged <- unique(c(hits$SNP[hits$SNP %in% ind_gws], ld[SNP_B %in% ind_gws, SNP_A]))
    # leads without LD information: any indicator GWS SNP within the clump window
    no_ld <- hits[!(in_1000g == TRUE & CHR <= 22L)]
    for (i in seq_len(nrow(no_ld))) {
      near <- ind$info[ind_hit & CHR == no_ld$CHR[i] & abs(BP - no_ld$BP[i]) <= flank]
      if (nrow(near)) hit_tagged <- c(hit_tagged, no_ld$SNP[i])
    }
    loci[, unique := vapply(strsplit(hits, ","), function(s) !any(s %in% hit_tagged), logical(1))]
    hits[, tags_indicator_hit := SNP %in% hit_tagged]

    # Q_SNP loci and the indicator driving each
    qd <- dt[Q_SNP_pval < qsnp_p][, .(SNP, CHR, BP, P = pmax(Q_SNP_pval, 10^-pmin(Q_LOG10P, 300)),
                                      LOG10P = Q_LOG10P, in_1000g)]
    ql <- data.table()
    if (nrow(qd)) {
      qh <- independent_hits(config, qd, bfile, cfg, file.path(work, paste0(f, "_qclump")), stage)
      zi <- ind$z[match(qh$SNP, ind$info$SNP), f_ind, drop = FALSE]
      ql <- qh[, .(factor = f, SNP, CHR, BP, Q_LOG10P = LOG10P)]
      ql[, driver := f_ind[apply(abs(zi), 1, function(r) if (all(is.na(r))) NA_integer_ else which.max(r))]]
      ql[, driver_abs_z := apply(abs(zi), 1, max, na.rm = TRUE)]
      ql[, indicators_gws := apply(abs(zi) > zthr, 1, function(r) paste(f_ind[which(r)], collapse = ","))]
      fac <- dt[match(ql$SNP, dt$SNP)]
      ql[, `:=`(factor_log10p = fac$LOG10P, factor_gws = fac$LOG10P > -log10(p_gws))]
    }

    fwrite(hits, file.path(out_dir, paste0(f, "_independent_hits.tsv")), sep = "\t")
    fwrite(loci, file.path(out_dir, paste0(f, "_loci.tsv")), sep = "\t")
    fwrite(ql, file.path(out_dir, paste0(f, "_qsnp_loci.tsv")), sep = "\t")
    regions[[f]] <- loci
    summ[[f]] <- data.table(factor = f, n_gws = nrow(gws), n_hits = nrow(hits), n_loci = nrow(loci),
                            n_heterogeneous = sum(loci$heterogeneous), n_unique = sum(loci$unique),
                            n_finemap_eligible = sum(loci$finemap_eligible), n_mhc = sum(loci$mhc),
                            n_qsnp_loci = nrow(ql))
    log_info(stage, sprintf("%s: %d loci (%d heterogeneous, %d unique vs indicators, %d MHC); %d Q_SNP loci",
                            f, nrow(loci), sum(loci$heterogeneous), sum(loci$unique), sum(loci$mhc), nrow(ql)))
  }

  r <- rbindlist(regions)
  if (nrow(r)) {
    setorder(r, factor, chr, start)
    r[, region_id := sprintf("%s_chr%d_%d_%d", factor, chr, start, end)]
    setcolorder(r, "region_id")
    # loci shared across factors (overlapping intervals)
    ov <- r[, .(region_id, factor, chr, start, end)]
    setkey(ov, chr, start, end)
    pairs <- foverlaps(ov, ov, nomatch = NULL)[factor != i.factor & region_id < i.region_id]
    fwrite(pairs[, .(region_a = region_id, region_b = i.region_id, chr)],
           file.path(out_dir, "cross_factor_overlap.tsv"), sep = "\t")
  }
  fwrite(r, file.path(out_dir, "regions.tsv"), sep = "\t")
  fwrite(rbindlist(summ), file.path(out_dir, "loci_summary.csv"))
  write_stage_manifest(stage, sex, config, list.files(out_dir, full.names = TRUE))
  invisible(r)
}
