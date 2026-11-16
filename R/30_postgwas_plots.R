# Post-GWAS plots (GWASTutorial "Visualization"): Q_SNP Manhattans, QQ plots (factor full and
# qfilter tiers, Q_SNP) and regional plots for every fine-mapped region with a credible set:
# -log10 P coloured by r2 with the credible-set lead (UK Biobank in-sample LD), PIP coloured by
# credible set, and GENCODE v19 protein-coding genes. Factor Manhattans are left to R/10_plots.R.

pg_manhattan <- function(d, title, sig = 5e-8) {
  d <- d[is.finite(lp)]
  d <- d[lp > 3 | runif(.N) < 0.1]
  setorder(d, CHR, BP)
  off <- d[, .(m = max(BP)), by = CHR][, .(CHR, off = cumsum(c(0, head(as.numeric(m), -1))))]
  d <- merge(d, off, by = "CHR")
  d[, x := BP + off]
  ax <- d[, .(x = mean(range(x))), by = CHR]
  ggplot2::ggplot(d, ggplot2::aes(x, lp, colour = factor(CHR %% 2))) +
    ggplot2::geom_point(size = 0.4, show.legend = FALSE) +
    ggplot2::geom_hline(yintercept = -log10(sig), linetype = "dashed", colour = "red") +
    ggplot2::scale_colour_manual(values = c("grey35", "steelblue")) +
    ggplot2::scale_x_continuous(breaks = ax$x, labels = ifelse(ax$CHR == 23, "X", ax$CHR)) +
    ggplot2::labs(x = NULL, y = expression(-log[10](P)), title = title) + ggplot2::theme_minimal(base_size = 10)
}

pg_qq <- function(lp, title) {
  lp <- sort(lp[is.finite(lp)], decreasing = TRUE)
  n <- length(lp)
  e <- -log10(seq_len(n) / (n + 1))
  keep <- lp > 2 | runif(n) < 0.02
  lam <- median(qchisq(10^-lp, 1, lower.tail = FALSE)) / qchisq(0.5, 1)
  ggplot2::ggplot(data.frame(e = e[keep], o = lp[keep]), ggplot2::aes(e, o)) + ggplot2::geom_point(size = 0.5) +
    ggplot2::geom_abline(colour = "red") + ggplot2::theme_minimal(base_size = 10) +
    ggplot2::labs(x = "expected -log10 P", y = "observed -log10 P",
                  title = sprintf("%s (lambda_GC %.3f, from P on the 1-df scale)", title, lam))
}

pg_regional <- function(config, sex, cs_rows, genes) {
  x <- cs_rows[1]
  pip <- fread(file.path(stage_dir(config, sex, "finemap"), x$factor, x$ld_id, "pip.tsv.gz"))
  ld <- read_region_ld(config, sex, x$ld_id)
  leads <- cs_rows$lead_snp
  r2 <- if (leads[1] %in% ld$snps) ld$R[pip$SNP, leads[1]]^2 else NA_real_
  pip[, r2 := pmin(1, r2)]
  pip[, r2_bin := cut(r2, c(-Inf, 0.2, 0.4, 0.6, 0.8, Inf), labels = c("<0.2", "0.2-0.4", "0.4-0.6", "0.6-0.8", ">0.8"))]
  pip[, csf := fifelse(is.na(cs), "none", paste0("CS", cs))]
  rng <- range(pip$BP)
  g <- genes[CHR == x$lead_chr & end >= rng[1] & start <= rng[2]]
  g[, y := (seq_len(.N) - 1) %% 4]
  long <- rbind(pip[, .(BP, v = LOG10P, panel = "-log10 P", colour = as.character(r2_bin))],
                pip[, .(BP, v = pip, panel = "PIP", colour = csf)])
  long[, panel := factor(panel, c("-log10 P", "PIP", "genes"))]
  p <- ggplot2::ggplot(long, ggplot2::aes(BP / 1e6, v, colour = colour)) + ggplot2::geom_point(size = 0.8) +
    ggplot2::facet_grid(panel ~ ., scales = "free_y", drop = FALSE) +
    ggplot2::labs(x = sprintf("chr%d (Mb, GRCh37)", x$lead_chr), y = NULL, colour = NULL,
                  title = sprintf("%s %s: %d credible set(s), lead %s", x$factor, x$ld_id, nrow(cs_rows), leads[1])) +
    ggplot2::theme_bw(base_size = 9)
  if (nrow(g)) {
    gd <- g[, .(x0 = start / 1e6, x1 = end / 1e6, y, lab = symbol, panel = factor("genes", levels(long$panel)))]
    p <- p + ggplot2::geom_segment(data = gd, ggplot2::aes(x = x0, xend = x1, y = y, yend = y), inherit.aes = FALSE, linewidth = 1.5, colour = "darkgreen") +
      ggplot2::geom_text(data = gd, ggplot2::aes(x = (x0 + x1) / 2, y = y + 0.4, label = lab), inherit.aes = FALSE, size = 2.2)
  }
  p
}

run_postgwas_plots <- function(config, sex) {
  stage <- "postgwas_plots"
  log_info(stage, sprintf("=== Post-GWAS plots: %s ===", sex))
  out_dir <- stage_dir(config, sex, stage); dir.create(file.path(out_dir, "regional"), recursive = TRUE, showWarnings = FALSE)
  files <- character(0)
  save <- function(p, f, w, h) { ggplot2::ggsave(f, p, width = w, height = h, dpi = 150); files <<- c(files, f) }
  for (f in postgwas_factors(config)) {
    d <- load_postgwas(config, sex, f, "full", cols = c("CHR", "BP", "LOG10P", "Q_LOG10P", "pass_qfilter"))
    save(pg_manhattan(d[, .(CHR, BP, lp = Q_LOG10P)], sprintf("%s: Q_SNP heterogeneity", f)),
         file.path(out_dir, paste0(f, "_qsnp_manhattan.png")), 12, 4)
    save(pg_qq(d$LOG10P, sprintf("%s factor (full)", f)), file.path(out_dir, paste0(f, "_qq_full.png")), 5, 5)
    save(pg_qq(d[pass_qfilter == TRUE, LOG10P], sprintf("%s factor (qfilter)", f)), file.path(out_dir, paste0(f, "_qq_qfilter.png")), 5, 5)
    save(pg_qq(d$Q_LOG10P, sprintf("%s Q_SNP", f)), file.path(out_dir, paste0(f, "_qq_qsnp.png")), 5, 5)
    clear_postgwas_cache()
  }
  csf <- file.path(stage_dir(config, sex, "finemap"), "all_cs.tsv")
  if (file.exists(csf) && file.size(csf) > 0) {
    cs <- fread(csf)
    cs[, ld_id := sub(paste0("^", factor, "_"), "", fm_id), by = factor]
    genes <- read_gencode_genes(config$refs$gencode_gtf)
    for (id in unique(cs$fm_id)) {
      p <- tryCatch(pg_regional(config, sex, cs[fm_id == id][order(-lead_pip)], genes), error = function(e) NULL)
      if (!is.null(p)) save(p, file.path(out_dir, "regional", paste0(id, ".png")), 8, 6)
    }
  } else log_info(stage, "no fine-mapping results yet: regional plots skipped")
  write_stage_manifest(stage, sex, config, files)
  invisible(files)
}
