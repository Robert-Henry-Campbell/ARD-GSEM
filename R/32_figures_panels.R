# Panel functions for the figures stage. Each fig_<name>(config, sex) returns a ggplot carrying
# attr "source_data" (the aggregate table behind the panel), or NULL when its inputs are not
# there yet. Individual-level data never enter a panel or a source-data file.

.fac_lab <- c(ARD = "ARD (general)", Cardiometabolic = "Cardiometabolic", UpperGI = "Upper GI", Osteoarthritis = "Osteoarthritis (excluded)")

# --- Figure 1: factor model ---------------------------------------------------------------------

fig_model_loadings <- function(config, sex) {
  if (!file.exists(file.path(config$paths$output_dir, sex, "cfa", "model_bifactor_full.rds"))) return(NULL)
  l <- marker_loadings(config, sex)
  ord <- l[factor == "ARD"][order(std), indicator]
  l[, indicator := factor(indicator, ord)]
  l[, factor := factor(factor, names(.fac_lab), .fac_lab)]
  p <- ggplot2::ggplot(l, ggplot2::aes(factor, indicator, fill = std)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.2) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", std)), size = 1.6, colour = INK[["primary"]]) +
    ggplot2::scale_fill_gradient(low = "#f7fbff", high = OKABE_ITO[["blue"]], limits = c(0, 1), name = "Standardised\nloading") +
    ggplot2::scale_x_discrete(labels = function(x) sub(" \\(", "\n(", x), position = "top") +
    ggplot2::labs(x = NULL, y = NULL) + theme_nature(base_family = figure_font(config)) +
    ggplot2::theme(axis.line = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank())
  with_source_data(p, l[, .(factor, indicator, standardised_loading = std, unstandardised_loading = lambda)])
}

fig_model_rg_vs_loading <- function(config, sex) {
  f <- fig_inputs(config, sex, "factor_ldsc", "rg_vs_loadings.csv"); if (is.null(f)) return(NULL)
  d <- fread(f[[1]])
  d[, factor := factor(trait_a, names(.fac_lab), .fac_lab)]
  ann <- d[, .(cor_r = cor(rg, std_loading), n = .N), by = factor]
  ann[, `:=`(x = 0.05, y = 1 - 0.07 * seq_len(.N), label = sprintf("%s: r = %.2f (n = %d)", factor, cor_r, n))]
  p <- ggplot2::ggplot(d, ggplot2::aes(std_loading, rg, colour = factor)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = INK[["muted"]], linewidth = 0.25) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = rg - 1.96 * rg_se, ymax = rg + 1.96 * rg_se), width = 0, linewidth = 0.25) +
    ggplot2::geom_point(size = 1) +
    ggplot2::geom_text(data = ann, ggplot2::aes(x = x, y = y, label = label), hjust = 0, size = 2, show.legend = FALSE) +
    ggplot2::scale_colour_manual(values = setNames(FACTOR_COLS[names(.fac_lab)], .fac_lab), name = NULL) +
    ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    ggplot2::labs(x = "Standardised loading (model)", y = "Genetic correlation with factor GWAS") +
    theme_nature(base_family = figure_font(config)) + ggplot2::theme(legend.position = "bottom")
  with_source_data(p, d[, .(factor = trait_a, indicator = trait_b, rg, rg_se, std_loading)])
}

fig_model_h2 <- function(config, sex) {
  p1 <- file.path(config$paths$output_dir, sex, "ldsc", "h2_qc.csv")
  if (!file.exists(p1) || !file.exists(file.path(config$paths$output_dir, sex, "cfa", "model_bifactor_full.rds"))) return(NULL)
  h <- fread(p1)
  tcol <- intersect(c("trait", "phenotype", "code"), names(h))[1]
  h2col <- intersect(c("h2", "h2_liab", "h2_liability"), names(h))[1]
  secol <- intersect(c("h2_se", "se", "h2_liab_se"), names(h))[1]
  if (is.na(tcol) || is.na(h2col)) return(NULL)
  l <- marker_loadings(config, sex)[factor == "ARD", .(indicator, std)]
  d <- merge(data.table(indicator = h[[tcol]], h2 = h[[h2col]], se = if (is.na(secol)) NA_real_ else h[[secol]]), l, by = "indicator")
  d[, indicator := factor(indicator, d[order(std), indicator])]
  p <- ggplot2::ggplot(d, ggplot2::aes(h2, indicator, colour = std)) +
    ggplot2::geom_errorbar(ggplot2::aes(xmin = h2 - 1.96 * se, xmax = h2 + 1.96 * se), width = 0, orientation = "y", linewidth = 0.25, na.rm = TRUE) +
    ggplot2::geom_point(size = 1) +
    ggplot2::scale_colour_gradient(low = "#9ecae1", high = OKABE_ITO[["blue"]], name = "ARD\nloading") +
    ggplot2::labs(x = "SNP heritability (liability scale)", y = NULL) + theme_nature(base_family = figure_font(config), grid = "y")
  with_source_data(p, d)
}

# --- Figure 2: factor GWAS ----------------------------------------------------------------------

.fig_manhattan_data <- function(d, keep_above = 2, sample_frac = 0.02, seed = 1) {
  set.seed(seed)
  d <- d[is.finite(lp)]
  d <- d[lp > keep_above | stats::runif(.N) < sample_frac]
  setorder(d, CHR, BP)
  off <- d[, .(m = max(BP)), by = CHR][order(CHR)][, .(CHR, off = cumsum(c(0, head(as.numeric(m), -1))))]
  d <- merge(d, off, by = "CHR")
  d[, x := BP + off]
  list(d = d, axis = d[, .(x = mean(range(x))), by = CHR][order(CHR)])
}

fig_gwas_manhattan_panel <- function(config, sex, f, title = f, qsnp = FALSE) {
  if (!file.exists(qc_rds_path(config, sex, f))) return(NULL)
  cols <- c("CHR", "BP", if (qsnp) "Q_LOG10P" else "LOG10P", "is_mhc")
  d <- load_postgwas(config, sex, f, "full", cols = cols)
  setnames(d, if (qsnp) "Q_LOG10P" else "LOG10P", "lp")
  m <- .fig_manhattan_data(d)
  dd <- m$d; dd[, col := fifelse(is_mhc, "mhc", fifelse(CHR %% 2 == 1, "odd", "even"))]
  p <- ggplot2::ggplot(dd, ggplot2::aes(x, lp, colour = col)) +
    ggplot2::geom_point(size = 0.25, shape = 16, show.legend = FALSE) +
    ggplot2::geom_hline(yintercept = -log10(5e-8), linetype = "dashed", colour = OKABE_ITO[["vermilion"]], linewidth = 0.25) +
    ggplot2::scale_colour_manual(values = MANHATTAN_COLS) +
    ggplot2::scale_x_continuous(breaks = m$axis$x, labels = ifelse(m$axis$CHR == 23, "X", m$axis$CHR), expand = c(0.01, 0)) +
    ggplot2::scale_y_continuous(expand = c(0, 0), limits = c(0, max(dd$lp) * 1.05)) +
    ggplot2::labs(x = NULL, y = if (qsnp) expression(-log[10](italic(P)[Q[SNP]])) else expression(-log[10](italic(P))), title = title) +
    theme_nature(base_family = figure_font(config)) + ggplot2::theme(axis.line.x = ggplot2::element_blank())
  loci <- file.path(stage_dir(config, sex, "loci"), paste0(f, if (qsnp) "_qsnp_loci.tsv" else "_loci.tsv"))
  with_source_data(p, if (file.exists(loci)) fread(loci) else dd[lp > 7.3, .(CHR, BP, lp)])
}
fig_gwas_manhattan_ARD <- function(config, sex) fig_gwas_manhattan_panel(config, sex, "ARD", "ARD (general factor)")
fig_gwas_manhattan_Cardiometabolic <- function(config, sex) fig_gwas_manhattan_panel(config, sex, "Cardiometabolic", "Cardiometabolic factor")
fig_gwas_manhattan_UpperGI <- function(config, sex) fig_gwas_manhattan_panel(config, sex, "UpperGI", "Upper GI factor")
fig_ed_qsnp_manhattan_ARD <- function(config, sex) fig_gwas_manhattan_panel(config, sex, "ARD", "ARD: Q_SNP heterogeneity", qsnp = TRUE)
fig_ed_qsnp_manhattan_Cardiometabolic <- function(config, sex) fig_gwas_manhattan_panel(config, sex, "Cardiometabolic", "Cardiometabolic: Q_SNP heterogeneity", qsnp = TRUE)
fig_ed_qsnp_manhattan_UpperGI <- function(config, sex) fig_gwas_manhattan_panel(config, sex, "UpperGI", "Upper GI: Q_SNP heterogeneity", qsnp = TRUE)

fig_gwas_miami_ard_qsnp <- function(config, sex) {
  if (!file.exists(qc_rds_path(config, sex, "ARD"))) return(NULL)
  d <- load_postgwas(config, sex, "ARD", "full", cols = c("CHR", "BP", "LOG10P", "Q_LOG10P", "is_mhc"))
  up <- .fig_manhattan_data(d[, .(CHR, BP, lp = LOG10P, is_mhc)])
  dn <- .fig_manhattan_data(d[, .(CHR, BP, lp = Q_LOG10P, is_mhc)])
  dd <- rbind(up$d[, .(CHR, BP, x, is_mhc, y = lp, track = "ARD factor")], dn$d[, .(CHR, BP, x, is_mhc, y = -pmin(lp, 60), track = "Q_SNP")])
  dd[, col := fifelse(is_mhc, "mhc", fifelse(CHR %% 2 == 1, "odd", "even"))]
  p <- ggplot2::ggplot(dd, ggplot2::aes(x, y, colour = col)) +
    ggplot2::geom_point(size = 0.25, shape = 16, show.legend = FALSE) +
    ggplot2::geom_hline(yintercept = c(-log10(5e-8), log10(5e-8)), linetype = "dashed", colour = OKABE_ITO[["vermilion"]], linewidth = 0.25) +
    ggplot2::geom_hline(yintercept = 0, colour = INK[["primary"]], linewidth = 0.25) +
    ggplot2::scale_colour_manual(values = MANHATTAN_COLS) +
    ggplot2::scale_x_continuous(breaks = up$axis$x, labels = ifelse(up$axis$CHR == 23, "X", up$axis$CHR), expand = c(0.01, 0)) +
    ggplot2::scale_y_continuous(labels = function(v) abs(v)) +
    ggplot2::labs(x = NULL, y = expression(-log[10](italic(P))*": factor ↑, "*Q[SNP]*" ↓"),
                  caption = "Q_SNP track capped at 60") +
    theme_nature(base_family = figure_font(config)) + ggplot2::theme(axis.line.x = ggplot2::element_blank())
  with_source_data(p, dd[abs(y) > 7.3, .(track, CHR, BP, neg_log10_p = abs(y))])
}

fig_ed_qq_panel <- function(config, sex, f) {
  if (!file.exists(qc_rds_path(config, sex, f))) return(NULL)
  d <- load_postgwas(config, sex, f, "full", cols = c("LOG10P", "pass_qfilter"))
  qq <- function(lp, tier) {
    lp <- sort(lp[is.finite(lp)], decreasing = TRUE); n <- length(lp)
    e <- -log10(seq_len(n) / (n + 1)); keep <- lp > 2 | stats::runif(n) < 0.01
    lam <- stats::median(stats::qchisq(10^-lp, 1, lower.tail = FALSE)) / stats::qchisq(0.5, 1)
    data.table(expected = e[keep], observed = lp[keep], tier = sprintf("%s (lambda = %.3f)", tier, lam))
  }
  set.seed(1)
  dd <- rbind(qq(d$LOG10P, "full"), qq(d[pass_qfilter == TRUE, LOG10P], "Q-filtered"))
  p <- ggplot2::ggplot(dd, ggplot2::aes(expected, observed, colour = tier)) + ggplot2::geom_point(size = 0.3, shape = 16) +
    ggplot2::geom_abline(colour = INK[["muted"]], linewidth = 0.25) +
    ggplot2::scale_colour_manual(values = unname(OKABE_ITO[c("blue", "orange")]), name = NULL) +
    ggplot2::labs(x = expression(Expected~-log[10](italic(P))), y = expression(Observed~-log[10](italic(P))), title = .fac_lab[[f]]) +
    theme_nature(base_family = figure_font(config)) + ggplot2::theme(legend.position = c(0.3, 0.9))
  with_source_data(p, dd[observed > 5])
}
fig_ed_qq_ARD <- function(config, sex) fig_ed_qq_panel(config, sex, "ARD")
fig_ed_qq_Cardiometabolic <- function(config, sex) fig_ed_qq_panel(config, sex, "Cardiometabolic")
fig_ed_qq_UpperGI <- function(config, sex) fig_ed_qq_panel(config, sex, "UpperGI")

# --- Figure 3: loci and heterogeneity -----------------------------------------------------------

fig_loci_counts <- function(config, sex) {
  f <- fig_inputs(config, sex, "loci", "loci_summary.csv"); if (is.null(f)) return(NULL)
  s <- fread(f[[1]])
  d <- s[, .(factor, `Factor-level` = n_loci - n_heterogeneous, Heterogeneous = n_heterogeneous)]
  long <- melt(d, id.vars = "factor", variable.name = "class", value.name = "n")
  long[, factor := factor(factor, names(.fac_lab), .fac_lab)]
  lab <- s[, .(factor = factor(factor, names(.fac_lab), .fac_lab), n_loci, lab = sprintf("%d unique", n_unique))]
  p <- ggplot2::ggplot(long, ggplot2::aes(factor, n, fill = class)) + ggplot2::geom_col(width = 0.6) +
    ggplot2::geom_text(data = lab, ggplot2::aes(factor, n_loci, label = lab), inherit.aes = FALSE, vjust = -0.4, size = 2) +
    ggplot2::scale_fill_manual(values = unname(OKABE_ITO[c("blue", "orange")]), name = NULL) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.15))) +
    ggplot2::labs(x = NULL, y = "Genome-wide significant loci") + theme_nature(base_family = figure_font(config), grid = "y") +
    ggplot2::theme(legend.position = "top")
  with_source_data(p, s)
}

fig_loci_genome_map <- function(config, sex) {
  fs <- postgwas_factors(config)
  paths <- file.path(stage_dir(config, sex, "loci"), paste0(fs, "_loci.tsv"))
  if (!all(file.exists(paths))) return(NULL)
  d <- rbindlist(lapply(paths, fread), fill = TRUE)
  g <- file.path(stage_dir(config, sex, "annotate"), "lead_nearest_genes.tsv")
  if (file.exists(g)) d <- merge(d, fread(g)[, .(lead_snp = SNP, factor, gene = nearest_tss_gene)], by = c("lead_snp", "factor"), all.x = TRUE)
  chrlen <- c(249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663, 146364022, 141213431, 135534747, 135006516, 133851895,
              115169878, 107349540, 102531392, 90354753, 81195210, 78077248, 59128983, 63025520, 48129895, 51304566, 155270560)
  off <- cumsum(c(0, head(chrlen, -1)))
  d[, x := lead_bp + off[chr]]
  d[, factor := factor(factor, names(.fac_lab), .fac_lab)]
  d[, class := fifelse(heterogeneous, "Heterogeneous", fifelse(unique, "Factor-level, unique", "Factor-level"))]
  top <- d[order(-lead_log10p)][, head(.SD, if (requireNamespace("ggrepel", quietly = TRUE)) 8 else 4), by = factor]
  p <- ggplot2::ggplot(d, ggplot2::aes(x, factor, colour = factor, shape = class)) +
    ggplot2::geom_vline(xintercept = off[-1], colour = INK[["grid"]], linewidth = 0.2) +
    ggplot2::geom_point(size = 1.4, position = ggplot2::position_jitter(height = 0.15, width = 0, seed = 1)) +
    ggplot2::scale_shape_manual(values = c("Factor-level" = 16, "Factor-level, unique" = 17, "Heterogeneous" = 1), name = NULL) +
    ggplot2::scale_colour_manual(values = setNames(FACTOR_COLS[names(.fac_lab)], .fac_lab), guide = "none") +
    ggplot2::scale_x_continuous(breaks = off + chrlen / 2, labels = c(1:22, "X"), expand = c(0.005, 0)) +
    ggplot2::labs(x = "Chromosome", y = NULL) + theme_nature(base_family = figure_font(config)) +
    ggplot2::theme(legend.position = "top", axis.line.x = ggplot2::element_blank())
  if ("gene" %in% names(d)) {
    p <- p + (if (requireNamespace("ggrepel", quietly = TRUE)) ggrepel::geom_text_repel(data = top, ggplot2::aes(label = gene), size = 1.6, max.overlaps = 30, segment.size = 0.2, show.legend = FALSE)
              else ggplot2::geom_text(data = top, ggplot2::aes(label = gene), size = 1.6, vjust = -1, show.legend = FALSE))
  }
  with_source_data(p, d[, intersect(c("factor", "chr", "start", "end", "lead_snp", "lead_bp", "lead_log10p", "heterogeneous", "unique", "mhc", "gene"), names(d)), with = FALSE])
}

fig_loci_qsnp_drivers <- function(config, sex) {
  fs <- postgwas_factors(config)
  paths <- file.path(stage_dir(config, sex, "loci"), paste0(fs, "_qsnp_loci.tsv"))
  if (!all(file.exists(paths))) return(NULL)
  d <- rbindlist(lapply(paths, fread), fill = TRUE)[, .N, by = .(factor, driver)]
  d <- d[order(factor, -N)][, head(.SD, 10), by = factor]
  d[, factor := factor(factor, names(.fac_lab), .fac_lab)]
  d[, driver := factor(driver, rev(unique(driver)))]
  p <- ggplot2::ggplot(d, ggplot2::aes(N, driver, fill = factor)) + ggplot2::geom_col(width = 0.7) +
    ggplot2::facet_wrap(~factor, scales = "free_y", ncol = 1) +
    ggplot2::scale_fill_manual(values = setNames(FACTOR_COLS[names(.fac_lab)], .fac_lab), guide = "none") +
    ggplot2::labs(x = expression(Q[SNP]~"loci driven by the diagnosis"), y = NULL) + theme_nature(base_family = figure_font(config), grid = "xy")
  with_source_data(p, d)
}

# --- Figure 4: enrichment -----------------------------------------------------------------------

fig_sldsc_baseline <- function(config, sex) {
  f <- fig_inputs(config, sex, "sldsc", "baselineLD_results.tsv"); if (is.null(f)) return(NULL)
  r <- fread(f[[1]])
  r <- r[binary == TRUE & annotation != "base" & !grepl("flanking|MAFbin", annotation)]
  keep <- r[enrichment_q < 0.05, unique(annotation)]
  d <- r[annotation %in% keep]
  d[, annotation := factor(annotation, d[factor == "ARD"][order(Enrichment), annotation])]
  d[, factor := factor(factor, names(.fac_lab), .fac_lab)]
  d[, sig := enrichment_q < 0.05]
  p <- ggplot2::ggplot(d, ggplot2::aes(Enrichment, annotation, colour = factor, alpha = sig)) +
    ggplot2::geom_vline(xintercept = 1, colour = INK[["muted"]], linewidth = 0.25) +
    ggplot2::geom_errorbar(ggplot2::aes(xmin = Enrichment - 1.96 * Enrichment_std_error, xmax = Enrichment + 1.96 * Enrichment_std_error), width = 0,
                           orientation = "y", linewidth = 0.25, position = ggplot2::position_dodge(width = 0.6)) +
    ggplot2::geom_point(size = 1, position = ggplot2::position_dodge(width = 0.6)) +
    ggplot2::scale_colour_manual(values = setNames(FACTOR_COLS[names(.fac_lab)], .fac_lab), name = NULL) +
    ggplot2::scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.35), guide = "none") +
    ggplot2::scale_x_continuous(trans = "log2") +
    ggplot2::labs(x = "Heritability enrichment (log scale)", y = NULL) + theme_nature(base_family = figure_font(config), grid = "y") +
    ggplot2::theme(legend.position = "bottom")
  with_source_data(p, d[, .(factor, annotation, Prop._SNPs, Prop._h2, Enrichment, Enrichment_std_error, Enrichment_p, enrichment_q)])
}

fig_sldsc_tissues <- function(config, sex) {
  f <- fig_inputs(config, sex, "sldsc", "gtex_ldsc_seg_results.tsv"); if (is.null(f)) return(NULL)
  ct <- fread(f[[1]])
  ct[, system := fifelse(grepl("^Brain", Name), "Brain", fifelse(grepl("^Artery|^Heart", Name), "Cardiovascular",
               fifelse(grepl("^Adipose|^Liver|^Pancreas|^Muscle", Name), "Metabolic",
               fifelse(grepl("Esophagus|Stomach|Colon|Small_Intestine", Name), "Digestive", "Other"))))]
  ct[, factor := factor(factor, names(.fac_lab), .fac_lab)]
  ct[, Name := factor(Name, ct[factor == levels(factor)[1]][order(system, Coefficient_P_value), Name])]
  ct[, sig := q_bh < 0.05]
  p <- ggplot2::ggplot(ct, ggplot2::aes(Name, -log10(Coefficient_P_value), fill = system, alpha = sig)) +
    ggplot2::geom_col(width = 0.8) + ggplot2::facet_wrap(~factor, ncol = 1) +
    ggplot2::geom_hline(yintercept = -log10(0.05 / 54), linetype = "dashed", colour = INK[["muted"]], linewidth = 0.25) +
    ggplot2::scale_fill_manual(values = unname(OKABE_ITO[c("purple", "vermilion", "green", "orange", "sky")]), name = NULL) +
    ggplot2::scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.35), guide = "none") +
    ggplot2::labs(x = NULL, y = expression(-log[10](italic(P))~"tissue-specific enrichment")) +
    theme_nature(base_family = figure_font(config), grid = "y") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5, size = 4), legend.position = "top")
  with_source_data(p, ct[, .(factor, tissue = Name, system, Coefficient, Coefficient_std_error, Coefficient_P_value, q_bh)])
}

fig_ed_sldsc_all_annotations <- function(config, sex) {
  f <- fig_inputs(config, sex, "sldsc", "baselineLD_results.tsv"); if (is.null(f)) return(NULL)
  r <- fread(f[[1]])[annotation != "base"]
  r[, factor := factor(factor, names(.fac_lab), .fac_lab)]
  r[, z := `Coefficient_z-score`]
  r[, annotation := factor(annotation, r[factor == levels(factor)[1]][order(z), annotation])]
  p <- ggplot2::ggplot(r, ggplot2::aes(z, annotation, colour = factor, shape = binary)) +
    ggplot2::geom_vline(xintercept = c(-1.96, 1.96), linetype = "dashed", colour = INK[["muted"]], linewidth = 0.25) +
    ggplot2::geom_point(size = 0.9) +
    ggplot2::scale_colour_manual(values = setNames(FACTOR_COLS[names(.fac_lab)], .fac_lab), name = NULL) +
    ggplot2::scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 17), labels = c(`TRUE` = "binary", `FALSE` = "continuous"), name = NULL) +
    ggplot2::labs(x = expression(tau~"z-score (baselineLD v2.2)"), y = NULL) + theme_nature(base_family = figure_font(config), grid = "y") +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 4), legend.position = "bottom")
  with_source_data(p, r[, .(factor, annotation, binary, Coefficient, `Coefficient_std_error`, z, tau_q, Enrichment, enrichment_q)])
}

# --- Extended data: EFA ------------------------------------------------------------------------

fig_ed_scree <- function(config, sex) {
  p1 <- file.path(config$paths$output_dir, sex, "efa", "scree.csv"); if (!file.exists(p1)) return(NULL)
  s <- fread(p1)
  xcol <- intersect(c("factor", "component", "k", "n"), names(s))[1]; ycol <- intersect(c("eigenvalue", "eigen", "value"), names(s))[1]
  if (is.na(xcol) || is.na(ycol)) return(NULL)
  d <- data.table(k = s[[xcol]], eigenvalue = s[[ycol]])
  p <- ggplot2::ggplot(d, ggplot2::aes(k, eigenvalue)) + ggplot2::geom_line(linewidth = NATURE_LWD, colour = OKABE_ITO[["blue"]]) +
    ggplot2::geom_point(size = 0.9, colour = OKABE_ITO[["blue"]]) +
    ggplot2::geom_hline(yintercept = 1, linetype = "dashed", colour = INK[["muted"]], linewidth = 0.25) +
    ggplot2::labs(x = "Factor", y = "Eigenvalue (odd-chromosome genetic correlation matrix)") + theme_nature(base_family = figure_font(config), grid = "y")
  with_source_data(p, d)
}

fig_ed_rg_heatmap <- function(config, sex) {
  p1 <- file.path(config$paths$output_dir, sex, "ldsc", "ldsc_full.rds"); if (!file.exists(p1)) return(NULL)
  x <- readRDS(p1)
  S <- x$S_Stand %||% x$S
  if (is.null(S)) return(NULL)
  S <- as.matrix(S); nm <- colnames(S) %||% rownames(S)
  R <- stats::cov2cor(S); dimnames(R) <- list(nm, nm)
  ord <- nm[stats::hclust(stats::as.dist(1 - R))$order]
  d <- as.data.table(as.table(R)); setnames(d, c("a", "b", "rg"))
  d[, a := factor(a, ord)]; d[, b := factor(b, ord)]
  p <- ggplot2::ggplot(d, ggplot2::aes(a, b, fill = rg)) + ggplot2::geom_tile() +
    ggplot2::scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", limits = c(-1, 1), name = expression(r[g])) +
    ggplot2::labs(x = NULL, y = NULL) + theme_nature(base_family = figure_font(config)) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5, size = 4), axis.text.y = ggplot2::element_text(size = 4),
                   axis.line = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank())
  with_source_data(p, d)
}
