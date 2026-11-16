# Panels for the PRS evaluation (figure 6, extended data 6-7) and the gene-level figure 5.
# All read aggregate tables written by prs_eval / finemap / coloc / gene_priority / magma.

.grp_lab <- c(SouthAsian = "South Asian", Black = "Black")
.fmt_thr <- function(x) ifelse(x >= 1e-3, sub("^0\\.", ".", format(x, scientific = FALSE, drop0trailing = TRUE)), sub("e-0", "e-", format(x, scientific = TRUE, digits = 1)))

.pe_results <- function(config, sex) {
  d <- stage_dir(config, sex, "prs_eval")
  fs <- list.files(d, pattern = "^results_.*\\.tsv$", full.names = TRUE)
  if (!length(fs)) { log_info("figures", "panel skipped: prs_eval results missing"); return(NULL) }
  r <- rbindlist(lapply(fs, fread), fill = TRUE)
  r[, group_label := .grp_lab[group]]
  r
}

# Primary weight set per (group, factor/comparator): the primary method; for C+T the threshold
# chosen in the other group (ct_selected_thresholds.tsv).
.pe_primary <- function(config, sex, r) {
  pm <- config$prs_eval$weights$primary_method %||% "ct"
  sel <- file.path(stage_dir(config, sex, "prs_eval"), "ct_selected_thresholds.tsv")
  keys <- unique(r[, .(group, name, kind, factor)])
  if (pm == "ct" && file.exists(sel)) {
    s <- fread(sel)
    keys <- merge(keys, s[, .(group, name, threshold)], by = c("group", "name"), all.x = TRUE)
    keys[, prs_id := sprintf("%s__ct__%s", name, threshold)]
  } else keys[, prs_id := sprintf("%s__%s", name, pm)]
  keys[prs_id %in% r$prs_id, .(group, name, kind, factor, prs_id)]
}

.prim_rows <- function(config, sex, r, outcome_pat, unit_val = "icd3") {
  p <- .pe_primary(config, sex, r)
  x <- merge(r[sample_set == "king" & covariate_set == "primary" & unit == unit_val & grepl(outcome_pat, outcome)], p, by = c("group", "name", "kind", "factor", "prs_id"))
  x
}

.forest <- function(d, y, x = "effect", xlab, colour = "group_label", vline = 1, config) {
  ggplot2::ggplot(d, ggplot2::aes(x = .data[[x]], y = .data[[y]], colour = .data[[colour]])) +
    ggplot2::geom_vline(xintercept = vline, colour = INK[["muted"]], linewidth = 0.25) +
    ggplot2::geom_errorbar(ggplot2::aes(xmin = ci_low, xmax = ci_high), width = 0, orientation = "y", linewidth = 0.3, position = ggplot2::position_dodge(width = 0.6)) +
    ggplot2::geom_point(size = 1.2, position = ggplot2::position_dodge(width = 0.6)) +
    ggplot2::scale_colour_manual(values = setNames(GROUP_COLS[names(.grp_lab)], .grp_lab), name = NULL) +
    ggplot2::scale_x_continuous(trans = "log2", breaks = function(l) { b <- pretty(l, 5); b[b > 0] }) +
    ggplot2::labs(x = xlab, y = NULL) + theme_nature(base_family = figure_font(config), grid = "y") +
    ggplot2::theme(legend.position = "bottom")
}

# Figure 6a: IRR per SD, factor scores on their own counts, both groups
fig_prs_irr <- function(config, sex) {
  r <- .pe_results(config, sex); if (is.null(r)) return(NULL)
  d <- .prim_rows(config, sex, r[kind == "factor"], "^count_")[outcome == paste0("count_", factor) & is.finite(effect)]
  if (!nrow(d)) return(NULL)
  d[, label := .fac_lab[factor]]
  d[, label := factor(label, rev(.fac_lab[c("ARD", "Cardiometabolic", "UpperGI")]))]
  p <- .forest(d, "label", xlab = "Incidence-rate ratio per SD of score (own indicator count)", config = config)
  with_source_data(p, d[, .(group, factor, prs_id, outcome, n, effect, ci_low, ci_high, p, lr_chisq, pseudo_r2)])
}

# Figure 6b: mean ARD count by score decile
fig_prs_decile <- function(config, sex) {
  d <- stage_dir(config, sex, "prs_eval")
  fs <- list.files(d, pattern = "^decile_.*\\.tsv$", full.names = TRUE)
  if (!length(fs)) return(NULL)
  r <- .pe_results(config, sex); if (is.null(r)) return(NULL)
  prim <- .pe_primary(config, sex, r[kind == "factor" & factor == "ARD"])
  dec <- rbindlist(lapply(fs, fread), fill = TRUE)
  dec <- merge(dec, prim[, .(group, prs_id)], by = c("group", "prs_id"))
  if (!nrow(dec)) return(NULL)
  dec[, group_label := .grp_lab[group]]
  p <- ggplot2::ggplot(dec, ggplot2::aes(decile, mean, colour = group_label)) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = ci_low, ymax = ci_high), width = 0, linewidth = 0.3) +
    ggplot2::geom_line(linewidth = NATURE_LWD) + ggplot2::geom_point(size = 1.2) +
    ggplot2::scale_x_continuous(breaks = 1:10) +
    ggplot2::scale_colour_manual(values = setNames(GROUP_COLS[names(.grp_lab)], .grp_lab), name = NULL) +
    ggplot2::labs(x = "ARD score decile", y = "Mean number of ARD diagnoses") + theme_nature(base_family = figure_font(config), grid = "y") +
    ggplot2::theme(legend.position = "bottom")
  with_source_data(p, dec[, .(group, prs_id, decile, n, mean, ci_low, ci_high)])
}

# Figure 6c: method comparison for the ARD score on the ARD count
fig_prs_methods <- function(config, sex) {
  f <- fig_inputs(config, sex, "prs_eval", "method_comparison.tsv"); if (is.null(f)) return(NULL)
  m <- fread(f[[1]])[kind == "factor" & factor == "ARD" & is.finite(irr)]
  if (!nrow(m)) return(NULL)
  m[, group_label := .grp_lab[group]]
  m[, thr := as.numeric(sub("^pT_", "", threshold))]
  m[, set := fifelse(method == "ct", paste0("C+T P<", .fmt_thr(thr)), c(prscs = "PRS-CS", ldpred2_auto = "LDpred2-auto", sbayesrc = "SBayesRC")[method])]
  m[, set := factor(set, rev(unique(set[order(method != "ct", -thr)])))]   # Bayesian methods first, then C+T from P<1 down
  setnames(m, c("irr"), c("effect"))
  p <- .forest(m, "set", xlab = "IRR per SD (ARD score -> ARD count)", config = config)
  with_source_data(p, m[, .(group, method, threshold, prs_id, n, effect, ci_low, ci_high, p, lr_chisq, pseudo_r2)])
}

# Figure 6d: comparator scores on the ARD count
fig_prs_comparators <- function(config, sex) {
  r <- .pe_results(config, sex); if (is.null(r)) return(NULL)
  d <- .prim_rows(config, sex, r, "^count_ARD$")[is.finite(effect)]
  if (!nrow(d) || !any(d$kind == "comparator")) return(NULL)
  ext <- read_external_manifest(config)
  d[, label := fifelse(kind == "factor", "ARD factor score (this study)", ext$label[match(name, ext$id)])]
  d[, label := factor(label, rev(unique(label[order(kind != "factor", -lr_chisq)])))]
  p <- .forest(d, "label", xlab = "IRR per SD (score -> ARD diagnosis count)", config = config)
  with_source_data(p, d[, .(group, kind, name, prs_id, n, effect, ci_low, ci_high, p, lr_chisq, pseudo_r2)])
}

# Figure 6e: FinnGen summary-level liability R2 per endpoint for the primary factor scores
fig_prs_finngen <- function(config, sex) {
  f <- fig_inputs(config, sex, "prs_eval", "finngen/results.tsv"); if (is.null(f)) return(NULL)
  r <- fread(f[[1]]); if (!nrow(r)) return(NULL)
  s <- r[kind == "factor"]
  best <- s[, .(m = mean(r2_liability)), by = .(factor, prs_id)][order(factor, -m)][, .SD[1], by = factor]
  d <- merge(s, best[, .(factor, prs_id)], by = c("factor", "prs_id"))
  mp <- config$prs_eval$finngen$endpoints %||% file.path(config$project$root, "manifest", "indicator_finngen_endpoints.tsv")
  if (!is_absolute_path(mp)) mp <- file.path(config$project$root, mp)
  if (file.exists(mp)) {   # "Obesity (E66)" instead of the FinnGen code
    em <- fread(mp)[, .(endpoint, ep_label = sprintf("%s (%s)", ifelse(nchar(endpoint_name) > 42, paste0(substr(endpoint_name, 1, 40), "\u2026"), endpoint_name), icd10))]
    em <- em[!duplicated(endpoint)]
    d <- merge(d, em, by = "endpoint", all.x = TRUE); d[is.na(ep_label), ep_label := endpoint]
  } else d[, ep_label := endpoint]
  d[, endpoint := factor(ep_label, d[factor == "ARD"][order(r2_liability), ep_label])]
  d[, flab := .fac_lab[factor]]
  p <- ggplot2::ggplot(d, ggplot2::aes(r2_liability, endpoint, colour = flab)) +
    ggplot2::geom_errorbar(ggplot2::aes(xmin = ci_low, xmax = ci_high), width = 0, orientation = "y", linewidth = 0.25, position = ggplot2::position_dodge(width = 0.5)) +
    ggplot2::geom_point(ggplot2::aes(shape = direction > 0), size = 1, position = ggplot2::position_dodge(width = 0.5)) +
    ggplot2::scale_colour_manual(values = setNames(FACTOR_COLS[names(.fac_lab)], .fac_lab), name = NULL) +
    ggplot2::scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), labels = c(`TRUE` = "risk-increasing", `FALSE` = "protective"), name = NULL) +
    ggplot2::labs(x = expression("Liability "*R^2*" in FinnGen (summary-level)"), y = NULL) +
    theme_nature(base_family = figure_font(config), grid = "y") + ggplot2::theme(legend.position = "bottom")
  with_source_data(p, d[, .(factor, prs_id, endpoint, icd10, n_cases, K, n_snps, r2_obs, r2_liability, ci_low, ci_high, direction)])
}

# Extended data 6: QC flow and within-group PCA
fig_ed_prs_qc_flow <- function(config, sex) {
  f <- fig_inputs(config, sex, "prs_eval", "sample_qc/sample_qc.json"); if (is.null(f)) return(NULL)
  j <- jsonlite::fromJSON(f[[1]])
  d <- rbindlist(lapply(names(j), function(g) data.table(group = g, step = names(j[[g]]), n = unlist(j[[g]]))))
  d <- d[!step %in% c("in_pgen", "n_final")]
  lab <- c(in_group = "Self-reported group", has_genetic_qc = "Genetic QC fields present", sex_match = "Reported sex = genetic sex",
           no_aneuploidy = "No sex-chromosome aneuploidy", no_het_missing_outlier = "Not a heterozygosity/missingness outlier",
           basic_keep = "After basic QC", unrelated_king = "Unrelated, KING 3rd degree (analysis set)",
           unrelated_gcta_list = "Unrelated, GCTA list (sensitivity)")
  d[, step := fifelse(step %in% names(lab), lab[step], step)]
  d[, step := factor(step, rev(unique(step)))]
  d[, group_label := .grp_lab[group]]
  p <- ggplot2::ggplot(d, ggplot2::aes(n, step, fill = group_label)) + ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.7), width = 0.6) +
    ggplot2::geom_text(ggplot2::aes(label = format(n, big.mark = ",")), position = ggplot2::position_dodge(width = 0.7), hjust = -0.1, size = 1.8) +
    ggplot2::scale_fill_manual(values = setNames(GROUP_COLS[names(.grp_lab)], .grp_lab), name = NULL) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0, 0.25))) +
    ggplot2::labs(x = "Participants", y = NULL) + theme_nature(base_family = figure_font(config), grid = "xy") + ggplot2::theme(legend.position = "bottom")
  with_source_data(p, d)
}

fig_ed_prs_pca <- function(config, sex) {
  d <- file.path(stage_dir(config, sex, "prs_eval"), "sample_qc", "work")
  fs <- list.files(d, pattern = "_pca\\.eigenvec$", full.names = TRUE)
  if (!length(fs)) return(NULL)
  ev <- rbindlist(lapply(fs, function(f) { x <- fread(f); data.table(group = sub("_pca\\.eigenvec$", "", basename(f)), PC1 = x$PC1, PC2 = x$PC2) }))
  ev[, group_label := .grp_lab[group]]
  p <- ggplot2::ggplot(ev, ggplot2::aes(PC1, PC2)) + ggplot2::geom_point(size = 0.2, alpha = 0.4, colour = OKABE_ITO[["blue"]]) +
    ggplot2::facet_wrap(~group_label, scales = "free") + theme_nature(base_family = figure_font(config))
  p  # individual-level coordinates: no source data by design
}

# Extended data 7: robustness
fig_ed_prs_gbd_vs_icd <- function(config, sex) {
  r <- .pe_results(config, sex); if (is.null(r)) return(NULL)
  d <- rbind(.prim_rows(config, sex, r[kind == "factor"], "^count_", "icd3"), .prim_rows(config, sex, r[kind == "factor"], "^count_", "gbd_cause"))
  d <- d[outcome %in% c(paste0("count_", factor), paste0("count_", factor, "_gbd")) & is.finite(effect)]
  if (!nrow(d)) return(NULL)
  d[, label := paste(.fac_lab[factor], fifelse(unit == "icd3", "ICD-10 code", "GBD cause"), sep = ": ")]
  d[, label := factor(label, rev(unique(label)))]
  p <- .forest(d, "label", xlab = "IRR per SD", config = config)
  with_source_data(p, d[, .(group, factor, unit, prs_id, outcome, n, effect, ci_low, ci_high, p)])
}

fig_ed_prs_mortality <- function(config, sex) {
  r <- .pe_results(config, sex); if (is.null(r)) return(NULL)
  d <- .prim_rows(config, sex, r, "^mortality$", "death")[is.finite(effect)]
  if (!nrow(d)) return(NULL)
  ext <- read_external_manifest(config)
  d[, label := fifelse(kind == "factor", paste(.fac_lab[factor], "score"), ext$label[match(name, ext$id)])]
  d[, label := factor(label, rev(unique(label)))]
  p <- .forest(d, "label", xlab = "HR per SD (all-cause mortality)", config = config)
  with_source_data(p, d[, .(group, kind, name, prs_id, n, n_events, effect, ci_low, ci_high, p, ph_p)])
}

fig_ed_prs_sensitivity <- function(config, sex) {
  r <- .pe_results(config, sex); if (is.null(r)) return(NULL)
  p0 <- .pe_primary(config, sex, r[kind == "factor" & factor == "ARD"])
  d <- merge(r[outcome == "count_ARD" & unit == "icd3" & is.finite(effect)], p0, by = c("group", "name", "kind", "factor", "prs_id"))
  if (!nrow(d)) return(NULL)
  d[, label := c("king primary" = "KING 3rd-degree unrelated, within-group PCs", "gcta_list primary" = "GCTA unrelated lists",
                 "king global_pcs" = "Global UKB PCs",
                 "king contact" = "Hospital-contact adjusted\n(log distinct diagnosis dates)")[paste(sample_set, covariate_set)]]
  d <- d[!is.na(label)]
  d[, label := factor(label, rev(unique(label)))]
  p <- .forest(d, "label", xlab = "IRR per SD (ARD score -> ARD count)", config = config)
  with_source_data(p, d[, .(group, sample_set, covariate_set, prs_id, n, effect, ci_low, ci_high, p)])
}

fig_ed_prs_thresholds <- function(config, sex) {
  r <- .pe_results(config, sex); if (is.null(r)) return(NULL)
  d <- r[kind == "factor" & method == "ct" & sample_set == "king" & covariate_set == "primary" & unit == "icd3" & outcome == paste0("count_", factor) & is.finite(lr_chisq)]
  if (!nrow(d)) return(NULL)
  d[, thr := as.numeric(sub("^pT_", "", threshold))]
  d[, flab := .fac_lab[factor]]
  p <- ggplot2::ggplot(d, ggplot2::aes(thr, lr_chisq, colour = .grp_lab[group])) + ggplot2::geom_line(linewidth = NATURE_LWD) + ggplot2::geom_point(size = 1) +
    ggplot2::facet_wrap(~flab, scales = "free_y", ncol = 3) +
    ggplot2::scale_x_log10(breaks = c(1e-6, 1e-4, 1e-2, 1), labels = c("1e-6", "1e-4", ".01", "1")) +
    ggplot2::scale_colour_manual(values = setNames(GROUP_COLS[names(.grp_lab)], .grp_lab), name = NULL) +
    ggplot2::labs(x = "C+T P threshold", y = expression("LR "*chi^2*" (score -> own count)")) +
    theme_nature(base_family = figure_font(config), grid = "y") +
    ggplot2::theme(legend.position = "bottom", panel.spacing.x = grid::unit(2, "mm"), strip.text = ggplot2::element_text(size = 6))
  with_source_data(p, d[, .(group, factor, threshold, prs_id, n, effect, ci_low, ci_high, lr_chisq, pseudo_r2)])
}

# Negative control (chapter-S injury count) next to the ARD count, each with and without the
# hospital-contact adjustment: HES-derived counts all rise with hospital contact.
fig_ed_prs_negative_control <- function(config, sex) {
  r <- .pe_results(config, sex); if (is.null(r)) return(NULL)
  p0 <- .pe_primary(config, sex, r[kind == "factor" & factor == "ARD"])
  d <- merge(r[sample_set == "king" & covariate_set %in% c("primary", "contact") & unit == "icd3" &
                 outcome %in% c("count_ARD", "injury_count") & is.finite(effect)],
             p0, by = c("group", "name", "kind", "factor", "prs_id"))
  if (!nrow(d)) return(NULL)
  lv <- c("ARD diagnoses", "ARD diagnoses,\nhospital-contact adjusted",
          "Injury diagnoses, chapter S\n(negative control)", "Injury diagnoses, chapter S\n(negative control),\nhospital-contact adjusted")
  d[, label := fifelse(outcome == "injury_count", fifelse(covariate_set == "contact", lv[4], lv[3]), fifelse(covariate_set == "contact", lv[2], lv[1]))]
  d[, label := factor(label, rev(lv))]
  p <- .forest(d, "label", xlab = "IRR per SD of the ARD score", config = config)
  with_source_data(p, d[, .(group, outcome, covariate_set, prs_id, n, effect, ci_low, ci_high, p)])
}

# --- Figure 5: genes (after finemap / coloc / gene_priority / magma) --------------------------

fig_finemap_cs <- function(config, sex) {
  f <- fig_inputs(config, sex, "finemap", "all_cs.tsv"); if (is.null(f)) return(NULL)
  cs <- fread(f[[1]]); if (!nrow(cs)) return(NULL)
  cs[, flab := .fac_lab[factor]]
  cs[, class := fifelse(factor_level, "factor-level", "heterogeneous")]
  p <- ggplot2::ggplot(cs, ggplot2::aes(size, lead_pip, colour = flab, shape = class)) + ggplot2::geom_point(size = 1, alpha = 0.8) +
    ggplot2::scale_x_log10() + ggplot2::scale_shape_manual(values = c("factor-level" = 16, "heterogeneous" = 1), name = NULL) +
    ggplot2::scale_colour_manual(values = setNames(FACTOR_COLS[names(.fac_lab)], .fac_lab), name = NULL) +
    ggplot2::labs(x = "Credible-set size (variants)", y = "Lead variant PIP") + theme_nature(base_family = figure_font(config), grid = "xy") +
    ggplot2::theme(legend.position = "bottom")
  with_source_data(p, cs[, .(factor, cs_id, lead_snp, lead_chr, lead_bp, lead_pip, size, purity_min, heterogeneous, factor_level)])
}

fig_coloc_heatmap <- function(config, sex) {
  f <- fig_inputs(config, sex, "coloc", "coloc_results.tsv"); if (is.null(f)) return(NULL)
  cl <- fread(f[[1]]); if (!nrow(cl)) return(NULL)
  d <- cl[, .(h4 = max(PP.H4.abf)), by = .(factor, gene_symbol, dataset = paste(study, tissue))][h4 >= 0.5]
  if (!nrow(d)) return(NULL)
  p <- ggplot2::ggplot(d, ggplot2::aes(dataset, gene_symbol, fill = h4)) + ggplot2::geom_tile(colour = "white", linewidth = 0.2) +
    ggplot2::facet_grid(factor ~ ., scales = "free_y", space = "free_y") +
    ggplot2::scale_fill_gradient(low = "#fee8c8", high = OKABE_ITO[["vermilion"]], limits = c(0.5, 1), name = "PP.H4") +
    ggplot2::labs(x = NULL, y = NULL) + theme_nature(base_family = figure_font(config)) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5, size = 4), axis.line = ggplot2::element_blank())
  with_source_data(p, d)
}

fig_gene_priority <- function(config, sex) {
  f <- fig_inputs(config, sex, "gene_priority", "cs_top_genes.tsv"); if (is.null(f)) return(NULL)
  g <- fread(f[[1]]); if (!nrow(g)) return(NULL)
  ev <- c("nearest_tss", "cs_overlap", "coding_pip", "magma_sig", "twas_sig", "eqtl_coloc_h4", "pqtl_coloc_h4")
  long <- melt(g[, c("factor", "cs_id", "symbol", intersect(ev, names(g))), with = FALSE], id.vars = c("factor", "cs_id", "symbol"), variable.name = "evidence")
  long[, has := fifelse(evidence %in% c("coding_pip"), !is.na(value) & value >= 0.1, fifelse(grepl("coloc", evidence), !is.na(value) & value >= 0.8, !is.na(value) & value == 1))]
  long[, label := paste0(symbol, " (", cs_id, ")")]
  p <- ggplot2::ggplot(long, ggplot2::aes(evidence, label, fill = has)) + ggplot2::geom_tile(colour = "white", linewidth = 0.2) +
    ggplot2::facet_grid(factor ~ ., scales = "free_y", space = "free_y") +
    ggplot2::scale_fill_manual(values = c(`TRUE` = OKABE_ITO[["blue"]], `FALSE` = "#f0f0f0"), guide = "none") +
    ggplot2::labs(x = NULL, y = NULL) + theme_nature(base_family = figure_font(config)) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1), axis.text.y = ggplot2::element_text(size = 4), axis.line = ggplot2::element_blank())
  with_source_data(p, long[, .(factor, cs_id, symbol, evidence, value, has)])
}

fig_magma_sets <- function(config, sex) {
  d <- stage_dir(config, sex, "magma")
  fs <- file.path(d, paste0(postgwas_factors(config), "_msigdb.gsa.out"))
  if (!all(file.exists(fs))) return(NULL)
  r <- rbindlist(lapply(seq_along(fs), function(i) { x <- magma_read_out(fs[i]); x[, factor := postgwas_factors(config)[i]]; x }), fill = TRUE)
  r[, q := p.adjust(P, "BH"), by = factor]
  top <- r[order(factor, P)][, head(.SD, 8), by = factor]
  top[, set := factor(FULL_NAME, rev(unique(FULL_NAME)))]
  top[, flab := .fac_lab[factor]]
  p <- ggplot2::ggplot(top, ggplot2::aes(-log10(P), set, fill = flab, alpha = q < 0.05)) + ggplot2::geom_col(width = 0.7) +
    ggplot2::facet_wrap(~flab, scales = "free_y", ncol = 1) +
    ggplot2::scale_fill_manual(values = setNames(FACTOR_COLS[names(.fac_lab)], .fac_lab), guide = "none") +
    ggplot2::scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.4), guide = "none") +
    ggplot2::labs(x = expression(-log[10](italic(P))~"(MAGMA gene-set)"), y = NULL) + theme_nature(base_family = figure_font(config), grid = "xy") +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 4))
  with_source_data(p, top[, .(factor, set = FULL_NAME, NGENES, BETA, SE, P, q)])
}
