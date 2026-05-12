suppressPackageStartupMessages({
  library(ggplot2)
})

plot_manhattan <- function(df, factor_name, sex, sig_line = 5e-8, sugg_line = 1e-5,
                            keep_threshold = 1e-3, bg_keep_frac = 0.1, seed = 1L) {
  d <- as.data.table(df)
  d <- d[!is.na(P) & !is.na(CHR) & !is.na(BP) & P > 0]
  if (nrow(d) == 0L) return(NULL)

  # Downsample non-significant SNPs to keep rendering tractable at 1M+ points.
  # Every SNP with P < keep_threshold is retained; a uniform random fraction
  # bg_keep_frac of the rest is kept. Signal & threshold lines are preserved.
  if (!is.null(bg_keep_frac) && bg_keep_frac < 1) {
    sig_mask <- d$P < keep_threshold
    set.seed(seed)
    bg_keep <- !sig_mask & runif(nrow(d)) < bg_keep_frac
    d <- d[sig_mask | bg_keep]
  }

  d[, CHR := as.character(CHR)]
  chr_levels <- unique(d$CHR)
  chr_levels <- chr_levels[order(suppressWarnings(as.integer(chr_levels)),
                                  chr_levels, na.last = TRUE)]
  d[, CHR := factor(CHR, levels = chr_levels)]
  d[, neg_log10_p := -log10(P)]
  d[, chr_int := as.integer(CHR)]

  chr_lens <- d[, .(max_bp = max(BP)), by = CHR][order(CHR)]
  chr_lens[, offset := cumsum(c(0, head(max_bp, -1)))]
  d <- merge(d, chr_lens[, .(CHR, offset)], by = "CHR", sort = FALSE)
  d[, x := BP + offset]
  mids <- d[, .(mid = mean(range(x))), by = CHR]

  alt_color <- c("#1f77b4", "#aec7e8")
  d[, color := alt_color[(chr_int %% 2L) + 1L]]

  ggplot(d, aes(x = x, y = neg_log10_p)) +
    geom_point(aes(color = color), size = 0.6, alpha = 0.75) +
    scale_color_identity() +
    scale_x_continuous(breaks = mids$mid, labels = mids$CHR, expand = c(0.01, 0)) +
    geom_hline(yintercept = -log10(sig_line), color = "red", linetype = "dashed") +
    geom_hline(yintercept = -log10(sugg_line), color = "blue", linetype = "dotted") +
    labs(x = "Chromosome",
         y = expression(-log[10](italic(p))),
         title = sprintf("Manhattan — %s factor (%s)", factor_name, sex)) +
    theme_minimal(base_size = 11) +
    theme(panel.grid.minor.x = element_blank(),
          panel.grid.major.x = element_blank())
}

plot_qq <- function(df, factor_name, sex) {
  d <- as.data.table(df)
  d <- d[!is.na(P) & P > 0 & P < 1]
  if (nrow(d) == 0L) return(NULL)
  observed <- -log10(sort(d$P))
  n <- length(observed)
  expected <- -log10(ppoints(n))
  lambda <- compute_lambda_gc(d$P)
  qq_df <- data.frame(expected = expected, observed = observed)

  ggplot(qq_df, aes(x = expected, y = observed)) +
    geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
    geom_point(size = 0.6, alpha = 0.75, color = "#1f77b4") +
    labs(x = expression(Expected~-log[10](italic(p))),
         y = expression(Observed~-log[10](italic(p))),
         title = sprintf("QQ — %s factor (%s)", factor_name, sex),
         subtitle = sprintf("λ_GC = %.3f, %s SNPs",
                            lambda, format(n, big.mark = ","))) +
    theme_minimal(base_size = 11)
}

plot_miami <- function(male_df, female_df, factor_name, sig_line = 5e-8,
                        keep_threshold = 1e-3, bg_keep_frac = 0.1, seed = 1L) {
  m <- as.data.table(male_df)[!is.na(P) & P > 0 & !is.na(CHR) & !is.na(BP)]
  f <- as.data.table(female_df)[!is.na(P) & P > 0 & !is.na(CHR) & !is.na(BP)]
  if (nrow(m) == 0L || nrow(f) == 0L) return(NULL)
  m[, sex := "male"]; f[, sex := "female"]
  d <- rbind(m, f, fill = TRUE)

  # Downsample non-significant SNPs identically across sexes; preserves all
  # significant signal in both panels.
  if (!is.null(bg_keep_frac) && bg_keep_frac < 1) {
    sig_mask <- d$P < keep_threshold
    set.seed(seed)
    bg_keep <- !sig_mask & runif(nrow(d)) < bg_keep_frac
    d <- d[sig_mask | bg_keep]
  }
  d[, CHR := as.character(CHR)]
  chr_levels <- unique(d$CHR)
  chr_levels <- chr_levels[order(suppressWarnings(as.integer(chr_levels)),
                                  chr_levels, na.last = TRUE)]
  d[, CHR := factor(CHR, levels = chr_levels)]
  chr_lens <- d[, .(max_bp = max(BP)), by = CHR][order(CHR)]
  chr_lens[, offset := cumsum(c(0, head(max_bp, -1)))]
  d <- merge(d, chr_lens[, .(CHR, offset)], by = "CHR", sort = FALSE)
  d[, x := BP + offset]
  d[, y := ifelse(sex == "male", -log10(P), log10(P))]
  d[, color := ifelse(as.integer(CHR) %% 2L == 1L, "#1f77b4", "#aec7e8")]
  mids <- d[, .(mid = mean(range(x))), by = CHR]

  ymax <- max(abs(d$y), na.rm = TRUE)
  ggplot(d, aes(x = x, y = y)) +
    geom_point(aes(color = color), size = 0.5, alpha = 0.7) +
    scale_color_identity() +
    scale_x_continuous(breaks = mids$mid, labels = mids$CHR, expand = c(0.01, 0)) +
    scale_y_continuous(
      breaks = scales::pretty_breaks(8),
      labels = function(y) sprintf("%.0f", abs(y)),
      limits = c(-ymax, ymax)) +
    geom_hline(yintercept = c(-log10(sig_line), log10(sig_line)),
               color = "red", linetype = "dashed") +
    geom_hline(yintercept = 0, color = "grey50") +
    annotate("text", x = -Inf, y = ymax * 0.95, hjust = -0.1, label = "male") +
    annotate("text", x = -Inf, y = -ymax * 0.95, hjust = -0.1, label = "female") +
    labs(x = "Chromosome",
         y = expression(-log[10](italic(p))),
         title = sprintf("Miami plot — %s factor (male up, female down)", factor_name)) +
    theme_minimal(base_size = 11) +
    theme(panel.grid.minor.x = element_blank(),
          panel.grid.major.x = element_blank())
}

plot_rg_heatmap <- function(S, sex, retained = NULL) {
  if (!is.null(retained)) S <- S[retained, retained, drop = FALSE]
  d <- sqrt(diag(S))
  rg <- S / outer(d, d)
  diag(rg) <- 1
  rg_df <- data.table(
    trait_y = rep(rownames(rg), times = ncol(rg)),
    trait_x = rep(colnames(rg), each = nrow(rg)),
    rg = as.vector(rg)
  )
  rg_df[, trait_x := factor(trait_x, levels = colnames(rg))]
  rg_df[, trait_y := factor(trait_y, levels = rev(rownames(rg)))]
  ggplot(rg_df, aes(x = trait_x, y = trait_y, fill = rg)) +
    geom_tile(color = "white") +
    geom_text(aes(label = sprintf("%.2f", rg)), size = 3) +
    scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
                         midpoint = 0, limits = c(-1, 1)) +
    coord_fixed() +
    labs(x = NULL, y = NULL, fill = "rg",
         title = sprintf("Genetic-correlation heatmap (%s)", sex)) +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

plot_loadings_bar <- function(loadings_df, sex) {
  d <- as.data.table(loadings_df)
  if (!all(c("trait", "factor", "loading") %in% names(d))) {
    stop("plot_loadings_bar: expected columns trait, factor, loading")
  }
  ggplot(d, aes(x = reorder(trait, loading), y = loading, fill = factor)) +
    geom_col() +
    coord_flip() +
    facet_wrap(~ factor, scales = "free_y") +
    labs(x = NULL, y = "Standardised loading",
         title = sprintf("Factor loadings (%s)", sex)) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "none")
}

plot_loading_diff_forest <- function(comparison_df) {
  d <- as.data.table(comparison_df)
  d[, ci_lo := diff - 1.96 * se_diff]
  d[, ci_hi := diff + 1.96 * se_diff]
  d[, label := sprintf("%s [%s]", trait, factor)]
  d[, sig := !is.na(p_diff) & p_diff < 0.05]

  ggplot(d, aes(x = diff, y = reorder(label, diff))) +
    geom_vline(xintercept = 0, color = "grey50", linetype = "dashed") +
    geom_pointrange(aes(xmin = ci_lo, xmax = ci_hi, color = sig)) +
    scale_color_manual(values = c(`TRUE` = "#b2182b", `FALSE` = "#3b3b3b"),
                       labels = c(`TRUE` = "p < 0.05", `FALSE` = "n.s."),
                       name = NULL) +
    labs(x = "Female loading - Male loading (unstandardised, with 95% CI)",
         y = NULL,
         title = "Sex differences in factor loadings (Wald Z-test)") +
    theme_minimal(base_size = 11)
}

dot_path_diagram <- function(apriori_model, loadings_df, sex) {
  factor_names <- parse_factor_names(apriori_model)
  ld <- as.data.table(loadings_df)
  if (!all(c("trait", "factor", "loading") %in% names(ld))) {
    stop("dot_path_diagram: loadings_df needs trait, factor, loading columns")
  }
  ld <- ld[factor %in% factor_names]

  dot <- character(0)
  dot <- c(dot,
           sprintf("digraph G_%s {", sex),
           "  rankdir=LR;",
           "  node [shape=ellipse, style=filled, fillcolor=\"#cfe0f5\"];",
           "  { rank=same;")
  for (f in factor_names) dot <- c(dot, sprintf("    \"%s\";", f))
  dot <- c(dot, "  }",
           "  node [shape=box, style=filled, fillcolor=\"#fcecc6\"];")
  traits <- unique(ld$trait)
  for (t in traits) dot <- c(dot, sprintf("  \"%s\";", t))
  for (i in seq_len(nrow(ld))) {
    dot <- c(dot,
             sprintf("  \"%s\" -> \"%s\" [label=\"%.2f\", fontsize=10];",
                     ld$factor[i], ld$trait[i], ld$loading[i]))
  }
  dot <- c(dot,
           sprintf("  labelloc=\"t\";"),
           sprintf("  label=\"GSEM a priori path diagram (%s)\";", sex),
           "}")
  paste(dot, collapse = "\n")
}

write_path_diagram <- function(apriori_model, loadings_df, sex, out_dir) {
  dot <- dot_path_diagram(apriori_model, loadings_df, sex)
  dot_path <- file.path(out_dir, sprintf("%s_path_diagram.dot", sex))
  writeLines(dot, dot_path)
  png_path <- file.path(out_dir, sprintf("%s_path_diagram.png", sex))
  rendered <- tryCatch({
    if (nzchar(Sys.which("dot"))) {
      ret <- system2("dot", c("-Tpng", shQuote(dot_path), "-o", shQuote(png_path)),
                     stdout = TRUE, stderr = TRUE)
      file.exists(png_path)
    } else FALSE
  }, error = function(e) FALSE)
  list(dot = dot_path, png = if (rendered) png_path else NA_character_)
}

save_plot <- function(plot, path, width = 10, height = 6, dpi = 150) {
  if (is.null(plot)) return(NA_character_)
  ggplot2::ggsave(path, plot, width = width, height = height, dpi = dpi,
                  bg = "white")
  path
}

run_plots <- function(config, sex) {
  log_info("plots", sprintf("=== Plots stage: %s ===", sex))
  out_dir <- file.path(config$paths$output_dir, sex, "plots")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  written <- character(0)

  gwas_path <- file.path(config$paths$output_dir, sex, "gwas",
                         paste0(sex, "_userGWAS_raw.rds"))
  if (file.exists(gwas_path)) {
    gwas_result <- readRDS(gwas_path)
    snp_path <- file.path(config$paths$output_dir, sex, "sumstats",
                          paste0(sex, "_snp_sumstats.rds"))
    snp_sumstats <- if (file.exists(snp_path)) readRDS(snp_path) else NULL

    sig_line  <- config$plots$manhattan_sig_line %||% 5e-8
    sugg_line <- config$plots$manhattan_sugg_line %||% 1e-5
    keep_thr  <- config$plots$manhattan_keep_threshold %||% 1e-3
    bg_keep   <- config$plots$manhattan_bg_keep_frac %||% 0.1
    for (f in names(gwas_result)) {
      g <- normalize_gwas_columns(gwas_result[[f]], snp_sumstats)
      man <- plot_manhattan(g, factor_name = f, sex = sex,
                             sig_line = sig_line, sugg_line = sugg_line,
                             keep_threshold = keep_thr,
                             bg_keep_frac = bg_keep)
      qq  <- plot_qq(g, factor_name = f, sex = sex)
      pm <- save_plot(man, file.path(out_dir, sprintf("%s_%s_manhattan.png", sex, f)),
                       width = 12, height = 4)
      pq <- save_plot(qq, file.path(out_dir, sprintf("%s_%s_qq.png", sex, f)),
                       width = 5, height = 5)
      written <- c(written, na.omit(c(pm, pq)))
    }
  } else {
    log_warn("plots", sprintf("No userGWAS results for %s; skipping Manhattan/QQ", sex))
  }

  ldsc_path <- file.path(config$paths$output_dir, sex, "ldsc", "ldsc_full.rds")
  h2_path <- file.path(config$paths$output_dir, sex, "ldsc", "h2_qc.csv")
  if (file.exists(ldsc_path) && file.exists(h2_path)) {
    ldsc_output <- readRDS(ldsc_path)
    retained <- fread(h2_path)[pass == TRUE]$trait
    rg_plot <- tryCatch(plot_rg_heatmap(ldsc_output$S, sex = sex, retained = retained),
                        error = function(e) {
                          log_warn("plots", sprintf("rg heatmap failed: %s", e$message))
                          NULL
                        })
    pr <- save_plot(rg_plot, file.path(out_dir, sprintf("%s_rg_heatmap.png", sex)),
                     width = 7, height = 6)
    written <- c(written, na.omit(pr))
  }

  read_loadings <- function(csv_path) {
    if (!file.exists(csv_path)) return(NULL)
    raw <- fread(csv_path)
    if (!all(c("lhs", "op", "rhs") %in% names(raw))) return(NULL)
    ld <- raw[op == "=~"]
    if (nrow(ld) == 0L) return(NULL)
    est_col <- intersect(c("STD_All", "Unstand_Est", "est", "Estimate"), names(ld))[1]
    if (is.na(est_col)) return(NULL)
    data.table(trait = ld$rhs, factor = ld$lhs, loading = ld[[est_col]])
  }

  cfa_dir <- file.path(config$paths$output_dir, sex, "cfa")
  ldgs_efa <- read_loadings(file.path(cfa_dir, "loadings.csv"))
  ldgs_ap  <- read_loadings(file.path(cfa_dir, "loadings_apriori.csv"))

  # Loadings bar: prefer a priori (= the model that drives the GWAS); fall back to EFA.
  ldgs_bar <- if (!is.null(ldgs_ap)) ldgs_ap else ldgs_efa
  if (!is.null(ldgs_bar)) {
    pl <- save_plot(plot_loadings_bar(ldgs_bar, sex = sex),
                     file.path(out_dir, sprintf("%s_loadings_bar.png", sex)),
                     width = 9, height = 6)
    written <- c(written, na.omit(pl))
  }

  if (!is.null(ldgs_ap)) {
    categories <- fread(file.path(config$paths$meta_dir, "icd10_categories.csv"))
    retained_traits <- unique(ldgs_ap$trait)
    if (length(retained_traits) > 0L) {
      apriori_model <- build_apriori_model(
        retained_traits, categories,
        min_indicators = config$cfa$min_indicators_per_factor)
      if (nchar(apriori_model) > 0L) {
        pd <- write_path_diagram(apriori_model, ldgs_ap, sex, out_dir)
        written <- c(written, pd$dot, na.omit(pd$png))
      }
    }
  }

  write_stage_manifest("plots", sex, config, written)
  log_info("plots", sprintf("Wrote %d plot file(s) for %s", length(written), sex))
  list(files = written)
}

run_comparison_plot <- function(config) {
  comp_path <- file.path(config$paths$output_dir, "comparison", "loading_diff.csv")
  if (!file.exists(comp_path)) return(invisible(NULL))
  comparison <- fread(comp_path)
  if (nrow(comparison) == 0L) return(invisible(NULL))
  out_dir <- file.path(config$paths$output_dir, "comparison")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  p <- plot_loading_diff_forest(comparison)
  save_plot(p, file.path(out_dir, "sex_diff_forest.png"), width = 8, height = 6)
}

run_miami <- function(config) {
  male_path <- file.path(config$paths$output_dir, "male", "gwas", "male_userGWAS_raw.rds")
  female_path <- file.path(config$paths$output_dir, "female", "gwas", "female_userGWAS_raw.rds")
  if (!file.exists(male_path) || !file.exists(female_path)) return(invisible(NULL))
  m_all <- readRDS(male_path)
  f_all <- readRDS(female_path)
  out_dir <- file.path(config$paths$output_dir, "comparison")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  written <- character(0)
  sig_line <- config$plots$manhattan_sig_line %||% 5e-8
  keep_thr <- config$plots$manhattan_keep_threshold %||% 1e-3
  bg_keep  <- config$plots$manhattan_bg_keep_frac %||% 0.1
  for (f in intersect(names(m_all), names(f_all))) {
    m <- normalize_gwas_columns(m_all[[f]], NULL)
    fm <- normalize_gwas_columns(f_all[[f]], NULL)
    p <- plot_miami(m, fm, factor_name = f, sig_line = sig_line,
                     keep_threshold = keep_thr, bg_keep_frac = bg_keep)
    out <- save_plot(p, file.path(out_dir, sprintf("miami_%s.png", f)),
                      width = 12, height = 6)
    written <- c(written, na.omit(out))
  }
  log_info("plots", sprintf("Wrote %d Miami plot(s)", length(written)))
  written
}
