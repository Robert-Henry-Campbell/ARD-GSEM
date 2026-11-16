suppressPackageStartupMessages({
  library(ggplot2)
})

# PLINK numeric codes 23-26 -> X, Y, XY, MT; strips a "chr" prefix.
chrom_labels <- function(chr) {
  chr <- sub("^chr", "", as.character(chr), ignore.case = TRUE)
  plink <- c(`23` = "X", `24` = "Y", `25` = "XY", `26` = "MT")
  hit <- chr %in% names(plink)
  chr[hit] <- plink[chr[hit]]
  chr
}

# -log10(P) that survives underflow: where P rounded to 0, recompute it from Z in log space.
neg_log10_p <- function(p, z = NULL) {
  out <- -log10(p)
  if (!is.null(z)) {
    uf <- !is.na(p) & p == 0 & is.finite(z)
    out[uf] <- -(pnorm(-abs(z[uf]), log.p = TRUE) + log(2)) / log(10)
  }
  out
}

# Greedy distance clumping on columns CHR, BP, lp: the strongest remaining SNP leads a
# locus and absorbs every SNP within +-window bp. SNPs inside the mhc window are first
# collapsed to their single strongest, so the extended MHC counts as one locus.
distance_clump <- function(d, window = 5e5, mhc = NULL) {
  d <- as.data.table(d)[order(-lp)]
  if (!is.null(mhc)) {
    in_mhc <- d$CHR == as.character(mhc$chr) & d$BP >= mhc$start & d$BP <= mhc$end
    d <- d[!in_mhc | !duplicated(in_mhc)]
  }
  lead <- taken <- logical(nrow(d))
  for (i in seq_len(nrow(d))) {
    if (taken[i]) next
    lead[i] <- TRUE
    taken <- taken | (d$CHR == d$CHR[i] & abs(d$BP - d$BP[i]) <= window)
  }
  d[lead]
}

# Linear up to `at`, then slope k above it: keeps one or two extreme peaks from
# flattening every other locus on the Manhattan.
compress_above_transform <- function(at, k) {
  scales::new_transform(
    "compress_above",
    transform = function(y) ifelse(y > at, at + (y - at) * k, y),
    inverse   = function(y) ifelse(y > at, at + (y - at) / k, y))
}

# "5 %*% 10^-8" for plotmath and "5×10⁻⁸" for plain text.
p_label <- function(p, plotmath = TRUE) {
  e <- floor(log10(p) + 1e-9)
  m <- signif(p / 10^e, 2)
  if (plotmath) {
    return(if (m == 1) sprintf("10^%d", e) else sprintf("%s %%*%% 10^%d", format(m), e))
  }
  sup <- chartr("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹", as.character(e))
  if (m == 1) paste0("10", sup) else paste0(format(m), "×10", sup)
}

parse_indicators <- function(model_str) {
  lines <- grep("=~", strsplit(model_str, "\n", fixed = TRUE)[[1L]], fixed = TRUE, value = TRUE)
  rhs <- strsplit(sub("^[^=]*=~", "", lines), "+", fixed = TRUE)
  setNames(lapply(rhs, function(x) sub("^.*\\*", "", trimws(x))),
           trimws(sub("\\s*=~.*$", "", lines)))
}

# df: normalize_gwas_columns() output. Every SNP with P < keep_threshold is drawn; the
# rest are thinned to one point per output pixel (width x height inches at dpi), which
# is visually lossless for opaque points but cuts ~8M SNPs to well under 1M.
plot_manhattan <- function(df, factor_name, sex, sig_line = 5e-8, sugg_line = 1e-5,
                            keep_threshold = 1e-3, qsnp_line = 5e-8,
                            mhc = list(chr = "6", start = 25e6, end = 35e6),
                            compress_above = 30, compress_frac = 0.25,
                            role = "factor", indicators = NULL,
                            width = 12, height = 4.5, dpi = 300) {
  ink <- c(primary = "#0b0b0b", secondary = "#52514e", muted = "#898781",
           grid = "#e1e0d9", band = "#f4f3ef")
  cols <- intersect(c("CHR", "BP", "P", "ES", "SE", "Q_SNP_pval"), names(df))
  d <- as.data.table(df)[, cols, with = FALSE]
  z <- if (all(c("ES", "SE") %in% cols)) d$ES / d$SE else NULL
  d[, lp := neg_log10_p(P, z)]
  n_no_est <- sum(!is.finite(d$lp))
  n_no_pos <- sum(is.finite(d$lp) & (is.na(d$CHR) | is.na(d$BP)))
  d <- d[is.finite(lp) & !is.na(CHR) & !is.na(BP)]
  if (nrow(d) == 0L) return(NULL)

  # Chromosomes laid end to end in genome order, with a small gap between them.
  d[, CHR := chrom_labels(CHR)]
  present <- unique(d$CHR)
  known <- c(as.character(1:22), "X", "Y", "XY", "MT")
  chr_order <- c(intersect(known, present), sort(setdiff(present, known)))
  chrs <- d[, .(len = max(BP)), by = CHR][match(chr_order, CHR)]
  gap <- 0.004 * sum(chrs$len)
  chrs[, `:=`(offset = cumsum(c(0, head(len, -1) + gap)), idx = .I)]
  d[chrs, on = "CHR", `:=`(x = BP + i.offset, cls = c("even", "odd")[i.idx %% 2L + 1L])]
  if (!is.null(mhc)) {
    mhc$chr <- chrom_labels(mhc$chr)
    d[CHR == mhc$chr & BP >= mhc$start & BP <= mhc$end, cls := "mhc"]
  }

  sig_lp <- -log10(sig_line)
  n_loci <- nrow(distance_clump(d[lp > sig_lp], mhc = mhc))
  het <- if ("Q_SNP_pval" %in% cols) {
    d[lp > sig_lp & !is.na(Q_SNP_pval) & Q_SNP_pval < qsnp_line]
  } else d[0L]
  lambda <- compute_lambda_gc(d$P)

  y_max <- max(d$lp, sig_lp + 1)
  compress <- is.finite(compress_above) && y_max > compress_above / (1 - compress_frac)
  y_disp_max <- if (compress) compress_above / (1 - compress_frac) else max(pretty(c(0, y_max)))

  px_w <- width * dpi
  px_h <- height * dpi
  keep_lp <- -log10(keep_threshold)
  bulk <- d[lp < keep_lp, .(x, lp, cls)]
  bulk[, `:=`(xb = floor(x / (diff(range(d$x)) / px_w)), yb = floor(lp / (y_disp_max / px_h)))]
  bulk <- unique(bulk, by = c("xb", "yb", "cls"))
  pts <- rbind(bulk[, .(x, lp, cls)], d[lp >= keep_lp, .(x, lp, cls)])

  if (compress) {
    k <- compress_frac / (1 - compress_frac) * compress_above / (y_max - compress_above)
    lower <- pretty(c(0, compress_above))
    lower <- lower[lower <= compress_above]
    # Upper ticks: round values at least a tenth of the panel apart on the display scale.
    vis_top <- 1.04 * y_disp_max
    upper <- numeric(0)
    last <- compress_above
    for (b in pretty(c(compress_above, compress_above + (vis_top - compress_above) / k), n = 6)) {
      if (b > compress_above && (b - last) * k >= 0.1 * y_disp_max &&
          compress_above + (b - compress_above) * k <= vis_top) {
        upper <- c(upper, b)
        last <- b
      }
    }
    y_scale <- scale_y_continuous(transform = compress_above_transform(compress_above, k),
                                  breaks = c(lower, upper),
                                  expand = expansion(mult = c(0.01, 0.04)))
  } else {
    y_scale <- scale_y_continuous(limits = c(0, y_disp_max), breaks = pretty(c(0, y_disp_max)),
                                  expand = expansion(mult = c(0.01, 0.04)))
  }

  subtitle <- paste(c(
    if (length(indicators) > 12L) sprintf("%d indicators", length(indicators))
    else if (length(indicators) > 0L) sprintf("indicators: %s", paste(indicators, collapse = ", ")),
    sprintf("λGC %.2f", lambda),
    sprintf("%d loci at P < %s", n_loci, p_label(sig_line, plotmath = FALSE))),
    collapse = "  ·  ")
  caption <- paste0(
    sprintf("Genomic SEM factor GWAS, %s  ·  %s SNPs", sex, format(nrow(d), big.mark = ",")),
    if (n_no_est > 0L) sprintf("; %s failed fits omitted", format(n_no_est, big.mark = ",")),
    if (n_no_pos > 0L) sprintf("; %s without a position omitted", format(n_no_pos, big.mark = ",")),
    "  ·  loci: greedy ±500 kb clumping, MHC counted once",
    if (compress) sprintf("  ·  shaded band: y-axis compressed above %g", compress_above))

  p <- ggplot(pts, aes(x = x, y = lp))
  if (compress) {
    p <- p + annotate("rect", xmin = -Inf, xmax = Inf, ymin = compress_above, ymax = Inf,
                      fill = ink[["band"]]) +
      geom_hline(yintercept = upper, colour = ink[["grid"]], linewidth = 0.3)
  }
  p <- p +
    geom_point(aes(colour = cls), shape = 16, size = 0.6, stroke = 0) +
    scale_colour_manual(
      values = c(odd = "#184f95", even = "#5598e7", mhc = "#eb6834"), breaks = "mhc",
      labels = sprintf("MHC (chr%s:%g–%g Mb)", mhc$chr, mhc$start / 1e6, mhc$end / 1e6),
      name = NULL, guide = guide_legend(order = 1, override.aes = list(size = 2.2)))
  if (nrow(het) > 0L) {
    p <- p +
      geom_point(data = het, aes(shape = "het"), colour = ink[["primary"]],
                 size = 1.3, stroke = 0.3) +
      scale_shape_manual(
        values = c(het = 5), name = NULL,
        labels = parse(text = sprintf('italic(Q)[SNP]~italic(P) < %s~"(%s SNPs)"',
                                      p_label(qsnp_line), format(nrow(het), big.mark = ","))),
        guide = guide_legend(order = 2, override.aes = list(size = 2, stroke = 0.5)))
  }
  p +
    geom_hline(yintercept = sig_lp, colour = ink[["secondary"]], linetype = "dashed",
               linewidth = 0.35) +
    geom_hline(yintercept = -log10(sugg_line), colour = ink[["muted"]], linetype = "dotted",
               linewidth = 0.35) +
    annotate("text", x = Inf, y = c(sig_lp, -log10(sugg_line)),
             label = c(p_label(sig_line), p_label(sugg_line)), parse = TRUE,
             hjust = -0.15, size = 2.6, colour = ink[c("secondary", "muted")]) +
    scale_x_continuous(breaks = chrs$offset + chrs$len / 2, labels = chrs$CHR,
                       expand = expansion(mult = 0.004)) +
    y_scale +
    coord_cartesian(clip = "off") +
    labs(x = "Chromosome", y = expression(-log[10](italic(P))),
         title = sprintf("%s %s", factor_name, role), subtitle = subtitle, caption = caption) +
    theme_minimal(base_size = 10) +
    theme(plot.title = element_text(face = "bold", colour = ink[["primary"]], size = 13),
          plot.subtitle = element_text(colour = ink[["secondary"]], size = 9.5),
          plot.caption = element_text(colour = ink[["muted"]], size = 7.5, hjust = 0),
          plot.title.position = "plot",
          plot.caption.position = "plot",
          axis.title = element_text(colour = ink[["secondary"]], size = 9.5),
          axis.text = element_text(colour = ink[["secondary"]], size = 8),
          panel.grid.major.x = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major.y = element_line(colour = ink[["grid"]], linewidth = 0.3),
          legend.position = "top",
          legend.justification = "right",
          legend.text = element_text(colour = ink[["secondary"]], size = 8.5),
          legend.margin = margin(0, 0, 0, 0),
          legend.box.spacing = unit(2, "pt"),
          plot.margin = margin(10, 46, 8, 10))
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

# One Manhattan per latent factor -> <output_dir>/<sex>/plots/<sex>_<factor>_manhattan.png.
# Reads <sex>_userGWAS_raw.rds unless run_plots() passes the already-loaded list; factor
# roles (general/specific) and indicators come from <sex>_model_spec.rds when present.
run_manhattans <- function(config, sex, gwas_result = NULL, snp_sumstats = NULL) {
  gwas_dir <- file.path(config$paths$output_dir, sex, "gwas")
  if (is.null(gwas_result)) {
    gwas_path <- file.path(gwas_dir, paste0(sex, "_userGWAS_raw.rds"))
    if (!file.exists(gwas_path)) {
      log_warn("plots", sprintf("No userGWAS results for %s; skipping Manhattans", sex))
      return(character(0))
    }
    log_info("plots", sprintf("Reading %s", gwas_path))
    gwas_result <- readRDS(gwas_path)
  }
  # userGWAS output normally carries CHR/BP; only fall back to the SNP sumstats if not.
  if (is.null(snp_sumstats) && !all(c("CHR", "BP") %in% names(gwas_result[[1L]]))) {
    snp_path <- file.path(config$paths$output_dir, sex, "sumstats",
                          paste0(sex, "_snp_sumstats.rds"))
    if (file.exists(snp_path)) snp_sumstats <- readRDS(snp_path)
  }
  spec_path <- file.path(gwas_dir, paste0(sex, "_model_spec.rds"))
  spec <- if (file.exists(spec_path)) readRDS(spec_path) else list()
  indicators <- parse_indicators(spec$measurement_model %||% spec$model %||% "")

  out_dir <- file.path(config$paths$output_dir, sex, "plots")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  mhc <- list(chr = config$gwas$mhc_chr %||% 6,
              start = config$gwas$mhc_start %||% 25e6,
              end = config$gwas$mhc_end %||% 35e6)
  written <- character(0)
  for (f in names(gwas_result)) {
    role <- if (is.null(spec$general_name)) "factor"
            else if (identical(f, spec$general_name)) "general factor"
            else "specific factor"
    g <- normalize_gwas_columns(gwas_result[[f]], snp_sumstats)
    man <- plot_manhattan(g, factor_name = f, sex = sex,
                          sig_line = config$plots$manhattan_sig_line %||% 5e-8,
                          sugg_line = config$plots$manhattan_sugg_line %||% 1e-5,
                          keep_threshold = config$plots$manhattan_keep_threshold %||% 1e-3,
                          qsnp_line = config$qc_filters$qsnp_pval_threshold %||% 5e-8,
                          mhc = mhc,
                          compress_above = config$plots$manhattan_compress_above %||% 30,
                          role = role, indicators = indicators[[f]],
                          width = 12, height = 4.5, dpi = 300)
    pm <- save_plot(man, file.path(out_dir, sprintf("%s_%s_manhattan.png", sex, f)),
                    width = 12, height = 4.5, dpi = 300)
    if (!is.na(pm)) log_info("plots", sprintf("Manhattan (%s %s): %s", sex, f, pm))
    written <- c(written, na.omit(pm))
  }
  written
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

    written <- c(written, run_manhattans(config, sex, gwas_result, snp_sumstats))
    for (f in names(gwas_result)) {
      g <- normalize_gwas_columns(gwas_result[[f]], snp_sumstats)
      qq  <- plot_qq(g, factor_name = f, sex = sex)
      pq <- save_plot(qq, file.path(out_dir, sprintf("%s_%s_qq.png", sex, f)),
                       width = 5, height = 5)
      written <- c(written, na.omit(pq))
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
    categories <- load_trait_categories(config, sex)
    retained_traits <- unique(ldgs_ap$trait)
    if (length(retained_traits) > 0L) {
      apriori_model <- build_apriori_model(
        retained_traits, categories,
        min_indicators = config$cfa$min_indicators_per_factor,
        code_format = if (identical(sex, "bothsex_meta")) "free" else "icd3")
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
