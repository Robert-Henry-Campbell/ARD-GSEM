# Publication figures (Nature format): a 7-pt sans theme, colour-blind-safe palettes, an
# exporter writing PDF (cairo) + 600-dpi PNG + a source-data CSV per figure, and a runner that
# assembles multi-panel figures from panel functions fig_<name>(config, sex) (R/32_figures_panels.R).
# Panels whose inputs do not exist yet return NULL and are recorded as skipped, so the stage can
# be rerun as later stages finish. Widths: 89 mm (single column), 120 mm, 183 mm (double);
# height <= 247 mm.

OKABE_ITO <- c(orange = "#E69F00", sky = "#56B4E9", green = "#009E73", yellow = "#F0E442", blue = "#0072B2",
               vermilion = "#D55E00", purple = "#CC79A7", black = "#000000")
FACTOR_COLS <- c(ARD = "#0072B2", Cardiometabolic = "#D55E00", UpperGI = "#009E73", Osteoarthritis = "#898781")
GROUP_COLS <- c(SouthAsian = "#E69F00", Black = "#CC79A7", FinnGen = "#56B4E9")
MANHATTAN_COLS <- c(odd = "#184f95", even = "#5598e7", mhc = "#eb6834")
INK <- c(primary = "#0b0b0b", secondary = "#52514e", muted = "#898781", grid = "#e1e0d9", band = "#f4f3ef")
NATURE_LWD <- 0.35
NATURE_WIDTHS <- c(89, 120, 183)

.fig_env <- new.env(parent = emptyenv())

figure_font <- function(config = NULL) {
  if (!is.null(.fig_env$font)) return(.fig_env$font)
  fam <- config$figures$font_family %||% "Arial"
  have <- tryCatch({
    out <- suppressWarnings(system2("fc-list", stdout = TRUE, stderr = FALSE))
    any(grepl(fam, out, fixed = TRUE))
  }, error = function(e) FALSE)
  if (!have) { log_warn("figures", sprintf("font '%s' not found by fc-list; using 'sans'", fam)); fam <- "sans" }
  .fig_env$font <- fam
  fam
}

theme_nature <- function(base_size = 7, base_family = "sans", grid = c("none", "y", "xy")) {
  grid <- match.arg(grid)
  gl <- ggplot2::element_line(colour = INK[["grid"]], linewidth = 0.25)
  ggplot2::theme_minimal(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      text = ggplot2::element_text(size = base_size, colour = INK[["primary"]]),
      axis.text = ggplot2::element_text(size = base_size - 1, colour = INK[["secondary"]]),
      axis.title = ggplot2::element_text(size = base_size),
      axis.line = ggplot2::element_line(colour = INK[["primary"]], linewidth = 0.25),
      axis.ticks = ggplot2::element_line(colour = INK[["primary"]], linewidth = 0.25),
      axis.ticks.length = grid::unit(1.5, "pt"),
      legend.text = ggplot2::element_text(size = base_size - 1),
      legend.title = ggplot2::element_text(size = base_size),
      legend.key.size = grid::unit(3, "mm"),
      legend.margin = ggplot2::margin(0, 0, 0, 0),
      legend.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(size = base_size, face = "bold", hjust = 0),
      strip.background = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(size = base_size, face = "bold"),
      plot.subtitle = ggplot2::element_text(size = base_size - 1, colour = INK[["secondary"]]),
      plot.caption = ggplot2::element_text(size = base_size - 2, colour = INK[["muted"]]),
      plot.tag = ggplot2::element_text(size = base_size + 1, face = "bold"),
      plot.tag.position = c(0, 1),
      plot.margin = ggplot2::margin(2, 2, 2, 2, "mm"),
      plot.background = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.x = if (grid == "xy") gl else ggplot2::element_blank(),
      panel.grid.major.y = if (grid %in% c("y", "xy")) gl else ggplot2::element_blank()
    )
}

# Attach the data behind a panel (a data.frame, or a named list of them) for the source-data file.
with_source_data <- function(plot, data) { attr(plot, "source_data") <- data; plot }

save_nature_fig <- function(plot, name, width_mm, height_mm, out_dir, dpi = 600,
                            source_data = attr(plot, "source_data")) {
  if (!width_mm %in% NATURE_WIDTHS) log_warn("figures", sprintf("%s: width %s mm is not a Nature column width (89/120/183)", name, width_mm))
  if (height_mm > 247) log_warn("figures", sprintf("%s: height %s mm exceeds 247 mm", name, height_mm))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  pdf_path <- file.path(out_dir, paste0(name, ".pdf")); png_path <- file.path(out_dir, paste0(name, ".png"))
  ggplot2::ggsave(pdf_path, plot, width = width_mm, height = height_mm, units = "mm", device = grDevices::cairo_pdf, bg = "white")
  ggplot2::ggsave(png_path, plot, width = width_mm, height = height_mm, units = "mm", dpi = dpi, bg = "white")
  files <- c(pdf_path, png_path)
  if (!is.null(source_data)) {
    if (is.data.frame(source_data)) source_data <- list(source_data)
    for (k in seq_along(source_data)) {
      tag <- names(source_data)[k]
      if (is.null(tag) || !nzchar(tag)) tag <- if (length(source_data) == 1L) "" else paste0("_", letters[k])
      else tag <- paste0("_", tag)
      p <- file.path(out_dir, paste0(name, "__source_data", tag, ".csv"))
      fwrite(as.data.table(source_data[[k]]), p); files <- c(files, p)
    }
  }
  files
}

# Existing stage outputs needed by a panel; NULL (with a log line) when any is missing.
fig_inputs <- function(config, sex, stage, files) {
  p <- file.path(stage_dir(config, sex, stage), files)
  miss <- !file.exists(p)
  if (any(miss)) { log_info("figures", sprintf("panel skipped: %s/%s missing", stage, paste(files[miss], collapse = ", "))); return(NULL) }
  setNames(as.list(p), files)
}

# entry: list(name, width_mm, height_mm, ncol, panels = character vector of fig_* function names)
assemble_figure <- function(entry, config, sex) {
  panels <- list(); status <- list()
  for (fn in entry$panels) {
    p <- if (exists(fn, mode = "function")) tryCatch(get(fn, mode = "function")(config, sex), error = function(e) {
      log_warn("figures", sprintf("%s: %s", fn, conditionMessage(e))); structure(list(), class = "fig_error", msg = conditionMessage(e)) })
    else structure(list(), class = "fig_error", msg = "no such panel function")
    if (inherits(p, "fig_error")) status[[fn]] <- paste("error:", attr(p, "msg"))
    else if (is.null(p)) status[[fn]] <- "skipped: inputs missing"
    else { panels[[fn]] <- p; status[[fn]] <- "ok" }
  }
  if (!length(panels)) return(list(plot = NULL, status = status, source_data = NULL))
  sd <- lapply(panels, function(p) attr(p, "source_data"))
  names(sd) <- paste0(letters[seq_along(panels)], "_", sub("^fig_", "", names(panels)))
  sd <- sd[!vapply(sd, is.null, logical(1))]
  if (length(panels) == 1L) {
    plot <- panels[[1]] + ggplot2::labs(tag = NULL)
  } else if (requireNamespace("patchwork", quietly = TRUE)) {
    design <- entry$design   # e.g. "AB\nCD\nEE" lets the last panel span the page; only valid when no panel was skipped
    if (!is.null(design) && length(unique(strsplit(gsub("[^A-Za-z]", "", design), "")[[1]])) != length(panels)) design <- NULL
    plot <- (if (is.null(design)) patchwork::wrap_plots(panels, ncol = entry$ncol %||% 2) else patchwork::wrap_plots(panels, design = design)) +
      patchwork::plot_annotation(tag_levels = "a") &
      ggplot2::theme(plot.tag = ggplot2::element_text(size = 8, face = "bold"))
  } else {
    log_warn("figures", sprintf("%s: patchwork not installed; panels saved separately", entry$name))
    plot <- panels   # list -> saved one by one by run_figures
  }
  list(plot = plot, status = status, source_data = sd)
}

run_figures <- function(config, sex) {
  stage <- "figures"
  log_info(stage, sprintf("=== Figures (Nature format): %s ===", sex))
  cfg <- config$figures %||% list()
  out_dir <- stage_dir(config, sex, stage)
  dir.create(file.path(out_dir, "extended"), recursive = TRUE, showWarnings = FALSE)
  .fig_env$font <- NULL
  only <- unlist(cfg$only %||% list())
  index <- list(); files <- character(0)
  for (set in c("main", "extended")) {
    for (entry in cfg[[set]] %||% list()) {
      if (length(only) && !entry$name %in% only) next
      dir <- if (set == "main") out_dir else file.path(out_dir, "extended")
      a <- assemble_figure(entry, config, sex)
      for (fn in names(a$status)) index[[length(index) + 1]] <- data.table(figure = entry$name, panel = fn, status = a$status[[fn]])
      if (is.null(a$plot)) { log_info(stage, sprintf("%s: nothing to draw yet", entry$name)); next }
      safe_save <- function(...) tryCatch(save_nature_fig(...), error = function(e) {
        log_warn(stage, sprintf("%s: drawing failed: %s", entry$name, conditionMessage(e))); character(0) })
      written <- if (is.list(a$plot) && !inherits(a$plot, "ggplot")) {
        unlist(lapply(names(a$plot), function(fn) safe_save(a$plot[[fn]], paste0(entry$name, "_", sub("^fig_", "", fn)),
                                                            entry$width_mm %||% 183, (entry$height_mm %||% 120) / max(1, length(a$plot) / (entry$ncol %||% 2)),
                                                            dir, dpi = cfg$dpi %||% 600, source_data = attr(a$plot[[fn]], "source_data"))))
      } else safe_save(a$plot, entry$name, entry$width_mm %||% 183, entry$height_mm %||% 120, dir,
                       dpi = cfg$dpi %||% 600, source_data = a$source_data)
      files <- c(files, written)
      log_info(stage, sprintf("%s: %d/%d panels drawn", entry$name, sum(unlist(a$status) == "ok"), length(a$status)))
    }
  }
  idx <- rbindlist(index)
  fwrite(idx, file.path(out_dir, "figure_index.tsv"), sep = "\t")
  write_stage_manifest(stage, sex, config, c(files, file.path(out_dir, "figure_index.tsv")))
  invisible(idx)
}
