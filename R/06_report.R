run_report <- function(config) {
  log_info("report", "=== Report stage ===")

  report_data <- list()

  for (sex in c("male", "female")) {
    h2_path <- file.path(config$paths$output_dir, sex, "ldsc", "h2_qc.csv")
    if (file.exists(h2_path)) report_data[[sex]]$h2 <- fread(h2_path)

    loadings_path <- file.path(config$paths$output_dir, sex, "efa", "loadings.csv")
    if (file.exists(loadings_path)) report_data[[sex]]$efa_loadings <- fread(loadings_path)

    fit_path <- file.path(config$paths$output_dir, sex, "cfa", "fit_comparison.csv")
    if (file.exists(fit_path)) report_data[[sex]]$cfa_fit <- fread(fit_path)

    cfa_loadings_path <- file.path(config$paths$output_dir, sex, "cfa", "loadings.csv")
    if (file.exists(cfa_loadings_path)) report_data[[sex]]$cfa_loadings <- fread(cfa_loadings_path)
  }

  comp_path <- file.path(config$paths$output_dir, "comparison", "loading_diff.csv")
  if (file.exists(comp_path)) report_data$comparison <- fread(comp_path)

  template_path <- "templates/report.Rmd"
  if (!file.exists(template_path)) {
    log_warn("report", "Report template not found; generating summary to console")
    print_summary(report_data)
    return(invisible(NULL))
  }

  output_path <- file.path(config$paths$output_dir, "report.html")
  tryCatch({
    rmarkdown::render(
      template_path,
      output_file = output_path,
      params = list(data = report_data, config = config),
      quiet = TRUE
    )
    file_size <- file.size(output_path)
    log_info("report", sprintf("Report written to %s (%s KB)",
                               output_path, round(file_size / 1024)))
  }, error = function(e) {
    log_error("report", sprintf("Report rendering failed: %s", e$message))
    log_info("report", "Falling back to console summary")
    print_summary(report_data)
  })
}

print_summary <- function(report_data) {
  cat("\n=== GSEM Pipeline Summary ===\n\n")
  for (sex in c("male", "female")) {
    if (!is.null(report_data[[sex]]$h2)) {
      cat(sprintf("--- %s h2 estimates ---\n", toupper(sex)))
      print(report_data[[sex]]$h2)
      cat("\n")
    }
    if (!is.null(report_data[[sex]]$cfa_fit)) {
      cat(sprintf("--- %s CFA fit ---\n", toupper(sex)))
      print(report_data[[sex]]$cfa_fit)
      cat("\n")
    }
  }
  if (!is.null(report_data$comparison)) {
    cat("--- Sex Comparison ---\n")
    print(report_data$comparison)
    cat("\n")
  }
}
