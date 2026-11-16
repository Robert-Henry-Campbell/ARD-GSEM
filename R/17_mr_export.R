# MR export stage: Q_SNP-filtered, clumped instruments per factor in the exposure format of
# the ardmr package (ARD-Mendelian-Randomization). MR itself runs in ardmr.
#
# Selection: qfilter tier (GenomicSEM QC: no Q_SNP hits or their LD partners), P < p_gws,
# autosomal, MHC excluded (configurable), effect-allele frequency known, F >= f_stat_min,
# ambiguous palindromes (A/T, C/G with EAF 0.42-0.58, ardmr's rule) dropped, then plink
# clumping r2 < clump_r2 within clump_kb on 1000G EUR (ardmr's defaults).
# Effects are per SD of the latent factor (std.lv); samplesize = N_hat_std.

palindrome <- function(a1, a2) {
  p <- paste0(toupper(a1), toupper(a2))
  p %in% c("AT", "TA", "CG", "GC")
}

ARDMR_EXPOSURE_COLS <- c("id.exposure", "SNP", "beta.exposure", "se.exposure",
                         "effect_allele.exposure", "other_allele.exposure", "pval.exposure",
                         "eaf.exposure", "samplesize.exposure", "Chr", "Pos")

# ardmr's own column check when the package source is available, else the same rule here.
validate_ardmr_exposure <- function(x, ardmr_dir) {
  f <- file.path(ardmr_dir %||% "", "R", "utils-data.R")
  if (nzchar(ardmr_dir %||% "") && file.exists(f)) {
    e <- new.env()
    sys.source(f, envir = e)
    e$assert_exposure(x)
    return("ardmr::assert_exposure")
  }
  req <- ARDMR_EXPOSURE_COLS[1:7]
  miss <- setdiff(req, names(x))
  if (length(miss)) stop("exposure missing columns: ", paste(miss, collapse = ", "))
  "built-in column check"
}

run_mr_export <- function(config, sex) {
  stage <- "mr_export"
  log_info(stage, sprintf("=== MR instrument export: %s ===", sex))
  pg <- config$postgwas %||% list(); cfg <- config$mr_export %||% list()
  out_dir <- stage_dir(config, sex, stage); work <- file.path(out_dir, "work")
  dir.create(work, recursive = TRUE, showWarnings = FALSE)
  p_gws <- as.numeric(pg$p_gws %||% 5e-8)
  reg <- read_factor_registry(config, sex)
  files <- character(0); summ <- list()

  for (f in postgwas_factors(config)) {
    dt <- load_postgwas(config, sex, f, "qfilter",
                        cols = c("SNP", "CHR", "BP", "A1", "A2", "EAF_A1", "beta", "SE", "P",
                                 "LOG10P", "is_mhc", "in_1000g"))
    n0 <- sum(dt$P < p_gws)
    cand <- dt[P < p_gws & CHR <= 22L & !is.na(EAF_A1)]
    if (isTRUE(cfg$exclude_mhc %||% TRUE)) cand <- cand[is_mhc == FALSE]
    cand[, F_stat := (beta / SE)^2]
    cand <- cand[F_stat >= as.numeric(cfg$f_stat_min %||% 10)]
    ambig <- palindrome(cand$A1, cand$A2) & cand$EAF_A1 >= 0.42 & cand$EAF_A1 <= 0.58
    cand <- cand[!ambig]
    not_in_panel <- sum(!cand$in_1000g)
    cand <- cand[in_1000g == TRUE]
    log_info(stage, sprintf("%s: %d qfilter GWS SNPs -> %d candidates (%d ambiguous palindromes, %d not in 1000G dropped)",
                            f, n0, nrow(cand), sum(ambig), not_in_panel))
    if (!nrow(cand)) {
      log_warn(stage, sprintf("%s: no instrument candidates", f))
      next
    }
    cand[, P := pmax(P, 10^-pmin(LOG10P, 300))]
    cl <- plink_clump(config, config$refs$ld_1000g_bfile, cand, p1 = p_gws, p2 = 1,
                      r2 = as.numeric(cfg$clump_r2 %||% 0.001), kb = as.integer(cfg$clump_kb %||% 10000L),
                      out_prefix = file.path(work, paste0(f, "_clump")), stage = stage)
    iv <- cand[SNP %in% cl$SNP][order(CHR, BP)]
    n_hat <- reg[[f]]$N_hat_std
    out <- iv[, .(id.exposure = sprintf("ARD_GSEM_%s_%s", sex, f), SNP, beta.exposure = beta,
                  se.exposure = SE, effect_allele.exposure = A1, other_allele.exposure = A2,
                  pval.exposure = P, eaf.exposure = EAF_A1, samplesize.exposure = round(n_hat),
                  Chr = CHR, Pos = BP)]
    setcolorder(out, ARDMR_EXPOSURE_COLS)
    check <- validate_ardmr_exposure(out, config$refs$ardmr_dir)
    path <- file.path(out_dir, paste0(f, "_instruments.tsv"))
    fwrite(out, path, sep = "\t")
    meta <- list(factor = f, n_instruments = nrow(out), validated_with = check,
                 scale = "per SD of the latent genetic factor (GenomicSEM std.lv = TRUE)",
                 samplesize = "N_hat_std (GenomicSEM wiki N_hat formula on the std.lv SEs)",
                 marker = reg[[f]]$marker, marker_lambda = reg[[f]]$marker_lambda,
                 N_hat_marker = reg[[f]]$N_hat_marker,
                 note_marker = "beta per SD of the marker indicator's liability = beta * marker_lambda",
                 filters = list(tier = "qfilter (Q_SNP p >= qsnp_p and not in LD with a Q_SNP hit)",
                                p = p_gws, autosomal = TRUE, exclude_mhc = isTRUE(cfg$exclude_mhc %||% TRUE),
                                f_stat_min = cfg$f_stat_min %||% 10,
                                ambiguous_palindromes = "dropped (EAF 0.42-0.58)",
                                clump = sprintf("plink r2 < %s within %s kb, 1000G EUR",
                                                cfg$clump_r2 %||% 0.001, cfg$clump_kb %||% 10000)),
                 eaf_source = "Pan-UKB EUR allele frequency of A1 (the effect allele)")
    jsonlite::write_json(meta, sub("\\.tsv$", ".meta.json", path), auto_unbox = TRUE, pretty = TRUE, digits = NA)
    files <- c(files, path)
    summ[[f]] <- data.table(factor = f, n_gws_qfilter = n0, n_candidates = nrow(cand),
                            n_instruments = nrow(out), median_F = median((out$beta.exposure / out$se.exposure)^2))
    log_info(stage, sprintf("%s: %d instruments -> %s (%s)", f, nrow(out), basename(path), check))
  }
  fwrite(rbindlist(summ), file.path(out_dir, "mr_export_summary.csv"))
  write_stage_manifest(stage, sex, config, list.files(out_dir, full.names = TRUE))
  invisible(files)
}
