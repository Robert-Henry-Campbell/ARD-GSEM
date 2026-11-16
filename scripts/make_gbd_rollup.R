#!/usr/bin/env Rscript
# Builds manifest/indicator_gbd_causes.tsv: one GBD cause ("unit") per ARD indicator ICD-10
# three-character code, from the GBD -> ICD-10 explode map in the ardmr package
# (ARD-Mendelian-Randomization/data/bothsex_all.rda; the ARD-filtered map lacks half the
# indicators). Rule: the GBD level-3 cause whose ICD-10 range includes the code (level-4 matches
# are lifted to their level-3 parent), preferring causes that cover the whole three-character
# code, then the number of matching sub-code ranges; a code with no level-3 cause stays its own
# unit. Four codes are fixed by hand (see OVERRIDES). The output is meant to be hand-edited;
# rerunning this script overwrites it, so edit the OVERRIDES here instead.
suppressMessages({library(data.table); library(jsonlite)})
root <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])), ".."))
rda <- file.path(root, "..", "ARD-Mendelian-Randomization", "data", "bothsex_all.rda")
out <- file.path(root, "manifest", "indicator_gbd_causes.tsv")
reg <- fromJSON(file.path(root, "output", "bothsex", "qc_filters", "factor_registry.json"))
ind <- sort(unique(unlist(lapply(reg, `[[`, "indicators"))))

OVERRIDES <- list(
  E11 = list(cause = "Diabetes mellitus", note = "the map also lists Chronic kidney disease via E11.2; diabetes chosen"),
  K63 = list(cause = "Other neoplasms", note = "reaches Other neoplasms via K63.5 (polyp of colon); kept per plan"),
  K52 = list(cause = "Inflammatory bowel disease", note = "GBD groups non-infective gastroenteritis/colitis with IBD (K50-K52)"),
  E66 = list(cause = "Endocrine, metabolic, blood, and immune disorders", note = "also under Adverse effects of medical treatment via E66.1; endocrine chosen"))

e <- new.env(); load(rda, envir = e); x <- as.data.table(get(ls(e)[1], e))
m <- unique(x[ICD10_explo %in% ind, .(icd10 = ICD10_explo, cause = `Cause Name`, level = as.integer(`Cause Hierarchy Level`),
                                      l3 = cause_level_3, cause_id = `Cause ID`, range = ICD10)])
tok <- function(r) trimws(unlist(strsplit(r, ",")))
m[, whole := mapply(function(k, r) any(grepl(paste0("(^|-)", k, "(\\.9+)?($|-)"), tok(r))), icd10, range)]
m[, n_tok := mapply(function(k, r) sum(startsWith(tok(r), k)), icd10, range)]
cand <- m[level >= 3 & !is.na(l3)]
cand[, unit := l3]
agg <- cand[, .(whole = any(whole), n_tok = sum(n_tok), best_level = max(level)), by = .(icd10, unit)]
pick <- agg[order(icd10, -whole, -n_tok, -best_level)][, .SD[1], by = icd10]
rows <- lapply(ind, function(k) {
  if (!is.null(OVERRIDES[[k]])) {
    return(data.table(icd10 = k, gbd_cause = OVERRIDES[[k]]$cause, gbd_level = 3L, whole_code = agg[icd10 == k & unit == OVERRIDES[[k]]$cause, whole][1],
                      note = paste("manual:", OVERRIDES[[k]]$note)))
  }
  p <- pick[icd10 == k]
  if (!nrow(p)) return(data.table(icd10 = k, gbd_cause = k, gbd_level = NA_integer_, whole_code = NA,
                                  note = "no level-3 GBD cause covers this code; kept as its own unit"))
  alt <- agg[icd10 == k & unit != p$unit, unit]
  data.table(icd10 = k, gbd_cause = p$unit, gbd_level = 3L, whole_code = p$whole,
             note = if (length(alt)) paste("other candidates:", paste(alt, collapse = "; ")) else "")
})
res <- rbindlist(rows)
res[, cause_id := cand$cause_id[match(paste(res$icd10, res$gbd_cause), paste(cand$icd10, cand$unit))]]
setcolorder(res, c("icd10", "gbd_cause", "gbd_level", "cause_id", "whole_code", "note"))
dir.create(dirname(out), showWarnings = FALSE)
fwrite(res, out, sep = "\t")
cat(sprintf("%s: %d indicators -> %d units\n", out, nrow(res), uniqueN(res$gbd_cause)))
print(res[, .(codes = paste(icd10, collapse = " "), n = .N), by = gbd_cause][order(-n)])
