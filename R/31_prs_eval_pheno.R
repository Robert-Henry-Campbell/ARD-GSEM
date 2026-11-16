# Phenotypes for the PRS evaluation (pure functions over the UKB application extract):
#  - pe_parse_icd: pipe-joined HES ICD-10 codes (41270) and first-diagnosis dates (41280.0.all),
#    aligned position by position, to one row per (participant, code)
#  - pe_build_pheno: per participant, the number of distinct ARD indicator ICD-10 three-character
#    codes ever recorded (count_<factor>; the primary unit), the same count after rolling codes up
#    to GBD causes (count_<factor>_gbd; sensitivity unit, manifest/indicator_gbd_causes.tsv), the
#    injury-chapter count (negative control), and all-cause mortality with age-scale follow-up

UKB_SENTINEL_DATES <- c("1900-01-01", "1901-01-01", "1902-02-02", "1903-03-03", "2037-07-07")

pe_parse_icd <- function(eid, codes, dates) {
  codes[is.na(codes)] <- ""; dates[is.na(dates)] <- ""
  cl <- strsplit(codes, "|", fixed = TRUE); dl <- strsplit(dates, "|", fixed = TRUE)
  nc <- lengths(cl); nd <- lengths(dl)
  bad <- nc > 0 & nc != nd
  if (any(bad)) {
    log_warn("prs_eval", sprintf("%d participants have unequal ICD-10 code and date counts; their dates are set to NA", sum(bad)))
    dl[bad] <- lapply(nc[bad], function(n) rep("", n))
  }
  out <- data.table(eid = rep(as.character(eid), nc), code = unlist(cl), date_chr = unlist(dl))
  out <- out[nzchar(code)]
  out[, date := as.IDate(as.Date(date_chr, format = "%Y-%m-%d"))]
  out[date_chr %in% UKB_SENTINEL_DATES, date := NA]
  out[, date_chr := NULL]
  out[, code3 := substr(code, 1, 3)]
  out[]
}

# Columns of the application extract needed downstream, restricted to `ids` (character eids).
pe_read_pheno_extract <- function(path, ids) {
  hdr <- names(fread(path, nrows = 0))
  pick <- function(f) { v <- hdr[startsWith(hdr, paste0(f, "."))]; if (length(v)) v[1] else NA_character_ }
  want <- c(eth = pick("21000"), sex = pick("31"), genetic_sex = pick("22001"), age = pick("21022"),
            assess = pick("53"), death_date = pick("40000"), death_cause = pick("40001"), icd = pick("41270"),
            icd_dates = pick("41280"), kinship = pick("22021"), aneuploidy = pick("22019"), het_miss = pick("22027"),
            townsend = pick("189"), bmi = pick("21001"))
  if (any(is.na(want[c("eth", "sex", "age", "assess", "death_date", "icd", "icd_dates")])))
    log_fatal("prs_eval", "phenotype extract lacks one of 21000/31/21022/53/40000/41270/41280")
  want <- want[!is.na(want)]
  pcs <- hdr[grepl("^22009\\.0\\.([1-9]|10)$", hdr)]
  d <- fread(path, select = c("eid", unname(want), pcs), colClasses = "character", na.strings = c(""))
  setnames(d, unname(want), names(want))
  if (length(pcs)) setnames(d, pcs, paste0("gpc", sub("^22009\\.0\\.", "", pcs)))
  d[eid %in% as.character(ids)]
}

# Self-reported ethnicity codes (field 21000, instance 0) -> group name or NA.
pe_assign_group <- function(code21000, groups) {
  v <- suppressWarnings(as.integer(round(as.numeric(code21000))))
  out <- rep(NA_character_, length(v))
  for (g in names(groups)) out[v %in% as.integer(unlist(groups[[g]]$codes))] <- g
  out
}

pe_read_gbd_map <- function(path) {
  m <- fread(path, sep = "\t", colClasses = "character")
  if (!all(c("icd10", "gbd_cause") %in% names(m))) log_fatal("prs_eval", sprintf("%s lacks icd10/gbd_cause columns", path))
  m[, .(icd10, gbd_cause)]
}

# ph: pe_read_pheno_extract output; indicators: named list factor -> ICD-10 three-character codes;
# gbd_map: pe_read_gbd_map output (NULL skips the roll-up); censor: list(hes, death) as "YYYY-MM-DD".
# The negative control counts chapter-S injuries only: chapter T (complications of medical care
# T80-T88, poisoning and adverse effects T36-T65) is coded during the hospital care that disease
# burden itself generates, so it tracks hospital contact rather than injury. n_hes_dates (distinct
# first-diagnosis dates) is the hospital-contact proxy used by the "contact" covariate set.
pe_build_pheno <- function(ph, indicators, gbd_map = NULL, censor = list(hes = "2022-10-31", death = "2022-12-19"),
                           injury_chapters = "S") {
  long <- pe_parse_icd(ph$eid, ph$icd, ph$icd_dates)
  out <- data.table(eid = as.character(ph$eid),
                    age = suppressWarnings(as.numeric(ph$age)),
                    sex = suppressWarnings(as.integer(round(as.numeric(ph$sex)))),
                    assess = as.IDate(as.Date(ph$assess, format = "%Y-%m-%d")),
                    death_date = as.IDate(as.Date(ph$death_date, format = "%Y-%m-%d")))
  out[, age2 := age^2]
  dcens <- as.IDate(censor$death)
  out[, death := as.integer(!is.na(death_date) & death_date <= dcens)]
  out[, fu_end := fifelse(death == 1L, death_date, dcens)]
  out[, t_death := as.numeric(fu_end - assess) / 365.25]
  out[, age_entry := age]
  out[, age_exit := age + pmax(t_death, 0.5 / 365.25)]
  add_count <- function(col, tab) { set(out, j = col, value = 0L); i <- match(tab$eid, out$eid); set(out, i = i[!is.na(i)], j = col, value = as.integer(tab$N[!is.na(i)])) }
  for (f in names(indicators)) {
    codes <- indicators[[f]]
    add_count(paste0("count_", f), long[code3 %in% codes, .(N = uniqueN(code3)), by = eid])
    if (!is.null(gbd_map)) {
      lg <- merge(long[code3 %in% codes, .(eid, code3)], gbd_map, by.x = "code3", by.y = "icd10")
      add_count(paste0("count_", f, "_gbd"), lg[, .(N = uniqueN(gbd_cause)), by = eid])
    }
  }
  add_count("injury_count", long[substr(code3, 1, 1) %in% injury_chapters, .(N = uniqueN(code3)), by = eid])
  add_count("n_hes_dates", long[!is.na(date), .(N = uniqueN(date)), by = eid])
  out[, any_hes := eid %in% long$eid]
  out[]
}
