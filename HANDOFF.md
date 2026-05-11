# GSEM Pipeline — Session Handoff

## What This Project Is

Sex-stratified Genomic SEM pipeline for Neale UKBB round-2 ICD-10 ARD GWAS summary statistics. Runs LDSC, EFA, CFA, and sex-comparison independently per sex on traits in `sumstats/male/` and `sumstats/female/`.

Root directory: `/mnt/sdg/robert/ardmr/GSEM/`

## What Has Been Done

### Completed

1. **Git repo initialized** (`git init`, `.gitignore` in place, no commits yet)
2. **renv environment set up** with all packages installed and locked:
   - GenomicSEM (from GitHub), lavaan, data.table, psych, GPArotation, furrr, future, ggplot2, rmarkdown, jsonlite, digest, testthat, and dependencies
   - `renv.lock` is populated and valid
   - Note: `semPlot` was dropped because its dependency `OpenMx` fails to compile on this server (Eigen template errors). Not critical — it's only for path diagrams.
3. **Full pipeline code written** — 8 R modules:
   - `R/utils.R` — logging system (dual console+file, 5 levels, JSON structured log), Neff calculation, variant parsing, rsID mapping, filtering, a priori model builder, EFA-to-lavaan converter, S matrix PSD check, loading diff formula, stage manifest system for `--resume`
   - `R/00_setup.R` — config loading, path validation, reference file checks, directory creation
   - `R/01_munge.R` — reads `.tsv.bgz` sumstats, joins to variants manifest for rsID mapping, computes Neff, filters, calls `GenomicSEM::munge()`
   - `R/02_ldsc.R` — runs `GenomicSEM::ldsc()`, saves S/V matrices, h² QC table, intercept matrix, filters traits by h²_z threshold
   - `R/03_efa.R` — parallel analysis for factor count, `psych::fa()` with geomin rotation, saves loadings
   - `R/04_cfa.R` — fits both EFA-derived and ICD-10-chapter a priori models via `GenomicSEM::usermodel()`
   - `R/05_comparison.R` — sex-difference loading comparison (Wald Z-test on unstandardized loadings), factor alignment check, Bonferroni + FDR correction
   - `R/06_report.R` — compiles all outputs, renders HTML report via rmarkdown
4. **Entry point**: `scripts/run_pipeline.R` with CLI args (`--sex`, `--mode`, `--stage`, `--threads`, `--resume`)
5. **Report template**: `templates/report.Rmd` (HTML output with h², EFA, CFA, comparison tables)
6. **111 unit tests written and all passing**:
   - `tests/testthat/test-neff.R` (6 tests)
   - `tests/testthat/test-munge-parsing.R` (8 tests)
   - `tests/testthat/test-munge-filtering.R` (7 tests)
   - `tests/testthat/test-ldsc-inputs.R` (7 tests)
   - `tests/testthat/test-efa.R` (5 tests)
   - `tests/testthat/test-cfa-model-building.R` (5 tests)
   - `tests/testthat/test-comparison.R` (11 tests)
   - `tests/testthat/test-config.R` (4 tests)
   - `tests/testthat/test-logging.R` (6 tests)
   - `tests/testthat/test-resume.R` (4 tests)
   - Run with: `Rscript -e 'testthat::test_dir("tests/testthat")'`
7. **Integration smoke test script**: `tests/test_smoke.R` (not yet runnable — needs reference data)
8. **Config**: `config/pipeline.yaml` with all paths, thresholds, parallelism settings

### NOT Done — Blocked

**The pipeline cannot run yet because the LD score reference data is missing.**

The `reference/` directory needs:
- `reference/eur_w_ld_chr/` — directory with 22 chromosome LD score files (`1.l2.ldscore.gz` through `22.l2.ldscore.gz`, plus `.l2.M`, `.l2.M_5_50` files per chromosome)
- `reference/w_hm3.snplist` — HapMap3 SNP list (~1.2M SNPs, columns: SNP, A1, A2)

All traditional download URLs are broken as of 2026-05-11:
- `https://data.broadinstitute.org/alkesgroup/LDSCORE/eur_w_ld_chr.tar.bz2` → 301 → 404
- `https://storage.googleapis.com/broad-alkesgroup-public/LDSCORE/eur_w_ld_chr.tar.bz2` → 404
- The Broad has migrated to a GCS bucket (`broad-alkesgroup-public`) but the old tar.bz2 archives appear to have been removed
- Zenodo record 7768714 → 504 timeout
- UT Austin Box link from GenomicSEM wiki → 404

**The user needs to obtain these files and place them at the paths above.** Possible sources:
- Check the [GenomicSEM wiki](https://github.com/GenomicSEM/GenomicSEM/wiki) for updated links
- Check [GenomicSEM GitHub issues](https://github.com/GenomicSEM/GenomicSEM/issues) — others may have reported the broken links
- Download from a colleague or institutional mirror
- The `eur_w_ld_chr.tar.bz2` archive is ~50 MB; `w_hm3.snplist.bz2` is ~2 MB

There is a junk 214-byte file that was already cleaned up from `reference/`.

## What To Do Next

Once reference data is in place:

1. **Verify reference data**:
   ```bash
   ls reference/eur_w_ld_chr/ | grep -c 'l2.ldscore.gz'  # should be 22
   wc -l reference/w_hm3.snplist  # should be ~1,217,312
   ```

2. **Run the smoke test**:
   ```bash
   cd /mnt/sdg/robert/ardmr/GSEM
   Rscript scripts/run_pipeline.R --sex both --mode smoke --threads 24
   ```

3. **Check outputs**:
   - Logs: `output/logs/latest.log`
   - h² tables: `output/{male,female}/ldsc/h2_qc.csv`
   - CFA fit: `output/{male,female}/cfa/fit_comparison.csv`
   - Sex comparison: `output/comparison/loading_diff.csv`
   - Report: `output/report.html`

4. **If issues arise**, run stages individually:
   ```bash
   Rscript scripts/run_pipeline.R --sex male --stage munge --threads 24
   Rscript scripts/run_pipeline.R --sex male --stage ldsc
   ```

5. **Make initial git commit** once smoke test passes.

## Key Design Decisions

- **Neff**: `4 * n_cases * n_controls / (n_cases + n_controls)` for all binary traits (Neale linear regression on 0/1)
- **A1/A2**: In `variant = chr:pos:ref:alt`, ALT = effect allele → A1 = alt, A2 = ref
- **rsID mapping**: Joins sumstats `variant` column to `manifest/variants.tsv.bgz` (13.8M rows) on `variant`, extracts `rsid`
- **h² threshold**: Z ≥ 2 (relaxed for sex-stratified halved sample sizes)
- **Observed-scale h²**: No population prevalence available → no liability-scale conversion
- **LDSC free intercept**: Accounts for within-sex sample overlap (all traits from same UKB cohort)
- **CFA**: Both EFA-derived AND a priori (ICD-10 chapter grouping) models fitted
- **Sex comparison**: Independent per-sex fits + Wald Z-test on unstandardized loading differences + Bonferroni/FDR correction
- **No semPlot**: Dropped due to OpenMx compile failure; path diagrams not available
- **All operations stay within `/mnt/sdg/robert/ardmr/GSEM/`** — never read or write outside this directory

## File Map

```
GSEM/
├── config/pipeline.yaml          ← single source of truth for all settings
├── R/
│   ├── utils.R                   ← logging, helpers, all pure functions (tested)
│   ├── 00_setup.R                ← config loading, validation, directory setup
│   ├── 01_munge.R                ← sumstats → .sumstats.gz (rsID mapped, Neff)
│   ├── 02_ldsc.R                 ← LDSC → S matrix, V matrix, h² QC
│   ├── 03_efa.R                  ← parallel analysis → factor loadings
│   ├── 04_cfa.R                  ← GenomicSEM::usermodel() for EFA + a priori models
│   ├── 05_comparison.R           ← sex-difference Wald tests + alignment check
│   └── 06_report.R               ← compile results → HTML report
├── scripts/
│   ├── run_pipeline.R            ← CLI entry point
│   ├── install_packages.R        ← renv install (already run, all packages installed)
│   └── download_reference.sh     ← BROKEN URLs — needs updated links
├── tests/
│   ├── testthat/                 ← 111 unit tests, ALL PASSING
│   └── test_smoke.R              ← end-to-end integration test
├── templates/report.Rmd          ← rmarkdown report template
├── manifest/                     ← neale manifests (.rda) + variants.tsv.bgz
├── meta/                         ← icd10_categories.csv, ARD .rda files
├── sumstats/{male,female}/       ← .tsv.bgz GWAS files (4-5 traits for smoke test)
├── reference/                    ← eur_w_ld_chr/ (22 chr) + w_hm3.snplist
├── renv.lock                     ← locked package versions
└── .Rprofile                     ← renv activation
```

## Server Info

- Host: doraemon6, Ubuntu, kernel 5.15.0-176
- 64 cores, 754 GB RAM, 3.6 TB free on /mnt/sdg
- R 4.4.2, renv 1.2.2
- Smoke test traits: D50, D64, E10, E11 (both sexes) + E14 (male only)
