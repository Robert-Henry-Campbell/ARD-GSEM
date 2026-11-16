# Comparator GWAS summary statistics

External summary statistics used to build comparator polygenic scores (same methods as the
factor scores: C+T, PRS-CS, LDpred2-auto, SBayesRC) for the PRS evaluation stage (`prs_eval`).
Files uploaded by Robert Campbell on 2026-11-10 (mvAge, lifegen, healthspan, mvARD) or fetched
from the EBI GWAS Catalog FTP on the same day (Deelen 2019, GrimAge). The machine-readable
manifest the pipeline reads is `comparator_sumstats.tsv` (one row per score, with column
mappings, genome build and sample size).

| Use | Paper | Sumstats link | File |
|---|---|---|---|
| Multifactor | Rosoff 2023 (mvAge) | Zenodo 7926323 | mvAge.summary.EUR.txt (489 MB), plus README |
| Parental lifespan | Timmers 2019 | Edinburgh DataShare, doi:10.7488/ds/2463 | lifegen_phase2_bothpl_alldr_2017_09_18.tsv.gz (461 MB) |
| Participant longevity | Deelen 2019 | GWAS Catalog GCST008598 (90th percentile) and GCST008599 (99th percentile) | Also METAANALYSIS_90_metaX.csv.gz and METAANALYSIS_99_metaX.csv.gz, formatted for S-PrediXcan |
| Healthspan | Zenin 2019 | Zenodo 1302861 | healthspan_summary.csv.gz (237 MB); p-values given only as -log10(p) |
| Clock | GrimAge 1 acceleration (GWAS by McCartney 2021) | Edinburgh DataShare, doi:10.7488/ds/2834 | GrimAge_EUR_summary_statistics.txt (646 MB, European) or GrimAge_MRMEGA_summary_statistics.txt (trans-ancestry) |
| Comparator | Dinh 2026 (mvARD) | Zenodo 15522815 | Link taken from the paper's data availability statement; Zenodo blocked my requests, so I could not list the files |

## Files in this directory

| File | Source | Notes |
|---|---|---|
| `rosoff 2023-mvAge.summary.EUR.txt` | Zenodo 7926323 | GenomicSEM common-factor GWAS (mvAge); `SNP CHR BP MAF effect_allele other_allele beta se Pvalue`; no N column (N_hat derived from MAF and SE as for our factors) |
| `timmers 2019 DS_10283_3209.zip` → `lifegen_phase2_bothpl_alldr_2017_09_18.tsv.gz`, `lifegen_README.txt` | Edinburgh DataShare doi:10.7488/ds/2463 | Parental lifespan (LifeGen phase 2); unpacked here |
| `Results_90th_percentile.txt.gz`, `Results_99th_percentile.txt.gz` | EBI GWAS Catalog GCST008598 / GCST008599 | Deelen 2019 longevity, cases above the 90th / 99th survival percentile |
| `zenin-2019-healthspan_summary.csv.gz` | Zenodo 1302861 | Healthspan Cox GWAS; `SNPID,chr,pos,EA,RA,EAF,beta,se,Z,-log10(p-value)` |
| `GCST90014288_buildGRCh37.tsv.gz` (+ `GCST90014288_README.txt`) | EBI GWAS Catalog GCST90014288 | GrimAge acceleration, European (McCartney 2021); GRCh37 |
| `dihn-2026-mvARD_factor.tsv.gz` | Zenodo 15522815 | GenomicSEM userGWAS output for the mvARD factor (`est`, `se_c`, `Z_Estimate`, `Pval_Estimate`, Q); no N column (N_hat derived) |

## Discovery ancestry and overlap with the evaluation sample

The scores are evaluated in UK Biobank participants of self-reported South Asian and Black
ethnicity, so a comparator whose discovery GWAS included those participants would be evaluated
in-sample. None does: every comparator is a European-ancestry GWAS (UK Biobank contributions are
the white-British / European subsets).

| Score | Discovery data (from the papers) |
|---|---|
| mvAge (Rosoff 2023) | Genomic SEM common factor of European GWAS of healthspan (Zenin 2019, UKB), parental lifespan (Timmers 2019), longevity (Deelen 2019), frailty (Atkins 2021, UKB) and epigenetic-age acceleration (McCartney 2021); Neff 1,958,774 |
| Parental lifespan (Timmers 2019) | UK Biobank European participants + LifeGen consortium, European |
| Longevity 90th / 99th (Deelen 2019) | European-ancestry meta-analysis of longevity cohorts |
| Healthspan (Zenin 2019) | UK Biobank white-British participants |
| GrimAge acceleration (McCartney 2021) | European-ancestry meta-analysis (`GrimAge_EUR` file; the MR-MEGA trans-ancestry file was not used) |
| mvARD (Dinh 2026, *GeroScience*) | Genomic SEM factor over heart attack (Hartiala 2021, N 471,717), high cholesterol and hypertension (Dönertaş 2021, UK Biobank European, N 484,598), stroke (Malik 2018 MEGASTROKE European, N 446,696) and type 2 diabetes (Xue 2018, N 659,316); inputs selected as "data from European population"; Neff 475,203; no polygenic score was built or validated in the paper |

Citations: Rosoff DB et al. 2023; Timmers PRHJ et al. *eLife* 2019; Deelen J et al. *Nat Commun*
2019; Zenin A et al. *Commun Biol* 2019; McCartney DL et al. *Genome Biol* 2021; Dinh et al. 2026.
