# Methods

## Contributing GWAS and trait selection

We modelled the joint genetic architecture of age-related disease using genome-wide
association summary statistics for ICD-10-coded diagnoses from the UK Biobank
European-ancestry cohort (Pan-UK Biobank; SAIGE logistic mixed models). Of the
diagnoses available, 171 had an effective sample size
(*N*eff = 4 / (1/*n*cases + 1/*n*controls)) of at least 5,000. We then retained those with
significant common-variant heritability (published Pan-UKB SNP-heritability *Z* ≥ 4,
re-checked on our own LD-score-regression estimates). Because the block jack-knife used to
estimate the sampling covariance of the genetic covariance matrix becomes unstable beyond
~43 traits, this threshold also set the practical ceiling on trait number. Forty-three
diagnoses met both criteria and entered structure discovery. Thirty-four of them loaded at
least 0.30 on a single genetic factor and form the indicators of the final model (see
*Model selection*).

## Estimation of the genetic covariance matrix

Summary statistics were restricted to HapMap3 variants with minor-allele frequency ≥ 1% and
imputation INFO ≥ 0.8, with the MHC removed. Genetic variances and covariances were then
estimated by multivariable LD-score regression using European LD scores, giving the genetic
covariance matrix (**S**) and its sampling covariance matrix (**V**). Case/control
heritabilities and co-heritabilities were placed on the liability scale using the
effective-sample-size convention (sample prevalence 0.5) and population prevalences equal
to the UK Biobank case fractions. Unknown — and potentially complete — sample overlap among
the UK Biobank phenotypes was absorbed into the cross-trait LD-score-regression intercepts,
and hence into **V**, so that downstream standard errors remained calibrated.

**S** was not positive definite and was smoothed to the nearest positive-definite matrix
before modelling. The largest change was for D12–K63, whose genetic correlation moved from
0.997 to 0.829; 23 of the 561 pairs among the 34 final indicators changed by more than
0.05.

## Split-sample discovery and confirmation of latent structure

To avoid confirming a factor structure on the same data used to derive it, we partitioned
the genome by chromosome parity and estimated two independent genetic covariance matrices,
one from odd and one from even chromosomes. Latent structure was **discovered** on the
odd-chromosome matrix by exploratory factor analysis (geomin rotation) and **confirmed** on
the held-out even-chromosome matrix by confirmatory factor analysis, both within the Genomic
SEM framework.

The factor-number criteria did not agree: parallel analysis suggested 17 factors, the
eigenvalue > 1 rule 12, BIC 10, Velicer's minimum average partial 3, and very simple
structure 1. We therefore took the domain structure from the four-factor solution. A
diagnosis was assigned to a domain only when its dominant standardised loading was at
least 0.40 and at least 1.3 times its next-largest loading. Three of the four factors were
specified as domain-specific factors: cardiometabolic, upper gastrointestinal and
osteoarthritis. The fourth (D12, G47, K63, M51, N20, Z88, Z90) was not retained as a domain,
and its diagnoses load on the general factor only.

## Model selection

On the held-out even-chromosome matrix we compared a single common ("umbrella") factor
with a bifactor model on the same 34 indicators. The bifactor model had a general
age-related-disease (ARD) factor loading on all 34 diagnoses and three orthogonal
domain-specific factors. It fitted substantially better than the single factor
(Δχ² = 6,759 on 15 df; CFI 0.938 vs 0.686; SRMR 0.086 vs 0.123). This is the only
like-for-like comparison. Correlated-factor and higher-order models were fitted only to a
15-diagnosis subset of domain indicators (CFI 0.914, SRMR 0.115). With three domains the
higher-order model is equivalent to the correlated-factor model, so neither is directly
comparable with the bifactor model.

Re-estimated on the full-genome matrix, the bifactor model fitted acceptably (CFI 0.933,
SRMR 0.073; χ² = 2,832 on 512 df), meeting conventional thresholds for acceptable fit
(CFI ≥ 0.90, SRMR < 0.10) but not for good fit (CFI ≥ 0.95, SRMR < 0.05). Standardised
general-factor loadings ranged from 0.31 to 0.87 (mean 0.61), and all specific loadings
were positive (minimum 0.15). Nine of the 43 diagnoses had been excluded beforehand because
they loaded below 0.30 on a single factor fitted to the odd-chromosome half. They included
non-melanoma skin cancer, breast cancer and hip osteoarthritis (M16).

Losing M16 left the osteoarthritis factor with only two indicators (M17 and Z96). In an
orthogonal bifactor model the loadings of a two-indicator specific factor are not
separately identified: the M17 loading could take any value between about 0.21 and 0.28
with identical fit. The osteoarthritis factor's GWAS is therefore not interpretable. It is
reported as under-identified and excluded from all post-GWAS analyses. Re-estimating the
model without it would have required repeating the multivariate GWAS.

## Multivariate GWAS of the latent factors

We performed a multivariate GWAS (GenomicSEM `userGWAS`) in which each SNP was regressed on
every latent factor of the full-sample bifactor model. It covered 7,982,480 variants with
minor-allele frequency ≥ 0.01 in the 1000 Genomes Phase 3 European reference. SNP effects
were prepared with GenomicSEM `sumstats` from the SAIGE log-odds ratios (`se.logit = TRUE`).
Estimation used DWLS with the measurement model fixed (`fix_measurement`), standard genomic
control and per-SNP smoothing checks. Factors were scaled by unit-variance identification
(`std.lv = TRUE`), so SNP effects are expressed per standard deviation of the latent genetic
factor.

For every SNP and factor we also computed the heterogeneity statistic *Q*SNP. It tests the
null hypothesis that a SNP's effects on a factor's indicators are fully mediated by the
factor; a large *Q*SNP identifies a variant acting on individual diagnoses rather than
through the shared factor. GenomicSEM assigns *Q*SNP (number of indicators − 1) degrees of
freedom, giving 33 for the general factor. In the bifactor model, however, four SNP paths
(one per factor) are estimated against the general factor's 34 indicators. We therefore
recomputed the general factor's *Q*SNP p-values on 30 degrees of freedom. The 33-df
p-values are conservative: 15,669 SNPs reach *P* < 5 × 10⁻⁸ at 33 df and 19,706 at 30 df.
The domain factors' *Q*SNP was used as computed.

## Post-GWAS quality control

We applied the GenomicSEM developers' recommendations to each interpretable factor (ARD,
cardiometabolic, upper gastrointestinal).

1. **Failed fits removed.** We removed SNPs whose model returned no estimate, an error, or
   lavaan's "solution has not been found" warning.
2. **Benign warnings retained.** Fits that warned only that the gradient was not near zero
   were retained and flagged. This warning arises because the residual variance of I25 sits
   on its bound, and re-fitting such SNPs left their estimates unchanged.
3. **Smoothing check.** We removed SNPs whose per-SNP smoothing changed the estimate by
   |*Z*smooth| > 1.96.
4. **Allele orientation.** Effect alleles were oriented to the GRCh37 reference genome and
   annotated with Pan-UKB European allele frequencies.

These steps define the **full** tier.

The **Q-filtered** tier additionally removes every SNP with *Q*SNP *P* < 5 × 10⁻⁸ together
with its LD partners. Partners are SNPs at *r*² > 0.1 within 1 Mb in 1000 Genomes Phase 3
Europeans; for heterogeneous SNPs on chromosome X or absent from the panel, a ±500 kb window
was removed instead. Analyses of the factors as such — LD-score regression, polygenic scores
and Mendelian-randomisation instruments — use the Q-filtered tier. Locus discovery uses the
full tier, so that LD blocks are not left with holes, and every locus is labelled by its
heterogeneity status.

The effective sample size of each factor (*N̂*) was estimated with the GenomicSEM formula,
*N̂* = mean[1 / (2·MAF·(1 − MAF)·SE²)] over SNPs with MAF ≥ 0.10. Under unit-variance
identification *N̂* refers to a fully heritable latent factor and is small by construction.
We therefore also report *N̂* rescaled to unit-loading identification on each factor's
highest-loading indicator (the marker). With λ the marker's unstandardised loading, effects
and standard errors are multiplied by λ and *N̂* is divided by λ². Z-statistics and
p-values are unaffected.

## Independent loci

Independent genome-wide significant SNPs (*P* < 5 × 10⁻⁸) were identified by clumping the
full tier in 1000 Genomes Phase 3 Europeans (plink 1.9; *r*² < 0.1 within 500 kb, 1 Mb on
chromosomes 6 and 8). Significant SNPs absent from the panel, and those on chromosome X,
were clumped by distance (±500 kb). Loci were defined as each independent SNP ±500 kb,
merged where overlapping; any locus touching the MHC (chr6:25–35 Mb) was treated as one
locus.

Each locus received two labels:
- **Heterogeneous:** its lead SNP was a *Q*SNP hit or an LD partner of one.
- **Unique:** none of its independent SNPs was, or was in LD (*r*² ≥ 0.1, within 1 Mb) with,
  a genome-wide significant SNP in any of that factor's indicator GWAS.

Clumps of *Q*SNP hits were tabulated separately, together with the indicator driving each.

## Summary-statistic release and Mendelian-randomisation instruments

Factor summary statistics were released in GWAS-VCF format: one multi-sample file and one
file per factor (GRCh37). Each record carries:
- the effect of the alternative allele and its standard error, −log₁₀ *P* and allele
  frequency;
- the effective sample size *N̂*;
- the *Q*SNP statistic and its p-value;
- flags for Q-filter failure and smoothing.

Instruments for Mendelian randomisation were taken from the Q-filtered tier using the
ardmr package's defaults. The criteria were genome-wide significance, autosomal, outside
the MHC and *F* ≥ 10. Ambiguous palindromic SNPs (A/T or C/G with allele frequency 0.42–0.58)
were dropped. The remainder were clumped at *r*² < 0.001 within 10 Mb in 1000 Genomes
Phase 3 Europeans. Instruments were written in ardmr's exposure format, with effects per
standard deviation of the factor and sample size *N̂*.

## LD-score regression of the factors

Each factor's Q-filtered and full summary statistics, and its *Q*SNP statistics, were
aligned to HapMap3 alleles and analysed by LD-score regression (GenomicSEM `ldsc`, European
LD scores, *N* = *N̂*). The LDSC intercept, attenuation ratio and heritability *Z*-statistic
do not depend on the sample size supplied. The heritability point estimate does depend on
it and is not interpretable on the unit-variance scale, so we report the former three.
Genetic correlations were estimated from the Q-filtered statistics between factors and
between each factor and each of the 34 indicators.

Partitioned heritability was estimated by stratified LDSC with the baselineLD v2.2 model
(1000 Genomes Phase 3 European LD scores and frequencies, HapMap3 regression weights
excluding the MHC). Heritability enrichment was assessed for the 83 binary functional
annotations (all except the all-SNP base annotation).
For the 13 continuous annotations, where enrichment is not defined, we used the
per-annotation coefficient τ. Enrichment and τ p-values were each corrected by
Benjamini–Hochberg.

Tissue enrichment followed the LDSC-SEG approach. Because the published LDSC-SEG tissue
LD scores could not be downloaded, we built GTEx v8 annotations ourselves. For each of 54
tissues we took the 10% of protein-coding genes with the highest expression specificity
(median TPM in that tissue divided by the sum over tissues), extended each gene body by
±100 kb (GENCODE v19), and computed LD scores in 1000 Genomes Phase 3 Europeans (1 cM
window). Each tissue annotation was tested jointly with an all-genes control annotation and
the baseline v1.2 model, with Benjamini–Hochberg correction over tissues.

## Polygenic scores

Score weights were computed from the Q-filtered tier by four methods: clumping and
thresholding (plink 1.9; autosomal SNPs present in 1000 Genomes Phase 3 Europeans clumped at
*r*² < 0.1 within 250 kb, ten p-value thresholds from 5 × 10⁻⁸ to 1), PRS-CS-auto, LDpred2-auto
(no tuning sample) and SBayesRC, each with the factor's effective sample size *N̂*.

Comparator scores were built with the same methods from published summary statistics of six
ageing-related traits (`extdata/comparator_sumstats.tsv`): the multivariate ageing factor mvAge
(Rosoff et al. 2023), parental lifespan (Timmers et al. 2019), longevity above the 90th and
99th survival percentiles (Deelen et al. 2019), healthspan (Zenin et al. 2019), GrimAge
acceleration (McCartney et al. 2021) and the multivariate age-related-disease factor mvARD
(Dinh et al. 2026). Every comparator was oriented so that a higher score means worse ageing or
more disease; the orientation of each file was verified at known risk alleles (LPA
rs10455872, 9p21 rs1333049, CHRNA3 rs1051730). Effective sample sizes were taken from the
files where given, otherwise estimated from MAF and standard errors as for the factors.

### Evaluation

Scores were evaluated in UK Biobank participants outside the discovery sample: self-reported
South Asian (field 21000 codes 3, 3001–3004; n = 9,414) and Black (4, 4001–4003; n = 7,641)
participants, who are not in the Pan-UKB European analysis, analysed separately. Participants
with missing genetic quality-control data, a mismatch between reported and genetic sex, or sex
chromosome aneuploidy were excluded; one of each pair related to the third degree was removed
(plink2 KING cut-off 0.0442 on LD-pruned array genotypes with long-range LD regions excluded),
and ten principal components were computed within each group. Scores were computed with
plink2 from imputed dosages (7.8 million variants, alleles aligned to the genotype reference)
and standardised within group.

The outcome was the number of distinct ARD indicator diagnoses (ICD-10 three-character codes
in the hospital inpatient record, field 41270), modelled by negative-binomial regression on the
score with age, age², sex and the ten within-group principal components; the effect is reported
as the incidence-rate ratio per standard deviation of the score with a likelihood-ratio test
against the covariates-only model. The Cardiometabolic and Upper GI scores were evaluated on
the counts of their own indicators. Two sensitivity analyses were pre-specified: the same model
after rolling the indicator codes up to Global Burden of Disease level-3 causes (34 codes → 22
units; `manifest/indicator_gbd_causes.tsv`), and all-cause mortality by Cox regression with age
as the timescale. Robustness checks repeated the primary model on the stricter GCTA-derived
unrelated subsets, with the global UK Biobank principal components, with an additional adjustment
for hospital contact (log of 1 + the number of distinct first-diagnosis dates in the inpatient
record), and on the count of injury diagnoses (ICD-10 chapter S) as a negative control. Chapter T
was not used for the control because complications of medical care (T80–T88) and poisonings and
adverse effects (T36–T65) are coded during the hospital care that disease burden itself generates,
so they track hospital contact rather than injury. Clumping-and-thresholding thresholds were chosen
for each weight set and group by the likelihood-ratio statistic in the *other* group. Comparator
scores were evaluated with the same models. An inverse-variance meta-analysis of the two groups
is given for completeness; the groups are reported separately throughout.

As a European-ancestry check without individual-level data, each score was also evaluated
against FinnGen R12 summary statistics of the endpoints matching the indicator diagnoses
(`manifest/indicator_finngen_endpoints.tsv`): the score–endpoint covariance was estimated from
the FinnGen per-allele effects and allele frequencies, the score variance from the 1000 Genomes
European reference, and the resulting R² converted to the liability scale, with
leave-one-chromosome-out jackknife intervals (the pseudo-validation approach of lassosum and
PUMAS).

## Key methods and software

- Genomic SEM 0.0.5 (Grotzinger et al., *Nat. Hum. Behav.* 2019), with lavaan 0.6-21 for
  model specification.
- LD-score regression (Bulik-Sullivan et al., *Nat. Genet.* 2015) for the genetic
  covariance and sampling-covariance matrices and for factor LDSC. Stratified LDSC and
  LDSC-SEG (Finucane et al., *Nat. Genet.* 2015, 2018) used a Python 3 port of ldsc whose
  LD-score computation we validated against the published baseline LD scores.
- plink 1.9 for LD and clumping.
- The *Q*SNP heterogeneity test as implemented in Genomic SEM's multivariate GWAS.

All summary statistics and LD references were of European ancestry throughout.
