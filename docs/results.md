# Results

*Draft of 11 November 2026. This draft covers the factor model, the factor GWAS and its
quality control, LD-score regression (including partitioned heritability and tissue
enrichment), loci, heterogeneity, Mendelian-randomisation instruments, clumping-and-thresholding
score weights and the polygenic-score evaluation (UK Biobank South Asian and Black participants;
FinnGen summary-level check; six comparator scores). Fine-mapping, conditional analysis, VEP
annotation and the eQTL/pQTL lookup have run and are being written up; colocalisation is
running. Still to run: MAGMA, TWAS, single-cell enrichment, gene prioritisation and the remaining
score methods (PRS-CS, LDpred2-auto, SBayesRC), after which the evaluation is repeated with those
weights. Methods are in [methods.md](methods.md).*

## Figures produced so far

Nature-format panels (PDF + 600-dpi PNG + source-data CSV) under `output/bothsex/figures/`:
- **Fig. 1** factor model: standardised loadings, SNP-heritability by loading, genetic
  correlation of each factor GWAS with its indicators against the model loading.
- **Fig. 2** factor GWAS: Manhattan plots for the three factors; ARD Miami plot with the
  *Q*SNP heterogeneity track.
- **Fig. 3** loci: counts by factor and heterogeneity class, *Q*SNP driver diagnoses, genome
  map of loci with nearest genes.
- **Fig. 4** enrichment: baselineLD binary annotations; GTEx tissue-specific expression.
- **Fig. 5** genes: credible-set sizes and lead PIPs by factor and heterogeneity class so far;
  the colocalisation, gene-prioritisation and MAGMA panels are drawn once those stages finish.
- **Fig. 6** polygenic-score evaluation: (a) IRR per SD of each factor score on its own diagnosis
  count in South Asian and Black participants; (b) mean ARD count by ARD-score decile; (c) IRR by
  C+T threshold; (d) comparator scores on the ARD count; (e) FinnGen summary-level liability *R*²
  per endpoint.
- **Extended data 1–5:** EFA scree, 43-trait genetic-correlation heatmap, QQ plots, *Q*SNP
  Manhattans, all baselineLD annotations.
- **Extended data 6:** PRS sample QC flow and within-group principal components.
- **Extended data 7:** ICD-10-code versus GBD-cause units, all-cause mortality, relatedness /
  covariate / hospital-contact sensitivities, the chapter-S injury negative control with and
  without the contact adjustment, and the likelihood-ratio statistic by C+T threshold.

## A general factor of age-related disease

Of 171 ICD-10 diagnoses with an effective sample size of at least 5,000, 43 had
significant SNP-heritability, and 34 of those loaded on a common genetic factor. A
bifactor model fitted acceptably both on held-out even chromosomes (CFI 0.938,
SRMR 0.086) and genome-wide (CFI 0.933, SRMR 0.073). It comprises a general
age-related-disease factor (ARD) on all 34 diagnoses plus three orthogonal domain
factors, and it fitted far better than a single common factor (Δχ² = 6,759, 15 df).

The general factor dominates the shared genetic variance: its explained common variance is
0.70. Its standardised loadings range from 0.31 to 0.87 (mean 0.61):
- **Highest:** absence of organs (Z90, 0.87), chest pain (R07, 0.85), drug-allergy status
  (Z88, 0.83), gastritis (K29, 0.83), history of medical treatment (Z92, 0.79) and COPD
  (J44, 0.79).
- **Lowest:** kidney stones (N20, 0.31) and atrial fibrillation (I48, 0.33).

The domain factors are as follows:
- **Cardiometabolic** is in effect a coronary-artery-disease factor. Its specific
  loadings are 0.78–0.85 for angina, myocardial infarction, chronic ischaemic heart disease
  and cardiac implants (I20, I21, I25, Z95), 0.59 for hyperlipidaemia (E78), and only 0.22
  and 0.15 for type 2 diabetes (E11) and obesity (E66).
- **Upper gastrointestinal** loads 0.43–0.77 on oesophageal and gastric diagnoses (K21,
  K22, K29, K44).
- **Osteoarthritis** retained only two indicators (M17, Z96) after hip osteoarthritis
  dropped out. Its loadings are not identified, so its GWAS is not interpreted and it is
  excluded from everything below.

Within their own indicators, the cardiometabolic and upper-gastrointestinal factors
account for 47% and 46% of the common genetic variance.

## Factor GWAS and quality control

The multivariate GWAS tested 7,982,480 variants. Its quality control is summarised below.
- **Failed fits:** 528 SNPs (0.007%) and were removed.
- **Gradient warnings:** 458,868 fits (5.7%) carried a gradient warning, traced to the I25
  residual variance on its bound; estimates were unchanged on refitting.
- **Smoothing:** no SNP's estimate changed by |*Z*smooth| > 1.96 under per-SNP smoothing.
- **Allele orientation:** all effect alleles mapped to the GRCh37 reference.

Removing *Q*SNP-heterogeneous SNPs and their LD partners shrank the general factor's
genome-wide significant SNPs by a third. The cardiometabolic factor lost 16%; the upper
gastrointestinal factor had a single *Q*SNP hit. Neither confounding nor stratification
contributes detectably: LD-score-regression intercepts were 0.99–1.00 with attenuation
ratios near zero. The general factor's λGC of 1.36 therefore reflects polygenic signal.

| | ARD | Cardiometabolic | Upper GI |
|---|---|---|---|
| *Q*SNP hits, *P* < 5 × 10⁻⁸ (df) | 19,706 (30) | 3,428 (8) | 1 (3) |
| LD partners removed | 97,524 | 20,452 | 99 |
| SNPs, Q-filtered tier | 7,864,726 | 7,958,072 | 7,981,852 |
| Genome-wide significant SNPs, full → Q-filtered | 2,089 → 1,367 | 3,203 → 2,694 | 1,063 → 1,063 |
| λGC, full tier | 1.36 | 1.14 | 1.13 |
| Mean χ², Q-filtered | 1.52 | 1.23 | 1.16 |
| LDSC intercept (SE), Q-filtered | 0.992 (0.009) | 0.988 (0.009) | 1.001 (0.008) |
| LDSC attenuation ratio | < 0 | < 0 | 0.01 |
| SNP-heritability *Z* | 28.6 | 12.1 | 13.2 |
| *N̂*, unit-variance scale | 16,845 | 7,255 | 4,782 |
| Marker indicator (loading λ) | Z90 (0.166) | I25 (0.296) | K22 (0.202) |
| *N̂*, marker scale | 610,290 | 82,771 | 117,330 |
| *Q*SNP λ, GenomicSEM df → effective df | 0.96 → 1.06 | 1.01 | 1.03 |

The general factor's *Q*SNP is miscalibrated at the 33 degrees of freedom GenomicSEM
assigns (λ 0.96, i.e. deflated). It is calibrated at the effective 30 df (λ 1.06), which
we used throughout.

## Genetic correlations

The three factor GWASs were close to orthogonal, as the model specifies. The genetic
correlations were −0.08 (SE 0.03) between ARD and cardiometabolic, −0.14 (0.04) between
ARD and upper GI, and 0.13 (0.05) between cardiometabolic and upper GI.

Each factor's genetic correlations with the diagnoses reproduced the measurement model.
Across indicators, the factor–indicator genetic correlation tracked the indicator's
standardised loading: r = 0.97 over the general factor's 34 indicators, 0.998 for
cardiometabolic and 0.999 for upper GI.

Heterogeneity itself is polygenic. LD-score regression on the *Q*SNP statistics gave a
heritability *Z* of 18.9 for the general factor (mean χ² 1.62) and 10.3 for
cardiometabolic, but only 3.1 for upper GI. Diagnosis-specific genetic effects are
therefore spread across the genome rather than confined to a few loci.

## Partitioned heritability and tissue enrichment

Stratified LDSC on the Q-filtered statistics placed each factor's heritability
disproportionately in evolutionarily conserved sequence.
- **Primate-conserved regions** hold 1.9% of SNPs but carried 20% of the general factor's
  SNP-heritability (10.2-fold enrichment, *q* = 2 × 10⁻⁶) and 35% of the upper-GI factor's
  (18-fold, *q* = 4 × 10⁻⁴). The cardiometabolic enrichment was smaller (8.8-fold,
  *q* = 0.04).
- **Active enhancers** separated the factors. For cardiometabolic, H3K27ac regions (39% of
  SNPs) carried 84% of the heritability (2.2-fold, *q* = 5 × 10⁻¹²) and super-enhancers
  2.6-fold. Super-enhancer enrichment was weaker for upper GI (2.1-fold) and ARD
  (1.4-fold).
- **Low-frequency variants:** for the general factor, the lowest-frequency common variants
  carried less heritability than expected (0.36–0.39-fold), a signature of negative
  selection.

Tissue-specific expression (LDSC-SEG, 54 GTEx tissues) separated the factors as well.
- **General factor:** all ten tissues significant at FDR < 5% were brain regions: frontal
  cortex, caudate, hippocampus, nucleus accumbens, cortex, anterior cingulate, amygdala,
  cerebellar hemisphere, putamen and cerebellum (*p* = 2 × 10⁻⁴ to 3 × 10⁻³). No other
  tissue reached significance. Together with its behavioural and neurodevelopmental loci,
  this suggests that genetic liability shared across age-related diseases acts partly
  through brain-expressed genes rather than through a single peripheral organ.
- **Cardiometabolic:** heritability concentrated in genes specific to the tibial, coronary
  and aortic arteries (*p* = 2–7 × 10⁻⁴), subcutaneous adipose tissue and ectocervix. This
  is the vascular pattern expected of coronary artery disease.
- **Upper GI:** no tissue was significant (top: liver, *p* = 0.008, *q* = 0.42).

## Genome-wide significant loci

| | ARD | Cardiometabolic | Upper GI |
|---|---|---|---|
| Independent SNPs | 54 | 96 | 19 |
| Loci | 49 | 52 | 11 |
| Heterogeneous (lead is a *Q*SNP hit or LD partner) | 10 | 8 | 0 |
| Not significant in any of the factor's indicators ("unique") | 28 | 13 | 5 |
| MHC | 1 | 0 | 1 |

Gene names below are the nearest protein-coding gene to the lead SNP. They locate loci
and are not causal assignments, which await fine-mapping and colocalisation.

- **Cardiometabolic** recovers the established coronary-artery-disease loci. These include
  LPA (*P* = 1 × 10⁻⁹⁷), 9p21/CDKN2B (2 × 10⁻⁹²), SORT1/PSRC1 (1 × 10⁻⁴⁰), APOE, LDLR,
  PHACTR1, TCF21, ADAMTS7, TRIB1, SH2B3, EDNRA, NOS3 and MIA3. They also include PCSK9,
  led by the missense variant rs11591147 (R46L); its minor T allele lowers the factor by
  0.57 SD.
  - Several of the strongest loci are heterogeneous: 9p21, SORT1, APOE, LDLR, PCSK9 and
    TRIB1. Their effects are larger on some cardiometabolic diagnoses than the factor
    predicts.
  - Thirteen loci were not significant in any of the factor's nine indicators, among them
    SMAD3, EHBP1, ZC3HC1, LIPC and an X-chromosome locus near FHL1.
- **Upper gastrointestinal** loci include FOXP1, BARX1, WT1 and CRLF1/TMEM59L, all
  previously associated with reflux disease and Barrett's oesophagus, plus an MHC locus.
  Five loci are unique, including MAPT (17q21.31).
- **General factor (ARD)** loci are mostly novel relative to its 34 indicators (28 of 49).
  - The strongest is near MLLT10 (*P* = 2 × 10⁻²⁰).
  - Others are near genes known for adiposity (TMEM18, SLC39A8) or for behavioural and
    neurodevelopmental traits (TTC12, FOXP2, NLGN1, CADM1, CELF4, TENM3).
  - Its heterogeneous loci are known single-disease signals that leak into the general
    factor. They are removed by the Q filter and kept out of factor-level conclusions:
    ABCG5/8 (gallstones), PITX2 (atrial fibrillation), PTPN22, IL2RA and the MHC
    (hypothyroidism), FTO (obesity), ARHGAP15 (diverticular disease) and 17q12 (asthma).

## Heterogeneity: which diagnoses depart from the factors

Clumping *Q*SNP hits gave 282 heterogeneity loci for the general factor, 81 for
cardiometabolic and 1 for upper GI. Each was attributed to the indicator with the largest
absolute *Z* at the lead SNP.

- **General factor:** heterogeneity is dominated by hypothyroidism (E03, 69 loci, including
  PTPN22, SH2B3 and FOXE1), atrial fibrillation (I48, 39; PITX2), type 2 diabetes (E11, 36;
  TCF7L2), asthma (J45, 31), diverticular disease (K57, 26), gallstones (K80, 20; ABCG8)
  and hyperlipidaemia (E78, 19).
- **Cardiometabolic:** heterogeneity comes from type 2 diabetes (48 loci; TCF7L2, CDKAL1,
  CCND2) and lipid loci (28; APOE, LDLR, SORT1, TRIB1, APOA5) whose effects on
  hyperlipidaemia exceed what the coronary-dominated factor implies. The FTO obesity locus
  is also heterogeneous.

## Mendelian-randomisation instruments

From the Q-filtered tier, after removing MHC variants and ambiguous palindromes and
clumping at *r*² < 0.001 within 10 Mb, we obtained 40 instruments for the general factor,
47 for cardiometabolic and 11 for upper GI. All have *F* > 29 (genome-wide significance).
The upper-GI count is low because 852 of its 1,063 significant SNPs lie in the MHC. The
instruments are released in the ardmr exposure format, with effects per SD of the factor.

## Polygenic-score weights

Clumping and thresholding of the Q-filtered tier retained about 369,000 approximately
independent SNPs per factor. Of these, 43 (ARD), 85 (cardiometabolic) and 30 (upper GI)
reach *P* < 5 × 10⁻⁸, and weights were defined at ten p-value thresholds.

Comparator scores were built with the same procedure from six published ageing-related
GWAS (mvAge, parental lifespan, longevity at the 90th and 99th percentiles, healthspan,
GrimAge acceleration and mvARD), each oriented so that a higher score means worse ageing.
Where the files carried no sample size, the *N̂* formula reproduced the published effective
sample sizes to within 0.3% (mvAge 1,961,102 vs 1,958,774; mvARD 476,662 vs 475,203).

## Polygenic-score evaluation

The scores were evaluated in UK Biobank participants outside the discovery sample: 8,764
self-reported South Asian and 7,033 Black participants after quality control and removal of
third-degree relatives (KING), of whom 6,144 and 4,478 had at least one ARD indicator diagnosis
and 551 and 409 had died. The two groups were analysed separately. The primary model was a
negative-binomial regression of the number of distinct ARD indicator diagnoses (ICD-10
three-character codes) on the standardised score, adjusted for age, age², sex and ten within-group
principal components. The clumping-and-thresholding *P* threshold for each group was the one with
the largest likelihood-ratio statistic in the *other* group, so no threshold was tuned on the
sample it is reported in (*P* < 10⁻³ for South Asian participants, *P* < 0.2 for Black
participants).

**Primary result.** One standard deviation of the ARD score was associated with 9.0% more ARD
diagnoses in South Asian participants (IRR 1.090, 95% CI 1.064–1.116, *P* = 2.4 × 10⁻¹²;
Cox–Snell pseudo-*R*² 0.56%) and 4.4% more in Black participants (IRR 1.044, 1.008–1.081,
*P* = 0.016; 0.08%). The association was present at every threshold from *P* < 10⁻⁴ upward in
South Asian participants (IRR 1.067–1.125) and smallest for the genome-wide-significant SNPs
alone (1.035, *P* = 0.006); in Black participants the IRRs were 0.99–1.05 across thresholds.
The Cardiometabolic score predicted the cardiometabolic diagnosis count in South Asian
participants (IRR 1.085, 1.052–1.118, *P* = 1.4 × 10⁻⁷) but not in Black participants (1.005,
0.968–1.043), and the Upper GI score predicted its count in neither group (1.027, 0.983–1.072;
1.012, 0.957–1.070).

**Sensitivity analyses.** Rolling the indicator codes up to 22 Global Burden of Disease causes
left the estimates unchanged (ARD score: IRR 1.087, 1.063–1.111 and 1.040, 1.007–1.074). The ARD
score was not associated with all-cause mortality over a median 13 years of follow-up (HR per SD
0.97, 0.89–1.05 in South Asian and 1.00, 0.89–1.12 in Black participants); the Cardiometabolic
score was, weakly, in South Asian participants (HR 1.095, 1.007–1.191, *P* = 0.035). Restricting
to the stricter GCTA-derived unrelated sets (IRR 1.094 and 1.023) or replacing the within-group
principal components with the global UK Biobank components (1.087 and 1.039) did not change the
primary estimates.

**Hospital contact and the negative control.** Any count derived from the inpatient record rises
with hospital contact. An initial negative control counting ICD-10 chapters S and T was itself
associated with the ARD score in South Asian participants (IRR up to 1.12), and this was driven by
chapter T (complications of medical care, T80–T88: IRR 1.16) and abolished by adjusting for the
number of distinct first-diagnosis dates (IRR 1.02, *P* = 0.5), whereas the ARD count remained
associated after the same adjustment (IRR 1.032, 1.017–1.047, *P* = 3.2 × 10⁻⁵ in South Asian
participants; 1.021, 0.999–1.044, *P* = 0.06 in Black participants). The negative control was
therefore defined as chapter-S injuries, which the ARD score did not predict in either group
(IRR 1.037, 0.965–1.115, *P* = 0.32; 1.053, 0.944–1.175, *P* = 0.36). The contact adjustment
over-corrects, because hospital contact is in part a consequence of the disease burden the score
predicts, so it bounds the association from below.

**Comparator scores.** On the same ARD count with the same models and threshold rule, the mvARD
score (Dinh 2026; a Genomic SEM factor of heart attack, high cholesterol, hypertension, stroke and
type 2 diabetes) was the strongest predictor in both groups (IRR 1.237, 1.206–1.269 in South Asian
and 1.215, 1.178–1.253 in Black participants) and the only comparator associated with mortality
(HR 1.26, 1.15–1.37 and 1.16, 1.05–1.29). The mvAge score (1.060, 1.034–1.086; 1.030, 1.000–1.060),
parental lifespan (1.048, 1.023–1.074; 1.020, 0.990–1.050), longevity (90th percentile 1.036 and
1.088; 99th 1.033 and 1.031) and healthspan (1.030; 1.017) were weaker than the ARD score in South
Asian participants, and GrimAge acceleration predicted nothing (1.012; 1.018). All comparator
discovery GWAS are European-ancestry (`extdata/README.md`), so none is evaluated in-sample. The
ARD factor therefore transfers to South Asian participants as well as or better than the published
ageing scores apart from mvARD, whose inputs are the largest cardiometabolic GWAS available and
whose diseases make up a large share of the indicator diagnoses.

**FinnGen summary-level check.** Against FinnGen R12 endpoints matching 27 of the 34 indicator
diagnoses, the ARD score explained a mean 0.3–0.5% of liability per endpoint (threshold
*P* < 10⁻³ and *P* < 0.2 respectively), in the expected direction for all 27 endpoints and highest
for obesity (0.8%), asthma (0.7%), type 2 diabetes (0.5%), knee osteoarthritis and intervertebral
disc disorders (0.5%). The Cardiometabolic score explained 0.8% on average over its six endpoints,
most for coronary atherosclerosis (1.5%), myocardial infarction (1.4%) and angina (1.3%). The
mvARD (0.4%) and mvAge (0.3%) scores were comparable to the ARD score on this European check.
These values are approximations (Finnish LD in the denominator, register-based endpoints) and are
reported per endpoint by necessity.

