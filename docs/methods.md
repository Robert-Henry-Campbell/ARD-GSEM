# Methods

## Contributing GWAS and trait selection

We modelled the joint genetic architecture of age-related disease using genome-wide
association summary statistics for ICD-10-coded diagnoses from the UK Biobank
European-ancestry cohort (Pan-UK Biobank). To ensure that every phenotype contributed
reliable genetic signal, we retained diagnoses with an effective sample size
(*N*eff = 4 / (1/*n*cases + 1/*n*controls)) of at least 5,000 and significant
common-variant heritability (stratified LD-score-regression *h*² *Z* ≥ 4). Because the
block jack-knife used to estimate the sampling covariance of the genetic covariance
matrix becomes unstable beyond ~43 traits, this heritability threshold also set the
practical ceiling on trait number. Forty-three diagnoses satisfied both criteria and
entered the model.

## Estimation of the genetic covariance matrix

Genetic variances and covariances were estimated by multivariable LD-score regression
using European LD scores computed from HapMap3 variants, yielding the genetic covariance
matrix (**S**) and its sampling covariance matrix (**V**). Case/control heritabilities and
co-heritabilities were placed on the liability scale using population prevalences derived
from the UK Biobank cohort. Unknown — and potentially complete — sample overlap among the
UK Biobank phenotypes was absorbed into the cross-trait LD-score-regression intercepts, and
hence into **V**, so that downstream standard errors remained calibrated.

## Split-sample discovery and confirmation of latent structure

To avoid confirming a factor structure on the same data used to derive it, we partitioned
the genome by chromosome parity and estimated two independent genetic covariance matrices,
one from odd and one from even chromosomes. Latent structure was **discovered** on the
odd-chromosome matrix by exploratory factor analysis and **confirmed** on the held-out
even-chromosome matrix by confirmatory factor analysis, both within the Genomic SEM
framework. The number of factors was chosen by convergent criteria — parallel analysis,
Velicer's minimum average partial, very-simple-structure, and the eigenvalue scree —
rather than any single rule. A diagnosis was assigned to a domain only when its dominant
standardised loading exceeded 0.40 and was at least 1.3-fold larger than its next-largest
loading; diagnoses whose within-factor genetic correlation lay near the genome-wide mean
were not forced into a factor.

## Model selection

On the held-out even-chromosome matrix, and on a common indicator set, we compared a single
common ("umbrella") factor, correlated domain factors, a higher-order factor, and a bifactor
model, evaluated jointly on the comparative fit index (CFI), the standardised root-mean-square
residual (SRMR), the Akaike information criterion, and nested χ²-difference tests. A bifactor
structure — a general age-related-disease factor loading on 34 diagnoses, together with three
orthogonal domain-specific factors (cardiometabolic, upper-gastrointestinal and osteoarthritis)
— provided the best and only adequate fit (held-out CFI = 0.938, SRMR = 0.086; full-sample
CFI = 0.933) with no anomalous specific-factor loadings, and significantly outperformed the
single-factor model (Δχ², *P* ≈ 0). The general factor loaded substantially across all four
domains (standardised loadings 0.27–0.86, mean 0.58), confirming a genuine pan-disease
dimension rather than one dominant organ system; nine diagnoses with negligible general
loadings (for example non-melanoma skin and breast cancer) were excluded. Parameters were
re-estimated on the full-sample covariance matrix for the genome-wide analysis.

## Multivariate GWAS of the latent factors

We performed a multivariate GWAS in which each SNP was regressed on every latent factor,
across ~7.98 million variants (1000 Genomes Phase 3 European reference, minor-allele
frequency > 0.005). Factor variances were fixed to unity, so SNP effects are expressed per
standard deviation of the latent factor. For every SNP and factor we additionally computed
the heterogeneity statistic *Q*SNP, which tests the null hypothesis that a SNP's effects
across a factor's indicator diagnoses are fully mediated by the common pathway; a large
*Q*SNP identifies a variant acting on individual diagnoses rather than on the shared factor.
Summary statistics for the general and domain factors were exported in GWAS-VCF format.

## Downstream use

The general age-related-disease factor provides a single, well-powered genetic index of
shared disease liability. For Mendelian-randomisation and polygenic analyses, instruments for
this factor are defined as variants reaching genome-wide significance for the factor and with
non-significant *Q*SNP — that is, variants acting through the shared liability rather than any
single diagnosis — with the domain-specific factors providing complementary, more granular
exposures.

## Key methods and software

Genomic SEM (Grotzinger et al., *Nat. Hum. Behav.* 2019) with lavaan for model specification;
LD-score regression (Bulik-Sullivan et al., *Nat. Genet.* 2015) for the genetic covariance and
sampling-covariance matrices; the *Q*SNP heterogeneity test as implemented in Genomic SEM's
multivariate GWAS. All summary statistics and LD scores were of European ancestry throughout.
