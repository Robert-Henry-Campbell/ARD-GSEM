#!/usr/bin/env bash
# Downloads the public reference data used by the post-GWAS stages into reference/postgwas/
# (gitignored). Data only — no code is executed. Idempotent and resumable: finished files are
# skipped, partial ones resume from *.part, archives are unpacked once (marked by .unpacked).
#
# Usage: bash scripts/download_postgwas_reference.sh [group ...]
#   groups: gene magma gtex ldsc fusion prscs sbayesrc eqtl vep scdrs   (default: all)
#           ldsc_seg is optional (currently 404; see group_ldsc_seg)
#   JOBS=4 sets how many files download in parallel within a group.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
REF="${REF:-$ROOT/reference/postgwas}"
JOBS="${JOBS:-4}"
mkdir -p "$REF"

log() { printf '%s [download] %s\n' "$(date '+%F %T')" "$*"; }

# fetch URL DEST: resumable download to DEST (skips if DEST exists and is non-empty)
fetch() {
  local url="$1" dest="$2"
  [[ -s "$dest" ]] && return 0
  mkdir -p "$(dirname "$dest")"
  log "GET $(basename "$dest")"
  curl -sS -L --fail --retry 8 --retry-delay 20 --retry-all-errors -C - -o "$dest.part" "$url"
  mv "$dest.part" "$dest"
}
export -f fetch log

# unpack ARCHIVE DIR: extract once (tar.gz/tgz/tar.bz2/tar/zip), leaving a marker
unpack() {
  local f="$1" dir="$2" mark="$1.unpacked"
  [[ -e "$mark" ]] && return 0
  mkdir -p "$dir"
  log "unpack $(basename "$f")"
  case "$f" in
    *.zip) unzip -q -o "$f" -d "$dir" ;;
    *.tar.bz2) tar -xjf "$f" -C "$dir" ;;
    *.tar.gz|*.tgz) tar -xzf "$f" -C "$dir" ;;
    *.tar) tar -xf "$f" -C "$dir" ;;
    *) log "don't know how to unpack $f"; return 1 ;;
  esac
  touch "$mark"
}

# fetch_many DIR URL...: parallel fetches into DIR, file named after the URL basename
fetch_many() {
  local dir="$1"; shift
  printf '%s\n' "$@" | xargs -P "$JOBS" -I{} bash -c 'u="{}"; b="${u%%\?*}"; fetch "$u" "'"$dir"'/$(basename "$b")"'
}

group_gene() {
  fetch https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz \
        "$REF/gene/gencode.v19.annotation.gtf.gz"
  fetch https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz \
        "$REF/chain/hg19ToHg38.over.chain.gz"
}

group_magma() {
  local d="$REF/magma"
  fetch https://vu.data.surfsara.nl/index.php/s/Pj2orwuF2JYyKxq/download "$d/NCBI37.3.zip"
  fetch https://vu.data.surfsara.nl/index.php/s/VZNByNwpD8qqINe/download "$d/g1000_eur.zip"
  fetch https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2024.1.Hs/msigdb.v2024.1.Hs.entrez.gmt \
        "$d/msigdb.v2024.1.Hs.entrez.gmt"
  unpack "$d/NCBI37.3.zip" "$d"
  unpack "$d/g1000_eur.zip" "$d"
}

group_gtex() {
  fetch https://storage.googleapis.com/adult-gtex/bulk-gex/v8/rna-seq/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz \
        "$REF/gtex/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz"
}

group_ldsc() {
  local d="$REF/ldsc" z=https://zenodo.org/records/10515792/files
  fetch_many "$d" \
    "$z/1000G_Phase3_baselineLD_v2.2_ldscores.tgz?download=1" \
    "$z/1000G_Phase3_baseline_v1.2_ldscores.tgz?download=1" \
    "$z/1000G_Phase3_weights_hm3_no_MHC.tgz?download=1" \
    "$z/1000G_Phase3_frq.tgz?download=1" \
    "$z/1000G_Phase3_plinkfiles.tgz?download=1" \
    "$z/hm3_no_MHC.list.txt?download=1"
  for f in "$d"/*.tgz; do unpack "$f" "$d"; done
}

# Optional, not in the default list: alkesgroup.broadinstitute.org returns 404 for these
# (the official LDSC-SEG scores now sit in a requester-pays bucket). The sldsc stage builds
# GTEx v8 / Tabula Sapiens specific-expression annotations instead.
group_ldsc_seg() {
  local d="$REF/ldsc_seg" b=https://alkesgroup.broadinstitute.org/LDSCORE/LDSC_SEG_ldscores
  fetch "$b/Multi_tissue_gene_expr.ldcts" "$d/Multi_tissue_gene_expr.ldcts"
  fetch "$b/Multi_tissue_gene_expr_1000Gv3_ldscores.tgz" "$d/Multi_tissue_gene_expr_1000Gv3_ldscores.tgz"
  gzip -t "$d/Multi_tissue_gene_expr_1000Gv3_ldscores.tgz" # the server sends no length; check it
  unpack "$d/Multi_tissue_gene_expr_1000Gv3_ldscores.tgz" "$d"
}

group_fusion() {
  local d="$REF/fusion"
  fetch https://alkesgroup.broadinstitute.org/FUSION/LDREF.tar.bz2 "$d/LDREF.tar.bz2"
  unpack "$d/LDREF.tar.bz2" "$d"
  mapfile -t urls < <(curl -sS -m 60 -L http://gusevlab.org/projects/fusion/ |
    grep -oE 'https://s3\.us-west-1\.amazonaws\.com/gtex\.v8\.fusion/ALL/GTExv8\.ALL\.[A-Za-z0-9_]+\.tar\.gz' | sort -u)
  log "FUSION GTEx v8 weight sets: ${#urls[@]}"
  fetch_many "$d/weights" "${urls[@]}"
  for f in "$d"/weights/*.tar.gz; do unpack "$f" "$d/weights"; done
}

group_prscs() {
  local d="$REF/prscs"
  fetch "https://www.dropbox.com/s/t9opx2ty6ucrpib/ldblk_ukbb_eur.tar.gz?dl=1" "$d/ldblk_ukbb_eur.tar.gz"
  unpack "$d/ldblk_ukbb_eur.tar.gz" "$d"
}

group_sbayesrc() {
  local d="$REF/sbayesrc" b=https://gctbhub.cloud.edu.au/data/SBayesRC/resources/v2.0
  fetch "$b/Annotation/annot_baseline2.2.zip" "$d/annot_baseline2.2.zip"
  fetch "$b/LD/Imputed/ukbEUR_Imputed.zip" "$d/ukbEUR_Imputed.zip"   # 57.9 GB
  unpack "$d/annot_baseline2.2.zip" "$d"
  unpack "$d/ukbEUR_Imputed.zip" "$d"
}

group_eqtl() {
  local d="$REF/eqtl" e=https://molgenis26.gcc.rug.nl/downloads/eqtlgen/cis-eqtl
  fetch_many "$d/eqtlgen" \
    "$e/2019-12-11-cis-eQTLsFDR-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz" \
    "$e/2018-07-18_SNP_AF_for_AlleleB_combined_allele_counts_and_MAF_pos_added.txt.gz" \
    "$e/README_cis"
  fetch https://raw.githubusercontent.com/eQTL-Catalogue/eQTL-Catalogue-resources/master/tabix/tabix_ftp_paths.tsv \
        "$d/eqtlcat/tabix_ftp_paths.tsv"
}

group_vep() {
  local d="$REF/vep" rel="${VEP_RELEASE:-115}"
  fetch "https://ftp.ensembl.org/pub/release-${rel}/variation/indexed_vep_cache/homo_sapiens_vep_${rel}_GRCh37.tar.gz" \
        "$d/homo_sapiens_vep_${rel}_GRCh37.tar.gz"
  unpack "$d/homo_sapiens_vep_${rel}_GRCh37.tar.gz" "$d/cache"
}

group_scdrs() {
  local d="$REF/scdrs"
  fetch https://ndownloader.figshare.com/files/40067134 "$d/TabulaSapiens.h5ad.zip"
  fetch https://ndownloader.figshare.com/files/40066912 "$d/Tabula_Sapiens_metadata.csv"
  unpack "$d/TabulaSapiens.h5ad.zip" "$d"
}

groups=("$@")
[[ ${#groups[@]} -eq 0 ]] && groups=(gene magma gtex ldsc fusion prscs eqtl scdrs vep sbayesrc)
for g in "${groups[@]}"; do
  log "=== group $g ==="
  "group_$g"
done
log "all requested groups done"
du -sh "$REF"/* 2>/dev/null
