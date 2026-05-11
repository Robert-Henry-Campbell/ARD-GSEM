#!/usr/bin/env bash
set -euo pipefail

REF="/mnt/sdg/robert/ardmr/GSEM/reference"
mkdir -p "$REF"

echo "Downloading EUR LD scores from Broad Institute..."
wget -q --show-progress -O "$REF/eur_w_ld_chr.tar.bz2" \
  https://data.broadinstitute.org/alkesgroup/LDSCORE/eur_w_ld_chr.tar.bz2

echo "Extracting LD scores..."
tar -xjf "$REF/eur_w_ld_chr.tar.bz2" -C "$REF"
rm "$REF/eur_w_ld_chr.tar.bz2"

echo "Downloading HapMap3 SNP list..."
wget -q --show-progress -O "$REF/w_hm3.snplist.bz2" \
  https://data.broadinstitute.org/alkesgroup/LDSCORE/w_hm3.snplist.bz2

echo "Extracting..."
bunzip2 "$REF/w_hm3.snplist.bz2"

echo ""
echo "=== Verification ==="
echo "LD score files: $(ls "$REF/eur_w_ld_chr/" | grep -c 'l2.ldscore.gz')"
echo "HM3 SNPs: $(wc -l < "$REF/w_hm3.snplist") lines"
echo "Done."
