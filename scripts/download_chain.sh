#!/usr/bin/env bash
# Download UCSC hg38ToHg19 chain file used for the bothsex_meta liftover step.
# The chain file is read by rtracklayer::import.chain in
# liftover_grch38_to_grch37() (see R/utils.R).
#
# rtracklayer must be installed separately (Bioconductor):
#   Rscript -e 'if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager"); BiocManager::install("rtracklayer"); renv::snapshot()'
set -euo pipefail

# Resolve script-relative reference/chains directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
DEST_DIR="$PROJECT_ROOT/reference/chains"
mkdir -p "$DEST_DIR"
cd "$DEST_DIR"

URL="https://hgdownload.cse.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz"
GZ="hg38ToHg19.over.chain.gz"
CHAIN="hg38ToHg19.over.chain"

if [[ -f "$CHAIN" ]]; then
  echo "Chain already present: $DEST_DIR/$CHAIN"
  exit 0
fi

echo "Downloading $URL -> $DEST_DIR/$GZ"
# Retry 5x with 30s backoff for transient UCSC failures.
wget -t 5 --waitretry=30 -O "$GZ" "$URL"

echo "Gunzipping $GZ"
gunzip -f "$GZ"

if [[ ! -f "$CHAIN" ]]; then
  echo "ERROR: chain file not present after gunzip" >&2
  exit 1
fi

echo "Done: $DEST_DIR/$CHAIN ($(du -h "$CHAIN" | cut -f1))"
