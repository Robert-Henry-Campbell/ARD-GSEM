#!/usr/bin/env bash
# Runs every post-GWAS stage that needs the tools from install_postgwas_tools.sh, in dependency
# order, each with its own log under logs/. Stops at the first failing stage. Safe to rerun:
# every stage skips work that already exists (weights, LD windows, scored chromosomes, ...).
#
#   bash scripts/run_post_install.sh            # all stages below
#   bash scripts/run_post_install.sh ld_panel finemap   # a subset, in the order given
set -uo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"
mkdir -p logs

# Synapse token for the UKB-PPP fetch: the client reads ~/.synapseConfig; export the variable too
if [[ -z "${SYNAPSE_AUTH_TOKEN:-}" && -f "$HOME/.synapseConfig" ]]; then
  SYNAPSE_AUTH_TOKEN="$(awk -F'= *' '/^authtoken/{print $2; exit}' "$HOME/.synapseConfig")"
  export SYNAPSE_AUTH_TOKEN
fi

for t in /mnt/sdg/robert/miniforge3/envs/gsem-cli/bin/plink2 /mnt/sdg/robert/miniforge3/envs/gsem-cli/bin/gcta64 \
         tools/magma/magma tools/vep_115.sif tools/fusion_twas/FUSION.assoc_test.R; do
  [[ -e "$t" ]] || { echo "missing $t: run scripts/install_postgwas_tools.sh first" >&2; exit 1; }
done

STAGES=("$@")
if [[ ${#STAGES[@]} -eq 0 ]]; then
  STAGES=(ld_panel finemap cojo annotate qtl_lookup coloc magma twas scdrs gene_priority prs prs_eval split_vcf figures factor_summary)
fi

for st in "${STAGES[@]}"; do
  log="logs/postgwas_${st}_$(date +%Y%m%d_%H%M).log"
  echo "$(date '+%F %T') stage $st -> $log"
  if ! Rscript scripts/run_pipeline.R --sex bothsex --stage "$st" --threads 16 > "$log" 2>&1; then
    echo "$(date '+%F %T') stage $st FAILED; see $log" >&2
    tail -20 "$log" >&2
    exit 1
  fi
done
echo "$(date '+%F %T') all stages done"
