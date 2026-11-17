#!/usr/bin/env bash
# Installs the external tools used by the post-GWAS stages. Idempotent: anything already present
# is skipped, so the script can be re-run after a failure or to pick up ANNOVAR later.
#
#   conda envs under $ENVS (separate envs avoid solver conflicts):
#     gsem-cli   plink2, gcta, bcftools/htslib (the system tabix/bgzip are broken), samtools,
#                ucsc-liftover
#     gsem-py    python 3.11 + gwaslab, synapseclient, pgenlib, pysam, PRS-CS
#     gsem-scdrs python 3.10 + scanpy + scdrs
#   $TOOLS:
#     magma/     official MAGMA v1.10 static binary (the bioconda "magma" package is a
#                different tool)
#     *.sif      apptainer images for VEP and SBayesRC (docker is not usable here)
#     fusion_twas/ + rlib/   FUSION and an R library for its scripts
#     annovar/   only if you supply annovar.latest.tar.gz (licence required)
#     ldsc/      patched copy of the GSEM-finn-ukb-mvp Python 3 LDSC port (S-LDSC, LDSC-SEG)
#   project renv: susieR, coloc, mockery, pgenlibr, bigsnpr (LDpred2)
#
# Usage (from the GSEM project root): bash scripts/install_postgwas_tools.sh
# ANNOVAR: register at https://www.openbioinformatics.org/annovar/annovar_download_form.php,
# put annovar.latest.tar.gz in tools/, then re-run.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TOOLS="${TOOLS:-$ROOT/tools}"
MAMBA="${MAMBA:-/mnt/sdg/robert/miniforge3/bin/mamba}"
ENVS="${ENVS:-/mnt/sdg/robert/miniforge3/envs}"
VEP_RELEASE="${VEP_RELEASE:-115}"
# / has little free space: keep every cache on /mnt/sdg.
CACHE="${CACHE:-/mnt/sdg/robert/.cache}"
export APPTAINER_CACHEDIR="${APPTAINER_CACHEDIR:-$CACHE/apptainer}"
export APPTAINER_TMPDIR="${APPTAINER_TMPDIR:-/mnt/sdg/robert/r-tmp}"
export PIP_CACHE_DIR="${PIP_CACHE_DIR:-$CACHE/pip}"
export CONDA_PKGS_DIRS="${CONDA_PKGS_DIRS:-$CACHE/conda-pkgs}"
export TMPDIR="${TMPDIR:-/mnt/sdg/robert/r-tmp}"
mkdir -p "$TOOLS" "$APPTAINER_CACHEDIR" "$APPTAINER_TMPDIR" "$PIP_CACHE_DIR" "$CONDA_PKGS_DIRS" "$TMPDIR"

log() { printf '%s [install] %s\n' "$(date '+%F %T')" "$*"; }

make_env() { # name, conda packages...
  local name="$1"; shift
  if [[ ! -d "$ENVS/$name/conda-meta" ]]; then
    log "creating conda env $name"
    "$MAMBA" create -y -p "$ENVS/$name" -c conda-forge -c bioconda --strict-channel-priority "$@"
  else
    log "conda env present: $name"
  fi
}

pip_need() { # env, import-name, pip-name
  local py="$ENVS/$1/bin/python"
  if ! "$py" -c "import $2" 2>/dev/null; then
    log "pip install $3 into $1"
    "$py" -m pip install --quiet "$3"
  fi
}

# --- 1. conda envs ----------------------------------------------------------------------------
make_env gsem-cli plink2 gcta bcftools htslib samtools ucsc-liftover
make_env gsem-py python=3.11 numpy scipy pandas h5py pysam matplotlib prscs pip
pip_need gsem-py gwaslab gwaslab
pip_need gsem-py synapseclient synapseclient
pip_need gsem-py pgenlib Pgenlib
make_env gsem-scdrs python=3.10 scanpy anndata statsmodels pip
pip_need gsem-scdrs scdrs scdrs

# --- 2. MAGMA v1.10 (official static Linux build) ---------------------------------------------
if [[ ! -x "$TOOLS/magma/magma" ]]; then
  log "downloading MAGMA"
  mkdir -p "$TOOLS/magma"
  curl -sS -L --fail -o "$TOOLS/magma/magma_static.zip" \
    https://vu.data.surfsara.nl/index.php/s/lxDgt2dNdNr6DYt/download
  unzip -q -o "$TOOLS/magma/magma_static.zip" -d "$TOOLS/magma"
  chmod +x "$TOOLS/magma/magma"
fi

# --- 3. apptainer images ----------------------------------------------------------------------
if [[ ! -s "$TOOLS/vep_${VEP_RELEASE}.sif" ]]; then
  log "pulling VEP release ${VEP_RELEASE} (must match the GRCh37 cache release)"
  apptainer pull "$TOOLS/vep_${VEP_RELEASE}.sif" "docker://ensemblorg/ensembl-vep:release_${VEP_RELEASE}.0"
fi
if [[ ! -s "$TOOLS/sbayesrc.sif" ]]; then
  # Container (Ubuntu 24.04) rather than the R package: the SBayesRC authors report wrong LD /
  # eigen results with Ubuntu 22.04's OpenBLAS 0.3.20, which is this server's OS.
  log "pulling SBayesRC"
  apptainer pull "$TOOLS/sbayesrc.sif" docker://zhiliz/sbayesrc
fi

# --- 4. FUSION (TWAS) -------------------------------------------------------------------------
if [[ ! -f "$TOOLS/fusion_twas/FUSION.assoc_test.R" ]]; then
  log "cloning FUSION"
  git clone --depth 1 https://github.com/gusevlab/fusion_twas "$TOOLS/fusion_twas"
fi
# FUSION's scripts run with `Rscript --vanilla` and R_LIBS=$TOOLS/rlib, outside the project renv.
mkdir -p "$TOOLS/rlib"
R_LIBS="$TOOLS/rlib" Rscript --vanilla -e '
  repos <- "https://cloud.r-project.org"
  need <- c("optparse", "RColorBrewer", "remotes", "Rcpp", "RcppEigen", "here")   # FUSION >= 2026-03 sources utils/ via here()
  miss <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss)) install.packages(miss, repos = repos)
  if (!requireNamespace("plink2R", quietly = TRUE))
    remotes::install_github("gabraham/plink2R/plink2R", upgrade = "never")
  stopifnot(requireNamespace("plink2R", quietly = TRUE))'

# --- 5. project renv packages -----------------------------------------------------------------
# susieR / coloc / mockery are linked from the shared renv cache; pgenlibr and bigsnpr build
# from CRAN.
( cd "$ROOT" && Rscript -e '
  need <- c("susieR", "coloc", "mockery", "pgenlibr", "bigsnpr", "patchwork", "ggrepel", "survival")
  miss <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss)) renv::install(miss, prompt = FALSE)
  stopifnot(all(vapply(need, requireNamespace, logical(1), quietly = TRUE)))' )

# --- 6. ANNOVAR (optional; licence required) --------------------------------------------------
if [[ -f "$TOOLS/annovar.latest.tar.gz" && ! -d "$TOOLS/annovar" ]]; then
  log "unpacking ANNOVAR"
  tar -xzf "$TOOLS/annovar.latest.tar.gz" -C "$TOOLS"
fi
if [[ -d "$TOOLS/annovar" ]]; then
  for db in refGene ensGene avsnp150 gnomad211_genome clinvar_20240917; do
    if ! ls "$TOOLS/annovar/humandb/hg19_${db}"* >/dev/null 2>&1; then
      log "ANNOVAR humandb: $db"
      ( cd "$TOOLS/annovar" && ./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar "$db" humandb/ )
    fi
  done
else
  log "ANNOVAR not supplied; the annotate stage will run VEP only"
fi

# --- 7. LDSC (Python 3 port) for the sldsc stage ---------------------------------------------
# A copy of the GSEM-finn-ukb-mvp project's Python 3 port of bulik/ldsc (the original is left
# untouched), patched where its --l2 and partitioned-h2 paths break under Python 3, pandas 2 and
# bitarray >= 3. Validated: chr22 LD scores reproduce the published baseline_v1.2 baseL2.
LDSC_SRC="${LDSC_SRC:-$ROOT/../GSEM-finn-ukb-mvp/repo/ldsc}"
if [[ ! -d "$TOOLS/ldsc" ]]; then
  log "copying the LDSC port from $LDSC_SRC"
  rsync -a --exclude .git "$LDSC_SRC/" "$TOOLS/ldsc/"
fi
python3 -I - "$TOOLS/ldsc" <<'EOF'
import pathlib, sys
d = pathlib.Path(sys.argv[1])
def patch(rel, old, new, count=1):
    p = d / rel
    s = p.read_text()
    n = s.count(old)
    if n == 0:
        return  # already patched
    assert n == count, f"{rel}: expected {count} x {old!r}, found {n}"
    p.write_text(s.replace(old, new))
# bitarray >= 3: decode() returns an iterator
patch("ldscore/ldscore.py", 'np.array(slice.decode(self._bedcode), dtype="float64")',
      'np.array(list(slice.decode(self._bedcode)), dtype="float64")')
# .M files were opened in binary mode
patch("ldsc.py", "+'.M','wb')", "+'.M','w')")
patch("ldsc.py", "+'.M_5_50','wb')", "+'.M_5_50','w')")
# np.c_ of mixed columns gives object dtype (float_format ignored, np.linalg.cond fails)
patch("ldsc.py", "    df.columns = new_colnames\n    if args.print_snps:",
      "    df.columns = new_colnames\n    for c in new_colnames[4:]:  # np.c_ of mixed columns gives object dtype\n"
      "        df[c] = df[c].astype(float)\n    if args.print_snps:")
# the chi^2 filters index with a numpy mask, not a Series
patch("ldscore/sumstats.py", "sumstats = sumstats.loc[ii.values, :]", "sumstats = sumstats.loc[ii, :]", count=2)
# annot() relied on fh leaking out of a Python 2 list comprehension
patch("ldscore/parse.py", "chrs = get_present_chrs(fh, num+1)\n        for i, fh in enumerate(fh_list):",
      "chrs = get_present_chrs(fh_list[-1], num+1)\n        for i, fh in enumerate(fh_list):")
# Python 3 traceback.format_exc() takes no exception argument
patch("ldsc.py", "traceback.format_exc(ex)", "traceback.format_exc()")
patch("ldscore/sumstats.py", "traceback.format_exc(ex)", "traceback.format_exc()")
EOF

# --- 8. PRS-CS (bioconda) for scipy >= 1.12 -------------------------------------------------------
# The bioconda PRS-CS imports scipy as sp and calls numpy aliases (sp.sqrt, sp.log, sp.array, ...)
# that scipy 1.12 removed; a copy in tools/prscs imports numpy as sp instead (config tools$prscs).
mkdir -p "$TOOLS/prscs"
for f in PRScs.py parse_genet.py mcmc_gtb.py gigrnd.py; do cp "$ENVS/gsem-py/bin/$f" "$TOOLS/prscs/"; done
sed -i 's/^import scipy as sp$/import numpy as sp  # scipy >= 1.12 dropped the numpy aliases (sqrt, log, array, ...)/' \
  "$TOOLS/prscs/parse_genet.py" "$TOOLS/prscs/mcmc_gtb.py"
# numpy >= 2: float() / %-formatting need 0-d values; PRS-CS passes one-element arrays
python3 -I - "$TOOLS/prscs" <<'PYEOF'
import pathlib, sys
d = pathlib.Path(sys.argv[1])
def patch(fn, old, new):
    p = d / fn; s = p.read_text()
    if new in s: return
    assert s.count(old) == 1, (fn, old); p.write_text(s.replace(old, new))
patch("gigrnd.py", "import math\n", "import math\nimport numpy as np\n")
patch("gigrnd.py", "    p = float(p); a = float(a); b = float(b)\n",
      "    p = float(np.ravel(p)[0]); a = float(np.ravel(a)[0]); b = float(np.ravel(b)[0])  # numpy >= 2\n")
patch("mcmc_gtb.py", "quad += sp.dot(sp.dot(beta[idx_blk].T, dinvt), beta[idx_blk])\n",
      "quad += sp.dot(sp.dot(beta[idx_blk].T, dinvt), beta[idx_blk]).item()\n")
patch("mcmc_gtb.py", "err = max(n/2.0*(1.0-2.0*sum(beta*beta_mrg)+quad), n/2.0*sum(beta**2/psi))",
      "err = max(n/2.0*(1.0-2.0*sp.sum(beta*beta_mrg)+quad), n/2.0*sp.sum(beta**2/psi))")
patch("mcmc_gtb.py", "phi = random.gamma(p*b+0.5, 1.0/(sum(delta)+w))", "phi = random.gamma(p*b+0.5, 1.0/(sp.sum(delta)+w))")
patch("mcmc_gtb.py", "% (chrom, snp, bp, a1, a2, beta))", "% (chrom, snp, bp, a1, a2, float(sp.ravel(beta)[0])))")
patch("mcmc_gtb.py", "% (snp, psi))", "% (snp, float(sp.ravel(psi)[0])))")
patch("mcmc_gtb.py", "% phi_est )", "% float(sp.ravel(phi_est)[0]) )")
PYEOF

# --- 9. scDRS 1.0.2 for numpy >= 2 ------------------------------------------------------------
# perform-downstream's Geary's C test casts with np.float_, which numpy 2.0 removed (it was np.float64)
SCDRS_PKG="$("$ENVS/gsem-scdrs/bin/python" -I -c 'import importlib.util as u; print(u.find_spec("scdrs").submodule_search_locations[0])')"
sed -i -E 's/np\.float_([,)])/np.float64\1/g' "$SCDRS_PKG/method.py"

# --- 10. report -------------------------------------------------------------------------------
log "versions:"
"$ENVS/gsem-cli/bin/plink2" --version | head -1
"$ENVS/gsem-cli/bin/bcftools" --version | head -1
ls "$ENVS/gsem-cli/bin" | grep -E '^(gcta|liftOver)' | tr '\n' ' '; echo
"$TOOLS/magma/magma" --version 2>&1 | head -1 || true
"$ENVS/gsem-py/bin/python" -c "import gwaslab, synapseclient, pgenlib; print('gwaslab', gwaslab.__version__, '| synapseclient', synapseclient.__version__)"
"$ENVS/gsem-scdrs/bin/python" -c "import scdrs; print('scdrs', scdrs.__version__)"
ls -la "$TOOLS"/*.sif
# the prs stage calls the SBayesRC R API inside the image (its entrypoint is a CLI wrapper)
apptainer exec "$TOOLS/sbayesrc.sif" Rscript --vanilla -e 'cat("SBayesRC", as.character(packageVersion("SBayesRC")), "\n")' ||
  log "WARNING: Rscript/SBayesRC not usable inside sbayesrc.sif; the prs stage's sbayesrc method will fail"
log "done"
