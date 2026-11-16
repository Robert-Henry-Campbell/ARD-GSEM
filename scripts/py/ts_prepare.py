#!/usr/bin/env python3
"""Prepare Tabula Sapiens for the scdrs stage: seeded subsample of at most --max-per-type cells
per cell_ontology_class (raw counts), a covariate file (const, n_genes) and per-cell-type mean
log1p(CPM) expression (genes x cell types) for MAGMA gene-property analysis.

Usage: ts_prepare.py --h5ad TabulaSapiens.h5ad --out-dir DIR [--max-per-type 1000] [--seed 20261108]
"""
import argparse
import os

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--h5ad", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--max-per-type", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=20261108)
    ap.add_argument("--celltype-col", default="cell_ontology_class")
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)
    full = ad.read_h5ad(a.h5ad, backed="r")
    obs = full.obs
    rng = np.random.default_rng(a.seed)
    keep = []
    for _, idx in obs.groupby(a.celltype_col, observed=True).indices.items():
        keep.extend(rng.choice(idx, size=min(len(idx), a.max_per_type), replace=False))
    keep = np.sort(np.array(keep))
    sub = full[keep].to_memory()
    counts = sub.layers["raw_counts"] if "raw_counts" in sub.layers else sub.X
    counts = sp.csr_matrix(counts)
    sub.X = counts
    sub.layers = {}
    if "gene_symbol" in sub.var.columns:
        sub.var_names = sub.var["gene_symbol"].astype(str)
    sub.var_names_make_unique()
    # the index inherits the name "gene_symbol", and after make_unique its values differ from that
    # column, which anndata refuses to write: drop the index name
    sub.var.index.name = None
    # every output is written to a temporary name and renamed, so a crash cannot leave a partial file
    # that a rerun mistakes for a finished one
    h5 = os.path.join(a.out_dir, "ts_subsample.h5ad")
    sub.write_h5ad(h5 + ".tmp.h5ad")
    os.replace(h5 + ".tmp.h5ad", h5)
    n_genes = np.asarray((counts > 0).sum(axis=1)).ravel()
    cov = os.path.join(a.out_dir, "ts_cov.tsv")
    pd.DataFrame({"index": sub.obs_names, "const": 1, "n_genes": n_genes}).to_csv(cov + ".tmp", sep="\t", index=False)
    os.replace(cov + ".tmp", cov)
    lib = np.asarray(counts.sum(axis=1)).ravel()
    norm = sp.diags(1e6 / np.maximum(lib, 1)) @ counts
    norm.data = np.log1p(norm.data)
    groups = sub.obs[a.celltype_col].astype(str).values
    means = {}
    for g in np.unique(groups):
        m = groups == g
        means[g] = np.asarray(norm[m].mean(axis=0)).ravel()
    mean_f = os.path.join(a.out_dir, "ts_celltype_mean_log1pcpm.tsv")
    pd.DataFrame(means, index=sub.var_names).to_csv(mean_f + ".tmp", sep="\t")
    os.replace(mean_f + ".tmp", mean_f)
    print(f"{sub.n_obs} cells, {len(means)} cell types")


if __name__ == "__main__":
    main()
