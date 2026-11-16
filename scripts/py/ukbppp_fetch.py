#!/usr/bin/env python3
"""Fetch UKB-PPP (Sun et al. 2023, Nature) European-discovery pQTL summary statistics from
Synapse (project syn51364943) for selected proteins, keeping only cis windows.

Called by the qtl_lookup stage (R/22_qtl_lookup.R) with a JSON spec:
  {"out_dir", "cache_dir", "keep_tars",
   "proteins": [{"gene", "chr", "windows": [[start, end], ...]}],   # GRCh37 windows
   "snps": [rsid, ...]}                                              # lookup SNPs
Writes <out_dir>/lookup.tsv (rows for the lookup SNPs) and <out_dir>/regions/<protein>.tsv.gz
(every variant in the protein's windows, for colocalisation).

Needs SYNAPSE_AUTH_TOKEN (a Synapse personal access token) and acceptance of the UKB-PPP
data-use terms on Synapse. Each protein's tar (~one file per chromosome) is downloaded, the
needed chromosome extracted, and the tar deleted unless keep_tars is true.
"""
import argparse
import gzip
import json
import os
import re
import sys
import tarfile

import pandas as pd
import synapseclient

DISCOVERY_FOLDER = "syn51365303"   # pGWAS summary statistics / European (discovery)
RSID_MAP_FOLDER = "syn51396727"    # Metadata / SNP RSID maps
NAME_RE = re.compile(r"^(?P<gene>.+?)_(?P<uniprot>[A-Z0-9]{6,10})_(?P<oid>OID\d+)_v1_(?P<panel>.+)\.tar$")


def log(msg):
    print(msg, flush=True)


class LazySynapse:
    """Logs in on first use, so a rerun whose listings and regions are all cached needs no token."""

    def __init__(self, cache_root, token):
        self._cache_root, self._token, self._syn = cache_root, token, None

    def __getattr__(self, name):
        if self._syn is None:
            if not self._token:
                sys.exit("SYNAPSE_AUTH_TOKEN is not set and a protein region is not cached yet")
            syn = synapseclient.Synapse(cache_root_dir=self._cache_root, silent=True)
            syn.login(authToken=self._token)
            self._syn = syn
        return getattr(self._syn, name)


def child_files(syn, folder, cache_file):
    if os.path.exists(cache_file):
        with open(cache_file) as fh:
            return json.load(fh)
    files = {c["name"]: c["id"] for c in syn.getChildren(folder, includeTypes=["file"])}
    with open(cache_file, "w") as fh:
        json.dump(files, fh)
    return files


def rsid_map(syn, chrom, cache, maps):
    """GRCh37 position and rsid for every UKB-PPP variant ID on one chromosome."""
    out = os.path.join(cache, f"rsid_map_chr{chrom}.parquet.pkl")
    if os.path.exists(out):
        return pd.read_pickle(out)
    name = next(n for n in maps if re.search(rf"_chr{chrom}_patched", n))
    path = syn.get(maps[name], downloadLocation=cache).path
    m = pd.read_csv(path, sep="\t", dtype=str)
    pos19 = next(c for c in m.columns if c.upper() in ("POS19", "POS_19", "GENPOS19"))
    m = m[["ID", "rsid", pos19]].rename(columns={pos19: "POS19"})
    m["POS19"] = pd.to_numeric(m["POS19"], errors="coerce")
    m.to_pickle(out)
    if os.path.exists(path):
        os.remove(path)
    return m


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--spec", required=True)
    args = ap.parse_args()
    with open(args.spec) as fh:
        spec = json.load(fh)
    out_dir, cache = spec["out_dir"], spec["cache_dir"]
    os.makedirs(os.path.join(out_dir, "regions"), exist_ok=True)
    os.makedirs(cache, exist_ok=True)
    syn = LazySynapse(os.path.join(cache, "synapse_cache"), os.environ.get("SYNAPSE_AUTH_TOKEN"))

    files = child_files(syn, DISCOVERY_FOLDER, os.path.join(cache, "discovery_listing.json"))
    maps = child_files(syn, RSID_MAP_FOLDER, os.path.join(cache, "rsid_map_listing.json"))
    wanted = {p["gene"]: p for p in spec["proteins"]}
    jobs = []
    for name, sid in sorted(files.items()):
        m = NAME_RE.match(name)
        if not m:
            continue
        genes = m.group("gene").split("_")
        hit = [g for g in genes if g in wanted]
        if hit:
            jobs.append((name, sid, wanted[hit[0]]))
    log(f"{len(jobs)} UKB-PPP protein files for {len(wanted)} cis genes")
    snps = set(spec["snps"])
    lookup_rows = []
    for k, (name, sid, prot) in enumerate(jobs, 1):
        protein = name[:-4]
        region_file = os.path.join(out_dir, "regions", f"{protein}.tsv.gz")
        if os.path.exists(region_file):
            reg = pd.read_csv(region_file, sep="\t")
        else:
            chrom = str(prot["chr"]).replace("23", "X")
            log(f"[{k}/{len(jobs)}] {protein} chr{chrom}")
            ent = syn.get(sid, downloadLocation=cache)
            with tarfile.open(ent.path) as tar:
                member = next(mm for mm in tar.getmembers() if re.search(rf"chr{chrom}_", os.path.basename(mm.name)))
                with tar.extractfile(member) as fh:
                    d = pd.read_csv(gzip.open(fh), sep=r"\s+", dtype={"ID": str})
            if not spec.get("keep_tars", False):
                os.remove(ent.path)
            rm = rsid_map(syn, chrom, cache, maps)
            d = d.merge(rm, on="ID", how="inner")
            keep = False
            for start, end in prot["windows"]:
                keep = keep | ((d["POS19"] >= start) & (d["POS19"] <= end))
            reg = d[keep][["rsid", "CHROM", "POS19", "ALLELE0", "ALLELE1", "A1FREQ", "N", "BETA", "SE", "LOG10P"]]
            reg.insert(0, "gene", prot["gene"])
            reg.insert(0, "protein", protein)
            tmp = region_file + ".tmp"
            reg.to_csv(tmp, sep="\t", index=False, compression="gzip")
            os.replace(tmp, region_file)
        lookup_rows.append(reg[reg["rsid"].isin(snps)])
    out = pd.concat(lookup_rows) if lookup_rows else pd.DataFrame(
        columns=["protein", "gene", "rsid", "CHROM", "POS19", "ALLELE0", "ALLELE1", "A1FREQ", "N", "BETA", "SE", "LOG10P"])
    out.to_csv(os.path.join(out_dir, "lookup.tsv"), sep="\t", index=False)
    log(f"wrote {len(out)} lookup rows")


if __name__ == "__main__":
    main()
