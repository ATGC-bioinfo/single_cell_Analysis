#!/usr/bin/env python3
import argparse, os
import sys
import requests
import numpy as np
import pandas as pd
import scanpy as sc
from io import StringIO

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

def ensure(p): os.makedirs(p, exist_ok=True)

def savefig(path):
    plt.savefig(path, dpi=200, bbox_inches="tight")
    plt.close()

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample_id", required=True)
    ap.add_argument("--in_h5ad", required=True)
    ap.add_argument("--out_h5ad", required=True)
    ap.add_argument("--figdir", required=True)
    ap.add_argument("--tabledir", required=True)
    args = ap.parse_args()

    ensure(args.figdir); ensure(args.tabledir)

    # Read input file with error handling
    try:
        adata = sc.read_h5ad(args.in_h5ad)
    except Exception as e:
        print(f"Error reading input file {args.in_h5ad}: {e}")
        sys.exit(1)

    # QC flags
    adata.var["mt"] = adata.var_names.str.startswith("MT-")

    # Download ribosomal genes with error handling
    ribo_url = "http://software.broadinstitute.org/gsea/msigdb/download_geneset.jsp?geneSetName=KEGG_RIBOSOME&fileType=txt"
    try:
        response = requests.get(ribo_url, timeout=30)
        response.raise_for_status()
        ribo_genes = pd.read_table(StringIO(response.text), skiprows=2, header=None)
        adata.var["ribo"] = adata.var_names.isin(ribo_genes[0].values)
    except Exception as e:
        print(f"Warning: Could not download ribosomal genes: {e}")
        print("Using default ribosomal gene prefix 'RPS' and 'RPL'")
        # Fallback to common ribosomal gene prefixes
        adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))

    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo"], percent_top=None, log1p=False, inplace=True)

    # filter genes
    sc.pp.filter_genes(adata, min_cells=3)

    # violin
    sc.pl.violin(
        adata,
        ["n_genes_by_counts", "total_counts", "pct_counts_mt", "pct_counts_ribo"],
        jitter=0.4,
        multi_panel=True,
        show=False,
    )
    savefig(os.path.join(args.figdir, f"{args.sample_id}_qc_violin.png"))

    # cell filters (same thresholds you used earlier)
    upper_lim = np.quantile(adata.obs["n_genes_by_counts"].values, 0.98)
    adata = adata[adata.obs["n_genes_by_counts"] < upper_lim].copy()
    adata = adata[adata.obs["pct_counts_mt"] < 20].copy()
    adata = adata[adata.obs["pct_counts_ribo"] < 2].copy()

    # keep raw counts layer for scVI
    adata.layers["counts"] = adata.X.copy()

    # Save outputs with error handling
    try:
        adata.obs.to_csv(os.path.join(args.tabledir, f"{args.sample_id}_obs_after_qc.csv"))
        adata.write_h5ad(args.out_h5ad)
    except Exception as e:
        print(f"Error saving outputs: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
