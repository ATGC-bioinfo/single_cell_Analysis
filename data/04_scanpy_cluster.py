#!/usr/bin/env python3
import argparse, os
import scanpy as sc

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

    adata = sc.read_h5ad(args.in_h5ad)

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.raw = adata

    sc.pp.highly_variable_genes(adata, n_top_genes=2000)
    sc.pl.highly_variable_genes(adata, show=False)
    savefig(os.path.join(args.figdir, f"{args.sample_id}_hvg.png"))

    adata = adata[:, adata.var["highly_variable"]].copy()

    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, svd_solver="arpack")
    sc.pl.pca_variance_ratio(adata, log=True, n_pcs=50, show=False)
    savefig(os.path.join(args.figdir, f"{args.sample_id}_pca_variance.png"))

    sc.pp.neighbors(adata, n_pcs=30)
    sc.tl.umap(adata)
    sc.tl.leiden(adata, resolution=0.5)

    sc.pl.umap(adata, color=["leiden"], show=False)
    savefig(os.path.join(args.figdir, f"{args.sample_id}_umap_scanpy_leiden05.png"))

    adata.write_h5ad(args.out_h5ad)

if __name__ == "__main__":
    main()
