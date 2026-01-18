#!/usr/bin/env python3
import argparse, os, math
import numpy as np
import pandas as pd
import scanpy as sc

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

def ensure(p): os.makedirs(p, exist_ok=True)

def savefig(path):
    plt.savefig(path, dpi=200, bbox_inches="tight")
    plt.close()

def volcano_plot(df: pd.DataFrame, out_png: str, title: str):
    # expects scanpy rank_genes_groups_df columns:
    # names, logfoldchanges, pvals_adj, scores, group
    x = df["logfoldchanges"].to_numpy()
    p = df["pvals_adj"].to_numpy()
    y = -np.log10(np.clip(p, 1e-300, 1.0))

    plt.figure()
    plt.scatter(x, y, s=6)
    plt.axvline(0.5, linestyle="--")
    plt.axvline(-0.5, linestyle="--")
    plt.axhline(-math.log10(0.05), linestyle="--")
    plt.title(title)
    plt.xlabel("logFC")
    plt.ylabel("-log10(adj p)")
    savefig(out_png)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample_id", required=True)
    ap.add_argument("--in_h5ad", required=True)
    ap.add_argument("--figdir", required=True)
    ap.add_argument("--tabledir", required=True)
    args = ap.parse_args()

    ensure(args.figdir); ensure(args.tabledir)

    adata = sc.read_h5ad(args.in_h5ad)

    # Scanpy DE on final leiden (after scVI step it’s resolution=1)
    sc.tl.rank_genes_groups(adata, "leiden")

    markers = sc.get.rank_genes_groups_df(adata, None)
    markers.to_csv(os.path.join(args.tabledir, f"{args.sample_id}_markers_scanpy_full.csv"), index=False)

    # filtered markers
    filt = markers[(markers["pvals_adj"] < 0.05) & (markers["logfoldchanges"] > 0.5)].copy()
    filt.to_csv(os.path.join(args.tabledir, f"{args.sample_id}_markers_scanpy_filtered.csv"), index=False)

    # Regenerate dendrogram to ensure it matches leiden groups
    sc.tl.dendrogram(adata, "leiden")
    
    # ✅ Heatmap (top markers per cluster)
    sc.pl.rank_genes_groups_heatmap(
        adata,
        n_genes=10,
        groupby="leiden",
        dendrogram=False,
        show=False
    )
    savefig(os.path.join(args.figdir, f"{args.sample_id}_heatmap_top_markers_per_cluster.png"))

    # ✅ Volcano plots per cluster
    for cl in sorted(markers["group"].unique(), key=lambda z: int(z) if str(z).isdigit() else str(z)):
        dfc = markers[markers["group"] == cl].copy()
        volcano_plot(
            dfc,
            os.path.join(args.figdir, f"{args.sample_id}_volcano_cluster_{cl}.png"),
            f"{args.sample_id} | Cluster {cl}"
        )

if __name__ == "__main__":
    main()
