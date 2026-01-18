#!/usr/bin/env python3
import argparse, os
import scanpy as sc
import scvi

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
    ap.add_argument("--model_dir", required=True)
    ap.add_argument("--figdir", required=True)
    ap.add_argument("--tabledir", required=True)
    ap.add_argument("--max_epochs", type=int, default=400, help="Maximum epochs for training")
    args = ap.parse_args()

    ensure(args.model_dir); ensure(args.figdir); ensure(args.tabledir)

    adata = sc.read_h5ad(args.in_h5ad)

    # scVI setup (expects counts layer from QC step)
    scvi.model.SCVI.setup_anndata(
        adata,
        layer="counts",
        categorical_covariate_keys=["Sample"],
        continuous_covariate_keys=["pct_counts_mt", "total_counts", "pct_counts_ribo"],
    )

    model = scvi.model.SCVI(adata)
    model.train(max_epochs=args.max_epochs)

    adata.obsm["X_scVI"] = model.get_latent_representation()
    adata.layers["scvi_normalized"] = model.get_normalized_expression(library_size=1e4)

    sc.pp.neighbors(adata, use_rep="X_scVI")
    sc.tl.umap(adata)
    sc.tl.leiden(adata, resolution=0.5)

    sc.pl.umap(adata, color=["leiden", "Sample"], frameon=False, show=False)
    savefig(os.path.join(args.figdir, f"{args.sample_id}_umap_scvi_leiden05.png"))

    # notebook later sets resolution=1 as well
    sc.tl.leiden(adata, resolution=1)

    # scVI DE (table output)
    markers_scvi = model.differential_expression(groupby="leiden")
    markers_scvi.to_csv(os.path.join(args.tabledir, f"{args.sample_id}_markers_scvi_full.csv"))

    # also save a filtered version commonly used in plots
    ms = markers_scvi.copy()
    if "is_de_fdr_0.05" in ms.columns:
        ms = ms[(ms["is_de_fdr_0.05"]) & (ms["lfc_mean"] > 0.5)]
    ms.to_csv(os.path.join(args.tabledir, f"{args.sample_id}_markers_scvi_filtered.csv"))

    model.save(args.model_dir, overwrite=True)

    adata.write_h5ad(args.out_h5ad)

if __name__ == "__main__":
    main()
