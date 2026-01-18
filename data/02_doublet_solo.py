#!/usr/bin/env python3
import argparse
import scanpy as sc
import scvi

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in_h5ad", required=True)
    ap.add_argument("--out_h5ad", required=True)
    ap.add_argument("--out_doublets", required=True)
    ap.add_argument("--max_epochs", type=int, default=400, help="Maximum epochs for training")
    args = ap.parse_args()

    adata = sc.read_h5ad(args.in_h5ad)

    ad_for_doublet = adata.copy()
    sc.pp.filter_genes(ad_for_doublet, min_cells=10)
    sc.pp.highly_variable_genes(ad_for_doublet, n_top_genes=2000, subset=True, flavor="seurat_v3")

    scvi.model.SCVI.setup_anndata(ad_for_doublet)
    vae = scvi.model.SCVI(ad_for_doublet)
    vae.train(max_epochs=args.max_epochs)

    solo = scvi.external.SOLO.from_scvi_model(vae)
    solo.train(max_epochs=args.max_epochs)

    df = solo.predict()
    df["prediction"] = solo.predict(soft=False)
    df.index = df.index.map(lambda x: x[:-2])  # keep same behavior as notebook
    df["dif"] = df["doublet"] - df["singlet"]

    doublets = df[(df["prediction"] == "doublet") & (df["dif"] > 1)].copy()
    doublets.to_csv(args.out_doublets)

    adata.obs["doublet"] = adata.obs_names.isin(doublets.index)
    adata = adata[~adata.obs["doublet"]].copy()
    adata.write_h5ad(args.out_h5ad)

if __name__ == "__main__":
    main()
