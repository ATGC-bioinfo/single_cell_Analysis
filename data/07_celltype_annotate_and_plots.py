#!/usr/bin/env python3
import argparse, os
import pandas as pd
import scanpy as sc

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

CELL_TYPE_MAP = {
    "0":"Macrophage","1":"Fibroblast","2":"CD4+ T-cell","3":"AT2","4":"AT1",
    "5":"CD8+ T-cell","6":"Endothelial cell","7":"Plasma cell","8":"Macrophage",
    "9":"AT2","10":"Fibroblast","11":"Fibroblast","12":"Macrophage","13":"Macrophage",
    "14":"Airway epithelial","15":"Airway epithelial","16":"Monocyte","17":"Airway epithelial",
    "18":"B-cell","19":"Aerocyte","20":"Airway epithelial","21":"Smooth muscle cell",
    "22":"Cycling T/NK","23":"Neuronal cell","24":"Denditic cell","25":"Pericyte",
    "26":"Fibroblast","27":"Erythroid-like","28":"Macrophage"
}

# ✅ canonical markers for dotplot (you can expand to match notebook exactly)
CANONICAL_MARKERS = {
    "AT2": ["SFTPC", "SFTPA1", "SLC34A2"],
    "AT1": ["AGER", "PDPN", "CAV1"],
    "Macrophage": ["LYZ", "CST3", "LGALS3"],
    "Monocyte": ["S100A8", "S100A9", "VCAN"],
    "Endothelial cell": ["PECAM1", "VWF", "KDR"],
    "Fibroblast": ["COL1A1", "COL1A2", "DCN"],
    "B-cell": ["MS4A1", "CD79A"],
    "Plasma cell": ["MZB1", "XBP1"],
    "CD4+ T-cell": ["CD3D", "IL7R", "LTB"],
    "CD8+ T-cell": ["CD3D", "NKG7", "GZMB"],
    "Airway epithelial": ["KRT18", "KRT19"],
}

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

    # annotate
    adata.obs["cell type"] = adata.obs["leiden"].map(CELL_TYPE_MAP).astype("category")

    # UMAP cell types
    sc.pl.umap(adata, color=["cell type"], frameon=False, show=False)
    savefig(os.path.join(args.figdir, f"{args.sample_id}_umap_celltypes.png"))

    # ✅ Cell-type proportion bar plot
    ct = (
        adata.obs["cell type"]
        .value_counts(dropna=False, normalize=True)
        .mul(100)
        .reset_index()
    )
    ct.columns = ["cell_type", "percent"]
    ct.to_csv(os.path.join(args.tabledir, f"{args.sample_id}_celltype_proportions.csv"), index=False)

    plt.figure()
    plt.bar(ct["cell_type"].astype(str), ct["percent"])
    plt.xticks(rotation=45, ha="right")
    plt.ylabel("Percent of cells")
    plt.title(f"{args.sample_id} cell-type proportions")
    savefig(os.path.join(args.figdir, f"{args.sample_id}_celltype_proportions_bar.png"))

    # ✅ DotPlot (canonical markers)
    # filter to genes present
    markers_present = {}
    gene_set = set(adata.var_names)
    for k, genes in CANONICAL_MARKERS.items():
        gg = [g for g in genes if g in gene_set]
        if len(gg) > 0:
            markers_present[k] = gg

    if len(markers_present) > 0:
        sc.pl.dotplot(
            adata,
            markers_present,
            groupby="cell type",
            show=False
        )
        savefig(os.path.join(args.figdir, f"{args.sample_id}_dotplot_canonical_markers.png"))

    adata.write_h5ad(args.out_h5ad)

if __name__ == "__main__":
    main()
