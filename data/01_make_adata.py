#!/usr/bin/env python3
import argparse
import scanpy as sc

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample_id", required=True)
    ap.add_argument("--counts_csv", required=True)
    ap.add_argument("--out_h5ad", required=True)
    args = ap.parse_args()

    adata = sc.read_csv(args.counts_csv).T
    adata.obs["Sample"] = args.sample_id
    adata.write_h5ad(args.out_h5ad)

if __name__ == "__main__":
    main()
