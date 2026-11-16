#!/usr/bin/env python
"""CLUSTER_EXPR: expression-only clustering (the standard single-cell approach).

kNN graph on the first --n_pcs principal components, then Leiden. Spot
locations are ignored. The Leiden resolution is tuned automatically to give
--n_clusters clusters (7 = L1-L6 + WM) so both methods are compared fairly.

Outputs: {id}_expr_clusters.csv, {id}_expr_umap.png, {id}_expr_spatial.png
"""
import argparse

import sppipe_utils as U  # noqa: I001
import scanpy as sc


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--input", required=True)
    ap.add_argument("--n_pcs", type=int, default=20)
    ap.add_argument("--n_clusters", type=int, default=7)
    ap.add_argument("--n_neighbors", type=int, default=15)
    args = ap.parse_args()

    adata = sc.read_h5ad(args.input)
    adata.obsm["X_expr"] = adata.obsm["X_pca"][:, : args.n_pcs]
    U.cluster_and_save(adata, "X_expr", "expr", args.n_clusters, args.n_neighbors)


if __name__ == "__main__":
    main()
