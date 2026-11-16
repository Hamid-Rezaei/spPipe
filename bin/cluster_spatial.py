#!/usr/bin/env python
"""CLUSTER_SPATIAL: spatially-aware clustering (BANKSY-style neighbour augmentation).

1. Build a spatial graph of the Visium hexagonal grid with squidpy
   (each spot's neighbours within --n_rings rings).
2. For each spot, average the PCs of its spatial neighbours.
3. Concatenate own PCs * sqrt(1 - lambda) with neighbour-mean PCs * sqrt(lambda)
   (Singhal et al. 2024, BANKSY). lambda = 0 is expression-only clustering;
   larger lambda puts more weight on the tissue neighbourhood.
4. kNN graph + Leiden on this joint representation, tuned to --n_clusters.

The same PCs as CLUSTER_EXPR are used, so the only difference is step 1-3.

Outputs: {id}_spatial_clusters.csv, {id}_spatial_umap.png, {id}_spatial_spatial.png
"""
import argparse

import sppipe_utils as U  # noqa: I001
import numpy as np
import scanpy as sc
import scipy.sparse as sp
import squidpy as sq


def neighbour_mean(adata, X, n_rings):
    sq.gr.spatial_neighbors(adata, coord_type="grid", n_neighs=6, n_rings=n_rings)
    W = adata.obsp["spatial_connectivities"].astype(float)
    W = W - sp.diags(W.diagonal())  # exclude the spot itself
    deg = np.asarray(W.sum(axis=1)).ravel()
    deg[deg == 0] = 1.0  # isolated spots keep a zero neighbour vector
    print(f"spatial graph: median {np.median(deg):.0f} neighbours per spot (n_rings={n_rings})")
    return sp.diags(1.0 / deg) @ W @ X


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--input", required=True)
    ap.add_argument("--n_pcs", type=int, default=20)
    ap.add_argument("--n_clusters", type=int, default=7)
    ap.add_argument("--n_neighbors", type=int, default=15)
    ap.add_argument("--n_rings", type=int, default=2)
    ap.add_argument("--lambda_", dest="lam", type=float, default=0.8)
    args = ap.parse_args()

    adata = sc.read_h5ad(args.input)
    X = adata.obsm["X_pca"][:, : args.n_pcs]
    Xn = neighbour_mean(adata, X, args.n_rings)
    adata.obsm["X_spatial"] = np.hstack([np.sqrt(1 - args.lam) * X, np.sqrt(args.lam) * Xn])
    U.cluster_and_save(adata, "X_spatial", "spatial", args.n_clusters, args.n_neighbors)


if __name__ == "__main__":
    main()
