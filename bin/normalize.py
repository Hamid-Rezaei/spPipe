#!/usr/bin/env python
"""NORMALIZE: library-size normalisation, log1p, highly variable genes, PCA.

Raw counts are kept in layers['counts']. PCA is computed on the HVGs and is
shared by both clustering methods, so they start from identical input.

Outputs: {id}_norm.h5ad, {id}_pca_variance.png
"""
import argparse

import sppipe_utils as U  # noqa: I001
import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--input", required=True)
    ap.add_argument("--n_hvg", type=int, default=3000)
    ap.add_argument("--n_pcs", type=int, default=50)
    args = ap.parse_args()
    U.set_seed()

    adata = sc.read_h5ad(args.input)
    sid = U.sample_id(adata)

    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=args.n_hvg, flavor="seurat")
    sc.pp.pca(adata, n_comps=args.n_pcs, mask_var="highly_variable",
              random_state=U.SEED, svd_solver="arpack")

    vr = adata.uns["pca"]["variance_ratio"]
    fig, ax = plt.subplots(figsize=(5, 3.5))
    ax.plot(np.arange(1, len(vr) + 1), vr, "o-", ms=3)
    ax.set_xlabel("Principal component")
    ax.set_ylabel("Variance ratio")
    ax.set_yscale("log")
    ax.set_title(f"{sid}: PCA variance ({args.n_hvg} HVGs)")
    fig.savefig(f"{sid}_pca_variance.png")
    plt.close(fig)

    print(f"{sid}: {adata.n_obs} spots, {int(adata.var.highly_variable.sum())} HVGs, "
          f"{args.n_pcs} PCs explain {vr.sum():.1%} of HVG variance")
    adata.write_h5ad(f"{sid}_norm.h5ad", compression="gzip")


if __name__ == "__main__":
    main()
