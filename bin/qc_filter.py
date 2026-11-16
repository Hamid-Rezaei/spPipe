#!/usr/bin/env python
"""QC_FILTER: compute per-spot QC metrics and remove low-quality spots/genes.

Spots are removed if total counts < --min_counts, detected genes < --min_genes,
or mitochondrial fraction > --max_pct_mt. Genes detected in fewer than
--min_cells spots are removed.

Outputs: {id}_qc.h5ad, {id}_qc_violin.png, {id}_qc_spatial.png, {id}_qc_summary.csv
"""
import argparse

import sppipe_utils as U  # noqa: I001  (sets matplotlib backend first)
import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--input", required=True)
    ap.add_argument("--min_counts", type=int, default=500)
    ap.add_argument("--min_genes", type=int, default=300)
    ap.add_argument("--max_pct_mt", type=float, default=30.0)
    ap.add_argument("--min_cells", type=int, default=10)
    args = ap.parse_args()

    adata = sc.read_h5ad(args.input)
    sid = U.sample_id(adata)
    n_spots0, n_genes0 = adata.shape

    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)

    # Plots are made before filtering so removed spots are visible.
    keys = ["total_counts", "n_genes_by_counts", "pct_counts_mt"]
    fig, axes = plt.subplots(1, 3, figsize=(11, 3.5))
    limits = [args.min_counts, args.min_genes, args.max_pct_mt]
    for ax, k, lim in zip(axes, keys, limits):
        sc.pl.violin(adata, k, jitter=0.4, ax=ax, show=False)
        ax.axhline(lim, color="red", ls="--", lw=1)
    fig.suptitle(f"{sid}: QC metrics (red = filter threshold)")
    fig.savefig(f"{sid}_qc_violin.png")
    plt.close(fig)
    U.spatial_plot(adata, keys + ["ground_truth"], f"{sid}_qc_spatial.png",
                   palette=[U.LAYER_COLORS[k] for k in U.LAYER_ORDER], ncols=4)

    keep = ((adata.obs["total_counts"] >= args.min_counts)
            & (adata.obs["n_genes_by_counts"] >= args.min_genes)
            & (adata.obs["pct_counts_mt"] <= args.max_pct_mt))
    adata = adata[keep].copy()
    sc.pp.filter_genes(adata, min_cells=args.min_cells)

    summary = pd.DataFrame([{
        "sample_id": sid,
        "spots_before": n_spots0, "spots_after": adata.n_obs,
        "genes_before": n_genes0, "genes_after": adata.n_vars,
        "median_counts": float(adata.obs["total_counts"].median()),
        "median_genes": float(adata.obs["n_genes_by_counts"].median()),
        "median_pct_mt": round(float(adata.obs["pct_counts_mt"].median()), 2),
    }])
    summary.to_csv(f"{sid}_qc_summary.csv", index=False)
    print(summary.T.to_string(header=False))

    adata.write_h5ad(f"{sid}_qc.h5ad", compression="gzip")


if __name__ == "__main__":
    main()
