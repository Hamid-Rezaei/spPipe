#!/usr/bin/env python
"""EVALUATE: score both clusterings against the manual layer annotation.

Metrics (computed on spots with a ground-truth label only):
  ARI  adjusted Rand index (0 = random, 1 = perfect agreement)
  NMI  normalised mutual information
Cluster IDs are arbitrary, so for plotting each cluster is matched one-to-one
to the layer it overlaps most (Hungarian algorithm) and drawn in that layer's
colour. Matching affects colours only, never the scores.

Outputs: {id}_metrics.csv, {id}_ari_bar.png, {id}_comparison_spatial.png,
         {id}_confusion.png
"""
import argparse

import sppipe_utils as U  # noqa: I001
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from matplotlib.patches import Patch
from scipy.optimize import linear_sum_assignment
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score

METHODS = {"expr": "Expression-only", "spatial": "Spatially-aware"}
EXTRA_COLORS = ["#a6cee3", "#b2df8a", "#fb9a99", "#cab2d6", "#fdbf6f", "#b15928"]


def match_clusters(truth: pd.Series, clusters: pd.Series) -> dict:
    """Map each cluster to a layer (or to itself if unmatched)."""
    tab = pd.crosstab(clusters, truth)
    rows, cols = linear_sum_assignment(-tab.to_numpy())
    mapping = {tab.index[r]: tab.columns[c] for r, c in zip(rows, cols)}
    for i, cl in enumerate(c for c in tab.index if c not in mapping):
        mapping[cl] = f"extra{i}"
    return mapping


def color_of(label: str) -> str:
    if label in U.LAYER_COLORS:
        return U.LAYER_COLORS[label]
    return EXTRA_COLORS[int(label.replace("extra", "")) % len(EXTRA_COLORS)]


def draw_map(ax, adata, labels, title):
    sid = U.sample_id(adata)
    lib = adata.uns["spatial"][sid]
    sf = lib["scalefactors"]["tissue_lowres_scalef"]
    xy = adata.obsm["spatial"] * sf
    ax.imshow(lib["images"]["lowres"])
    colors = [color_of(lab) if isinstance(lab, str) else "#ffffff00" for lab in labels]
    ax.scatter(xy[:, 0], xy[:, 1], c=colors, s=6, linewidths=0)
    ax.set_xlim(xy[:, 0].min() - 10, xy[:, 0].max() + 10)
    ax.set_ylim(xy[:, 1].max() + 10, xy[:, 1].min() - 10)
    ax.set_title(title, fontsize=10)
    ax.axis("off")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--input", required=True, help="{id}_norm.h5ad")
    ap.add_argument("--expr", required=True, help="{id}_expr_clusters.csv")
    ap.add_argument("--spatial", required=True, help="{id}_spatial_clusters.csv")
    args = ap.parse_args()

    adata = sc.read_h5ad(args.input)
    sid = U.sample_id(adata)
    truth = adata.obs["ground_truth"].astype("object")
    labeled = truth.notna()

    rows, matched = [], {}
    for method, path in [("expr", args.expr), ("spatial", args.spatial)]:
        cl = pd.read_csv(path, index_col="barcode")["cluster"].astype(str).reindex(adata.obs_names)
        if cl.isna().any():
            raise ValueError(f"{path}: {int(cl.isna().sum())} spots without a cluster")
        t, c = truth[labeled], cl[labeled]
        rows.append({
            "sample_id": sid, "method": method,
            "ARI": round(adjusted_rand_score(t, c), 4),
            "NMI": round(normalized_mutual_info_score(t, c), 4),
            "n_clusters": int(cl.nunique()), "n_spots_scored": int(labeled.sum()),
        })
        mapping = match_clusters(t, c)
        matched[method] = cl.map(mapping)
    metrics = pd.DataFrame(rows)
    metrics.to_csv(f"{sid}_metrics.csv", index=False)
    print(metrics.to_string(index=False))

    # Bar plot of ARI / NMI
    fig, ax = plt.subplots(figsize=(4.5, 3.5))
    m = metrics.set_index("method")[["ARI", "NMI"]].rename(index=METHODS)
    m.T.plot.bar(ax=ax, color=["#9e9e9e", "#2b8cbe"], rot=0, width=0.7)
    for cont in ax.containers:
        ax.bar_label(cont, fmt="%.2f", fontsize=8)
    ax.set_ylim(0, 1)
    ax.set_ylabel("score vs manual layers")
    ax.set_title(f"{sid}: clustering accuracy")
    ax.legend(fontsize=8, frameon=False)
    fig.savefig(f"{sid}_ari_bar.png")
    plt.close(fig)

    # Side-by-side maps: truth / expression-only / spatially-aware
    ari = metrics.set_index("method")["ARI"]
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.8))
    draw_map(axes[0], adata, truth.to_numpy(), "Manual annotation")
    for ax, meth in zip(axes[1:], ["expr", "spatial"]):
        draw_map(ax, adata, matched[meth].to_numpy(), f"{METHODS[meth]} (ARI={ari[meth]:.2f})")
    handles = [Patch(color=U.LAYER_COLORS[k], label=k) for k in U.LAYER_ORDER]
    fig.legend(handles=handles, loc="center right", frameon=False, title="Layer / matched cluster")
    fig.suptitle(f"{sid}: cortical layers vs clusters")
    fig.subplots_adjust(right=0.88, wspace=0.03)
    fig.savefig(f"{sid}_comparison_spatial.png")
    plt.close(fig)

    # Confusion matrices (fraction of each layer assigned to each matched cluster)
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    for ax, meth in zip(axes, ["expr", "spatial"]):
        tab = pd.crosstab(truth[labeled], matched[meth][labeled], normalize="index")
        tab = tab.reindex(index=U.LAYER_ORDER,
                          columns=[c for c in U.LAYER_ORDER if c in tab.columns]
                          + sorted(c for c in tab.columns if c not in U.LAYER_ORDER))
        im = ax.imshow(tab.to_numpy(), cmap="Blues", vmin=0, vmax=1)
        ax.set_xticks(range(tab.shape[1]), tab.columns, rotation=45, ha="right", fontsize=8)
        ax.set_yticks(range(tab.shape[0]), tab.index, fontsize=8)
        ax.set_xlabel("matched cluster")
        ax.set_ylabel("manual layer")
        ax.set_title(METHODS[meth], fontsize=10)
        for (i, j), v in np.ndenumerate(tab.to_numpy()):
            if v >= 0.1:
                ax.text(j, i, f"{v:.2f}", ha="center", va="center", fontsize=6,
                        color="white" if v > 0.6 else "black")
    fig.colorbar(im, ax=axes, shrink=0.8, label="fraction of layer")
    fig.suptitle(f"{sid}: layer-to-cluster agreement")
    fig.savefig(f"{sid}_confusion.png")
    plt.close(fig)


if __name__ == "__main__":
    main()
