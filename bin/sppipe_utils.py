"""Shared helpers for the spPipe scripts (imported from bin/)."""
import os

os.environ.setdefault("MPLCONFIGDIR", os.path.join(os.getcwd(), ".mplconfig"))
os.environ.setdefault("NUMBA_CACHE_DIR", os.path.join(os.getcwd(), ".numba_cache"))

import random
import warnings

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import scanpy as sc  # noqa: E402

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)

SEED = 0
LAYER_ORDER = ["L1", "L2", "L3", "L4", "L5", "L6", "WM"]
# One fixed colour per cortical layer, reused for matched clusters.
LAYER_COLORS = {
    "L1": "#e41a1c", "L2": "#ff7f00", "L3": "#ffd92f", "L4": "#4daf4a",
    "L5": "#377eb8", "L6": "#984ea3", "WM": "#666666",
}

sc.settings.verbosity = 1
sc.settings.figdir = "."
plt.rcParams.update({"figure.dpi": 110, "savefig.dpi": 120, "savefig.bbox": "tight"})


def set_seed(seed: int = SEED) -> None:
    random.seed(seed)
    np.random.seed(seed)


def sample_id(adata) -> str:
    return next(iter(adata.uns["spatial"]))


def spatial_plot(adata, color, path, title=None, ncols=4, palette=None, wspace=0.1, **kw):
    """Spot plot over the low-res H&E image, saved to `path`."""
    sid = sample_id(adata)
    sc.pl.spatial(adata, color=color, library_id=sid, img_key="lowres",
                  title=title, ncols=ncols, palette=palette, wspace=wspace, show=False,
                  frameon=False, **kw)
    plt.savefig(path)
    plt.close("all")


def leiden_n_clusters(adata, target, key, neighbors_key=None, lo=0.02, hi=2.5,
                      max_iter=30, seed=SEED):
    """Binary-search the Leiden resolution so that the number of clusters
    equals `target` (or is as close as possible). Returns the resolution used."""
    def run(res):
        sc.tl.leiden(adata, resolution=res, key_added=key, random_state=seed,
                     flavor="igraph", n_iterations=2, directed=False,
                     neighbors_key=neighbors_key)
        return adata.obs[key].nunique()

    best = None  # (|n - target|, res)
    for _ in range(max_iter):
        mid = (lo + hi) / 2
        n = run(mid)
        if best is None or abs(n - target) < best[0]:
            best = (abs(n - target), mid)
        if n == target:
            break
        if n < target:
            lo = mid
        else:
            hi = mid
    res = best[1]
    n = run(res)
    print(f"Leiden ({key}): resolution={res:.4f} -> {n} clusters (target {target})")
    return res


def cluster_and_save(adata, use_rep, method, n_clusters, n_neighbors=15):
    """Shared by CLUSTER_EXPR and CLUSTER_SPATIAL so both methods differ only
    in the representation (`use_rep`) they cluster on.

    Writes {id}_{method}_clusters.csv, {id}_{method}_umap.png,
    {id}_{method}_spatial.png."""
    sid = sample_id(adata)
    set_seed()
    key = f"cluster_{method}"
    sc.pp.neighbors(adata, use_rep=use_rep, n_neighbors=n_neighbors, random_state=SEED)
    sc.tl.umap(adata, random_state=SEED)
    res = leiden_n_clusters(adata, n_clusters, key)

    out = adata.obs[[key]].rename(columns={key: "cluster"})
    out.index.name = "barcode"
    out.to_csv(f"{sid}_{method}_clusters.csv")

    adata.uns["ground_truth_colors"] = [LAYER_COLORS[k] for k in LAYER_ORDER]
    sc.pl.umap(adata, color=[key, "ground_truth"], show=False, frameon=False,
               wspace=0.3, title=[f"{method} clusters (res={res:.2f})", "ground truth"])
    plt.savefig(f"{sid}_{method}_umap.png")
    plt.close("all")
    spatial_plot(adata, [key, "ground_truth"], f"{sid}_{method}_spatial.png",
                 title=[f"{method} clusters", "ground truth"], wspace=0.15)
    return res
