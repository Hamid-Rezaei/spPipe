#!/usr/bin/env python
"""SPATIAL_GENES: genes that explain the spatial domains.

1. Moran's I (squidpy) on the highly variable genes (MT- genes excluded): how strongly each gene's
   expression is spatially autocorrelated (I near 1 = smooth spatial pattern,
   near 0 = random). Top --n_top genes are reported and mapped.
2. Known DLPFC layer markers (Maynard et al. 2021) mapped on the tissue.
3. Marker genes of each spatially-aware domain (Wilcoxon test, one domain vs rest).

Outputs: {id}_moran_top.csv, {id}_svg_spatial.png, {id}_markers_spatial.png,
         {id}_domain_markers.csv, {id}_domain_markers_dotplot.png
"""
import argparse

import sppipe_utils as U  # noqa: I001
import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
import squidpy as sq

# Layer markers reported in Maynard et al. 2021 (Fig. 2-3)
KNOWN_MARKERS = ["MOBP", "PCP4", "SNAP25", "HPCAL1", "AQP4", "CCK", "FABP7", "KRT17"]


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--input", required=True, help="{id}_norm.h5ad")
    ap.add_argument("--spatial", required=True, help="{id}_spatial_clusters.csv")
    ap.add_argument("--n_top", type=int, default=20)
    args = ap.parse_args()
    U.set_seed()

    adata = sc.read_h5ad(args.input)
    sid = U.sample_id(adata)

    # 1. Moran's I on HVGs (analytic p-values, no permutations -> deterministic)
    sq.gr.spatial_neighbors(adata, coord_type="grid", n_neighs=6, n_rings=1)
    # Mitochondrial genes are spatially patterned for technical reasons
    # (tissue quality), not cortical layering, so they are not ranked.
    hvgs = [g for g in adata.var_names[adata.var["highly_variable"]] if not g.startswith("MT-")]
    sq.gr.spatial_autocorr(adata, mode="moran", genes=hvgs, n_perms=None, n_jobs=1)
    moran = adata.uns["moranI"].sort_values("I", ascending=False)
    top = moran.head(args.n_top)[["I", "pval_norm", "pval_norm_fdr_bh"]]
    top.index.name = "gene"
    top.round(6).to_csv(f"{sid}_moran_top.csv")
    print(f"Top spatially variable genes:\n{top.head(10).round(3).to_string()}")

    U.spatial_plot(adata, top.index[:8].tolist(), f"{sid}_svg_spatial.png",
                   title=[f"{g} (I={i:.2f})" for g, i in top["I"].head(8).items()],
                   cmap="magma", ncols=4)

    # 2. Known layer markers
    markers = [g for g in KNOWN_MARKERS if g in adata.var_names]
    U.spatial_plot(adata, markers, f"{sid}_markers_spatial.png", cmap="magma", ncols=4)

    # 3. Marker genes of each spatial domain
    cl = pd.read_csv(args.spatial, index_col="barcode")["cluster"].astype(str)
    adata.obs["domain"] = pd.Categorical(cl.reindex(adata.obs_names))
    sc.tl.rank_genes_groups(adata, "domain", method="wilcoxon", use_raw=False)
    sc.tl.dendrogram(adata, groupby="domain", use_rep="X_pca")
    df = sc.get.rank_genes_groups_df(adata, group=None)
    df = df.groupby("group", observed=True).head(10)
    df.round(6).to_csv(f"{sid}_domain_markers.csv", index=False)

    sc.pl.rank_genes_groups_dotplot(adata, n_genes=4, show=False, standard_scale="var",
                                    title=f"{sid}: top markers per spatial domain")
    plt.savefig(f"{sid}_domain_markers_dotplot.png")
    plt.close("all")


if __name__ == "__main__":
    main()
