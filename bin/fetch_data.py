#!/usr/bin/env python
"""FETCH_DATA: download one DLPFC Visium sample and build an AnnData object.

Inputs (downloaded, or read from --local_dir if given):
  {id}_filtered_feature_bc_matrix.h5   counts matrix (spots x genes)
  {id}_tissue_lowres_image.png         H&E image
  tissue_positions_list.txt            spot coordinates (no header)
  scalefactors_json.json               pixel scale factors
  cluster_labels_{id}.csv              manual layer annotation (ground_truth)

Output:
  {id}_raw.h5ad with obsm['spatial'], uns['spatial'][id] (image + scale
  factors) and obs['ground_truth'] (L1..L6, WM; NaN if unlabeled).
"""
import argparse
import json
import shutil
import time
import urllib.request
from pathlib import Path

import sppipe_utils  # noqa: F401,I001  (sets writable matplotlib/numba cache dirs)
import matplotlib.image as mpimg
import numpy as np
import pandas as pd
import scanpy as sc

H5_BASE = "https://spatial-dlpfc.s3.us-east-2.amazonaws.com/h5"
IMG_BASE = "https://spatial-dlpfc.s3.us-east-2.amazonaws.com/images"
SPOT_BASE = "https://raw.githubusercontent.com/LieberInstitute/HumanPilot/master/10X"
LABEL_BASE = ("https://raw.githubusercontent.com/LieberInstitute/HumanPilot/"
              "master/outputs/SpatialDE_clustering")

LAYER_NAMES = {f"Layer_{i}": f"L{i}" for i in range(1, 7)}
LAYER_NAMES["WM"] = "WM"
LAYER_ORDER = ["L1", "L2", "L3", "L4", "L5", "L6", "WM"]


def download(url: str, dest: Path, retries: int = 3) -> None:
    for attempt in range(1, retries + 1):
        try:
            print(f"Downloading {url}", flush=True)
            with urllib.request.urlopen(url, timeout=120) as r, open(dest, "wb") as f:
                shutil.copyfileobj(r, f)
            return
        except Exception as e:  # network errors: retry with backoff
            print(f"  attempt {attempt} failed: {e}", flush=True)
            if attempt == retries:
                raise
            time.sleep(5 * attempt)


def get_files(sid: str, workdir: Path, local_dir: Path | None) -> dict:
    files = {
        "h5": (f"{sid}_filtered_feature_bc_matrix.h5", f"{H5_BASE}/{sid}_filtered_feature_bc_matrix.h5"),
        "img": (f"{sid}_tissue_lowres_image.png", f"{IMG_BASE}/{sid}_tissue_lowres_image.png"),
        "pos": ("tissue_positions_list.txt", f"{SPOT_BASE}/{sid}/tissue_positions_list.txt"),
        "sf": ("scalefactors_json.json", f"{SPOT_BASE}/{sid}/scalefactors_json.json"),
        "lab": (f"cluster_labels_{sid}.csv", f"{LABEL_BASE}/cluster_labels_{sid}.csv"),
    }
    out = {}
    for key, (name, url) in files.items():
        dest = workdir / name
        if local_dir is not None and (local_dir / name).exists():
            shutil.copy(local_dir / name, dest)
        else:
            download(url, dest)
        out[key] = dest
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sample_id", required=True)
    ap.add_argument("--outdir", default=".")
    ap.add_argument("--local_dir", default=None, help="use files from here instead of downloading")
    args = ap.parse_args()

    sid = args.sample_id
    outdir = Path(args.outdir)
    workdir = outdir / f"{sid}_download"
    workdir.mkdir(parents=True, exist_ok=True)
    f = get_files(sid, workdir, Path(args.local_dir) if args.local_dir else None)

    # Expression matrix
    adata = sc.read_10x_h5(f["h5"])
    adata.var_names_make_unique()

    # Spot coordinates (no header in this file version)
    pos = pd.read_csv(f["pos"], header=None,
                      names=["barcode", "in_tissue", "array_row", "array_col",
                             "pxl_row_in_fullres", "pxl_col_in_fullres"]).set_index("barcode")
    missing = adata.obs_names.difference(pos.index)
    if len(missing):
        raise ValueError(f"{len(missing)} barcodes have no coordinates")
    pos = pos.loc[adata.obs_names]
    adata.obs[["in_tissue", "array_row", "array_col"]] = pos[["in_tissue", "array_row", "array_col"]]
    adata = adata[adata.obs["in_tissue"] == 1].copy()
    pos = pos.loc[adata.obs_names]
    # scanpy convention: obsm['spatial'] = (x, y) = (pxl_col, pxl_row) in full-res pixels
    adata.obsm["spatial"] = pos[["pxl_col_in_fullres", "pxl_row_in_fullres"]].to_numpy(dtype=float)

    # Image + scale factors (layout expected by sc.pl.spatial / squidpy)
    with open(f["sf"]) as fh:
        scalefactors = json.load(fh)
    adata.uns["spatial"] = {sid: {
        "images": {"lowres": mpimg.imread(f["img"])},
        "scalefactors": scalefactors,
        "metadata": {"source": "spatialLIBD / Maynard et al. 2021"},
    }}

    # Manual layer annotation
    lab = pd.read_csv(f["lab"])
    lab["barcode"] = lab["key"].str.replace(f"{sid}_", "", n=1, regex=False)
    truth = lab.set_index("barcode")["ground_truth"].map(LAYER_NAMES)
    adata.obs["ground_truth"] = pd.Categorical(
        truth.reindex(adata.obs_names).to_numpy(), categories=LAYER_ORDER)
    adata.obs["sample_id"] = sid

    n_lab = int(adata.obs["ground_truth"].notna().sum())
    print(f"{sid}: {adata.n_obs} spots, {adata.n_vars} genes, "
          f"{n_lab} labeled, {adata.n_obs - n_lab} unlabeled")
    print(adata.obs["ground_truth"].value_counts(dropna=False).to_string())

    adata.write_h5ad(outdir / f"{sid}_raw.h5ad", compression="gzip")
    shutil.rmtree(workdir)


if __name__ == "__main__":
    main()
