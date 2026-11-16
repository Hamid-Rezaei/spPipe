#!/usr/bin/env python
"""REPORT: combine all samples into one self-contained HTML report.

Reads every {id}_metrics.csv, {id}_qc_summary.csv, {id}_moran_top.csv and the
per-sample figures found in --indir.

Outputs: summary_metrics.csv, summary_ari.png, report.html
"""
import argparse
import base64
import datetime as dt
import io
from pathlib import Path

import sppipe_utils as U  # noqa: I001,F401  (sets matplotlib backend)
import matplotlib.pyplot as plt
import pandas as pd
from jinja2 import Template
from PIL import Image

FIGURES = [  # (file suffix, caption, pipeline step)
    ("comparison_spatial", "Manual layers vs the two clustering methods (clusters coloured by their best-matching layer)", "EVALUATE"),
    ("ari_bar", "ARI / NMI of each method against the manual annotation", "EVALUATE"),
    ("confusion", "Fraction of each manual layer assigned to each matched cluster", "EVALUATE"),
    ("expr_umap", "UMAP of expression-only clusters next to the manual layers", "CLUSTER_EXPR"),
    ("spatial_umap", "UMAP of spatially-aware clusters next to the manual layers", "CLUSTER_SPATIAL"),
    ("markers_spatial", "Known DLPFC layer markers on the tissue", "SPATIAL_GENES"),
    ("svg_spatial", "Top spatially variable genes (Moran's I)", "SPATIAL_GENES"),
    ("domain_markers_dotplot", "Top marker genes of each spatially-aware domain", "SPATIAL_GENES"),
    ("qc_spatial", "QC metrics and manual annotation on the tissue", "QC_FILTER"),
    ("qc_violin", "QC metric distributions with filter thresholds", "QC_FILTER"),
    ("pca_variance", "Variance explained by each principal component", "NORMALIZE"),
]

TEMPLATE = Template("""<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>spPipe report</title>
<style>
 body{font-family:system-ui,-apple-system,Segoe UI,Helvetica,Arial,sans-serif;max-width:1100px;margin:0 auto;padding:16px;color:#222;background:#fff;line-height:1.45}
 h1{margin-bottom:0} .sub{color:#666;margin-top:4px}
 table{border-collapse:collapse;margin:8px 0 16px;font-size:14px}
 th,td{border:1px solid #ddd;padding:4px 10px;text-align:right} th{background:#f4f4f4}
 td:first-child,th:first-child{text-align:left}
 figure{margin:12px 0 24px} img{max-width:100%;height:auto;border:1px solid #eee}
 figcaption{font-size:13px;color:#555} .step{font-family:monospace;background:#eef;padding:1px 5px;border-radius:3px}
 details{margin:6px 0 20px} summary{cursor:pointer;font-weight:600}
 .better{color:#1a7f37;font-weight:600} .worse{color:#b42318;font-weight:600}
</style></head><body>
<h1>spPipe report</h1>
<p class="sub">Does spatial information help recover cortical layers in human DLPFC (10x Visium)? Generated {{ date }}.</p>

<h2>Summary</h2>
<table><tr><th>Sample</th><th>Spots scored</th><th>ARI expression-only</th><th>ARI spatially-aware</th><th>&Delta;ARI</th><th>NMI expression-only</th><th>NMI spatially-aware</th></tr>
{% for r in wide %}<tr><td>{{ r.sample_id }}</td><td>{{ r.n_spots_scored }}</td><td>{{ "%.3f"|format(r.ARI_expr) }}</td><td>{{ "%.3f"|format(r.ARI_spatial) }}</td>
<td class="{{ 'better' if r.dARI > 0 else 'worse' }}">{{ "%+.3f"|format(r.dARI) }}</td><td>{{ "%.3f"|format(r.NMI_expr) }}</td><td>{{ "%.3f"|format(r.NMI_spatial) }}</td></tr>{% endfor %}
</table>
<p>Spatially-aware clustering scored higher ARI than expression-only clustering in <b>{{ n_better }} of {{ wide|length }}</b> sample(s) (mean &Delta;ARI = {{ "%+.3f"|format(mean_d) }}).</p>
<figure><img src="data:image/jpeg;base64,{{ summary_png }}" alt="ARI per sample"><figcaption>ARI per sample and method (<span class="step">REPORT</span>)</figcaption></figure>

<h2>QC</h2>
<table><tr>{% for c in qc.columns %}<th>{{ c }}</th>{% endfor %}</tr>
{% for _, r in qc.iterrows() %}<tr>{% for c in qc.columns %}<td>{{ r[c] }}</td>{% endfor %}</tr>{% endfor %}</table>

{% for s in samples %}
<h2>Sample {{ s.id }}</h2>
<details><summary>Top spatially variable genes (Moran's I)</summary>
<table><tr><th>gene</th><th>Moran's I</th><th>FDR</th></tr>
{% for g, r in s.moran.iterrows() %}<tr><td>{{ g }}</td><td>{{ "%.3f"|format(r.I) }}</td><td>{{ "%.2e"|format(r.pval_norm_fdr_bh) }}</td></tr>{% endfor %}
</table></details>
{% for f in s.figs %}<figure><img src="data:image/jpeg;base64,{{ f.b64 }}" alt="{{ f.caption }}">
<figcaption>{{ f.caption }} (<span class="step">{{ f.step }}</span>, <code>{{ f.file }}</code>)</figcaption></figure>{% endfor %}
{% endfor %}
</body></html>
""")


def b64(path: Path, max_width: int = 1400) -> str:
    """Embed a figure as a downsized JPEG so the report stays a few MB.
    Full-resolution PNGs remain available as separate pipeline outputs."""
    img = Image.open(path).convert("RGB")
    if img.width > max_width:
        img = img.resize((max_width, round(img.height * max_width / img.width)), Image.LANCZOS)
    buf = io.BytesIO()
    img.save(buf, format="JPEG", quality=85, optimize=True)
    return base64.b64encode(buf.getvalue()).decode()


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--indir", default=".")
    args = ap.parse_args()
    indir = Path(args.indir)

    metrics = pd.concat([pd.read_csv(p) for p in sorted(indir.glob("*_metrics.csv"))
                         if p.name != "summary_metrics.csv"], ignore_index=True)
    metrics["sample_id"] = metrics["sample_id"].astype(str)
    metrics = metrics.sort_values(["sample_id", "method"])
    metrics.to_csv("summary_metrics.csv", index=False)
    qc = pd.concat([pd.read_csv(p) for p in sorted(indir.glob("*_qc_summary.csv"))], ignore_index=True)

    wide = metrics.pivot(index="sample_id", columns="method", values=["ARI", "NMI"])
    wide.columns = [f"{a}_{b}" for a, b in wide.columns]
    wide["n_spots_scored"] = metrics.groupby("sample_id")["n_spots_scored"].first()
    wide["dARI"] = wide["ARI_spatial"] - wide["ARI_expr"]
    wide = wide.reset_index()

    fig, ax = plt.subplots(figsize=(max(4, 1.2 * len(wide) + 2), 3.5))
    x = range(len(wide))
    ax.bar([i - 0.2 for i in x], wide["ARI_expr"], 0.4, label="Expression-only", color="#9e9e9e")
    ax.bar([i + 0.2 for i in x], wide["ARI_spatial"], 0.4, label="Spatially-aware", color="#2b8cbe")
    ax.set_xticks(list(x), wide["sample_id"])
    ax.set_ylabel("ARI vs manual layers")
    ax.set_ylim(0, max(0.6, float(wide[["ARI_expr", "ARI_spatial"]].max().max()) + 0.1))
    ax.legend(frameon=False, fontsize=8)
    ax.set_title("Layer recovery per sample")
    fig.savefig("summary_ari.png")
    plt.close(fig)

    samples = []
    for sid in wide["sample_id"]:
        figs = []
        for suffix, caption, step in FIGURES:
            p = indir / f"{sid}_{suffix}.png"
            if p.exists():
                figs.append({"b64": b64(p), "caption": caption, "step": step, "file": p.name})
        moran_p = indir / f"{sid}_moran_top.csv"
        moran = pd.read_csv(moran_p, index_col="gene") if moran_p.exists() else pd.DataFrame()
        samples.append({"id": sid, "figs": figs, "moran": moran})

    html = TEMPLATE.render(
        date=dt.date.today().isoformat(), wide=wide.to_dict("records"),
        n_better=int((wide["dARI"] > 0).sum()), mean_d=float(wide["dARI"].mean()),
        summary_png=b64(Path("summary_ari.png")), qc=qc, samples=samples)
    Path("report.html").write_text(html)
    print(wide.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
