#!/usr/bin/env python3
"""
Coarse infected Myeloid / CD8 on meRN metabolic embedding: Leiden at chosen resolutions,
then heatmap of fine types (paper_names) × cluster.

  - Myeloid: default resolution 0.675
  - CD8: default resolution 1.0, excluding CD8-Tcell_terminal

Heatmap values are **row-normalized** (each cluster row sums to 1): fraction of cells in that
cluster from each fine type. Raw counts + proportions also written as CSV.

CLI:
  python coarse_leiden_celltype_heatmaps.py

Jupyter:
  from coarse_leiden_celltype_heatmaps import run_heatmaps
  run_heatmaps()
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns

DEFAULT_MERN = "/Users/yashkulkarni/Downloads/spleen_mern.h5ad"
REP = "X_mern_metabolic_rep0"
COARSE = "coarse_celltypes"
PAPER = "paper_names"
N_NEIGHBORS = 15
RANDOM_STATE = 0

RES_MYELOID_DEFAULT = 0.675
RES_CD8_DEFAULT = 1.0
CD8_EXCLUDE_TERMINAL = ["CD8-Tcell_terminal"]


def load_mern_infected(path: str) -> sc.AnnData:
    ad = sc.read_h5ad(path)
    h = ad.obs["health"].astype(str).replace({"1": "infected", "2": "infected"})
    ad.obs["health"] = h
    return ad[ad.obs["health"] == "infected"].copy()


def _argv_after_script() -> list[str]:
    av = sys.argv[1:]
    if not av:
        return av
    if "ipykernel_launcher" in sys.argv[0]:
        return []
    return av


def _leiden_on_subset(
    sub: sc.AnnData, resolution: float, key_added: str
) -> sc.AnnData:
    if REP not in sub.obsm:
        raise SystemExit(f"Missing {REP}; have {list(sub.obsm.keys())}")
    ad = sub.copy()
    sc.pp.neighbors(ad, use_rep=REP, n_neighbors=N_NEIGHBORS, random_state=RANDOM_STATE)
    sc.tl.leiden(
        ad,
        resolution=float(resolution),
        key_added=key_added,
        random_state=RANDOM_STATE,
    )
    return ad


def _cluster_label_order(index: pd.Index) -> list[str]:
    labs = [str(x) for x in index.unique()]

    def sort_key(x: str):
        try:
            return (0, int(x))
        except ValueError:
            return (1, x)

    return sorted(labs, key=sort_key)


def _heatmap_from_crosstab(
    ct: pd.DataFrame,
    title: str,
    outpath: Path,
    dpi: int,
) -> None:
    # rows = cluster, cols = fine type; row-normalize
    row_sum = ct.sum(axis=1).replace(0, np.nan)
    prop = ct.div(row_sum, axis=0).fillna(0.0)

    # consistent ordering
    prop = prop.reindex(index=_cluster_label_order(prop.index))
    prop = prop.reindex(columns=sorted(prop.columns, key=str))

    n_r, n_c = prop.shape
    w = min(22, max(8, 0.35 * n_c + 4))
    h = min(16, max(5, 0.35 * n_r + 3))
    fig, ax = plt.subplots(figsize=(w, h))
    sns.heatmap(
        prop,
        ax=ax,
        cmap="mako",
        vmin=0,
        vmax=1,
        linewidths=0.5,
        linecolor="#e5e5e5",
        cbar_kws={"label": "Fraction of cluster"},
        annot=n_r * n_c <= 200,
        fmt=".2f" if n_r * n_c <= 200 else ".1f",
        annot_kws={"size": 7},
    )
    ax.set_xlabel(PAPER)
    ax.set_ylabel("Leiden cluster (metabolic embedding)")
    ax.set_title(title, fontsize=11)
    plt.xticks(rotation=45, ha="right")
    plt.yticks(rotation=0)
    plt.tight_layout()
    fig.savefig(outpath, dpi=dpi, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def run_heatmaps(
    mern: str = DEFAULT_MERN,
    res_myeloid: float = RES_MYELOID_DEFAULT,
    res_cd8: float = RES_CD8_DEFAULT,
    outdir: str | Path = "coarse_leiden_celltype_figs",
    dpi: int = 200,
) -> None:
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)

    m = load_mern_infected(mern)
    if PAPER not in m.obs.columns:
        raise SystemExit(f"Missing obs['{PAPER}']")

    # --- Myeloid ---
    mye = m[m.obs[COARSE].astype(str) == "Myeloid"].copy()
    if mye.n_obs < 50:
        raise SystemExit(f"Too few myeloid: {mye.n_obs}")
    key_m = "leiden_coarse_mye_metab"
    mye = _leiden_on_subset(mye, res_myeloid, key_m)
    ct_m = pd.crosstab(mye.obs[key_m].astype(str), mye.obs[PAPER].astype(str))
    ct_m.to_csv(out / f"crosstab_counts_myeloid_r{res_myeloid}.csv")
    ct_m.div(ct_m.sum(axis=1).replace(0, np.nan), axis=0).fillna(0).to_csv(
        out / f"crosstab_prop_within_cluster_myeloid_r{res_myeloid}.csv"
    )
    _heatmap_from_crosstab(
        ct_m,
        title=(
            f"Infected coarse Myeloid  |  Leiden on {REP}  |  res={res_myeloid}\n"
            "Rows: clusters; colors: fraction of cluster from each fine type"
        ),
        outpath=out / f"heatmap_myeloid_leiden_r{res_myeloid}.png",
        dpi=dpi,
    )
    print(f"Myeloid n={mye.n_obs}  K={mye.obs[key_m].nunique()}  -> heatmap + CSV")

    # --- CD8 (no terminal) ---
    cd8 = m[m.obs[COARSE].astype(str) == "CD8-Tcell"].copy()
    cd8 = cd8[~cd8.obs[PAPER].isin(CD8_EXCLUDE_TERMINAL)].copy()
    if cd8.n_obs < 50:
        raise SystemExit(f"Too few CD8 after dropping terminal: {cd8.n_obs}")
    key_c = "leiden_coarse_cd8_metab"
    cd8 = _leiden_on_subset(cd8, res_cd8, key_c)
    ct_c = pd.crosstab(cd8.obs[key_c].astype(str), cd8.obs[PAPER].astype(str))
    ct_c.to_csv(out / f"crosstab_counts_cd8_no_terminal_r{res_cd8}.csv")
    ct_c.div(ct_c.sum(axis=1).replace(0, np.nan), axis=0).fillna(0).to_csv(
        out / f"crosstab_prop_within_cluster_cd8_no_terminal_r{res_cd8}.csv"
    )
    _heatmap_from_crosstab(
        ct_c,
        title=(
            f"Infected coarse CD8-Tcell (excl. terminal)  |  Leiden on {REP}  |  res={res_cd8}\n"
            "Rows: clusters; colors: fraction of cluster from each fine type"
        ),
        outpath=out / f"heatmap_cd8_no_terminal_leiden_r{res_cd8}.png",
        dpi=dpi,
    )
    print(f"CD8 (no terminal) n={cd8.n_obs}  K={cd8.obs[key_c].nunique()}  -> heatmap + CSV")
    print("Wrote:", out.resolve())


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, allow_abbrev=False)
    ap.add_argument("--mern", default=DEFAULT_MERN)
    ap.add_argument("--res-myeloid", type=float, default=RES_MYELOID_DEFAULT)
    ap.add_argument("--res-cd8", type=float, default=RES_CD8_DEFAULT)
    ap.add_argument("--outdir", default="coarse_leiden_celltype_figs")
    ap.add_argument("--dpi", type=int, default=200)
    args, _ = ap.parse_known_args(_argv_after_script())
    run_heatmaps(
        mern=args.mern,
        res_myeloid=args.res_myeloid,
        res_cd8=args.res_cd8,
        outdir=args.outdir,
        dpi=args.dpi,
    )


if __name__ == "__main__":
    main()
