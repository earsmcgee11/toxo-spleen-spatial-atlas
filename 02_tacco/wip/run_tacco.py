import scanpy as sc
import pandas as pd
import numpy as np
import tacco as tc

# ── 1. Load data ─────────────────────────────────────────────────────────────
print("Loading single-cell data...")
sc_adata = sc.read_h5ad("/Users/yashkulkarni/cellxgene_data/sc_spleen_clean.h5ad")
print(f"  SC: {sc_adata.shape}")

print("Loading Visium data...")
st_adata = sc.read_h5ad("/Users/yashkulkarni/cellxgene_data/1_18_25_visium_annotated.h5ad")
print(f"  ST: {st_adata.shape}")

# ── 2. Rename tissue regions to figure 2 names ──────────────────────────────
region_rename = {
    "WP-A: TZ": "TZ",
    "WP-B: BZ": "BZ",
    "WP-D: GC": "2° Fol",
    "WP-C: MZ": "MZ",
    "RP-A: GZMK": "RP-A-PC",
    "RP-B: H_MK": "RP-B-MK",
    "RP-C: NGP": "RP-C-Neut",
    "RP-D: F480": "RP-D",
    "RP-E: Rhag": "RP-E-RBC",
}
st_adata.obs["tissue_regions"] = st_adata.obs["paper_clusters"].map(region_rename)
print(f"\nVisium tissue regions:")
print(st_adata.obs["tissue_regions"].value_counts())

# ── 3. Prep: filter genes, use raw counts ────────────────────────────────────
exclude_prefixes = ("Rps", "Rpl", "mt-", "ERCC")

st_adata = st_adata[:, [g for g in st_adata.var_names if not g.startswith(exclude_prefixes)]].copy()
st_adata.X = st_adata.layers["counts"].copy()

sc_adata = sc_adata[:, [g for g in sc_adata.var_names if not g.startswith(exclude_prefixes)]].copy()

sc.pp.filter_genes(sc_adata, min_cells=5)
print(f"\nAfter filtering: SC={sc_adata.shape}, ST={st_adata.shape}")

# Normalize (TACCO expects normalized counts)
sc.pp.normalize_total(st_adata)
sc.pp.normalize_total(sc_adata)

# ── 4. Split by condition and run TACCO separately ───────────────────────────
results = []

for condition, sc_tp, st_health in [("infected", "3wk", "infected"),
                                     ("uninfected", "0wk", "uninfected")]:
    print(f"\n{'='*60}")
    print(f"  Running TACCO for: {condition}")
    print(f"{'='*60}")

    sc_sub = sc_adata[sc_adata.obs["timepoint"] == sc_tp].copy()
    st_sub = st_adata[st_adata.obs["health"] == st_health].copy()

    print(f"  SC cells: {sc_sub.n_obs}, ST spots: {st_sub.n_obs}")
    print(f"  ST regions: {st_sub.obs['tissue_regions'].value_counts().to_dict()}")

    tc.tl.annotate(
        sc_sub, st_sub,
        "tissue_regions",
        result_key="tissue_regions",
        normalize_to="reference",
        bisections=10,
        bisection_divisor=9,
        assume_valid_counts=True,
        lamb=0,
        multi_center=None,
        min_log2foldchange=0,
        min_genes_per_cell=1,
    )

    result_df = sc_sub.obsm["tissue_regions"].copy()
    results.append(result_df)
    print(f"  Done. Mean scores:")
    print(result_df.mean().to_string())

# ── 5. Combine and export ────────────────────────────────────────────────────
tissue_df = pd.concat(results)
tissue_df = tissue_df.loc[sc_adata.obs_names]  # reorder to original cell order
tissue_df.index.name = "cell_barcode"

out_path = "/Users/yashkulkarni/cellxgene_data/tissue_region_annotations.csv"
tissue_df.to_csv(out_path)
print(f"\nSaved {tissue_df.shape[0]} cells x {tissue_df.shape[1]} regions to: {out_path}")
print(f"\nFirst 5 rows:")
print(tissue_df.head().to_string())

# ── 6. Validation plot: UMAP colored by tissue region scores ─────────────────
import matplotlib.pyplot as plt

sc_full = sc.read_h5ad("/Users/yashkulkarni/cellxgene_data/sc_spleen_clean.h5ad")
for col in tissue_df.columns:
    sc_full.obs[col] = tissue_df[col].values

fig, axes = plt.subplots(3, 3, figsize=(15, 15))
axes = axes.flatten()

region_order = ["BZ", "2° Fol", "RP-B-MK", "MZ", "RP-D", "RP-A-PC",
                "RP-C-Neut", "RP-E-RBC", "TZ"]

for i, region in enumerate(region_order):
    sc.pl.embedding(
        sc_full,
        basis="X_umap_scanvi",
        color=region,
        ax=axes[i],
        show=False,
        frameon=False,
        alpha=0.9,
        cmap="YlGnBu",
        title=region,
        colorbar_loc="right",
        vmin=0,
        vmax=1,
    )

plt.tight_layout()
plt.savefig("/Users/yashkulkarni/cellxgene_data/tissue_regions_umap_validation.png",
            dpi=200, bbox_inches="tight")
plt.savefig("/Users/yashkulkarni/cellxgene_data/tissue_regions_umap_validation.pdf",
            bbox_inches="tight")
print("\nSaved validation plot to tissue_regions_umap_validation.png/pdf")
