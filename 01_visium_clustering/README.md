# Visium clustering / niches (Fig. 2, S3)

[`visium_harmony_leiden_res1.ipynb`](visium_harmony_leiden_res1.ipynb) — Visium QC, Harmony batch correction, Leiden niches at resolution 1.0 (`leiden_res1`), spatial maps, and marker genes.

`leiden_res1` is the clustering closest to the final niche labels stored downstream as `paper_clusters` (after biological renaming for the figures).

```bash
export VISIUM_FILTERED_H5AD=/path/to/filtered_spot_adata.h5ad
export VISIUM_CLUSTER_OUT=outputs
```
