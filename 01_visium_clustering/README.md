# Visium clustering / niches (Fig. 2, S3)

## Main notebook

- [`visium_harmony_leiden_res1.ipynb`](visium_harmony_leiden_res1.ipynb) — Visium QC, Harmony batch correction, Leiden niches at resolution 1.0 (`leiden_res1`), spatial maps, and marker genes.

`leiden_res1` is the clustering closest to the final niche labels stored downstream as `paper_clusters` (after biological renaming for the figures).

Set the data path before running (or edit the first code cell):

```bash
export VISIUM_FILTERED_H5AD=/path/to/filtered_spot_adata.h5ad
export VISIUM_CLUSTER_OUT=outputs
```

## Additional notebooks

| File | Role |
|------|------|
| `wip/1_16_24_newvisium cluster.ipynb` | Earlier exploratory clustering notebook |
| `wip/abundance_visium_clustering_12226.ipynb` | Abundance-based clustering vs `paper_clusters` (related to Fig. 3) |
| `wip/visium_de_071425.ipynb` | Downstream differential expression |
| `wip/mz_marker_genes_0726.ipynb` | Marginal-zone marker helper |
