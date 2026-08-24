# Figure → code map

| Paper figure | Folder | Main code | Notes |
|--------------|--------|-----------|-------|
| Fig. 2 / S3 | `01_visium_clustering` | `visium_harmony_leiden_res1.ipynb` | Harmony + Leiden (`leiden_res1`); niches later stored as `paper_clusters` |
| Fig. 2D | `02_tacco` | `tacco_visium_to_scrna.ipynb` | TACCO parameters from Methods; not a fully executed run in this repo |
| Fig. 3 / S4 | `03_cell2location` | `cell2location_visium_sd.ipynb` | Standard Visium (not HD) |
| Fig. 3B / S4 | `04_cd8_myeloid_spatial` | `cd8_myeloid_spatial_coupling.ipynb` | Pearson (Fig. 3B) + bivariate Moran’s I (Fig. S4D) |
| Fig. 4D / Table 1 | `05_cellchat` | `cellchat_fig4_table1.Rmd`, `fig4d_communication_rewiring.ipynb` | CellChat exports + circle plot |
| S2 / Fig. 2C / S3E | `06_expression_panels` | `expression_marker_heatmaps.ipynb` | Marker / IFN / focus-gene panels |

Not included here: MeRN / metabolic analyses, Visium HD exploratory work, TCR repertoire, or Fig. 5 IL-27 CRISPR / flow cytometry.
