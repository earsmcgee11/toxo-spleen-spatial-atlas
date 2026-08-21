# Figure → code map

Some exploratory notebooks still use provisional figure letters; prefer the manuscript figure numbering below when reading the main notebooks.

| Paper figure | Biology | Folder | Main code | Notes |
|--------------|---------|--------|-----------|-------|
| Fig. 2 / S3 | Visium Leiden niches, cluster markers | `01_visium_clustering` | `visium_harmony_leiden_res1.ipynb` (`leiden_res1`) | Niches later stored as `paper_clusters` |
| Fig. 2D | TACCO spatial→sc scores | `02_tacco` | `wip/run_tacco.py` | Notebook version forthcoming |
| Fig. 3 / S4 | cell2location abundances | `03_cell2location` | notebooks in `wip/` | Standard Visium (not HD) |
| Fig. 3B / S4 | CD8–myeloid spatial coupling | `04_cd8_myeloid_spatial` | notebooks in `wip/` | Bivariate Moran’s I on cell2location abundances |
| Fig. 4 / Table 1 | CellChat L–R networks | `05_cellchat` | forthcoming | |
| S2 / Fig. 2C / related | Marker & IFN expression panels | `06_expression_panels` | notebooks in `wip/` | Retag labels to final figure order as needed |

## Out of scope for this repository

- MeRN / metabolic analyses  
- Visium HD exploratory work  
- TCR repertoire  
- Fig. 5 IL-27 CRISPR / flow cytometry (wet-lab)
