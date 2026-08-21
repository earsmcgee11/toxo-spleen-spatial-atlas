# toxo-spleen-spatial-atlas

Analysis code supporting a single-cell RNA-seq and Visium spatial transcriptomics atlas of mouse spleen during chronic *Toxoplasma gondii* infection (infected and uninfected), used to map CD8 T cell and myeloid niches and nominate regulators of effector differentiation (including IL-27).

## Citation

Manuscript under review (not yet published).

## Abstract (from manuscript)

The spleen is an important site for maintaining CD8 T cell responses during chronic infection, but a comprehensive analysis of cell types and molecules that regulate CD8 T cell differentiation in the spleen is lacking. In a well-controlled *Toxoplasma gondii* chronic infection model, armed effector T cells are continuously produced in the spleen, making this a useful model to address this gap. We used single cell RNAseq and spatial transcriptomics to provide an atlas of spleen from both infected and uninfected mice. We identified candidate regulators of T cell differentiation, including IL27: a cytokine whose role in chronic *T. gondii* infection was previously uncharacterized. Using CRISPR knockout of IL27R on parasite-specific CD8 T cells, we provided experimental evidence that IL27 promotes an ongoing effector response by sustaining a proliferative intermediate population. Given the relatively unperturbed splenic architecture and cell composition in this infection setting, these data should be broadly useful for studies of spleen biology.

## What this repository contains

Computational workflows for the atlas analyses above, including:

- Visium spatial niche clustering and marker panels  
- Integration of scRNA-seq with Visium (TACCO; cell2location deconvolution)  
- CD8–myeloid spatial coupling  
- Expression heatmaps / related figure panels  
- CellChat ligand–receptor analysis (see `05_cellchat/`; forthcoming)

Notebooks may still contain local absolute paths or provisional figure labels while the manuscript is finalized. Large data objects (`.h5ad`, Loupe, raw sequencing) are **not** included here — see `data/README.md`.

## Repository layout

| Folder | Paper focus |
|--------|-------------|
| `01_visium_clustering/` | Visium niches / Leiden clustering (Fig. 2, S3) |
| `02_tacco/` | Transfer of spatial annotations to scRNA-seq (Fig. 2D) |
| `03_cell2location/` | Spatial deconvolution and abundance heatmaps (Fig. 3, S4) |
| `04_cd8_myeloid_spatial/` | CD8–myeloid spatial coupling / Moran-style analyses (Fig. 3B, S4) |
| `05_cellchat/` | Ligand–receptor networks (Fig. 4, Table 1) |
| `06_expression_panels/` | Marker / IFN / heatmaps and dotplots (S2, Fig. 2C, related panels) |

See `docs/FIGURE_MAP.md` for a figure-to-code map.

## Authors

Rachel Coombs, Yash Kulkarni (#), Can Ergen (#), Mari Brady, Leon Han, Eli Benuck, Aaron Streets\*, Nir Yosef\*, Ellen A. Robey\* (# contributed equally).

## Data availability (from manuscript)

Processed scRNA-seq and Visium data from this study will be deposited in GEO and available through CELLxGENE upon publication. Loupe Browser files for Visium samples will be available at Zenodo upon publication.
