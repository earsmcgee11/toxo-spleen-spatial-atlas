# A single-cell and spatial transcriptomics spleen atlas reveals candidate regulators of CD8 T cell differentiation

Code associated with the manuscript above (under review).

**Authors:** Rachel Coombs, Yash Kulkarni (#), Can Ergen (#), Mari Brady, Leon Han, Eli Benuck, Aaron Streets*, Nir Yosef*, Ellen A. Robey* (# contributed equally).

## Summary (from manuscript)

Authors combine single-cell and spatial transcriptomics to map the cellular organization of the spleen during chronic infection. Their atlas identifies candidate signals regulating CD8 T cell differentiation and reveals a role for IL-27 in sustaining ongoing effector responses.

## Abstract (from manuscript)

The spleen is an important site for maintaining CD8 T cell responses during chronic infection, but a comprehensive analysis of cell types and molecules that regulate CD8 T cell differentiation in the spleen is lacking. In a well-controlled *Toxoplasma gondii* chronic infection model, armed effector T cells are continuously produced in the spleen, making this a useful model to address this gap. We used single cell RNAseq and spatial transcriptomics to provide an atlas of spleen from both infected and uninfected mice. We identified candidate regulators of T cell differentiation, including IL27: a cytokine whose role in chronic *T. gondii* infection was previously uncharacterized. Using CRISPR knockout of IL27R on parasite-specific CD8 T cells, we provided experimental evidence that IL27 promotes an ongoing effector response by sustaining a proliferative intermediate population. Given the relatively unperturbed splenic architecture and cell composition in this infection setting, these data should be broadly useful for studies of spleen biology.

## Repository layout

| Folder | Related figures |
|--------|-----------------|
| `01_visium_clustering/` | Fig. 2, S3 |
| `02_tacco/` | Fig. 2D |
| `03_cell2location/` | Fig. 3, S4 |
| `04_cd8_myeloid_spatial/` | Fig. 3B, S4 |
| `05_cellchat/` | Fig. 4D, Table 1 |
| `06_expression_panels/` | S2, Fig. 2C, S3E |

See `docs/FIGURE_MAP.md`. Large data objects (`.h5ad`, Loupe, raw sequencing) are not included — see `data/README.md`.

## Data availability (from manuscript)

Processed scRNA-seq and Visium data from this study will be deposited in GEO and available through CELLxGENE upon publication. Loupe Browser files for Visium samples will be available at Zenodo upon publication.
