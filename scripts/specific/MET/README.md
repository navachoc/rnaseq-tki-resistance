# scripts/specific/MET

Figures for the MET design (parental vs resistant populations). The analysis itself runs in `scripts/simple/`; this folder only draws figures from its tables.

| Script | What it does |
|--------|--------------|
| `MET_01_figura14.R` | PCA of the top 500 most variable genes (with and without sample names), MA plot of resistant vs parental, and GO BP and KEGG dot plots of the resistance signature |

## Run order

1. `scripts/simple/rna_01.1` → `rna_02.1` → `rna_03.1` → `rna_04.1` with the MET dataset set in the CONFIG block.
2. Optional: `scripts/utilidades/met_egfr_offon_genes_regenerar.R`, which extends the signature and redraws figures and tables.
3. `MET_01_figura14.R`.

## Things to know

- `MET_01` recomputes the TMM matrix for the PCA with the same filter and parameters as `rna_01.1`, and reads the LFC and enrichment tables from `results/MET/`.
- `met_egfr_offon_genes_regenerar.R` writes files with the same names as the pipeline (`04_*`, `01_pca_top500.png`, `02_lfc_scatter.png`). Run it after the pipeline and before `MET_01`.
- Label placement in the PCA and MA plot uses a fixed `ggrepel` seed, so the PNGs are reproducible.
- `MET_01` does not drop Ensembl-ID rows or gene names containing a dot before the PCA, while `rna_01.1` and `scripts/utilidades/clean_figures.R` do. The explained-variance percentages of the PCA therefore differ slightly between them.
