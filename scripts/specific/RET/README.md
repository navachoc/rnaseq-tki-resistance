# scripts/specific/RET

Analysis of the RET design: each model has an untreated parental population, a parental population exposed to the drug, and several resistant populations selected with one of two drugs. The generic pipeline in `scripts/simple/` handles two groups, so this folder holds the "treated parental" variants (`b` scripts) and the figures built on top of them.

## Pipeline variants

| Script | What it does |
|--------|--------------|
| `rna_02b.1_lfc.R` | Fold changes per drug and for all resistant populations together. Four explicit comparisons instead of the five-class scheme |
| `rna_03b.1_key_genes.R` | Resistance signature from the untreated vs resistant comparison. A `Resistance_specific` flag marks genes that also pass the threshold in treated vs resistant with the same sign, which separates resistance from the acute drug response |
| `rna_03c.1_wes_panel_combined.R` | Heatmaps and strip plots of a curated gene panel, embedded in the script |
| `rna_04b.1_pathways.R` | GO and KEGG over-representation and networks for the signature. No GSEA, because no single ranking metric represents the four comparisons |
| `rna_08.1_fusions_ret.R` | Fusion calls (STAR-Fusion) per sample: driver trajectory, most striking fusions per model, fusions shared between models |

## Figures

| Script | What it does |
|--------|--------------|
| `RET_08_figura1_heatmaps_individuales.R` | PCA and heatmap (top 500 most variable genes, z-score) per model, plus a combined panel |
| `RET_08b_figura25_pca.R` | Same computation with a PCA with and without sample names, and the combined figure |
| `RET_09_figura2_lfc_heatmaps.R` | LFC heatmaps per model, GO/KEGG dot plot and a dot plot of fusions |

## Run order

1. `scripts/simple/rna_01.1_qc.R`, once per model.
2. `rna_02b.1` → `rna_03b.1` → `rna_04b.1`, once per model (set `CELL_LINE` in the CONFIG block).
3. `rna_03c.1` and `rna_08.1`.
4. `RET_08`, `RET_08b` and `RET_09`. `RET_09` reads the tables written by `rna_04b.1` and `rna_08.1`.

## Things to know

- `RET_08` and `RET_08b` share the same computation and differ only in how the PCA is drawn.
- The gene lists highlighted in `RET_09` (bypass genes, selected pathways, fusions) are curated by hand and defined at the top of the script.
