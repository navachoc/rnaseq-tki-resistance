# scripts/specific/EGFR

Figures and gene lists for the EGFR design (one parental and one resistant population). The analysis itself runs in `scripts/simple/`; these scripts read its tables and do not redo TMM or fold changes.

| Script | What it does |
|--------|--------------|
| `egfr_figura3_panels.R` | Four panels: heatmap of the top 500 most variable genes, MA plot with labelled key genes, bar plots of the conservative signature (top 15 up and top 15 down), and a KEGG up/down dot plot |
| `signature_amplicon_egfr1.R` | Filters the gene lists of a genomic amplicon (focal and broader region) by expression, using the LFC master table. Writes two CSV files |

## Run order

1. `scripts/simple/rna_01.1` → `rna_02.1` → `rna_03.1` → `rna_04.1` with the EGFR dataset set in the CONFIG block.
2. `egfr_figura3_panels.R` and `signature_amplicon_egfr1.R`, in any order.

## Things to know

- With only two samples, a per-gene z-score is degenerate (every row collapses into two flat colors). The heatmap falls back to raw log2CPM with a sequential palette when the matrix has two columns or fewer.
- `LOGFC_THRESHOLD` in `egfr_figura3_panels.R` must match the value used when running `rna_02.1_lfc.R`.
- The amplicon gene lists come from an exome sequencing analysis and are fixed inside `signature_amplicon_egfr1.R`.
