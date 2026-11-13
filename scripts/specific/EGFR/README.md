# scripts/specific/EGFR

Summary figures for the EGFR design (one parental and one resistant population). The analysis itself runs in `scripts/simple/`; these scripts read its tables and do not redo TMM or fold changes.

| Script | What it does |
|--------|--------------|
| `egfr_summary_panels.R` | Four panels: heatmap of the top 500 most variable genes, MA plot with labelled key genes, bar plots of the conservative signature (top 15 up and top 15 down), and a KEGG up/down dot plot |

## Run order

1. `scripts/simple/rna_01.1` → `rna_02.1` → `rna_03.1` → `rna_04.1` with the EGFR dataset set in the CONFIG block.
2. `egfr_summary_panels.R`.

## Things to know

- With only two samples, a per-gene z-score is degenerate (every row collapses into two flat colors). The heatmap falls back to raw log2CPM with a sequential palette when the matrix has two columns or fewer.
- `LOGFC_THRESHOLD` in `egfr_summary_panels.R` must match the value used when running `rna_02.1_lfc.R`.
