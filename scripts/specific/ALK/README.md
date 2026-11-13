# scripts/specific/ALK

ALK-adapted version of the pipeline. There is one folder per cell line, each with numbered steps, plus combined scripts at the root of this folder that draw figures across all lines. The per-line scripts started as copies of `scripts/simple/` and then diverged to fit this design (several resistant populations per line, one line with two drugs), so they are kept separate.

## Per-line steps (`ALK1/`, `ALK2/`, `ALK3/`)

| Step | Script | What it does |
|------|--------|--------------|
| 01 | `ALK_01_<line>_qc.R` | TMM normalization, PCA, sample distances, top-variable heatmap |
| 02 | `ALK_02_<line>_lfc.R` | Fold changes of each resistant population and gene classification |
| 03 | `ALK_03_<line>_key_genes.R` | Resistance signature (full and conservative) and heatmaps |
| 04 | `ALK_04_<line>_pathways.R` | GO and KEGG over-representation, GSEA, STRING network |
| 05 | `ALK_05_<line>_between_lines.R` | Compares the resistant populations with each other |
| 08.1 | `ALK_08.1_<line>_fusions.R` | STAR-Fusion calls per sample |

The line with two drugs processes both drugs together in each script.

## Combined scripts

| Script | What it does |
|--------|--------------|
| `ALK_09_figura4_combinada.R` | PCA per line, and a joint heatmap of the most variable genes |
| `ALK_09b_kegg_vs_untreated_exploratorio.R` | Exploratory. KEGG and GO with three alternative Up/Down definitions, to compare which pathways each one gives. Writes tables that `ALK_09c` reads |
| `ALK_09c_figura2_kegg_heatmap.R` | LFC heatmap across lines and KEGG dot plots |
| `ALK_10_figura3_fusiones.R` | Fusion figure, built from the tables of the `08.1` scripts |
| `ALK_11_figura1_heatmaps_individuales.R` | One heatmap per line (top 500 most variable genes, z-score) |
| `ALK_12_figura2_lfc_heatmaps.R` | One LFC heatmap per line against the untreated parental, plus a combined panel |

## Run order

1. Steps `01` → `05` and `08.1` for each line.
2. `ALK_09b`, then `ALK_09c`.
3. `ALK_09`, `ALK_10`, `ALK_11` and `ALK_12` read the tables written by steps 01, 02 and 08.1.

## Things to know

- `ALK_09c` also needs `results/ALK/combined/tables/ALK_09c_kegg_top5_final.csv`, a curated ranking of KEGG pathways. It is supplied by hand and no script in this repository generates it, so this figure cannot be rebuilt from scratch.
- The gene lists forced into the heatmaps of `ALK_09c` and `ALK_12` are curated by hand and defined at the top of each script.
- Output file names keep the numbering of the figures they were first made for (`fig1_`, `fig2_`, `fig3_`).
