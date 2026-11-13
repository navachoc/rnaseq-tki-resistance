# =============================================================================
# ALK_11_individual_heatmaps.R
# Purpose: 3 individual heatmaps per cell model (ALK1/ALK2/ALK3),
#   top 500 most variable genes of each model, per-gene Z-score in the
#   colour gradient. Per-model alternative to the joint TMM heatmap of the
#   23 samples (Panel A of ALK_09_combined_figure.R), which is not modified
#   (ALK_09 still also generates the PCA of Panel B-D).
#
# With 500 genes per heatmap it makes no sense for the row names to be
# legible (global view instead of gene-by-gene reference, same criterion as
# the WES figure panels script) -- gene names are omitted
# (show_row_names = FALSE) and the image height is fixed instead of
# proportional to the number of genes.
#
# TMM (not LFC) per model, per-gene Z-score computed AFTER selecting
# the top500 by variance (the ranking uses raw TMM, not Z-score, like
# Panel A of ALK_09) -- top500 computed PER MODEL, with its own columns
# only (unlike Panel A of ALK_09, which merges the 3 models into
# a single matrix before ranking).
#
# Columns per model = all the samples of that model (Untreated + Treated
# + resistant lines), same criterion as Panel A of ALK_09 (not only
# resistant).
#
# Does NOT reanalyse anything -- reads the same TMM matrices (01_mat_tmm.csv) as
# ALK_09_combined_figure.R.
#
# Inputs:
#   results/ALK/ALK1/tables/01_mat_tmm.csv
#   results/ALK/ALK2/tables/01_mat_tmm.csv
#   results/ALK/ALK3/tables/01_mat_tmm.csv
#
# Outputs (results/ALK/combined/figures/):
#   ALK_11_fig1_heatmap_ALK1.png
#   ALK_11_fig1_heatmap_ALK2.png
#   ALK_11_fig1_heatmap_ALK3.png
#   ALK_11_fig1_heatmap_combinado.png   (the 3 side by side, quick reference)
#
# Labels, axes and legends in English.
#
# Usage: source(here::here("scripts/specific/ALK/ALK_11_individual_heatmaps.R"))
#   or Rscript from any folder of the repo.
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

OUT_DIR <- here::here("results", "ALK", "combined")
dir.create(file.path(OUT_DIR, "figures"), recursive = TRUE, showWarnings = FALSE)

# Palette shared by the whole figure -- identical to ALK_09_combined_figure.R
# (one gradient colour per technique: black-grey-gold for RNA-seq).
GRADIENTE_RNA    <- c("#000000", "#808080", "#f2a900")
COLOR_CONDICION  <- c(Parental = "#e6a817", Resistant = "#013c58")
COLOR_DRUG_FIRMA <- c(No = "#95A5A6", Lorlatinib = "#00B894", Alectinib = "#F4978E")

TOP_N_VARIANZA <- 500

# FIXED Z-score scale, IDENTICAL in the 3 models (-2 to 2). With a cap
# of its own per model (min(3, max(abs(z)))) each heatmap would have a
# different legend (-2/2, -3/3, -4/4) and they would not be comparable at
# a glance.
Z_MAX        <- 2
COL_FUN_RNA  <- circlize::colorRamp2(c(-Z_MAX, 0, Z_MAX), GRADIENTE_RNA)

# FIXED, identical dimensions for the 3 individual PNGs -- the width does not
# scale with the number of columns of each model (ALK3 has 11 samples and
# ALK1 5), the 3 files have the same size even though the number of samples
# differs (cells are compressed more in ALK3, an expected effect of fixing
# the image size).
ANCHO_HEATMAP  <- 8
ALTO_HEATMAP   <- 10

# =============================================================================
# 1. Read TMM per model
# =============================================================================
cat("1. Reading TMM per model...\n")

leer_tmm <- function(cell_line) {
  ruta <- here::here("results", "ALK", cell_line, "tables", "01_mat_tmm.csv")
  df <- read.csv(ruta, check.names = FALSE)
  rownames(df) <- df$Gene
  as.matrix(df[, -1, drop = FALSE])
}

tmm_alk1 <- leer_tmm("ALK1")
tmm_alk2 <- leer_tmm("ALK2")
tmm_alk3  <- leer_tmm("ALK3")

# Same columns per model as ALK_09_combined_figure.R (Untreated +
# Treated + resistant lines).
cols_alk1  <- c("ALK1", "ALK1_L", "ALK1_L1", "ALK1_L2", "ALK1_L3")
cols_alk2  <- c("ALK2", "ALK2_L", "ALK2_L1", "ALK2_L2", "ALK2_L3", "ALK2_L4", "ALK2_L5")
cols_alk3_l <- c("ALK3", "ALK3_L", "ALK3_L2", "ALK3_L3", "ALK3_L7", "ALK3_L9", "ALK3_L10")
cols_alk3_a <- c("ALK3_A", "ALK3_A2", "ALK3_A3", "ALK3_A4")  # "ALK3" (Untreated) is already in cols_alk3_l
cols_alk3   <- c(cols_alk3_l, cols_alk3_a)

derivar_condicion <- function(col) ifelse(grepl("_[LA][0-9]+$", col), "Resistant", "Parental")
derivar_drug <- function(col) {
  case_when(
    !grepl("_", col) ~ "No",
    grepl("_A", col)  ~ "Alectinib",
    TRUE              ~ "Lorlatinib"
  )
}

genes_validos <- function(mat) {
  g <- rownames(mat)
  g[!grepl("^ENSG", g) & !grepl("\\.", g)]
}

# =============================================================================
# 2. One heatmap per model -- top500 most variable WITHIN that model,
#    per-gene Z-score, without row names (global view).
# =============================================================================
cat("2. Generating individual heatmaps (top", TOP_N_VARIANZA, "most variable per model)...\n")

# mostrar_leyenda = FALSE is used only for the copies that go inside the
# combined panel (section 3), where the legend is single and shared instead
# of one per panel (a single horizontal legend).
construir_heatmap_modelo <- function(modelo, mat, cols, top_n, mostrar_leyenda = TRUE) {
  genes_ok   <- intersect(genes_validos(mat), rownames(mat))
  mat_modelo <- mat[genes_ok, cols, drop = FALSE]

  var_g     <- apply(mat_modelo, 1, var, na.rm = TRUE)
  genes_top <- names(sort(var_g, decreasing = TRUE))[1:min(top_n, length(var_g))]
  mat_top   <- mat_modelo[genes_top, , drop = FALSE]
  mat_z     <- t(scale(t(mat_top)))  # per-gene (row) Z-score, on TMM

  col_condicion <- setNames(derivar_condicion(colnames(mat_z)), colnames(mat_z))
  col_drug      <- setNames(derivar_drug(colnames(mat_z)), colnames(mat_z))

  ann_col <- HeatmapAnnotation(
    Condition = col_condicion,
    Drug      = col_drug,
    col = list(Condition = COLOR_CONDICION, Drug = COLOR_DRUG_FIRMA),
    show_annotation_name = FALSE,
    show_legend = mostrar_leyenda,
    simple_anno_size = unit(3, "mm")
  )

  Heatmap(
    mat_z,
    name = "Z-score",
    col = COL_FUN_RNA,
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    show_row_names = FALSE,
    show_heatmap_legend = mostrar_leyenda,
    top_annotation = ann_col,
    column_names_gp = gpar(fontsize = 9),
    rect_gp = gpar(col = NA),
    heatmap_legend_param = list(title = "Z-score\n(TMM per gene)", at = c(-Z_MAX, 0, Z_MAX)),
    column_title = paste0(modelo, " -- top ", length(genes_top), " most variable genes (TMM, Z-score)")
  )
}

ht_alk1 <- construir_heatmap_modelo("ALK1", tmm_alk1, cols_alk1, TOP_N_VARIANZA)
ht_alk2 <- construir_heatmap_modelo("ALK2", tmm_alk2, cols_alk2, TOP_N_VARIANZA)
ht_alk3  <- construir_heatmap_modelo("ALK3",  tmm_alk3,  cols_alk3,  TOP_N_VARIANZA)

heatmaps <- list(ALK1 = ht_alk1, ALK2 = ht_alk2, ALK3 = ht_alk3)

for (modelo in names(heatmaps)) {
  out <- file.path(OUT_DIR, "figures", paste0("ALK_11_fig1_heatmap_", modelo, ".png"))
  png(out, width = ANCHO_HEATMAP, height = ALTO_HEATMAP, units = "in", res = 300)
  draw(heatmaps[[modelo]], heatmap_legend_side = "right", annotation_legend_side = "right")
  dev.off()
  cat("  ", out, "\n")
}

# =============================================================================
# 3. Combined panel -- the 3 heatmaps side by side, quick reference
#    (grid viewports, not patchwork -- ComplexHeatmap is not a ggplot grob)
# =============================================================================
cat("3. Generating combined panel...\n")

# Copies WITHOUT their own legend (a single horizontal legend) -- the 3
# share Z_MAX/COL_FUN_RNA and the same Condition/Drug colours, so a single
# shared legend is valid for all 3.
ht_alk1_comb <- construir_heatmap_modelo("ALK1", tmm_alk1, cols_alk1, TOP_N_VARIANZA, mostrar_leyenda = FALSE)
ht_alk2_comb <- construir_heatmap_modelo("ALK2", tmm_alk2, cols_alk2, TOP_N_VARIANZA, mostrar_leyenda = FALSE)
ht_alk3_comb  <- construir_heatmap_modelo("ALK3",  tmm_alk3,  cols_alk3,  TOP_N_VARIANZA, mostrar_leyenda = FALSE)

leyenda_compartida <- packLegend(
  Legend(col_fun = COL_FUN_RNA, title = "Z-score (TMM per gene)", direction = "horizontal", legend_width = unit(4, "cm")),
  Legend(labels = names(COLOR_CONDICION), legend_gp = gpar(fill = COLOR_CONDICION), title = "Condition", direction = "horizontal"),
  Legend(labels = names(COLOR_DRUG_FIRMA), legend_gp = gpar(fill = COLOR_DRUG_FIRMA), title = "Drug", direction = "horizontal"),
  direction = "horizontal"
)

out_combinado <- file.path(OUT_DIR, "figures", "ALK_11_fig1_heatmap_combinado.png")

# The 3 panels with the SAME width and height -- widths = 1/1/1 instead of
# proportional to the number of columns of each model. Bottom row reserved
# for the single legend.
png(out_combinado, width = 3 * ANCHO_HEATMAP, height = ALTO_HEATMAP + 1.5, units = "in", res = 300)
grid.newpage()
pushViewport(viewport(layout = grid.layout(2, 3,
  heights = unit(c(ALTO_HEATMAP, 1.5), c("inches", "inches")),
  widths  = unit(c(1, 1, 1), "null"))))

pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
draw(ht_alk1_comb, newpage = FALSE)
upViewport()

pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
draw(ht_alk2_comb, newpage = FALSE)
upViewport()

pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 3))
draw(ht_alk3_comb, newpage = FALSE)
upViewport()

pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1:3))
draw(leyenda_compartida, just = "center")
upViewport()

dev.off()
cat("  ", out_combinado, "\n")

cat("\n=== ALK_11_individual_heatmaps.R completed (3 individual PNGs + 1 combined) ===\n")
