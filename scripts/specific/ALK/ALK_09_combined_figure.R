# =============================================================================
# ALK_09_combined_figure.R
# Purpose: combined bulk RNA-seq figure of the 3 ALK models.
#   Panel A: joint clustering heatmap, top 50 most variable genes on the
#   joint TMM matrix ALK1+ALK2+ALK3, per-gene Z-score.
#   Panel B-D: PCA per cell line, top 500 most variable genes of each one on
#   TMM with per-gene Z-score, with and without names, stacked with patchwork.
#   The combined KEGG + LFC heatmap is in ALK_09c_kegg_heatmap.R and
#   the fusions in ALK_10_fusions_figure.R.
# Inputs: results/ALK/{ALK1,ALK2,ALK3}/tables/01_mat_tmm.csv
# Outputs (results/ALK/combined/figures/):
#   ALK_09_fig1_panelA_heatmap_variables.png,
#   ALK_09_fig1_panelBCD_pca_lineas[_sin_nombres].png,
#   ALK_09_fig1_pca_{ALK1,ALK2,ALK3}[_sin_nombres].png
# Usage: source(here::here("scripts/specific/ALK/ALK_09_combined_figure.R"))
#   or Rscript from any folder of the repo.
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(ggrepel)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

OUT_DIR <- here::here("results", "ALK", "combined")
dir.create(file.path(OUT_DIR, "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT_DIR, "tables"),  recursive = TRUE, showWarnings = FALSE)

# Palette shared by the whole figure (one gradient colour per technique and
# categorical colours unified with the rest of the figures).
GRADIENTE_RNA    <- c("#000000", "#808080", "#f2a900")  # black-grey-gold
COLOR_CONDICION  <- c(Parental = "#e6a817", Resistant = "#013c58")
COLOR_DRUG_FIRMA <- c(No = "#95A5A6", Lorlatinib = "#00B894", Alectinib = "#F4978E")

# =============================================================================
# COMMON INFRASTRUCTURE -- TMM and per-sample LFC (vs Untreated of its own
# model/drug). Used by Panel A and by the PCAs.
# =============================================================================
cat("1. Reading TMM and computing LFC per individual sample (vs Untreated)...\n")

leer_tmm <- function(cell_line) {
  ruta <- here::here("results", "ALK", cell_line, "tables", "01_mat_tmm.csv")
  df <- read.csv(ruta, check.names = FALSE)
  rownames(df) <- df$Gene
  as.matrix(df[, -1, drop = FALSE])
}

tmm_alk1 <- leer_tmm("ALK1")
tmm_alk2 <- leer_tmm("ALK2")
tmm_alk3  <- leer_tmm("ALK3")

common_genes_full <- Reduce(intersect, list(rownames(tmm_alk1), rownames(tmm_alk2), rownames(tmm_alk3)))
common_genes_full <- common_genes_full[!grepl("^ENSG", common_genes_full) & !grepl("\\.", common_genes_full)]

lfc_vs_untreated <- function(mat, genes, untreated_col, target_cols) {
  mat[genes, target_cols, drop = FALSE] - mat[genes, untreated_col]
}

cols_alk1  <- c("ALK1", "ALK1_L", "ALK1_L1", "ALK1_L2", "ALK1_L3")
cols_alk2  <- c("ALK2", "ALK2_L", "ALK2_L1", "ALK2_L2", "ALK2_L3", "ALK2_L4", "ALK2_L5")
cols_alk3_l <- c("ALK3", "ALK3_L", "ALK3_L2", "ALK3_L3", "ALK3_L7", "ALK3_L9", "ALK3_L10")
cols_alk3_a <- c("ALK3_A", "ALK3_A2", "ALK3_A3", "ALK3_A4")  # "ALK3" (Untreated) is already in cols_alk3_l

full_alk1  <- lfc_vs_untreated(tmm_alk1, common_genes_full, "ALK1", cols_alk1)
full_alk2  <- lfc_vs_untreated(tmm_alk2, common_genes_full, "ALK2", cols_alk2)
full_alk3_l <- lfc_vs_untreated(tmm_alk3,  common_genes_full, "ALK3",  cols_alk3_l)
full_alk3_a <- lfc_vs_untreated(tmm_alk3,  common_genes_full, "ALK3",  cols_alk3_a)


# Condition (Parental/Resistant) and Drug per column, derived from the sample
# name (actual pattern of a resistant line: underscore + L/A + number at the
# end -- "_L1", "_A10"..., not just any final digit, since the cell line
# names ALK1/ALK2/ALK3 also end in a digit).
derivar_condicion <- function(col) ifelse(grepl("_[LA][0-9]+$", col), "Resistant", "Parental")
derivar_drug <- function(col) {
  case_when(
    !grepl("_", col) ~ "No",
    grepl("_A", col)  ~ "Alectinib",
    TRUE              ~ "Lorlatinib"
  )
}

# =============================================================================
# PANEL A -- joint heatmap: the TMM matrices (not LFC) of the
# 3 models are merged into a single matrix (23 samples), a per-gene Z-score
# is computed on that joint matrix, and the top50 genes with the highest
# VARIANCE of that joint TMM matrix are selected (ranking is done BEFORE
# standardising -- once in Z-score all genes have variance 1 by
# definition, so the order has to come from the raw TMM values).
# The starting point is joint TMM, not LFC versus Untreated, and the selection
# is made on the joint matrix, not per cell line.
#
# The columns carry a Model annotation (ALK1/ALK2/ALK3) with its own
# legend, stacked next to Condition/Drug (ComplexHeatmap stacks all the
# legends of heatmap_legend_side/annotation_legend_side in one column).
# TOP_N_VARIANZA = 50 and the colour scale is capped at |Z| <= 3
# (Z_MAX_CONJUNTO).
# =============================================================================
cat("2. Generating Panel A (joint TMM heatmap + Z-score, top50 most variable)...\n")

TOP_N_VARIANZA <- 50

joint_tmm <- cbind(
  tmm_alk1[common_genes_full, cols_alk1,  drop = FALSE],
  tmm_alk2[common_genes_full, cols_alk2,  drop = FALSE],
  tmm_alk3[common_genes_full,  cols_alk3_l, drop = FALSE],
  tmm_alk3[common_genes_full,  cols_alk3_a, drop = FALSE]
)

var_genes <- apply(joint_tmm, 1, var, na.rm = TRUE)
genes_variables <- names(sort(var_genes, decreasing = TRUE))[1:TOP_N_VARIANZA]

mat_variable   <- joint_tmm[genes_variables, , drop = FALSE]
mat_variable_z <- t(scale(t(mat_variable)))  # per-gene (row) Z-score, on TMM

col_condicion_full <- setNames(derivar_condicion(colnames(mat_variable_z)), colnames(mat_variable_z))
col_drug_full      <- setNames(derivar_drug(colnames(mat_variable_z)), colnames(mat_variable_z))
col_modelo_full    <- setNames(
  c(rep("ALK1", ncol(tmm_alk1[common_genes_full, cols_alk1, drop = FALSE])),
    rep("ALK2", ncol(tmm_alk2[common_genes_full, cols_alk2, drop = FALSE])),
    rep("ALK3",  ncol(tmm_alk3[common_genes_full, cols_alk3_l, drop = FALSE])),
    rep("ALK3",  ncol(tmm_alk3[common_genes_full, cols_alk3_a, drop = FALSE]))),
  colnames(joint_tmm)
)

COLOR_MODELO_TMM <- c(ALK1 = "#89CFF0", ALK2 = "#DDA0DD", ALK3 = "#D9A441")

ann_col_conjunto <- HeatmapAnnotation(
  Model     = col_modelo_full[colnames(mat_variable_z)],
  Condition = col_condicion_full[colnames(mat_variable_z)],
  Drug      = col_drug_full[colnames(mat_variable_z)],
  col = list(Model = COLOR_MODELO_TMM, Condition = COLOR_CONDICION, Drug = COLOR_DRUG_FIRMA),
  show_annotation_name = FALSE,
  simple_anno_size = unit(3, "mm")
)

Z_MAX_CONJUNTO   <- min(3, max(abs(mat_variable_z), na.rm = TRUE))
COL_FUN_CONJUNTO <- circlize::colorRamp2(c(-Z_MAX_CONJUNTO, 0, Z_MAX_CONJUNTO), GRADIENTE_RNA)

ht_conjunto <- Heatmap(
  mat_variable_z,
  name = "Z-score",
  col = COL_FUN_CONJUNTO,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  top_annotation = ann_col_conjunto,
  row_names_gp = gpar(fontsize = 6.5),
  column_names_gp = gpar(fontsize = 8),
  rect_gp = gpar(col = "grey30", lwd = 1),
  heatmap_legend_param = list(title = "Z-score\n(TMM por gen)"),
  column_title = paste0("Top ", TOP_N_VARIANZA, " genes mas variables — TMM conjunto Z-score (ALK1 + ALK2 + ALK3)")
)

out_png_fig1a <- file.path(OUT_DIR, "figures", "ALK_09_fig1_panelA_heatmap_variables.png")
png(out_png_fig1a, width = 11, height = max(8, TOP_N_VARIANZA * 0.16) + 1, units = "in", res = 300)
draw(ht_conjunto, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()
cat("  ", out_png_fig1a, "\n\n")

# =============================================================================
# PANEL B-D -- PCA per cell line (ALK1 / ALK2 / ALK3), each with its own
# top TOP_N_PCA most variable genes (computed ONLY with its own samples,
# which prevents the PCA from being dominated by identity differences between
# lines instead of by the response to treatment/resistance within each
# one). Sample labels in black bold, stacked with patchwork.
#
# PCA instead of UMAP: with 5-11 samples per line UMAP is not reliable (it
# needs many more points for the notion of "neighbour" to make statistical
# sense), PCA is the standard choice for such a small n.
#
# geom_text_repel (ggrepel) instead of geom_text spreads the labels so that
# they do not overlap and pushes them away from the panel edge. It is
# combined with expansion(mult=0.18) on both axes and clip="off" so that no
# label gets cut off even if the point is near the limit.
#
# Fixed shapes per drug in the 3 panels (scale_shape_manual): circle =
# untreated parental, triangle = Lorlatinib, square = Alectinib. With the
# automatic assignment (alphabetical order of the Drug levels present in
# each model) the same drug came out with a different shape depending on the
# panel. Black thick panel.border, because theme_minimal() does not draw a
# panel border by default.
# =============================================================================
cat("3. Generating Panel B-D (PCA per cell line)...\n")

# The PCA is computed as in the RET analysis and on the SAME
# matrix as the individual heatmap of ALK_11: top 500 most variable genes of
# each line on TMM, per-gene Z-score and prcomp without rescaling (each gene
# already has mean 0 and sd 1).
TOP_N_PCA <- 500
FORMA_DRUG_PCA <- c(No = 16, Lorlatinib = 17, Alectinib = 15)

genes_validos <- function(mat) {
  g <- rownames(mat)
  g[!grepl("^ENSG", g) & !grepl("\\.", g)]
}

construir_pca_linea <- function(nombre_modelo, mat_tmm, cols, top_n) {
  genes_ok   <- intersect(genes_validos(mat_tmm), rownames(mat_tmm))
  mat_modelo <- mat_tmm[genes_ok, cols, drop = FALSE]

  var_g     <- apply(mat_modelo, 1, var, na.rm = TRUE)
  genes_top <- names(sort(var_g, decreasing = TRUE))[1:min(top_n, length(var_g))]
  mat_top   <- t(scale(t(mat_modelo[genes_top, , drop = FALSE])))  # per-gene Z-score

  pca <- prcomp(t(mat_top), center = TRUE, scale. = FALSE)
  var_exp <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)

  df_pca <- data.frame(
    Sample    = colnames(mat_top),
    PC1       = pca$x[, 1],
    PC2       = pca$x[, 2],
    Condition = derivar_condicion(colnames(mat_top)),
    Drug      = derivar_drug(colnames(mat_top))
  )

  cat("  ", nombre_modelo, ": PC1", var_exp[1], "% | PC2", var_exp[2], "%\n")

  # Same approach as the RET PCA: return the version with names (con) and
  # without names (sin), the latter to label by hand during layout.
  p_base <- ggplot(df_pca, aes(x = PC1, y = PC2, color = Condition, shape = Drug)) +
    geom_point(size = 4) +
    scale_color_manual(values = COLOR_CONDICION) +
    scale_shape_manual(values = FORMA_DRUG_PCA, drop = TRUE) +
    scale_x_continuous(expand = expansion(mult = 0.18)) +
    scale_y_continuous(expand = expansion(mult = 0.18)) +
    coord_cartesian(clip = "off") +
    labs(
      title = paste0(nombre_modelo, " — PCA (top ", length(genes_top), " most variable genes, TMM Z-score)"),
      x = paste0("PC1 (", var_exp[1], "%)"),
      y = paste0("PC2 (", var_exp[2], "%)")
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 12),
      plot.margin = margin(12, 16, 12, 16),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      panel.background = element_rect(fill = "white"),
      plot.background  = element_rect(fill = "white")
    )

  list(
    con = p_base +
      geom_text_repel(aes(label = Sample), color = "black", fontface = "bold", size = 3.3,
                       max.overlaps = Inf, box.padding = 0.5, point.padding = 0.3,
                       min.segment.length = 0, segment.size = 0.3, seed = 1, show.legend = FALSE),
    sin = p_base
  )
}

pcas_rna <- list(
  ALK1 = construir_pca_linea("ALK1", tmm_alk1, cols_alk1, TOP_N_PCA),
  ALK2 = construir_pca_linea("ALK2", tmm_alk2, cols_alk2, TOP_N_PCA),
  ALK3  = construir_pca_linea("ALK3",  tmm_alk3,  c(cols_alk3_l, cols_alk3_a), TOP_N_PCA)
)

titulo_pca <- "PCA by cell line — top 500 most variable genes (TMM, Z-score)"

for (version in c("con", "sin")) {
  p_pca_apilado <- wrap_plots(lapply(pcas_rna, function(x) x[[version]]), ncol = 1) +
    plot_annotation(title = titulo_pca,
                    theme = theme(plot.title = element_text(face = "bold", size = 14)))
  sufijo <- if (version == "sin") "_sin_nombres" else ""
  out_png_fig1bcd <- file.path(OUT_DIR, "figures", paste0("ALK_09_fig1_panelBCD_pca_lineas", sufijo, ".png"))
  ggsave(out_png_fig1bcd, p_pca_apilado, width = 8, height = 15, dpi = 300)
  cat("  ", out_png_fig1bcd, "\n")
}

# Individual PNG per line, with and without names
for (m in names(pcas_rna)) {
  out_m     <- file.path(OUT_DIR, "figures", paste0("ALK_09_fig1_pca_", m, ".png"))
  out_m_sin <- file.path(OUT_DIR, "figures", paste0("ALK_09_fig1_pca_", m, "_sin_nombres.png"))
  ggsave(out_m,     pcas_rna[[m]]$con, width = 7, height = 6.5, dpi = 300)
  ggsave(out_m_sin, pcas_rna[[m]]$sin, width = 7, height = 6.5, dpi = 300)
  cat("  ", out_m, "\n  ", out_m_sin, "\n")
}
cat("\n")

cat("=== ALK_09_combined_figure.R — Panel A + B-D completed ===\n")
