# =============================================================================
# RET_08b_pca.R
#
# Purpose: RNA-seq PCA of RET1 and RET2 with and without sample names, and a
# full figure (PCA + heatmaps) using the PCA without names, so that samples
# can be labelled manually (same approach as MET_01_summary_figures.R).
#
# The computation is identical to RET_08_individual_heatmaps.R: top
# 500 most variable genes of each model on TMM, per-gene Z-score, prcomp on
# the samples and a ComplexHeatmap heatmap of the same matrix. Same colors,
# shapes, legends and layout.
#
# Inputs:
#   results/RET/{RET1,RET2}/tables/01_mat_tmm.csv
#   results/RET/{RET1,RET2}/tables/01_metadata_unit.csv
#
# Outputs (results/RET/combined/figures/):
#   RET_08b_fig25_pca_{RET1,RET2}.png              -- with names   
#   RET_08b_fig25_pca_{RET1,RET2}_sin_nombres.png  -- without names
#   RET_08b_fig25_panel_combinado_sin_nombres.png   -- PCA without names + heatmaps
#
# Dependencies (run first): scripts/simple/rna_01.1_qc.R with
# DATASET = "RET1" and DATASET = "RET2".
#
# Usage: paths are resolved with here::here() from the repository root
# (.here file), so it works from any working directory:
#   Rscript scripts/specific/RET/RET_08b_pca.R
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

OUT_DIR <- here::here("results", "RET", "combined", "figures")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

MODELOS        <- c("RET1", "RET2")
TOP_N_VARIANZA <- 500

COLOR_CONDICION     <- c(Parental = "#e6a817", Resistente = "#013c58")
FORMA_INHIBIDOR_PCA <- c(No = 16, BLU667 = 17, LOXO292 = 15)
COLOR_DRUG          <- c(No = "#95A5A6", BLU667 = "#6C5CE7", LOXO292 = "#00B894")

GRADIENTE_RNA <- c("#000000", "#808080", "#f2a900")
Z_MAX         <- 2
COL_FUN_RNA   <- circlize::colorRamp2(c(-Z_MAX, 0, Z_MAX), GRADIENTE_RNA)

ANCHO_HEATMAP <- 8
ALTO_HEATMAP  <- 10
ALTO_PCA      <- 6.5

# ── Data (same functions as RET_08) ──────────────────────────────────────
leer_tmm <- function(modelo) {
  df <- read.csv(here::here("results", "RET", modelo, "tables", "01_mat_tmm.csv"),
                 check.names = FALSE)
  rownames(df) <- df$Gene
  as.matrix(df[, -1, drop = FALSE])
}

leer_metadata <- function(modelo) {
  meta <- read.csv(here::here("results", "RET", modelo, "tables", "01_metadata_unit.csv"),
                   stringsAsFactors = FALSE)
  meta$Condition3 <- ifelse(meta$Condition == "Resistance", "Resistente", "Parental")
  meta
}

genes_validos <- function(mat) {
  g <- rownames(mat)
  g[!grepl("^ENSG", g) & !grepl("\\.", g)]
}

seleccionar_top_z <- function(mat, meta, top_n) {
  cols <- meta$Sample[meta$Sample %in% colnames(mat)]
  genes_ok   <- intersect(genes_validos(mat), rownames(mat))
  mat_modelo <- mat[genes_ok, cols, drop = FALSE]

  var_g     <- apply(mat_modelo, 1, var, na.rm = TRUE)
  genes_top <- names(sort(var_g, decreasing = TRUE))[1:min(top_n, length(var_g))]
  mat_z     <- t(scale(t(mat_modelo[genes_top, , drop = FALSE])))

  list(mat_z = mat_z, meta_ord = meta[match(colnames(mat_z), meta$Sample), ],
       n_genes = length(genes_top))
}

# ── PCA: returns the version with and without names ───────────────────────────────
construir_pca <- function(sel, modelo) {
  pca     <- prcomp(t(sel$mat_z), center = TRUE, scale. = FALSE)
  var_exp <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)

  df_pca <- data.frame(
    Sample    = colnames(sel$mat_z),
    PC1       = pca$x[, 1],
    PC2       = pca$x[, 2],
    Condition = sel$meta_ord$Condition3,
    Drug      = sel$meta_ord$Drug
  )
  cat(modelo, ": PC1", var_exp[1], "% | PC2", var_exp[2], "%\n")

  p_base <- ggplot(df_pca, aes(x = PC1, y = PC2, color = Condition, shape = Drug)) +
    geom_point(size = 4) +
    scale_color_manual(values = COLOR_CONDICION, name = "Condición") +
    scale_shape_manual(values = FORMA_INHIBIDOR_PCA, drop = TRUE, name = "Inhibidor") +
    scale_x_continuous(expand = expansion(mult = 0.18)) +
    scale_y_continuous(expand = expansion(mult = 0.18)) +
    coord_cartesian(clip = "off") +
    labs(
      title = paste0(modelo, " — PCA (top ", sel$n_genes, " most variable genes, TMM Z-score)"),
      x = paste0("PC1 (", var_exp[1], "%)"),
      y = paste0("PC2 (", var_exp[2], "%)")
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title       = element_text(face = "bold", size = 12),
      plot.margin      = margin(12, 16, 12, 16),
      panel.border     = element_rect(color = "black", fill = NA, linewidth = 1),
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

# ── Heatmap (same as construir_heatmap_modelo() in RET_08, without its own
#    legend because in the combined figure it is shared at the bottom) ───────
construir_heatmap <- function(sel, modelo) {
  ann_col <- HeatmapAnnotation(
    `Condición` = setNames(sel$meta_ord$Condition3, colnames(sel$mat_z)),
    Inhibidor   = setNames(sel$meta_ord$Drug, colnames(sel$mat_z)),
    col = list(`Condición` = COLOR_CONDICION, Inhibidor = COLOR_DRUG),
    show_annotation_name = FALSE,
    show_legend = FALSE,
    simple_anno_size = unit(3, "mm")
  )

  Heatmap(
    sel$mat_z,
    name = "Z-score",
    col = COL_FUN_RNA,
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    show_row_names = FALSE,
    show_heatmap_legend = FALSE,
    top_annotation = ann_col,
    column_names_gp = gpar(fontsize = 9),
    rect_gp = gpar(col = NA),
    column_title = paste0(modelo, " -- top ", sel$n_genes, " most variable genes (TMM, Z-score)")
  )
}

# ── Individual PCAs ─────────────────────────────────────────────────────────
plots_pca_sin <- list()
heatmaps      <- list()

for (modelo in MODELOS) {
  sel <- seleccionar_top_z(leer_tmm(modelo), leer_metadata(modelo), TOP_N_VARIANZA)
  p   <- construir_pca(sel, modelo)
  plots_pca_sin[[modelo]] <- p$sin
  heatmaps[[modelo]]      <- construir_heatmap(sel, modelo)

  out_con <- file.path(OUT_DIR, paste0("RET_08b_fig25_pca_", modelo, ".png"))
  out_sin <- file.path(OUT_DIR, paste0("RET_08b_fig25_pca_", modelo, "_sin_nombres.png"))
  ggsave(out_con, p$con, width = 7, height = 6.5, dpi = 300)
  ggsave(out_sin, p$sin, width = 7, height = 6.5, dpi = 300)
  cat("  ", out_con, "\n  ", out_sin, "\n")
}

# ── Combined figure: PCA without names on top, heatmaps below and the shared
#    heatmap legend at the bottom (same layout as RET_08) ────────────────────
leyenda_compartida <- packLegend(
  Legend(col_fun = COL_FUN_RNA, title = "Z-score (TMM per gene)", direction = "horizontal",
         legend_width = unit(4, "cm")),
  Legend(labels = names(COLOR_CONDICION), legend_gp = gpar(fill = COLOR_CONDICION),
         title = "Condición", direction = "horizontal"),
  Legend(labels = names(COLOR_DRUG), legend_gp = gpar(fill = COLOR_DRUG),
         title = "Inhibidor", direction = "horizontal"),
  direction = "horizontal"
)

out_combinado <- file.path(OUT_DIR, "RET_08b_fig25_panel_combinado_sin_nombres.png")
png(out_combinado, width = 2 * ANCHO_HEATMAP, height = ALTO_PCA + ALTO_HEATMAP + 1.5,
    units = "in", res = 300)
grid.newpage()
pushViewport(viewport(layout = grid.layout(3, 2,
  heights = unit(c(ALTO_PCA, ALTO_HEATMAP, 1.5), "inches"),
  widths  = unit(c(1, 1), "null"))))

for (i in seq_along(MODELOS)) {
  pushViewport(viewport(layout.pos.row = 1, layout.pos.col = i))
  print(plots_pca_sin[[MODELOS[i]]], newpage = FALSE)
  upViewport()

  pushViewport(viewport(layout.pos.row = 2, layout.pos.col = i))
  draw(heatmaps[[MODELOS[i]]], newpage = FALSE)
  upViewport()
}

pushViewport(viewport(layout.pos.row = 3, layout.pos.col = 1:2))
draw(leyenda_compartida, just = "center")
upViewport()

dev.off()
cat("  ", out_combinado, "\n")
