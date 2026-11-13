# =============================================================================
# RET_08_figura1_heatmaps_individuales.R
#
# Purpose: PCA and expression heatmap per RET model (RET1, RET2). RET version
# of ALK_11_figura1_heatmaps_individuales.R with 2 columns instead of 3.
#   - Heatmap: top 500 most variable genes WITHIN each model (TMM),
#     per-gene Z-score, fixed -2/2 scale in both models.
#   - PCA: computed on the SAME top 500 Z-score matrix as the heatmap
#     (seleccionar_top_z()), not on its own top-N, so that both panels are
#     strictly comparable. Visual design (color = Condition, shape = Inhibitor,
#     geom_text_repel, panel border) follows the ALK RNA-seq PCA
#     (ALK_09_figura4_combinada.R).
#   - Combined panel: PCA (top row) + heatmap (bottom row), 2 rows x 2 columns.
# No reanalysis: reads the TMM matrices and metadata produced by
# rna_01.1_qc.R for RET.
#
# Per-sample metadata: 01_metadata_unit.csv is used (columns
# Model/Sample/Cell_Line/Condition/Treatment/Drug already resolved) instead of
# inferring the condition from the column name by regex. Mapping (RET1; RET2
# is identical with BR1-5 and LR1/LR3/LR4, there is no LR2):
#   RET1          -> Condition=Parental,   Treatment=No,  Drug=No
#   RET1_B250     -> Condition=Parental,   Treatment=Yes, Drug=BLU667   (acute)
#   RET1_L250     -> Condition=Parental,   Treatment=Yes, Drug=LOXO292  (acute)
#   RET1_BR1/2/3  -> Condition=Resistance, Treatment=Yes, Drug=BLU667
#   RET1_LR1/2/3  -> Condition=Resistance, Treatment=Yes, Drug=LOXO292
# Condition3 collapses Condition into 2 levels (Parental/Resistente): the
# acutely treated parental (B250/L250) counts as Parental and the acute
# treatment is still carried by the Inhibidor annotation (No/BLU667/LOXO292).
# Same criterion as ALK (Parental/Resistant, Drug separate).
#
# Labels and axes are in English; titles and categories of the color legends
# are in Spanish ("Condición", "Inhibidor"). Palettes are defined inline,
# without sourcing config/colors.R (same approach as the equivalent ALK
# scripts). COLOR_DRUG does NOT use drug_colors from config/colors.R (ALK
# scheme) but the BLU-667 purple / LOXO-292 green colors used in the other
# RET scripts (the RET WES oncoplot, rna_08.1_fusions_ret.R).
#
# Inputs:
#   results/RET/RET1/tables/01_mat_tmm.csv
#   results/RET/RET1/tables/01_metadata_unit.csv
#   results/RET/RET2/tables/01_mat_tmm.csv
#   results/RET/RET2/tables/01_metadata_unit.csv
#
# Outputs (results/RET/combined/figures/):
#   RET_08_fig1_pca_RET1.png / _RET2.png       (PCA)
#   RET_08_fig1_heatmap_RET1.png / _RET2.png   (heatmap)
#   RET_08_fig1_panel_combinado.png             (PCA + heatmap, 2 x 2)
#
# Dependencies (run first): scripts/simple/rna_01.1_qc.R with
# DATASET = "RET1" and DATASET = "RET2".
#
# Usage: paths are resolved with here::here() from the repository root
# (.here file), so it works from any working directory:
#   Rscript scripts/specific/RET/RET_08_figura1_heatmaps_individuales.R
# =============================================================================

suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
  library(ggplot2)  # PCA
  library(ggrepel)  # geom_text_repel for the PCA
})

OUT_DIR <- here::here("results", "RET", "combined")
dir.create(file.path(OUT_DIR, "figures"), recursive = TRUE, showWarnings = FALSE)

# Shared palette, identical to ALK_09/ALK_11 (black-grey-gold for RNA-seq,
# technique identity shared across all figures).
GRADIENTE_RNA <- c("#000000", "#808080", "#f2a900")

# Condition3 (Parental/Resistente), 2 levels. The treated/untreated distinction
# of the parental is still available via Inhibidor (shape in the PCA, color
# strip in the heatmap: No/BLU667/LOXO292), without duplicating it as its own
# Condition level. Same criterion as ALK (Parental/Resistant, Drug separate).
COLOR_CONDICION <- c(Parental = "#e6a817", Resistente = "#013c58")

# COLOR_DRUG aligned with the RET WES oncoplot and rna_08.1_fusions_
# ret.R (fusions): BLU-667 purple / LOXO-292 green in ALL RET figures, not the
# salmon/cyan values of drug_colors in config/colors.R (those belong to the
# ALK/Lorlatinib-Alectinib scheme and are not reused here).
COLOR_DRUG  <- c(No = "#95A5A6", BLU667 = "#6C5CE7", LOXO292 = "#00B894")

TOP_N_VARIANZA <- 500

# FIXED Z-score scale -2/2, identical in both models (same as ALK_11, for
# visual comparability between panels).
Z_MAX       <- 2
COL_FUN_RNA <- circlize::colorRamp2(c(-Z_MAX, 0, Z_MAX), GRADIENTE_RNA)

# FIXED, identical dimensions for the 2 individual PNGs (same as ALK_11).
ANCHO_HEATMAP <- 8
ALTO_HEATMAP  <- 10

# =============================================================================
# 1. Read TMM + metadata per model
# =============================================================================
cat("1. Reading TMM and metadata per model...\n")

leer_tmm <- function(modelo) {
  ruta <- here::here("results", "RET", modelo, "tables", "01_mat_tmm.csv")
  df <- read.csv(ruta, check.names = FALSE)
  rownames(df) <- df$Gene
  as.matrix(df[, -1, drop = FALSE])
}

leer_metadata <- function(modelo) {
  ruta <- here::here("results", "RET", modelo, "tables", "01_metadata_unit.csv")
  meta <- read.csv(ruta, stringsAsFactors = FALSE)
  # Parental/Resistente only: the acutely treated parental counts as
  # "Parental" like the untreated one; the acute treatment distinction is
  # still carried by the Inhibidor annotation.
  meta$Condition3 <- ifelse(meta$Condition == "Resistance", "Resistente", "Parental")
  meta
}

tmm_ret1 <- leer_tmm("RET1")
tmm_ret2  <- leer_tmm("RET2")
meta_ret1 <- leer_metadata("RET1")
meta_ret2  <- leer_metadata("RET2")

genes_validos <- function(mat) {
  g <- rownames(mat)
  g[!grepl("^ENSG", g) & !grepl("\\.", g)]
}

# =============================================================================
# 2. One heatmap per model: top 500 most variable genes WITHIN that model,
#    per-gene Z-score, no row names (global view); same criterion as ALK_11.
# =============================================================================
cat("2. Generating individual heatmaps (top", TOP_N_VARIANZA, "most variable per model)...\n")

# Shared function so that the PCA (section 2b) reuses EXACTLY the same top 500
# Z-score matrix as the heatmap, instead of recomputing its own top-N (which
# could select a slightly different gene set and make PCA and heatmap not
# comparable).
seleccionar_top_z <- function(mat, meta, top_n) {
  cols <- meta$Sample
  cols <- cols[cols %in% colnames(mat)]
  genes_ok   <- intersect(genes_validos(mat), rownames(mat))
  mat_modelo <- mat[genes_ok, cols, drop = FALSE]

  var_g     <- apply(mat_modelo, 1, var, na.rm = TRUE)
  genes_top <- names(sort(var_g, decreasing = TRUE))[1:min(top_n, length(var_g))]
  mat_top   <- mat_modelo[genes_top, , drop = FALSE]
  mat_z     <- t(scale(t(mat_top)))  # per-gene (row) Z-score on TMM

  meta_ord <- meta[match(colnames(mat_z), meta$Sample), ]
  list(mat_z = mat_z, meta_ord = meta_ord, n_genes = length(genes_top))
}

construir_heatmap_modelo <- function(modelo, mat, meta, top_n, mostrar_leyenda = TRUE) {
  sel      <- seleccionar_top_z(mat, meta, top_n)
  mat_z    <- sel$mat_z
  meta_ord <- sel$meta_ord

  # Color legend in Spanish (title + categories): "Condición" and "Inhibidor"
  # (same term as the ALK legends, not "Fármaco"). Argument names are in
  # backticks because of the accent; the name passed here is what
  # ComplexHeatmap uses as the default legend title (show_annotation_name=FALSE
  # only hides the side label of the color strip, not the legend title).
  ann_col <- HeatmapAnnotation(
    `Condición` = setNames(meta_ord$Condition3, colnames(mat_z)),
    Inhibidor   = setNames(meta_ord$Drug, colnames(mat_z)),
    col = list(`Condición` = COLOR_CONDICION, Inhibidor = COLOR_DRUG),
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
    column_title = paste0(modelo, " -- top ", sel$n_genes, " most variable genes (TMM, Z-score)")
  )
}

ht_ret1 <- construir_heatmap_modelo("RET1", tmm_ret1, meta_ret1, TOP_N_VARIANZA)
ht_ret2  <- construir_heatmap_modelo("RET2",  tmm_ret2,  meta_ret2,  TOP_N_VARIANZA)

heatmaps <- list(RET1 = ht_ret1, RET2 = ht_ret2)

for (modelo in names(heatmaps)) {
  out <- file.path(OUT_DIR, "figures", paste0("RET_08_fig1_heatmap_", modelo, ".png"))
  png(out, width = ANCHO_HEATMAP, height = ALTO_HEATMAP, units = "in", res = 300)
  draw(heatmaps[[modelo]], heatmap_legend_side = "right", annotation_legend_side = "right")
  dev.off()
  cat("  ", out, "\n")
}

# =============================================================================
# 2b. PCA per model, COMPUTED ON THE SAME top 500 Z-score MATRIX used by the
#     heatmap (seleccionar_top_z(), section 2), instead of its own top-N.
#     Visual design (color=Condition, shape=Inhibitor, geom_text_repel, bold
#     panel border) follows the ALK RNA-seq PCA (ALK_09_figura4_combinada.R;
#     ALK_11 has no PCA).
#
#     Color = Condition (Parental/Resistente, 2 levels, see Condition3 in
#     section 1). Shape = Inhibitor (No/BLU667/LOXO292), same scheme as ALK
#     (Drug as shape, not as color).
# =============================================================================
cat("2b. Generating PCA per model (same top 500 Z-score matrix)...\n")

FORMA_INHIBIDOR_PCA <- c(No = 16, BLU667 = 17, LOXO292 = 15)

construir_pca_modelo_rna <- function(modelo, mat, meta, top_n) {
  sel      <- seleccionar_top_z(mat, meta, top_n)
  mat_z    <- sel$mat_z
  meta_ord <- sel$meta_ord

  # prcomp on the SAMPLES (t(mat_z)); each row of mat_z is already a Z-scored
  # gene (mean 0, sd 1 by construction), so scale.=FALSE here (already
  # standardized gene by gene before the PCA; rescaling per sample column
  # makes no sense for this matrix).
  pca <- prcomp(t(mat_z), center = TRUE, scale. = FALSE)
  var_exp <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)

  df_pca <- data.frame(
    Sample    = colnames(mat_z),
    PC1       = pca$x[, 1],
    PC2       = pca$x[, 2],
    Condition = meta_ord$Condition3,
    Drug      = meta_ord$Drug
  )

  ggplot(df_pca, aes(x = PC1, y = PC2, color = Condition, shape = Drug)) +
    geom_point(size = 4) +
    geom_text_repel(aes(label = Sample), color = "black", fontface = "bold", size = 3.3,
                     max.overlaps = Inf, box.padding = 0.5, point.padding = 0.3,
                     min.segment.length = 0, segment.size = 0.3, seed = 1, show.legend = FALSE) +
    # Color legend in Spanish: "Condición"/"Inhibidor".
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
      plot.title = element_text(face = "bold", size = 12),
      plot.margin = margin(12, 16, 12, 16),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      panel.background = element_rect(fill = "white"),
      plot.background  = element_rect(fill = "white")
    )
}

p_pca_ret1 <- construir_pca_modelo_rna("RET1", tmm_ret1, meta_ret1, TOP_N_VARIANZA)
p_pca_ret2  <- construir_pca_modelo_rna("RET2",  tmm_ret2,  meta_ret2,  TOP_N_VARIANZA)

plots_pca <- list(RET1 = p_pca_ret1, RET2 = p_pca_ret2)

for (modelo in names(plots_pca)) {
  out <- file.path(OUT_DIR, "figures", paste0("RET_08_fig1_pca_", modelo, ".png"))
  ggsave(out, plots_pca[[modelo]], width = 7, height = 6.5, dpi = 300)
  cat("  ", out, "\n")
}

# =============================================================================
# 3. Combined panel: PCA (top row) + heatmap (bottom row), 2 rows x 2 columns
#    (RET1, RET2). Mixes ggplot (PCA) + ComplexHeatmap (heatmap) in the SAME
#    PNG via grid viewports. ComplexHeatmap is not a ggplot/patchwork grob, so
#    patchwork's `+`/`/` cannot be used here; each ggplot is drawn with
#    print(p, newpage=FALSE) inside its own viewport, just like draw() for the
#    Heatmap objects.
# =============================================================================
cat("3. Generating combined panel (PCA + heatmap)...\n")

ht_ret1_comb <- construir_heatmap_modelo("RET1", tmm_ret1, meta_ret1, TOP_N_VARIANZA, mostrar_leyenda = FALSE)
ht_ret2_comb  <- construir_heatmap_modelo("RET2",  tmm_ret2,  meta_ret2,  TOP_N_VARIANZA, mostrar_leyenda = FALSE)

leyenda_compartida <- packLegend(
  Legend(col_fun = COL_FUN_RNA, title = "Z-score (TMM per gene)", direction = "horizontal", legend_width = unit(4, "cm")),
  Legend(labels = names(COLOR_CONDICION), legend_gp = gpar(fill = COLOR_CONDICION), title = "Condición", direction = "horizontal"),
  Legend(labels = names(COLOR_DRUG), legend_gp = gpar(fill = COLOR_DRUG), title = "Inhibidor", direction = "horizontal"),
  direction = "horizontal"
)

out_combinado <- file.path(OUT_DIR, "figures", "RET_08_fig1_panel_combinado.png")

# PCA row (ALTO_PCA, with its own color/shape legend included in each ggplot)
# + heatmap row (ALTO_HEATMAP) + a final strip for the heatmap's shared
# legend. The PCA legend is not shared via packLegend (mixing ggplot legends
# with ComplexHeatmap ones in one Legend object is not trivial): each PCA
# keeps its own Condición/Inhibidor legend at the right of the panel, and the
# heatmap keeps its own (Z-score/Condición/Inhibidor) in the shared strip
# below.
ALTO_PCA <- 6.5
png(out_combinado, width = 2 * ANCHO_HEATMAP, height = ALTO_PCA + ALTO_HEATMAP + 1.5, units = "in", res = 300)
grid.newpage()
pushViewport(viewport(layout = grid.layout(3, 2,
  heights = unit(c(ALTO_PCA, ALTO_HEATMAP, 1.5), c("inches", "inches", "inches")),
  widths  = unit(c(1, 1), "null"))))

pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
print(plots_pca[["RET1"]], newpage = FALSE)
upViewport()

pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
print(plots_pca[["RET2"]], newpage = FALSE)
upViewport()

pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1))
draw(ht_ret1_comb, newpage = FALSE)
upViewport()

pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 2))
draw(ht_ret2_comb, newpage = FALSE)
upViewport()

pushViewport(viewport(layout.pos.row = 3, layout.pos.col = 1:2))
draw(leyenda_compartida, just = "center")
upViewport()

dev.off()
cat("  ", out_combinado, "\n")

cat("\n=== RET_08_figura1_heatmaps_individuales.R finished (4 individual PNGs: 2 PCA + 2 heatmap; + 1 combined PCA+heatmap panel) ===\n")
