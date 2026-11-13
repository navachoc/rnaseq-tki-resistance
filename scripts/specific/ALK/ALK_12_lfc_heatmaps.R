# =============================================================================
# ALK_12_lfc_heatmaps.R
# Purpose: LFC heatmaps per model (ALK1/ALK2/ALK3), 3 individual PNGs
#   + 1 combined panel (the 3 side by side, quick reference), same format
#   as ALK_11_individual_heatmaps.R. A combined heatmap with the
#   3 models as columns would be fragmented (the top-LFC genes of one model
#   are flat in the other 2).
#
# Criteria:
#  1. Comparator: ALWAYS Resistant (or treated parental) vs Untreated
#     Parental, not vs Treated (which is the primary comparator of
#     `02_lfc_master_table.csv`).
#  2. 1vs1 LFC: the LFC is computed per INDIVIDUAL LINE against the untreated
#     parental of its model (not an average of all the resistant lines
#     of the model, which is what `02_lfc_master_table.csv` does). The heatmap
#     shows one column per line, with its real line-by-line LFC.
#  3. Columns: acutely treated parental(s) + resistant lines
#     (structure in `COLUMNAS`, with explicit label/condition/drug per
#     column). The UNTREATED parental is not a column (it is only used as
#     the reference to subtract, its LFC against itself would be 0 in all
#     rows, explained in the figure caption). The treated parental does NOT
#     take part in the Top LFC ranking (seleccionar_top_lfc() only looks at
#     the columns Condition=="Resistant"), it only provides visual context of
#     the acute pre-resistance response.
#  4. Gene selection per panel (TOP_N_LFC = 20 per direction, 40 top-LFC
#     + 10 WES/receptors):
#     - 40 top-LFC: for each gene, the most extreme value in magnitude (|LFC|)
#       among ALL the resistant lines of THAT model decides its ranking.
#       Top 20 with positive sign (Up) + top 20 with negative sign (Down).
#       Genes with a bare ENSG without symbol (`^ENSG[0-9]+$`) are EXCLUDED
#       from the ranking -- readthrough genes with a compound symbol (e.g.
#       "GASK1B.AS1", "P3R3URF.PIK3R3") DO count, they have a real transcript.
#     - 10 WES/receptors (fixed, the same in the 3 panels, manually curated
#       list): ALK, APOBEC3A, MET, SMO, KRAS, EGFR, STK11, KEAP1,
#       BRCA2, CDK12. EML4 is not included (it is only the fusion partner of ALK,
#       redundant signal). APOBEC3A fits the APOBEC mutagenesis hypothesis
#       explored in WES (elevated APOBEC score in ALK1 L1-L3 although
#       without significance at whole-exome level; significant enrichment
#       in ALK3_L9), its expression LFC is the RNA complement
#       of that signal.
#     - If a model lacks genes of the WES/receptor panel in its TMM matrix
#       (APOBEC3A is not quantified in ALK1/ALK2),
#       `n_extra_up_down()` extends the top-LFC by that same amount
#       (split Up/Down) so that the 3 panels add up to 50 genes.
#     - Union of both blocks per panel (may be <50 if some gene of the
#       block of 10 was already in the top-LFC of that model).
#  5. Rows: a single block with cluster_rows=TRUE (no row_split between
#     Top-LFC and WES/Receptor); the Category strip is kept as an
#     informative annotation, without a physical gap. Columns not clustered,
#     in the order of COLUMNAS.
#  6. Model/Condition/Drug annotations above each column, same palette
#     as ALK_11 (Model with the same colours as the final WES oncoplot).
#  7. FIXED colour scale from -LFC_MAX to LFC_MAX (LFC_MAX = 8). Genes of the
#     Top LFC block with |LFC| > LFC_MAX saturate at the extreme colour, an
#     expected effect explained in the figure caption.
#
# Does NOT reanalyse anything -- reads the same TMM matrices (01_mat_tmm.csv, already in
# log2) as ALK_09_combined_figure.R and ALK_11_individual_heatmaps.R.
#
# Inputs:
#   results/ALK/ALK1/tables/01_mat_tmm.csv
#   results/ALK/ALK2/tables/01_mat_tmm.csv
#   results/ALK/ALK3/tables/01_mat_tmm.csv
#
# Outputs (results/ALK/combined/figures/):
#   ALK_12_fig2_lfc_heatmap_ALK1.png
#   ALK_12_fig2_lfc_heatmap_ALK2.png
#   ALK_12_fig2_lfc_heatmap_ALK3.png
#   ALK_12_fig2_lfc_heatmap_combinado.png   (the 3 side by side, quick reference)
#   ALK_12_fig2_lfc_heatmap_<model>_gene_order.csv (row order after clustering)
#   ALK_12_fig2_lfc_heatmap_pie.txt (figure caption text)
#
# Labels, axes and legends in English.
#
# Usage: source(here::here("scripts/specific/ALK/ALK_12_lfc_heatmaps.R"))
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

# Same RNA palette as ALK_09/ALK_11 (black-grey-gold) -- colour identity
# shared across all RNA-seq figures.
GRADIENTE_RNA <- c("#000000", "#808080", "#f2a900")

# Colours of the category strip (Top LFC vs WES/Receptor panel) -- a new
# pair, different from the Founder/Resistance of the WES oncoplot (they mean
# different things, do not reuse the same pair to avoid confusion).
COLOR_CATEGORY <- c("Top LFC (Res vs Unt, per line)" = "#2C7FB8",
                     "WES/Receptor panel"             = "#D95F02")

# Condition/Drug annotation above each column -- same palette as
# ALK_11_individual_heatmaps.R (shared colour identity).
# The untreated parental is not shown; the columns are the treated
# parental(s) (Condition = "Parental") and the resistant lines
# (Condition = "Resistant").
COLOR_CONDICION  <- c(Parental = "#e6a817", Resistant = "#013c58")
COLOR_DRUG_FIRMA <- c(No = "#95A5A6", Lorlatinib = "#00B894", Alectinib = "#F4978E")
# Model also as an annotation (constant per panel, for consistency with
# the rest of the figures) -- same colours as the Model strip of the final
# WES oncoplot, colour identity shared between WES and RNA.
COLOR_MODELO <- c(ALK1 = "#89CFF0", ALK2 = "#DDA0DD", ALK3 = "#D9A441")

TOP_N_LFC <- 20  # per direction (Up/Down) -- 40 top-LFC genes in total, except deficit (see section 2)

# Manually curated list. APOBEC3A instead of EML4 (fusion partner of
# ALK, redundant signal) -- fits the APOBEC mutagenesis hypothesis
# explored in WES for ALK1/ALK3 (elevated APOBEC score in ALK1 L1-L3
# although without significance at whole-exome level, and significant
# enrichment in ALK3_L9) -- its expression LFC here is the RNA complement
# of that WES signal.
GENES_WES_RECEPTORES <- c("ALK", "APOBEC3A", "MET", "SMO", "KRAS", "EGFR",
                           "STK11", "KEAP1", "BRCA2", "CDK12")

# FIXED, identical dimensions for the 3 individual PNGs (same criterion
# as ALK_11) -- ~50 genes with legible names need more height than the 500
# unnamed ones of ALK_11.
ANCHO_HEATMAP <- 8
ALTO_HEATMAP  <- 13

# =============================================================================
# 1. Read TMM per model + resistant line / parental columns
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

# Untreated parental = LFC reference (never a column of the heatmap,
# only used for subtraction). Columns that ARE shown: acutely treated
# parental(s) + resistant lines. ALK3 has 2
# treated parentals (one per drug, Lorlatinib/Alectinib) because it has
# 2 drugs; ALK1/ALK2 only have Lorlatinib.
UNTREATED <- list(ALK1 = "ALK1", ALK2 = "ALK2", ALK3 = "ALK3")

COLUMNAS <- list(
  ALK1 = data.frame(
    raw       = c("ALK1_L", "ALK1_L1", "ALK1_L2", "ALK1_L3"),
    label     = c("Treated", "L1", "L2", "L3"),
    condition = c("Parental", "Resistant", "Resistant", "Resistant"),
    drug      = c("Lorlatinib", "Lorlatinib", "Lorlatinib", "Lorlatinib"),
    stringsAsFactors = FALSE
  ),
  ALK2 = data.frame(
    raw       = c("ALK2_L", paste0("ALK2_L", 1:5)),
    label     = c("Treated", paste0("L", 1:5)),
    condition = c("Parental", rep("Resistant", 5)),
    drug      = rep("Lorlatinib", 6),
    stringsAsFactors = FALSE
  ),
  ALK3 = data.frame(
    raw       = c("ALK3_L", "ALK3_A", "ALK3_L2", "ALK3_L3", "ALK3_L7", "ALK3_L9", "ALK3_L10", "ALK3_A2", "ALK3_A3", "ALK3_A4"),
    label     = c("Treated_L", "Treated_A", "L2", "L3", "L7", "L9", "L10", "A2", "A3", "A4"),
    condition = c("Parental", "Parental", rep("Resistant", 8)),
    drug      = c("Lorlatinib", "Alectinib", "Lorlatinib", "Lorlatinib", "Lorlatinib", "Lorlatinib", "Lorlatinib", "Alectinib", "Alectinib", "Alectinib"),
    stringsAsFactors = FALSE
  )
)
MATS <- list(ALK1 = tmm_alk1, ALK2 = tmm_alk2, ALK3 = tmm_alk3)

# Genes with a bare ENSG without symbol -- excluded from the top-LFC RANKING
# (no annotated transcript). Different from the ALK_11 filter (which also
# excluded any gene with a dot in its name) -- here readthrough genes with a
# compound symbol DO count, only the bare Ensembl ID is discarded.
es_ensg_puro <- function(g) grepl("^ENSG[0-9]+$", g)

# =============================================================================
# 2. 1vs1 LFC (each column of COLUMNAS[[modelo]] - untreated parental of the
#    same model, including the treated one(s)) + selection of the
#    top-LFC genes (20 Up + 20 Down by default) by maximum magnitude AMONG THE
#    RESISTANT LINES (the treated parental does not take part in the ranking,
#    it is only shown as an additional context column).
# =============================================================================
cat("2. Computing 1vs1 LFC (column - Untreated parental) and selecting genes...\n")

calcular_lfc_1v1 <- function(modelo) {
  mat  <- MATS[[modelo]]
  u    <- mat[, UNTREATED[[modelo]]]
  cols <- COLUMNAS[[modelo]]
  lfc  <- sapply(cols$raw, function(cl) mat[, cl] - u)
  rownames(lfc) <- rownames(mat)
  colnames(lfc) <- cols$label
  lfc
}

# n_up/n_down make it possible to compensate the deficit of a model when it
# lacks genes of the WES/receptor panel in its TMM matrix (ALK1 and ALK2
# lack APOBEC3A, not quantified in those 2 matrices, so the
# gap is filled with 1 extra top-LFC gene so that the 3 panels have 50 genes).
seleccionar_top_lfc <- function(lfc_mat, modelo, n_up = TOP_N_LFC, n_down = TOP_N_LFC) {
  cols_resistentes <- COLUMNAS[[modelo]]$label[COLUMNAS[[modelo]]$condition == "Resistant"]
  lfc_res <- lfc_mat[, cols_resistentes, drop = FALSE]
  genes_ok <- rownames(lfc_res)[!es_ensg_puro(rownames(lfc_res))]
  lfc_ok   <- lfc_res[genes_ok, , drop = FALSE]
  # For each gene, the value with the largest magnitude among its resistant
  # lines (with its sign) -- that value decides the Up/Down ranking.
  extremo <- apply(lfc_ok, 1, function(v) v[which.max(abs(v))])
  up   <- sort(extremo[extremo > 0], decreasing = TRUE)
  down <- sort(extremo[extremo < 0], decreasing = FALSE)
  c(names(head(up, n_up)), names(head(down, n_down)))
}

lfc_alk1 <- calcular_lfc_1v1("ALK1")
lfc_alk2 <- calcular_lfc_1v1("ALK2")
lfc_alk3  <- calcular_lfc_1v1("ALK3")

# Deficit per model (WES/receptor panel genes absent from its TMM
# matrix) split half Up / half Down (ceiling/floor) -- with deficit=1 this
# adds exactly 1 extra gene.
n_extra_up_down <- function(lfc_mat) {
  deficit <- length(GENES_WES_RECEPTORES) - length(intersect(GENES_WES_RECEPTORES, rownames(lfc_mat)))
  c(up = ceiling(deficit / 2), down = floor(deficit / 2))
}
extra_alk1 <- n_extra_up_down(lfc_alk1)
extra_alk2 <- n_extra_up_down(lfc_alk2)
extra_alk3  <- n_extra_up_down(lfc_alk3)

top40_alk1 <- seleccionar_top_lfc(lfc_alk1, "ALK1", TOP_N_LFC + extra_alk1["up"], TOP_N_LFC + extra_alk1["down"])
top40_alk2 <- seleccionar_top_lfc(lfc_alk2, "ALK2", TOP_N_LFC + extra_alk2["up"], TOP_N_LFC + extra_alk2["down"])
top40_alk3  <- seleccionar_top_lfc(lfc_alk3,  "ALK3",  TOP_N_LFC + extra_alk3["up"],  TOP_N_LFC + extra_alk3["down"])

cat("   Top-LFC: ALK1", length(top40_alk1), "(+", extra_alk1["up"] + extra_alk1["down"], "for deficit) | ALK2", length(top40_alk2),
    "(+", extra_alk2["up"] + extra_alk2["down"], "for deficit) | ALK3", length(top40_alk3), "genes\n\n")

# =============================================================================
# 3. Heatmap per model -- union of top40 + 10 WES/receptors, rows
#    clustered WITHOUT splitting by block (a single block with
#    cluster_rows=TRUE, no row_split/row_gap between Top-LFC and WES/Receptor;
#    the Category colour strip is kept as an informative annotation, without a
#    physical gap). Columns =
#    treated parental(s) + resistant lines (see COLUMNAS), not
#    clustered -- they stay in the order of COLUMNAS[[modelo]].
# =============================================================================
cat("3. Generating individual heatmaps...\n")

construir_heatmap_modelo <- function(modelo, lfc_mat, top40, mostrar_leyenda = TRUE) {
  wes_disponibles <- intersect(GENES_WES_RECEPTORES, rownames(lfc_mat))
  faltan <- setdiff(GENES_WES_RECEPTORES, rownames(lfc_mat))
  if (length(faltan) > 0) {
    cat("   ", modelo, "- WES/receptor genes not found in the TMM matrix:", paste(faltan, collapse = ", "), "\n")
  }

  genes_finales <- union(top40, wes_disponibles)
  mat_final <- lfc_mat[genes_finales, , drop = FALSE]

  categoria <- factor(
    ifelse(genes_finales %in% top40, "Top LFC (Res vs Unt, per line)", "WES/Receptor panel"),
    levels = names(COLOR_CATEGORY)
  )

  left_ann <- rowAnnotation(
    Category = categoria,
    col = list(Category = COLOR_CATEGORY),
    show_annotation_name = FALSE,
    show_legend = mostrar_leyenda,
    simple_anno_size = unit(3, "mm")
  )

  cols_meta     <- COLUMNAS[[modelo]]
  idx_meta      <- match(colnames(mat_final), cols_meta$label)
  col_modelo    <- setNames(rep(modelo, ncol(mat_final)), colnames(mat_final))
  col_condicion <- setNames(cols_meta$condition[idx_meta], colnames(mat_final))
  col_drug      <- setNames(cols_meta$drug[idx_meta], colnames(mat_final))
  ann_top <- HeatmapAnnotation(
    Model     = col_modelo,
    Condition = col_condicion,
    Drug      = col_drug,
    col = list(Model = COLOR_MODELO, Condition = COLOR_CONDICION, Drug = COLOR_DRUG_FIRMA),
    show_annotation_name = FALSE,
    show_legend = mostrar_leyenda,
    simple_anno_size = unit(3, "mm")
  )

  ht <- Heatmap(
    mat_final,
    name = "LFC",
    col = COL_FUN_LFC,
    cluster_rows = TRUE,
    cluster_columns = FALSE,
    show_row_names = TRUE,
    row_names_gp = gpar(fontsize = 7, fontface = "italic"),
    column_names_gp = gpar(fontsize = 10, fontface = "bold"),
    column_names_rot = 0,
    top_annotation = ann_top,
    left_annotation = left_ann,
    show_heatmap_legend = mostrar_leyenda,
    rect_gp = gpar(col = "grey85", lwd = 0.5),
    heatmap_legend_param = list(title = "LFC\n(vs Untreated parental)",
                                 at = c(-LFC_MAX, 0, LFC_MAX)),
    column_title = paste0(modelo, " -- Treated parental + Resistant lines vs Untreated parental (", length(genes_finales), " genes)")
  )
  # genes_finales/top40 are also returned (not only the Heatmap object) so
  # the row order can be exported AFTER clustering (with
  # cluster_rows=TRUE the final order is only known after building the
  # Heatmap, row_order()).
  list(ht = ht, genes = genes_finales, top40 = top40, categoria = categoria)
}

# FIXED colour scale -8 to 8 instead of a dynamic scale (rounded to the
# real maximum, ~13), so that the 3 panels are comparable.
LFC_MAX <- 8
COL_FUN_LFC <- circlize::colorRamp2(c(-LFC_MAX, 0, LFC_MAX), GRADIENTE_RNA)

res_alk1 <- construir_heatmap_modelo("ALK1", lfc_alk1, top40_alk1)
res_alk2 <- construir_heatmap_modelo("ALK2", lfc_alk2, top40_alk2)
res_alk3  <- construir_heatmap_modelo("ALK3",  lfc_alk3,  top40_alk3)

resultados <- list(ALK1 = res_alk1, ALK2 = res_alk2, ALK3 = res_alk3)

# Figure caption: explains why the untreated
# parental is not a column of the heatmap, and that the treated
# parental IS (unlike the resistant lines, it does not take part in
# the Top LFC ranking, it is only shown as a context column).
PIE_FIGURA <- paste0(
  "LFC de cada columna (parental tratada de forma aguda + lineas resistentes) frente a la parental sin tratar del mismo modelo ",
  "(comparacion 1:1 por muestra, no promedio del grupo). La parental sin tratar no se muestra como columna: por definicion su LFC ",
  "frente a si misma seria 0 en todas las filas. La parental tratada SI se muestra (columna 'Treated', o 'Treated_L'/'Treated_A' en ",
  "ALK3 segun farmaco) pero no participa en la seleccion de genes Top LFC, solo aporta contexto de la respuesta aguda pre-resistencia. ",
  "Escala de color fija de -", LFC_MAX, " a ", LFC_MAX, ": los genes con |LFC|>", LFC_MAX, " (criterio de seleccion del bloque Top LFC) saturan al color extremo. ",
  "Categoria 'Top LFC': ", TOP_N_LFC, " genes mas sobreexpresados + ", TOP_N_LFC, " mas infraexpresados por magnitud de LFC entre las lineas resistentes de ese modelo ",
  "(con 1 gen extra en ALK1 y ALK2 para compensar que APOBEC3A no esta cuantificado en sus matrices TMM). ",
  "Categoria 'WES/Receptor panel': panel fijo de ", length(GENES_WES_RECEPTORES), " genes (", paste(GENES_WES_RECEPTORES, collapse = ", "), ")."
)
cat("Figure caption:\n ", PIE_FIGURA, "\n\n")
writeLines(PIE_FIGURA, file.path(OUT_DIR, "figures", "ALK_12_fig2_lfc_heatmap_pie.txt"))

for (modelo in names(resultados)) {
  out <- file.path(OUT_DIR, "figures", paste0("ALK_12_fig2_lfc_heatmap_", modelo, ".png"))
  png(out, width = ANCHO_HEATMAP, height = ALTO_HEATMAP, units = "in", res = 300)
  draw(resultados[[modelo]]$ht, heatmap_legend_side = "right", annotation_legend_side = "right")
  dev.off()
  cat("  ", out, "\n")
}

# Row order AFTER clustering -- row_order() is only known once the
# Heatmap is built (cluster_rows=TRUE), it is exported to CSV per model
# so that it does not have to be read off the PNG.
for (modelo in names(resultados)) {
  res <- resultados[[modelo]]
  orden_idx <- unlist(row_order(res$ht))
  genes_orden <- res$genes[orden_idx]
  categoria_orden <- as.character(res$categoria[orden_idx])
  df_orden <- data.frame(Rank = seq_along(genes_orden), Gene = genes_orden, Category = categoria_orden)
  out_orden <- file.path(OUT_DIR, "figures", paste0("ALK_12_fig2_lfc_heatmap_", modelo, "_gene_order.csv"))
  write.csv(df_orden, out_orden, row.names = FALSE)
  cat("  ", out_orden, "\n")
}

# =============================================================================
# 4. Combined panel -- the 3 heatmaps side by side, quick reference
#    (grid viewports, not patchwork -- ComplexHeatmap is not a ggplot grob,
#    same mechanism as ALK_11_individual_heatmaps.R)
# =============================================================================
cat("4. Generating combined panel...\n")

ht_alk1_comb <- construir_heatmap_modelo("ALK1", lfc_alk1, top40_alk1, mostrar_leyenda = FALSE)$ht
ht_alk2_comb <- construir_heatmap_modelo("ALK2", lfc_alk2, top40_alk2, mostrar_leyenda = FALSE)$ht
ht_alk3_comb  <- construir_heatmap_modelo("ALK3",  lfc_alk3,  top40_alk3,  mostrar_leyenda = FALSE)$ht

leyenda_compartida <- packLegend(
  Legend(col_fun = COL_FUN_LFC, title = "LFC (Resistant vs Untreated parental)",
         direction = "horizontal", legend_width = unit(4, "cm"), at = c(-LFC_MAX, 0, LFC_MAX)),
  Legend(labels = names(COLOR_CATEGORY), legend_gp = gpar(fill = COLOR_CATEGORY), title = "Category", direction = "horizontal"),
  Legend(labels = names(COLOR_MODELO), legend_gp = gpar(fill = COLOR_MODELO), title = "Model", direction = "horizontal"),
  Legend(labels = names(COLOR_CONDICION), legend_gp = gpar(fill = COLOR_CONDICION), title = "Condition", direction = "horizontal"),
  Legend(labels = names(COLOR_DRUG_FIRMA), legend_gp = gpar(fill = COLOR_DRUG_FIRMA), title = "Drug", direction = "horizontal"),
  direction = "horizontal"
)

out_combinado <- file.path(OUT_DIR, "figures", "ALK_12_fig2_lfc_heatmap_combinado.png")

# The 3 panels with the SAME width (same criterion as ALK_11) -- bottom row
# reserved for the single shared legend.
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

cat("\n=== ALK_12_lfc_heatmaps.R completed (3 individual PNGs + 1 combined) ===\n")
