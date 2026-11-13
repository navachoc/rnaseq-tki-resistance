# =============================================================================
# RET_09_figura2_lfc_heatmaps.R
#
# Purpose: differential expression figure for the RET models (RET1, RET2),
# RET version of the LFC heatmap in ALK_12_figura2_lfc_heatmaps.R. Three blocks:
#
#   1. LFC heatmap per model (sections 1-4). One column per displayed sample
#      (acutely treated parental + each resistant line), with 1 vs 1 LFC
#      against the untreated parental of the same model, computed from the TMM
#      matrix (same method as ALK_12). The untreated parental is not a column
#      (its LFC would be 0). Genes: 25 Up + 25 Down (TOP_N_LFC) by maximum LFC
#      magnitude across the resistant lines, excluding ENSG IDs without a
#      symbol, plus a minimal fixed panel of bypass genes (GENES_FORZADOS, see
#      section 2) that takes slots from the Up quota; 50 genes per model in
#      total. Fixed scale -LFC_MAX/+LFC_MAX (LFC_MAX = 6, same threshold as
#      ALK_12), black-grey-gold gradient (GRADIENTE_RNA).
#   2. GO BP + KEGG functional enrichment (section 5), on the FULL differential
#      expression signature (774 genes RET1 / 794 genes RET2), not on the 50
#      heatmap genes. Curated selection of 6 pathways per model
#      (RUTAS_GO_CURADAS), one per theme to avoid redundant GO terms, read
#      from the ORA tables already generated (clusterProfiler is not rerun).
#      RET2 has no significant Down pathway (neither GO BP nor KEGG), so it
#      only has Up pathways.
#   3. Fusion dotplot (section 6): Untreated -> Treated -> Resistant
#      trajectory of 2 fixed fusions per model (FUSIONES_PANEL_B): position 1 =
#      RET driver fusion (CCDC6--RET in RET1, KIF5B--RET in RET2), position 2
#      = FN1--ENSG00000228618 (RET1) / MYH9--CACNG2 (RET2). A fixed list is
#      used instead of an FFPM ranking because the ranking left CCDC6--RET out
#      in RET1. Size = JunctionReadCount, color = FFPM with GRADIENTE_RNA;
#      FFPM = 0 where the fusion is not detected.
#   Combined figure (section 7): LFC heatmaps + shared legend + enrichment +
#   fusion dotplot in a single PNG. The Heatmap objects are built once in
#   memory and reused, without intermediate .rds files.
#
# Inputs:
#   results/RET/{RET1,RET2}/tables/01_mat_tmm.csv
#   results/RET/{RET1,RET2}/tables/01_metadata_unit.csv
#   results/RET/{RET1,RET2}/tables/04_go_bp_{up,down}regulated.csv
#   results/RET/{RET1,RET2}/tables/04_kegg_{up,down}regulated.csv
#   results/RET/combined/tables/08_fusions_all_samples.csv (the STAR-Fusion
#     XLSX files are not reparsed here)
#
# Outputs (results/RET/combined/):
#   figures/RET_09_fig2_lfc_heatmap_{RET1,RET2}.png   (LFC heatmap per model)
#   figures/RET_09_fig2_lfc_heatmap_combinado.png      (LFC heatmaps side by side)
#   figures/RET_09_fig2_lfc_heatmap_{RET1,RET2}_gene_order.csv
#                                                      (row order after clustering)
#   figures/RET_09_fig2_lfc_heatmap_pie.txt            (descriptive text of the heatmap)
#   figures/RET_09_fig2_go_enrichment.png, tables/RET_09_fig2_go_enrichment.csv
#   figures/RET_09_fig2_fusions_dotplot.png, tables/RET_09_fig2_fusions_trajectory.csv
#   figures/RET_09_fig2_panel_combinado.png            (combined figure)
#
# Labels, axes and heatmap titles in English; color legend (title +
# categories) in Spanish where applicable.
#
# Dependencies (run first):
#   scripts/simple/rna_01.1_qc.R with DATASET = "RET1" and "RET2" (01_* tables)
#   scripts/specific/RET/rna_04b.1_pathways.R for RET1 and RET2 (04_* tables)
#   scripts/specific/RET/rna_08.1_fusions_ret.R (08_fusions_all_samples.csv)
#
# Usage: paths are resolved with here::here() from the repository root
# (.here file), so it works from any working directory:
#   Rscript scripts/specific/RET/RET_09_figura2_lfc_heatmaps.R
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
  library(ggplot2)
})

OUT_DIR <- here::here("results", "RET", "combined")
dir.create(file.path(OUT_DIR, "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT_DIR, "tables"),  recursive = TRUE, showWarnings = FALSE)

MODELOS <- c("RET1", "RET2")

GRADIENTE_RNA <- c("#000000", "#808080", "#f2a900")
COLOR_MODEL   <- c(RET1 = "#A55EEA", RET2 = "#F368E0")
COLOR_CONDICION <- c(Parental = "#e6a817", Resistente = "#013c58")
COLOR_DRUG      <- c(BLU667 = "#6C5CE7", LOXO292 = "#00B894")  # same hex as RET_08/rna_08.1_fusions_ret.R ("No" does not apply here, the untreated parental is not a column)

TOP_N_LFC <- 25  # per direction (Up/Down); 50 genes per model, including GENES_FORZADOS (section 2)

ANCHO_HEATMAP <- 8
ALTO_HEATMAP  <- 13  # same as ALK_12 (heatmap with readable gene names, not an unnamed top 500 as in RET_08)

LFC_MAX <- 6  # same threshold as ALK_12
COL_FUN_LFC <- circlize::colorRamp2(c(-LFC_MAX, 0, LFC_MAX), GRADIENTE_RNA)

# =============================================================================
# 1. Read TMM + metadata per model (LFC heatmap): same files as
#    RET_08_figura1_heatmaps_individuales.R (01_mat_tmm.csv / 01_metadata_
#    unit.csv); nothing is reanalyzed.
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
  read.csv(ruta, stringsAsFactors = FALSE)
}

tmm_ret1  <- leer_tmm("RET1")
tmm_ret2   <- leer_tmm("RET2")
meta_ret1 <- leer_metadata("RET1")
meta_ret2  <- leer_metadata("RET2")

MATS <- list(RET1 = tmm_ret1, RET2 = tmm_ret2)
METAS <- list(RET1 = meta_ret1, RET2 = meta_ret2)

# Bare ENSG IDs without a symbol are excluded from the top-LFC RANKING (same
# criterion as ALK_12): readthroughs with a compound symbol are valid, only the
# bare Ensembl ID is discarded.
es_ensg_puro <- function(g) grepl("^ENSG[0-9]+$", g)

# =============================================================================
# 2. 1 vs 1 LFC (each displayed column - untreated parental of the same model)
#    + selection of the 50 top-LFC genes (25 Up + 25 Down) by maximum
#    magnitude ACROSS THE RESISTANT LINES; same method as ALK_12 section 2.
#
# Minimal fixed panel (GENES_FORZADOS): the pure top-LFC ranking leaves out
# known bypass genes whose change is real and consistent but of moderate
# magnitude, whereas pseudogenes/ZNF/LINC without functional annotation
# dominate the ranking by raw magnitude. Genes verified against
# 01_mat_tmm.csv:
#  - RET2: *ALK* gains expression consistently in the 8 resistant lines
#    (LFC +0.7 to +3.5, highest in BR1), possible RTK bypass, absent from the
#    pure top 50.
#  - RET1: *AXL* and *FGFR1* both rise specifically in the 3 BLU-667
#    resistant lines (AXL +1.9 to +2.1, FGFR1 +1.1 to +1.4) and stay flat in
#    the 3 LOXO-292 ones, the cleanest bypass signature of the whole model,
#    also absent from the pure top 50.
# The 3 forced genes replace the 3 weakest members of the Up ranking (top-25
# boundary, no clear functional identity): PGM5 and BIVM.ERCC5 in RET1 (ranks
# 24-25 of the original Up), SNX32 in RET2 (rank 25 of the original Up, per
# RET_09_fig2_lfc_heatmap_RET2_gene_order.csv without forced genes). The total
# per model stays at 50 genes.
# =============================================================================
cat("2. Computing 1 vs 1 LFC (column - Untreated parental) and selecting genes...\n")

GENES_FORZADOS <- list(
  RET1 = c("AXL", "FGFR1"),
  RET2  = c("ALK")
)

calcular_lfc_1v1 <- function(modelo) {
  mat  <- MATS[[modelo]]
  meta <- METAS[[modelo]]
  untreated_sample <- meta$Sample[meta$Condition == "Parental" & meta$Treatment == "No"]
  # Displayed columns: acutely treated parental(s) + resistant lines. The
  # UNTREATED parental is never a column (its LFC against itself would be 0 in
  # all rows, same criterion as ALK).
  cols_meta <- meta[meta$Sample != untreated_sample, ]
  u   <- mat[, untreated_sample]
  lfc <- sapply(cols_meta$Sample, function(s) mat[, s] - u)
  rownames(lfc) <- rownames(mat)
  colnames(lfc) <- cols_meta$Sample
  list(lfc = lfc, cols_meta = cols_meta)
}

res_ret1 <- calcular_lfc_1v1("RET1")
res_ret2  <- calcular_lfc_1v1("RET2")
LFC_MATS  <- list(RET1 = res_ret1$lfc, RET2 = res_ret2$lfc)
COLS_META <- list(RET1 = res_ret1$cols_meta, RET2 = res_ret2$cols_meta)

seleccionar_top_lfc <- function(lfc_mat, cols_meta, forzados = character(0), n_up = TOP_N_LFC, n_down = TOP_N_LFC) {
  cols_resistentes <- cols_meta$Sample[cols_meta$Condition == "Resistance"]
  lfc_res  <- lfc_mat[, cols_resistentes, drop = FALSE]
  genes_ok <- rownames(lfc_res)[!es_ensg_puro(rownames(lfc_res))]
  # Forced genes are removed from the ranking pool (they must not compete or
  # appear duplicated); their slot in the Up quota is subtracted below (the 3
  # current forced genes are all Up).
  genes_ok <- setdiff(genes_ok, forzados)
  lfc_ok   <- lfc_res[genes_ok, , drop = FALSE]
  # For each gene, the value with the largest magnitude among its resistant
  # lines (keeping its sign); that value decides the Up/Down ranking.
  extremo <- apply(lfc_ok, 1, function(v) v[which.max(abs(v))])
  up   <- sort(extremo[extremo > 0], decreasing = TRUE)
  down <- sort(extremo[extremo < 0], decreasing = FALSE)
  c(forzados, names(head(up, n_up - length(forzados))), names(head(down, n_down)))
}

TOP_LFC <- list(
  RET1 = seleccionar_top_lfc(LFC_MATS[["RET1"]], COLS_META[["RET1"]], GENES_FORZADOS[["RET1"]]),
  RET2  = seleccionar_top_lfc(LFC_MATS[["RET2"]],  COLS_META[["RET2"]],  GENES_FORZADOS[["RET2"]])
)

cat("   Top-LFC: RET1", length(TOP_LFC[["RET1"]]), "genes | RET2", length(TOP_LFC[["RET2"]]), "genes\n\n")

# =============================================================================
# 3. LFC heatmap per model: the 50 top-LFC genes of that model, rows
#    clustered (cluster_rows=TRUE works fine here: each matrix belongs to a
#    single model with several real columns/lines, no NA). Columns = treated
#    parental(s) + resistant lines, NOT clustered (kept in the order of
#    COLS_META[[modelo]]). No gene category strip (left_annotation): the whole
#    block is a single category (top-LFC).
# =============================================================================
cat("3. Generating individual LFC heatmaps...\n")

construir_heatmap_modelo <- function(modelo, lfc_mat, cols_meta, top_lfc, mostrar_leyenda = TRUE) {
  mat_final <- lfc_mat[top_lfc, , drop = FALSE]

  idx_meta      <- match(colnames(mat_final), cols_meta$Sample)
  col_modelo    <- setNames(rep(modelo, ncol(mat_final)), colnames(mat_final))
  col_condicion <- setNames(ifelse(cols_meta$Condition[idx_meta] == "Resistance", "Resistente", "Parental"), colnames(mat_final))
  col_drug      <- setNames(cols_meta$Drug[idx_meta], colnames(mat_final))
  # Column label more readable than the raw sample name: drops the "MODEL_"
  # prefix (e.g. "RET1_BR1" -> "BR1"); acutely treated parentals
  # ("RET1_B250" -> "Treated_BLU667", according to Drug).
  etiqueta_col <- sub(paste0("^", modelo, "_"), "", colnames(mat_final))
  etiqueta_col <- ifelse(col_condicion == "Parental", paste0("Treated_", col_drug), etiqueta_col)
  colnames(mat_final) <- etiqueta_col
  names(col_modelo) <- names(col_condicion) <- names(col_drug) <- etiqueta_col

  ann_top <- HeatmapAnnotation(
    Modelo    = col_modelo,
    `Condición` = col_condicion,
    Inhibidor = col_drug,
    col = list(Modelo = COLOR_MODEL, `Condición` = COLOR_CONDICION, Inhibidor = COLOR_DRUG),
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
    show_heatmap_legend = mostrar_leyenda,
    rect_gp = gpar(col = "grey85", lwd = 0.5),
    heatmap_legend_param = list(title = "LFC\n(Resistente vs Parental sin tratar,\npor línea)",
                                 at = c(-LFC_MAX, 0, LFC_MAX)),
    column_title = paste0(modelo, " -- Treated parental + Resistant lines vs Untreated parental (", length(top_lfc), " genes)")
  )
  list(ht = ht, genes = top_lfc)
}

resultados <- list(
  RET1 = construir_heatmap_modelo("RET1", LFC_MATS[["RET1"]], COLS_META[["RET1"]], TOP_LFC[["RET1"]]),
  RET2  = construir_heatmap_modelo("RET2",  LFC_MATS[["RET2"]],  COLS_META[["RET2"]],  TOP_LFC[["RET2"]])
)

for (modelo in names(resultados)) {
  out <- file.path(OUT_DIR, "figures", paste0("RET_09_fig2_lfc_heatmap_", modelo, ".png"))
  png(out, width = ANCHO_HEATMAP, height = ALTO_HEATMAP, units = "in", res = 300)
  draw(resultados[[modelo]]$ht, heatmap_legend_side = "right", annotation_legend_side = "right")
  dev.off()
  cat("  ", out, "\n")
}

# Row order AFTER clustering, same criterion as ALK_12 (row_order() is only
# known once the Heatmap is built).
for (modelo in names(resultados)) {
  res <- resultados[[modelo]]
  orden_idx <- unlist(row_order(res$ht))
  genes_orden <- res$genes[orden_idx]
  df_orden <- data.frame(Rank = seq_along(genes_orden), Gene = genes_orden)
  out_orden <- file.path(OUT_DIR, "figures", paste0("RET_09_fig2_lfc_heatmap_", modelo, "_gene_order.csv"))
  write.csv(df_orden, out_orden, row.names = FALSE)
  cat("  ", out_orden, "\n")
}

# Descriptive text of the LFC heatmap (same spirit as ALK_12): explains why the
# untreated parental is not a column and that the LFC is 1 vs 1 per line.
PIE_FIGURA <- paste0(
  "LFC de cada columna (parental tratada de forma aguda + lineas resistentes) frente a la parental sin tratar del mismo modelo ",
  "(comparacion 1:1 por muestra, no promedio del grupo -- igual criterio que en ALK). La parental sin tratar no se muestra como ",
  "columna: por definicion su LFC frente a si misma seria 0 en todas las filas. La parental tratada SI se muestra (columna ",
  "'Treated_BLU667'/'Treated_LOXO292' segun farmaco) pero no participa en la seleccion de genes Top LFC, solo aporta contexto de la ",
  "respuesta aguda pre-resistencia. Escala de color fija de -", LFC_MAX, " a ", LFC_MAX, ": los genes con |LFC|>", LFC_MAX,
  " saturan al color extremo. Los ", TOP_N_LFC * 2, " genes de cada panel (", TOP_N_LFC, " mas sobreexpresados + ", TOP_N_LFC,
  " mas infraexpresados por magnitud de LFC entre las lineas resistentes de ese modelo) son mayoritariamente por ranking de LFC, ",
  "salvo un panel fijo minimo de genes de bypass conocidos (AXL y FGFR1 en RET1, ALK en RET2), ",
  "verificados a mano por su patron consistente pese a no alcanzar el top-LFC por magnitud bruta."
)
cat("Figure footnote (LFC heatmap):\n ", PIE_FIGURA, "\n\n")
writeLines(PIE_FIGURA, file.path(OUT_DIR, "figures", "RET_09_fig2_lfc_heatmap_pie.txt"))

# =============================================================================
# 4. Combined LFC heatmap (the 2 heatmaps side by side, quick reference);
#    same mechanism as ALK_12 section 4, adapted to 2 models instead of 3.
# =============================================================================
cat("4. Generating combined LFC heatmap (RET1 + RET2 side by side)...\n")

ht_ret1_comb <- construir_heatmap_modelo("RET1", LFC_MATS[["RET1"]], COLS_META[["RET1"]], TOP_LFC[["RET1"]], mostrar_leyenda = FALSE)$ht
ht_ret2_comb  <- construir_heatmap_modelo("RET2",  LFC_MATS[["RET2"]],  COLS_META[["RET2"]],  TOP_LFC[["RET2"]],  mostrar_leyenda = FALSE)$ht

leyenda_panel_a <- packLegend(
  Legend(col_fun = COL_FUN_LFC, title = "LFC (Resistente vs Parental sin tratar)",
         direction = "horizontal", legend_width = unit(4, "cm"), at = c(-LFC_MAX, 0, LFC_MAX)),
  Legend(labels = names(COLOR_MODEL), legend_gp = gpar(fill = COLOR_MODEL), title = "Modelo", direction = "horizontal"),
  Legend(labels = names(COLOR_CONDICION), legend_gp = gpar(fill = COLOR_CONDICION), title = "Condición", direction = "horizontal"),
  Legend(labels = names(COLOR_DRUG), legend_gp = gpar(fill = COLOR_DRUG), title = "Inhibidor", direction = "horizontal"),
  direction = "horizontal"
)

out_panel_a_combinado <- file.path(OUT_DIR, "figures", "RET_09_fig2_lfc_heatmap_combinado.png")
png(out_panel_a_combinado, width = 2 * ANCHO_HEATMAP, height = ALTO_HEATMAP + 1.5, units = "in", res = 300)
grid.newpage()
pushViewport(viewport(layout = grid.layout(2, 2,
  heights = unit(c(ALTO_HEATMAP, 1.5), c("inches", "inches")),
  widths  = unit(c(1, 1), "null"))))

pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
draw(ht_ret1_comb, newpage = FALSE)
upViewport()

pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
draw(ht_ret2_comb, newpage = FALSE)
upViewport()

pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1:2))
draw(leyenda_panel_a, just = "center")
upViewport()

dev.off()
cat("  ", out_panel_a_combinado, "\n\n")

# =============================================================================
# 5. Functional enrichment (GO BP + KEGG) per model, on the FULL differential
#    expression signature (NOT the 50 top-LFC heatmap genes; it is a different
#    analysis).
#    CURATED selection of 6 non-redundant pathways per model (the raw
#    clusterProfiler results include nearly identical GO BP terms, e.g.
#    "extracellular matrix organization" / "extracellular structure
#    organization" / "external encapsulating structure organization" are the
#    same signal under 3 different IDs; one per theme is chosen), read
#    directly from the ORA tables already generated by the RNA-seq pipeline
#    (04_go_bp_*.csv / 04_kegg_*.csv), without rerunning clusterProfiler here.
#    RET1: 3 Down (GO BP x2 + KEGG x1) + 3 Up (GO BP x3). RET2: 6 Up, NO Down
#    pathway: confirmed that none is significant (neither GO BP nor KEGG) in
#    the RET2 ORA tables; an Up/Down symmetry that does not exist in the data
#    is not forced.
# =============================================================================
cat("5. Computing GO BP + KEGG functional enrichment...\n")

RUTAS_GO_CURADAS <- list(
  RET1 = data.frame(
    ID        = c("GO:0097530", "GO:0050900", "hsa04060",
                   "GO:0030198", "GO:0006805", "GO:0007411"),
    Direccion = c("Down", "Down", "Down", "Up", "Up", "Up"),
    Fuente    = c("GO_BP", "GO_BP", "KEGG", "GO_BP", "GO_BP", "GO_BP")
  ),
  RET2 = data.frame(
    ID        = c("GO:0050919", "GO:0007157", "GO:0030198",
                   "GO:0007411", "GO:0006936", "GO:0001822"),
    Direccion = rep("Up", 6),
    Fuente    = rep("GO_BP", 6)
  )
)

leer_ora <- function(modelo, direccion, fuente) {
  archivo <- if (fuente == "GO_BP") paste0("04_go_bp_", tolower(direccion), "regulated.csv")
             else paste0("04_kegg_", tolower(direccion), "regulated.csv")
  ruta <- here::here("results", "RET", modelo, "tables", archivo)
  read.csv(ruta, stringsAsFactors = FALSE)
}

go_panel <- bind_rows(lapply(MODELOS, function(modelo) {
  curada <- RUTAS_GO_CURADAS[[modelo]]
  bind_rows(lapply(seq_len(nrow(curada)), function(i) {
    tabla <- leer_ora(modelo, curada$Direccion[i], curada$Fuente[i])
    fila  <- tabla[tabla$ID == curada$ID[i], c("ID", "Description", "GeneRatio", "FoldEnrichment", "p.adjust", "Count")]
    if (nrow(fila) != 1) {
      stop("Curated pathway ", curada$ID[i], " (", modelo, "/", curada$Direccion[i], "/", curada$Fuente[i],
           ") does not appear exactly once in its ORA table; check RUTAS_GO_CURADAS.")
    }
    fila$Modelo    <- modelo
    fila$Direccion <- curada$Direccion[i]
    fila$Fuente    <- curada$Fuente[i]
    fila
  }))
}))

go_panel <- go_panel %>%
  mutate(
    Modelo    = factor(Modelo, levels = MODELOS),
    Direccion = factor(Direccion, levels = c("Up", "Down"))
  ) %>%
  arrange(Modelo, desc(FoldEnrichment)) %>%
  mutate(Description = factor(Description, levels = rev(unique(Description))))

write.csv(go_panel, file.path(OUT_DIR, "tables", "RET_09_fig2_go_enrichment.csv"), row.names = FALSE)

cat("   GO pathways -- RET1:", paste(go_panel$Description[go_panel$Modelo == "RET1"], collapse = " | "),
    "\n   RET2:", paste(go_panel$Description[go_panel$Modelo == "RET2"], collapse = " | "), "\n\n")

# Dotplot, same style as ALK_09c_figura2_kegg_heatmap.R and the KEGG
# enrichment of ALK_09_figura3_egfr1_egfri.R: upward triangle (pch 24) =
# Upregulated, downward triangle (pch 25) = Downregulated; shape encodes
# direction. Color = p.adjust (terracotta = most significant / steel blue =
# least significant gradient), NOT direction, so color and shape are not
# redundant. Size = Count (number of genes in the term). One facet per model
# (scales="free_y", RET2 has no Down pathways), x axis = FoldEnrichment, as in
# ALK/EGFR.
ALTO_GO_PANEL <- 4.5
COLOR_KEGG_SIG   <- "#D6604D"  # most significant (terracotta); same hex as ALK/EGFR
COLOR_KEGG_NOSIG <- "#4393C3"  # least significant (steel blue); same hex as ALK/EGFR

p_go_panel <- ggplot(go_panel, aes(x = FoldEnrichment, y = Description,
                                    size = Count, fill = p.adjust, shape = Direccion)) +
  geom_point(color = "grey20", stroke = 0.3) +
  facet_wrap(~ Modelo, scales = "free_y") +
  scale_shape_manual(values = c(Up = 24, Down = 25), name = "Dirección",
                      labels = c(Up = "Sobreexpresado", Down = "Reprimido")) +
  scale_fill_gradient(low = COLOR_KEGG_SIG, high = COLOR_KEGG_NOSIG, trans = "log10",
                       name = "p ajustado", guide = guide_colorbar(reverse = TRUE)) +
  scale_size_continuous(range = c(3, 10), name = "Genes",
                         guide = guide_legend(override.aes = list(shape = 24, fill = "grey40"))) +
  labs(
    title = "Functional enrichment (GO BP + KEGG) on the full differential expression signature",
    x = "Fold enrichment", y = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.background = element_rect(fill = "white"),
    plot.background  = element_rect(fill = "white")
  )

out_panel_go <- file.path(OUT_DIR, "figures", "RET_09_fig2_go_enrichment.png")
ggsave(out_panel_go, p_go_panel, width = 2 * ANCHO_HEATMAP, height = ALTO_GO_PANEL, dpi = 300)
cat("  ", out_panel_go, "\n\n")

# =============================================================================
# 6. Fusion dotplot: trajectory of 2 fixed fusions per model: position 1 =
#    RET driver fusion, position 2 = specific gene (FN1 in RET1, MYH9 in
#    RET2). It is not an automatic FFPM ranking (see header). Dotplot style:
#    geom_point, size=JunctionReadCount, color=FFPM.
#    `08_fusions_all_samples.csv` is read (already aggregated from the
#    STAR-Fusion XLSX files by rna_08.1_fusions_ret.R); no XLSX is reparsed.
# =============================================================================
cat("6. Computing RET fusion + specific gene dotplot per model...\n")

RUTA_FUSIONES_CSV <- file.path(OUT_DIR, "tables", "08_fusions_all_samples.csv")
if (!file.exists(RUTA_FUSIONES_CSV)) {
  stop("Missing ", RUTA_FUSIONES_CSV, "; run scripts/specific/RET/rna_08.1_fusions_ret.R first.")
}
fusions_all <- read.csv(RUTA_FUSIONES_CSV, check.names = FALSE, stringsAsFactors = FALSE)

# Curated, fixed selection, not top-N by FFPM: position 1 is ALWAYS the RET
# driver fusion, position 2 a model-specific gene; 2 fusions per model to
# avoid overloading the figure.
FUSIONES_PANEL_B <- list(
  RET1 = c("CCDC6--RET", "FN1--ENSG00000228618"),
  RET2  = c("KIF5B--RET", "MYH9--CACNG2")
)

# Stage order per model: Untreated -> (Treated + Resistant) per drug; same as
# build_stage_order() in rna_08.1_fusions_ret.R, but reusing METAS[[modelo]]
# (already loaded in section 1 of this script, same 01_metadata_unit.csv file)
# instead of rereading metadata_*.csv.
orden_etapas_fusiones <- function(modelo) {
  m <- METAS[[modelo]]
  untreated <- m$Sample[m$Condition == "Parental" & m$Treatment == "No"]
  drugs     <- setdiff(unique(m$Drug), "No")
  ordenado  <- untreated
  for (d in drugs) {
    tratada    <- m$Sample[m$Condition == "Parental" & m$Drug == d]
    resistente <- sort(m$Sample[m$Condition == "Resistance" & m$Drug == d])
    ordenado   <- c(ordenado, tratada, resistente)
  }
  ordenado
}

ETAPAS_FUSIONES <- setNames(lapply(MODELOS, orden_etapas_fusiones), MODELOS)

TOP_FUSIONES <- bind_rows(lapply(MODELOS, function(modelo) {
  data.frame(`#FusionName` = FUSIONES_PANEL_B[[modelo]],
             Rank = seq_along(FUSIONES_PANEL_B[[modelo]]),
             Cell_Line = modelo,
             check.names = FALSE)
}))

# Trajectory of each of the 2 curated fusions: FFPM=0 in samples where it is
# not detected (real loss of the fusion), not NA; same criterion as
# rna_08.1_fusions_ret.R section 4.
driver_traj <- bind_rows(lapply(seq_len(nrow(TOP_FUSIONES)), function(i) {
  modelo  <- TOP_FUSIONES$Cell_Line[i]
  fname   <- TOP_FUSIONES$`#FusionName`[i]
  rank    <- TOP_FUSIONES$Rank[i]
  samples <- ETAPAS_FUSIONES[[modelo]]
  detectada <- fusions_all %>%
    filter(Cell_Line == modelo, `#FusionName` == fname) %>%
    group_by(Sample) %>%
    summarise(FFPM = sum(FFPM, na.rm = TRUE),
              JunctionReadCount = sum(JunctionReadCount, na.rm = TRUE),
              .groups = "drop")
  data.frame(Sample = samples) %>%
    left_join(detectada, by = "Sample") %>%
    mutate(
      FFPM = ifelse(is.na(FFPM), 0, FFPM),
      JunctionReadCount = ifelse(is.na(JunctionReadCount), 0, JunctionReadCount),
      Cell_Line = modelo,
      FusionName = fname,
      Rank = rank,
      Sample = factor(Sample, levels = samples)
    )
}))

# Facet label ("RET1: CCDC6--RET") ordered first by model (MODELOS order,
# RET1 before RET2) and then by Rank (1=RET, 2=specific gene, within each
# model); this level order is what lets facet_wrap(dir="v") place the 2 RET1
# facets in the left column and the 2 RET2 facets in the right column, with
# RET always in the top row.
orden_facetas <- TOP_FUSIONES %>%
  mutate(Cell_Line = factor(Cell_Line, levels = MODELOS)) %>%
  arrange(Cell_Line, Rank) %>%
  mutate(Facet_Label = paste0(Cell_Line, ": ", `#FusionName`))

driver_traj <- driver_traj %>%
  mutate(Facet_Label = factor(paste0(Cell_Line, ": ", FusionName),
                               levels = orden_facetas$Facet_Label))

write.csv(driver_traj, file.path(OUT_DIR, "tables", "RET_09_fig2_fusions_trajectory.csv"), row.names = FALSE)

cat("   Fusions -- RET1:", paste(TOP_FUSIONES$`#FusionName`[TOP_FUSIONES$Cell_Line == "RET1"], collapse = ", "),
    "| RET2:", paste(TOP_FUSIONES$`#FusionName`[TOP_FUSIONES$Cell_Line == "RET2"], collapse = ", "), "\n\n")

# Fusion dotplot: size=JunctionReadCount, color=FFPM. Color scale
# GRADIENTE_RNA (black-grey-gold), the same as the LFC heatmap, to unify the
# RNA scale across the whole figure. 2 rows x 2 columns grid (facet_wrap
# dir="v"): left column = RET + FN1 of RET1, right column = RET + MYH9 of
# RET2, RET always on top in both columns.
ALTO_FUSION_DOT <- 6.5  # 2 rows of facets

p_driver_dot <- ggplot(driver_traj, aes(x = Sample, y = 1,
                                          size = JunctionReadCount, color = FFPM)) +
  geom_point() +
  facet_wrap(~ Facet_Label, scales = "free_x", ncol = 2, dir = "v") +
  scale_color_gradientn(colors = GRADIENTE_RNA) +
  scale_size_continuous(range = c(1, 12)) +
  labs(
    title = "Fusiones RET + gen especifico -- evolución Untreated -> Treated -> Resistant",
    subtitle = "Izquierda: RET1 (CCDC6--RET, FN1)  |  Derecha: RET2 (KIF5B--RET, MYH9)",
    x = NULL, y = NULL, color = "FFPM", size = "Lecturas\nde unión"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.title  = element_text(face = "bold"),
    panel.grid.major.y = element_blank(),
    panel.background = element_rect(fill = "white"),
    plot.background  = element_rect(fill = "white")
  )

out_panel_b <- file.path(OUT_DIR, "figures", "RET_09_fig2_fusions_dotplot.png")
ggsave(out_panel_b, p_driver_dot, width = 2 * ANCHO_HEATMAP, height = ALTO_FUSION_DOT, dpi = 300)
cat("  ", out_panel_b, "\n\n")

# =============================================================================
# 7. Combined figure: LFC heatmaps (RET1/RET2 side by side, objects
#    ht_ret1_comb/ht_ret2_comb from section 4) + functional enrichment
#    (ggplot object p_go_panel from section 5) + RET + specific gene dotplot
#    (ggplot object p_driver_dot from section 6), in a single PNG, IN ADDITION
#    to keeping the individual PNGs (sections 3-6).
#
#    Layout: row 1 = the 2 LFC heatmaps (RET1 | RET2), row 2 = shared heatmap
#    legend, row 3 = functional enrichment, row 4 = RET + specific gene
#    dotplot (it already carries its own embedded color+size legend, so no
#    separate legend row is needed as for the heatmaps). The ggplots are
#    inserted with `print(..., newpage = FALSE)` inside their viewport, the
#    same grid mechanism the rest of the script uses to mix Heatmap objects.
# =============================================================================
cat("7. Generating combined figure (LFC heatmaps + GO/KEGG + RET + specific gene dotplot)...\n")

ANCHO_FIG_COMBINADA <- 2 * ANCHO_HEATMAP  # same as the width of the combined LFC heatmap (section 4)

out_png_combinado <- file.path(OUT_DIR, "figures", "RET_09_fig2_panel_combinado.png")
png(out_png_combinado, width = ANCHO_FIG_COMBINADA, height = ALTO_HEATMAP + 1.5 + ALTO_GO_PANEL + ALTO_FUSION_DOT, units = "in", res = 300)
grid.newpage()
pushViewport(viewport(layout = grid.layout(4, 2,
  heights = unit(c(ALTO_HEATMAP, 1.5, ALTO_GO_PANEL, ALTO_FUSION_DOT), c("inches", "inches", "inches", "inches")),
  widths  = unit(c(1, 1), "null"))))

pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
draw(ht_ret1_comb, newpage = FALSE)
upViewport()

pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
draw(ht_ret2_comb, newpage = FALSE)
upViewport()

pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1:2))
draw(leyenda_panel_a, just = "center")
upViewport()

pushViewport(viewport(layout.pos.row = 3, layout.pos.col = 1:2))
print(p_go_panel, newpage = FALSE)
upViewport()

pushViewport(viewport(layout.pos.row = 4, layout.pos.col = 1:2))
print(p_driver_dot, newpage = FALSE)
upViewport()

dev.off()
cat("  ", out_png_combinado, "\n")

cat("\n=== RET_09_figura2_lfc_heatmaps.R finished (2 individual LFC heatmaps + combined + row order CSV + footnote",
    "-- GO BP + KEGG functional enrichment (ggplot) -- RET + specific gene dotplot (ggplot)",
    "-- combined figure) ===\n")
