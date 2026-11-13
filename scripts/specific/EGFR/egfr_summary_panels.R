# =============================================================================
# egfr_summary_panels.R
#
# Purpose: EGFR expression panels (EGFR1 parental vs EGFR1_R resistant, bulk
# RNA-seq). Reads the CSV tables produced by the pipeline scripts (it does not
# redo TMM or LFC, to avoid duplicating logic or introducing inconsistencies).
# Ordered from general to specific:
#
# Panel A — Heatmap of the top 500 most variable genes (TMM log2CPM), same code
#   as rna_01.1_qc.R step 6. With only 2 samples (EGFR1/EGFR1_R) the per-gene
#   z-score is degenerate (each row ends up as a solid 2-colour block with no
#   real gradient), so with ncol(mat) <= 2 raw log2CPM is used with a
#   black-grey-gold gradient instead of z-score.
# Panel B — MA plot (Mean_Expr_Par vs LFC_Res_vs_Par), same code as
#   rna_02.1_lfc.R step "FIGURA 2" (2-class branch, Upregulated/
#   Downregulated/Other), key genes labelled from egfr_resistance_core with a
#   fixed label list. LOGFC_THRESHOLD = 2.0 must match the value used when
#   running rna_02.1_lfc.R for EGFR1_EGFRi (default value).
# Panel C — EGFR conservative signature, top 15 Up + top 15 Down by LFC, plus
#   the genes forced into the signature. The base is
#   EGFR_resistance_signature_conservative.csv (the SAME signature used to
#   project onto scRNA-seq), not egfr_resistance_core from gene_lists.R (a
#   hand-curated mechanism list that does not contain MDM2/TP53). The
#   signature already has its own Forced column (FGFR3/TACC3/MDM2/RB1, forced
#   by rna_03.1_key_genes.R regardless of the threshold). Discarded
#   alternatives: top15/15 within egfr_resistance_core (only 13 downregulated
#   genes) and top15/15 of the whole genome (dominated by genes with no
#   mechanistic anchor, CEACAM5/KLK6/ALPK2...). Barplot of LFC (Res vs Par),
#   sorted, gold/black colours to match the Panel A gradient (local change in
#   this script, NOT in colors.R).
# Panel C2 — EGFR resistance mechanisms (egfr_resistance_core, gene_lists.R),
#   additional panel to C: mechanism categories from Pros et al. /
#   Quintanal-Villalonga 2020 (SCLC/LuSCC transformation, amplification, EMT,
#   bypass, EGFR pathway), not verified gene by gene against the original
#   text. All genes in the list with an available LFC (no top15/15), standard
#   orange/blue colours (classification_colors).
# Panel D — KEGG enrichment (ORA, up/down combined on a single Y axis), same
#   visual style as ALK_09c_kegg_heatmap.R: upward triangle (pch 24) =
#   Upregulated, downward triangle (pch 25) = Downregulated, size = Count,
#   colour = p.adjust (terracotta #D6604D -> steel blue #4393C3 gradient).
#   Without the ALK "Convergence" column (EGFR1_EGFRi has a single comparator).
#   Top 5 Up + top 5 Down by p.adjust, read from the rna_04.1_pathways.R
#   tables (enrichKEGG is not repeated). If they are missing, the panel is
#   skipped with a console warning.
#
# Axes, titles and subtitles in English; colour legend (title + categories) in
# Spanish: Clasificacion/Direccion (Sobreexpresado/Infraexpresado/Otro),
# p ajustado. Z-score/log2CPM/Genes are not translated.
#
# Inputs (results/EGFR/EGFR1/tables/):
#   01_mat_tmm.csv, 01_metadata_unit.csv          (rna_01.1_qc.R)
#   02_lfc_master_table.csv                       (rna_02.1_lfc.R)
#   EGFR_resistance_signature_conservative.csv    (rna_03.1_key_genes.R)
#   04_kegg_upregulated.csv / 04_kegg_downregulated.csv (rna_04.1_pathways.R,
#     optional: if missing, panel D is skipped with a warning)
#   config/colors.R, config/gene_lists.R
#
# Outputs (results/EGFR/EGFR1/figures/):
#   fig3_panelA_heatmap_top500.png
#   fig3_panelB_maplot.png
#   fig3_panelC_signature_barplot.png
#   fig3_panelC2_mechanisms_barplot.png
#   fig3_panelD_kegg_updown.png
#
# Dependencies (run first, scripts/simple/): rna_01.1_qc.R
# (DATASET="EGFR1_EGFRi", TOP_N_HEATMAP=500), rna_02.1_lfc.R and
# rna_03.1_key_genes.R (DATASET="EGFR1_EGFRi"); for panel D also
# rna_04.1_pathways.R (DATASET="EGFR1_EGFRi", RECEPTOR="EGFR", CELL_LINE="EGFR1").
#
# Usage: paths are resolved with here::here() from the repository root
# (.here file), so it works from any working directory:
#   Rscript scripts/specific/EGFR/egfr_summary_panels.R
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
})

CONFIG_DIR  <- here::here("config")
RESULTS_DIR <- here::here("results", "EGFR", "EGFR1")
TAB_DIR     <- file.path(RESULTS_DIR, "tables")
FIG_DIR     <- file.path(RESULTS_DIR, "figures")

CELL_LINE <- "placeholder"
DRUG      <- "placeholder"
source(file.path(CONFIG_DIR, "colors.R"))
source(file.path(CONFIG_DIR, "gene_lists.R"))

TOP_N_HEATMAP <- 500   # same value as rna_01.1_qc.R for EGFR1_EGFRi
LOGFC_THRESHOLD <- 2.0   # same default value as rna_02.1_lfc.R

IN_MAT       <- file.path(TAB_DIR, "01_mat_tmm.csv")
IN_META      <- file.path(TAB_DIR, "01_metadata_unit.csv")
IN_LFC       <- file.path(TAB_DIR, "02_lfc_master_table.csv")
IN_SIG_CONS  <- file.path(TAB_DIR, "EGFR_resistance_signature_conservative.csv")
IN_KEGG_UP   <- file.path(TAB_DIR, "04_kegg_upregulated.csv")
IN_KEGG_DOWN <- file.path(TAB_DIR, "04_kegg_downregulated.csv")

for (f in c(IN_MAT, IN_META, IN_LFC, IN_SIG_CONS)) {
  if (!file.exists(f)) stop(
    "Missing ", f, " -- run rna_01.1_qc.R, rna_02.1_lfc.R and rna_03.1_key_genes.R (DATASET='EGFR1_EGFRi') first."
  )
}

save_pheatmap <- function(ph, filepath, width, height, res = 300,
                           legend_title = NULL, legend_range = NULL, legend_colors = NULL) {
  png(filepath, width = width, height = height, units = "in", res = res)
  if (!is.null(ph$gtable)) {
    grid::grid.draw(ph$gtable)
    # Hand-drawn colour scale, placed below the Treatment/Condition annotation
    # legend, with a title. pheatmap does not allow reordering its native
    # legends -- by default it places the colour scale and the annotation
    # legend in separate columns, side by side rather than stacked -- so the
    # native legend is disabled (legend = FALSE in the pheatmap call) and
    # replaced by this one, positioned below where Treatment/Condition falls.
    # Approximate position (leg_y_top/leg_y_bot); it may need fine tuning
    # after rendering if the number of annotation categories changes.
    if (!is.null(legend_title) && !is.null(legend_range) && !is.null(legend_colors)) {
      leg_x     <- 0.965
      leg_y_top <- 0.55
      leg_y_bot <- 0.30

      grid::grid.raster(
        matrix(rev(legend_colors), ncol = 1),
        x = leg_x, y = (leg_y_top + leg_y_bot) / 2,
        width = grid::unit(0.02, "npc"),
        height = grid::unit(leg_y_top - leg_y_bot, "npc"),
        interpolate = TRUE
      )
      grid::grid.text(legend_title, x = leg_x, y = leg_y_top + 0.03,
                       gp = grid::gpar(fontsize = 11, fontface = "bold"))
      grid::grid.text(sprintf("%.0f", legend_range[2]), x = leg_x + 0.02, y = leg_y_top,
                       just = "left", gp = grid::gpar(fontsize = 9))
      grid::grid.text(sprintf("%.0f", mean(legend_range)), x = leg_x + 0.02,
                       y = (leg_y_top + leg_y_bot) / 2, just = "left", gp = grid::gpar(fontsize = 9))
      grid::grid.text(sprintf("%.0f", legend_range[1]), x = leg_x + 0.02, y = leg_y_bot,
                       just = "left", gp = grid::gpar(fontsize = 9))
    }
  } else {
    grid::grid.newpage()
    grid::grid.draw(ph)
  }
  dev.off()
}

mat_df  <- read.csv(IN_MAT, check.names = FALSE)
meta_df <- read.csv(IN_META, check.names = FALSE)

# =============================================================================
# Panel A — Heatmap of the top 500 most variable genes
# =============================================================================
cat("Panel A — heatmap top", TOP_N_HEATMAP, "variable genes...\n")

mat <- as.matrix(mat_df[, -1, drop = FALSE])
rownames(mat) <- mat_df$Gene

metadata_unit <- meta_df[, -1, drop = FALSE]
rownames(metadata_unit) <- meta_df$Model
mat <- mat[, rownames(metadata_unit), drop = FALSE]

annotation <- metadata_unit[, !colnames(metadata_unit) %in% c("Sample", "Model", "Cell_Line"), drop = FALSE]
ann_colors <- filter_annotation_colors(annotation)

# Local Spanish translation of the pheatmap annotation legend --
# Treatment/Condition and their values (Yes/No/Parental/Resistance) come from
# 01_metadata_unit.csv and from the global `annotation_colors` list in
# config/colors.R, SHARED with MET/RET/ALK -- it is translated only here
# (columns + values + ann_colors of this script), without touching colors.R,
# so the other pipelines that already use those English names are unaffected.
VALORES_ANOT_ES  <- c(Yes = "Sí", No = "No", Parental = "Parental", Resistance = "Resistencia")
COLUMNAS_ANOT_ES <- c(Treatment = "Tratamiento", Condition = "Condición")
recode_es <- function(x) ifelse(x %in% names(VALORES_ANOT_ES), VALORES_ANOT_ES[x], x)

for (col in intersect(names(ann_colors), colnames(annotation))) {
  names(ann_colors[[col]]) <- recode_es(names(ann_colors[[col]]))
  annotation[[col]]        <- recode_es(as.character(annotation[[col]]))
}
names(ann_colors)    <- ifelse(names(ann_colors)    %in% names(COLUMNAS_ANOT_ES), COLUMNAS_ANOT_ES[names(ann_colors)],    names(ann_colors))
colnames(annotation) <- ifelse(colnames(annotation) %in% names(COLUMNAS_ANOT_ES), COLUMNAS_ANOT_ES[colnames(annotation)], colnames(annotation))

gene_var <- apply(mat, 1, var)
n_hm     <- min(TOP_N_HEATMAP, nrow(mat))
top_idx  <- order(gene_var, decreasing = TRUE)[seq_len(n_hm)]
fsize    <- ifelse(n_hm > 50, 6, 8)

# With only 2 samples (EGFR1/EGFR1_R) the per-gene z-score is degenerate -- each
# row has a single sign, a solid 2-colour block with no real gradient
# ("flag"). With ncol(mat) <= 2 raw log2CPM and a magnitude gradient are used
# instead of z-score (same criterion as rna_01.1_qc.R).
usar_zscore <- ncol(mat) > 2

if (usar_zscore) {
  cal_zscore <- function(x) if (sd(x) == 0) rep(0, length(x)) else (x - mean(x)) / sd(x)
  plot_mat     <- t(apply(mat[top_idx, ], 1, cal_zscore))
  color_top    <- heatmap_zscore
  titulo_leyenda <- "Z-score"
} else {
  plot_mat     <- mat[top_idx, ]
  color_top    <- colorRampPalette(c("#000000", "#808080", "#f2a900"))(100)
  titulo_leyenda <- "log2CPM"
}

# Treatment/Condition annotation bars above the heatmap WITH their category
# legend. annotation_names_col = FALSE removes the "Treatment"/"Condition"
# name attached to the right side of each colour bar.
ph_top <- pheatmap(plot_mat,
         color                = color_top,
         show_rownames        = (n_hm <= 100),
         show_colnames        = TRUE,
         fontsize_row         = fsize,
         annotation_col       = annotation,
         annotation_colors    = ann_colors,
         annotation_names_col = FALSE,
         legend               = FALSE,
         cluster_rows         = TRUE,
         cluster_cols         = TRUE,
         clustering_method    = "ward.D2",
         main              = paste0("Top ", n_hm, " variable genes — EGFR1 / EGFR1_R",
                                     if (!usar_zscore) " (log2CPM crudo, n<=2 muestras)" else ""),
         silent            = TRUE)

out_a <- file.path(FIG_DIR, "fig3_panelA_heatmap_top500.png")
# 10x9: with only 2 columns (EGFR1/EGFR1_R) a wider format leaves empty
# horizontal space; 10x9 leaves room for the dendrogram + legends on the right.
save_pheatmap(ph_top, out_a, width = 10, height = 9,
              legend_title = titulo_leyenda,
              legend_range = range(plot_mat, na.rm = TRUE),
              legend_colors = color_top)
cat("  ", out_a, "\n")

# =============================================================================
# Panel B — MA plot
# =============================================================================
cat("Panel B — MA plot...\n")

lfc_table <- read.csv(IN_LFC, check.names = FALSE)

key_genes <- egfr_resistance_core
top_up   <- lfc_table %>% filter(Classification == "Upregulated") %>%
  arrange(desc(LFC_Res_vs_Par)) %>% slice_head(n = 15) %>% pull(Gene)
top_down <- lfc_table %>% filter(Classification == "Downregulated") %>%
  arrange(LFC_Res_vs_Par) %>% slice_head(n = 15) %>% pull(Gene)
genes_to_label <- unique(c(
  intersect(key_genes, lfc_table$Gene[lfc_table$Classification != "Other"]),
  top_up, top_down
))
# Genes left unlabelled so that the large-font names do not overlap. The
# points are still drawn.
GENES_SIN_ETIQUETA <- c("ZNF804A", "ALDH3A1", "LINC00261", "AGR3", "TP63")
genes_to_label <- setdiff(genes_to_label, GENES_SIN_ETIQUETA)
# Fixed label list. With a finite max.overlaps, ggrepel dropped different
# labels depending on spacing (FGFR3 among them); fixing the list and
# max.overlaps = Inf always yields the same labels.
GENES_ETIQUETA <- c("CEACAM5", "CYP2B7P", "SFTA3.1", "SLC34A2", "SGPP2", "NKX2.1",
                    "HOPX", "CEACAM6", "SOX2", "KRT40", "KRT4", "FGFR3",
                    "PLCB4", "KLK10", "SIRPB1", "XG", "ROBO4", "CXorf49", "VWA5A",
                    "KLK5", "KLK6", "ALPK2")
genes_to_label <- intersect(genes_to_label, GENES_ETIQUETA)

plot_df <- lfc_table %>%
  mutate(
    Classification = factor(Classification, levels = c("Upregulated", "Downregulated", "Other")),
    label = ifelse(Gene %in% genes_to_label & Classification != "Other", Gene, "")
  ) %>%
  arrange(desc(Classification))

p_panelB <- ggplot(plot_df,
                    aes(Mean_Expr_Par, LFC_Res_vs_Par)) +
  geom_point(aes(color = Classification), alpha = 0.4, size = 0.7) +
  geom_hline(yintercept = c(-LOGFC_THRESHOLD, LOGFC_THRESHOLD),
             linetype = "dashed", color = "gray50", linewidth = 0.5) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.3) +
  # Names in black bold, like the MET MA plot (MET_01_summary_figures.R);
  # size 5 so they are legible.
  geom_text_repel(aes(label = label), size = 5, fontface = "bold", color = "black",
                  max.overlaps = Inf, segment.size = 0.2, box.padding = 0.5,
                  force = 2, max.iter = 20000, seed = 42) +
  # Legend in Spanish -- axes/titles in English.
  scale_color_manual(values = classification_colors, name = "Clasificación",
                     labels = c(Upregulated = "Sobreexpresado",
                                Downregulated = "Infraexpresado",
                                Other = "Otro")) +
  labs(
    title = "MA plot, EGFR1_R vs EGFR1",
    subtitle = paste0("|LFC| threshold = ", LOGFC_THRESHOLD, "  |  genes clave: egfr_resistance_core"),
    x = "Mean Parental log2CPM",
    y = "LFC Resistance vs Parental"
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(size = 8, color = "gray40"),
        legend.background = element_rect(fill = "white", color = NA),
        legend.key         = element_rect(fill = "white", color = NA))

out_b <- file.path(FIG_DIR, "fig3_panelB_maplot.png")
ggsave(out_b, p_panelB, width = 9, height = 7, dpi = 300)
cat("  ", out_b, "\n")

# =============================================================================
# Panel C — EGFR conservative signature (top 15 Up + top 15 Down + forced genes)
# =============================================================================
cat("Panel C — conservative signature (top 15 up + top 15 down)...\n")

# The base is NOT egfr_resistance_core (gene_lists.R, a hand-curated mechanism
# list that never contained MDM2/TP53), but
# EGFR_resistance_signature_conservative.csv, the SAME signature used to
# project onto scRNA-seq (rna_03.1_key_genes.R). That table already has its
# own Forced column -- FGFR3/TACC3/MDM2/RB1 are forced there by the pipeline
# itself (regardless of threshold), so MDM2 is recovered just by reading this
# file, with no new hand-made lists.
sig_cons <- read.csv(IN_SIG_CONS)

sig_up   <- sig_cons %>% filter(Classification == "Upregulated") %>%
  arrange(desc(LFC_Res_vs_Par)) %>% slice_head(n = 15)
sig_down <- sig_cons %>% filter(Classification == "Downregulated") %>%
  arrange(LFC_Res_vs_Par) %>% slice_head(n = 15)
sig_forced <- sig_cons %>% filter(Forced == TRUE)

df_c1 <- bind_rows(
    sig_up, sig_down,
    sig_forced %>% filter(!Gene %in% c(sig_up$Gene, sig_down$Gene))
  ) %>%
  distinct(Gene, .keep_all = TRUE) %>%
  mutate(Gene = factor(Gene, levels = Gene[order(LFC_Res_vs_Par)]))

genes_forzados_extra <- setdiff(sig_forced$Gene, c(sig_up$Gene, sig_down$Gene))
cat("  Genes: ", nrow(df_c1), " (15 up + 15 down + ", length(genes_forzados_extra),
    " forced genes outside the top 30 [", paste(genes_forzados_extra, collapse = ", "), "])\n", sep = "")

# Gold/black colours instead of the standard orange/blue of
# classification_colors -- so the barplot matches the Panel A heatmap palette
# (same black-grey-gold gradient). Local change to this panel, NOT in
# colors.R -- the rest of the project figures (ALK, Panel B MA plot, etc.)
# keep the standard orange/blue.
COLORES_PANEL_C <- c("Upregulated" = "#f2a900", "Downregulated" = "#000000")

p_panelC <- ggplot(df_c1, aes(x = LFC_Res_vs_Par, y = Gene, fill = Classification)) +
  geom_col(width = 0.7) +
  geom_vline(xintercept = 0, color = "gray40", linewidth = 0.3) +
  scale_fill_manual(values = COLORES_PANEL_C, name = "Dirección",
                     labels = c(Upregulated = "Sobreexpresado", Downregulated = "Infraexpresado")) +
  labs(
    title = "EGFR conservative resistance signature",
    subtitle = paste0("Top 15 Up + top 15 Down by LFC (EGFR_resistance_signature_conservative.csv)",
                       "  |  + genes included regardless of threshold: ", paste(sig_forced$Gene, collapse = ", ")),
    x = "log2FC (Res vs Par)", y = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 8, color = "gray40"),
    axis.text.y = element_text(face = "italic", size = 9),
    legend.position = "none"
  )

out_c <- file.path(FIG_DIR, "fig3_panelC_signature_barplot.png")
ggsave(out_c, p_panelC, width = 8, height = max(6, 1.5 + 0.28 * nrow(df_c1)), dpi = 300)
cat("  ", out_c, "\n")

# =============================================================================
# Panel C2 — EGFR resistance mechanisms (egfr_resistance_core)
# =============================================================================
cat("Panel C2 — egfr_resistance_core (Pros et al. / Quintanal-Villalonga 2020)...\n")

# Additional panel, does NOT replace Panel C -- same LFC table already read
# for Panel B, but filtered to the hand-curated mechanism genes
# (egfr_resistance_core, gene_lists.R: SCLC/LuSCC transformation,
# amplification, EMT, bypass, EGFR pathway), no top15/15 -- all genes with an
# available LFC are shown. Colours: standard orange/blue
# (classification_colors), see below.
df_c2 <- lfc_table %>%
  filter(Gene %in% egfr_resistance_core, !is.na(LFC_Res_vs_Par)) %>%
  distinct(Gene, .keep_all = TRUE) %>%
  mutate(
    Classification = if_else(LFC_Res_vs_Par >= 0, "Upregulated", "Downregulated"),
    Gene = factor(Gene, levels = Gene[order(LFC_Res_vs_Par)])
  )

genes_core_faltantes <- setdiff(egfr_resistance_core, df_c2$Gene)
if (length(genes_core_faltantes) > 0) {
  cat("  WARNING: egfr_resistance_core genes without LFC in the table: ",
      paste(genes_core_faltantes, collapse = ", "), "\n", sep = "")
}

p_panelC2 <- ggplot(df_c2, aes(x = LFC_Res_vs_Par, y = Gene, fill = Classification)) +
  geom_col(width = 0.7) +
  geom_vline(xintercept = 0, color = "gray40", linewidth = 0.3) +
  # Standard project orange/blue (classification_colors), different from the
  # gold/black of Panel C.
  scale_fill_manual(values = classification_colors[c("Upregulated", "Downregulated")], name = "Dirección",
                     labels = c(Upregulated = "Sobreexpresado", Downregulated = "Infraexpresado")) +
  labs(
    title = "EGFR resistance mechanisms (Pros et al. / Quintanal-Villalonga 2020)",
    subtitle = paste0("n = ", nrow(df_c2), " genes (egfr_resistance_core: SCLC/LuSCC transformation, ",
                       "amplification, EMT, bypass, EGFR pathway)  |  LFC EGFR1_R vs EGFR1"),
    x = "log2FC (Res vs Par)", y = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 8, color = "gray40"),
    axis.text.y = element_text(face = "italic", size = 9),
    legend.position = "none"
  )

out_c2 <- file.path(FIG_DIR, "fig3_panelC2_mechanisms_barplot.png")
ggsave(out_c2, p_panelC2, width = 8, height = max(6, 1.5 + 0.28 * nrow(df_c2)), dpi = 300)
cat("  ", out_c2, "\n")

# =============================================================================
# Panel D — KEGG enrichment (up/down), same style as ALK
# (ALK_09c_kegg_heatmap.R): upward triangle (pch 24) = Upregulated,
# downward triangle (pch 25) = Downregulated, size = Count. As in ALK, colour
# encodes p.adjust (terracotta -> steel blue gradient) and shape encodes
# direction; colour and shape are not redundant.
# Top 5 Up + top 5 Down by p.adjust, on a single Y axis (not faceted), like
# the ALK KEGG panels.
# =============================================================================
cat("Panel D — KEGG up/down...\n")

if (!file.exists(IN_KEGG_UP)) {
  cat("  WARNING: missing ", IN_KEGG_UP, " -- run rna_04.1_pathways.R (DATASET='EGFR1_EGFRi', ",
      "RECEPTOR='EGFR', CELL_LINE='EGFR1') first. Panel D skipped.\n", sep = "")
} else {
  kegg_up <- read.csv(IN_KEGG_UP, check.names = FALSE) %>% mutate(Direction = "Upregulated")
  kegg_down <- if (file.exists(IN_KEGG_DOWN)) {
    read.csv(IN_KEGG_DOWN, check.names = FALSE) %>% mutate(Direction = "Downregulated")
  } else {
    cat("  WARNING: missing ", IN_KEGG_DOWN, ", panel D with Upregulated only.\n", sep = "")
    NULL
  }

  kegg_top <- bind_rows(kegg_up, kegg_down) %>%
    group_by(Direction) %>%
    slice_min(p.adjust, n = 5, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(Direction = factor(Direction, levels = c("Upregulated", "Downregulated")))

  # A single combined Y axis (Up+Down together, ordered by FoldEnrichment),
  # same pattern as plot_panel_kegg_combinado() in ALK_09c.
  orden_desc <- kegg_top %>% arrange(FoldEnrichment) %>% pull(Description) %>% unique()
  kegg_top <- kegg_top %>% mutate(Description = factor(Description, levels = orden_desc))

  # Colour = p.adjust (terracotta->steel blue gradient, as in ALK), NOT
  # direction -- the shape (up/down triangle) already encodes direction.
  # Without the ALK "Convergence" text colouring of the Y axis (colouring the
  # pathway name by how many LFC comparators it was significant in) -- it does
  # not apply here, EGFR1_EGFRi has a single comparator (2 classes, Resistance vs
  # Parental), so there is no convergence to encode.
  COLOR_KEGG_SIG   <- "#D6604D"  # most significant (terracotta)
  COLOR_KEGG_NOSIG <- "#4393C3"  # least significant (steel blue)

  p_panelD <- ggplot(kegg_top, aes(x = FoldEnrichment, y = Description,
                                    size = Count, fill = p.adjust, shape = Direction)) +
    geom_point(color = "grey20", stroke = 0.3) +
    scale_shape_manual(values = c(Upregulated = 24, Downregulated = 25), name = "Dirección",
                       labels = c(Upregulated = "Sobreexpresado", Downregulated = "Infraexpresado")) +
    scale_fill_gradient(low = COLOR_KEGG_SIG, high = COLOR_KEGG_NOSIG, trans = "log10",
                         name = "p ajustado", guide = guide_colorbar(reverse = TRUE)) +
    scale_size_continuous(range = c(3, 10), name = "Genes",
                           guide = guide_legend(override.aes = list(shape = 24, fill = "grey40"))) +
    labs(
      title = "KEGG pathway enrichment (ORA)",
      subtitle = "Top 5 Up + top 5 Down (by p.adjust)  |  candidate seed for scRNA-seq signature screening (as in ALK)",
      x = "Fold enrichment", y = NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(size = 8.5, color = "gray40"),
      axis.text.y = element_text(size = 9),
      legend.background = element_rect(fill = "white", color = NA),
      legend.key         = element_rect(fill = "white", color = NA)
    )

  out_d <- file.path(FIG_DIR, "fig3_panelD_kegg_updown.png")
  ggsave(out_d, p_panelD, width = 10, height = max(6, 1.5 + 0.35 * nrow(kegg_top)), dpi = 300)
  cat("  ", out_d, "\n")
}

cat("\n=== egfr_summary_panels.R finished (panels A, B, C, C2, D generated here) ===\n")
