# ============================================================
# ALK_03_alk2_key_genes.R
# Purpose: key genes and data-driven resistance signature
#   (Upregulated/Downregulated on LFC_Res_vs_Trt) in ALK2. Cross-reference with
#   alk_resistance_core and with the curated WES panel (wes_panel_alk).
#   Dedicated per-model version, derived from scripts/simple/rna_03.1_key_genes.R.
#   2-class scheme (see ALK_02_alk2_lfc.R). The signature keeps the column
#   Ya_Elevado_en_Treated (inherited from the master table) as informative
#   data, not as a separate class. ALK and EML4 (driver fusion) are forced
#   into the signature and the figures.
# Inputs (results/ALK/ALK2/tables/): 01_mat_tmm.csv, 01_metadata_unit.csv,
#   02_lfc_master_table.csv; config/colors.R, config/gene_lists.R
# Outputs (results/ALK/ALK2/): tables/ALK_resistance_signature.csv,
#   tables/ALK_resistance_signature_conservative.csv (Mean_Expr_Unt > -5),
#   tables/ALK_signature_wes_panel_overlap.csv, figures/03_*.png,
#   report/03_key_genes_report.md
# Usage: source(here::here("scripts/specific/ALK/ALK2/ALK_03_alk2_key_genes.R"))
#   or Rscript from any folder of the repo.
# ============================================================
CELL_LINE       <- "ALK2"
DRUG            <- "Lorlatinib"
RECEPTOR        <- "ALK"
# ============================================================
RESULTS_DIR <- here::here("results", RECEPTOR, CELL_LINE)
CONFIG_DIR  <- here::here("config")
source(file.path(CONFIG_DIR, "colors.R"))
source(file.path(CONFIG_DIR, "gene_lists.R"))

OUT_PREFIX <- paste0(CELL_LINE, "_", DRUG)
sig_prefix <- RECEPTOR

# Genes forced into the signature and the figures of this script, whether or not they
# pass the LFC filter — ALK and EML4 are the two genes of the driver fusion.
FORCED_GENES <- c("ALK", "EML4")

# LIBRARIES ------------------------------------------------------------------
library(ggplot2)
library(pheatmap)
library(dplyr)

save_pheatmap <- function(ph, filepath, width, height, res = 300) {
  png(filepath, width = width, height = height, units = "in", res = res)
  if (!is.null(ph$gtable)) {
    grid::grid.draw(ph$gtable)
  } else {
    grid::grid.newpage()
    grid::grid.draw(ph)
  }
  dev.off()
}

cal_zscore <- function(x) if (sd(x) == 0) rep(0, length(x)) else (x - mean(x)) / sd(x)

# RNA-seq gradient (low -> high: black -> grey -> orange/gold), same
# as rna_03.1_key_genes.R (col_fun_firma).
GRADIENTE_RNA <- c("#000000", "#808080", "#f2a900")
heatmap_expr  <- colorRampPalette(GRADIENTE_RNA)(100)

is_clean_gene <- function(g) !grepl("\\.", g) & !grepl("^[0-9]", g)

cat("\n================================================================\n")
cat("  KEY GENES —", OUT_PREFIX, "\n")
cat("================================================================\n\n")

# LOAD TABLES --------------------------------------------------------------
mat_file  <- file.path(RESULTS_DIR, "tables", "01_mat_tmm.csv")
meta_file <- file.path(RESULTS_DIR, "tables", "01_metadata_unit.csv")
lfc_file  <- file.path(RESULTS_DIR, "tables", "02_lfc_master_table.csv")

if (!file.exists(mat_file))  stop("01_mat_tmm.csv not found — run ALK_01_alk2_qc.R")
if (!file.exists(lfc_file))  stop("02_lfc_master_table.csv not found — run ALK_02_alk2_lfc.R")

mat_df        <- read.csv(mat_file, check.names = FALSE)
rownames(mat_df) <- mat_df$Gene
mat           <- as.matrix(mat_df[, -1])

meta_df       <- read.csv(meta_file, check.names = FALSE)
rownames(meta_df) <- meta_df$Model
metadata_unit <- meta_df[, -1]

lfc_table <- read.csv(lfc_file, check.names = FALSE)

parental_unt <- rownames(metadata_unit)[
  metadata_unit$Condition == "Parental" & metadata_unit$Treatment == "No"]
parental_trt <- rownames(metadata_unit)[
  metadata_unit$Condition == "Parental" & metadata_unit$Treatment == "Yes"]
resistant <- rownames(metadata_unit)[metadata_unit$Condition == "Resistance"]

cat("Tables loaded: 01_mat_tmm | 01_metadata_unit | 02_lfc_master_table\n\n")

# GENE SELECTION ---------------------------------------------------------
cat("1. Selecting key genes...\n")

up_genes   <- lfc_table %>% filter(Classification == "Upregulated") %>%
  arrange(desc(LFC_Res_vs_Trt))
down_genes <- lfc_table %>% filter(Classification == "Downregulated") %>%
  arrange(LFC_Res_vs_Trt)

cat("  Upregulated:  ", nrow(up_genes), "\n", sep = "")
cat("  Downregulated:", nrow(down_genes), "\n\n")

# CROSS-REFERENCE WITH alk_resistance_core -----------------------------------------------
cat("2. Cross-referencing with alk_resistance_core...\n")

core_up   <- intersect(up_genes$Gene,   alk_resistance_core)
core_down <- intersect(down_genes$Gene, alk_resistance_core)
cat("  alk_resistance_core — Upregulated:  ", paste(core_up,   collapse = ", "), "\n")
cat("  alk_resistance_core — Downregulated:", paste(core_down, collapse = ", "), "\n\n")

# CROSS-REFERENCE WITH wes_panel_alk (curated SNV/CNV/LOH panel of the same model) -----
cat("2b. Cross-referencing with wes_panel_alk (curated WES panel)...\n")

wes_up   <- intersect(up_genes$Gene,   wes_panel_alk)
wes_down <- intersect(down_genes$Gene, wes_panel_alk)
cat("  wes_panel_alk — Upregulated:  ", paste(wes_up,   collapse = ", "), "\n")
cat("  wes_panel_alk — Downregulated:", paste(wes_down, collapse = ", "), "\n\n")

annotation_col <- metadata_unit[,
  !colnames(metadata_unit) %in% c("Sample", "Model", "Cell_Line"), drop = FALSE]
ann_colors <- filter_annotation_colors(annotation_col)

# RESISTANCE SIGNATURE (OUTPUT → scRNA-seq / oncoprint) -----------------------
cat("3. Generating resistance signature (data-driven)...\n")

signature <- bind_rows(up_genes, down_genes) %>%
  mutate(
    In_alk_resistance_core = Gene %in% alk_resistance_core,
    In_wes_panel_alk       = Gene %in% wes_panel_alk
  )

# Forced genes (ALK/EML4, driver fusion) — add them to the signature even if their
# real classification is Other (they do not exceed the LFC threshold)
forced_df <- lfc_table %>%
  filter(Gene %in% FORCED_GENES, !Gene %in% signature$Gene) %>%
  mutate(
    In_alk_resistance_core = Gene %in% alk_resistance_core,
    In_wes_panel_alk       = Gene %in% wes_panel_alk
  )
if (nrow(forced_df) > 0) {
  cat("  Forced genes added to the signature (outside threshold):",
      paste0(forced_df$Gene, " [", forced_df$Classification, "]", collapse = ", "), "\n")
  signature <- bind_rows(signature, forced_df)
}

write.csv(
  signature %>% filter(In_wes_panel_alk) %>%
    dplyr::select(Gene, Classification, LFC_Res_vs_Trt, LFC_Res_vs_Unt, LFC_Trt_vs_Unt, Ya_Elevado_en_Treated),
  file.path(RESULTS_DIR, "tables", paste0(sig_prefix, "_signature_wes_panel_overlap.csv")),
  row.names = FALSE
)
cat("  →", paste0(sig_prefix, "_signature_wes_panel_overlap.csv"), "(",
    length(c(wes_up, wes_down)), "WES panel genes present in the signature)\n\n")

write.csv(signature,
          file.path(RESULTS_DIR, "tables", paste0(sig_prefix, "_resistance_signature.csv")),
          row.names = FALSE)
cat("  →", paste0(sig_prefix, "_resistance_signature.csv"), "(", nrow(signature), "genes: ",
    sum(signature$Classification == "Upregulated"), "up /",
    sum(signature$Classification == "Downregulated"), "down, ",
    sum(signature$Ya_Elevado_en_Treated), "already elevated in Treated)\n")

# Conservative: expressed above noise in the untreated parental (baseline)
signature_cons <- signature %>% filter(Mean_Expr_Unt > -5)
write.csv(signature_cons,
          file.path(RESULTS_DIR, "tables",
                    paste0(sig_prefix, "_resistance_signature_conservative.csv")),
          row.names = FALSE)
cat("  →", paste0(sig_prefix, "_resistance_signature_conservative.csv"), "(", nrow(signature_cons),
    "genes: ", sum(signature_cons$Classification == "Upregulated"), "up /",
    sum(signature_cons$Classification == "Downregulated"), "down)\n\n")

# FIGURE: SIGNATURE HEATMAP -------------------------------------------------------
cat("4. Resistance signature heatmap...\n")

sig_genes <- signature$Gene
sig_genes <- sig_genes[sig_genes %in% rownames(mat) & is_clean_gene(sig_genes)]

if (length(sig_genes) >= 3) {
  sig_up_all   <- signature %>% filter(Classification == "Upregulated") %>%
    arrange(desc(LFC_Res_vs_Trt)) %>% pull(Gene)
  sig_down_all <- signature %>% filter(Classification == "Downregulated") %>%
    arrange(LFC_Res_vs_Trt) %>% pull(Gene)

  top25_up   <- head(sig_up_all[sig_up_all %in% sig_genes], 25)
  top25_down <- head(sig_down_all[sig_down_all %in% sig_genes], 25)
  forced_in_mat <- FORCED_GENES[FORCED_GENES %in% rownames(mat)]
  sig_ordered <- unique(c(forced_in_mat, top25_up, top25_down))

  z_sig <- t(apply(mat[sig_ordered, ], 1, cal_zscore))
  ann_sig_colors <- filter_annotation_colors(annotation_col)

  ph_sig <- pheatmap(z_sig,
           color             = heatmap_expr,
           annotation_col    = annotation_col,
           annotation_colors = ann_sig_colors,
           cluster_rows      = TRUE,
           cluster_cols      = TRUE,
           clustering_method = "ward.D2",
           show_rownames     = TRUE,
           show_colnames     = TRUE,
           fontsize_row      = 8,
           fontsize_col      = 9,
           main = paste0(sig_prefix, " resistance signature — top50 — ", OUT_PREFIX,
                         " (", length(top25_up), " up / ",
                         length(top25_down), " down / +",
                         length(forced_in_mat), " forzados: ",
                         paste(forced_in_mat, collapse = ","), ")"),
           silent = TRUE)

  grid::grid.newpage()
  grid::grid.draw(ph_sig$gtable)
  save_pheatmap(ph_sig,
                file.path(RESULTS_DIR, "figures", "03_heatmap_signature.png"),
                width = 12, height = 14)
  cat("  → 03_heatmap_signature.png (top50:",
      length(top25_up), "up /", length(top25_down), "down)\n\n")
} else {
  cat("  Fewer than 3 genes in the signature — heatmap skipped\n\n")
}

# FIGURE: FILTERED WES PANEL HEATMAP ------------------------------------------
# Genes of wes_panel_alk that DID pass the LFC filter (Classification !=
# "Other" in ALK_02), + forced ALK/EML4 -- visual version of wes_up/wes_down
# and of ALK_signature_wes_panel_overlap.csv computed above (step 2b).
cat("4b. Filtered WES panel heatmap (wes_panel_alk ∩ signature)...\n")

genes_wes_filtrado <- unique(c(wes_up, wes_down))
genes_wes_filtrado_forced <- union(FORCED_GENES, genes_wes_filtrado)
genes_wes_filtrado_forced <- genes_wes_filtrado_forced[genes_wes_filtrado_forced %in% rownames(mat)]
cat("  WES panel genes passing the filter:", length(genes_wes_filtrado),
    "(", paste(genes_wes_filtrado, collapse = ", "), ") + forced:",
    paste(FORCED_GENES, collapse = ", "), "\n")

if (length(genes_wes_filtrado_forced) >= 2) {
  z_wes <- t(apply(mat[genes_wes_filtrado_forced, , drop = FALSE], 1, cal_zscore))

  ph_wes <- pheatmap(z_wes,
           color             = heatmap_expr,
           show_rownames     = TRUE,
           show_colnames     = TRUE,
           fontsize_row      = 8,
           annotation_col    = annotation_col,
           annotation_colors = ann_colors,
           cluster_rows      = TRUE,
           cluster_cols      = TRUE,
           clustering_method = "ward.D2",
           main = paste0("WES curated panel (wes_panel_alk) — filtrado |LFC| >= threshold — ",
                         OUT_PREFIX, " (", length(genes_wes_filtrado_forced), " genes)"),
           silent = TRUE)

  grid::grid.newpage()
  grid::grid.draw(ph_wes$gtable)
  save_pheatmap(ph_wes,
                file.path(RESULTS_DIR, "figures", "03_heatmap_wes_panel_filtrado.png"),
                width = 10, height = 8)
  cat("  → 03_heatmap_wes_panel_filtrado.png (", length(genes_wes_filtrado_forced), "genes)\n\n")
} else {
  cat("  Fewer than 2 WES panel genes pass the filter — heatmap skipped\n\n")
}

# GENE SELECTION FOR STRIPPLOT / LINEPLOT (top LFC + forced ALK/EML4) ----
top5_up   <- up_genes   %>% filter(Gene %in% rownames(mat)) %>% slice_head(n = 5) %>% pull(Gene)
top5_down <- down_genes %>% filter(Gene %in% rownames(mat)) %>% slice_head(n = 5) %>% pull(Gene)
vis_genes <- unique(c(FORCED_GENES[FORCED_GENES %in% rownames(mat)], top5_up, top5_down))

cat("  Stripplot/lineplot genes (top LFC + forced):\n")
cat("    Forced:            ", paste(FORCED_GENES[FORCED_GENES %in% rownames(mat)], collapse = ", "), "\n")
cat("    Top5 Upregulated:  ", paste(top5_up,   collapse = ", "), "\n")
cat("    Top5 Downregulated:", paste(top5_down, collapse = ", "), "\n\n")

# FIGURE: STRIPPLOT ------------------------------------------------------------
if (length(vis_genes) > 0) {
  cat("5. Key genes stripplot...\n")

  sample_order <- c(parental_unt, parental_trt, resistant)
  strip_df <- lapply(vis_genes, function(g) {
    data.frame(
      Gene      = g,
      Sample    = factor(colnames(mat), levels = sample_order),
      Expr      = as.numeric(mat[g, ]),
      Condition = metadata_unit$Condition,
      Treatment = metadata_unit$Treatment
    )
  }) %>% bind_rows() %>%
    mutate(Gene = factor(Gene, levels = vis_genes))

  p_strip <- ggplot(strip_df, aes(x = Sample, y = Expr, color = Condition, shape = Treatment)) +
    geom_point(size = 3, alpha = 0.9) +
    scale_color_manual(values = condition_colors) +
    scale_shape_manual(values = c(Yes = 17, No = 16)) +
    facet_wrap(~ Gene, scales = "free_y", ncol = 5) +
    labs(
      title    = paste("Expression per sample —", OUT_PREFIX),
      subtitle = paste0("Forzados: ", paste(FORCED_GENES[FORCED_GENES %in% rownames(mat)], collapse = ", "),
                       " | Top Upregulated: ", paste(top5_up, collapse = ", "),
                       " | Top Downregulated: ", paste(top5_down, collapse = ", ")),
      x = NULL, y = "log2CPM (TMM)"
    ) +
    theme_minimal(base_size = 10) +
    theme(
      legend.position  = "right",
      axis.text.x      = element_text(angle = 40, hjust = 1, size = 7),
      strip.text       = element_text(face = "bold.italic", size = 8),
      plot.title       = element_text(face = "bold"),
      plot.subtitle    = element_text(size = 7, color = "gray50"),
      panel.background = element_rect(fill = "white"),
      plot.background  = element_rect(fill = "white")
    )

  n_cols <- min(5, length(vis_genes))
  n_rows <- ceiling(length(vis_genes) / n_cols)
  ggsave(file.path(RESULTS_DIR, "figures", "03_stripplot_key_genes.png"),
         p_strip, width = n_cols * 2.5, height = n_rows * 3, dpi = 300)
  cat("  → 03_stripplot_key_genes.png (", length(vis_genes), "genes)\n\n")
}

# FIGURE: LINEPLOT --------------------------------------------------------------
if (length(vis_genes) > 0) {
  cat("6. Key genes lineplot (trajectory vs Treated)...\n")

  treated_mat <- mat[vis_genes, parental_trt, drop = FALSE]
  treated_expr <- setNames(
    if (length(parental_trt) > 1) rowMeans(treated_mat) else as.numeric(treated_mat),
    vis_genes
  )

  line_df <- lapply(vis_genes, function(g) {
    data.frame(
      Gene         = g,
      Sample       = resistant,
      Expr         = as.numeric(mat[g, resistant]),
      Treated_Expr = treated_expr[g]
    )
  }) %>% bind_rows() %>%
    mutate(
      Gene   = factor(Gene, levels = vis_genes),
      Sample = factor(Sample, levels = resistant)
    )

  p_line <- ggplot(line_df, aes(x = Sample, y = Expr, group = Gene)) +
    geom_hline(aes(yintercept = Treated_Expr), linetype = "dashed",
               color = "gray60", linewidth = 0.5) +
    geom_line(linewidth = 0.6, alpha = 0.5, color = "#013c58") +
    geom_point(size = 3, alpha = 0.9, color = "#013c58") +
    facet_wrap(~ Gene, scales = "free_y", ncol = 5) +
    labs(
      title    = paste("Key gene expression per resistant model —", OUT_PREFIX),
      subtitle = "Dashed line = Treated parental level",
      x = NULL, y = "log2CPM (TMM)"
    ) +
    theme_minimal(base_size = 10) +
    theme(
      axis.text.x      = element_text(angle = 40, hjust = 1, size = 7),
      strip.text       = element_text(face = "bold.italic", size = 9),
      plot.title       = element_text(face = "bold"),
      panel.background = element_rect(fill = "white"),
      plot.background  = element_rect(fill = "white")
    )

  n_cols <- min(5, length(vis_genes))
  n_rows <- ceiling(length(vis_genes) / n_cols)
  ggsave(file.path(RESULTS_DIR, "figures", "03_lineplot_key_genes.png"),
         p_line, width = n_cols * 2.5, height = n_rows * 3, dpi = 300)
  cat("  → 03_lineplot_key_genes.png (", length(vis_genes), "genes)\n\n")
}

# REPORT --------------------------------------------------------------------
cat("7. Generating report...\n")

report_text <- paste0(
  "# Key Genes Report — ", OUT_PREFIX, "\n\n",
  "**Firma:** data-driven (sin lista de genes conocida) — Upregulated + Downregulated ",
  "(Ya_Elevado_en_Treated se conserva como columna informativa, no como clase)\n\n",
  "## Top 15 Upregulated\n",
  paste(capture.output(print(
    as.data.frame(slice_head(up_genes, n = 15) %>%
                    dplyr::select(Gene, LFC_Res_vs_Trt, Ya_Elevado_en_Treated)),
    row.names = FALSE)), collapse = "\n"), "\n\n",
  "## alk_resistance_core detectados\n",
  "- Upregulated:  ", paste(core_up,   collapse = ", "), "\n",
  "- Downregulated:", paste(core_down, collapse = ", "), "\n\n",
  "## wes_panel_alk detectados (mismo panel curado de SNV/CNV/LOH)\n",
  "- Upregulated:  ", paste(wes_up,   collapse = ", "), "\n",
  "- Downregulated:", paste(wes_down, collapse = ", "), "\n",
  "- Ver también: ", paste0(sig_prefix, "_signature_wes_panel_overlap.csv"), "\n\n",
  paste0("## Firma ", sig_prefix, "_resistance_signature_conservative.csv\n"),
  "- N genes: ", nrow(signature_cons), " (",
  sum(signature_cons$Classification == "Upregulated"), " up / ",
  sum(signature_cons$Classification == "Downregulated"), " down)\n\n",
  "## Genes forzados (fusión driver)\n",
  "- ", paste(FORCED_GENES, collapse = ", "), " — siempre incluidos en la firma, el heatmap y ",
  "el stripplot/lineplot, pasen o no el filtro de LFC\n",
  if (nrow(forced_df) > 0)
    paste0("- Clasificación real (fuera de threshold): ",
           paste0(forced_df$Gene, " [", forced_df$Classification, "]", collapse = ", "), "\n")
  else "- Ambos superan threshold de forma natural (ya incluidos sin forzar)\n",
  "\n",
  "## Genes stripplot/lineplot (forzados + top LFC)\n",
  "- Forzados:          ", paste(FORCED_GENES[FORCED_GENES %in% rownames(mat)], collapse = ", "), "\n",
  "- Top5 Upregulated:  ", paste(top5_up,   collapse = ", "), "\n",
  "- Top5 Downregulated:", paste(top5_down, collapse = ", "), "\n\n",
  "## Figuras\n",
  "- 03_heatmap_signature.png   — firma completa (top25 up + top25 down + ",
  paste(FORCED_GENES, collapse = "/"), " forzados; gradiente negro-gris-dorado)\n",
  "- 03_heatmap_wes_panel_filtrado.png — wes_panel_alk ∩ firma (", length(genes_wes_filtrado),
  " genes) + ", paste(FORCED_GENES, collapse = "/"), " forzados\n",
  "- 03_stripplot_key_genes.png — expresión por muestra\n",
  "- 03_lineplot_key_genes.png  — trayectoria por modelo resistente vs nivel Treated\n"
)

writeLines(report_text, file.path(RESULTS_DIR, "report", "03_key_genes_report.md"))
cat("  Report: report/03_key_genes_report.md\n\n")

cat("================================================================\n")
cat("  KEY GENES completed —", OUT_PREFIX, "\n")
cat("================================================================\n\n")
