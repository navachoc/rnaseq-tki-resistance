# ============================================================
# ALK_03_alk3_key_genes.R
# Purpose: key genes and data-driven resistance signature
#   (Upregulated/Downregulated on LFC_Res_vs_Trt) in COMBINED ALK3.
#   Cross-reference with alk_resistance_core and with the curated WES panel
#   (wes_panel_alk). Derived from ALK_03_alk2_key_genes.R.
#   The input table (02_lfc_master_table.csv) is long (1 row per
#   gene x drug); the signature keeps the Drug column to know in which arm
#   (Lorlatinib/Alectinib) each gene originated. 2-class scheme (see
#   ALK_02_alk3_lfc.R); Ya_Elevado_en_Treated is kept as informative
#   data, not as a separate class.
# Inputs (results/ALK/ALK3/tables/): 01_mat_tmm.csv, 01_metadata_unit.csv,
#   02_lfc_master_table.csv; config/colors.R, config/gene_lists.R
# Outputs (results/ALK/ALK3/): tables/ALK_resistance_signature*.csv,
#   tables/ALK_signature_wes_panel_overlap.csv, figures/03_*.png,
#   report/03_key_genes_report.md
# Usage: source(here::here("scripts/specific/ALK/ALK3/ALK_03_alk3_key_genes.R"))
#   or Rscript from any folder of the repo.
# ============================================================
CELL_LINE       <- "ALK3"
DRUG            <- "Combined"
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
DRUGS <- c("Lorlatinib", "Alectinib")

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

GRADIENTE_RNA <- c("#000000", "#808080", "#f2a900")
heatmap_expr  <- colorRampPalette(GRADIENTE_RNA)(100)

# Own drug colours for ALK3 -- LOCAL override, colors.R is not touched
# (shared with RET). Green=Lorlatinib, purple=Alectinib, same as GX_cluster in MET.
DRUG_COLORS_ALK <- c(Lorlatinib = "#00B894", Alectinib = "#6C5CE7", No = "#DFE6E9")

is_clean_gene <- function(g) !grepl("\\.", g) & !grepl("^[0-9]", g)

cat("\n================================================================\n")
cat("  KEY GENES —", OUT_PREFIX, "\n")
cat("================================================================\n\n")

# LOAD TABLES --------------------------------------------------------------
mat_file  <- file.path(RESULTS_DIR, "tables", "01_mat_tmm.csv")
meta_file <- file.path(RESULTS_DIR, "tables", "01_metadata_unit.csv")
lfc_file  <- file.path(RESULTS_DIR, "tables", "02_lfc_master_table.csv")

if (!file.exists(mat_file))  stop("01_mat_tmm.csv not found — run ALK_01_alk3_qc.R")
if (!file.exists(lfc_file))  stop("02_lfc_master_table.csv not found — run ALK_02_alk3_lfc.R")

mat_df        <- read.csv(mat_file, check.names = FALSE)
rownames(mat_df) <- mat_df$Gene
mat           <- as.matrix(mat_df[, -1])

meta_df       <- read.csv(meta_file, check.names = FALSE)
rownames(meta_df) <- meta_df$Model
metadata_unit <- meta_df[, -1]

lfc_table <- read.csv(lfc_file, check.names = FALSE)

parental_unt <- rownames(metadata_unit)[
  metadata_unit$Condition == "Parental" & metadata_unit$Treatment == "No"]
trt_by_drug <- setNames(lapply(DRUGS, function(d)
  rownames(metadata_unit)[metadata_unit$Condition == "Parental" & metadata_unit$Drug == d]), DRUGS)
resistant_by_drug <- setNames(lapply(DRUGS, function(d)
  rownames(metadata_unit)[metadata_unit$Condition == "Resistance" & metadata_unit$Drug == d]), DRUGS)
resistant_all <- unlist(resistant_by_drug, use.names = FALSE)
drug_of_model <- setNames(rep(names(resistant_by_drug), lengths(resistant_by_drug)), resistant_all)

cat("Tables loaded: 01_mat_tmm | 01_metadata_unit | 02_lfc_master_table\n\n")

# GENE SELECTION (both drugs, long table with Drug column) ----------
cat("1. Selecting key genes...\n")

up_genes   <- lfc_table %>% filter(Classification == "Upregulated") %>%
  arrange(desc(LFC_Res_vs_Trt))
down_genes <- lfc_table %>% filter(Classification == "Downregulated") %>%
  arrange(LFC_Res_vs_Trt)

cat("  Upregulated (both drugs):      ", nrow(up_genes), "\n", sep = "")
cat("  Downregulated (both drugs):    ", nrow(down_genes), "\n\n")

# CROSS-REFERENCE WITH alk_resistance_core -----------------------------------------------
cat("2. Cross-referencing with alk_resistance_core...\n")

core_up   <- intersect(up_genes$Gene,   alk_resistance_core)
core_down <- intersect(down_genes$Gene, alk_resistance_core)
cat("  alk_resistance_core — Upregulated:  ", paste(unique(core_up),   collapse = ", "), "\n")
cat("  alk_resistance_core — Downregulated:", paste(unique(core_down), collapse = ", "), "\n\n")

# CROSS-REFERENCE WITH wes_panel_alk (curated SNV/CNV/LOH panel of the same model) -----
cat("2b. Cross-referencing with wes_panel_alk (curated WES panel)...\n")

wes_up   <- intersect(up_genes$Gene,   wes_panel_alk)
wes_down <- intersect(down_genes$Gene, wes_panel_alk)
cat("  wes_panel_alk — Upregulated:  ", paste(unique(wes_up),   collapse = ", "), "\n")
cat("  wes_panel_alk — Downregulated:", paste(unique(wes_down), collapse = ", "), "\n\n")

# Annotation: Drug + Condition (Treatment is omitted, already implicit in Drug).
# Explicit Drug/Condition order -- pheatmap stacks the first column next
# to the heatmap and the following one outwards; Condition goes last to end up
# at the very top of the annotation bar.
annotation_col <- metadata_unit[, c("Drug", "Condition"), drop = FALSE]
ann_colors <- filter_annotation_colors(annotation_col)
ann_colors$Drug <- DRUG_COLORS_ALK

# RESISTANCE SIGNATURE (OUTPUT → scRNA-seq / oncoprint) -----------------------
cat("3. Generating resistance signature (data-driven)...\n")

# Upregulated + Downregulated, with the Drug column kept -- the same gene
# can appear twice (once per drug) if it is classified in both arms.
signature <- bind_rows(up_genes, down_genes) %>%
  mutate(
    In_alk_resistance_core = Gene %in% alk_resistance_core,
    In_wes_panel_alk       = Gene %in% wes_panel_alk
  )

# Forced genes (ALK/EML4) — add to the signature for each drug where they are
# not already present, keeping their real classification
forced_df <- lfc_table %>%
  filter(Gene %in% FORCED_GENES) %>%
  anti_join(signature, by = c("Gene", "Drug")) %>%
  mutate(
    In_alk_resistance_core = Gene %in% alk_resistance_core,
    In_wes_panel_alk       = Gene %in% wes_panel_alk
  )
if (nrow(forced_df) > 0) {
  cat("  Forced genes added to the signature (outside threshold):",
      paste0(forced_df$Gene, " [", forced_df$Drug, ": ", forced_df$Classification, "]", collapse = ", "), "\n")
  signature <- bind_rows(signature, forced_df)
}

write.csv(
  signature %>% filter(In_wes_panel_alk) %>%
    dplyr::select(Gene, Drug, Classification, LFC_Res_vs_Trt, LFC_Res_vs_Unt, LFC_Trt_vs_Unt, Ya_Elevado_en_Treated),
  file.path(RESULTS_DIR, "tables", paste0(sig_prefix, "_signature_wes_panel_overlap.csv")),
  row.names = FALSE
)
cat("  →", paste0(sig_prefix, "_signature_wes_panel_overlap.csv"), "(",
    length(unique(c(wes_up, wes_down))), "WES panel genes present in the signature, in 1 or both drugs)\n\n")

write.csv(signature,
          file.path(RESULTS_DIR, "tables", paste0(sig_prefix, "_resistance_signature.csv")),
          row.names = FALSE)
cat("  →", paste0(sig_prefix, "_resistance_signature.csv"), "(", nrow(signature), "gene x drug rows: ",
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
    "rows: ", sum(signature_cons$Classification == "Upregulated"), "up /",
    sum(signature_cons$Classification == "Downregulated"), "down)\n\n")

# FIGURE: SIGNATURE HEATMAP (all samples, annotated by Drug) ----------------
cat("4. Resistance signature heatmap...\n")

# Top25 up/down per drug, union -- same criterion as ALK_02_alk3_lfc.R
top_genes_drug <- function(drug, n = 25) {
  sig_d <- signature %>% filter(Drug == drug)
  up   <- sig_d %>% filter(Classification == "Upregulated")   %>% arrange(desc(LFC_Res_vs_Trt)) %>% slice_head(n = n) %>% pull(Gene)
  down <- sig_d %>% filter(Classification == "Downregulated") %>% arrange(LFC_Res_vs_Trt)       %>% slice_head(n = n) %>% pull(Gene)
  list(up = up, down = down)
}
lorla_sig <- top_genes_drug("Lorlatinib")
alec_sig  <- top_genes_drug("Alectinib")

sig_genes <- unique(c(lorla_sig$up, lorla_sig$down, alec_sig$up, alec_sig$down))
sig_genes <- sig_genes[sig_genes %in% rownames(mat) & is_clean_gene(sig_genes)]

if (length(sig_genes) >= 3) {
  forced_in_mat <- FORCED_GENES[FORCED_GENES %in% rownames(mat)]
  sig_ordered <- unique(c(forced_in_mat, sig_genes))

  z_sig <- t(apply(mat[sig_ordered, ], 1, cal_zscore))

  ph_sig <- pheatmap(z_sig,
           color             = heatmap_expr,
           annotation_col    = annotation_col,
           annotation_colors = ann_colors,
           cluster_rows      = TRUE,
           cluster_cols      = TRUE,
           clustering_method = "ward.D2",
           show_rownames     = TRUE,
           show_colnames     = TRUE,
           fontsize_row      = 7,
           fontsize_col      = 9,
           main = paste0(sig_prefix, " resistance signature — top50 por fármaco (unión) — ", OUT_PREFIX,
                         " (", length(sig_genes), " genes + ",
                         length(forced_in_mat), " forzados: ",
                         paste(forced_in_mat, collapse = ","), ")"),
           silent = TRUE)

  grid::grid.newpage()
  grid::grid.draw(ph_sig$gtable)
  save_pheatmap(ph_sig,
                file.path(RESULTS_DIR, "figures", "03_heatmap_signature.png"),
                width = 13, height = 15)
  cat("  → 03_heatmap_signature.png (", length(sig_genes), "genes, union Lorlatinib+Alectinib)\n\n")
} else {
  cat("  Fewer than 3 genes in the signature — heatmap skipped\n\n")
}

# FIGURE: FILTERED WES PANEL HEATMAP ------------------------------------------
cat("4b. Filtered WES panel heatmap (wes_panel_alk ∩ signature)...\n")

genes_wes_filtrado <- unique(c(wes_up, wes_down))
genes_wes_filtrado_forced <- union(FORCED_GENES, genes_wes_filtrado)
genes_wes_filtrado_forced <- genes_wes_filtrado_forced[genes_wes_filtrado_forced %in% rownames(mat)]
cat("  WES panel genes passing the filter (in 1 or both drugs):", length(genes_wes_filtrado),
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
top5_up   <- up_genes   %>% filter(Gene %in% rownames(mat)) %>% distinct(Gene, .keep_all = TRUE) %>%
  slice_head(n = 5) %>% pull(Gene)
top5_down <- down_genes %>% filter(Gene %in% rownames(mat)) %>% distinct(Gene, .keep_all = TRUE) %>%
  slice_head(n = 5) %>% pull(Gene)
vis_genes <- unique(c(FORCED_GENES[FORCED_GENES %in% rownames(mat)], top5_up, top5_down))

cat("  Stripplot/lineplot genes (top LFC + forced, both drugs combined):\n")
cat("    Forced:              ", paste(FORCED_GENES[FORCED_GENES %in% rownames(mat)], collapse = ", "), "\n")
cat("    Top5 Upregulated:    ", paste(top5_up,   collapse = ", "), "\n")
cat("    Top5 Downregulated:  ", paste(top5_down, collapse = ", "), "\n\n")

# FIGURE: STRIPPLOT ------------------------------------------------------------
if (length(vis_genes) > 0) {
  cat("5. Key genes stripplot...\n")

  sample_order <- c(parental_unt, trt_by_drug$Lorlatinib, trt_by_drug$Alectinib, resistant_all)
  strip_df <- lapply(vis_genes, function(g) {
    data.frame(
      Gene      = g,
      Sample    = factor(colnames(mat), levels = sample_order),
      Expr      = as.numeric(mat[g, ]),
      Condition = metadata_unit$Condition,
      Drug      = metadata_unit$Drug
    )
  }) %>% bind_rows() %>%
    mutate(Gene = factor(Gene, levels = vis_genes))

  p_strip <- ggplot(strip_df, aes(x = Sample, y = Expr, color = Condition, shape = Drug)) +
    geom_point(size = 3, alpha = 0.9) +
    scale_color_manual(values = condition_colors) +
    scale_shape_manual(values = c(Lorlatinib = 17, Alectinib = 15, No = 16)) +
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
      axis.text.x      = element_text(angle = 40, hjust = 1, size = 6),
      strip.text       = element_text(face = "bold.italic", size = 8),
      plot.title       = element_text(face = "bold"),
      plot.subtitle    = element_text(size = 7, color = "gray50"),
      panel.background = element_rect(fill = "white"),
      plot.background  = element_rect(fill = "white")
    )

  n_cols <- min(5, length(vis_genes))
  n_rows <- ceiling(length(vis_genes) / n_cols)
  ggsave(file.path(RESULTS_DIR, "figures", "03_stripplot_key_genes.png"),
         p_strip, width = n_cols * 2.7, height = n_rows * 3, dpi = 300)
  cat("  → 03_stripplot_key_genes.png (", length(vis_genes), "genes)\n\n")
}

# FIGURE: LINEPLOT --------------------------------------------------------------
# Trajectory per resistant line versus ITS OWN treated parental -- two dashed
# baselines (one per drug, colour=Drug), since Lorlatinib and Alectinib
# compare against different references (ALK3_L / ALK3_A).
if (length(vis_genes) > 0) {
  cat("6. Key genes lineplot (trajectory vs Treated, per drug)...\n")

  baseline_df <- bind_rows(lapply(DRUGS, function(d) {
    trt <- trt_by_drug[[d]]
    m   <- mat[vis_genes, trt, drop = FALSE]
    te  <- if (length(trt) > 1) rowMeans(m) else as.numeric(m)
    data.frame(Gene = vis_genes, Drug = d, Treated_Expr = te)
  })) %>% mutate(Gene = factor(Gene, levels = vis_genes))

  line_df <- lapply(resistant_all, function(s) {
    data.frame(
      Gene   = vis_genes,
      Sample = s,
      Drug   = drug_of_model[[s]],
      Expr   = as.numeric(mat[vis_genes, s])
    )
  }) %>% bind_rows() %>%
    mutate(
      Gene   = factor(Gene, levels = vis_genes),
      Sample = factor(Sample, levels = resistant_all)
    )

  p_line <- ggplot(line_df, aes(x = Sample, y = Expr, group = Drug, color = Drug)) +
    geom_hline(data = baseline_df, aes(yintercept = Treated_Expr, color = Drug),
               linetype = "dashed", linewidth = 0.5, inherit.aes = FALSE) +
    geom_line(linewidth = 0.6, alpha = 0.6) +
    geom_point(size = 3, alpha = 0.9) +
    facet_wrap(~ Gene, scales = "free_y", ncol = 5) +
    scale_color_manual(values = DRUG_COLORS_ALK) +
    labs(
      title    = paste("Key gene expression per resistant model —", OUT_PREFIX),
      subtitle = "Líneas punteadas = nivel Treated de cada fármaco (ALK3_L / ALK3_A)",
      x = NULL, y = "log2CPM (TMM)"
    ) +
    theme_minimal(base_size = 10) +
    theme(
      axis.text.x      = element_text(angle = 40, hjust = 1, size = 6),
      strip.text       = element_text(face = "bold.italic", size = 9),
      plot.title       = element_text(face = "bold"),
      plot.subtitle    = element_text(size = 7, color = "gray50"),
      panel.background = element_rect(fill = "white"),
      plot.background  = element_rect(fill = "white")
    )

  n_cols <- min(5, length(vis_genes))
  n_rows <- ceiling(length(vis_genes) / n_cols)
  ggsave(file.path(RESULTS_DIR, "figures", "03_lineplot_key_genes.png"),
         p_line, width = n_cols * 2.7, height = n_rows * 3, dpi = 300)
  cat("  → 03_lineplot_key_genes.png (", length(vis_genes), "genes)\n\n")
}

# REPORT --------------------------------------------------------------------
cat("7. Generating report...\n")

report_text <- paste0(
  "# Key Genes Report — ", OUT_PREFIX, "\n\n",
  "**Firma:** data-driven (sin lista de genes conocida) — Upregulated + Downregulated, ",
  "calculada por separado para Lorlatinib y Alectinib (tabla larga, columna Drug)\n\n",
  "## Top 15 Upregulated (ambos fármacos, ordenado por LFC)\n",
  paste(capture.output(print(
    as.data.frame(up_genes %>% slice_head(n = 15) %>%
                    dplyr::select(Gene, Drug, LFC_Res_vs_Trt, Ya_Elevado_en_Treated)),
    row.names = FALSE)), collapse = "\n"), "\n\n",
  "## alk_resistance_core detectados\n",
  "- Upregulated:  ", paste(unique(core_up),   collapse = ", "), "\n",
  "- Downregulated:", paste(unique(core_down), collapse = ", "), "\n\n",
  "## wes_panel_alk detectados (mismo panel curado de SNV/CNV/LOH)\n",
  "- Upregulated:  ", paste(unique(wes_up),   collapse = ", "), "\n",
  "- Downregulated:", paste(unique(wes_down), collapse = ", "), "\n",
  "- Ver también: ", paste0(sig_prefix, "_signature_wes_panel_overlap.csv"), " y 02_convergencia_farmacos.csv\n\n",
  paste0("## Firma ", sig_prefix, "_resistance_signature_conservative.csv\n"),
  "- N filas: ", nrow(signature_cons), " (",
  sum(signature_cons$Classification == "Upregulated"), " up / ",
  sum(signature_cons$Classification == "Downregulated"), " down)\n\n",
  "## Genes forzados (fusión driver)\n",
  "- ", paste(FORCED_GENES, collapse = ", "), " — siempre incluidos en la firma (ambos fármacos), el heatmap y ",
  "el stripplot/lineplot, pasen o no el filtro de LFC\n",
  if (nrow(forced_df) > 0)
    paste0("- Clasificación real (fuera de threshold): ",
           paste0(forced_df$Gene, " [", forced_df$Drug, ": ", forced_df$Classification, "]", collapse = ", "), "\n")
  else "- Ambos superan threshold de forma natural en ambos fármacos (ya incluidos sin forzar)\n",
  "\n",
  "## Genes stripplot/lineplot (forzados + top LFC combinado)\n",
  "- Forzados:          ", paste(FORCED_GENES[FORCED_GENES %in% rownames(mat)], collapse = ", "), "\n",
  "- Top5 Upregulated:  ", paste(top5_up,   collapse = ", "), "\n",
  "- Top5 Downregulated:", paste(top5_down, collapse = ", "), "\n\n",
  "## Figuras\n",
  "- 03_heatmap_signature.png   — firma completa (unión top25 up + top25 down por fármaco + ",
  paste(FORCED_GENES, collapse = "/"), " forzados; gradiente negro-gris-dorado; anotado por Drug)\n",
  "- 03_heatmap_wes_panel_filtrado.png — wes_panel_alk ∩ firma (", length(genes_wes_filtrado),
  " genes) + ", paste(FORCED_GENES, collapse = "/"), " forzados\n",
  "- 03_stripplot_key_genes.png — expresión por muestra (forma = fármaco)\n",
  "- 03_lineplot_key_genes.png  — trayectoria por modelo resistente vs nivel Treated de SU fármaco\n"
)

writeLines(report_text, file.path(RESULTS_DIR, "report", "03_key_genes_report.md"))
cat("  Report: report/03_key_genes_report.md\n\n")

cat("================================================================\n")
cat("  KEY GENES completed —", OUT_PREFIX, "\n")
cat("================================================================\n\n")
