# ============================================================
# rna_02.1_lfc.R
# log2 fold change (LFC) calculation on TMM log2CPM and gene classification.
# Simplified scheme (Upregulated/Downregulated/Other) for the datasets in
# TWO_CLASS_DATASETS (no treated parental); 5-class scheme (Acquired,
# Pre_existing, Maintained, Downregulated, Other) for the rest.
#
# Inputs:    results/<RECEPTOR>/<CELL_LINE>/tables/01_mat_tmm.csv,
#              01_metadata_unit.csv (output of rna_01.1_qc.R)
#            config/colors.R, config/gene_lists.R
# Outputs:   results/<RECEPTOR>/<CELL_LINE>/tables/02_lfc_master_table.csv,
#              02_classified_genes.csv (input of rna_03 and rna_04)
#            results/<RECEPTOR>/<CELL_LINE>/figures/02_*.png
#            results/<RECEPTOR>/<CELL_LINE>/report/02_lfc_report.{md,html}
# Usage:     set the CONFIG block and, from any folder of the repo,
#            source(here::here("scripts/simple/rna_02.1_lfc.R"))
# ============================================================
# CONFIG — change here for each run
# ============================================================
DATASET         <- "MET1_METi"   # "ALK12" | "ALK3" | "MET1_METi" | "EGFR1_EGFRi"
CELL_LINE       <- "MET1"        # "ALK2" | "ALK1" | "ALK3" | "MET1" | "EGFR1"
DRUG            <- "MET_inhibitor" # "Lorlatinib" | "Alectinib" | "MET_inhibitor" | "Gefitinib_Afatinib"
RECEPTOR        <- "MET"         # "MET" | "ALK" | "EGFR"
LOGFC_THRESHOLD <- 2.0
# ============================================================
RESULTS_DIR <- here::here("results", RECEPTOR, CELL_LINE)
CONFIG_DIR  <- here::here("config")
source(file.path(CONFIG_DIR, "colors.R"))
source(file.path(CONFIG_DIR, "gene_lists.R"))

# Category names of the 5-class scheme (datasets outside TWO_CLASS_DATASETS)
CLASS_LABELS <- c(down = "Downregulated", pre = "Pre_existing",
                   maint = "Maintained", acq = "Acquired")

# Datasets with the simplified scheme (Resistance vs Parental, no treated/untreated distinction)
TWO_CLASS_DATASETS <- c("MET1_METi", "EGFR1_EGFRi")

# LIBRARIES ------------------------------------------------------------------
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(dplyr)
library(tibble)

save_pheatmap <- function(ph, filepath, width, height, res = 300) {
  png(filepath, width = width, height = height, units = "in", res = res)
  if (methods::is(ph, "Heatmap") || methods::is(ph, "HeatmapList")) {
    ComplexHeatmap::draw(ph)
  } else if (!is.null(ph$gtable)) {
    grid::grid.draw(ph$gtable)
  } else {
    grid::grid.newpage()
    grid::grid.draw(ph)
  }
  dev.off()
}

detect_sep <- function(file) {
  if (grepl(";", readLines(file, n = 1))) ";" else ","
}

cat("\n================================================================\n")
cat("  LFC —", OUT_PREFIX, "\n")
cat("================================================================\n\n")

# LOAD TABLES FROM rna_01 ----------------------------------------------------
mat_file  <- file.path(RESULTS_DIR, "tables", "01_mat_tmm.csv")
meta_file <- file.path(RESULTS_DIR, "tables", "01_metadata_unit.csv")

if (!file.exists(mat_file)) stop("01_mat_tmm.csv not found — run rna_01.1_qc.R first")

mat_df        <- read.csv(mat_file, check.names = FALSE)
rownames(mat_df) <- mat_df$Gene
mat           <- as.matrix(mat_df[, -1])

meta_df       <- read.csv(meta_file, check.names = FALSE)
rownames(meta_df) <- meta_df$Model
metadata_unit <- meta_df[, -1]

if (DATASET == "MET1_METi" && !"GX_cluster" %in% colnames(metadata_unit)) {
  rg_map <- c(
    MET1    = "Parental",
    MET1_C1 = "Group1_EGFR_bypass", MET1_PR1 = "Group1_EGFR_bypass",
    MET1_PR2 = "Group1_EGFR_bypass", MET1_V8  = "Group1_EGFR_bypass",
    MET1_PR3 = "Group2_NF2_YAP",    MET1_R = "Group2_NF2_YAP"
  )
  metadata_unit$GX_cluster <- rg_map[rownames(metadata_unit)]
}

cat("Tables 01_mat_tmm.csv | 01_metadata_unit.csv loaded\n\n")

# IDENTIFY GROUPS ---------------------------------------------------------
cat("1. Identifying groups...\n")

parental  <- rownames(metadata_unit)[metadata_unit$Condition == "Parental"]
resistant <- rownames(metadata_unit)[metadata_unit$Condition == "Resistance"]

cat("  Parental: ", paste(parental,  collapse = ", "), "\n")
cat("  Resistant:", paste(resistant, collapse = ", "), "\n\n")

if (length(parental)  == 0) stop("No parental sample found.")
if (length(resistant) == 0) stop("No resistant models found.")

if (!DATASET %in% TWO_CLASS_DATASETS) {
  parental_unt <- rownames(metadata_unit)[
    metadata_unit$Condition == "Parental" & metadata_unit$Treatment == "No"
  ]
  parental_trt <- rownames(metadata_unit)[
    metadata_unit$Condition == "Parental" & metadata_unit$Treatment == "Yes"
  ]
  if (length(parental_unt) == 0) stop("Untreated parental not found.")
  if (length(parental_trt) == 0) stop("Treated parental not found.")
  cat("  Parental untreated:", parental_unt, "\n")
  cat("  Parental treated:  ", parental_trt, "\n\n")
}

# LFC -------------------------------------------------------------------------
cat("2. Computing LFC...\n")

expr_par <- if (length(parental) > 1) rowMeans(mat[, parental]) else mat[, parental]
expr_res <- if (length(resistant) > 1) rowMeans(mat[, resistant]) else mat[, resistant]

lfc_res_vs_par <- expr_res - expr_par   # main comparison

if (!DATASET %in% TWO_CLASS_DATASETS) {
  expr_unt       <- if (length(parental_unt) > 1) rowMeans(mat[, parental_unt]) else mat[, parental_unt]
  expr_trt       <- if (length(parental_trt) > 1) rowMeans(mat[, parental_trt]) else mat[, parental_trt]
  lfc_res_vs_trt <- expr_res - expr_trt
  lfc_res_vs_unt <- expr_res - expr_unt
  lfc_trt_vs_unt <- expr_trt - expr_unt

  cat("  lfc_res_vs_trt: range [",
      round(min(lfc_res_vs_trt), 2), ",",
      round(max(lfc_res_vs_trt), 2), "]\n")
  cat("  lfc_res_vs_unt: range [",
      round(min(lfc_res_vs_unt), 2), ",",
      round(max(lfc_res_vs_unt), 2), "]\n\n")
} else {
  cat("  lfc_res_vs_par: range [",
      round(min(lfc_res_vs_par), 2), ",",
      round(max(lfc_res_vs_par), 2), "]\n\n")
}

# LFC EXPLORATION -----------------------------------------------
lfc_primary <- if (DATASET %in% TWO_CLASS_DATASETS) lfc_res_vs_par else lfc_res_vs_trt

cat("  LFC distribution:\n")
q <- quantile(lfc_primary, probs = c(0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99))
for (nm in names(q)) cat("    ", nm, ":", round(q[nm], 2), "\n")
cat("    Genes |LFC| ≥", LOGFC_THRESHOLD, ":",
    sum(abs(lfc_primary) >= LOGFC_THRESHOLD), "\n\n")

p_hist <- ggplot(data.frame(LFC = lfc_primary), aes(LFC)) +
  geom_histogram(bins = 80, fill = "#74b9ff", color = "white", linewidth = 0.2) +
  geom_vline(xintercept = c(-LOGFC_THRESHOLD, LOGFC_THRESHOLD),
             linetype = "dashed", color = "#e17055") +
  labs(title = paste("LFC distribution — Resistance vs Parental —", OUT_PREFIX),
       x = "log2 Fold Change", y = "Genes") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"),
        panel.background = element_rect(fill = "white"),
        plot.background  = element_rect(fill = "white"))
print(p_hist)
ggsave(file.path(RESULTS_DIR, "figures", "02_lfc_histogram.png"),
       p_hist, width = 8, height = 5, dpi = 300)

# CLASSIFICATION ------------------------------------------------------------------
cat("3. Classifying genes...\n")

if (DATASET %in% TWO_CLASS_DATASETS) {
  # 2 classes: no treated sample — Upregulated/Downregulated in Resistance vs Parental
  classification <- rep("Other", nrow(mat))
  classification[lfc_res_vs_par >  LOGFC_THRESHOLD] <- "Upregulated"
  classification[lfc_res_vs_par < -LOGFC_THRESHOLD] <- "Downregulated"
} else {
  # 5 classes: with treated sample
  classify_gene_pattern <- function(lfc_res_vs_trt, lfc_res_vs_unt, lfc_trt_vs_unt, threshold, labels) {
    class <- rep("Other", length(lfc_res_vs_trt))
    class[lfc_res_vs_unt < -threshold] <- labels["down"]
    class[lfc_res_vs_unt > threshold &
          abs(lfc_res_vs_trt) < threshold &
          lfc_trt_vs_unt > 0.5] <- labels["pre"]
    class[lfc_trt_vs_unt > threshold &
          abs(lfc_res_vs_trt) < threshold] <- labels["maint"]
    class[lfc_res_vs_trt > threshold] <- labels["acq"]
    class
  }
  classification <- classify_gene_pattern(
    lfc_res_vs_trt, lfc_res_vs_unt, lfc_trt_vs_unt, LOGFC_THRESHOLD, CLASS_LABELS
  )
}

tab <- table(classification)
cat("  Distribution:\n")
for (cls in names(tab)) {
  pct <- round(100 * tab[cls] / sum(tab), 1)
  cat("    ", cls, ":", tab[cls], "(", pct, "%)\n")
}
cat("\n")

# MASTER TABLE --------------------------------------------------------------
if (DATASET %in% TWO_CLASS_DATASETS) {
  lfc_table <- tibble(
    Gene           = rownames(mat),
    LFC_Res_vs_Par = lfc_res_vs_par,
    Mean_Expr_Par  = expr_par,
    Mean_Expr_Res  = expr_res,
    Classification = classification
  )
} else {
  lfc_table <- tibble(
    Gene            = rownames(mat),
    LFC_Res_vs_Trt  = lfc_res_vs_trt,
    LFC_Res_vs_Unt  = lfc_res_vs_unt,
    LFC_Trt_vs_Unt  = lfc_trt_vs_unt,
    Mean_Expr_Unt   = expr_unt,
    Mean_Expr_Trt   = expr_trt,
    Mean_Expr_Res   = expr_res,
    Classification  = classification
  )
}

write.csv(lfc_table, file.path(RESULTS_DIR, "tables", "02_lfc_master_table.csv"), row.names = FALSE)

classified_genes <- lfc_table %>% filter(Classification != "Other")
write.csv(classified_genes, file.path(RESULTS_DIR, "tables", "02_classified_genes.csv"), row.names = FALSE)
cat("  Tables saved: 02_lfc_master_table.csv | 02_classified_genes.csv\n\n")

# FIGURE 1: CLASSIFICATION BARPLOT --------------------------------------------
cat("4. Figures...\n")

class_levels <- if (DATASET %in% TWO_CLASS_DATASETS) {
  c("Upregulated", "Downregulated", "Other")
} else {
  c(unname(CLASS_LABELS["acq"]), unname(CLASS_LABELS["pre"]), unname(CLASS_LABELS["maint"]),
    "Downregulated", "Other")
}

class_summary <- lfc_table %>%
  count(Classification) %>%
  mutate(
    pct   = round(100 * n / sum(n), 1),
    label = paste0(n, "\n(", pct, "%)"),
    Classification = factor(Classification, levels = class_levels)
  )

p_bar <- ggplot(class_summary, aes(x = Classification, y = n, fill = Classification)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = label), vjust = -0.3, size = 3.5) +
  scale_fill_manual(values = classification_colors) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(
    title = paste("Gene classification —", OUT_PREFIX),
    subtitle = paste("|LFC| threshold =", LOGFC_THRESHOLD),
    x = NULL, y = "Number of genes"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none",
        panel.grid.major.x = element_blank(),
        plot.title = element_text(face = "bold"),
        panel.background = element_rect(fill = "white"),
        plot.background  = element_rect(fill = "white"))

print(p_bar)  # viewer only — count summary

# FIGURE 2: MA PLOT colored by class --------------------------------------
if (DATASET %in% TWO_CLASS_DATASETS) {
  # Genes to label: (1) genes from key lists, (2) top LFC per class
  key_genes <- if (DATASET == "MET1_METi") {
    unique(c(met_resistance_core, gimenez_common, gimenez_egfr_bypass, gimenez_group2))
  } else {
    egfr_resistance_core
  }
  top_up   <- lfc_table %>% filter(Classification == "Upregulated") %>%
    arrange(desc(LFC_Res_vs_Par)) %>% slice_head(n = 15) %>% pull(Gene)
  top_down <- lfc_table %>% filter(Classification == "Downregulated") %>%
    arrange(LFC_Res_vs_Par) %>% slice_head(n = 15) %>% pull(Gene)
  genes_to_label <- unique(c(
    intersect(key_genes, lfc_table$Gene[lfc_table$Classification != "Other"]),
    top_up, top_down
  ))

  plot_df <- lfc_table %>%
    mutate(
      Classification = factor(Classification, levels = c("Upregulated", "Downregulated", "Other")),
      label = ifelse(Gene %in% genes_to_label & Classification != "Other", Gene, "")
    ) %>%
    arrange(desc(Classification))

  p_scatter <- ggplot(plot_df,
                      aes(Mean_Expr_Par, LFC_Res_vs_Par, color = Classification)) +
    geom_point(alpha = 0.4, size = 0.7) +
    geom_hline(yintercept = c(-LOGFC_THRESHOLD, LOGFC_THRESHOLD),
               linetype = "dashed", color = "gray50", linewidth = 0.5) +
    geom_hline(yintercept = 0, color = "black", linewidth = 0.3) +
    geom_text_repel(aes(label = label), size = 2.5, max.overlaps = 20,
                    segment.size = 0.2, box.padding = 0.3, seed = 1) +
    scale_color_manual(values = classification_colors) +
    labs(
      title = paste("MA plot — Resistance vs Parental —", OUT_PREFIX),
      x = "Mean Parental log2CPM",
      y = "LFC Resistance vs Parental"
    ) +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold"),
          panel.background = element_rect(fill = "white"),
          plot.background  = element_rect(fill = "white"))
} else {
  plot_df <- lfc_table %>%
    mutate(
      Classification = factor(Classification, levels = class_levels),
      label = ifelse(Classification %in% c(unname(CLASS_LABELS["acq"]), "Downregulated") &
                       (abs(LFC_Res_vs_Trt) > LOGFC_THRESHOLD * 1.5 |
                        abs(LFC_Res_vs_Unt) > LOGFC_THRESHOLD * 1.5), Gene, "")
    ) %>%
    arrange(Classification)

  p_scatter <- ggplot(plot_df, aes(LFC_Res_vs_Unt, LFC_Res_vs_Trt, color = Classification)) +
    geom_point(alpha = 0.4, size = 0.8) +
    geom_hline(yintercept = c(-LOGFC_THRESHOLD, LOGFC_THRESHOLD),
               linetype = "dashed", color = "gray50", linewidth = 0.5) +
    geom_vline(xintercept = c(-LOGFC_THRESHOLD, LOGFC_THRESHOLD),
               linetype = "dashed", color = "gray50", linewidth = 0.5) +
    geom_text_repel(aes(label = label), size = 2.5, max.overlaps = 20,
                    segment.size = 0.2, box.padding = 0.3, seed = 1) +
    scale_color_manual(values = classification_colors) +
    labs(
      title = paste("LFC scatter —", OUT_PREFIX),
      x = "LFC Resistant vs Untreated",
      y = "LFC Resistant vs Treated (primary)"
    ) +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold"),
          panel.background = element_rect(fill = "white"),
          plot.background  = element_rect(fill = "white"))
}

ggsave(file.path(RESULTS_DIR, "figures", "02_lfc_scatter.png"),
       p_scatter, width = 9, height = 7, dpi = 300)
cat("  → 02_lfc_scatter.png\n")

# FIGURE 3: CLASSIFIED GENES HEATMAP ----------------------------------------
# Subsample: top 100 Up + top 100 Down by LFC so that the heatmap stays readable.
classified_idx <- which(classification != "Other")
cat("  Total classified genes:", length(classified_idx), "\n")

if (length(classified_idx) >= 2) {
  cal_zscore <- function(x) if (sd(x) == 0) rep(0, length(x)) else (x - mean(x)) / sd(x)

  up_label  <- if (DATASET %in% TWO_CLASS_DATASETS) "Upregulated" else unname(CLASS_LABELS["acq"])
  lfc_sort  <- if (DATASET %in% TWO_CLASS_DATASETS) lfc_res_vs_par else lfc_res_vs_trt
  up_idx    <- which(classification == up_label)
  down_idx  <- which(classification == "Downregulated")
  top_up    <- up_idx[order(lfc_sort[up_idx],   decreasing = TRUE)][seq_len(min(100, length(up_idx)))]
  top_down  <- down_idx[order(lfc_sort[down_idx])][seq_len(min(100, length(down_idx)))]
  hm_idx   <- c(top_up, top_down)

  cat("  Heatmap:", length(hm_idx), "genes (top", length(top_up), "Up + top", length(top_down), "Down)\n")

  z_mat <- t(apply(mat[hm_idx, ], 1, cal_zscore))

  if (DATASET == "MET1_METi") {
    annotation <- metadata_unit[, c("Condition", "GX_cluster"), drop = FALSE]
  } else {
    annotation <- metadata_unit[, !colnames(metadata_unit) %in% c("Sample", "Model", "Cell_Line"), drop = FALSE]
  }
  ann_colors  <- filter_annotation_colors(annotation)

  ph_lfc <- pheatmap(z_mat,
           color             = heatmap_zscore,
           show_rownames     = FALSE,
           show_colnames     = TRUE,
           annotation_col    = annotation,
           annotation_colors = ann_colors,
           cluster_rows      = TRUE,
           cluster_cols      = TRUE,
           clustering_method = "ward.D2",
           main = paste0("Classified genes (top 100 Up + top 100 Down) — ", OUT_PREFIX),
           silent = TRUE)

  grid::grid.newpage()
  grid::grid.draw(ph_lfc$gtable)

  save_pheatmap(ph_lfc,
                file.path(RESULTS_DIR, "figures", "02_heatmap_classified.png"),
                width = 10, height = 8)
  cat("  → 02_heatmap_classified.png\n\n")
} else {
  cat("  Fewer than 2 classified genes — heatmap skipped\n\n")
}

# Tables 02_lfc_master_table.csv and 02_classified_genes.csv are already saved above.
# rna_03 and rna_04 load them directly — no .RData is needed.
cat("Tables ready for rna_03 and rna_04.\n\n")

# REPORT --------------------------------------------------------------------
n_up   <- sum(classification == "Upregulated")
n_acq  <- sum(classification == unname(CLASS_LABELS["acq"]))
n_down <- sum(classification == "Downregulated")

if (DATASET %in% TWO_CLASS_DATASETS) {
  top10_acq <- lfc_table %>%
    filter(Classification == "Upregulated") %>%
    arrange(desc(LFC_Res_vs_Par)) %>%
    slice_head(n = 10) %>%
    dplyr::select(Gene, LFC_Res_vs_Par)

  top10_down <- lfc_table %>%
    filter(Classification == "Downregulated") %>%
    arrange(LFC_Res_vs_Par) %>%
    slice_head(n = 10) %>%
    dplyr::select(Gene, LFC_Res_vs_Par)

  report_text <- paste0(
    "# LFC Report — ", OUT_PREFIX, "\n\n",
    "**Threshold:** |LFC| ≥ ", LOGFC_THRESHOLD, " (simplified classification: no treated sample)\n\n",
    "## Classification\n",
    "| Class | N | % |\n",
    "|-------|---|---|\n",
    "| Upregulated | ", n_up, " | ", round(100*n_up/nrow(lfc_table),1), "% |\n",
    "| Downregulated | ", n_down, " | ", round(100*n_down/nrow(lfc_table),1), "% |\n",
    "| Other | ", nrow(lfc_table)-n_up-n_down, " | ",
    round(100*(nrow(lfc_table)-n_up-n_down)/nrow(lfc_table),1), "% |\n\n",
    "## Top 10 Upregulated (by LFC_Res_vs_Par)\n",
    paste(capture.output(print(as.data.frame(top10_acq), row.names = FALSE)), collapse = "\n"), "\n\n",
    "## Top 10 Downregulated (by LFC_Res_vs_Par)\n",
    paste(capture.output(print(as.data.frame(top10_down), row.names = FALSE)), collapse = "\n"), "\n\n",
    "## Figures\n",
    "- 02_lfc_histogram.png\n",
    "- 02_lfc_scatter.png (colored MA plot, key genes labelled)\n",
    "- 02_heatmap_classified.png (z-score per sample, top 100 Up + top 100 Down of the classified genes)\n"
  )
} else {
  acq_lbl  <- unname(CLASS_LABELS["acq"])
  pre_lbl  <- unname(CLASS_LABELS["pre"])
  maint_lbl <- unname(CLASS_LABELS["maint"])

  n_pre  <- sum(classification == pre_lbl)
  n_main <- sum(classification == maint_lbl)
  ratio  <- if (n_pre > 0) round(n_acq / n_pre, 2) else "Inf"

  top10_acq <- lfc_table %>%
    filter(Classification == acq_lbl) %>%
    arrange(desc(LFC_Res_vs_Trt)) %>%
    slice_head(n = 10) %>%
    dplyr::select(Gene, LFC_Res_vs_Trt)

  top10_down <- lfc_table %>%
    filter(Classification == "Downregulated") %>%
    arrange(LFC_Res_vs_Unt) %>%
    slice_head(n = 10) %>%
    dplyr::select(Gene, LFC_Res_vs_Unt)

  report_text <- paste0(
    "# LFC Report — ", OUT_PREFIX, "\n\n",
    "**Threshold:** |LFC| ≥ ", LOGFC_THRESHOLD, "\n\n",
    "## Classification\n",
    "| Class | N | % |\n",
    "|-------|---|---|\n",
    "| ", acq_lbl, " | ", n_acq, " | ", round(100*n_acq/nrow(lfc_table),1), "% |\n",
    "| ", pre_lbl, " | ", n_pre, " | ", round(100*n_pre/nrow(lfc_table),1), "% |\n",
    "| Downregulated | ", n_down, " | ", round(100*n_down/nrow(lfc_table),1), "% |\n",
    "| ", maint_lbl, " | ", n_main, " | ", round(100*n_main/nrow(lfc_table),1), "% |\n\n",
    "**Ratio ", acq_lbl, "/", pre_lbl, ":** ", ratio,
    " → [INTERPRET: does de novo resistance or pre-existing heterogeneity predominate?]\n\n",
    "## Top 10 ", acq_lbl, " (by LFC_Res_vs_Trt)\n",
    paste(capture.output(print(as.data.frame(top10_acq), row.names = FALSE)), collapse = "\n"), "\n\n",
    "## Top 10 Downregulated (by LFC_Res_vs_Unt)\n",
    paste(capture.output(print(as.data.frame(top10_down), row.names = FALSE)), collapse = "\n"), "\n\n",
    "## Figures\n",
    "- 02_lfc_scatter.png (LFC scatter)\n",
    "- 02_heatmap_classified.png (z-score per sample, ordered by LFC)\n"
  )
}
writeLines(report_text, file.path(RESULTS_DIR, "report", "02_lfc_report.md"))
cat("Report: report/02_lfc_report.md\n\n")

HAVE_REPORT <- requireNamespace("rmarkdown", quietly = TRUE) && rmarkdown::pandoc_available()
if (!HAVE_REPORT) cat("rmarkdown or pandoc not available: skipping the HTML report
")
if (HAVE_REPORT) {
  cat("Generating HTML report...\n")
  rmarkdown::render(here::here("scripts", "simple", "rna_02_lfc_report.Rmd"),
                    output_dir  = file.path(RESULTS_DIR, "report"),
                    output_file = "02_lfc_report.html",
                    quiet       = TRUE)
  cat("  → report/02_lfc_report.html\n\n")
}

cat("================================================================\n")
cat("  LFC completed —", OUT_PREFIX, "\n")
cat("================================================================\n\n")

