# =============================================================================
# met_egfr_offon_regenerate.R
#
# WARNING: this script writes files with the SAME NAMES as the simple pipeline
# (scripts/simple/rna_01.1_qc.R ... rna_03.1_key_genes.R / rna_04.1_pathways.R),
# e.g. figures/01_pca_top500.png, figures/02_lfc_scatter.png,
# figures/03_heatmap_signature.png and tables/figures 04_go_bp_* / 04_kegg_*,
# and it will OVERWRITE them. Run it only when you want the extended-signature
# versions to replace the standard pipeline outputs.
#
# Purpose: extended conservative signature with "off in parental -> on in
# resistant" genes for MET (MET1) and EGFR (EGFR1), and regeneration of the
# enrichment figures and tables with it.
#
# Standalone -- it does not modify rna_03.1_key_genes.R or
# rna_04.1_pathways.R; it directly reads the tables those scripts already
# generated (02_lfc_master_table, resistance_signature.csv, 01_mat_tmm.csv,
# 01_metadata_unit.csv) and the raw counts, as is done for RET in the
# dedicated RET figure-regeneration scripts.
#
# Rationale: the conservative signature (Mean_Expr_Par > -5 | Forced)
# excludes genes with 0 reads in the parental even when they are clearly and
# consistently activated in the resistant models. This script identifies
# them, adds them to an extended conservative signature (OffOn column) and
# regenerates with it the signature heatmap (black->grey->gold gradient) and
# the GO BP / KEGG ORA (default clusterProfiler colour scale). It also
# regenerates the PCA (sample labels in black) and the MA plot.
#
# "off->on" criterion: Classification == "Upregulated" & Mean_Expr_Par <= -5 &
# raw reads in the parental == 0 & Mean_Expr_Res > 2 (log2CPM, ~4 CPM,
# substantial expression). It excludes genes sitting at the noise floor in
# both parental and resistant (high LFC but a numerical artefact of two
# near-zero values).
#
# Inputs (per dataset, results/MET/MET1/ and results/EGFR/EGFR1/):
#   tables/<MET|EGFR>_resistance_signature.csv
#   tables/<MET|EGFR>_resistance_signature_conservative.csv
#   tables/01_mat_tmm.csv, tables/01_metadata_unit.csv
#   tables/02_lfc_master_table.csv
#   data/counts_MET1_METi.csv, data/counts_EGFR1_EGFRi.csv
#   config/colors.R, config/gene_lists.R
#
# Outputs (per dataset):
#   tables/<MET|EGFR>_resistance_signature_conservative_offon.csv (new, does
#     not overwrite the original signature)
#   figures/03_heatmap_signature.png
#   figures/04_go_bp_{upregulated,downregulated}.png + tables/*.csv
#   figures/04_kegg_{upregulated,downregulated}.png + tables/*.csv
#   figures/01_pca_top500.png, figures/02_lfc_scatter.png
#   (they overwrite the same-named pipeline files; the MET1 04_* tables are
#   read by scripts/specific/MET/MET_01_summary_figures.R)
#
# Dependencies (run first, scripts/simple/): rna_01.1_qc.R,
# rna_02.1_lfc.R and rna_03.1_key_genes.R with DATASET = "MET1_METi" and
# DATASET = "EGFR1_EGFRi".
#
# Usage: paths are resolved with here::here() from the repository root
# (.here file), so it works from any working directory:
#   Rscript scripts/utils/met_egfr_offon_regenerate.R
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(clusterProfiler)
  library(org.Hs.eg.db)
})

CONFIG_DIR  <- here::here("config")

# colors.R expects CELL_LINE/DRUG to be defined (it uses them for its own
# OUT_PREFIX, which is not used here; each dataset defines its own
# cfg$out_prefix)
CELL_LINE <- "placeholder"
DRUG      <- "placeholder"
source(file.path(CONFIG_DIR, "colors.R"))
source(file.path(CONFIG_DIR, "gene_lists.R"))

OFFON_RES_THRESHOLD <- 2   # log2CPM en resistentes para considerar "on" real
TOP_N_PATHWAYS      <- 20
MIN_GENES_SET       <- 10

# New gradient (low -> high: black -> grey -> gold), same as the RET
# gradient figures. Only used in the signature heatmap.
GRADIENTE_RNA    <- c("#000000", "#808080", "#f2a900")
col_fun_heatmap  <- colorRampPalette(GRADIENTE_RNA)(100)

detect_sep <- function(file) if (grepl(";", readLines(file, n = 1))) ";" else ","
cal_zscore <- function(x) if (sd(x) == 0) rep(0, length(x)) else (x - mean(x)) / sd(x)
is_clean_gene <- function(g) !grepl("\\.", g) & !grepl("^[0-9]", g)

save_pheatmap <- function(ph, filepath, width, height, res = 300) {
  png(filepath, width = width, height = height, units = "in", res = res)
  grid::grid.draw(ph$gtable)
  dev.off()
}

# =============================================================================
# PER-DATASET CONFIGURATION
# =============================================================================
DATASETS <- list(
  MET = list(
    sig_prefix   = "MET",
    results_dir  = here::here("results", "MET", "MET1"),
    counts_file  = here::here("data", "counts_MET1_METi.csv"),
    parental_col = "MET1",
    out_prefix   = "MET1_METi",
    gx_cluster   = TRUE
  ),
  EGFR = list(
    sig_prefix   = "EGFR",
    results_dir  = here::here("results", "EGFR", "EGFR1"),
    counts_file  = here::here("data", "counts_EGFR1_EGFRi.csv"),
    parental_col = "EGFR1",
    out_prefix   = "EGFR1_EGFR",
    gx_cluster   = FALSE
  )
)

# =============================================================================
# FUNCTION: identify off->on genes and build the extended conservative signature
# =============================================================================
build_offon_signature <- function(cfg) {
  cat("\n=== ", cfg$sig_prefix, " — detecting off->on genes ===\n", sep = "")

  sig_full <- read.csv(file.path(cfg$results_dir, "tables",
                                  paste0(cfg$sig_prefix, "_resistance_signature.csv")))
  sig_cons <- read.csv(file.path(cfg$results_dir, "tables",
                                  paste0(cfg$sig_prefix, "_resistance_signature_conservative.csv")))

  counts <- read.csv(cfg$counts_file, sep = detect_sep(cfg$counts_file), check.names = FALSE)

  # Count-matrix columns are named by Sample; the parental is configured by Model
  meta_unit    <- read.csv(file.path(cfg$results_dir, "tables", "01_metadata_unit.csv"),
                           check.names = FALSE)
  parental_col <- meta_unit$Sample[meta_unit[[1]] == cfg$parental_col]

  candidatos <- sig_full %>%
    filter(Classification == "Upregulated", Mean_Expr_Par <= -5, Mean_Expr_Res > OFFON_RES_THRESHOLD) %>%
    inner_join(counts[, c("gene_name", parental_col)], by = c("Gene" = "gene_name")) %>%
    filter(.data[[parental_col]] == 0) %>%
    dplyr::select(-all_of(parental_col))

  cat("  off->on genes found:", nrow(candidatos), "\n")
  if (nrow(candidatos) > 0) cat("  ", paste(candidatos$Gene, collapse = ", "), "\n")

  ya_incluidos <- intersect(candidatos$Gene, sig_cons$Gene)
  if (length(ya_incluidos) > 0) {
    cat("  (", length(ya_incluidos), "already in the conservative signature, ignored)\n")
    candidatos <- candidatos %>% filter(!Gene %in% ya_incluidos)
  }

  if (!"OffOn" %in% colnames(sig_cons)) sig_cons$OffOn <- FALSE

  if (nrow(candidatos) > 0) {
    candidatos$OffOn <- TRUE
    sig_extendida <- bind_rows(sig_cons, candidatos)
    cat("  Extended signature:", nrow(sig_extendida), "genes (", nrow(sig_cons), "+", nrow(candidatos), "new )\n")
  } else {
    cat("  No new off->on genes, the conservative signature already includes them\n")
    sig_extendida <- sig_cons
  }

  out_csv <- file.path(cfg$results_dir, "tables",
                        paste0(cfg$sig_prefix, "_resistance_signature_conservative_offon.csv"))
  write.csv(sig_extendida, out_csv, row.names = FALSE)
  cat("  Saved:", out_csv, "\n")

  sig_extendida
}

# =============================================================================
# FUNCTION: signature heatmap (top 15 up + top 15 down), new gradient
# (the variables top25_up/top25_down keep their name but take 15 genes)
# =============================================================================
plot_signature_heatmap <- function(cfg, sig_extendida) {
  cat("\n=== ", cfg$sig_prefix, " — signature heatmap (black-grey-gold gradient) ===\n", sep = "")

  mat_df <- read.csv(file.path(cfg$results_dir, "tables", "01_mat_tmm.csv"), check.names = FALSE)
  rownames(mat_df) <- mat_df$Gene
  mat <- as.matrix(mat_df[, -1])

  meta_df <- read.csv(file.path(cfg$results_dir, "tables", "01_metadata_unit.csv"), check.names = FALSE)
  rownames(meta_df) <- meta_df$Model
  metadata_unit <- meta_df[, -1, drop = FALSE]

  if (cfg$gx_cluster && !"GX_cluster" %in% colnames(metadata_unit)) {
    rg_map <- c(
      MET1 = "Parental",
      MET1_C1 = "Group1_EGFR_bypass", MET1_PR1 = "Group1_EGFR_bypass",
      MET1_PR2 = "Group1_EGFR_bypass", MET1_V8  = "Group1_EGFR_bypass",
      MET1_PR3 = "Group2_NF2_YAP",     MET1_R = "Group2_NF2_YAP"
    )
    metadata_unit$GX_cluster <- rg_map[rownames(metadata_unit)]
  }
  ann_sig <- if (cfg$gx_cluster) metadata_unit[, c("Condition", "GX_cluster"), drop = FALSE] else
                                 metadata_unit[, "Condition", drop = FALSE]
  ann_colors <- filter_annotation_colors(ann_sig)

  sig_genes <- sig_extendida$Gene[sig_extendida$Gene %in% rownames(mat) & is_clean_gene(sig_extendida$Gene)]

  sig_up_all   <- sig_extendida %>% filter(Classification == "Upregulated") %>%
    arrange(desc(LFC_Res_vs_Par)) %>% pull(Gene)
  sig_down_all <- sig_extendida %>% filter(Classification == "Downregulated") %>%
    arrange(LFC_Res_vs_Par) %>% pull(Gene)

  top25_up   <- head(sig_up_all[sig_up_all %in% sig_genes], 15)
  top25_down <- head(sig_down_all[sig_down_all %in% sig_genes], 15)
  sig_ordered <- c(top25_up, top25_down)

  z_sig <- t(apply(mat[sig_ordered, , drop = FALSE], 1, cal_zscore))

  ph_sig <- pheatmap(z_sig,
           color             = col_fun_heatmap,
           annotation_col    = ann_sig,
           annotation_colors = ann_colors,
           cluster_rows      = TRUE,
           cluster_cols      = TRUE,
           clustering_method = "ward.D2",
           show_rownames     = TRUE,
           show_colnames     = TRUE,
           fontsize_row      = 11,
           fontsize_col      = 10,
           main = paste0(cfg$sig_prefix, " resistance signature (+ off->on) — top30 — ", cfg$out_prefix,
                         " (", length(top25_up), " up / ", length(top25_down), " down)"),
           silent = TRUE)

  out_png <- file.path(cfg$results_dir, "figures", "03_heatmap_signature.png")
  save_pheatmap(ph_sig, out_png, width = 12, height = 14)
  cat("  Saved:", out_png, "\n")
}

# =============================================================================
# FUNCTION: ORA GO BP + KEGG (Upregulated / Downregulated), colours unchanged
# =============================================================================
run_enrichment <- function(gene_symbols, label, results_dir, out_prefix, top_n) {
  cat("  ORA:", label, "(", length(gene_symbols), "genes)\n")
  if (length(gene_symbols) < 5) { cat("    Fewer than 5 genes, skipped\n"); return(invisible(NULL)) }

  entrez <- bitr(gene_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  if (nrow(entrez) == 0) { cat("    No Entrez conversion\n"); return(invisible(NULL)) }
  prefix <- tolower(gsub(" ", "_", label))

  tryCatch({
    ego <- enrichGO(gene = entrez$ENTREZID, OrgDb = org.Hs.eg.db, ont = "BP",
                    pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = TRUE)
    if (!is.null(ego) && nrow(ego@result) > 0) {
      p_go <- dotplot(ego, showCategory = min(top_n, nrow(ego@result))) +
        labs(title = paste("GO BP —", label, "—", out_prefix)) +
        theme(plot.title = element_text(face = "bold"),
              panel.background = element_rect(fill = "white"),
              plot.background  = element_rect(fill = "white"))
      ggsave(file.path(results_dir, "figures", paste0("04_go_bp_", prefix, ".png")),
             p_go, width = 10, height = 8, dpi = 300)
      write.csv(ego@result, file.path(results_dir, "tables", paste0("04_go_bp_", prefix, ".csv")),
                row.names = FALSE)
      cat("    ->  04_go_bp_", prefix, ".png\n", sep = "")
    } else cat("    GO BP: no results\n")
  }, error = function(e) cat("    GO BP error:", conditionMessage(e), "\n"))

  tryCatch({
    ekegg <- enrichKEGG(gene = entrez$ENTREZID, organism = "hsa",
                        pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2)
    if (!is.null(ekegg) && nrow(ekegg@result) > 0) {
      p_kegg <- dotplot(ekegg, showCategory = min(top_n, nrow(ekegg@result))) +
        labs(title = paste("KEGG —", label, "—", out_prefix)) +
        theme(plot.title = element_text(face = "bold"),
              panel.background = element_rect(fill = "white"),
              plot.background  = element_rect(fill = "white"))
      ggsave(file.path(results_dir, "figures", paste0("04_kegg_", prefix, ".png")),
             p_kegg, width = 10, height = 8, dpi = 300)
      write.csv(ekegg@result, file.path(results_dir, "tables", paste0("04_kegg_", prefix, ".csv")),
                row.names = FALSE)
      cat("    ->  04_kegg_", prefix, ".png\n", sep = "")
    } else cat("    KEGG: no results\n")
  }, error = function(e) cat("    KEGG error:", conditionMessage(e), "\n"))
}

regenerate_pathways <- function(cfg, sig_extendida) {
  cat("\n=== ", cfg$sig_prefix, " — GO BP / KEGG (extended signature, colours unchanged) ===\n", sep = "")
  genes_up   <- sig_extendida %>% filter(Classification == "Upregulated") %>% pull(Gene)
  genes_down <- sig_extendida %>% filter(Classification == "Downregulated") %>% pull(Gene)

  run_enrichment(genes_up, "Upregulated", cfg$results_dir, cfg$out_prefix, TOP_N_PATHWAYS)
  if (length(genes_down) >= MIN_GENES_SET) {
    run_enrichment(genes_down, "Downregulated", cfg$results_dir, cfg$out_prefix, TOP_N_PATHWAYS)
  } else {
    cat("  Downregulated: only", length(genes_down), "genes, skipped\n")
  }
}

# =============================================================================
# FUNCTION: PCA (01_pca_top500.png) — dpi 300, sample labels in black
# (in the original, rna_01.1_qc.R, the labels inherit the Condition colour,
# which makes them hard to read over a point of the same colour)
# =============================================================================
plot_pca <- function(cfg) {
  cat("\n=== ", cfg$sig_prefix, " — PCA (dpi 300, black labels) ===\n", sep = "")

  mat_df <- read.csv(file.path(cfg$results_dir, "tables", "01_mat_tmm.csv"), check.names = FALSE)
  rownames(mat_df) <- mat_df$Gene
  mat <- as.matrix(mat_df[, -1])

  meta_df <- read.csv(file.path(cfg$results_dir, "tables", "01_metadata_unit.csv"), check.names = FALSE)
  rownames(meta_df) <- meta_df$Model
  metadata_unit <- meta_df[, -1, drop = FALSE]

  TOP_N_PCA <- 500
  gene_var  <- apply(mat, 1, var)
  n_pca     <- min(TOP_N_PCA, nrow(mat))
  top_genes <- names(sort(gene_var, decreasing = TRUE))[seq_len(n_pca)]
  pca       <- prcomp(t(mat[top_genes, ]), scale. = TRUE)
  var_pc1   <- round(summary(pca)$importance[2, 1] * 100, 1)
  var_pc2   <- round(summary(pca)$importance[2, 2] * 100, 1)
  pca_df    <- data.frame(pca$x, metadata_unit)
  pca_df$Sample <- rownames(pca_df)

  p_pca <- ggplot(pca_df, aes(PC1, PC2, color = Condition, shape = Treatment, label = Sample)) +
    geom_point(size = 4, alpha = 0.85) +
    geom_text_repel(color = "black", size = 3, max.overlaps = 20,
                     box.padding = 0.3, segment.size = 0.2, seed = 1) +
    scale_color_manual(values = condition_colors) +
    scale_shape_manual(values = c(Yes = 17, No = 16)) +
    labs(
      title = paste0("PCA — ", cfg$out_prefix, " (top ", n_pca, " genes)"),
      x = paste0("PC1 (", var_pc1, "%)"),
      y = paste0("PC2 (", var_pc2, "%)")
    ) +
    theme_minimal(base_size = 12) +
    theme(plot.title       = element_text(face = "bold"),
          panel.grid       = element_line(color = "gray90", linetype = "dashed"),
          panel.background = element_rect(fill = "white"),
          plot.background  = element_rect(fill = "white"))

  out_png <- file.path(cfg$results_dir, "figures", paste0("01_pca_top", n_pca, ".png"))
  ggsave(out_png, p_pca, width = 10, height = 7, dpi = 300)
  cat("  Saved:", out_png, "\n")
}

# =============================================================================
# FUNCTION: MA plot (02_lfc_scatter.png) — dpi 300
# Gene labels (not model labels), kept in the Classification colour as in the
# original (rna_02.1_lfc.R).
# =============================================================================
plot_ma <- function(cfg) {
  cat("\n=== ", cfg$sig_prefix, " — MA plot (dpi 300) ===\n", sep = "")

  lfc_table <- read.csv(file.path(cfg$results_dir, "tables", "02_lfc_master_table.csv"), check.names = FALSE)

  key_genes <- if (cfg$sig_prefix == "MET") {
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

  p_scatter <- ggplot(plot_df, aes(Mean_Expr_Par, LFC_Res_vs_Par, color = Classification)) +
    geom_point(alpha = 0.4, size = 0.7) +
    geom_hline(yintercept = c(-2, 2), linetype = "dashed", color = "gray50", linewidth = 0.5) +
    geom_hline(yintercept = 0, color = "black", linewidth = 0.3) +
    geom_text_repel(aes(label = label), size = 2.5, max.overlaps = 20,
                     segment.size = 0.2, box.padding = 0.3, seed = 1) +
    scale_color_manual(values = classification_colors) +
    guides(color = guide_legend(override.aes = list(size = 5, alpha = 1))) +
    labs(
      title = paste("MA plot — Resistance vs Parental —", cfg$out_prefix),
      x = "Mean Parental log2CPM",
      y = "LFC Resistance vs Parental"
    ) +
    theme_minimal(base_size = 12) +
    theme(plot.title       = element_text(face = "bold"),
          panel.background = element_rect(fill = "white"),
          plot.background  = element_rect(fill = "white"))

  out_png <- file.path(cfg$results_dir, "figures", "02_lfc_scatter.png")
  ggsave(out_png, p_scatter, width = 9, height = 7, dpi = 300)
  cat("  Saved:", out_png, "\n")
}

# =============================================================================
# EXECUTION — MET and EGFR
# =============================================================================
for (nombre in names(DATASETS)) {
  cfg <- DATASETS[[nombre]]
  sig_extendida <- build_offon_signature(cfg)
  plot_signature_heatmap(cfg, sig_extendida)
  regenerate_pathways(cfg, sig_extendida)
  plot_pca(cfg)
  plot_ma(cfg)
}

cat("\n=== met_egfr_offon_regenerate.R finished ===\n")
