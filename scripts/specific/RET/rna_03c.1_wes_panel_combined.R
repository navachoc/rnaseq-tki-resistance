# ============================================================
# rna_03c.1_wes_panel_combined.R
# Curated WES (CNV) + scDNA (point SNV) + recurrent SNV panel: RET
# Everything related to this panel lives HERE (embedded data + heatmaps +
# stripplots); rna_03b.1_key_genes.R does not touch it, it only handles the
# d) signature.
#
# Generates, for each model (RET1, RET2):
#   - Panel heatmap (filtered by Model %in% c(model, "Ambos"))
#   - CNV genes stripplot (CORE_BYPASS) and SNV genes stripplot
#     (scDNA_panel/EXTRA_SNV)
# and a final combined heatmap (both models, FULL unfiltered panel, all
# samples together). Everything is saved in results/RET/combined/figures/.
# RNAseq_bulk pipeline
#
# Inputs: results/RET/{RET1,RET2}/tables/01_mat_tmm.csv and
#   01_metadata_unit.csv; config/colors.R. Embedded gene panel.
# Outputs (results/RET/combined/figures/): 03c_heatmap_wes_panel_<model>.png,
#   03c_stripplot_cnv_panel_<model>.png, 03c_stripplot_snv_panel_<model>.png,
#   03c_heatmap_wes_panel_combined.png
# Dependencies (run first): scripts/simple/rna_01.1_qc.R with
#   DATASET = "RET1" and DATASET = "RET2".
# Usage: paths are resolved with here::here() from the repository root
# (.here file), so it works from any working directory:
#   Rscript scripts/specific/RET/rna_03c.1_wes_panel_combined.R
# ============================================================
# CONFIG
# ============================================================
RECEPTOR   <- "RET"
CELL_LINES <- c("RET1", "RET2")
CONFIG_DIR <- here::here("config")
# ============================================================

DRUG      <- RECEPTOR  # placeholder: colors.R computes OUT_PREFIX from DRUG
CELL_LINE <- paste(CELL_LINES, collapse = "_")
source(file.path(CONFIG_DIR, "colors.R"))
OUT_PREFIX <- paste0(CELL_LINE, "_", RECEPTOR)  # overrides the one from colors.R

library(ggplot2)
library(pheatmap)
library(dplyr)

# CURATED WES/scDNA/SNV PANEL: embedded object, no external CSV --------------
# Model: RET1 | RET2 | Ambos.
ret_wes_gene_panel <- data.frame(
  Gene = c("EGFR", "MAP2K1", "MAP2K1", "PTEN", "PTEN", "RET", "PIK3CA",
           "STK11", "ALK", "TP53", "KRAS", "KRAS", "CDKN2A", "NRAS",
           "MET", "BRAF",
           "ITSN1", "FAM221A", "ARMCX4", "PGM2L1", "EMP2", "CROCC",
           "TPSD1", "KRAS", "NRAS", "AKAP1", "DDAH1", "PLEKHG1", "PRDX6",
           "HERC2", "TSC22D1",
           "TRPC4", "COL7A1", "DUS1L", "SEMA6A", "VWDE", "ABCA6",
           "IBA57", "ADAM12"),
  Model = c("RET1", "RET1", "RET2", "RET1", "RET2", "RET2", "RET2",
            "RET2", "RET2", "RET1", "RET1", "RET2", "RET2", "RET2",
            "Both", "Both",
            "RET1", "RET1", "RET1", "RET1", "RET1", "RET2",
            "RET2", "RET2", "RET2", "RET2", "RET2", "RET2", "RET2",
            "RET2", "RET2",
            "RET1", "RET1", "RET1", "RET1", "RET1", "RET2",
            "RET2", "RET2"),
  stringsAsFactors = FALSE
)

# CNV (WES, amp/del) vs SNV (point scDNA + recurrent) split; independent of
# ret_wes_gene_panel, used only to separate the 2 stripplots.
cnv_panel_genes <- c("EGFR", "MAP2K1", "PTEN", "RET", "PIK3CA", "STK11",
                     "ALK", "TP53", "KRAS", "CDKN2A", "NRAS", "MET", "BRAF")
snv_panel_genes <- c("ITSN1", "FAM221A", "ARMCX4", "PGM2L1", "EMP2", "CROCC",
                     "TPSD1", "KRAS", "NRAS", "AKAP1", "DDAH1", "PLEKHG1",
                     "PRDX6", "HERC2", "TSC22D1",
                     "TRPC4", "COL7A1", "DUS1L", "SEMA6A", "VWDE", "ABCA6",
                     "IBA57", "ADAM12")

save_pheatmap <- function(ph, filepath, width, height, res = 300) {
  png(filepath, width = width, height = height, units = "in", res = res)
  grid::grid.draw(ph$gtable)
  dev.off()
}
cal_zscore <- function(x) if (sd(x) == 0) rep(0, length(x)) else (x - mean(x)) / sd(x)

# Column annotation without Treatment (Drug already distinguishes No/BLU667/LOXO292)
make_annotation <- function(metadata_unit, extra_cols = character(0)) {
  cols <- c(extra_cols, "Condition", "Drug")
  ann <- metadata_unit[, cols, drop = FALSE]
  ann_colors <- filter_annotation_colors(ann)
  ann_colors <- setNames(lapply(names(ann_colors), function(nm) {
    vals <- ann_colors[[nm]]
    vals[names(vals) %in% unique(ann[[nm]])]
  }), names(ann_colors))

  # RET-specific colors (same as the RET WES heatmap and oncoprint scripts,
  # COLOR_DRUG/COLOR_MODEL).
  if (!is.null(ann_colors$Drug)) {
    ann_colors$Drug <- c(
      "BLU667"  = "#6C5CE7",
      "LOXO292" = "#00B894",
      "No"      = "#95A5A6"
    )[names(ann_colors$Drug)]
  }
  if (!is.null(ann_colors$Cell_Line)) {
    ann_colors$Cell_Line <- c(
      "RET1" = "#7570B3",
      "RET2"  = "#E7298A"
    )[names(ann_colors$Cell_Line)]
  }
  ann_colors
}

make_stripplot <- function(genes, mat, metadata_unit, sample_order, drugs,
                            title, subtitle, filepath) {
  genes <- unique(genes[genes %in% rownames(mat)])
  if (length(genes) < 1) return(invisible(NULL))

  df <- lapply(genes, function(g) {
    data.frame(
      Gene      = g,
      Sample    = factor(colnames(mat), levels = sample_order),
      Expr      = as.numeric(mat[g, ]),
      Condition = metadata_unit$Condition,
      Drug      = metadata_unit$Drug
    )
  }) %>% bind_rows() %>%
    mutate(Gene = factor(Gene, levels = genes))

  p <- ggplot(df, aes(x = Sample, y = Expr, color = Condition, shape = Drug)) +
    geom_point(size = 3, alpha = 0.9) +
    scale_color_manual(values = condition_colors) +
    scale_shape_manual(values = c(No = 15, setNames(c(16, 17, 18, 3, 4)[seq_along(drugs)], drugs))) +
    facet_wrap(~ Gene, scales = "free_y", ncol = 5) +
    labs(title = title, subtitle = subtitle, x = NULL, y = "log2CPM (TMM)") +
    theme_minimal(base_size = 10) +
    theme(
      legend.position   = "right",
      axis.text.x       = element_text(angle = 40, hjust = 1, size = 7),
      strip.text        = element_text(face = "bold.italic", size = 8),
      plot.title        = element_text(face = "bold"),
      plot.subtitle     = element_text(size = 7, color = "gray50"),
      panel.background  = element_rect(fill = "white"),
      plot.background   = element_rect(fill = "white")
    )

  n_cols <- min(5, length(genes))
  n_rows <- ceiling(length(genes) / n_cols)
  ggsave(filepath, p, width = n_cols * 2.5, height = n_rows * 3, dpi = 300)
  invisible(genes)
}

cat("\n================================================================\n")
cat("  PANEL WES —", OUT_PREFIX, "\n")
cat("================================================================\n\n")

# LOAD MATRICES AND METADATA FOR BOTH MODELS ---------------------------------
load_model <- function(cell_line) {
  mat_file  <- here::here("results", RECEPTOR, cell_line, "tables", "01_mat_tmm.csv")
  meta_file <- here::here("results", RECEPTOR, cell_line, "tables", "01_metadata_unit.csv")

  mat_df           <- read.csv(mat_file, check.names = FALSE)
  rownames(mat_df) <- mat_df$Gene
  mat              <- as.matrix(mat_df[, -1])

  meta_df           <- read.csv(meta_file, check.names = FALSE)
  rownames(meta_df) <- meta_df$Model
  metadata_unit     <- meta_df[, -1]

  list(mat = mat, metadata = metadata_unit)
}

models <- lapply(CELL_LINES, load_model)
names(models) <- CELL_LINES

# HEATMAP + STRIPPLOTS PER MODEL ---------------------------------------------
for (cl in CELL_LINES) {
  cat("1.", cl, "— heatmap + CNV/SNV stripplots...\n")

  mat           <- models[[cl]]$mat
  metadata_unit <- models[[cl]]$metadata
  RESULTS_DIR   <- here::here("results", RECEPTOR, "combined")
  dir.create(file.path(RESULTS_DIR, "figures"), recursive = TRUE, showWarnings = FALSE)

  untreated <- rownames(metadata_unit)[metadata_unit$Condition == "Parental" & metadata_unit$Treatment == "No"]
  drugs     <- setdiff(unique(metadata_unit$Drug), "No")
  treated_order   <- unlist(lapply(drugs, function(d)
    rownames(metadata_unit)[metadata_unit$Condition == "Parental" &
                            metadata_unit$Treatment == "Yes" & metadata_unit$Drug == d]))
  resistant_order <- unlist(lapply(drugs, function(d)
    rownames(metadata_unit)[metadata_unit$Condition == "Resistance" & metadata_unit$Drug == d]))
  sample_order <- c(untreated, treated_order, resistant_order)

  annotation_col <- metadata_unit[, c("Condition", "Drug"), drop = FALSE]
  ann_colors     <- make_annotation(metadata_unit)

  wes_panel <- ret_wes_gene_panel %>%
    filter(Model %in% c(cl, "Both")) %>%
    distinct(Gene)

  wes_genes     <- wes_panel$Gene[wes_panel$Gene %in% rownames(mat)]
  missing_genes <- setdiff(wes_panel$Gene, rownames(mat))
  if (length(missing_genes) > 0) {
    cat("   genes not found:", paste(missing_genes, collapse = ", "), "\n")
  }

  if (length(wes_genes) >= 3) {
    z_wes <- t(apply(mat[wes_genes, ], 1, cal_zscore))

    ph_wes <- pheatmap(z_wes,
             annotation_col    = annotation_col,
             annotation_colors = ann_colors,
             cluster_rows      = TRUE,
             cluster_cols      = TRUE,
             clustering_method = "ward.D2",
             show_rownames     = TRUE,
             show_colnames     = TRUE,
             fontsize_row      = 8,
             fontsize_col      = 9,
             main = paste0("Panel WES/scDNA/SNV — ", cl, "_", RECEPTOR, " (", length(wes_genes), " genes)"),
             silent = TRUE)

    save_pheatmap(ph_wes, file.path(RESULTS_DIR, "figures", paste0("03c_heatmap_wes_panel_", cl, ".png")),
                  width = 12, height = 14)
    cat("   → 03c_heatmap_wes_panel_", cl, ".png (", length(wes_genes), "genes)\n", sep = "")
  }

  cnv_genes <- wes_panel$Gene[wes_panel$Gene %in% cnv_panel_genes]
  snv_genes <- wes_panel$Gene[wes_panel$Gene %in% snv_panel_genes]

  make_stripplot(cnv_genes, mat, metadata_unit, sample_order, drugs,
                 paste0("Panel CNV (WES) — ", cl, "_", RECEPTOR),
                 "Genes with amplification/deletion (CORE_BYPASS)",
                 file.path(RESULTS_DIR, "figures", paste0("03c_stripplot_cnv_panel_", cl, ".png")))
  cat("   → 03c_stripplot_cnv_panel_", cl, ".png (", length(cnv_genes), "genes)\n", sep = "")

  make_stripplot(snv_genes, mat, metadata_unit, sample_order, drugs,
                 paste0("SNV panel (scDNA + recurrent) — ", cl, "_", RECEPTOR),
                 "Genes with missense/stop-gained/recurrent SNV",
                 file.path(RESULTS_DIR, "figures", paste0("03c_stripplot_snv_panel_", cl, ".png")))
  cat("   → 03c_stripplot_snv_panel_", cl, ".png (", length(snv_genes), "genes)\n\n", sep = "")
}

# COMBINED HEATMAP (both models, full panel, not filtered by Model) -----
cat("2. Combined heatmap (", paste(CELL_LINES, collapse = " + "), ", full panel)...\n")

common_genes <- Reduce(intersect, lapply(models, function(m) rownames(m$mat)))
mat_combined <- do.call(cbind, lapply(models, function(m) m$mat[common_genes, , drop = FALSE]))

meta_combined <- do.call(rbind, lapply(names(models), function(cl) {
  models[[cl]]$metadata[, c("Cell_Line", "Condition", "Drug")]
}))
meta_combined <- meta_combined[colnames(mat_combined), , drop = FALSE]

ann_colors_combined <- make_annotation(meta_combined, extra_cols = "Cell_Line")

wes_panel_full <- ret_wes_gene_panel %>%
  distinct(Gene)

wes_genes_full     <- wes_panel_full$Gene[wes_panel_full$Gene %in% common_genes]
missing_genes_full <- setdiff(wes_panel_full$Gene, common_genes)
if (length(missing_genes_full) > 0) {
  cat("   genes not found in the combined matrix:",
      paste(missing_genes_full, collapse = ", "), "\n")
}

if (length(wes_genes_full) >= 3) {
  z_wes_full <- t(apply(mat_combined[wes_genes_full, ], 1, cal_zscore))

  ph_combined <- pheatmap(z_wes_full,
           annotation_col    = meta_combined,
           annotation_colors = ann_colors_combined,
           cluster_rows      = TRUE,
           cluster_cols      = TRUE,
           clustering_method = "ward.D2",
           show_rownames     = TRUE,
           show_colnames     = TRUE,
           fontsize_row      = 8,
           fontsize_col      = 8,
           main = paste0("Combined WES/scDNA/SNV panel — ", paste(CELL_LINES, collapse = " + "),
                         " (", length(wes_genes_full), " genes, ", ncol(mat_combined), " samples)"),
           silent = TRUE)

  dir.create(here::here("results", RECEPTOR, "combined", "figures"), recursive = TRUE, showWarnings = FALSE)
  save_pheatmap(ph_combined,
                here::here("results", RECEPTOR, "combined", "figures", "03c_heatmap_wes_panel_combined.png"),
                width = 16, height = 14)
  cat("   → 03c_heatmap_wes_panel_combined.png (", length(wes_genes_full), "genes,",
      ncol(mat_combined), "samples)\n\n")
}

cat("================================================================\n")
cat("  PANEL WES finished —", OUT_PREFIX, "\n")
cat("================================================================\n\n")
