# ============================================================
# rna_01.1_qc.R
# Exploratory quality control of bulk RNA-seq: edgeR TMM normalization
# (log2CPM), expression distribution, PCA, sample distances and heatmap of
# the most variable genes. Generic for any dataset/receptor.
#
# Inputs:    data/counts_<DATASET>.csv, data/metadata_<DATASET>.csv
#            config/colors.R, config/gene_lists.R
# Outputs:   results/<RECEPTOR>/<CELL_LINE>/figures/01_*.png
#            results/<RECEPTOR>/<CELL_LINE>/tables/01_mat_tmm.csv,
#              01_metadata_unit.csv, 01_qc_stats.csv (input of rna_02)
#            results/<RECEPTOR>/<CELL_LINE>/report/01_qc_report.{md,html}
# Usage:     set the CONFIG block and, from any folder of the repo,
#            source(here::here("scripts/simple/rna_01.1_qc.R"))
# ============================================================
# CONFIG — change here for each run
# ============================================================
DATASET         <- "RET2"    # "ALK12" | "ALK3" | "MET1_METi" | "EGFR1_EGFRi" | "RET1" | "RET2"
CELL_LINE       <- "RET2"         # "ALK2" | "ALK1" | "ALK3" | "MET1" | "EGFR1" | "RET1" | "RET2"
DRUG            <- "RET_inhibitor" # "Lorlatinib" | "Alectinib" | "MET_inhibitor" | "Gefitinib_Afatinib" | "BLU667" | "LOXO292" | "RET_inhibitor"
RECEPTOR        <- "RET"          # "MET" | "ALK" | "EGFR" | "RET"
MIN_COUNT       <- 10
TOP_N_PCA       <- 500   # 50 | 100 | 500 | Inf
TOP_N_HEATMAP   <- 30   # 50 | 100 | 500
# ============================================================
DATA_DIR    <- here::here("data")
RESULTS_DIR <- here::here("results", RECEPTOR, CELL_LINE)
CONFIG_DIR  <- here::here("config")
source(file.path(CONFIG_DIR, "colors.R"))
source(file.path(CONFIG_DIR, "gene_lists.R"))
dir.create(file.path(RESULTS_DIR, "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(RESULTS_DIR, "tables"),  recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(RESULTS_DIR, "report"),  recursive = TRUE, showWarnings = FALSE)

# LIBRARIES ------------------------------------------------------------------
library(edgeR)
library(ggplot2)
library(pheatmap)
library(ggrepel)
library(RColorBrewer)
library(dplyr)

# FUNCTIONS ------------------------------------------------------------------
detect_sep <- function(file) {
  if (grepl(";", readLines(file, n = 1))) ";" else ","
}
cal_zscore <- function(x) if (sd(x) == 0) rep(0, length(x)) else (x - mean(x)) / sd(x)

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

cat("\n================================================================\n")
cat("  QC —", OUT_PREFIX, "\n")
cat("================================================================\n\n")

# STEP 1: LOAD DATA -------------------------------------------------------
cat("1. Loading data...\n")

counts_file <- file.path(DATA_DIR, paste0("counts_", DATASET, ".csv"))
meta_file   <- file.path(DATA_DIR, paste0("metadata_", DATASET, ".csv"))

counts   <- read.csv(counts_file, sep = detect_sep(counts_file), check.names = FALSE)
metadata <- read.csv(meta_file,   sep = detect_sep(meta_file))

rownames(counts) <- make.names(counts[, 2], unique = TRUE)
counts <- counts[, -c(1:3)]

sample_to_model <- setNames(metadata$Model, metadata$Sample)
colnames(counts) <- sample_to_model[colnames(counts)]
rownames(metadata) <- metadata$Model

if (DATASET %in% c("MET1_METi", "EGFR1_EGFRi")) {
  metadata_unit <- metadata
} else if (DATASET == "ALK12") {
  metadata_unit <- metadata %>% filter(Cell_Line == CELL_LINE)
} else if (DATASET %in% c("RET1", "RET2")) {
  # A single cell line per file — no split by drug (BLU667+LOXO292 together)
  metadata_unit <- metadata
} else {
  metadata_unit <- metadata %>% filter(Drug %in% c("No", DRUG))
}
counts_filt <- counts[, rownames(metadata_unit)]
genes_before <- nrow(counts_filt)
cat("  Samples:", ncol(counts_filt), "| Raw genes:", genes_before, "\n\n")

# STEP 2: TMM NORMALIZATION --------------------------------------------------
cat("2. Normalizing (edgeR TMM)...\n")

keep <- rowSums(counts_filt) >= MIN_COUNT
dge  <- DGEList(counts = counts_filt[keep, ])
dge  <- calcNormFactors(dge, method = "TMM")
mat  <- cpm(dge, log = TRUE, prior.count = 1)

retention_pct <- round(100 * nrow(mat) / genes_before, 1)
cat("  Genes after MIN_COUNT filter =", MIN_COUNT, ":", nrow(mat),
    "(", retention_pct, "% retained)\n\n")

# Reads per sample, medians and TMM factors
total_reads <- colSums(counts_filt)
med_expr    <- round(apply(mat, 2, median), 2)
tmm_factors <- setNames(round(dge$samples$norm.factors, 3), rownames(dge$samples))
cat("  Reads per sample:\n")
for (s in names(total_reads)) {
  cat("    ", s, ":", total_reads[s],
      "| median log2CPM:", med_expr[s],
      "| TMM factor:", tmm_factors[s], "\n")
}
cat("\n")

annotation <- metadata_unit[, !colnames(metadata_unit) %in% c("Sample", "Model", "Cell_Line"), drop = FALSE]
ann_colors <- filter_annotation_colors(annotation)

# STEP 3: DISTRIBUTION -------------------------------------------------------
cat("3. Expression distribution...\n")

df_dist <- data.frame(
  Expression = as.vector(mat),
  Sample     = rep(colnames(mat), each = nrow(mat)),
  Condition  = rep(metadata_unit$Condition, each = nrow(mat))
)

p_dist <- ggplot(df_dist, aes(Sample, Expression, fill = Condition)) +
  geom_violin(alpha = 0.7) +
  geom_boxplot(width = 0.1, alpha = 0.5, outlier.shape = NA) +
  scale_fill_manual(values = annotation_colors$Condition) +
  labs(title = paste("log2CPM distribution —", OUT_PREFIX),
       y = "log2CPM", x = NULL) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x   = element_text(angle = 45, hjust = 1),
        plot.title     = element_text(face = "bold"),
        panel.background = element_rect(fill = "white"),
        plot.background  = element_rect(fill = "white"))

print(p_dist)

# TMM factors (viewer only — normalization check)
tmm_df <- data.frame(
  Sample  = names(tmm_factors),
  Factor  = tmm_factors,
  Condition = metadata_unit$Condition
)
p_tmm <- ggplot(tmm_df, aes(Sample, Factor, fill = Condition)) +
  geom_col(width = 0.6) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray40") +
  scale_fill_manual(values = annotation_colors$Condition) +
  labs(title = paste("TMM normalization factors —", OUT_PREFIX),
       y = "Factor", x = NULL) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x    = element_text(angle = 45, hjust = 1),
        plot.title      = element_text(face = "bold"),
        panel.background = element_rect(fill = "white"),
        plot.background  = element_rect(fill = "white"))
print(p_tmm)

# STEP 4: PCA ----------------------------------------------------------------
cat("4. PCA (top", TOP_N_PCA, "variable genes)...\n")

gene_var  <- apply(mat, 1, var)
n_pca     <- min(TOP_N_PCA, nrow(mat))
top_genes <- names(sort(gene_var, decreasing = TRUE))[seq_len(n_pca)]
pca       <- prcomp(t(mat[top_genes, ]), scale. = TRUE)
var_pc1   <- round(summary(pca)$importance[2, 1] * 100, 1)
var_pc2   <- round(summary(pca)$importance[2, 2] * 100, 1)
pca_df    <- data.frame(pca$x, metadata_unit)
pca_df$Sample <- rownames(pca_df)

cat("  PC1:", var_pc1, "% | PC2:", var_pc2, "%\n")

p_pca <- ggplot(pca_df, aes(PC1, PC2, color = Condition, shape = Treatment, label = Sample)) +
  geom_point(size = 4, alpha = 0.85) +
  geom_text_repel(size = 3, max.overlaps = 20, box.padding = 0.3, segment.size = 0.2,
                  seed = 1) +
  scale_color_manual(values = condition_colors) +
  scale_shape_manual(values = c(Yes = 17, No = 16)) +
  labs(
    title = paste0("PCA — ", OUT_PREFIX, " (top ", n_pca, " genes)"),
    x = paste0("PC1 (", var_pc1, "%)"),
    y = paste0("PC2 (", var_pc2, "%)")
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title       = element_text(face = "bold"),
        panel.grid       = element_line(color = "gray90", linetype = "dashed"),
        panel.background = element_rect(fill = "white"),
        plot.background  = element_rect(fill = "white"))

print(p_pca)
ggsave(file.path(RESULTS_DIR, "figures", paste0("01_pca_top", n_pca, ".png")),
       p_pca, width = 10, height = 7, dpi = 300)
cat("  → 01_pca_top", n_pca, ".png\n\n", sep = "")

# STEP 5: SAMPLE DISTANCES ------------------------------------------
cat("5. Euclidean distance heatmap...\n")

dist_mat <- as.matrix(dist(t(mat)))

ph_dist <- pheatmap(dist_mat,
         clustering_distance_rows = as.dist(dist_mat),
         clustering_distance_cols = as.dist(dist_mat),
         color          = heatmap_distances,
         annotation_col = annotation,
         annotation_colors = ann_colors,
         main           = paste("Sample distances —", OUT_PREFIX),
         silent         = TRUE)

grid::grid.newpage()
grid::grid.draw(ph_dist$gtable)
save_pheatmap(ph_dist,
              file.path(RESULTS_DIR, "figures", "01_sample_distances.png"),
              width = 8, height = 7)
cat("  → 01_sample_distances.png\n\n")

# STEP 6: VARIABLE GENES HEATMAP --------------------------------------------
cat("6. Heatmap top", TOP_N_HEATMAP, "variable genes...\n")

top_idx <- order(gene_var, decreasing = TRUE)[seq_len(TOP_N_HEATMAP)]
fsize   <- ifelse(TOP_N_HEATMAP > 50, 6, 8)

# With only 2 samples (e.g. EGFR1_EGFRi: EGFR1/EGFR1_R) the per-gene z-score is
# degenerate -- each row ends up with a single sign, a solid block of 2
# colors with no real gradient (a "flag"), not a graded signature. With
# ncol(mat) <= 2, raw log2CPM and a magnitude gradient are used instead of
# the z-score, so that the color reflects actual expression.
usar_zscore <- ncol(mat) > 2

if (usar_zscore) {
  plot_mat  <- t(apply(mat[top_idx, ], 1, cal_zscore))
  color_top <- heatmap_zscore
} else {
  plot_mat  <- mat[top_idx, ]
  color_top <- colorRampPalette(c("#000000", "#808080", "#f2a900"))(100)
}

ph_top <- pheatmap(plot_mat,
         color             = color_top,
         show_rownames     = (TOP_N_HEATMAP <= 100),
         show_colnames     = TRUE,
         fontsize_row      = fsize,
         annotation_col    = annotation,
         annotation_colors = ann_colors,
         cluster_rows      = TRUE,
         cluster_cols      = TRUE,
         clustering_method = "ward.D2",
         main              = paste0("Top ", TOP_N_HEATMAP, " variable genes — ", OUT_PREFIX,
                                     if (!usar_zscore) " (raw log2CPM, n<=2 samples)" else ""),
         silent            = TRUE)

grid::grid.newpage()
grid::grid.draw(ph_top$gtable)
save_pheatmap(ph_top,
              file.path(RESULTS_DIR, "figures", paste0("01_heatmap_top", TOP_N_HEATMAP, ".png")),
              width = 14, height = 9)
cat("  → 01_heatmap_top", TOP_N_HEATMAP, ".png\n\n", sep = "")

# INTERMEDIATE TABLES (input for rna_02) -------------------------------------
mat_out <- as.data.frame(mat)
mat_out <- cbind(Gene = rownames(mat_out), mat_out)
write.csv(mat_out, file.path(RESULTS_DIR, "tables", "01_mat_tmm.csv"), row.names = FALSE)
write.csv(cbind(Model = rownames(metadata_unit), metadata_unit),
          file.path(RESULTS_DIR, "tables", "01_metadata_unit.csv"), row.names = FALSE)

qc_stats <- data.frame(
  Sample         = names(total_reads),
  Condition      = metadata_unit[names(total_reads), "Condition"],
  Reads          = as.integer(total_reads),
  Median_log2CPM = as.numeric(med_expr),
  TMM_factor     = as.numeric(tmm_factors),
  Genes_input    = genes_before,
  Genes_filtered = nrow(mat),
  Retention_pct  = retention_pct,
  PC1_var        = var_pc1,
  PC2_var        = var_pc2,
  stringsAsFactors = FALSE
)
write.csv(qc_stats, file.path(RESULTS_DIR, "tables", "01_qc_stats.csv"), row.names = FALSE)
cat("Tables saved: 01_mat_tmm.csv | 01_metadata_unit.csv | 01_qc_stats.csv\n\n")

# REPORT --------------------------------------------------------------------
core_list     <- if (DATASET == "MET1_METi") met_resistance_core else alk_resistance_core
core_in_top_heatmap <- intersect(rownames(mat)[top_idx], core_list)

sink(file.path(RESULTS_DIR, "report", "01_qc_report.md"))
tryCatch({
  cat("================================================================\n")
  cat("QC REPORT —", OUT_PREFIX, "\n")
cat("Date:", as.character(Sys.time()), "\n")
cat("================================================================\n\n")
cat("PROCESSING:\n")
cat("  • Initial genes:", genes_before, "\n")
cat("  • Genes after filter (MIN_COUNT =", MIN_COUNT, "):", nrow(mat), "\n")
cat("  • Retained:", retention_pct, "%\n")
cat("  • Samples:", ncol(mat), "\n\n")
cat("QUALITY PER SAMPLE:\n")
cat(sprintf("  %-20s %12s %18s %12s\n", "Sample", "Reads", "Median log2CPM", "TMM factor"))
for (s in names(total_reads)) {
  cat(sprintf("  %-20s %12d %18.2f %12.3f\n",
              s, total_reads[s], med_expr[s], tmm_factors[s]))
}
cat("\n")
cat("PCA (top", TOP_N_PCA, "genes):\n")
cat("  • PC1:", var_pc1, "% | PC2:", var_pc2, "%\n\n")
cat("CORE GENES IN TOP", TOP_N_HEATMAP, ":\n")
cat(" ", if (length(core_in_top_heatmap) > 0) paste(core_in_top_heatmap, collapse = ", ") else "None", "\n\n")
cat("FOR INTERPRETATION:\n")
cat("  1. 01_pca_top", TOP_N_PCA, ".png      — Separation of parental vs resistance?\n", sep = "")
cat("  2. 01_sample_distances.png — Coherent clustering?\n")
cat("  3. 01_heatmap_top", TOP_N_HEATMAP, ".png  — Clear expression patterns?\n", sep = "")
}, finally = sink())

cat("Report: report/01_qc_report.md\n\n")

if (requireNamespace("rmarkdown", quietly = TRUE)) {
  cat("Generating HTML report...\n")
  rmarkdown::render(here::here("scripts", "simple", "rna_01_qc_report.Rmd"),
                    output_dir  = file.path(RESULTS_DIR, "report"),
                    output_file = "01_qc_report.html",
                    quiet       = TRUE)
  cat("  → report/01_qc_report.html\n\n")
}

cat("================================================================\n")
cat("  QC completed —", OUT_PREFIX, "\n")
cat("================================================================\n\n")

