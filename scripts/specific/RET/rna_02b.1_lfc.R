# ============================================================
# rna_02b.1_lfc.R
# LFC per drug + combined: family with treated parental (RET)
# No 5-category classification: 4 explicit LFC comparisons
# RNAseq_bulk RET resistance pipeline
#
# Inputs: results/RET/<CELL_LINE>/tables/01_mat_tmm.csv and
#   01_metadata_unit.csv; config/colors.R, config/gene_lists.R
# Outputs (results/RET/<CELL_LINE>/): tables/02b_lfc_by_drug.csv,
#   tables/02b_lfc_combined.csv, figures/02b_*.png, report/02b_lfc_report.md
# Dependencies (run first): scripts/simple/rna_01.1_qc.R with the same
#   DATASET. One run per model (CELL_LINE = "RET1" | "RET2").
# Usage: paths are resolved with here::here() from the repository root
# (.here file), so it works from any working directory:
#   Rscript scripts/specific/RET/rna_02b.1_lfc.R
# ============================================================
# CONFIG: change here for each run
# ============================================================
DATASET         <- "RET1"   # "RET1" | "RET2"
CELL_LINE       <- "RET1"
RECEPTOR        <- "RET"
LOGFC_THRESHOLD <- 2.0
TOP_N_HEATMAP   <- 25   # top N Upregulated + top N Downregulated (combined signature)
# ============================================================
RESULTS_DIR <- here::here("results", RECEPTOR, CELL_LINE)
CONFIG_DIR  <- here::here("config")
DRUG <- RECEPTOR  # placeholder: colors.R computes OUT_PREFIX from DRUG, but there is no single drug here
source(file.path(CONFIG_DIR, "colors.R"))
source(file.path(CONFIG_DIR, "gene_lists.R"))
OUT_PREFIX <- paste0(CELL_LINE, "_", RECEPTOR)  # overrides the one from colors.R

# LIBRARIES ------------------------------------------------------------------
library(ggplot2)
library(pheatmap)
library(dplyr)
library(tidyr)
library(tibble)
library(VennDiagram)
library(grid)

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

classify_lfc <- function(lfc, threshold) {
  cls <- rep("Other", length(lfc))
  cls[lfc >  threshold] <- "Upregulated"
  cls[lfc < -threshold] <- "Downregulated"
  cls
}

cat("\n================================================================\n")
cat("  LFC (treated parental) —", OUT_PREFIX, "\n")
cat("================================================================\n\n")

# LOAD TABLES FROM rna_01 ----------------------------------------------------
mat_file  <- file.path(RESULTS_DIR, "tables", "01_mat_tmm.csv")
meta_file <- file.path(RESULTS_DIR, "tables", "01_metadata_unit.csv")

if (!file.exists(mat_file)) stop("01_mat_tmm.csv not found; run rna_01.1_qc.R first")

mat_df <- read.csv(mat_file, check.names = FALSE)
rownames(mat_df) <- mat_df$Gene
mat <- as.matrix(mat_df[, -1])

meta_df <- read.csv(meta_file, check.names = FALSE)
rownames(meta_df) <- meta_df$Model
metadata_unit <- meta_df[, -1]

cat("Tables 01_mat_tmm.csv | 01_metadata_unit.csv loaded —",
    nrow(mat), "genes,", ncol(mat), "samples\n\n")

# IDENTIFY GROUPS ----------------------------------------------------------
cat("1. Identifying groups...\n")

untreated <- rownames(metadata_unit)[
  metadata_unit$Condition == "Parental" & metadata_unit$Treatment == "No"]
drugs <- setdiff(unique(metadata_unit$Drug), "No")

if (length(untreated) == 0) stop("Untreated parental sample not found.")
if (length(drugs) == 0)     stop("No drugs found in the Drug column.")

cat("  Untreated: ", paste(untreated, collapse = ", "), "\n")
cat("  Drugs: ", paste(drugs, collapse = ", "), "\n\n")

expr_unt <- if (length(untreated) > 1) rowMeans(mat[, untreated]) else mat[, untreated]

# LFC PER DRUG (b, c, d) ---------------------------------------------------
cat("2. Computing LFC per drug (b: Unt-vs-Trt, c: Trt-vs-Res, d: Unt-vs-Res)...\n")

per_drug_expr <- list()  # stores expr_trt_D / expr_res_D to reuse in the combined analysis

lfc_by_drug <- lapply(drugs, function(D) {
  treated_D   <- rownames(metadata_unit)[
    metadata_unit$Condition == "Parental"  & metadata_unit$Treatment == "Yes" & metadata_unit$Drug == D]
  resistant_D <- rownames(metadata_unit)[
    metadata_unit$Condition == "Resistance" & metadata_unit$Drug == D]

  if (length(treated_D) == 0)   stop("Treated parental not found for ", D)
  if (length(resistant_D) == 0) stop("Resistant samples not found for ", D)

  expr_trt_D <- if (length(treated_D) > 1) rowMeans(mat[, treated_D]) else mat[, treated_D]
  expr_res_D <- if (length(resistant_D) > 1) rowMeans(mat[, resistant_D]) else mat[, resistant_D]

  per_drug_expr[[D]] <<- list(expr_trt = expr_trt_D, expr_res = expr_res_D)

  lfc_unt_trt <- expr_trt_D - expr_unt
  lfc_trt_res <- expr_res_D - expr_trt_D
  lfc_unt_res <- expr_res_D - expr_unt

  cat("  ", D, "— treated:", paste(treated_D, collapse = ","),
      "| resistant:", paste(resistant_D, collapse = ","), "\n")

  tibble(
    Gene            = rownames(mat),
    Drug            = D,
    LFC_Unt_vs_Trt  = lfc_unt_trt,
    LFC_Trt_vs_Res  = lfc_trt_res,
    LFC_Unt_vs_Res  = lfc_unt_res,
    Mean_Expr_Unt   = expr_unt,
    Mean_Expr_Trt   = expr_trt_D,
    Mean_Expr_Res   = expr_res_D,
    Class_UntTrt    = classify_lfc(lfc_unt_trt, LOGFC_THRESHOLD),
    Class_TrtRes    = classify_lfc(lfc_trt_res, LOGFC_THRESHOLD),
    Class_UntRes    = classify_lfc(lfc_unt_res, LOGFC_THRESHOLD)
  )
}) %>% bind_rows()

cat("\n")

write.csv(lfc_by_drug, file.path(RESULTS_DIR, "tables", "02b_lfc_by_drug.csv"), row.names = FALSE)
cat("  → 02b_lfc_by_drug.csv (", nrow(lfc_by_drug), "rows —", length(drugs), "drugs x", nrow(mat), "genes)\n\n")

# COMBINED LFC (a) -------------------------------------------------------------
cat("3. Computing combined LFC (a: Unt-vs-Rest, drugs combined 50/50)...\n")

trt_means <- sapply(per_drug_expr, function(x) x$expr_trt)  # genes x drugs matrix
res_means <- sapply(per_drug_expr, function(x) x$expr_res)

treated_comb   <- if (length(drugs) > 1) rowMeans(trt_means) else trt_means[, 1]
resistant_comb <- if (length(drugs) > 1) rowMeans(res_means) else res_means[, 1]
rest_comb      <- (treated_comb + resistant_comb) / 2
lfc_unt_rest   <- rest_comb - expr_unt

lfc_combined <- tibble(
  Gene                 = rownames(mat),
  Mean_Expr_Unt        = expr_unt,
  Mean_Expr_Trt_comb   = treated_comb,
  Mean_Expr_Res_comb   = resistant_comb,
  Mean_Expr_Rest_comb  = rest_comb,
  LFC_Unt_vs_Rest      = lfc_unt_rest,
  Classification       = classify_lfc(lfc_unt_rest, LOGFC_THRESHOLD)
)

write.csv(lfc_combined, file.path(RESULTS_DIR, "tables", "02b_lfc_combined.csv"), row.names = FALSE)

tab_comb <- table(lfc_combined$Classification)
cat("  → 02b_lfc_combined.csv (", nrow(lfc_combined), "genes: ",
    tab_comb["Upregulated"], "up /", tab_comb["Downregulated"], "down)\n\n")

# FIGURES: HISTOGRAMS PER DRUG --------------------------------------------
cat("4. Histograms per drug...\n")

for (D in drugs) {
  df_d <- lfc_by_drug %>% filter(Drug == D) %>%
    select(Gene, LFC_Unt_vs_Trt, LFC_Trt_vs_Res, LFC_Unt_vs_Res) %>%
    pivot_longer(cols = -Gene, names_to = "Comparison", values_to = "LFC") %>%
    mutate(Comparison = factor(Comparison,
             levels = c("LFC_Unt_vs_Trt", "LFC_Trt_vs_Res", "LFC_Unt_vs_Res"),
             labels = c("Untreated vs Treated", "Treated vs Resistance", "Untreated vs Resistance")))

  p_hist_d <- ggplot(df_d, aes(LFC)) +
    geom_histogram(bins = 80, fill = "#74b9ff", color = "white", linewidth = 0.2) +
    geom_vline(xintercept = c(-LOGFC_THRESHOLD, LOGFC_THRESHOLD),
               linetype = "dashed", color = "#e17055") +
    facet_wrap(~ Comparison, ncol = 3) +
    labs(title = paste("LFC distribution —", D, "—", OUT_PREFIX),
         x = "log2 Fold Change", y = "Genes") +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold"),
          panel.background = element_rect(fill = "white"),
          plot.background  = element_rect(fill = "white"))

  ggsave(file.path(RESULTS_DIR, "figures", paste0("02b_lfc_histogram_", D, ".png")),
         p_hist_d, width = 12, height = 5, dpi = 300)
  cat("  →", paste0("02b_lfc_histogram_", D, ".png"), "\n")
}
cat("\n")

# FIGURE: COMBINED HISTOGRAM --------------------------------------------------
cat("5. Combined histogram...\n")

p_hist_comb <- ggplot(lfc_combined, aes(LFC_Unt_vs_Rest)) +
  geom_histogram(bins = 80, fill = "#74b9ff", color = "white", linewidth = 0.2) +
  geom_vline(xintercept = c(-LOGFC_THRESHOLD, LOGFC_THRESHOLD),
             linetype = "dashed", color = "#e17055") +
  labs(title = paste("LFC distribution — Untreated vs Rest (joint) —", OUT_PREFIX),
       x = "log2 Fold Change", y = "Genes") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"),
        panel.background = element_rect(fill = "white"),
        plot.background  = element_rect(fill = "white"))

ggsave(file.path(RESULTS_DIR, "figures", "02b_lfc_histogram_combined.png"),
       p_hist_comb, width = 8, height = 5, dpi = 300)
cat("  → 02b_lfc_histogram_combined.png\n\n")

# FIGURES: UP/DOWN VENN PER DRUG --------------------------------------------
cat("6. Venn diagrams (Up/Down, b vs c vs d) per drug...\n")

draw_venn3 <- function(set_b, set_c, set_d, labels, title, filepath) {
  n12  <- length(intersect(set_b, set_c))
  n23  <- length(intersect(set_c, set_d))
  n13  <- length(intersect(set_b, set_d))
  n123 <- length(intersect(intersect(set_b, set_c), set_d))

  venn.plot <- draw.triple.venn(
    area1 = length(set_b), area2 = length(set_c), area3 = length(set_d),
    n12 = n12, n23 = n23, n13 = n13, n123 = n123,
    category = labels,
    fill = c("#e67e22", "#9b59b6", "#3498db"),
    alpha = 0.4, cex = 1.3, cat.cex = 1.1, cat.fontface = "bold",
    margin = 0.1, ind = FALSE
  )
  png(filepath, width = 700, height = 650, res = 120)
  grid::grid.newpage()
  grid::grid.draw(venn.plot)
  grid::grid.text(title, y = 0.97, gp = grid::gpar(fontface = "bold", fontsize = 13))
  dev.off()
}

venn_labels <- c("Unt vs Trt", "Trt vs Res", "Unt vs Res")

for (D in drugs) {
  df_d <- lfc_by_drug %>% filter(Drug == D)

  up_b <- df_d$Gene[df_d$Class_UntTrt == "Upregulated"]
  up_c <- df_d$Gene[df_d$Class_TrtRes == "Upregulated"]
  up_d <- df_d$Gene[df_d$Class_UntRes == "Upregulated"]

  down_b <- df_d$Gene[df_d$Class_UntTrt == "Downregulated"]
  down_c <- df_d$Gene[df_d$Class_TrtRes == "Downregulated"]
  down_d <- df_d$Gene[df_d$Class_UntRes == "Downregulated"]

  tryCatch({
    draw_venn3(up_b, up_c, up_d, venn_labels,
               paste("Upregulated —", D),
               file.path(RESULTS_DIR, "figures", paste0("02b_venn_up_", D, ".png")))
    cat("  →", paste0("02b_venn_up_", D, ".png"), "\n")
  }, error = function(e) cat("  ⚠ Venn Up", D, "error:", conditionMessage(e), "\n"))

  tryCatch({
    draw_venn3(down_b, down_c, down_d, venn_labels,
               paste("Downregulated —", D),
               file.path(RESULTS_DIR, "figures", paste0("02b_venn_down_", D, ".png")))
    cat("  →", paste0("02b_venn_down_", D, ".png"), "\n")
  }, error = function(e) cat("  ⚠ Venn Down", D, "error:", conditionMessage(e), "\n"))
}
cat("\n")

# FIGURE: CLASSIFIED HEATMAP (combined, default pheatmap colors) -----
cat("7. Combined signature heatmap (default pheatmap colors)...\n")

is_clean_gene <- function(g) !grepl("\\.", g) & !grepl("^[0-9]", g)

sig_up   <- lfc_combined %>% filter(Classification == "Upregulated", is_clean_gene(Gene)) %>%
  arrange(desc(LFC_Unt_vs_Rest)) %>% slice_head(n = TOP_N_HEATMAP) %>% pull(Gene)
sig_down <- lfc_combined %>% filter(Classification == "Downregulated", is_clean_gene(Gene)) %>%
  arrange(LFC_Unt_vs_Rest) %>% slice_head(n = TOP_N_HEATMAP) %>% pull(Gene)
sig_genes <- c(sig_up, sig_down)

if (length(sig_genes) >= 3) {
  z_mat <- t(apply(mat[sig_genes, ], 1, cal_zscore))

  annotation <- metadata_unit[, c("Condition", "Drug"), drop = FALSE]
  ann_colors <- filter_annotation_colors(annotation)
  # Restrict each palette to the values present (prevents Lorlatinib/Alectinib
  # from other ALK units from appearing in the Drug legend of RET1/RET2)
  ann_colors <- setNames(lapply(names(ann_colors), function(nm) {
    vals <- ann_colors[[nm]]
    vals[names(vals) %in% unique(annotation[[nm]])]
  }), names(ann_colors))
  # RET-specific colors (same as the RET WES heatmap and oncoprint scripts,
  # COLOR_DRUG).
  if (!is.null(ann_colors$Drug)) {
    ann_colors$Drug <- c(
      "BLU667"  = "#6C5CE7",
      "LOXO292" = "#00B894",
      "No"      = "#95A5A6"
    )[names(ann_colors$Drug)]
  }

  ph_sig <- pheatmap::pheatmap(z_mat,
           # no `color=`: the default pheatmap palette is kept
           annotation_col    = annotation,
           annotation_colors = ann_colors,
           cluster_rows      = TRUE,
           cluster_cols      = TRUE,
           clustering_method = "ward.D2",
           show_rownames     = TRUE,
           show_colnames     = TRUE,
           fontsize_row      = 8,
           fontsize_col      = 9,
           main = paste0("Joint signature — top ", length(sig_up), " up / ",
                         length(sig_down), " down — ", OUT_PREFIX),
           silent = TRUE)

  save_pheatmap(ph_sig,
                file.path(RESULTS_DIR, "figures", "02b_heatmap_classified.png"),
                width = 11, height = 13)
  cat("  → 02b_heatmap_classified.png (", length(sig_up), "up /", length(sig_down), "down)\n\n")
} else {
  cat("  Fewer than 3 classified genes — heatmap skipped\n\n")
}

# REPORT ----------------------------------------------------------------------
cat("8. Generating report...\n")

report_text <- paste0(
  "# LFC Report (treated parental) — ", OUT_PREFIX, "\n\n",
  "**Drugs:** ", paste(drugs, collapse = ", "), "\n",
  "**Threshold:** |LFC| ≥ ", LOGFC_THRESHOLD, "\n\n",
  "## Comparisons per drug\n",
  "- b) Untreated vs Treated — acute response to the drug\n",
  "- c) Treated vs Resistance — additional change on reaching resistance\n",
  "- d) Untreated vs Resistance — total change\n\n",
  "## Joint comparison (a)\n",
  "Untreated vs Rest (50/50 mean of Treated_comb and Resistant_comb, each drug weighted equally)\n",
  "- Upregulated: ", tab_comb["Upregulated"], "\n",
  "- Downregulated: ", tab_comb["Downregulated"], "\n",
  "- Other: ", tab_comb["Other"], "\n\n",
  "## Tables\n",
  "- 02b_lfc_by_drug.csv (", nrow(lfc_by_drug), " rows)\n",
  "- 02b_lfc_combined.csv (", nrow(lfc_combined), " rows)\n\n",
  "## Figures\n",
  paste(paste0("- 02b_lfc_histogram_", drugs, ".png\n"), collapse = ""),
  "- 02b_lfc_histogram_combined.png\n",
  paste(paste0("- 02b_venn_up_", drugs, ".png\n"), collapse = ""),
  paste(paste0("- 02b_venn_down_", drugs, ".png\n"), collapse = ""),
  "- 02b_heatmap_classified.png (default pheatmap colors)\n"
)
writeLines(report_text, file.path(RESULTS_DIR, "report", "02b_lfc_report.md"))
cat("  Report: report/02b_lfc_report.md\n\n")

cat("================================================================\n")
cat("  LFC (treated parental) finished —", OUT_PREFIX, "\n")
cat("================================================================\n\n")
