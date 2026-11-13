# ============================================================
# rna_03b.1_key_genes.R
# Key genes + resistance signature: RET (RET1/RET2)
# "Treated parental" family; complements rna_02b.1_lfc.R
#
# Signature based on d) Untreated vs Resistance (full baseline-to-resistant
# comparison, consistent with the rest of the pipeline: MET1/EGFR1_EGFRi also
# use Resistance vs Parental without subtracting anything).
#
# To minimize the acute drug effect, a Resistance_specific column is added:
# TRUE only if the gene ALSO exceeds the threshold in c) Treated vs Resistance
# with the same sign (same LOGFC_THRESHOLD, no relaxed criterion). A gene that
# only changes in d) but not in c) is, mathematically (d = b + c), a change
# inherited almost entirely from b) Untreated vs Treated (acute response), not
# from resistance itself.
#
# No forced genes and no signature for scRNA-seq (there is no single-cell data
# for RET).
# RNAseq_bulk pipeline
#
# Inputs: results/RET/<CELL_LINE>/tables/01_mat_tmm.csv,
#   01_metadata_unit.csv, 02b_lfc_combined.csv; config/colors.R,
#   config/gene_lists.R
# Outputs (results/RET/<CELL_LINE>/): tables/03b_lfc_combined_4comparisons.csv,
#   tables/RET_resistance_signature.csv,
#   tables/RET_resistance_signature_conservative.csv, figures/03b_*.png,
#   report/03b_key_genes_report.md
# Dependencies (run first): rna_02b.1_lfc.R for the same CELL_LINE.
#   One run per model (CELL_LINE = "RET1" | "RET2").
# Usage: paths are resolved with here::here() from the repository root
# (.here file), so it works from any working directory:
#   Rscript scripts/specific/RET/rna_03b.1_key_genes.R
# ============================================================
# CONFIG: change here for each run
# ============================================================
DATASET   <- "RET1"   # "RET1" | "RET2"
CELL_LINE <- "RET1"
RECEPTOR  <- "RET"
LOGFC_THRESHOLD <- 2.0
TOP_N_HEATMAP   <- 15   # number of up/down genes in the signature heatmap (top15+15=30)
TOP_N_VIS       <- 5    # number of up/down genes in the stripplot
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

cal_zscore    <- function(x) if (sd(x) == 0) rep(0, length(x)) else (x - mean(x)) / sd(x)
is_clean_gene <- function(g) !grepl("\\.", g) & !grepl("^[0-9]", g)

cat("\n================================================================\n")
cat("  KEY GENES —", OUT_PREFIX, "\n")
cat("================================================================\n\n")

# LOAD TABLES ---------------------------------------------------------------
mat_file      <- file.path(RESULTS_DIR, "tables", "01_mat_tmm.csv")
meta_file     <- file.path(RESULTS_DIR, "tables", "01_metadata_unit.csv")
lfc_comb_file <- file.path(RESULTS_DIR, "tables", "02b_lfc_combined.csv")

if (!file.exists(mat_file))      stop("01_mat_tmm.csv not found; run rna_01")
if (!file.exists(lfc_comb_file)) stop("02b_lfc_combined.csv not found; run rna_02b")

mat_df           <- read.csv(mat_file, check.names = FALSE)
rownames(mat_df) <- mat_df$Gene
mat              <- as.matrix(mat_df[, -1])

meta_df           <- read.csv(meta_file, check.names = FALSE)
rownames(meta_df) <- meta_df$Model
metadata_unit     <- meta_df[, -1]

lfc_combined <- read.csv(lfc_comb_file, check.names = FALSE)

untreated <- rownames(metadata_unit)[metadata_unit$Condition == "Parental" & metadata_unit$Treatment == "No"]
drugs     <- setdiff(unique(metadata_unit$Drug), "No")

treated_order  <- unlist(lapply(drugs, function(d)
  rownames(metadata_unit)[metadata_unit$Condition == "Parental" &
                          metadata_unit$Treatment == "Yes" & metadata_unit$Drug == d]))
resistant_order <- unlist(lapply(drugs, function(d)
  rownames(metadata_unit)[metadata_unit$Condition == "Resistance" & metadata_unit$Drug == d]))
sample_order <- c(untreated, treated_order, resistant_order)

cat("Tables loaded: 01_mat_tmm | 01_metadata_unit | 02b_lfc_combined\n")
cat("Untreated:", untreated, "| Drugs:", paste(drugs, collapse = ", "), "\n\n")

annotation_col <- metadata_unit[, c("Condition", "Drug"), drop = FALSE]
ann_colors     <- filter_annotation_colors(annotation_col)
# Restrict each palette to the values present (prevents Lorlatinib/Alectinib
# from other ALK units from appearing in the Drug legend of RET1/RET2)
ann_colors <- setNames(lapply(names(ann_colors), function(nm) {
  vals <- ann_colors[[nm]]
  vals[names(vals) %in% unique(annotation_col[[nm]])]
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

# THE 4 LFC COMPARISONS (b/c/d derived from columns already in 02b_lfc_combined)
# b/c/d are derived from Mean_Expr_Trt_comb / Mean_Expr_Res_comb, with the same
# 50/50 weighting already computed in 02b, without recomputing or inventing new
# categories.
cat("1. Computing the 4 comparisons (a/b/c/d)...\n")

lfc4 <- lfc_combined %>%
  mutate(
    LFC_Unt_vs_Trt_comb = Mean_Expr_Trt_comb - Mean_Expr_Unt,
    LFC_Trt_vs_Res_comb = Mean_Expr_Res_comb - Mean_Expr_Trt_comb,
    LFC_Unt_vs_Res_comb = Mean_Expr_Res_comb - Mean_Expr_Unt,
    Class_a = Classification,
    Class_b = case_when(
      LFC_Unt_vs_Trt_comb >  LOGFC_THRESHOLD ~ "Upregulated",
      LFC_Unt_vs_Trt_comb < -LOGFC_THRESHOLD ~ "Downregulated",
      TRUE ~ "Other"),
    Class_c = case_when(
      LFC_Trt_vs_Res_comb >  LOGFC_THRESHOLD ~ "Upregulated",
      LFC_Trt_vs_Res_comb < -LOGFC_THRESHOLD ~ "Downregulated",
      TRUE ~ "Other"),
    Class_d = case_when(
      LFC_Unt_vs_Res_comb >  LOGFC_THRESHOLD ~ "Upregulated",
      LFC_Unt_vs_Res_comb < -LOGFC_THRESHOLD ~ "Downregulated",
      TRUE ~ "Other"),
    # TRUE only if it exceeds the threshold in c) with the SAME sign as d);
    # same LOGFC_THRESHOLD as the rest of the script, no relaxed criterion.
    Resistance_specific = Class_d != "Other" & Class_c == Class_d
  )

write.csv(lfc4 %>% dplyr::select(Gene, Mean_Expr_Unt, Mean_Expr_Trt_comb, Mean_Expr_Res_comb,
                                  LFC_Unt_vs_Trt_comb, Class_b,
                                  LFC_Trt_vs_Res_comb, Class_c,
                                  LFC_Unt_vs_Res_comb, Class_d, Resistance_specific,
                                  LFC_Unt_vs_Rest, Class_a),
          file.path(RESULTS_DIR, "tables", "03b_lfc_combined_4comparisons.csv"),
          row.names = FALSE)
cat("  → 03b_lfc_combined_4comparisons.csv\n\n")

# KEY GENE SELECTION (based on d) --------------------------------------
cat("2. Selecting key genes (d) classification, LFC_Unt_vs_Res_comb)...\n")

up_genes   <- lfc4 %>% filter(Class_d == "Upregulated") %>%
  arrange(desc(LFC_Unt_vs_Res_comb))
down_genes <- lfc4 %>% filter(Class_d == "Downregulated") %>%
  arrange(LFC_Unt_vs_Res_comb)

cat("  Upregulated:  ", nrow(up_genes), " (", sum(up_genes$Resistance_specific), "Resistance_specific)\n")
cat("  Downregulated:", nrow(down_genes), " (", sum(down_genes$Resistance_specific), "Resistance_specific)\n\n")

# RESISTANCE SIGNATURE (data-driven, no forced genes) ----------------------
cat("3. Generating resistance signature...\n")

signature <- bind_rows(up_genes, down_genes) %>%
  dplyr::select(Gene, Mean_Expr_Unt, Mean_Expr_Trt_comb, Mean_Expr_Res_comb,
                LFC_Unt_vs_Res_comb, Class_d, Resistance_specific,
                LFC_Trt_vs_Res_comb) %>%
  rename(Classification = Class_d)

write.csv(signature,
          file.path(RESULTS_DIR, "tables", paste0(RECEPTOR, "_resistance_signature.csv")),
          row.names = FALSE)
cat("  →", paste0(RECEPTOR, "_resistance_signature.csv"), "(", nrow(signature), "genes: ",
    sum(signature$Classification == "Upregulated"), "up /",
    sum(signature$Classification == "Downregulated"), "down /",
    sum(signature$Resistance_specific), "Resistance_specific)\n")

# Conservative: expressed above noise in the untreated parental (baseline)
signature_cons <- signature %>% filter(Mean_Expr_Unt > -5)
write.csv(signature_cons,
          file.path(RESULTS_DIR, "tables",
                    paste0(RECEPTOR, "_resistance_signature_conservative.csv")),
          row.names = FALSE)
cat("  →", paste0(RECEPTOR, "_resistance_signature_conservative.csv"), "(", nrow(signature_cons),
    "genes: ", sum(signature_cons$Classification == "Upregulated"), "up /",
    sum(signature_cons$Classification == "Downregulated"), "down /",
    sum(signature_cons$Resistance_specific), "Resistance_specific)\n\n")

# HEATMAP: SIGNATURE d) top30 --------------------------------------------------
cat("4. Heatmap d) signature Untreated vs Resistance (top", 2 * TOP_N_HEATMAP, ")...\n")

sig_up_all   <- signature_cons %>% filter(Classification == "Upregulated") %>%
  arrange(desc(LFC_Unt_vs_Res_comb)) %>% pull(Gene)
sig_down_all <- signature_cons %>% filter(Classification == "Downregulated") %>%
  arrange(LFC_Unt_vs_Res_comb) %>% pull(Gene)

sig_up_all   <- sig_up_all[sig_up_all %in% rownames(mat) & is_clean_gene(sig_up_all)]
sig_down_all <- sig_down_all[sig_down_all %in% rownames(mat) & is_clean_gene(sig_down_all)]

top_up   <- head(sig_up_all, TOP_N_HEATMAP)
top_down <- head(sig_down_all, TOP_N_HEATMAP)
sig_ordered <- c(top_up, top_down)

if (length(sig_ordered) >= 3) {
  z_sig <- t(apply(mat[sig_ordered, ], 1, cal_zscore))

  ph_d <- pheatmap::pheatmap(z_sig,
           annotation_col    = annotation_col,
           annotation_colors = ann_colors,
           cluster_rows      = TRUE,
           cluster_cols      = TRUE,
           clustering_method = "ward.D2",
           show_rownames     = TRUE,
           show_colnames     = TRUE,
           fontsize_row      = 8,
           fontsize_col      = 9,
           main = paste0("d) Untreated vs Resistance — ", OUT_PREFIX,
                         " (", length(top_up), " up / ", length(top_down), " down)"),
           silent = TRUE)

  save_pheatmap(ph_d, file.path(RESULTS_DIR, "figures", "03b_heatmap_d_signature.png"),
                width = 12, height = 14)
  cat("  → 03b_heatmap_d_signature.png (", length(top_up), "up /", length(top_down), "down)\n\n")
} else {
  cat("  Fewer than 3 genes in the signature — heatmap skipped\n\n")
}

# GENE SELECTION FOR THE SIGNATURE STRIPPLOT (d) -------------------------------
cat("5. Selecting genes for the signature stripplot (top", TOP_N_VIS, "up/down by LFC d)...\n")

vis_up   <- head(up_genes$Gene[up_genes$Gene %in% rownames(mat)], TOP_N_VIS)
vis_down <- head(down_genes$Gene[down_genes$Gene %in% rownames(mat)], TOP_N_VIS)
vis_genes <- unique(c(vis_up, vis_down))

cat("  Top", TOP_N_VIS, "Up:  ", paste(vis_up,   collapse = ", "), "\n")
cat("  Top", TOP_N_VIS, "Down:", paste(vis_down, collapse = ", "), "\n\n")

# FIGURE: SIGNATURE STRIPPLOT (d) --------------------------------------------------
cat("6. Signature stripplot (d)...\n")

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
  scale_shape_manual(values = c(No = 15, setNames(c(16, 17, 18, 3, 4)[seq_along(drugs)], drugs))) +
  facet_wrap(~ Gene, scales = "free_y", ncol = 5) +
  labs(
    title    = paste("Expression per sample —", OUT_PREFIX),
    subtitle = paste0("Top", TOP_N_VIS, " Up: ", paste(vis_up, collapse = ", "),
                       " | Top", TOP_N_VIS, " Down: ", paste(vis_down, collapse = ", ")),
    x = NULL, y = "log2CPM (TMM)"
  ) +
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

n_cols <- min(5, length(vis_genes))
n_rows <- ceiling(length(vis_genes) / n_cols)
ggsave(file.path(RESULTS_DIR, "figures", "03b_stripplot_key_genes.png"),
       p_strip, width = n_cols * 2.5, height = n_rows * 3, dpi = 300)
cat("  → 03b_stripplot_key_genes.png (", length(vis_genes), "genes)\n\n")

# REPORT ----------------------------------------------------------------------
cat("7. Generating report...\n")

report_text <- paste0(
  "# Key Genes Report — ", OUT_PREFIX, "\n\n",
  "**Signature:** data-driven, based on d) Untreated vs Resistance (full comparison ",
  "baseline→resistant, consistent with MET1/EGFR1_EGFRi). Column `Resistance_specific` = TRUE ",
  "if the gene also passes the threshold in c) Treated vs Resistance with the same sign ",
  "(same LOGFC_THRESHOLD=", LOGFC_THRESHOLD, ") — filters out genes whose change in resistance ",
  "is almost entirely inherited from the acute response to the drug (b).\n\n",
  "**Note:** no signature is generated for scRNA-seq — there is no RET single-cell data yet.\n\n",
  "## Top 15 Upregulated (d)\n",
  paste(capture.output(print(
    as.data.frame(slice_head(up_genes, n = 15) %>%
                    dplyr::select(Gene, LFC_Unt_vs_Res_comb, Resistance_specific)),
    row.names = FALSE)), collapse = "\n"), "\n\n",
  "## Top 15 Downregulated (d)\n",
  paste(capture.output(print(
    as.data.frame(slice_head(down_genes, n = 15) %>%
                    dplyr::select(Gene, LFC_Unt_vs_Res_comb, Resistance_specific)),
    row.names = FALSE)), collapse = "\n"), "\n\n",
  "## Signature ", RECEPTOR, "_resistance_signature_conservative.csv\n",
  "- N genes: ", nrow(signature_cons), " (",
  sum(signature_cons$Classification == "Upregulated"), " up / ",
  sum(signature_cons$Classification == "Downregulated"), " down / ",
  sum(signature_cons$Resistance_specific), " Resistance_specific)\n\n",
  "## Genes stripplot (top", TOP_N_VIS, "by LFC d, no known list)\n",
  "- Top", TOP_N_VIS, " Up:  ", paste(vis_up,   collapse = ", "), "\n",
  "- Top", TOP_N_VIS, " Down:", paste(vis_down, collapse = ", "), "\n\n",
  "## Figures\n",
  "- 03b_heatmap_d_signature.png — signature d) top", TOP_N_HEATMAP, " up + top", TOP_N_HEATMAP,
  " down, annotated with Resistance_specific (corroboration in c)\n",
  "- 03b_stripplot_key_genes.png — expression per sample (top", TOP_N_VIS, " up/down of d)\n\n",
  "**WES/scDNA/SNV panel:** see rna_03c.1_wes_panel_combined.R (heatmaps + CNV/SNV stripplots)\n"
)

writeLines(report_text, file.path(RESULTS_DIR, "report", "03b_key_genes_report.md"))
cat("  Report: report/03b_key_genes_report.md\n\n")

cat("================================================================\n")
cat("  KEY GENES finished —", OUT_PREFIX, "\n")
cat("================================================================\n\n")
