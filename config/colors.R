# colors.R
# Color palettes shared by the RNAseq_bulk scripts (all receptors and
# datasets).
# Usage: source(here::here("config", "colors.R")). If CELL_LINE and DRUG are
# defined before source(), OUT_PREFIX is also created.
# Unified system: condition_colors and annotation_colors$Condition are identical.
# Same criterion as scRNA-seq and scDNA-seq for comparability across analyses.

library(RColorBrewer)

# ------------------------------------------------------------------
# METADATA VARIABLES — a single color system
# ------------------------------------------------------------------

# Condition (main biological variable)
#   Parental   = gold      #e6a817  (same as Treatment=No)
#   Resistance = navy blue #013c58  (same as Treatment=Yes)
# Note: avoids a clash with the z-score gradient (blue-red)
condition_colors <- c(
  "Parental"   = "#e6a817",
  "Resistance" = "#013c58"
)

# Treatment (whether the sample was exposed to the drug)
#   Yes = salmon     #ef745c
#   No  = turquoise  #12CBC4
treatment_colors <- c(
  "Yes" = "#ef745c",
  "No"  = "#12CBC4"
)

# Cell lines
cellline_colors <- c(
  "ALK2" = "#6C5CE7",
  "ALK1" = "#00B894",
  "ALK3"  = "#FD79A8",
  "MET1"  = "#E17055",
  "RET1" = "#A55EEA",
  "RET2"  = "#F368E0"
)

# Drugs
drug_colors <- c(
  "Lorlatinib" = "#E17055",
  "Alectinib"  = "#0984E3",
  "BLU667"     = "#FA8072",
  "LOXO292"    = "#00CEC9",
  "No"         = "#DFE6E9"
)

# Gene classification
# ALK (5 classes): Acquired, Pre_existing, Maintained, Downregulated, Other
# RET (5 classes, same calculations, neutral names): Resistance_specific,
#   Present_in_Treated, Treatment_response, Downregulated, Other
# MET/MET1 (3 classes): Upregulated, Downregulated, Other
classification_colors <- c(
  Upregulated         = "#e67e22",
  Acquired            = "#e67e22",
  Resistance_specific = "#e67e22",
  Pre_existing        = "#e74c3c",
  Present_in_Treated  = "#e74c3c",
  Downregulated       = "#3498db",
  Maintained          = "#95a5a6",
  Treatment_response  = "#95a5a6",
  Other               = "#ecf0f1"
)

# ------------------------------------------------------------------
# HEATMAP GRADIENTS
# ------------------------------------------------------------------
heatmap_zscore    <- colorRampPalette(c("#3498DB", "white", "#E74C3C"))(100)
heatmap_distances <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
heatmap_standard  <- colorRampPalette(c("#2E86AB", "white", "#A23B72"))(100)

# ------------------------------------------------------------------
# GENERAL CATEGORICAL PALETTE
# ------------------------------------------------------------------
categorical_palette <- c(
  "#E74C3C", "#3498DB", "#2ECC71", "#F39C12", "#9B59B6",
  "#1ABC9C", "#34495E", "#E67E22", "#95A5A6", "#16A085"
)

# ------------------------------------------------------------------
# PHEATMAP ANNOTATIONS — synchronized with the vectors above
# ------------------------------------------------------------------
annotation_colors <- list(
  Condition        = condition_colors,
  Treatment        = treatment_colors,
  Cell_Line        = c(ALK2      = "#6C5CE7", ALK1      = "#00B894",
                       ALK3       = "#FD79A8", MET1       = "#E17055",
                       EGFR1       = "#74B9FF",
                       RET1      = "#A55EEA", RET2       = "#F368E0"),
  Drug             = c(Lorlatinib = "#E17055", Alectinib  = "#0984E3",
                       BLU667     = "#FA8072", LOXO292    = "#00CEC9",
                       No         = "#DFE6E9"),
  GX_cluster = c(Group1_EGFR_bypass = "#6C5CE7", Group2_NF2_YAP = "#00B894",
                 Parental           = "#e7d7c1")
)

# Filter annotation_colors to the columns present in an annotation data.frame
filter_annotation_colors <- function(annotation_df) {
  annotation_colors[names(annotation_colors) %in% colnames(annotation_df)]
}

if (exists("CELL_LINE") && exists("DRUG")) OUT_PREFIX <- paste0(CELL_LINE, "_", DRUG)

cat("colors.R loaded\n")
