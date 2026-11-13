# ============================================================
# signature_amplicon_egfr1.R
# Chr12 amplicon signatures (focal / non-focal) filtered by EGFR1 RNA-seq
# Reads 02_lfc_master_table.csv — does not recompute TMM.
# Input: results/EGFR/EGFR1/tables/02_lfc_master_table.csv
# Output: EGFR1_amplicon_focal.csv · EGFR1_amplicon_no_focal.csv
#   → results/EGFR/EGFR1/tables/
# Dependencies (run first): scripts/simple/rna_02.1_lfc.R with
#   DATASET = "EGFR1_EGFRi".
# ============================================================
# Usage: paths are resolved with here::here() from the repository root
# (.here file), so it works from any working directory:
#   Rscript scripts/specific/EGFR/signature_amplicon_egfr1.R
# ============================================================

RESULTS_DIR <- here::here("results", "EGFR", "EGFR1", "tables")
EXPR_CUTOFF <- -5   # minimum log2CPM in parental (same as rna_03.1)

library(dplyr)

# Focal zone: high amplification specific to EGFR1_R (log2FC ~4.1–4.3)
# Boundaries: SLC35E3 (5') → CCT2 (3') — WES CNV vs blood
focal_genes <- c(
  "SLC35E3", "MDM2", "CPM", "CPSF6", "LYZ", "YEATS4", "FRS2", "CCT2"
)

# Non-focal zone: the whole broad amplicon (log2FC ~2.0–2.7 in both samples)
# Boundaries: GRIP1 (5') → TRHDE (3') — GRIP1 included as the region boundary gene
no_focal_genes <- c(
  "GRIP1",
  "CAND1", "DYRK2", "IFNG", "IL26", "IL22", "MDM1",
  "RAP1B", "NUP107",
  "SLC35E3", "MDM2", "CPM", "CPSF6", "LYZ", "YEATS4", "FRS2", "CCT2",
  "LRRC10", "BEST3", "RAB3IP", "CNOT2", "KCNMB4",
  "PTPRB", "PTPRR", "TSPAN8", "LGR5", "ZFC3H1", "THAP2", "TMEM19",
  "RAB21", "TBC1D15", "TPH2", "TRHDE"
)

lfc_table <- read.csv(file.path(RESULTS_DIR, "02_lfc_master_table.csv"),
                      check.names = FALSE)

cat("\n================================================================\n")
cat("  Chr12 amplicon signatures — EGFR1 (EGFR)\n")
cat("  Genes in RNA-seq:", nrow(lfc_table), "\n")
cat("================================================================\n\n")

filter_signature <- function(gene_list, sig_name) {
  found <- lfc_table %>%
    filter(Gene %in% gene_list, Mean_Expr_Par > EXPR_CUTOFF) %>%
    arrange(match(Gene, gene_list)) %>%
    select(Gene)

  absent   <- setdiff(gene_list, lfc_table$Gene)
  low_expr <- setdiff(intersect(gene_list, lfc_table$Gene), found$Gene)

  cat(sprintf("  %s\n", sig_name))
  cat(sprintf("    Genes in list:            %d\n", length(gene_list)))
  cat(sprintf("    Expressed in RNA-seq:     %d\n", nrow(found)))
  if (length(absent)   > 0) cat("    Absent:     ", paste(absent,   collapse=", "), "\n")
  if (length(low_expr) > 0) cat("    Low expr.:  ", paste(low_expr, collapse=", "), "\n")
  cat("\n")

  found
}

sig_focal    <- filter_signature(focal_genes,    "FOCAL    (SLC35E3 → CCT2)")
sig_no_focal <- filter_signature(no_focal_genes, "NO FOCAL (GRIP1   → TRHDE)")

write.csv(sig_focal,    file.path(RESULTS_DIR, "EGFR1_amplicon_focal.csv"),    row.names=FALSE)
write.csv(sig_no_focal, file.path(RESULTS_DIR, "EGFR1_amplicon_no_focal.csv"), row.names=FALSE)

cat("  Saved:\n")
cat("   →", file.path(RESULTS_DIR, "EGFR1_amplicon_focal.csv"),    "\n")
cat("   →", file.path(RESULTS_DIR, "EGFR1_amplicon_no_focal.csv"), "\n")
cat("================================================================\n\n")
