# ============================================================
# rna_03.1_key_genes.R
# Key genes and transcriptomic resistance signature: cross-referencing of the
# classified genes with reference lists (config/gene_lists.R), full and
# conservative signatures, heatmaps, stripplot and lineplot of selected genes.
# Generic: 2-class scheme (TWO_CLASS_DATASETS) or 5-class scheme.
#
# Inputs:    results/<RECEPTOR>/<CELL_LINE>/tables/01_mat_tmm.csv,
#              01_metadata_unit.csv, 02_lfc_master_table.csv
#            data/counts_<DATASET>.csv (2-class scheme only, off->on genes)
#            config/colors.R, config/gene_lists.R
# Outputs:   results/<RECEPTOR>/<CELL_LINE>/tables/<PREFIX>_resistance_signature.csv,
#              <PREFIX>_resistance_signature_conservative.csv (input of rna_04)
#            results/<RECEPTOR>/<CELL_LINE>/figures/03_*.png
#            results/<RECEPTOR>/<CELL_LINE>/report/03_key_genes_report.{md,html}
# Usage:     set the CONFIG block and, from any folder of the repo,
#            source(here::here("scripts/simple/rna_03.1_key_genes.R"))
# ============================================================
# CONFIG — change here for each run
# ============================================================
DATASET         <- "EGFR1_EGFRi"      # "ALK12" | "ALK3" | "MET1_METi" | "EGFR1_EGFRi"
CELL_LINE       <- "EGFR1"           # "ALK2" | "ALK1" | "ALK3" | "MET1" | "EGFR1"
DRUG            <- "Gefitinib_Afatinib"  # "Lorlatinib" | "Alectinib" | "MET_inhibitor" | "Gefitinib_Afatinib"
RECEPTOR        <- "EGFR"            # "MET" | "ALK" | "EGFR"
LOGFC_THRESHOLD <- 2.0

# Genes forced into the signature regardless of LFC
FORCED_GENES <- if (exists("DATASET") && DATASET == "EGFR1_EGFRi") {
  c("RB1", "FGFR3", "MDM2", "TACC3")          # RB1 inactivated, EGFR1_R FGFR3-TACC3 fusion
} else {
  c("NF2", "NRG1")            # Giménez-Xavier 2017 — MET1_METi
}

# Fixed genes for stripplot and lineplot
GENES_FIXED <- if (exists("DATASET") && DATASET == "EGFR1_EGFRi") {
  c("EGFR", "RB1", "FGFR3", "TACC3", "MYCL", "MDM2", "TP53")
} else {
  c("MET", "EGFR", "NF2", "NRG1")     # MET1_METi
}
# ============================================================
RESULTS_DIR <- here::here("results", RECEPTOR, CELL_LINE)
CONFIG_DIR  <- here::here("config")
source(file.path(CONFIG_DIR, "colors.R"))
source(file.path(CONFIG_DIR, "gene_lists.R"))

# Category names of the 5-class scheme (datasets outside TWO_CLASS_DATASETS;
# must match rna_02.1_lfc.R)
CLASS_LABELS <- c(down = "Downregulated", pre = "Pre_existing",
                   maint = "Maintained", acq = "Acquired")

# Datasets with the simplified scheme (Resistance vs Parental, no treated/untreated distinction)
TWO_CLASS_DATASETS <- c("MET1_METi", "EGFR1_EGFRi")

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

cal_zscore <- function(x) if (sd(x) == 0) rep(0, length(x)) else (x - mean(x)) / sd(x)

# Gradient for the signature heatmap (03_heatmap_signature.png): black -> gray
# -> gold. heatmap_zscore (colors.R) is left untouched and is still used for
# the rest of the heatmaps (Giménez-Xavier, etc.)
GRADIENTE_FIRMA <- c("#000000", "#808080", "#f2a900")
col_fun_firma   <- colorRampPalette(GRADIENTE_FIRMA)(100)

cat("\n================================================================\n")
cat("  KEY GENES —", OUT_PREFIX, "\n")
cat("================================================================\n\n")

# LOAD TABLES --------------------------------------------------------------
mat_file  <- file.path(RESULTS_DIR, "tables", "01_mat_tmm.csv")
meta_file <- file.path(RESULTS_DIR, "tables", "01_metadata_unit.csv")
lfc_file  <- file.path(RESULTS_DIR, "tables", "02_lfc_master_table.csv")

if (!file.exists(mat_file))  stop("01_mat_tmm.csv not found — run rna_01")
if (!file.exists(lfc_file))  stop("02_lfc_master_table.csv not found — run rna_02")

mat_df        <- read.csv(mat_file, check.names = FALSE)
rownames(mat_df) <- mat_df$Gene
mat           <- as.matrix(mat_df[, -1])

meta_df       <- read.csv(meta_file, check.names = FALSE)
rownames(meta_df) <- meta_df$Model
metadata_unit <- meta_df[, -1]

lfc_table <- read.csv(lfc_file, check.names = FALSE)

if (DATASET == "MET1_METi" && !"GX_cluster" %in% colnames(metadata_unit)) {
  rg_map <- c(
    MET1    = "Parental",
    MET1_C1 = "Group1_EGFR_bypass", MET1_PR1 = "Group1_EGFR_bypass",
    MET1_PR2 = "Group1_EGFR_bypass", MET1_V8  = "Group1_EGFR_bypass",
    MET1_PR3 = "Group2_NF2_YAP",    MET1_R = "Group2_NF2_YAP"
  )
  metadata_unit$GX_cluster <- rg_map[rownames(metadata_unit)]
}

parental  <- rownames(metadata_unit)[metadata_unit$Condition == "Parental"]
resistant <- rownames(metadata_unit)[metadata_unit$Condition == "Resistance"]

cat("Tables loaded: 01_mat_tmm | 01_metadata_unit | 02_lfc_master_table\n\n")

# GENE SELECTION ---------------------------------------------------------
cat("1. Selecting key genes...\n")

if (DATASET %in% TWO_CLASS_DATASETS) {
  up_genes   <- lfc_table %>% filter(Classification == "Upregulated") %>%
    arrange(desc(LFC_Res_vs_Par))
  down_genes <- lfc_table %>% filter(Classification == "Downregulated") %>%
    arrange(LFC_Res_vs_Par)

  cat("  Upregulated:  ", nrow(up_genes), "\n")
  cat("  Downregulated:", nrow(down_genes), "\n\n")

} else {
  up_genes   <- lfc_table %>% filter(Classification == unname(CLASS_LABELS["acq"])) %>%
    arrange(desc(LFC_Res_vs_Trt))
  pre_genes  <- lfc_table %>% filter(Classification == unname(CLASS_LABELS["pre"])) %>%
    arrange(desc(abs(LFC_Res_vs_Trt)))
  down_genes <- lfc_table %>% filter(Classification == "Downregulated") %>%
    arrange(LFC_Res_vs_Unt)

  cat("  ", CLASS_LABELS["acq"], ":    ", nrow(up_genes), "\n", sep = "")
  cat("  ", CLASS_LABELS["pre"], ":", nrow(pre_genes), "\n", sep = "")
  cat("  Downregulated:", nrow(down_genes), "\n\n")
}

# CROSS-REFERENCING WITH KNOWN LISTS ------------------------------------------------
cat("2. Cross-referencing with known gene lists...\n")

if (DATASET == "MET1_METi") {
  lists_to_check <- list(
    met_resistance_core  = met_resistance_core,
    gimenez_common       = gimenez_common,
    gimenez_egfr_bypass  = gimenez_egfr_bypass,
    gimenez_group2       = gimenez_group2,
    met_core             = met_core,
    met_bypass           = met_bypass
  )
  for (nm in names(lists_to_check)) {
    hits_up   <- intersect(up_genes$Gene,   lists_to_check[[nm]])
    hits_down <- intersect(down_genes$Gene, lists_to_check[[nm]])
    if (length(hits_up) > 0 || length(hits_down) > 0) {
      cat("  [", nm, "]\n")
      if (length(hits_up)   > 0) cat("    ↑ Up:  ", paste(hits_up,   collapse = ", "), "\n")
      if (length(hits_down) > 0) cat("    ↓ Down:", paste(hits_down, collapse = ", "), "\n")
    }
  }
} else if (DATASET == "EGFR1_EGFRi") {
  hits_up   <- intersect(up_genes$Gene,   egfr_resistance_core)
  hits_down <- intersect(down_genes$Gene, egfr_resistance_core)
  cat("  [egfr_resistance_core]\n")
  if (length(hits_up)   > 0) cat("    ↑ Up:  ", paste(hits_up,   collapse = ", "), "\n")
  if (length(hits_down) > 0) cat("    ↓ Down:", paste(hits_down, collapse = ", "), "\n")
} else {
  core_acq  <- intersect(up_genes$Gene,   alk_resistance_core)
  core_pre  <- intersect(pre_genes$Gene,  alk_resistance_core)
  core_down <- intersect(down_genes$Gene, alk_resistance_core)
  cat("  alk_resistance_core — Acquired:", paste(core_acq,  collapse = ", "), "\n")
  cat("  alk_resistance_core — Pre:     ", paste(core_pre,  collapse = ", "), "\n")
  cat("  alk_resistance_core — Down:    ", paste(core_down, collapse = ", "), "\n")
}
cat("\n")

# Annotation for heatmaps
if (DATASET == "MET1_METi") {
  annotation_col <- metadata_unit[, c("Condition", "GX_cluster"), drop = FALSE]
} else {
  annotation_col <- metadata_unit[,
    !colnames(metadata_unit) %in% c("Sample", "Model", "Cell_Line"), drop = FALSE]
}
ann_colors <- filter_annotation_colors(annotation_col)

is_clean_gene <- function(g) !grepl("\\.", g) & !grepl("^[0-9]", g)

# FIGURE 1: GIMÉNEZ-XAVIER HEATMAP -------------------------------------------
if (DATASET == "MET1_METi") {
  cat("3. Giménez-Xavier heatmap (Figure 1)...\n")

  gx_genes <- unique(c(gimenez_common, gimenez_egfr_bypass, gimenez_group2))
  gx_genes <- gx_genes[gx_genes %in% rownames(mat)]

  if (length(gx_genes) >= 3) {
    z_gx <- t(apply(mat[gx_genes, ], 1, cal_zscore))
    ph_gx <- pheatmap(z_gx,
             color             = heatmap_zscore,
             annotation_col    = annotation_col,
             annotation_colors = ann_colors,
             cluster_rows      = TRUE,
             cluster_cols      = TRUE,
             clustering_method = "ward.D2",
             cutree_rows       = 3,
             show_rownames     = TRUE,
             show_colnames     = TRUE,
             fontsize_row      = 8,
             fontsize_col      = 9,
             main = paste0("Giménez-Xavier 2017 gene lists — ", OUT_PREFIX),
             silent = TRUE)

    grid::grid.newpage()
    grid::grid.draw(ph_gx$gtable)
    save_pheatmap(ph_gx,
                  file.path(RESULTS_DIR, "figures", "03_heatmap_gimenez.png"),
                  width = 14, height = 9)
    cat("  → 03_heatmap_gimenez.png\n\n")
  } else {
    cat("  Fewer than 3 Giménez-Xavier genes in the matrix — skipped\n\n")
  }
}

# RESISTANCE SIGNATURE (OUTPUT → scRNA-seq) ----------------------------------
cat("4. Generating resistance signature...\n")

sig_prefix <- if (DATASET == "MET1_METi") "MET" else RECEPTOR

if (DATASET %in% TWO_CLASS_DATASETS) {
  # "Off in parental -> on in resistant" genes: 0 raw reads in the
  # parental and substantial expression (log2CPM > 2) in the resistant models.
  # The conservative filter below (Mean_Expr_Par > -5) would exclude them
  # because they sit at the noise floor in the parental, but unlike the rest
  # of the genes at that noise floor, these are clearly and consistently
  # switched on in the resistant models (not noise against noise).
  off_on_genes <- character(0)
  counts_file  <- here::here("data", paste0("counts_", DATASET, ".csv"))
  if (file.exists(counts_file)) {
    sep_counts <- if (grepl(";", readLines(counts_file, n = 1))) ";" else ","
    counts_raw <- read.csv(counts_file, sep = sep_counts, check.names = FALSE)

    # Count-matrix columns are named by Sample, the parental by Model
    parental_col <- metadata_unit[parental, "Sample"]

    off_on_genes <- up_genes %>%
      filter(Mean_Expr_Par <= -5, Mean_Expr_Res > 2) %>%
      inner_join(counts_raw[, c("gene_name", parental_col)], by = c("Gene" = "gene_name")) %>%
      filter(.data[[parental_col]] == 0) %>%
      pull(Gene)

    if (length(off_on_genes) > 0) {
      cat("  Off->on genes (0 reads in parental, log2CPM > 2 in resistant):\n")
      cat("   ", paste(off_on_genes, collapse = ", "), "\n")
    }
  } else {
    cat("  ", counts_file, "not found, off->on gene detection skipped\n")
  }

  add_flags <- function(df) {
    df %>%
      dplyr::select(Gene, LFC_Res_vs_Par, Mean_Expr_Par, Mean_Expr_Res,
                    Classification) %>%
      mutate(
        In_met_core       = Gene %in% met_resistance_core,
        In_gimenez_common = Gene %in% gimenez_common,
        In_gimenez_egfr   = Gene %in% gimenez_egfr_bypass,
        In_gimenez_group2 = Gene %in% gimenez_group2,
        Forced            = Gene %in% FORCED_GENES,
        OffOn             = Gene %in% off_on_genes
      )
  }

  # Forced genes: add FORCED_GENES even if they do not pass the threshold
  forced_df <- lfc_table %>%
    filter(Gene %in% FORCED_GENES, !Gene %in% up_genes$Gene,
           !Gene %in% down_genes$Gene) %>%
    mutate(Classification = ifelse(LFC_Res_vs_Par >= 0, "Upregulated", "Downregulated"))

  if (nrow(forced_df) > 0) {
    cat("  Forced genes added (outside threshold):",
        paste(forced_df$Gene, collapse = ", "), "\n")
  }

  signature <- bind_rows(
    add_flags(up_genes),
    add_flags(down_genes),
    add_flags(forced_df)
  )

  sig_prefix <- if (DATASET == "MET1_METi") "MET" else RECEPTOR
  write.csv(signature,
            file.path(RESULTS_DIR, "tables",
                      paste0(sig_prefix, "_resistance_signature.csv")),
            row.names = FALSE)
  cat("  →", paste0(sig_prefix, "_resistance_signature.csv"), "(", nrow(signature), "genes: ",
      sum(signature$Classification == "Upregulated"), "up /",
      sum(signature$Classification == "Downregulated"), "down)\n")

  signature_cons <- signature %>% filter(Mean_Expr_Par > -5 | Forced | OffOn)
  write.csv(signature_cons,
            file.path(RESULTS_DIR, "tables",
                      paste0(sig_prefix, "_resistance_signature_conservative.csv")),
            row.names = FALSE)
  cat("  →", paste0(sig_prefix, "_resistance_signature_conservative.csv"), "(", nrow(signature_cons),
      "genes: ", sum(signature_cons$Classification == "Upregulated"), "up /",
      sum(signature_cons$Classification == "Downregulated"), "down)\n\n")
} else {
  # 5 classes — no known gene list, data-driven signature
  # (Acquired + Pre_existing + Downregulated; Maintained is left out of the
  # signature because it reflects drug response shared with the treated
  # parental, not something specific to resistance)
  signature <- bind_rows(up_genes, pre_genes, down_genes)

  write.csv(signature,
            file.path(RESULTS_DIR, "tables", paste0(sig_prefix, "_resistance_signature.csv")),
            row.names = FALSE)
  cat("  →", paste0(sig_prefix, "_resistance_signature.csv"), "(", nrow(signature), "genes: ",
      sum(signature$Classification == unname(CLASS_LABELS["acq"])), CLASS_LABELS["acq"], "/",
      sum(signature$Classification == unname(CLASS_LABELS["pre"])), CLASS_LABELS["pre"], "/",
      sum(signature$Classification == "Downregulated"), "down)\n")

  # Conservative: expressed above noise in the untreated parental (baseline)
  signature_cons <- signature %>% filter(Mean_Expr_Unt > -5)
  write.csv(signature_cons,
            file.path(RESULTS_DIR, "tables",
                      paste0(sig_prefix, "_resistance_signature_conservative.csv")),
            row.names = FALSE)
  cat("  →", paste0(sig_prefix, "_resistance_signature_conservative.csv"), "(", nrow(signature_cons),
      "genes: ", sum(signature_cons$Classification == unname(CLASS_LABELS["acq"])), CLASS_LABELS["acq"], "/",
      sum(signature_cons$Classification == unname(CLASS_LABELS["pre"])), CLASS_LABELS["pre"], "/",
      sum(signature_cons$Classification == "Downregulated"), "down)\n\n")
}

# FIGURE 3: SIGNATURE HEATMAP ----------------------------------------------------
if (DATASET %in% TWO_CLASS_DATASETS) {
  cat("5. Resistance signature heatmap (Figure 2)...\n")

  sig_genes <- signature_cons$Gene
  sig_genes <- sig_genes[sig_genes %in% rownames(mat) & is_clean_gene(sig_genes)]

  if (length(sig_genes) >= 3) {
    # Top 50: top 25 Up (LFC desc) + top 25 Down (LFC asc), FORCED_GENES included
    sig_up_all   <- signature_cons %>% filter(Classification == "Upregulated") %>%
      arrange(desc(LFC_Res_vs_Par)) %>% pull(Gene)
    sig_down_all <- signature_cons %>% filter(Classification == "Downregulated") %>%
      arrange(LFC_Res_vs_Par) %>% pull(Gene)

    # Forced genes count within the limit of 25 per group
    forced_up   <- intersect(FORCED_GENES, sig_up_all[sig_up_all %in% sig_genes])
    forced_down <- intersect(FORCED_GENES, sig_down_all[sig_down_all %in% sig_genes])

    top25_up   <- unique(c(forced_up,
                           head(sig_up_all[sig_up_all %in% sig_genes &
                                           !sig_up_all %in% forced_up],
                                25 - length(forced_up))))
    top25_down <- unique(c(forced_down,
                           head(sig_down_all[sig_down_all %in% sig_genes &
                                             !sig_down_all %in% forced_down],
                                25 - length(forced_down))))
    sig_ordered <- c(top25_up, top25_down)

    z_sig <- t(apply(mat[sig_ordered, ], 1, cal_zscore))

    ann_sig        <- if (DATASET == "EGFR1_EGFRi") {
      metadata_unit[, "Condition", drop = FALSE]
    } else annotation_col
    ann_sig_colors <- filter_annotation_colors(ann_sig)

    ph_sig <- pheatmap(z_sig,
             color             = col_fun_firma,
             annotation_col    = ann_sig,
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
                           length(top25_down), " down)"),
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
} else {
  cat("5. Resistance signature heatmap (Figure 2)...\n")

  sig_genes <- signature$Gene
  sig_genes <- sig_genes[sig_genes %in% rownames(mat) & is_clean_gene(sig_genes)]

  if (length(sig_genes) >= 3) {
    # Top 25 Acquired (by LFC_Res_vs_Trt) + top 25 Downregulated (by LFC_Res_vs_Unt)
    sig_up_all   <- signature %>% filter(Classification == unname(CLASS_LABELS["acq"])) %>%
      arrange(desc(LFC_Res_vs_Trt)) %>% pull(Gene)
    sig_down_all <- signature %>% filter(Classification == "Downregulated") %>%
      arrange(LFC_Res_vs_Unt) %>% pull(Gene)

    top25_up   <- head(sig_up_all[sig_up_all %in% sig_genes], 25)
    top25_down <- head(sig_down_all[sig_down_all %in% sig_genes], 25)
    sig_ordered <- c(top25_up, top25_down)

    z_sig <- t(apply(mat[sig_ordered, ], 1, cal_zscore))
    ann_sig_colors <- filter_annotation_colors(annotation_col)

    ph_sig <- pheatmap(z_sig,
             color             = heatmap_zscore,
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
                           length(top25_down), " down)"),
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
}

# AUTOMATIC GENE SELECTION FOR STRIPPLOT / LINEPLOT -----------------------
# Biologically fixed genes + top 3 Up and top 3 Down with the highest absolute
# LFC, cross-referenced with known-mechanism lists
if (DATASET %in% TWO_CLASS_DATASETS) {
  known_mech <- if (DATASET == "MET1_METi") {
    unique(c(met_resistance_core, met_bypass, gimenez_common,
             gimenez_egfr_bypass, gimenez_group2))
  } else {
    egfr_resistance_core
  }

  top3_up <- up_genes %>%
    filter(Gene %in% known_mech, Gene %in% rownames(mat)) %>%
    filter(!Gene %in% GENES_FIXED) %>%
    slice_head(n = 3) %>% pull(Gene)

  top3_down <- down_genes %>%
    filter(Gene %in% known_mech, Gene %in% rownames(mat)) %>%
    filter(!Gene %in% GENES_FIXED) %>%
    slice_head(n = 3) %>% pull(Gene)

  # If there are not enough in known lists, fill up with the highest absolute LFC
  if (length(top3_up) < 3) {
    extra <- up_genes %>%
      filter(Gene %in% rownames(mat), !Gene %in% GENES_FIXED,
             !Gene %in% top3_up) %>%
      slice_head(n = 3 - length(top3_up)) %>% pull(Gene)
    top3_up <- c(top3_up, extra)
  }
  if (length(top3_down) < 3) {
    extra <- down_genes %>%
      filter(Gene %in% rownames(mat), !Gene %in% GENES_FIXED,
             !Gene %in% top3_down) %>%
      slice_head(n = 3 - length(top3_down)) %>% pull(Gene)
    top3_down <- c(top3_down, extra)
  }

  vis_genes <- unique(c(
    GENES_FIXED[GENES_FIXED %in% rownames(mat)],
    top3_up, top3_down
  ))
  vis_genes <- head(vis_genes, 10)

  cat("  Stripplot/lineplot genes:\n")
  cat("    Fixed:   ", paste(GENES_FIXED[GENES_FIXED %in% rownames(mat)],
                             collapse = ", "), "\n")
  cat("    Top3 Up: ", paste(top3_up,   collapse = ", "), "\n")
  cat("    Top3 Down:", paste(top3_down, collapse = ", "), "\n\n")
} else {
  # 5-class scheme — no cross-reference with a known list in this selection,
  # top genes purely by LFC
  parental_unt <- rownames(metadata_unit)[
    metadata_unit$Condition == "Parental" & metadata_unit$Treatment == "No"]
  parental_trt <- rownames(metadata_unit)[
    metadata_unit$Condition == "Parental" & metadata_unit$Treatment == "Yes"]

  top5_up   <- up_genes   %>% filter(Gene %in% rownames(mat)) %>% slice_head(n = 5) %>% pull(Gene)
  top5_down <- down_genes %>% filter(Gene %in% rownames(mat)) %>% slice_head(n = 5) %>% pull(Gene)
  vis_genes <- unique(c(top5_up, top5_down))

  cat("  Stripplot/lineplot genes (top LFC, no known list):\n")
  cat("    Top5 Acquired:     ", paste(top5_up,   collapse = ", "), "\n")
  cat("    Top5 Downregulated:", paste(top5_down, collapse = ", "), "\n\n")
}

# FIGURE 4: STRIPPLOT --------------------------------------------------------
if (DATASET %in% TWO_CLASS_DATASETS) {
  cat("6. Key genes stripplot (Figure 3)...\n")

  group_var    <- if (DATASET == "MET1_METi") metadata_unit$GX_cluster else metadata_unit$Cell_Line
  group_colors <- if (DATASET == "MET1_METi") annotation_colors$GX_cluster else annotation_colors$Cell_Line

  if (DATASET == "EGFR1_EGFRi") {
    # EGFR1_EGFRi: 2 samples (EGFR1 parental, EGFR1_R resistant) — Sample on x,
    # color by Condition.
    sample_order <- c(parental, resistant)
    strip_df <- lapply(vis_genes, function(g) {
      data.frame(
        Gene      = g,
        Sample    = factor(colnames(mat), levels = sample_order),
        Expr      = as.numeric(mat[g, ]),
        Cell_Line = metadata_unit$Cell_Line,
        Condition = metadata_unit$Condition
      )
    }) %>% bind_rows() %>%
      mutate(Gene = factor(Gene, levels = vis_genes))

    p_strip <- ggplot(strip_df, aes(x = Sample, y = Expr,
                                     color = Condition, shape = Cell_Line)) +
      geom_point(size = 3, alpha = 0.9) +
      scale_color_manual(values = condition_colors) +
      scale_shape_manual(values = c(EGFR1 = 17)) +
      facet_wrap(~ Gene, scales = "free_y", ncol = 5) +
      labs(
        title    = paste("Expression per sample —", OUT_PREFIX),
        subtitle = paste("Fixed:", paste(GENES_FIXED, collapse = ", "),
                         "| Top3 Up:", paste(top3_up, collapse = ", "),
                         "| Top3 Down:", paste(top3_down, collapse = ", ")),
        x = NULL, y = "log2CPM (TMM)"
      )
  } else {
    strip_df <- lapply(vis_genes, function(g) {
      data.frame(
        Gene      = g,
        Sample    = colnames(mat),
        Expr      = as.numeric(mat[g, ]),
        Group     = group_var,
        Condition = metadata_unit$Condition
      )
    }) %>% bind_rows() %>%
      mutate(Gene = factor(Gene, levels = vis_genes))

    p_strip <- ggplot(strip_df, aes(x = Group, y = Expr, color = Group)) +
      geom_jitter(width = 0.15, size = 2.5, alpha = 0.9) +
      stat_summary(fun = mean, geom = "crossbar", width = 0.4,
                   linewidth = 0.5, color = "black") +
      scale_color_manual(values = group_colors) +
      facet_wrap(~ Gene, scales = "free_y", ncol = 5) +
      labs(
        title    = paste("Expression by GX cluster —", OUT_PREFIX),
        subtitle = paste("Fixed:", paste(GENES_FIXED, collapse = ", "),
                         "| Top3 Up:", paste(top3_up, collapse = ", "),
                         "| Top3 Down:", paste(top3_down, collapse = ", ")),
        x = NULL, y = "log2CPM (TMM)"
      )
  }

  p_strip <- p_strip +
    theme_minimal(base_size = 10) +
    theme(
      legend.position  = if (DATASET == "EGFR1_EGFRi") "right" else "none",
      axis.text.x      = element_text(angle = 30, hjust = 1, size = 7),
      strip.text       = element_text(face = "bold", size = 8),
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
} else if (length(vis_genes) > 0) {
  cat("6. Key genes stripplot (Figure 3)...\n")

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
      subtitle = paste0("Top ", CLASS_LABELS["acq"], ": ", paste(top5_up, collapse = ", "),
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

# FIGURE 5: LINEPLOT ---------------------------------------------------------
if (DATASET %in% TWO_CLASS_DATASETS) {
  cat("7. Key genes lineplot (Figure 4)...\n")

  if (DATASET == "EGFR1_EGFRi") {
    # Paired lineplot: Parental → Resistance, one line per pair, color by Condition
    pair_map <- setNames(metadata_unit$Cell_Line, rownames(metadata_unit))
    line_df <- lapply(vis_genes, function(g) {
      do.call(rbind, lapply(rownames(metadata_unit), function(s) {
        data.frame(
          Gene      = g,
          Sample    = s,
          Expr      = mat[g, s],
          Cell_Line = pair_map[s],
          Condition = metadata_unit[s, "Condition"]
        )
      }))
    }) %>% bind_rows() %>%
      mutate(
        Gene      = factor(Gene, levels = vis_genes),
        Condition = factor(Condition, levels = c("Parental", "Resistance"))
      )

    cellline_pal <- c(EGFR1 = "#444444")

    p_line <- ggplot(line_df, aes(x = Condition, y = Expr, group = Cell_Line)) +
      geom_line(aes(color = Cell_Line, linetype = Cell_Line),
                linewidth = 0.8, alpha = 0.7) +
      geom_point(aes(fill = Condition, shape = Cell_Line),
                 size = 3.5, alpha = 0.95, color = "white", stroke = 0.3) +
      scale_color_manual(values = cellline_pal, guide = "none") +
      scale_fill_manual(values = condition_colors) +
      scale_shape_manual(values = c(EGFR1 = 24)) +
      scale_linetype_manual(values = c(EGFR1 = "dashed")) +
      facet_wrap(~ Gene, scales = "free_y", ncol = 5) +
      labs(
        title    = paste("Key gene expression per model —", OUT_PREFIX),
        subtitle = paste("Each line = one patient pair —", OUT_PREFIX),
        x = NULL, y = "log2CPM (TMM)",
        fill = "Condition", shape = "Cell Line", linetype = "Cell Line"
      )
  } else {
    parental_expr <- setNames(
      as.numeric(mat[vis_genes, parental[1], drop = FALSE]),
      vis_genes
    )
    line_df <- lapply(vis_genes, function(g) {
      data.frame(
        Gene          = g,
        Sample        = colnames(mat),
        Expr          = as.numeric(mat[g, ]),
        Group         = group_var,
        Condition     = metadata_unit$Condition,
        Parental_Expr = parental_expr[g]
      )
    }) %>% bind_rows() %>%
      filter(Condition == "Resistance") %>%
      mutate(
        Gene   = factor(Gene, levels = vis_genes),
        Sample = factor(Sample, levels = resistant)
      )

    p_line <- ggplot(line_df, aes(x = Sample, y = Expr, color = Group,
                                   group = Gene)) +
      geom_hline(aes(yintercept = Parental_Expr), linetype = "dashed",
                 color = "gray60", linewidth = 0.5) +
      geom_line(linewidth = 0.6, alpha = 0.5) +
      geom_point(size = 3, alpha = 0.9) +
      scale_color_manual(values = group_colors) +
      facet_wrap(~ Gene, scales = "free_y", ncol = 5) +
      labs(
        title    = paste("Key gene expression per model —", OUT_PREFIX),
        subtitle = "Dashed line = Parental level",
        x = NULL, y = "log2CPM (TMM)", color = "GX cluster"
      )
  }

  p_line <- p_line +
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
} else if (length(vis_genes) > 0) {
  cat("7. Key genes lineplot (Figure 4)...\n")

  # Trajectory per resistant model — dashed line = Treated level (the reference
  # against which the classification is defined). With several treated samples the
  # mean is used; with a single one, its value directly.
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
cat("8. Generating report...\n")

if (DATASET %in% TWO_CLASS_DATASETS) {
  known_hits_up   <- intersect(up_genes$Gene,   known_mech)
  known_hits_down <- intersect(down_genes$Gene, known_mech)

  report_text <- paste0(
    "# Key Genes Report — ", OUT_PREFIX, "\n\n",
    "## Top 15 Upregulated\n",
    paste(capture.output(print(
      as.data.frame(slice_head(up_genes, n = 15) %>%
                      dplyr::select(Gene, LFC_Res_vs_Par)),
      row.names = FALSE)), collapse = "\n"), "\n\n",
    "## Top 15 Downregulated\n",
    paste(capture.output(print(
      as.data.frame(slice_head(down_genes, n = 15) %>%
                      dplyr::select(Gene, LFC_Res_vs_Par)),
      row.names = FALSE)), collapse = "\n"), "\n\n",
    "## Forced genes (NF2/NRG1)\n",
    if (nrow(forced_df) > 0)
      paste(forced_df$Gene, "LFC =",
            round(forced_df$LFC_Res_vs_Par, 2), collapse = "; ")
    else "All pass the threshold", "\n\n",
    "## Known genes detected (Upregulated)\n",
    if (length(known_hits_up) > 0)
      paste(known_hits_up, collapse = ", ") else "none", "\n\n",
    "## Known genes detected (Downregulated)\n",
    if (length(known_hits_down) > 0)
      paste(known_hits_down, collapse = ", ") else "none", "\n\n",
    "## Genes stripplot/lineplot\n",
    "- Fixed:    ", paste(GENES_FIXED, collapse = ", "), "\n",
    "- Top3 Up:  ", paste(top3_up,    collapse = ", "), "\n",
    "- Top3 Down:", paste(top3_down,  collapse = ", "), "\n\n",
    paste0("## Signature ", sig_prefix, "_resistance_signature_conservative.csv\n"),
    "- N genes: ", nrow(signature_cons), " (",
    sum(signature_cons$Classification == "Upregulated"), " up / ",
    sum(signature_cons$Classification == "Downregulated"), " down)\n",
    "- Forced genes included: ",
    paste(FORCED_GENES[FORCED_GENES %in% signature_cons$Gene],
          collapse = ", "), "\n\n",
    "## Figures\n",
    "- 03_heatmap_gimenez.png     — exact Giménez-Xavier 2017 genes\n",
    "- 03_heatmap_signature.png   — full signature with forced NF2/NRG1\n",
    "- 03_stripplot_key_genes.png — expression per GX cluster\n",
    "- 03_lineplot_key_genes.png  — trajectory per model vs parental\n"
  )
} else {
  acq_lbl <- unname(CLASS_LABELS["acq"])
  pre_lbl <- unname(CLASS_LABELS["pre"])

  report_text <- paste0(
    "# Key Genes Report — ", OUT_PREFIX, "\n\n",
    "**Signature:** data-driven (no known gene list) — ", acq_lbl, " + ", pre_lbl, " + Downregulated\n\n",
    "## Top 15 ", acq_lbl, "\n",
    paste(capture.output(print(
      as.data.frame(slice_head(up_genes, n = 15) %>%
                      dplyr::select(Gene, LFC_Res_vs_Trt)),
      row.names = FALSE)), collapse = "\n"), "\n\n",
    "## Top 10 ", pre_lbl, "\n",
    paste(capture.output(print(
      as.data.frame(slice_head(pre_genes, n = 10) %>%
                      dplyr::select(Gene, LFC_Res_vs_Trt, LFC_Res_vs_Unt)),
      row.names = FALSE)), collapse = "\n"), "\n\n",
    "## alk_resistance_core detected (generic reference, not specific to ", RECEPTOR, ")\n",
    "- ", acq_lbl, ": ", paste(core_acq,  collapse = ", "), "\n",
    "- ", pre_lbl, ": ", paste(core_pre,  collapse = ", "), "\n",
    "- Downregulated:", paste(core_down, collapse = ", "), "\n\n",
    paste0("## Signature ", sig_prefix, "_resistance_signature_conservative.csv\n"),
    "- N genes: ", nrow(signature_cons), " (",
    sum(signature_cons$Classification == acq_lbl), " ", acq_lbl, " / ",
    sum(signature_cons$Classification == pre_lbl), " ", pre_lbl, " / ",
    sum(signature_cons$Classification == "Downregulated"), " down)\n\n",
    "## Stripplot/lineplot genes (top LFC, no known list)\n",
    "- Top5 ", acq_lbl, ":     ", paste(top5_up,   collapse = ", "), "\n",
    "- Top5 Downregulated:", paste(top5_down, collapse = ", "), "\n\n",
    "## Figures\n",
    "- 03_heatmap_signature.png   — full signature (top25 up + top25 down)\n",
    "- 03_stripplot_key_genes.png — expression per sample\n",
    "- 03_lineplot_key_genes.png  — trajectory per resistant model vs Treated level\n"
  )
}

writeLines(report_text, file.path(RESULTS_DIR, "report", "03_key_genes_report.md"))
cat("  Report: report/03_key_genes_report.md\n\n")

HAVE_REPORT <- requireNamespace("rmarkdown", quietly = TRUE) && rmarkdown::pandoc_available()
if (!HAVE_REPORT) cat("rmarkdown or pandoc not available: skipping the HTML report
")
if (HAVE_REPORT) {
  cat("Generating HTML report...\n")
  rmarkdown::render(here::here("scripts", "simple", "rna_03_key_genes_report.Rmd"),
                    output_dir  = file.path(RESULTS_DIR, "report"),
                    output_file = "03_key_genes_report.html",
                    quiet       = TRUE)
  cat("  → report/03_key_genes_report.html\n\n")
}

cat("================================================================\n")
cat("  KEY GENES completed —", OUT_PREFIX, "\n")
cat("================================================================\n\n")
