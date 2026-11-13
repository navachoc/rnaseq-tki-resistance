# =============================================================================
# MET_01_summary_figures.R
#
# Purpose: MET1 (MET) expression figures: PCA, MA plot and GO BP / KEGG
# functional enrichment of the resistance signature.
#
#   Sections 3-4. PCA of parental MET1 and its six resistant models (C1, V8,
#     PR1-MET1R). Reproduces the PCA of scripts/simple/rna_01.1_qc.R with
#     DATASET = "MET1_METi" (same data, MIN_COUNT filter, TMM normalisation
#     and top 500 most variable genes), in two versions: with model names next
#     to each point and without names. Variance with one decimal on the axes.
#     Colour by Condition and shape by Treatment, as in rna_01.1_qc.R.
#   Section 5. MA plot resistant vs parental from the rna_02.1_lfc.R table.
#     Same thresholds and colours as 02_lfc_scatter.png, with a fixed list of
#     genes labelled in black.
#   Section 6. Functional enrichment of the signature, one GO BP plot and one
#     KEGG plot, with upregulated and repressed genes (same style as RET_09
#     and ALK_09c). Legends in English.
#   PCA and MA plot without legend or title.
#
# Inputs:
#   data/counts_MET1_METi.csv, data/metadata_MET1_METi.csv
#   config/gene_lists.R
#   results/MET/MET1/tables/02_lfc_master_table.csv
#   results/MET/MET1/tables/04_go_bp_{up,down}regulated.csv
#   results/MET/MET1/tables/04_kegg_{up,down}regulated.csv
#
# Outputs (results/MET/MET1/):
#   figures/MET_01_fig14_pca_with_names.png, MET_01_fig14_pca_no_names.png
#   figures/MET_01_fig14_ma_plot.png
#   figures/MET_01_fig15b_{go,kegg}.png, tables/MET_01_fig15b_{go,kegg}.csv
#
# Dependencies (run first):
#   scripts/simple/rna_02.1_lfc.R with DATASET = "MET1_METi" (02_* table)
#   scripts/utils/met_egfr_offon_regenerate.R (04_* tables)
#
# Usage: paths are resolved with here::here() from the repository root
# (.here file), so it works from any working directory:
#   Rscript scripts/specific/MET/MET_01_summary_figures.R
# =============================================================================

suppressPackageStartupMessages({
  library(edgeR)
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
})

DATASET   <- "MET1_METi"
MIN_COUNT <- 10
TOP_N_PCA <- 500

DATA_DIR <- here::here("data")
OUT_DIR  <- here::here("results", "MET", "MET1", "figures")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Same values as condition_colors in config/colors.R. That file is not
# source()d because it depends on CELL_LINE (same approach as RET_08).
condition_colors <- c(Parental = "#e6a817", Resistance = "#013c58")

detect_sep <- function(file) {
  if (grepl(";", readLines(file, n = 1))) ";" else ","
}

# 1. Data ---------------------------------------------------------------------
counts_file <- file.path(DATA_DIR, paste0("counts_", DATASET, ".csv"))
meta_file   <- file.path(DATA_DIR, paste0("metadata_", DATASET, ".csv"))

counts   <- read.csv(counts_file, sep = detect_sep(counts_file), check.names = FALSE)
metadata <- read.csv(meta_file,   sep = detect_sep(meta_file))

rownames(counts) <- make.names(counts[, 2], unique = TRUE)
counts <- counts[, -c(1:3)]

sample_to_model <- setNames(metadata$Model, metadata$Sample)
colnames(counts) <- sample_to_model[colnames(counts)]
rownames(metadata) <- metadata$Model
counts <- counts[, rownames(metadata)]

# 2. TMM normalisation ---------------------------------------------------------
keep <- rowSums(counts) >= MIN_COUNT
dge  <- DGEList(counts = counts[keep, ])
dge  <- calcNormFactors(dge, method = "TMM")
mat  <- cpm(dge, log = TRUE, prior.count = 1)
cat("Genes after filter (MIN_COUNT =", MIN_COUNT, "):", nrow(mat), "\n")

# 3. PCA -------------------------------------------------------------------------
gene_var  <- apply(mat, 1, var)
n_pca     <- min(TOP_N_PCA, nrow(mat))
top_genes <- names(sort(gene_var, decreasing = TRUE))[seq_len(n_pca)]
pca       <- prcomp(t(mat[top_genes, ]), scale. = TRUE)
var_pc1   <- round(summary(pca)$importance[2, 1] * 100, 1)
var_pc2   <- round(summary(pca)$importance[2, 2] * 100, 1)
cat("PCA on", n_pca, "genes | PC1:", var_pc1, "% | PC2:", var_pc2, "%\n")

pca_df <- data.frame(pca$x, metadata)
pca_df$Sample <- rownames(pca_df)

# 4. Figure ----------------------------------------------------------------------
p_base <- ggplot(pca_df, aes(PC1, PC2, color = Condition, shape = Treatment)) +
  geom_point(size = 4, alpha = 0.85) +
  scale_color_manual(values = condition_colors) +
  scale_shape_manual(values = c(Yes = 17, No = 16)) +
  labs(x = paste0("PC1 (", var_pc1, "%)"), y = paste0("PC2 (", var_pc2, "%)")) +
  theme_minimal(base_size = 12) +
  theme(legend.position  = "none",
        panel.grid       = element_line(color = "gray90", linetype = "dashed"),
        panel.border     = element_rect(color = "black", fill = NA),
        panel.background = element_rect(fill = "white"),
        plot.background  = element_rect(fill = "white"))

p_con_nombres <- p_base +
  geom_text_repel(aes(label = Sample), size = 3, color = "black",
                  max.overlaps = 20, box.padding = 0.3, segment.size = 0.2, seed = 1)

ggsave(file.path(OUT_DIR, "MET_01_fig14_pca_with_names.png"),
       p_con_nombres, width = 10, height = 7, dpi = 300)
ggsave(file.path(OUT_DIR, "MET_01_fig14_pca_no_names.png"),
       p_base, width = 10, height = 7, dpi = 300)

cat("Saved in", OUT_DIR, ":\n  MET_01_fig14_pca_with_names.png\n  MET_01_fig14_pca_no_names.png\n")

# 5. MA plot -------------------------------------------------------------------
LOGFC_THRESHOLD <- 2.0   # same threshold as rna_02.1_lfc.R
TOP_N_LABEL     <- 15    # genes labelled per direction, same criterion

# Same values as classification_colors in config/colors.R (MET classes)
classification_colors <- c(Upregulated = "#e67e22", Downregulated = "#3498db",
                           Other = "#ecf0f1")

source(here::here("config", "gene_lists.R"))  # MET and Gimenez key genes

lfc_file <- here::here("results", "MET", "MET1", "tables", "02_lfc_master_table.csv")
if (!file.exists(lfc_file)) {
  stop("Missing ", lfc_file, " -- first run scripts/simple/rna_02.1_lfc.R with DATASET = \"MET1_METi\".")
}
lfc_table <- read.csv(lfc_file, stringsAsFactors = FALSE)

key_genes <- unique(c(met_resistance_core, gimenez_common, gimenez_egfr_bypass, gimenez_group2))
top_up    <- lfc_table %>% filter(Classification == "Upregulated") %>%
  arrange(desc(LFC_Res_vs_Par)) %>% slice_head(n = TOP_N_LABEL) %>% pull(Gene)
top_down  <- lfc_table %>% filter(Classification == "Downregulated") %>%
  arrange(LFC_Res_vs_Par) %>% slice_head(n = TOP_N_LABEL) %>% pull(Gene)
genes_to_label <- unique(c(
  intersect(key_genes, lfc_table$Gene[lfc_table$Classification != "Other"]),
  top_up, top_down
))
# Fixed label list, as in the EGFR MA plot (egfr_summary_panels.R), so that a
# large font can be used. With a finite max.overlaps ggrepel drops different
# labels depending on spacing; fixing the list and max.overlaps = Inf always
# yields the same labels.
GENES_ETIQUETA <- c("ABI3BP", "C9orf84", "DDR1", "ADORA2A", "NNMT", "SAMD11", "TLN2",
                    "CELF2", "ANKRD1", "RSPH10B2", "SARDH", "AGR2", "C5orf63",
                    "SLC24A5", "RSPH10B", "TTN", "CAPN12",
                    "SPP1", "FAM151A", "SLCO1B3")
genes_to_label <- intersect(genes_to_label, GENES_ETIQUETA)
cat("MA plot:", sum(lfc_table$Classification == "Upregulated"), "upregulated,",
    sum(lfc_table$Classification == "Downregulated"), "downregulated,",
    length(genes_to_label), "labelled genes\n")

ma_df <- lfc_table %>%
  mutate(
    Classification = factor(Classification, levels = c("Upregulated", "Downregulated", "Other")),
    label = ifelse(Gene %in% genes_to_label & Classification != "Other", Gene, "")
  ) %>%
  arrange(desc(Classification))

p_ma <- ggplot(ma_df, aes(Mean_Expr_Par, LFC_Res_vs_Par)) +
  geom_point(aes(color = Classification), alpha = 0.4, size = 0.7) +
  geom_hline(yintercept = c(-LOGFC_THRESHOLD, LOGFC_THRESHOLD),
             linetype = "dashed", color = "gray50", linewidth = 0.5) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.3) +
  # Size 6 and more spacing between labels so they are legible.
  geom_text_repel(aes(label = label), size = 6, fontface = "bold", color = "black",
                  max.overlaps = Inf, segment.size = 0.2, box.padding = 0.5,
                  force = 2, max.iter = 20000, seed = 42) +
  scale_color_manual(values = classification_colors) +
  labs(x = "Mean Parental log2CPM", y = "LFC Resistance vs Parental") +
  theme_minimal(base_size = 12) +
  theme(legend.position  = "none",
        panel.background = element_rect(fill = "white"),
        plot.background  = element_rect(fill = "white"))

ggsave(file.path(OUT_DIR, "MET_01_fig14_ma_plot.png"), p_ma, width = 9, height = 7, dpi = 300)
cat("Saved in", OUT_DIR, ":\n  MET_01_fig14_ma_plot.png\n")

# 6. Functional enrichment, one GO BP plot and one KEGG plot --------------------
# Same style as the RET_09 and ALK_09c enrichment: upward triangle =
# upregulated, downward triangle = repressed, size = number of genes,
# colour = adjusted p, X axis = fold enrichment. Reads the ORA tables already
# generated by scripts/utils/met_egfr_offon_regenerate.R on the
# extended signature (04_go_bp_*regulated.csv / 04_kegg_*regulated.csv),
# without rerunning clusterProfiler.
#
# Pathway selection: in each plot and direction, up to N_VIAS pathways with
# p.adjust < 0.05 are included. Mandatory ones first (VIAS_FIJAS), then
# filled by increasing p.adjust, skipping those that share more than
# JACCARD_MAX of their genes with an already chosen one (avoids repeating the
# same signal under different IDs, e.g. "extracellular matrix /
# extracellular structure / external encapsulating structure organization",
# which have the same 31 genes).
N_VIAS      <- 5
JACCARD_MAX <- 0.5
P_ADJ_MAX   <- 0.05
MIN_GENES   <- 5     # pathways with fewer genes inflate the fold enrichment and
                     # stretch the X axis (e.g. "preganglionic parasympathetic
                     # fiber development", 4 genes, FE 27.9)
VIAS_FIJAS  <- list(
  GO_BP = list(Up = c("GO:0016055",   # Wnt signaling pathway
                      "GO:0030198",   # extracellular matrix organization
                      "GO:0070371"),  # ERK1 and ERK2 cascade
               Down = character(0)),
  KEGG  = list(Up = c("hsa04015",     # Rap1 signaling pathway
                      "hsa04390"),    # Hippo signaling pathway
               Down = character(0))
)

# Explicitly excluded pathways: very high fold enrichment with few genes,
# they separated from the rest on the X axis. They are removed without
# adding others in their place.
VIAS_EXCLUIDAS <- c("GO:0048483",   # autonomic nervous system development (5 genes, FE 11.6)
                    "GO:0001570")   # vasculogenesis (7 genes, FE 9.4)

COLOR_SIG   <- "#D6604D"  # most significant, same hex values as RET_09 and ALK
COLOR_NOSIG <- "#4393C3"  # least significant

tablas_dir <- here::here("results", "MET", "MET1", "tables")

# excluir (exclude): IDs that cannot be included (the repressed set excludes
# the pathways already chosen as upregulated, so that a pathway does not
# appear in both directions).
seleccionar_vias <- function(fuente, direccion, excluir = character(0)) {
  archivo <- paste0(if (fuente == "GO_BP") "04_go_bp_" else "04_kegg_",
                    tolower(direccion), "regulated.csv")
  ruta <- file.path(tablas_dir, archivo)
  if (!file.exists(ruta)) {
    stop("Missing ", ruta, " -- first run scripts/utils/met_egfr_offon_regenerate.R.")
  }
  tabla <- read.csv(ruta, stringsAsFactors = FALSE) %>%
    filter(p.adjust < P_ADJ_MAX, Count >= MIN_GENES, !ID %in% excluir) %>%
    arrange(p.adjust)
  genes <- setNames(strsplit(tabla$geneID, "/"), tabla$ID)

  fijas <- VIAS_FIJAS[[fuente]][[direccion]]
  faltan <- setdiff(fijas, tabla$ID)
  if (length(faltan) > 0) {
    stop("Mandatory pathway without p.adjust < ", P_ADJ_MAX, " in ", archivo, ": ", paste(faltan, collapse = ", "))
  }
  elegidas <- fijas
  for (id in tabla$ID) {
    if (length(elegidas) >= N_VIAS) break
    if (id %in% elegidas) next
    jaccard <- sapply(elegidas, function(e) {
      length(intersect(genes[[id]], genes[[e]])) / length(union(genes[[id]], genes[[e]]))
    })
    if (length(jaccard) > 0 && any(jaccard > JACCARD_MAX)) next
    elegidas <- c(elegidas, id)
  }

  # VIAS_EXCLUIDAS are removed after selection, without filling their gap
  tabla %>%
    filter(ID %in% elegidas, !ID %in% VIAS_EXCLUIDAS) %>%
    select(ID, Description, FoldEnrichment, p.adjust, Count) %>%
    mutate(Direccion = direccion, Fuente = fuente)
}

dibujar_enriquecimiento <- function(fuente) {
  vias_up   <- seleccionar_vias(fuente, "Up")
  vias_down <- seleccionar_vias(fuente, "Down", excluir = vias_up$ID)
  df <- bind_rows(vias_up, vias_down) %>%
    mutate(Direccion = factor(Direccion, levels = c("Up", "Down"))) %>%
    arrange(Direccion, desc(FoldEnrichment)) %>%
    mutate(Description = factor(Description, levels = rev(unique(Description))))

  cat("\n", fuente, "-- selected pathways:\n")
  print(df %>% mutate(p.adjust = signif(p.adjust, 3), FoldEnrichment = round(FoldEnrichment, 2)) %>%
          select(Direccion, ID, Description, Count, p.adjust, FoldEnrichment), row.names = FALSE)

  p <- ggplot(df, aes(x = FoldEnrichment, y = Description,
                      size = Count, fill = p.adjust, shape = Direccion)) +
    geom_point(color = "grey20", stroke = 0.3) +
    scale_shape_manual(values = c(Up = 24, Down = 25), name = "Direction",
                       labels = c(Up = "Upregulated", Down = "Downregulated"), drop = FALSE) +
    scale_fill_gradient(low = COLOR_SIG, high = COLOR_NOSIG, trans = "log10",
                        name = "adjusted p-value", guide = guide_colorbar(reverse = TRUE)) +
    scale_size_continuous(range = c(3, 10), name = "Genes",
                          guide = guide_legend(override.aes = list(shape = 24, fill = "grey40"))) +
    labs(x = "Fold enrichment", y = NULL) +
    theme_minimal(base_size = 12) +
    theme(panel.border     = element_rect(color = "black", fill = NA),
          panel.background = element_rect(fill = "white"),
          plot.background  = element_rect(fill = "white"))

  sufijo <- if (fuente == "GO_BP") "go" else "kegg"
  write.csv(df, here::here("results", "MET", "MET1", "tables", paste0("MET_01_fig15b_", sufijo, ".csv")),
            row.names = FALSE)
  ggsave(file.path(OUT_DIR, paste0("MET_01_fig15b_", sufijo, ".png")), p, width = 8, height = 5, dpi = 300)
  cat("Saved in", OUT_DIR, ":\n  MET_01_fig15b_", sufijo, ".png\n", sep = "")
}

dibujar_enriquecimiento("GO_BP")
dibujar_enriquecimiento("KEGG")
