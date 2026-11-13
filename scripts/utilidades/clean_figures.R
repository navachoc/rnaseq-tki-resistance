# =============================================================================
# clean_figures.R
# Purpose: minimal, anonymized versions of the summary figures, for showcasing
#   the pipeline output without sample or cell line names.
#
#   Per dataset:
#     <ID>_pca.png                PCA of the top 500 most variable genes
#     <ID>_signature_heatmap.png  heatmap of the genes in the conservative
#                                 signature, split into up and down (no gene
#                                 names)
#     <ID>_ma.png                 MA plot of the fold changes with the
#                                 |log2FC| threshold (no gene labels, no counts)
#   Datasets with a single parental and a single resistant sample (EGFR) only
#   get the MA plot. Datasets without an LFC master table (RET) get no MA plot.
#
#   Samples are renamed "<RECEPTOR>-<n> <R|P>": receptor, a running number
#   (parental samples first) and R (resistant) or P (parental). Datasets are
#   named "<RECEPTOR>_<k>". The key between anonymous labels and real sample
#   names is written to results/clean/_label_key.csv and must not be published.
#
# Method: same as the per-target figure scripts. For the PCA, genes starting
#   with "ENSG" or containing "." are dropped, the 500 most variable genes are
#   taken from the TMM matrix, z-scored per gene, and the PCA is run on that
#   matrix. The signature heatmap shows every gene of the conservative
#   signature, z-scored per gene across samples, with the genes of each
#   direction clustered within their own block. A gene listed for more than one
#   drug is drawn once; genes with opposite directions between drugs are left
#   out.
#
# Inputs (written by the pipeline, for each dataset):
#   results/<RECEPTOR>/<MODEL>/tables/01_mat_tmm.csv
#   results/<RECEPTOR>/<MODEL>/tables/01_metadata_unit.csv
#   results/<RECEPTOR>/<MODEL>/tables/02_lfc_master_table.csv   (MA plot only)
#   results/<RECEPTOR>/<MODEL>/tables/<RECEPTOR>_resistance_signature_conservative.csv
# Outputs: results/clean/
#
# Run: source(here::here("scripts/utilidades/clean_figures.R"))
#   or Rscript from any folder inside the repository.
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

OUT_DIR <- here::here("results", "clean")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

COLOR_CONDITION <- c(Parental = "#e6a817", Resistant = "#013c58")
GRADIENT_RNA    <- c("#000000", "#808080", "#f2a900")
COLOR_UP        <- "#e67e22"
COLOR_DOWN      <- "#2980b9"
COLOR_NS        <- "#b2bec3"

TOP_N           <- 500   # most variable genes for the PCA
Z_MAX           <- 2     # z-score scale of the heatmap, -Z_MAX to +Z_MAX
LOGFC_THRESHOLD <- 2     # same threshold as the pipeline

# One entry per dataset. figs: which figures to draw. lfc and base: columns of
# 02_lfc_master_table.csv used for the MA plot (log2 fold change, and the
# expression of the reference group).
DATASETS <- list(
  list(receptor = "MET",  model = "MET1",  index = 1, figs = c("pca", "signature", "ma"),
       lfc = "LFC_Res_vs_Par", base = "Mean_Expr_Par"),
  list(receptor = "EGFR", model = "EGFR1",  index = 1, figs = "ma",
       lfc = "LFC_Res_vs_Par", base = "Mean_Expr_Par"),
  list(receptor = "RET",  model = "RET1", index = 1, figs = c("pca", "signature"),
       lfc = NA, base = NA),
  list(receptor = "RET",  model = "RET2",  index = 2, figs = c("pca", "signature"),
       lfc = NA, base = NA),
  list(receptor = "ALK",  model = "ALK1", index = 1, figs = c("pca", "signature", "ma"),
       lfc = "LFC_Res_vs_Unt", base = "Mean_Expr_Unt"),
  list(receptor = "ALK",  model = "ALK2", index = 2, figs = c("pca", "signature", "ma"),
       lfc = "LFC_Res_vs_Unt", base = "Mean_Expr_Unt"),
  list(receptor = "ALK",  model = "ALK3",  index = 3, figs = c("pca", "signature", "ma"),
       lfc = "LFC_Res_vs_Unt", base = "Mean_Expr_Unt")
)

dataset_id <- function(ds) paste0(ds$receptor, "_", ds$index)
tables_dir <- function(ds) here::here("results", ds$receptor, ds$model, "tables")

theme_clean_panel <- function() {
  theme_classic(base_size = 11) +
    theme(panel.border = element_rect(fill = NA, color = "black", linewidth = 0.8),
          axis.line = element_blank())
}

# =============================================================================
# 1. Read and anonymize
# =============================================================================
read_tmm <- function(ds) {
  df <- read.csv(file.path(tables_dir(ds), "01_mat_tmm.csv"), check.names = FALSE)
  m <- as.matrix(df[, -1, drop = FALSE])
  rownames(m) <- df$Gene
  m
}

# The first column of the sample sheet holds the names used as matrix columns.
read_samples <- function(ds) {
  meta <- read.csv(file.path(tables_dir(ds), "01_metadata_unit.csv"),
                   check.names = FALSE, stringsAsFactors = FALSE)
  data.frame(sample = meta[[1]],
             condition = ifelse(meta$Condition == "Resistance", "Resistant", "Parental"),
             stringsAsFactors = FALSE)
}

anonymize <- function(receptor, condition) {
  ord <- order(factor(condition, levels = names(COLOR_CONDITION)), seq_along(condition))
  lab <- character(length(condition))
  lab[ord] <- paste0(receptor, "-", seq_along(ord), " ",
                     ifelse(condition[ord] == "Resistant", "R", "P"))
  lab
}

load_dataset <- function(ds) {
  mat <- read_tmm(ds)
  smp <- read_samples(ds)
  smp <- smp[match(colnames(mat), smp$sample), ]
  stopifnot(!anyNA(smp$sample))
  smp$label <- anonymize(ds$receptor, smp$condition)
  key <- data.frame(dataset = dataset_id(ds), label = smp$label, sample = smp$sample,
                    condition = smp$condition, stringsAsFactors = FALSE)
  colnames(mat) <- smp$label
  list(mat = mat, smp = smp, key = key)
}

top_z <- function(mat, n = TOP_N) {
  g <- rownames(mat)
  g <- g[!grepl("^ENSG", g) & !grepl("\\.", g)]
  m <- mat[g, , drop = FALSE]
  v <- apply(m, 1, var, na.rm = TRUE)
  top <- names(sort(v, decreasing = TRUE))[seq_len(min(n, length(v)))]
  raw <- m[top, , drop = FALSE]
  list(raw = raw, z = t(scale(t(raw))))
}

# =============================================================================
# 2. PCA
# =============================================================================
plot_pca <- function(d) {
  top <- top_z(d$mat)
  pca <- prcomp(t(top$z), center = TRUE, scale. = FALSE)
  ve  <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
  df  <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], label = rownames(pca$x),
                    stringsAsFactors = FALSE)
  df$condition <- d$smp$condition[match(df$label, d$smp$label)]

  ggplot(df, aes(PC1, PC2, color = condition, label = label)) +
    geom_point(size = 3.5) +
    geom_text_repel(size = 3, seed = 1, max.overlaps = Inf, show.legend = FALSE,
                    box.padding = 0.6, point.padding = 0.5, force = 2) +
    scale_color_manual(values = COLOR_CONDITION, name = "Condition") +
    labs(x = paste0("PC1 (", ve[1], "%)"), y = paste0("PC2 (", ve[2], "%)")) +
    theme_clean_panel()
}

# =============================================================================
# 3. Heatmap of the conservative signature
# =============================================================================
build_signature_heatmap <- function(d, ds) {
  f <- file.path(tables_dir(ds), paste0(ds$receptor, "_resistance_signature_conservative.csv"))
  sig <- read.csv(f, stringsAsFactors = FALSE)
  sig <- sig[sig$Classification %in% c("Upregulated", "Downregulated"), c("Gene", "Classification")]

  # A gene can be listed once per drug. Keep it once, and leave out genes whose
  # direction differs between drugs, because they have no single direction.
  n_signature <- length(unique(sig$Gene))
  n_directions <- tapply(sig$Classification, sig$Gene, function(v) length(unique(v)))
  ambiguous <- names(n_directions)[n_directions > 1]
  sig <- sig[!sig$Gene %in% ambiguous & !duplicated(sig$Gene), ]
  if (length(ambiguous) > 0) {
    cat("   genes left out for opposite directions between drugs:", length(ambiguous), "\n")
  }
  sig <- sig[sig$Gene %in% rownames(d$mat), ]
  raw <- d$mat[sig$Gene, , drop = FALSE]
  ok  <- apply(raw, 1, var, na.rm = TRUE) > 0
  sig <- sig[ok, ]
  raw <- raw[ok, , drop = FALSE]
  cat("   signature genes drawn:", nrow(sig), "of", n_signature, "\n")

  z <- t(scale(t(raw)))
  direction <- factor(ifelse(sig$Classification == "Upregulated", "Up", "Down"),
                      levels = c("Up", "Down"))
  smp <- d$smp[match(colnames(z), d$smp$label), ]

  ann_col <- HeatmapAnnotation(
    Condition = smp$condition,
    col = list(Condition = COLOR_CONDITION),
    show_annotation_name = FALSE,
    simple_anno_size = unit(3, "mm")
  )
  ann_row <- rowAnnotation(
    Direction = direction,
    col = list(Direction = c(Up = COLOR_UP, Down = COLOR_DOWN)),
    show_annotation_name = FALSE,
    simple_anno_size = unit(3, "mm")
  )

  Heatmap(
    z,
    name = "Z-score\n(per gene)",
    col = circlize::colorRamp2(c(-Z_MAX, 0, Z_MAX), GRADIENT_RNA),
    row_split = direction,
    cluster_row_slices = FALSE,
    cluster_rows = TRUE,
    show_row_dend = FALSE,
    show_row_names = FALSE,
    row_title = NULL,
    row_gap = unit(1.5, "mm"),
    cluster_columns = TRUE,
    top_annotation = ann_col,
    left_annotation = ann_row,
    column_names_gp = gpar(fontsize = 10),
    rect_gp = gpar(col = NA)
  )
}

# =============================================================================
# 4. MA plot
# =============================================================================
plot_ma <- function(ds) {
  f <- file.path(tables_dir(ds), "02_lfc_master_table.csv")
  if (is.na(ds$lfc) || !file.exists(f)) return(NULL)
  df <- read.csv(f)
  df$M <- df[[ds$lfc]]
  df$A <- (df[[ds$base]] + df$Mean_Expr_Res) / 2
  df <- df[is.finite(df$M) & is.finite(df$A), ]

  df$group <- factor(ifelse(df$M >= LOGFC_THRESHOLD, "Up",
                     ifelse(df$M <= -LOGFC_THRESHOLD, "Down", "Not selected")),
                     levels = c("Not selected", "Down", "Up"))
  df <- df[order(df$group), ]  # selected genes are drawn on top

  # The legend carries no gene counts: the signature table applies extra
  # criteria, so its size differs from the number of genes beyond the threshold.
  ggplot(df, aes(A, M, color = group)) +
    geom_point(size = 0.6, alpha = 0.6) +
    geom_hline(yintercept = c(-LOGFC_THRESHOLD, LOGFC_THRESHOLD),
               linetype = "dashed", color = "grey40", linewidth = 0.4) +
    scale_color_manual(values = c(Up = COLOR_UP, Down = COLOR_DOWN, "Not selected" = COLOR_NS),
                       breaks = c("Up", "Down"),
                       name = NULL) +
    guides(color = guide_legend(override.aes = list(size = 2.5, alpha = 1))) +
    labs(x = "Mean expression (log2CPM)", y = "log2 fold change") +
    theme_clean_panel()
}

# =============================================================================
# 5. Run
# =============================================================================
keys <- list()
for (ds in DATASETS) {
  id <- dataset_id(ds)
  cat("\n==", id, "\n")
  tryCatch({
    d <- NULL
    if (any(c("pca", "signature") %in% ds$figs)) {
      d <- load_dataset(ds)
      keys[[id]] <- d$key
    }

    if ("pca" %in% ds$figs) {
      f <- file.path(OUT_DIR, paste0(id, "_pca.png"))
      ggsave(f, plot_pca(d), width = 5.5, height = 4.2, dpi = 300)
      cat("  ", f, "\n")
    }

    if ("signature" %in% ds$figs) {
      f <- file.path(OUT_DIR, paste0(id, "_signature_heatmap.png"))
      png(f, width = 6.5, height = 8, units = "in", res = 300)
      draw(build_signature_heatmap(d, ds), heatmap_legend_side = "right",
           annotation_legend_side = "right")
      dev.off()
      cat("  ", f, "\n")
    }

    if ("ma" %in% ds$figs) {
      f <- file.path(OUT_DIR, paste0(id, "_ma.png"))
      ggsave(f, plot_ma(ds), width = 5.5, height = 4.2, dpi = 300)
      cat("  ", f, "\n")
    }
  }, error = function(e) {
    if (dev.cur() > 1) dev.off()
    cat("   ERROR in", id, ":", conditionMessage(e), "\n")
  })
}

key_file <- file.path(OUT_DIR, "_label_key.csv")
write.csv(bind_rows(keys), key_file, row.names = FALSE)
cat("\n", key_file, "(key to the anonymous labels: do not publish)\n")

cat("\n=== clean_figures.R completed ===\n")
