# =============================================================================
# ALK_09b_kegg_vs_untreated_exploratory.R
# Purpose: exploratory script, NOT part of the main pipeline. It does not
#   modify the outputs of ALK_02/ALK_04/ALK_09; all its files carry their own
#   prefix/suffix.
#
# Runs KEGG + GO Biological Process with THREE different Up/Down
# classification bases, to compare which pathways/processes emerge depending
# on the comparator:
#
#   1) vs_no_tratada -- LFC_Res_vs_Unt (total change of the resistant line
#      versus the untreated parental, the metric of the per-sample LFC
#      heatmaps)
#   2) vs_tratada     -- LFC_Res_vs_Trt (the main comparator of ALK_02/
#      ALK_04, isolates the resistance-specific effect from the acute drug
#      effect)
#   3) combinado      -- intersection of the two above (a gene counts as
#      Up/Down only if it passes the threshold in BOTH, same direction in
#      both). It combines the total change with the filter against the
#      treated parental, which removes genes whose "total" change is actually
#      inherited from the acute drug effect (significant in vs_no_tratada but
#      flat in vs_tratada). NOTE: this does NOT resolve heterogeneity between
#      resistant populations of the same model (e.g. KRAS/BCAT1 in ALK2,
#      flat in the pooled Res_vs_Trt although the change in L4/L5 is real);
#      that still requires manually forcing genes, it is an aggregation
#      problem, not a matter of which comparator is used.
#
# Inputs: results/ALK/{ALK1,ALK2,ALK3}/tables/02_lfc_master_table.csv
# Outputs: 3 KEGG figures + 3 GO BP figures (6 in total), each with the
#   4 facets (ALK1, ALK2, ALK3-Lorlatinib, ALK3-Alectinib) as subpanels:
#   results/ALK/combined/figures/ALK_09b_{kegg,go}_<mode>.png
#   results/ALK/combined/tables/ALK_09b_{kegg,go}_top8_<mode>.csv
#   Full CSVs (not only top8) per facet/direction/mode in
#   results/ALK/<model>/tables/04_{kegg,go}_{up,down}regulated_<mode>_<facet>.csv
# Usage: source(here::here("scripts/specific/ALK/ALK_09b_kegg_vs_untreated_exploratory.R"))
#   or Rscript from any folder of the repo.
# =============================================================================

suppressMessages({
  library(dplyr)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ggplot2)
  library(patchwork)
})

OUT_DIR <- here::here("results", "ALK", "combined")
dir.create(file.path(OUT_DIR, "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT_DIR, "tables"),  recursive = TRUE, showWarnings = FALSE)

LOGFC_THRESHOLD <- 2.0  # same threshold as ALK_02_*_lfc.R / classify_lfc()
TOP_N           <- 8    # same top-N as the main KEGG dotplot, for both KEGG and GO

COLOR_KEGG_SIG   <- "#D6604D"
COLOR_KEGG_NOSIG <- "#4393C3"

FACETAS <- list(
  ALK1 = list(ruta = here::here("results", "ALK", "ALK1", "tables", "02_lfc_master_table.csv"), drug = NULL),
  ALK2 = list(ruta = here::here("results", "ALK", "ALK2", "tables", "02_lfc_master_table.csv"), drug = NULL),
  `ALK3 - Lorlatinib` = list(ruta = here::here("results", "ALK", "ALK3", "tables", "02_lfc_master_table.csv"), drug = "Lorlatinib"),
  `ALK3 - Alectinib`  = list(ruta = here::here("results", "ALK", "ALK3", "tables", "02_lfc_master_table.csv"), drug = "Alectinib")
)

RUTAS_IMPORTANTES_ACTUAL <- c(
  "Cell cycle", "PI3K-Akt signaling pathway", "MAPK signaling pathway",
  "Ras signaling pathway", "Small cell lung cancer", "p53 signaling pathway",
  "ABC transporters", "Focal adhesion", "Wnt signaling pathway",
  "Hippo signaling pathway", "TGF-beta signaling pathway", "Rap1 signaling pathway",
  "Calcium signaling pathway", "ECM-receptor interaction", "Integrin signaling",
  "Cadherin signaling"
)

LFC_MODOS <- c("vs_no_tratada", "vs_tratada", "combinado")

# Classify Up/Down according to the mode. "combinado" requires passing the
# threshold in both LFC_Res_vs_Unt AND LFC_Res_vs_Trt, same direction in both.
clasificar_genes <- function(df, modo, threshold) {
  up_unt   <- df$LFC_Res_vs_Unt >  threshold
  down_unt <- df$LFC_Res_vs_Unt < -threshold
  up_trt   <- df$LFC_Res_vs_Trt >  threshold
  down_trt <- df$LFC_Res_vs_Trt < -threshold

  if (modo == "vs_no_tratada") {
    list(up = df$Gene[up_unt], down = df$Gene[down_unt])
  } else if (modo == "vs_tratada") {
    list(up = df$Gene[up_trt], down = df$Gene[down_trt])
  } else if (modo == "combinado") {
    list(up = df$Gene[up_unt & up_trt], down = df$Gene[down_unt & down_trt])
  } else {
    stop("Unknown mode: ", modo)
  }
}

run_kegg <- function(gene_symbols, label) {
  cat("  KEGG:", label, "(", length(gene_symbols), "genes)\n")
  if (length(gene_symbols) < 5) { cat("   Fewer than 5 genes, skipping\n"); return(NULL) }
  entrez <- bitr(gene_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  if (nrow(entrez) == 0) { cat("   No Entrez conversion\n"); return(NULL) }
  ekegg <- tryCatch(
    enrichKEGG(gene = entrez$ENTREZID, organism = "hsa",
               pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2),
    error = function(e) { cat("   Error enrichKEGG:", conditionMessage(e), "\n"); NULL }
  )
  if (is.null(ekegg) || nrow(ekegg@result) == 0) { cat("   No significant pathways\n"); return(NULL) }
  ekegg@result
}

run_go <- function(gene_symbols, label) {
  cat("  GO BP:", label, "(", length(gene_symbols), "genes)\n")
  if (length(gene_symbols) < 5) { cat("   Fewer than 5 genes, skipping\n"); return(NULL) }
  entrez <- bitr(gene_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  if (nrow(entrez) == 0) { cat("   No Entrez conversion\n"); return(NULL) }
  ego <- tryCatch(
    enrichGO(gene = entrez$ENTREZID, OrgDb = org.Hs.eg.db, ont = "BP",
             pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = TRUE),
    error = function(e) { cat("   Error enrichGO:", conditionMessage(e), "\n"); NULL }
  )
  if (is.null(ego) || nrow(ego@result) == 0) { cat("   No significant processes\n"); return(NULL) }
  ego@result
}

negrita_actuales <- function(labels) {
  lapply(labels, function(l) if (l %in% RUTAS_IMPORTANTES_ACTUAL) bquote(bold(.(l))) else l)
}

plot_panel <- function(df_sub, titulo) {
  rango_color <- range(df_sub$p.adjust)
  rango_count <- range(df_sub$Count)
  ggplot(df_sub, aes(x = FoldEnrichment, y = reorder(Description, FoldEnrichment),
                      size = Count, fill = p.adjust, shape = Direction)) +
    geom_point(color = "grey20", stroke = 0.3) +
    scale_fill_gradient(low = COLOR_KEGG_SIG, high = COLOR_KEGG_NOSIG, trans = "log10",
                         name = "p.adjust", limits = rango_color,
                         guide = guide_colorbar(reverse = TRUE, order = 3)) +
    scale_size_continuous(range = c(2, 8), name = "Genes", limits = rango_count,
                           guide = guide_legend(override.aes = list(shape = 24, fill = "grey40"), order = 2)) +
    scale_shape_manual(values = c(Upregulated = 24, Downregulated = 25), name = "Direction",
                        guide = guide_legend(order = 1)) +
    scale_y_discrete(labels = negrita_actuales) +
    labs(title = titulo, x = "Fold enrichment", y = NULL) +
    theme_minimal(base_size = 10) +
    theme(
      axis.text.y  = element_text(size = 7.5),
      plot.title   = element_text(face = "bold", size = 11),
      legend.position = "right",
      legend.box   = "vertical",
      panel.background = element_rect(fill = "white"),
      plot.background  = element_rect(fill = "white")
    )
}

# Run one analysis type (KEGG or GO) for the 4 facets with a given
# classification mode, save the full CSVs, and return the combined top8
# (Up+Down per facet) ready for the dotplot.
correr_analisis <- function(tipo, modo, run_fun) {
  cat("\n===", tipo, "-", modo, "===\n")
  bind_rows(lapply(names(FACETAS), function(f) {
    cat(f, "\n")
    e  <- FACETAS[[f]]
    df <- read.csv(e$ruta, stringsAsFactors = FALSE)
    if (!is.null(e$drug)) df <- df %>% filter(Drug == e$drug)

    genes <- clasificar_genes(df, modo, LOGFC_THRESHOLD)
    cat("  Upregulated:", length(genes$up), "| Downregulated:", length(genes$down), "\n")

    res_up   <- run_fun(genes$up,   paste(f, modo, "Upregulated"))
    res_down <- run_fun(genes$down, paste(f, modo, "Downregulated"))

    slug <- tolower(gsub("[^A-Za-z0-9]+", "_", f))
    if (!is.null(res_up)) {
      write.csv(res_up, file.path(dirname(e$ruta), paste0("04_", tipo, "_upregulated_", modo, "_", slug, ".csv")), row.names = FALSE)
    }
    if (!is.null(res_down)) {
      write.csv(res_down, file.path(dirname(e$ruta), paste0("04_", tipo, "_downregulated_", modo, "_", slug, ".csv")), row.names = FALSE)
    }

    bind_rows(
      if (!is.null(res_up))   res_up   %>% arrange(p.adjust) %>% slice_head(n = TOP_N) %>%
        transmute(Faceta = f, Direction = "Upregulated", Description, p.adjust, Count, FoldEnrichment) else NULL,
      if (!is.null(res_down)) res_down %>% arrange(p.adjust) %>% slice_head(n = TOP_N) %>%
        transmute(Faceta = f, Direction = "Downregulated", Description, p.adjust, Count, FoldEnrichment) else NULL
    )
  }))
}

generar_figura <- function(top_df, tipo, modo, titulo_tipo) {
  if (nrow(top_df) == 0) {
    cat("  No data for", tipo, modo, "-- figure not generated\n")
    return(invisible(NULL))
  }
  top_df <- top_df %>%
    mutate(
      Faceta    = factor(Faceta, levels = names(FACETAS)),
      Direction = factor(Direction, levels = c("Upregulated", "Downregulated"))
    )

  out_csv <- file.path(OUT_DIR, "tables", paste0("ALK_09b_", tipo, "_top8_", modo, ".csv"))
  write.csv(top_df, out_csv, row.names = FALSE)
  cat("  ", out_csv, "\n")

  facetas_con_datos <- names(FACETAS)[names(FACETAS) %in% unique(top_df$Faceta)]
  paneles <- lapply(facetas_con_datos, function(f) {
    plot_panel(top_df %>% filter(Faceta == f), as.character(f))
  })

  p <- wrap_plots(paneles, ncol = 2) +
    plot_annotation(
      title = paste0(titulo_tipo, " (top 8, Up+Down) — ", modo, " — exploratorio, ALK1 + ALK2 + ALK3"),
      caption = paste0(
        "EXPLORATORIO -- clasificacion Up/Down por '", modo, "'. ",
        "Circulo = Upregulated, triangulo = Downregulated. En negrita: rutas ya en RUTAS_IMPORTANTES (ALK_09, solo aplica a KEGG)."
      ),
      theme = theme(
        plot.title   = element_text(face = "bold", size = 14),
        plot.caption = element_text(hjust = 0, size = 8.5, color = "grey35")
      )
    )

  out_png <- file.path(OUT_DIR, "figures", paste0("ALK_09b_", tipo, "_", modo, ".png"))
  ggsave(out_png, p, width = 15, height = 15, dpi = 300)
  cat("  ", out_png, "\n")
}

for (modo in LFC_MODOS) {
  kegg_top <- correr_analisis("kegg", modo, run_kegg)
  generar_figura(kegg_top, "kegg", modo, "KEGG enrichment")

  go_top <- correr_analisis("go", modo, run_go)
  generar_figura(go_top, "go", modo, "GO Biological Process")
}

cat("\nDone. 3 KEGG figures + 3 GO figures in", file.path(OUT_DIR, "figures"), "\n")
cat("Full CSVs (beyond the top8) in tables/ of each model, suffix _<mode>_<facet>.csv\n")
