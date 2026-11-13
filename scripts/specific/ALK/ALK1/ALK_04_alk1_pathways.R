# ============================================================
# ALK_04_alk1_pathways.R
# Purpose: functional enrichment in ALK1. ORA (GO BP and KEGG) of the
#   Upregulated/Downregulated genes, GSEA (GO BP and KEGG) ranked by
#   LFC_Res_vs_Trt, cnetplot and STRINGdb network (v12.0, score >= 400).
#   Dedicated per-model version, derived from scripts/simple/rna_04.1_pathways.R.
# Inputs: results/ALK/ALK1/tables/02_lfc_master_table.csv;
#   config/colors.R, config/gene_lists.R
# Outputs (results/ALK/ALK1/): tables/04_*.csv, figures/04_*.png,
#   report/04_pathways_report.md
# Usage: source(here::here("scripts/specific/ALK/ALK1/ALK_04_alk1_pathways.R"))
#   or Rscript from any folder of the repo.
# ============================================================
CELL_LINE        <- "ALK1"
DRUG             <- "Lorlatinib"
RECEPTOR         <- "ALK"
TOP_N_PATHWAYS   <- 20
MIN_GENES_SET    <- 10
# ============================================================
RESULTS_DIR <- here::here("results", RECEPTOR, CELL_LINE)
CONFIG_DIR  <- here::here("config")
source(file.path(CONFIG_DIR, "colors.R"))
source(file.path(CONFIG_DIR, "gene_lists.R"))

OUT_PREFIX <- paste0(CELL_LINE, "_", DRUG)

# LIBRARIES ------------------------------------------------------------------
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(dplyr)

cat("\n================================================================\n")
cat("  PATHWAYS —", OUT_PREFIX, "\n")
cat("================================================================\n\n")

# LOAD TABLES --------------------------------------------------------------
lfc_file <- file.path(RESULTS_DIR, "tables", "02_lfc_master_table.csv")
if (!file.exists(lfc_file)) stop("02_lfc_master_table.csv not found — run ALK_02_alk1_lfc.R")
lfc_table <- read.csv(lfc_file, check.names = FALSE)
cat("Table 02_lfc_master_table.csv loaded\n\n")

# GENE LISTS --------------------------------------------------------------
genes_up   <- lfc_table %>% filter(Classification == "Upregulated") %>% pull(Gene)
genes_down <- lfc_table %>% filter(Classification == "Downregulated") %>% pull(Gene)
cat("  Upregulated:", length(genes_up),
    "| Downregulated:", length(genes_down), "\n\n", sep = "")

# ORA FUNCTION ----------------------------------------------------------------
run_enrichment <- function(gene_symbols, label, results_dir, top_n) {
  cat("  ORA:", label, "(", length(gene_symbols), "genes)\n")
  if (length(gene_symbols) < 5) {
    cat("  ⚠ Fewer than 5 genes — skipping\n\n"); return(NULL)
  }

  entrez <- bitr(gene_symbols, fromType = "SYMBOL", toType = "ENTREZID",
                 OrgDb = org.Hs.eg.db)
  if (nrow(entrez) == 0) { cat("  ⚠ No Entrez conversion\n\n"); return(NULL) }

  results <- list()
  prefix  <- tolower(gsub(" ", "_", label))

  tryCatch({
    ego <- enrichGO(gene = entrez$ENTREZID, OrgDb = org.Hs.eg.db, ont = "BP",
                    pAdjustMethod = "BH", pvalueCutoff = 0.05,
                    qvalueCutoff = 0.2, readable = TRUE)
    if (!is.null(ego) && nrow(ego@result) > 0) {
      p_go <- dotplot(ego, showCategory = min(top_n, nrow(ego@result))) +
        labs(title = paste("GO BP —", label, "—", OUT_PREFIX)) +
        theme(plot.title = element_text(face = "bold"),
              panel.background = element_rect(fill = "white"),
              plot.background  = element_rect(fill = "white"))
      ggsave(file.path(results_dir, "figures", paste0("04_go_bp_", prefix, ".png")),
             p_go, width = 10, height = 8, dpi = 300)
      cat("    → 04_go_bp_", prefix, ".png\n", sep = "")
      write.csv(ego@result,
                file.path(results_dir, "tables", paste0("04_go_bp_", prefix, ".csv")),
                row.names = FALSE)
      results$go <- ego
    } else { cat("    ⚠ GO BP: no results\n") }
  }, error = function(e) cat("    ⚠ GO BP error:", conditionMessage(e), "\n"))

  tryCatch({
    ekegg <- enrichKEGG(gene = entrez$ENTREZID, organism = "hsa",
                        pAdjustMethod = "BH", pvalueCutoff = 0.05,
                        qvalueCutoff = 0.2)
    if (!is.null(ekegg) && nrow(ekegg@result) > 0) {
      p_kegg <- dotplot(ekegg, showCategory = min(top_n, nrow(ekegg@result))) +
        labs(title = paste("KEGG —", label, "—", OUT_PREFIX)) +
        theme(plot.title = element_text(face = "bold"),
              panel.background = element_rect(fill = "white"),
              plot.background  = element_rect(fill = "white"))
      ggsave(file.path(results_dir, "figures", paste0("04_kegg_", prefix, ".png")),
             p_kegg, width = 10, height = 8, dpi = 300)
      cat("    → 04_kegg_", prefix, ".png\n", sep = "")
      write.csv(ekegg@result,
                file.path(results_dir, "tables", paste0("04_kegg_", prefix, ".csv")),
                row.names = FALSE)
      results$kegg <- ekegg
    } else { cat("    ⚠ KEGG: no results\n") }
  }, error = function(e) cat("    ⚠ KEGG error:", conditionMessage(e), "\n"))

  cat("\n")
  results
}

# PRIMARY ORA ---------------------------------------------------------------
up_label <- "Upregulated"
cat("1. ORA —", up_label, "\n")
res_up <- run_enrichment(genes_up, up_label, RESULTS_DIR, TOP_N_PATHWAYS)

cat("2. ORA — Downregulated\n")
res_down <- NULL
if (length(genes_down) >= MIN_GENES_SET) {
  res_down <- run_enrichment(genes_down, "Downregulated", RESULTS_DIR, TOP_N_PATHWAYS)
} else {
  cat("  ⚠ Only", length(genes_down), "genes — skipping\n\n")
}

# GSEA -----------------------------------------------------------------------
cat("3. GSEA — ranking by LFC_Res_vs_Trt\n")

gene_list <- lfc_table %>%
  filter(!is.na(LFC_Res_vs_Trt)) %>%
  arrange(desc(LFC_Res_vs_Trt))

entrez_all <- bitr(gene_list$Gene, fromType = "SYMBOL", toType = "ENTREZID",
                   OrgDb = org.Hs.eg.db)
gene_list_ranked <- gene_list %>%
  inner_join(entrez_all, by = c("Gene" = "SYMBOL")) %>%
  arrange(desc(LFC_Res_vs_Trt))

ranked_vec <- setNames(gene_list_ranked$LFC_Res_vs_Trt, gene_list_ranked$ENTREZID)
ranked_vec <- ranked_vec[!duplicated(names(ranked_vec))]

res_gsea_go   <- NULL
res_gsea_kegg <- NULL

tryCatch({
  set.seed(1)
  gsea_go <- gseGO(geneList     = ranked_vec,
                   OrgDb        = org.Hs.eg.db,
                   ont          = "BP",
                   minGSSize    = 10,
                   maxGSSize    = 500,
                   pvalueCutoff = 0.05,
                   verbose      = FALSE)

  if (!is.null(gsea_go) && nrow(gsea_go@result) > 0) {
    gsea_go <- setReadable(gsea_go, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")

    p_gsea <- dotplot(gsea_go, showCategory = min(TOP_N_PATHWAYS,
                                                   nrow(gsea_go@result)),
                      split = ".sign") +
      facet_grid(. ~ .sign) +
      labs(title = paste("GSEA GO BP —", OUT_PREFIX)) +
      theme(plot.title       = element_text(face = "bold"),
            panel.background = element_rect(fill = "white"),
            plot.background  = element_rect(fill = "white"))

    ggsave(file.path(RESULTS_DIR, "figures", "04_gsea_go_bp.png"),
           p_gsea, width = 14, height = 8, dpi = 300)
    cat("  → 04_gsea_go_bp.png\n")

    write.csv(gsea_go@result,
              file.path(RESULTS_DIR, "tables", "04_gsea_go_bp.csv"),
              row.names = FALSE)
    res_gsea_go <- gsea_go
  } else { cat("  ⚠ GSEA GO: no significant results\n") }
}, error = function(e) cat("  ⚠ GSEA GO error:", conditionMessage(e), "\n"))

tryCatch({
  set.seed(1)
  gsea_kegg <- gseKEGG(geneList     = ranked_vec,
                        organism     = "hsa",
                        minGSSize    = 10,
                        maxGSSize    = 500,
                        pvalueCutoff = 0.05,
                        verbose      = FALSE)

  if (!is.null(gsea_kegg) && nrow(gsea_kegg@result) > 0) {
    p_gsea_k <- dotplot(gsea_kegg,
                        showCategory = min(TOP_N_PATHWAYS, nrow(gsea_kegg@result)),
                        split = ".sign") +
      facet_grid(. ~ .sign) +
      labs(title = paste("GSEA KEGG —", OUT_PREFIX)) +
      theme(plot.title       = element_text(face = "bold"),
            panel.background = element_rect(fill = "white"),
            plot.background  = element_rect(fill = "white"))

    ggsave(file.path(RESULTS_DIR, "figures", "04_gsea_kegg.png"),
           p_gsea_k, width = 14, height = 8, dpi = 300)
    cat("  → 04_gsea_kegg.png\n")

    write.csv(gsea_kegg@result,
              file.path(RESULTS_DIR, "tables", "04_gsea_kegg.csv"),
              row.names = FALSE)
    res_gsea_kegg <- gsea_kegg
  } else { cat("  ⚠ GSEA KEGG: no significant results\n") }
}, error = function(e) cat("  ⚠ GSEA KEGG error:", conditionMessage(e), "\n"))
cat("\n")

# CNETPLOT (enrichplot) ------------------------------------------------------
cat("4. Network — cnetplot (GO BP)...\n")

if (!is.null(res_up$go)) {
  tryCatch({
    p_cnet <- cnetplot(res_up$go,
                       showCategory  = 5,
                       color_gene    = "steelblue",
                       color_category = classification_colors["Upregulated"]) +
      labs(title = paste("Gene-concept network —", up_label, "—", OUT_PREFIX)) +
      theme(plot.title = element_text(face = "bold"),
            plot.background = element_rect(fill = "white"))
    ggsave(file.path(RESULTS_DIR, "figures", "04_cnetplot_upregulated.png"),
           p_cnet, width = 12, height = 10, dpi = 300)
    cat("  → 04_cnetplot_upregulated.png\n")
  }, error = function(e) cat("  ⚠ cnetplot error:", conditionMessage(e), "\n"))
} else {
  cat("  ⚠ No GO results for cnetplot\n")
}

if (!is.null(res_down$go)) {
  tryCatch({
    p_cnet_dn <- cnetplot(res_down$go,
                          showCategory   = 5,
                          color_gene     = "steelblue",
                          color_category = classification_colors["Downregulated"]) +
      labs(title = paste("Gene-concept network — Downregulated —", OUT_PREFIX)) +
      theme(plot.title = element_text(face = "bold"),
            plot.background = element_rect(fill = "white"))
    ggsave(file.path(RESULTS_DIR, "figures", "04_cnetplot_downregulated.png"),
           p_cnet_dn, width = 12, height = 10, dpi = 300)
    cat("  → 04_cnetplot_downregulated.png\n")
  }, error = function(e) cat("  ⚠ cnetplot error:", conditionMessage(e), "\n"))
}
cat("\n")

# STRINGDB NETWORK -----------------------------------------------------------
cat("5. Network — STRINGdb (", up_label, ")...\n")

tryCatch({
  library(STRINGdb)

  string_db <- STRINGdb$new(
    version    = "12.0",
    species    = 9606,
    score_threshold = 400,
    network_type = "functional",
    input_directory = file.path(RESULTS_DIR, "tables")
  )

  sig_df <- data.frame(Gene = genes_up)
  sig_mapped <- string_db$map(sig_df, "Gene", removeUnmappedRows = TRUE)

  if (nrow(sig_mapped) >= 3) {
    png(file.path(RESULTS_DIR, "figures", "04_stringdb_upregulated.png"),
        width = 2400, height = 2000, res = 200)
    string_db$plot_network(sig_mapped$STRING_id)
    dev.off()
    cat("  → 04_stringdb_upregulated.png\n")

    enrichment_string <- string_db$get_enrichment(sig_mapped$STRING_id)
    if (!is.null(enrichment_string) && nrow(enrichment_string) > 0) {
      write.csv(enrichment_string,
                file.path(RESULTS_DIR, "tables", "04_stringdb_enrichment.csv"),
                row.names = FALSE)
      cat("  → 04_stringdb_enrichment.csv\n")
    }
  } else {
    cat("  ⚠ Fewer than 3 genes mapped in STRINGdb\n")
  }
}, error = function(e) cat("  ⚠ STRINGdb error:", conditionMessage(e), "\n"))
cat("\n")

# REPORT --------------------------------------------------------------------
extract_top <- function(res, type = "go", n = 10) {
  obj <- res[[type]]
  if (is.null(obj) || nrow(obj@result) == 0) return("Sin resultados significativos")
  top <- head(obj@result %>% arrange(p.adjust) %>%
                dplyr::select(Description, p.adjust, Count), n)
  paste(capture.output(print(top, row.names = FALSE)), collapse = "\n")
}

extract_top_gsea <- function(res, n = 10) {
  if (is.null(res) || nrow(res@result) == 0) return("Sin resultados significativos")
  top <- head(res@result %>% arrange(p.adjust) %>%
                dplyr::select(Description, NES, p.adjust), n)
  paste(capture.output(print(top, row.names = FALSE)), collapse = "\n")
}

report_text <- paste0(
  "# Pathways Report — ", OUT_PREFIX, "\n\n",
  "## ORA GO BP — ", up_label, " (top 10)\n",
  extract_top(res_up, "go", 10), "\n\n",
  "## ORA KEGG — ", up_label, " (top 10)\n",
  extract_top(res_up, "kegg", 10), "\n\n",
  "## GSEA GO BP (top 10 por p.adjust)\n",
  extract_top_gsea(res_gsea_go, 10), "\n\n",
  "## GSEA KEGG (top 10 por p.adjust)\n",
  extract_top_gsea(res_gsea_kegg, 10), "\n\n",
  "## Figuras\n",
  "- 04_go_bp_", tolower(gsub(" ", "_", up_label)), ".png\n",
  "- 04_kegg_", tolower(gsub(" ", "_", up_label)), ".png\n",
  if (!is.null(res_down)) "- 04_go_bp_downregulated.png\n- 04_kegg_downregulated.png\n" else "",
  "- 04_gsea_go_bp.png\n",
  "- 04_gsea_kegg.png\n",
  "- 04_cnetplot_upregulated.png\n",
  "- 04_cnetplot_downregulated.png\n",
  "- 04_stringdb_upregulated.png\n"
)
writeLines(report_text, file.path(RESULTS_DIR, "report", "04_pathways_report.md"))
cat("Report: report/04_pathways_report.md\n\n")

cat("================================================================\n")
cat("  PATHWAYS completed —", OUT_PREFIX, "\n")
cat("================================================================\n\n")
