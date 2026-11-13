# ============================================================
# rna_04.1_pathways.R
# Functional enrichment: ORA (GO BP/KEGG) of the up/down genes, GSEA
# (GO BP/KEGG) on the full LFC ranking, cnetplot and STRINGdb network.
# Generic: 2-class scheme (TWO_CLASS_DATASETS) or 5-class scheme.
#
# Inputs:    results/<RECEPTOR>/<CELL_LINE>/tables/02_lfc_master_table.csv
#            <PREFIX>_resistance_signature[_conservative].csv (2-class scheme)
#            01_mat_tmm.csv, 01_metadata_unit.csv (MET1_METi only)
#            config/colors.R, config/gene_lists.R
# Outputs:   results/<RECEPTOR>/<CELL_LINE>/tables/04_*.csv
#            results/<RECEPTOR>/<CELL_LINE>/figures/04_*.png
#            results/<RECEPTOR>/<CELL_LINE>/report/04_pathways_report.{md,html}
# Usage:     set the CONFIG block and, from any folder of the repo,
#            source(here::here("scripts/simple/rna_04.1_pathways.R"))
#            (requires an internet connection for KEGG and STRINGdb)
# ============================================================
# CONFIG — change here for each run
# ============================================================
DATASET         <- "EGFR1_EGFRi"       # "ALK12" | "ALK3" | "MET1_METi" | "EGFR1_EGFRi"
CELL_LINE       <- "EGFR1"           # "ALK2" | "ALK1" | "ALK3" | "MET1" | "EGFR1"
DRUG            <- "Gefitinib_Afatinib"  # "Lorlatinib" | "Alectinib" | "MET_inhibitor" | "Gefitinib_Afatinib"
RECEPTOR        <- "EGFR"           # "MET" | "ALK" | "EGFR"
LOGFC_THRESHOLD <- 2.0
TOP_N_PATHWAYS  <- 20
MIN_GENES_SET   <- 10
USE_CONSERVATIVE <- TRUE   # TRUE = conservative signature
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
if (!file.exists(lfc_file)) stop("02_lfc_master_table.csv not found — run rna_02")
lfc_table <- read.csv(lfc_file, check.names = FALSE)
cat("Table 02_lfc_master_table.csv loaded\n\n")

# GENE LISTS ------------------------------------------------------------
if (DATASET %in% TWO_CLASS_DATASETS) {
  sig_prefix <- if (DATASET == "MET1_METi") "MET" else RECEPTOR
  sig_file <- if (USE_CONSERVATIVE) {
    paste0(sig_prefix, "_resistance_signature_conservative.csv")
  } else {
    paste0(sig_prefix, "_resistance_signature.csv")
  }
  sig_path <- file.path(RESULTS_DIR, "tables", sig_file)
  if (!file.exists(sig_path)) stop("Signature not found: ", sig_path,
                                   "\nRun rna_03.1_key_genes.R first")
  signature  <- read.csv(sig_path)
  genes_up   <- signature %>% filter(Classification == "Upregulated") %>% pull(Gene)
  genes_down <- signature %>% filter(Classification == "Downregulated") %>% pull(Gene)

  cat("  Signature:", sig_file, "\n")
  cat("  Upregulated:", length(genes_up), "| Downregulated:", length(genes_down), "\n\n")

  if (DATASET == "MET1_METi") {
    # GX_cluster per sample — to compare Group1 vs Group2 genes (04_genes_by_gx_cluster.csv)
    meta_file <- file.path(RESULTS_DIR, "tables", "01_metadata_unit.csv")
    meta_df   <- read.csv(meta_file, check.names = FALSE)
    rownames(meta_df) <- meta_df$Model

    if (!"GX_cluster" %in% colnames(meta_df)) {
      rg_map <- c(
        MET1     = "Parental",
        MET1_C1  = "Group1_EGFR_bypass", MET1_PR1 = "Group1_EGFR_bypass",
        MET1_PR2 = "Group1_EGFR_bypass", MET1_V8  = "Group1_EGFR_bypass",
        MET1_PR3 = "Group2_NF2_YAP",     MET1_R = "Group2_NF2_YAP"
      )
      meta_df$GX_cluster <- rg_map[rownames(meta_df)]
    }

    mat_df <- read.csv(file.path(RESULTS_DIR, "tables", "01_mat_tmm.csv"),
                       check.names = FALSE)
    rownames(mat_df) <- mat_df$Gene
    mat <- as.matrix(mat_df[, -1])

    g1_samples  <- rownames(meta_df)[meta_df$GX_cluster == "Group1_EGFR_bypass"]
    g2_samples  <- rownames(meta_df)[meta_df$GX_cluster == "Group2_NF2_YAP"]
    par_sample  <- rownames(meta_df)[meta_df$GX_cluster == "Parental"]
  }

} else {
  genes_acquired    <- lfc_table %>%
    filter(Classification == unname(CLASS_LABELS["acq"]), LFC_Res_vs_Trt > LOGFC_THRESHOLD) %>%
    pull(Gene)
  genes_preexisting <- lfc_table %>%
    filter(Classification == unname(CLASS_LABELS["pre"])) %>% pull(Gene)
  genes_down        <- lfc_table %>%
    filter(Classification == "Downregulated") %>% pull(Gene)
  genes_up          <- genes_acquired   # alias — the rest of the script uses genes_up generically
  cat("  ", CLASS_LABELS["acq"], ":", length(genes_acquired),
      "| ", CLASS_LABELS["pre"], ":", length(genes_preexisting),
      "| Downregulated:", length(genes_down), "\n\n", sep = "")
}

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
up_label <- if (DATASET %in% TWO_CLASS_DATASETS) "Upregulated" else unname(CLASS_LABELS["acq"])
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
# Ranking: all genes by LFC from highest to lowest (LFC_Res_vs_Par in the
# 2-class scheme, LFC_Res_vs_Trt in the 5-class scheme)
cat("3. GSEA — LFC ranking\n")

lfc_col <- if (DATASET %in% TWO_CLASS_DATASETS) "LFC_Res_vs_Par" else "LFC_Res_vs_Trt"

gene_list <- lfc_table %>%
  filter(!is.na(.data[[lfc_col]])) %>%
  arrange(desc(.data[[lfc_col]]))

entrez_all <- bitr(gene_list$Gene, fromType = "SYMBOL", toType = "ENTREZID",
                   OrgDb = org.Hs.eg.db)
gene_list_ranked <- gene_list %>%
  inner_join(entrez_all, by = c("Gene" = "SYMBOL")) %>%
  arrange(desc(.data[[lfc_col]]))

ranked_vec <- setNames(gene_list_ranked[[lfc_col]], gene_list_ranked$ENTREZID)
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

# GROUP1 VS GROUP2 GENES (MET1_METi only; table, no Venn figure) -----------
if (DATASET == "MET1_METi") {
  cat("4. Venn Group1 vs Group2...\n")

  # LFC per group: mean of the samples of each group vs parental
  expr_par <- rowMeans(mat[, par_sample, drop = FALSE])

  lfc_g1 <- if (length(g1_samples) > 1) rowMeans(mat[, g1_samples]) - expr_par else
               mat[, g1_samples] - expr_par
  lfc_g2 <- if (length(g2_samples) > 1) rowMeans(mat[, g2_samples]) - expr_par else
               mat[, g2_samples] - expr_par

  up_g1 <- names(lfc_g1)[lfc_g1 >  LOGFC_THRESHOLD]
  up_g2 <- names(lfc_g2)[lfc_g2 >  LOGFC_THRESHOLD]
  dn_g1 <- names(lfc_g1)[lfc_g1 < -LOGFC_THRESHOLD]
  dn_g2 <- names(lfc_g2)[lfc_g2 < -LOGFC_THRESHOLD]

  cat("  Upregulated — Group1:", length(up_g1), "| Group2:", length(up_g2),
      "| Shared:", length(intersect(up_g1, up_g2)), "\n")
  cat("  Downregulated — Group1:", length(dn_g1), "| Group2:", length(dn_g2),
      "| Shared:", length(intersect(dn_g1, dn_g2)), "\n\n")

  # Save group-specific genes
  write.csv(
    data.frame(
      Gene      = unique(c(up_g1, up_g2, dn_g1, dn_g2)),
      Up_G1     = unique(c(up_g1, up_g2, dn_g1, dn_g2)) %in% up_g1,
      Up_G2     = unique(c(up_g1, up_g2, dn_g1, dn_g2)) %in% up_g2,
      Down_G1   = unique(c(up_g1, up_g2, dn_g1, dn_g2)) %in% dn_g1,
      Down_G2   = unique(c(up_g1, up_g2, dn_g1, dn_g2)) %in% dn_g2
    ),
    file.path(RESULTS_DIR, "tables", "04_genes_by_gx_cluster.csv"),
    row.names = FALSE
  )
}

# CNETPLOT (enrichplot) ------------------------------------------------------
cat("5. Network — cnetplot (GO BP)...\n")

if (!is.null(res_up$go)) {
  tryCatch({
    p_cnet <- cnetplot(res_up$go,
                       showCategory  = 5,
                       color_gene    = "steelblue",
                       color_category = classification_colors["Upregulated"]) +
      labs(title = paste("Gene-concept network — Upregulated —", OUT_PREFIX)) +
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
cat("6. Network — STRINGdb (full signature)...\n")

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

    # STRINGdb enrichment
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
  if (is.null(obj) || nrow(obj@result) == 0) return("No significant results")
  top <- head(obj@result %>% arrange(p.adjust) %>%
                dplyr::select(Description, p.adjust, Count), n)
  paste(capture.output(print(top, row.names = FALSE)), collapse = "\n")
}

extract_top_gsea <- function(res, n = 10) {
  if (is.null(res) || nrow(res@result) == 0) return("No significant results")
  top <- head(res@result %>% arrange(p.adjust) %>%
                dplyr::select(Description, NES, p.adjust), n)
  paste(capture.output(print(top, row.names = FALSE)), collapse = "\n")
}

report_text <- paste0(
  "# Pathways Report — ", OUT_PREFIX, "\n\n",
  if (DATASET == "MET1_METi")
    paste0("**Signature:** ",
           if (USE_CONSERVATIVE) "conservative (Mean_Expr_Par > -5 + forced NF2/NRG1)"
           else "full", "\n\n")
  else "",
  "## ORA GO BP — Upregulated (top 10)\n",
  extract_top(res_up, "go", 10), "\n\n",
  "## ORA KEGG — Upregulated (top 10)\n",
  extract_top(res_up, "kegg", 10), "\n\n",
  "## GSEA GO BP (top 10 by p.adjust)\n",
  extract_top_gsea(res_gsea_go, 10), "\n\n",
  "## GSEA KEGG (top 10 by p.adjust)\n",
  extract_top_gsea(res_gsea_kegg, 10), "\n\n",
  "## Figures\n",
  "- 04_go_bp_upregulated.png\n",
  "- 04_kegg_upregulated.png\n",
  if (!is.null(res_down)) "- 04_go_bp_downregulated.png\n- 04_kegg_downregulated.png\n" else "",
  "- 04_gsea_go_bp.png\n",
  "- 04_gsea_kegg.png\n",
  "- 04_cnetplot_upregulated.png\n",
  "- 04_cnetplot_downregulated.png\n",
  "- 04_stringdb_upregulated.png\n"
)
writeLines(report_text, file.path(RESULTS_DIR, "report", "04_pathways_report.md"))
cat("Report: report/04_pathways_report.md\n\n")

if (requireNamespace("rmarkdown", quietly = TRUE)) {
  cat("Generating HTML report...\n")
  rmarkdown::render(here::here("scripts", "simple", "rna_04_pathways_report.Rmd"),
                    output_dir  = file.path(RESULTS_DIR, "report"),
                    output_file = "04_pathways_report.html",
                    quiet       = TRUE)
  cat("  → report/04_pathways_report.html\n\n")
}

cat("================================================================\n")
cat("  PATHWAYS completed —", OUT_PREFIX, "\n")
cat("================================================================\n\n")
