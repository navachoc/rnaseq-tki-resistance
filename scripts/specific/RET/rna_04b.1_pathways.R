# ============================================================
# rna_04b.1_pathways.R
# Functional enrichment: ORA (GO/KEGG) + networks, RET family
# (RET1, RET2; treated parental per drug). Consumes the signature from
# rna_03b.1_key_genes.R (comparison d, Untreated-vs-Resistance).
# No GSEA: there is no single ranking metric that represents the 4 RET LFC
# comparisons well; the analysis is already driven by the discrete up/down
# signature from 03b.
# RNAseq_bulk pipeline
#
# Inputs: results/RET/<CELL_LINE>/tables/RET_resistance_signature_conservative.csv
#   (or RET_resistance_signature.csv if USE_CONSERVATIVE = FALSE);
#   config/colors.R, config/gene_lists.R
# Outputs (results/RET/<CELL_LINE>/): tables/04_go_bp_*.csv, tables/04_kegg_*.csv,
#   tables/04_stringdb_enrichment.csv, figures/04_*.png,
#   report/04_pathways_report.md
# Dependencies (run first): rna_03b.1_key_genes.R for the same CELL_LINE.
#   One run per model (CELL_LINE = "RET1" | "RET2").
#   Its 04_* tables are read by RET_09_lfc_heatmaps.R.
# Usage: paths are resolved with here::here() from the repository root
# (.here file), so it works from any working directory:
#   Rscript scripts/specific/RET/rna_04b.1_pathways.R
# ============================================================
# CONFIG: change here for each run
# ============================================================
DATASET          <- "RET1"   # "RET1" | "RET2"
CELL_LINE        <- "RET1"
RECEPTOR         <- "RET"
LOGFC_THRESHOLD  <- 2.0
TOP_N_PATHWAYS   <- 20
MIN_GENES_SET    <- 10
USE_CONSERVATIVE <- TRUE   # TRUE = conservative signature
# ============================================================
RESULTS_DIR <- here::here("results", RECEPTOR, CELL_LINE)
CONFIG_DIR  <- here::here("config")
DRUG <- RECEPTOR  # placeholder: same pattern as 02b/03b, there is no single drug
source(file.path(CONFIG_DIR, "colors.R"))
source(file.path(CONFIG_DIR, "gene_lists.R"))
OUT_PREFIX <- paste0(CELL_LINE, "_", RECEPTOR)  # overrides the one from colors.R

# LIBRARIES ------------------------------------------------------------------
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(dplyr)
library(grid)

cat("\n================================================================\n")
cat("  PATHWAYS —", OUT_PREFIX, "\n")
cat("================================================================\n\n")

# SIGNATURE (rna_03b.1_key_genes.R) -----------------------------------------------
sig_file <- if (USE_CONSERVATIVE) {
  paste0(RECEPTOR, "_resistance_signature_conservative.csv")
} else {
  paste0(RECEPTOR, "_resistance_signature.csv")
}
sig_path <- file.path(RESULTS_DIR, "tables", sig_file)
if (!file.exists(sig_path)) stop("Signature not found: ", sig_path,
                                 "\nRun rna_03b.1_key_genes.R first")
signature  <- read.csv(sig_path)
genes_up   <- signature %>% filter(Classification == "Upregulated") %>% pull(Gene)
genes_down <- signature %>% filter(Classification == "Downregulated") %>% pull(Gene)

cat("  Signature:", sig_file, "\n")
cat("  Upregulated:", length(genes_up), "| Downregulated:", length(genes_down), "\n\n")

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
      write.csv(ego@result,
                file.path(results_dir, "tables", paste0("04_go_bp_", prefix, ".csv")),
                row.names = FALSE)
      n_sig <- sum(ego@result$p.adjust <= 0.05)
      # nrow(ego@result) > 0 does not imply significant terms: enrichGO()
      # keeps in @result everything that passes the raw p-value before FDR
      # correction; some terms may not survive the BH correction.
      # dotplot()/cnetplot() filter by p.adjust internally and return an empty
      # figure (without warning) if n_sig == 0, so it must be checked first.
      if (n_sig > 0) {
        p_go <- dotplot(ego, showCategory = min(top_n, nrow(ego@result))) +
          labs(title = paste("GO BP —", label, "—", OUT_PREFIX)) +
          theme(plot.title = element_text(face = "bold"),
                panel.background = element_rect(fill = "white"),
                plot.background  = element_rect(fill = "white"))
        ggsave(file.path(results_dir, "figures", paste0("04_go_bp_", prefix, ".png")),
               p_go, width = 10, height = 8, dpi = 300)
        cat("    → 04_go_bp_", prefix, ".png (", n_sig, "/", nrow(ego@result),
            " terms with p.adjust ≤ 0.05)\n", sep = "")
        results$go <- ego
      } else {
        cat("    ⚠ GO BP: 0/", nrow(ego@result),
            " terms with p.adjust ≤ 0.05 — no figure (table still saved)\n", sep = "")
      }
    } else { cat("    ⚠ GO BP: no results\n") }
  }, error = function(e) cat("    ⚠ GO BP error:", conditionMessage(e), "\n"))

  tryCatch({
    ekegg <- enrichKEGG(gene = entrez$ENTREZID, organism = "hsa",
                        pAdjustMethod = "BH", pvalueCutoff = 0.05,
                        qvalueCutoff = 0.2)
    if (!is.null(ekegg) && nrow(ekegg@result) > 0) {
      write.csv(ekegg@result,
                file.path(results_dir, "tables", paste0("04_kegg_", prefix, ".csv")),
                row.names = FALSE)
      n_sig <- sum(ekegg@result$p.adjust <= 0.05)
      if (n_sig > 0) {
        p_kegg <- dotplot(ekegg, showCategory = min(top_n, nrow(ekegg@result))) +
          labs(title = paste("KEGG —", label, "—", OUT_PREFIX)) +
          theme(plot.title = element_text(face = "bold"),
                panel.background = element_rect(fill = "white"),
                plot.background  = element_rect(fill = "white"))
        ggsave(file.path(results_dir, "figures", paste0("04_kegg_", prefix, ".png")),
               p_kegg, width = 10, height = 8, dpi = 300)
        cat("    → 04_kegg_", prefix, ".png (", n_sig, "/", nrow(ekegg@result),
            " terms with p.adjust ≤ 0.05)\n", sep = "")
        results$kegg <- ekegg
      } else {
        cat("    ⚠ KEGG: 0/", nrow(ekegg@result),
            " terms with p.adjust ≤ 0.05 — no figure (table still saved)\n", sep = "")
      }
    } else { cat("    ⚠ KEGG: no results\n") }
  }, error = function(e) cat("    ⚠ KEGG error:", conditionMessage(e), "\n"))

  cat("\n")
  results
}

# ORA --------------------------------------------------------------------
cat("1. ORA — Upregulated\n")
res_up <- run_enrichment(genes_up, "Upregulated", RESULTS_DIR, TOP_N_PATHWAYS)

cat("2. ORA — Downregulated\n")
res_down <- NULL
if (length(genes_down) >= MIN_GENES_SET) {
  res_down <- run_enrichment(genes_down, "Downregulated", RESULTS_DIR, TOP_N_PATHWAYS)
} else {
  cat("  ⚠ Only", length(genes_down), "genes — skipping\n\n")
}

# CNETPLOT (enrichplot) ------------------------------------------------------
cat("3. Network — cnetplot (GO BP)...\n")

if (!is.null(res_up$go)) {
  tryCatch({
    p_cnet <- cnetplot(res_up$go,
                       showCategory  = 5,
                       color_item    = "steelblue",
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
                          color_item     = "steelblue",
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
cat("4. Network — STRINGdb (full signature)...\n")

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

report_text <- paste0(
  "# Pathways Report — ", OUT_PREFIX, "\n\n",
  "**Signature:** ",
  if (USE_CONSERVATIVE) "conservative (Mean_Expr_Unt > -5)" else "full", "\n\n",
  "## ORA GO BP — Upregulated (top 10)\n",
  extract_top(res_up, "go", 10), "\n\n",
  "## ORA KEGG — Upregulated (top 10)\n",
  extract_top(res_up, "kegg", 10), "\n\n",
  if (!is.null(res_down)) paste0(
    "## ORA GO BP — Downregulated (top 10)\n",
    extract_top(res_down, "go", 10), "\n\n",
    "## ORA KEGG — Downregulated (top 10)\n",
    extract_top(res_down, "kegg", 10), "\n\n"
  ) else "",
  "## Figures\n",
  "- 04_go_bp_upregulated.png\n",
  "- 04_kegg_upregulated.png\n",
  if (!is.null(res_down)) "- 04_go_bp_downregulated.png\n- 04_kegg_downregulated.png\n" else "",
  "- 04_cnetplot_upregulated.png\n",
  if (!is.null(res_down$go)) "- 04_cnetplot_downregulated.png\n" else "",
  "- 04_stringdb_upregulated.png\n"
)
writeLines(report_text, file.path(RESULTS_DIR, "report", "04_pathways_report.md"))
cat("Report: report/04_pathways_report.md\n\n")

cat("================================================================\n")
cat("  PATHWAYS finished —", OUT_PREFIX, "\n")
cat("================================================================\n\n")
