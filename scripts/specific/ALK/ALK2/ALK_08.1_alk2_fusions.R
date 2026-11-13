# ============================================================
# ALK_08.1_alk2_fusions.R
# Purpose: RNA fusions (STAR-Fusion) in ALK2, all samples
#   (untreated/treated parental + resistant L1-L5). Trajectory of the driver
#   fusion EML4-ALK v1 and selection of the most striking fusions.
#   Single-model dedicated version, derived from
#   scripts/simple/rna_08.1_fusions_ret.R (there combined RET1+RET2).
#
# Input format: 1 TSV per sample, standard STAR-Fusion output
# (1 row per fusion, 27 columns), WITHOUT the 2-row split format of
# FusionInspector used by rna_08.1_fusions_ret.R for RET. There is also an
# ALK2_all_fusions.xlsx with one sheet per sample (same content as the
# TSVs); the TSVs are used because they are simpler to read.
# Note: LargeAnchorSupport arrives as "YES_LDAS", not "YES" as in RET,
# hence grepl("YES", ...) is used instead of an exact comparison.
#
# Inputs: data/fusions/{Sample}_star-fusion.fusion_predictions.abridged.coding_effect.tsv,
#   data/metadata_ALK12.csv, config/colors.R
# Outputs (results/ALK/ALK2/): tables/08_fusions_all_samples.csv,
#   tables/08_fusions_driver_trajectory.csv, tables/08_fusions_top_striking.csv,
#   figures/08_fusions_*.png, report/08_fusions_report.md
# Usage: source(here::here("scripts/specific/ALK/ALK2/ALK_08.1_alk2_fusions.R"))
#   or Rscript from any folder of the repo.
# ============================================================
# CONFIG
# ============================================================
RECEPTOR      <- "ALK"
CELL_LINE     <- "ALK2"
FUSIONS_DIR   <- here::here("data", "fusions")
CONFIG_DIR    <- here::here("config")
OUT_DIR       <- here::here("results", RECEPTOR, CELL_LINE)
TOP_N         <- 10   # most striking fusions to report
TOP_N_HEATMAP <- 3    # fusions in the heatmap
# Expected driver fusion (to confirm at RNA level) — EML4-ALK variant 1
driver_fusion <- "EML4--ALK"
# ============================================================

DRUG <- "Lorlatinib"
source(file.path(CONFIG_DIR, "colors.R"))
OUT_PREFIX <- paste0(CELL_LINE, "_", RECEPTOR)

library(dplyr)
library(tidyr)
library(ggplot2)
library(pheatmap)

dir.create(file.path(OUT_DIR, "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT_DIR, "tables"),  recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT_DIR, "report"),  recursive = TRUE, showWarnings = FALSE)

# RNA-seq gradient (low -> high: black -> grey -> orange/gold)
GRADIENTE_RNA <- c("#000000", "#808080", "#f2a900")
heatmap_expr  <- colorRampPalette(GRADIENTE_RNA)(100)

cat("\n================================================================\n")
cat("  ALK FUSIONS —", OUT_PREFIX, "\n")
cat("================================================================\n\n")

# ------------------------------------------------------------------
# 1. PARSER — flat STAR-Fusion TSV, 1 row per fusion (27 columns,
#    without the FusionInspector continuation as in RET)
# ------------------------------------------------------------------
numeric_cols <- c("JunctionReadCount", "SpanningFragCount", "est_J", "est_S",
                   "FFPM", "LeftBreakEntropy", "RightBreakEntropy")

read_fusion_file <- function(path) {
  out <- read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  if (nrow(out) == 0) return(NULL)
  out[numeric_cols] <- lapply(out[numeric_cols], function(x) suppressWarnings(as.numeric(x)))
  out
}

# ------------------------------------------------------------------
# 2. LOAD ALL SAMPLES + METADATA
# ------------------------------------------------------------------
cat("1. Loading fusions of", CELL_LINE, "...\n")

meta_all <- read.csv(here::here("data", paste0("metadata_", "ALK12", ".csv")), sep = ";") %>%
  filter(Cell_Line == CELL_LINE)
rownames(meta_all) <- meta_all$Sample

# All the figures/tables of this script are labelled by Model (ALK2, ALK2_L,
# ALK2_L1...), not by the sequencing ID (ABA...) -- meta_by_model allows
# indexing metadata rows directly by Model.
meta_by_model <- meta_all
rownames(meta_by_model) <- meta_by_model$Model

fusion_files <- list.files(
  FUSIONS_DIR,
  pattern = "_star-fusion\\.fusion_predictions\\.abridged\\.coding_effect\\.tsv$",
  full.names = TRUE
)
names(fusion_files) <- sub("_star-fusion\\.fusion_predictions\\.abridged\\.coding_effect\\.tsv$",
                            "", basename(fusion_files))
fusion_files <- fusion_files[names(fusion_files) %in% meta_all$Sample]
cat("  Fusion files found:", length(fusion_files), "of", nrow(meta_all), "samples in metadata\n\n")

if (length(fusion_files) == 0) {
  stop("No files data/fusions/*_star-fusion.fusion_predictions.abridged.coding_effect.tsv ",
       "for the ALK2 samples (", paste(meta_all$Sample, collapse = ", "),
       "). Place them in data/fusions/ before running this script.")
}

fusions_all <- bind_rows(lapply(names(fusion_files), function(sample) {
  df <- read_fusion_file(fusion_files[[sample]])
  if (is.null(df) || nrow(df) == 0) return(NULL)
  df$Sample <- sample
  df
}))

fusions_all <- fusions_all %>%
  mutate(
    LeftGeneSymbol  = sub("\\^.*", "", LeftGene),
    RightGeneSymbol = sub("\\^.*", "", RightGene)
  ) %>%
  left_join(meta_all[, c("Sample", "Model", "Cell_Line", "Condition", "Treatment")],
            by = "Sample") %>%
  # From here on "Sample" holds the model name (ALK2,
  # ALK2_L, ALK2_L1...) instead of the sequencing ID -- the
  # original ID is kept in Sample_ID for traceability.
  rename(Sample_ID = Sample) %>%
  rename(Sample = Model)

cat("  Total fusions (all samples):", nrow(fusions_all), "\n")
cat("  Unique fusions (#FusionName):", length(unique(fusions_all$`#FusionName`)), "\n\n")

write.csv(fusions_all, file.path(OUT_DIR, "tables", "08_fusions_all_samples.csv"), row.names = FALSE)
cat("  → tables/08_fusions_all_samples.csv\n\n")

# ------------------------------------------------------------------
# 3. STAGE ORDER (Untreated -> Treated -> Resistant)
# ------------------------------------------------------------------
# In Model, not in Sample (sequencing ID) -- see the note in the loading of fusions_all
sample_order <- c(
  meta_all$Model[meta_all$Condition == "Parental" & meta_all$Treatment == "No"],
  meta_all$Model[meta_all$Condition == "Parental" & meta_all$Treatment == "Yes"],
  sort(meta_all$Model[meta_all$Condition == "Resistance"])
)

# ------------------------------------------------------------------
# 4. DRIVER FUSION EVOLUTION (fill with FFPM=0 if not detected
#    in a sample -> real loss of the fusion transcript)
# ------------------------------------------------------------------
cat("2. Driver fusion trajectory (", driver_fusion, ")...\n")

detected <- fusions_all %>%
  filter(`#FusionName` == driver_fusion) %>%
  group_by(Sample) %>%
  summarise(FFPM = sum(FFPM, na.rm = TRUE),
            JunctionReadCount = sum(JunctionReadCount, na.rm = TRUE),
            .groups = "drop")

driver_traj <- data.frame(Sample = sample_order) %>%
  left_join(detected[, c("Sample", "FFPM", "JunctionReadCount")], by = "Sample") %>%
  mutate(
    FFPM              = ifelse(is.na(FFPM), 0, FFPM),
    JunctionReadCount = ifelse(is.na(JunctionReadCount), 0, JunctionReadCount),
    FusionName = driver_fusion,
    Condition  = meta_by_model[as.character(Sample), "Condition"],
    Sample = factor(Sample, levels = sample_order)
  )

if (all(driver_traj$FFPM == 0)) {
  cat("  ⚠ ", driver_fusion, " not detected in any sample — check the exact fusion name",
      "in 08_fusions_all_samples.csv (variant v1 may be reported with another breakpoint/gene order)\n\n")
}

write.csv(driver_traj, file.path(OUT_DIR, "tables", "08_fusions_driver_trajectory.csv"), row.names = FALSE)

driver_theme <- theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    plot.title  = element_text(face = "bold"),
    legend.position = "none",
    panel.background = element_rect(fill = "white"),
    plot.background  = element_rect(fill = "white")
  )
driver_labs <- labs(
  subtitle = paste0(CELL_LINE, ": ", driver_fusion),
  x = NULL, y = "FFPM (fusion fragments per million)"
)

p_driver_line <- ggplot(driver_traj, aes(x = Sample, y = FFPM, group = 1)) +
  geom_line(linewidth = 0.9, alpha = 0.85, color = cellline_colors[[CELL_LINE]]) +
  geom_point(size = 2.5, color = cellline_colors[[CELL_LINE]]) +
  labs(title = paste("EML4-ALK v1 driver fusion — evolution Untreated -> Treated -> Resistant (lineplot) —", CELL_LINE)) +
  driver_labs + driver_theme

ggsave(file.path(OUT_DIR, "figures", "08_fusions_driver_lineplot.png"),
       p_driver_line, width = 9, height = 5.5, dpi = 300)
cat("  → figures/08_fusions_driver_lineplot.png\n")

p_driver_dot <- ggplot(driver_traj, aes(x = Sample, y = 1,
                                          size = JunctionReadCount, color = FFPM)) +
  geom_point() +
  scale_color_gradientn(colors = GRADIENTE_RNA) +
  scale_size_continuous(range = c(1, 12)) +
  labs(
    title = paste("EML4-ALK v1 driver fusion — evolution Untreated -> Treated -> Resistant (dotplot) —", CELL_LINE),
    subtitle = paste0(CELL_LINE, ": ", driver_fusion),
    x = NULL, y = NULL, color = "FFPM", size = "Junction\nreads"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.title  = element_text(face = "bold"),
    panel.grid.major.y = element_blank(),
    panel.background = element_rect(fill = "white"),
    plot.background  = element_rect(fill = "white")
  )

ggsave(file.path(OUT_DIR, "figures", "08_fusions_driver_dotplot.png"),
       p_driver_dot, width = 9, height = 4.5, dpi = 300)
cat("  → figures/08_fusions_driver_dotplot.png\n\n")

# ------------------------------------------------------------------
# 5. FUSION HEATMAP (top N by FFPM, stage order)
# ------------------------------------------------------------------
cat("3. Fusion heatmap (top", TOP_N_HEATMAP, "FFPM)...\n")

fusions_agg <- fusions_all %>%
  group_by(Sample, `#FusionName`) %>%
  summarise(FFPM = sum(FFPM, na.rm = TRUE), .groups = "drop")

top_fusions <- fusions_agg %>%
  filter(Sample %in% sample_order) %>%
  group_by(`#FusionName`) %>%
  summarise(max_FFPM = max(FFPM, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(max_FFPM)) %>%
  slice_head(n = TOP_N_HEATMAP) %>%
  pull(`#FusionName`)

if (length(top_fusions) >= 2) {
  grid <- expand.grid(Sample = sample_order, FusionName = top_fusions, stringsAsFactors = FALSE) %>%
    left_join(fusions_agg, by = c("Sample" = "Sample", "FusionName" = "#FusionName")) %>%
    mutate(FFPM = ifelse(is.na(FFPM), 0, FFPM))
  mat_log <- grid %>%
    pivot_wider(names_from = Sample, values_from = FFPM) %>%
    tibble::column_to_rownames("FusionName") %>%
    as.matrix()
  mat_log <- log2(mat_log[, sample_order, drop = FALSE] + 1)

  ann_col <- meta_by_model[sample_order, c("Condition", "Treatment"), drop = FALSE]
  ann_colors <- filter_annotation_colors(ann_col)

  pheatmap(mat_log,
           annotation_col = ann_col,
           annotation_colors = ann_colors,
           color = heatmap_expr,
           cluster_rows = TRUE,
           cluster_cols = FALSE,
           show_rownames = TRUE,
           show_colnames = TRUE,
           fontsize_row = 9,
           main = paste0("RNA fusions (top ", TOP_N_HEATMAP, " FFPM) — ", CELL_LINE),
           filename = file.path(OUT_DIR, "figures", paste0("08_fusions_heatmap_", CELL_LINE, ".png")),
           width = 9, height = max(4, length(top_fusions) * 0.5))
  cat("  →", paste0("figures/08_fusions_heatmap_", CELL_LINE, ".png"), "\n\n")
} else {
  cat("  ⚠ Fewer than 2 fusions — heatmap skipped\n\n")
}

# ------------------------------------------------------------------
# 6. MOST STRIKING FUSIONS
#    Criterion: robust support (LargeAnchorSupport=YES), high FFPM,
#    evidence in known databases (ChimerKB/COSMIC/Mitelman in 'annots')
# ------------------------------------------------------------------
cat("4. Most striking fusions...\n")

top_striking <- fusions_all %>%
  mutate(known_db = grepl("ChimerKB|Cosmic|COSMIC|Mitelman|ChimerSeq|FusionGDB", annots, ignore.case = TRUE)) %>%
  group_by(`#FusionName`) %>%
  summarise(
    n_samples          = n_distinct(Sample),
    samples            = paste(sort(unique(Sample)), collapse = ", "),
    mean_FFPM          = round(mean(FFPM, na.rm = TRUE), 2),
    max_FFPM           = round(max(FFPM, na.rm = TRUE), 2),
    max_JunctionReads  = max(JunctionReadCount, na.rm = TRUE),
    any_LargeAnchor    = any(grepl("YES", LargeAnchorSupport), na.rm = TRUE),
    prot_effect        = paste(unique(na.omit(PROT_FUSION_TYPE)), collapse = "/"),
    known_db           = any(known_db),
    .groups = "drop"
  ) %>%
  mutate(is_driver = `#FusionName` == driver_fusion) %>%
  filter(any_LargeAnchor | is_driver) %>%
  arrange(desc(is_driver), desc(mean_FFPM)) %>%
  slice_head(n = TOP_N)

write.csv(top_striking, file.path(OUT_DIR, "tables", "08_fusions_top_striking.csv"), row.names = FALSE)
cat("  → tables/08_fusions_top_striking.csv (top", TOP_N, ")\n\n")

# ------------------------------------------------------------------
# 7. REPORT
# ------------------------------------------------------------------
driver_check <- driver_traj %>%
  summarise(min_FFPM = round(min(FFPM), 2), max_FFPM = round(max(FFPM), 2))

report_text <- paste0(
  "# RNA fusions (STAR-Fusion) — ALK: ", CELL_LINE, "\n\n",
  "## Source\n",
  "- ", length(fusion_files), " tsv files in `data/fusions/` (STAR-Fusion, abridged.coding_effect output)\n",
  "- Samples identified by Model (ALK2, ALK2_L, ALK2_L1-L5); original sequencing ID in `Sample_ID`\n",
  "- ", nrow(fusions_all), " total fusion events, ", length(unique(fusions_all$`#FusionName`)), " unique fusions\n\n",
  "## Driver fusion — confirmation at RNA level\n",
  "- ", CELL_LINE, " (EML4-ALK v1): FFPM range [", driver_check$min_FFPM, ", ", driver_check$max_FFPM, "]\n\n",
  "## Most striking fusions (top ", TOP_N, ")\n",
  "See `tables/08_fusions_top_striking.csv`. Criterion: LargeAnchorSupport=YES and/or driver fusion, ",
  "ordered by mean FFPM; `known_db` is flagged if present in ChimerKB/COSMIC/Mitelman.\n\n",
  "## Interpretation\n",
  "[INTERPRET]\n",
  "- Is the driver fusion maintained in all 5 resistant lines (L1-L5) or lost in any (driver switch)?\n",
  "- Does loss/maintenance of the fusion coincide with the IC50 group (L1/L4/L5 shifted vs L2/L3 unshifted)?\n",
  "- Does any new fusion appear only in resistant lines (acquired) vs. already present in the parental?\n\n",
  "## Figures\n",
  "- 08_fusions_driver_lineplot.png\n",
  "- 08_fusions_driver_dotplot.png (size=JunctionReadCount, color=FFPM, black-grey-gold gradient)\n",
  "- 08_fusions_heatmap_", CELL_LINE, ".png (top ", TOP_N_HEATMAP, "; black-grey-gold gradient)\n"
)
writeLines(report_text, file.path(OUT_DIR, "report", "08_fusions_report.md"))
cat("Report: ALK/", CELL_LINE, "/report/08_fusions_report.md\n\n", sep = "")

cat("================================================================\n")
cat("  ALK FUSIONS completed —", OUT_PREFIX, "\n")
cat("================================================================\n\n")
