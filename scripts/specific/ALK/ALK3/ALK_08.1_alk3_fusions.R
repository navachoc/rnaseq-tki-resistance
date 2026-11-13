# ============================================================
# ALK_08.1_alk3_fusions.R
# Purpose: RNA fusions (STAR-Fusion) in ALK3, all samples
#   (untreated parental + 2 treated parentals + 8 resistant). Trajectory
#   of the driver fusion EML4-ALK (variant not specified, Pros et al. 2020)
#   and selection of the most striking fusions.
#
# Unlike the rest of the ALK3 scripts (01-05, processed in a single
# "combined" block), this one DOES split by drug in the figures, same
# criterion as rna_08.1_fusions_ret.R splits by model (RET1/RET2): individual
# figures per drug (Lorlatinib/Alectinib) + one combined figure with
# the union of the most striking fusions of both.
#
# Input format: 1 TSV per sample (11 ALK3 samples), same format as
# ALK1/ALK2 (standard STAR-Fusion, NOT the FusionInspector split format
# used by RET). LargeAnchorSupport arrives as "YES_LDAS".
# Derived from ALK_08.1_alk2_fusions.R and rna_08.1_fusions_ret.R (facet
# splitting pattern + combined figure).
#
# Inputs: data/fusions/{Sample}_star-fusion.fusion_predictions.abridged.coding_effect.tsv,
#   data/metadata_ALK3.csv, config/colors.R
# Outputs (results/ALK/ALK3/): tables/08_*.csv, figures/08_*.png,
#   report/08_fusions_report.md
# Usage: source(here::here("scripts/specific/ALK/ALK3/ALK_08.1_alk3_fusions.R"))
#   or Rscript from any folder of the repo.
# ============================================================
# CONFIG
# ============================================================
RECEPTOR      <- "ALK"
CELL_LINE     <- "ALK3"
DATASET       <- "ALK3"
FUSIONS_DIR   <- here::here("data", "fusions")
CONFIG_DIR    <- here::here("config")
OUT_DIR       <- here::here("results", RECEPTOR, CELL_LINE)
TOP_N         <- 10   # most striking fusions to report per drug
TOP_N_HEATMAP <- 3    # fusions per drug in the heatmaps (individual and combined)
DRUGS         <- c("Lorlatinib", "Alectinib")
# Expected driver fusion (to confirm at RNA level) — same fusion in both
# drug arms, since they start from the same ALK3 parental clone
driver_fusion <- "EML4--ALK"
# ============================================================

DRUG <- "Combined"
source(file.path(CONFIG_DIR, "colors.R"))
OUT_PREFIX <- paste0(CELL_LINE, "_", RECEPTOR)

library(dplyr)
library(tidyr)
library(ggplot2)
library(pheatmap)

dir.create(file.path(OUT_DIR, "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT_DIR, "tables"),  recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT_DIR, "report"),  recursive = TRUE, showWarnings = FALSE)

GRADIENTE_RNA <- c("#000000", "#808080", "#f2a900")
heatmap_expr  <- colorRampPalette(GRADIENTE_RNA)(100)

# Own drug colours for ALK3 -- LOCAL override, colors.R is not touched
# (shared with RET, which uses BLU667/LOXO292). Green=Lorlatinib, purple=Alectinib,
# same as GX_cluster (Group2_NF2_YAP/Group1_EGFR_bypass) in MET.
DRUG_COLORS_ALK <- c(Lorlatinib = "#00B894", Alectinib = "#6C5CE7", No = "#DFE6E9")

# Alectinib-specific heatmap substitution: EML4-ALK is already shown in
# the 3 driver trajectory figures (line/bar/dot) -- repeating it in the
# Alectinib top-N heatmap is redundant. It is replaced by LINC02842-ASPH,
# an additional real fusion detected in ALK3_A3 (ABA47810, LargeAnchorSupport=
# YES_LDAS), not present in any Lorlatinib sample.
FUSION_SUSTITUCION <- list(
  Alectinib = list(quitar = "EML4--ALK", poner = "LINC02842--ASPH")
)

cat("\n================================================================\n")
cat("  ALK FUSIONS —", OUT_PREFIX, "(split by drug)\n")
cat("================================================================\n\n")

# ------------------------------------------------------------------
# 1. PARSER — flat STAR-Fusion TSV, 1 row per fusion (27 columns)
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
cat("1. Loading fusions of", CELL_LINE, "(", length(DRUGS), "drugs)...\n")

meta_all <- read.csv(here::here("data", paste0("metadata_", DATASET, ".csv")), sep = ";")
rownames(meta_all) <- meta_all$Sample

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
       "for the ALK3 samples (", paste(meta_all$Sample, collapse = ", "),
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
  left_join(meta_all[, c("Sample", "Model", "Cell_Line", "Condition", "Treatment", "Drug")],
            by = "Sample") %>%
  # From here on "Sample" holds the model name (ALK3,
  # ALK3_L, ALK3_A, ALK3_L2...) instead of the sequencing ID -- the
  # original ID is kept in Sample_ID for traceability.
  rename(Sample_ID = Sample) %>%
  rename(Sample = Model)

cat("  Total fusions (all samples):", nrow(fusions_all), "\n")
cat("  Unique fusions (#FusionName):", length(unique(fusions_all$`#FusionName`)), "\n\n")

write.csv(fusions_all, file.path(OUT_DIR, "tables", "08_fusions_all_samples.csv"), row.names = FALSE)
cat("  → tables/08_fusions_all_samples.csv\n\n")

# ------------------------------------------------------------------
# 3. STAGE ORDER PER DRUG (Untreated -> Treated -> Resistant),
#    same untreated parental shared by both arms
# ------------------------------------------------------------------
build_stage_order <- function(drug) {
  untreated <- meta_all$Model[meta_all$Condition == "Parental" & meta_all$Treatment == "No"]
  treated   <- meta_all$Model[meta_all$Condition == "Parental" & meta_all$Drug == drug]
  resistant <- meta_all$Model[meta_all$Condition == "Resistance" & meta_all$Drug == drug]
  # alphabetical sort() puts L10 before L2 -- numeric order by the suffix
  resistant <- resistant[order(as.numeric(sub(".*?(\\d+)$", "\\1", resistant)))]
  c(untreated, treated, resistant)
}
stage_orders <- setNames(lapply(DRUGS, build_stage_order), DRUGS)

# ------------------------------------------------------------------
# 4. DRIVER FUSION EVOLUTION PER DRUG
#    (fill with FFPM=0 in samples where it is not detected -> real loss)
# ------------------------------------------------------------------
cat("2. Driver fusion trajectory (", driver_fusion, ") per drug...\n")

driver_traj <- bind_rows(lapply(DRUGS, function(drug) {
  samples <- stage_orders[[drug]]
  detected <- fusions_all %>%
    filter(`#FusionName` == driver_fusion) %>%
    group_by(Sample) %>%
    summarise(FFPM = sum(FFPM, na.rm = TRUE),
              JunctionReadCount = sum(JunctionReadCount, na.rm = TRUE),
              .groups = "drop")

  data.frame(Sample = samples) %>%
    left_join(detected[, c("Sample", "FFPM", "JunctionReadCount")], by = "Sample") %>%
    mutate(
      FFPM              = ifelse(is.na(FFPM), 0, FFPM),
      JunctionReadCount = ifelse(is.na(JunctionReadCount), 0, JunctionReadCount),
      Drug       = drug,
      FusionName = driver_fusion,
      Condition  = meta_by_model[Sample, "Condition"],
      Sample     = factor(Sample, levels = samples)
    )
}))

if (all(driver_traj$FFPM == 0)) {
  cat("  ⚠ ", driver_fusion, " not detected in any sample — check the exact fusion name",
      "in 08_fusions_all_samples.csv\n\n")
}

write.csv(driver_traj, file.path(OUT_DIR, "tables", "08_fusions_driver_trajectory.csv"), row.names = FALSE)

driver_theme <- theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    plot.title  = element_text(face = "bold"),
    legend.position = "none",
    strip.text = element_text(face = "bold"),
    panel.background = element_rect(fill = "white"),
    plot.background  = element_rect(fill = "white")
  )
driver_labs <- labs(
  subtitle = paste0("ALK3: ", driver_fusion, " (same fusion in both drug branches)"),
  x = NULL, y = "FFPM (fusion fragments per million)"
)

p_driver_line <- ggplot(driver_traj, aes(x = Sample, y = FFPM, group = Drug, color = Drug)) +
  geom_line(linewidth = 0.9, alpha = 0.85) +
  geom_point(size = 2.5) +
  facet_wrap(~ Drug, scales = "free_x", ncol = 1) +
  scale_color_manual(values = DRUG_COLORS_ALK) +
  labs(title = paste("EML4-ALK driver fusion — evolution Untreated -> Treated -> Resistant (lineplot) —", CELL_LINE)) +
  driver_labs + driver_theme

ggsave(file.path(OUT_DIR, "figures", "08_fusions_driver_lineplot.png"),
       p_driver_line, width = 11, height = 8, dpi = 300)
cat("  → figures/08_fusions_driver_lineplot.png\n")

p_driver_bar <- ggplot(driver_traj, aes(x = Sample, y = FFPM, fill = Drug)) +
  geom_col(alpha = 0.9) +
  facet_wrap(~ Drug, scales = "free_x", ncol = 1) +
  scale_fill_manual(values = DRUG_COLORS_ALK) +
  labs(title = paste("EML4-ALK driver fusion — evolution Untreated -> Treated -> Resistant (barplot) —", CELL_LINE)) +
  driver_labs + driver_theme

ggsave(file.path(OUT_DIR, "figures", "08_fusions_driver_barplot.png"),
       p_driver_bar, width = 11, height = 8, dpi = 300)
cat("  → figures/08_fusions_driver_barplot.png\n")

p_driver_dot <- ggplot(driver_traj, aes(x = Sample, y = 1,
                                          size = JunctionReadCount, color = FFPM)) +
  geom_point() +
  facet_wrap(~ Drug, scales = "free_x", ncol = 1) +
  scale_color_gradientn(colors = GRADIENTE_RNA) +
  scale_size_continuous(range = c(1, 12)) +
  labs(
    title = paste("EML4-ALK driver fusion — evolution Untreated -> Treated -> Resistant (dotplot) —", CELL_LINE),
    subtitle = paste0("ALK3: ", driver_fusion),
    x = NULL, y = NULL, color = "FFPM", size = "Junction\nreads"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.title  = element_text(face = "bold"),
    strip.text = element_text(face = "bold"),
    panel.grid.major.y = element_blank(),
    panel.background = element_rect(fill = "white"),
    plot.background  = element_rect(fill = "white")
  )

ggsave(file.path(OUT_DIR, "figures", "08_fusions_driver_dotplot.png"),
       p_driver_dot, width = 11, height = 5.5, dpi = 300)
cat("  → figures/08_fusions_driver_dotplot.png\n\n")

# ------------------------------------------------------------------
# 5. FUSION HEATMAP PER DRUG (top N by FFPM, stage order)
#    + combined heatmap with the union of the top N of each drug
# ------------------------------------------------------------------
cat("3. Fusion heatmaps per drug (top", TOP_N_HEATMAP, "FFPM)...\n")

fusions_agg <- fusions_all %>%
  group_by(Sample, `#FusionName`) %>%
  summarise(FFPM = sum(FFPM, na.rm = TRUE), .groups = "drop")

get_top_fusions <- function(drug, n = TOP_N_HEATMAP) {
  ranked <- fusions_agg %>%
    filter(Sample %in% stage_orders[[drug]]) %>%
    group_by(`#FusionName`) %>%
    summarise(max_FFPM = max(FFPM, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(max_FFPM)) %>%
    pull(`#FusionName`)

  # See FUSION_SUSTITUCION above -- removes a fusion redundant with the
  # driver trajectory and forces another one in its place
  sust <- FUSION_SUSTITUCION[[drug]]
  if (!is.null(sust)) {
    ranked <- setdiff(ranked, sust$quitar)
    ranked <- unique(c(sust$poner, ranked))
  }
  head(ranked, n)
}

build_fusion_matrix <- function(samples, top_fusions) {
  grid <- expand.grid(Sample = samples, FusionName = top_fusions, stringsAsFactors = FALSE) %>%
    left_join(fusions_agg, by = c("Sample" = "Sample", "FusionName" = "#FusionName")) %>%
    mutate(FFPM = ifelse(is.na(FFPM), 0, FFPM))
  mat <- grid %>%
    pivot_wider(names_from = Sample, values_from = FFPM) %>%
    tibble::column_to_rownames("FusionName") %>%
    as.matrix()
  log2(mat[, samples, drop = FALSE] + 1)
}

plot_fusion_heatmap <- function(drug, n = TOP_N_HEATMAP) {
  samples     <- stage_orders[[drug]]
  top_fusions <- get_top_fusions(drug, n)

  if (length(top_fusions) < 2) {
    cat("  ⚠", drug, "— fewer than 2 fusions, heatmap skipped\n")
    return(invisible(NULL))
  }

  mat_log <- build_fusion_matrix(samples, top_fusions)
  # Treatment/Condition -- Condition last so it ends up at the top of the bar
  ann_col <- meta_by_model[samples, c("Treatment", "Condition"), drop = FALSE]
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
           main = paste0("RNA fusions (top ", n, " FFPM) — ", drug, " — ", CELL_LINE),
           filename = file.path(OUT_DIR, "figures", paste0("08_fusions_heatmap_", drug, ".png")),
           width = 9, height = max(4, length(top_fusions) * 0.5))
  cat("  →", paste0("figures/08_fusions_heatmap_", drug, ".png"), "\n")
}
invisible(lapply(DRUGS, plot_fusion_heatmap))

# Combined heatmap: union of the top N of Lorlatinib and the top N of Alectinib,
# all samples of both drugs, stage order per drug
top_fusions_combined <- Reduce(union, lapply(DRUGS, get_top_fusions))
combined_samples <- unlist(stage_orders[DRUGS], use.names = FALSE)
combined_samples <- combined_samples[!duplicated(combined_samples)]

if (length(top_fusions_combined) >= 2) {
  mat_log_combined <- build_fusion_matrix(combined_samples, top_fusions_combined)
  # Drug/Treatment/Condition -- Condition last so it ends up at the top of the bar
  ann_col_combined <- meta_by_model[combined_samples, c("Drug", "Treatment", "Condition"), drop = FALSE]
  ann_colors_combined <- filter_annotation_colors(ann_col_combined)
  ann_colors_combined$Drug <- DRUG_COLORS_ALK

  combined_fig <- paste0("08_fusions_heatmap_combined_top", TOP_N_HEATMAP, ".png")
  pheatmap(mat_log_combined,
           annotation_col = ann_col_combined,
           annotation_colors = ann_colors_combined,
           color = heatmap_expr,
           cluster_rows = TRUE,
           cluster_cols = FALSE,
           show_rownames = TRUE,
           show_colnames = TRUE,
           fontsize_row = 9,
           fontsize_col = 8,
           main = paste0("RNA fusions — union top ", TOP_N_HEATMAP, " per drug (Lorlatinib + Alectinib) — ", CELL_LINE),
           filename = file.path(OUT_DIR, "figures", combined_fig),
           width = 12, height = max(4, length(top_fusions_combined) * 0.5))
  cat("  →", file.path("figures", combined_fig), "(", length(top_fusions_combined), "fusions )\n\n")
}

# ------------------------------------------------------------------
# 6. MOST STRIKING FUSIONS PER DRUG
#    Criterion: robust support (LargeAnchorSupport=YES_LDAS), high FFPM,
#    evidence in known databases (ChimerKB/COSMIC/Mitelman in 'annots')
# ------------------------------------------------------------------
cat("4. Most striking fusions per drug...\n")

top_striking <- fusions_all %>%
  mutate(known_db = grepl("ChimerKB|Cosmic|COSMIC|Mitelman|ChimerSeq|FusionGDB", annots, ignore.case = TRUE)) %>%
  group_by(Drug, `#FusionName`) %>%
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
  filter(Drug %in% DRUGS) %>%
  mutate(is_driver = `#FusionName` == driver_fusion) %>%
  filter(any_LargeAnchor | is_driver) %>%
  arrange(Drug, desc(is_driver), desc(mean_FFPM)) %>%
  group_by(Drug) %>%
  slice_head(n = TOP_N) %>%
  ungroup()

write.csv(top_striking, file.path(OUT_DIR, "tables", "08_fusions_top_striking.csv"), row.names = FALSE)
cat("  → tables/08_fusions_top_striking.csv (top", TOP_N, "per drug)\n\n")

# ------------------------------------------------------------------
# 7. COINCIDENCES BETWEEN DRUGS
# ------------------------------------------------------------------
cat("5. Fusions shared between Lorlatinib and Alectinib...\n")

by_drug <- fusions_all %>%
  filter(Drug %in% DRUGS) %>%
  group_by(`#FusionName`, Drug) %>%
  summarise(n_samples = n_distinct(Sample), mean_FFPM = round(mean(FFPM, na.rm = TRUE), 2), .groups = "drop")

shared <- by_drug %>%
  group_by(`#FusionName`) %>%
  filter(n_distinct(Drug) == length(DRUGS)) %>%
  ungroup() %>%
  pivot_wider(names_from = Drug, values_from = c(n_samples, mean_FFPM)) %>%
  arrange(desc(pmax(!!!rlang::syms(paste0("mean_FFPM_", DRUGS)))))

write.csv(shared, file.path(OUT_DIR, "tables", "08_fusions_shared_drugs.csv"), row.names = FALSE)
cat("  → tables/08_fusions_shared_drugs.csv (", nrow(shared), "shared fusions)\n\n")

# ------------------------------------------------------------------
# 8. REPORT
# ------------------------------------------------------------------
driver_check <- driver_traj %>%
  group_by(Drug) %>%
  summarise(min_FFPM = round(min(FFPM), 2), max_FFPM = round(max(FFPM), 2), .groups = "drop")

report_text <- paste0(
  "# RNA fusions (STAR-Fusion) — ALK: ALK3 (Lorlatinib + Alectinib)\n\n",
  "## Source\n",
  "- ", length(fusion_files), " tsv files in `data/fusions/` (STAR-Fusion, abridged.coding_effect output)\n",
  "- Samples identified by Model (ALK3, ALK3_L, ALK3_A, ALK3_L2-L10, ALK3_A2-A4); original sequencing ID in `Sample_ID`\n",
  "- ", nrow(fusions_all), " total fusion events, ", length(unique(fusions_all$`#FusionName`)), " unique fusions\n\n",
  "## Driver fusion — confirmation at RNA level, per drug\n",
  paste(capture.output(print(as.data.frame(driver_check), row.names = FALSE)), collapse = "\n"), "\n",
  "- Same fusion expected in both branches (same ALK3 parental clone): ", driver_fusion, "\n\n",
  "## Most striking fusions per drug (top ", TOP_N, ")\n",
  "See `tables/08_fusions_top_striking.csv`. Criterion: LargeAnchorSupport=YES_LDAS and/or driver fusion, ",
  "ordered by mean FFPM; `known_db` is flagged if present in ChimerKB/COSMIC/Mitelman.\n\n",
  "## Overlaps between drugs\n",
  "See `tables/08_fusions_shared_drugs.csv` — ", nrow(shared), " fusions detected in both branches (Lorlatinib and Alectinib).\n\n",
  "## Interpretation\n",
  "[INTERPRET]\n",
  "- Is the driver fusion maintained in all 8 resistant lines or lost in any (driver switch)?\n",
  "- Are the fusions shared between drugs technical noise (read-through) or biological candidates?\n",
  "- Does any new fusion appear in only one drug (acquired, specific to that inhibitor generation)?\n\n",
  "## Figures\n",
  "- 08_fusions_driver_lineplot.png (faceted by drug)\n",
  "- 08_fusions_driver_barplot.png (faceted by drug)\n",
  "- 08_fusions_driver_dotplot.png (size=JunctionReadCount, color=FFPM, faceted by drug)\n",
  "- 08_fusions_heatmap_Lorlatinib.png (top ", TOP_N_HEATMAP, ")\n",
  "- 08_fusions_heatmap_Alectinib.png (top ", TOP_N_HEATMAP, ")\n",
  "- 08_fusions_heatmap_combined_top", TOP_N_HEATMAP, ".png (union Lorlatinib + Alectinib)\n"
)
writeLines(report_text, file.path(OUT_DIR, "report", "08_fusions_report.md"))
cat("Report: ALK/ALK3/report/08_fusions_report.md\n\n")

cat("================================================================\n")
cat("  ALK FUSIONS completed —", OUT_PREFIX, "\n")
cat("================================================================\n\n")
