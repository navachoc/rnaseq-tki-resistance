# ============================================================
# rna_08.1_fusions_ret.R
# RNA fusions (STAR-Fusion): RET, RET1 + RET2, all samples
# Evolution per model/drug, most striking fusions per model, and overlaps
# between both models.
# Source: data/fusions/*.xlsx (STAR-Fusion + FusionInspector output,
# 1 file per sample, format split into 2 rows per fusion).
# RNAseq_bulk pipeline, RET
#
# Inputs: data/fusions/*.xlsx, data/metadata_{RET1,RET2}.csv;
#   config/colors.R
# Outputs (results/RET/combined/): tables/08_fusions_all_samples.csv,
#   tables/08_fusions_driver_trajectory.csv, tables/08_fusions_top_striking.csv,
#   tables/08_fusions_shared_models.csv, figures/08_fusions_*.png,
#   report/08_fusions_report.md
# Dependencies: none within the pipeline (reads the STAR-Fusion output).
#   Its table 08_fusions_all_samples.csv is read by RET_09_lfc_heatmaps.R.
# Usage: paths are resolved with here::here() from the repository root
# (.here file), so it works from any working directory:
#   Rscript scripts/specific/RET/rna_08.1_fusions_ret.R
# ============================================================
# CONFIG
# ============================================================
RECEPTOR    <- "RET"
CELL_LINES  <- c("RET1", "RET2")
FUSIONS_DIR <- here::here("data", "fusions")
CONFIG_DIR  <- here::here("config")
OUT_DIR     <- here::here("results", RECEPTOR, "combined")
TOP_N       <- 10   # most striking fusions to report per model
TOP_N_HEATMAP <- 3  # fusions per model in the heatmaps (individual and combined)
# Expected driver fusions (to be confirmed at the RNA level)
driver_fusion <- c(RET1 = "CCDC6--RET", RET2 = "KIF5B--RET")
# ============================================================

DRUG      <- RECEPTOR
CELL_LINE <- paste(CELL_LINES, collapse = "_")
source(file.path(CONFIG_DIR, "colors.R"))
OUT_PREFIX <- paste0(CELL_LINE, "_", RECEPTOR)

library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(pheatmap)

dir.create(file.path(OUT_DIR, "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT_DIR, "tables"),  recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT_DIR, "report"),  recursive = TRUE, showWarnings = FALSE)

cat("\n================================================================\n")
cat("  RET FUSIONS —", OUT_PREFIX, "\n")
cat("================================================================\n\n")

# ------------------------------------------------------------------
# 1. PARSER: each xlsx has 1 fusion per 2 physical rows:
#    main row (24 STAR-Fusion cols) + continuation row
#    (11 protein-effect annotation cols, in cols 2-12)
# ------------------------------------------------------------------
numeric_cols <- c("JunctionReadCount", "SpanningFragCount", "est_J", "est_S",
                   "FAR_left", "FAR_right", "LeftBreakEntropy", "RightBreakEntropy",
                   "FFPM", "microh_brkpt_dist", "num_microh_near_brkpt")

read_fusion_file <- function(path) {
  raw <- suppressMessages(read_excel(path, col_names = FALSE, .name_repair = "minimal"))
  main_header <- as.character(raw[2, 1:24])
  cont_header <- as.character(raw[3, 2:12])
  data <- raw[4:nrow(raw), ]

  main_idx <- which(!is.na(data[[1]]))
  if (length(main_idx) == 0) return(NULL)

  main_df <- as.data.frame(data[main_idx, 1:24])
  colnames(main_df) <- main_header

  cont_idx <- pmin(main_idx + 1, nrow(data))
  cont_df  <- as.data.frame(data[cont_idx, 2:12])
  colnames(cont_df) <- cont_header

  out <- cbind(main_df, cont_df)
  out[numeric_cols] <- lapply(out[numeric_cols], function(x) suppressWarnings(as.numeric(x)))
  out
}

# ------------------------------------------------------------------
# 2. LOAD ALL SAMPLES + METADATA
# ------------------------------------------------------------------
cat("1. Loading fusions from", length(CELL_LINES), "cell lines...\n")

meta_all <- do.call(rbind, lapply(CELL_LINES, function(cl) {
  read.csv(here::here("data", paste0("metadata_", cl, ".csv")), sep = ";")
}))
rownames(meta_all) <- meta_all$Sample

fusion_files <- list.files(FUSIONS_DIR, pattern = "\\.xlsx$", full.names = TRUE)
names(fusion_files) <- tools::file_path_sans_ext(basename(fusion_files))
fusion_files <- fusion_files[names(fusion_files) %in% meta_all$Sample]
cat("  Fusion files found:", length(fusion_files), "of", nrow(meta_all), "samples in metadata\n\n")

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
            by = "Sample")

cat("  Total fusions (all samples):", nrow(fusions_all), "\n")
cat("  Unique fusions (#FusionName):", length(unique(fusions_all$`#FusionName`)), "\n\n")

write.csv(fusions_all, file.path(OUT_DIR, "tables", "08_fusions_all_samples.csv"), row.names = FALSE)
cat("  → tables/08_fusions_all_samples.csv\n\n")

# ------------------------------------------------------------------
# 3. STAGE ORDER PER MODEL/DRUG (Untreated -> Treated -> Resistant)
#    same criterion as rna_03b (untreated parental -> acutely treated -> resistant)
# ------------------------------------------------------------------
build_stage_order <- function(cl) {
  m <- meta_all[meta_all$Cell_Line == cl, ]
  untreated <- m$Sample[m$Condition == "Parental" & m$Treatment == "No"]
  drugs     <- setdiff(unique(m$Drug), "No")
  ordered   <- untreated
  for (d in drugs) {
    treated   <- m$Sample[m$Condition == "Parental" & m$Drug == d]
    resistant <- sort(m$Sample[m$Condition == "Resistance" & m$Drug == d])
    ordered   <- c(ordered, treated, resistant)
  }
  ordered
}
stage_orders <- setNames(lapply(CELL_LINES, build_stage_order), CELL_LINES)

# ------------------------------------------------------------------
# 4. DRIVER FUSION EVOLUTION PER MODEL/DRUG
#    (fill with FFPM=0 in samples where it is not detected -> real loss)
# ------------------------------------------------------------------
cat("2. Driver fusion trajectory per model...\n")

driver_traj <- bind_rows(lapply(CELL_LINES, function(cl) {
  samples <- stage_orders[[cl]]
  fname   <- driver_fusion[[cl]]
  # Sum in case the same gene pair is reported with several breakpoints
  detected <- fusions_all %>%
    filter(Cell_Line == cl, `#FusionName` == fname) %>%
    group_by(Sample) %>%
    summarise(FFPM = sum(FFPM, na.rm = TRUE),
              JunctionReadCount = sum(JunctionReadCount, na.rm = TRUE),
              .groups = "drop")

  data.frame(Sample = samples) %>%
    left_join(detected[, c("Sample", "FFPM", "JunctionReadCount")], by = "Sample") %>%
    mutate(
      FFPM              = ifelse(is.na(FFPM), 0, FFPM),
      JunctionReadCount = ifelse(is.na(JunctionReadCount), 0, JunctionReadCount),
      Cell_Line = cl,
      FusionName = fname,
      Sample = factor(Sample, levels = samples)
    ) %>%
    left_join(meta_all[, c("Sample", "Drug", "Condition")], by = "Sample")
}))

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
  subtitle = paste0("RET1: ", driver_fusion["RET1"], "  |  RET2: ", driver_fusion["RET2"]),
  x = NULL, y = "FFPM (fusion fragments per million)"
)

p_driver_line <- ggplot(driver_traj, aes(x = Sample, y = FFPM, group = Cell_Line, color = Cell_Line)) +
  geom_line(linewidth = 0.9, alpha = 0.85) +
  geom_point(size = 2.5) +
  facet_wrap(~ Cell_Line, scales = "free_x", ncol = 1) +
  scale_color_manual(values = cellline_colors) +
  labs(title = "RET driver fusion — evolution Untreated -> Treated -> Resistant (lineplot)") +
  driver_labs + driver_theme

ggsave(file.path(OUT_DIR, "figures", "08_fusions_driver_lineplot.png"),
       p_driver_line, width = 11, height = 8, dpi = 300)
cat("  → figures/08_fusions_driver_lineplot.png\n")

p_driver_bar <- ggplot(driver_traj, aes(x = Sample, y = FFPM, fill = Cell_Line)) +
  geom_col(alpha = 0.9) +
  facet_wrap(~ Cell_Line, scales = "free_x", ncol = 1) +
  scale_fill_manual(values = cellline_colors) +
  labs(title = "RET driver fusion — evolution Untreated -> Treated -> Resistant (barplot)") +
  driver_labs + driver_theme

ggsave(file.path(OUT_DIR, "figures", "08_fusions_driver_barplot.png"),
       p_driver_bar, width = 11, height = 8, dpi = 300)
cat("  → figures/08_fusions_driver_barplot.png\n")

# Dotplot: combines two support metrics in one compact panel
# (size = JunctionReadCount, color = FFPM); more informative than line/bar alone
p_driver_dot <- ggplot(driver_traj, aes(x = Sample, y = 1,
                                          size = JunctionReadCount, color = FFPM)) +
  geom_point() +
  facet_wrap(~ Cell_Line, scales = "free_x", ncol = 1) +
  scale_color_gradientn(colors = c("#3498db", "#ecf0f1", "#922b21")) +
  scale_size_continuous(range = c(1, 12)) +
  labs(
    title = "RET driver fusion — evolution Untreated -> Treated -> Resistant (dotplot)",
    subtitle = paste0("RET1: ", driver_fusion["RET1"], "  |  RET2: ", driver_fusion["RET2"]),
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
       p_driver_dot, width = 11, height = 5.5, dpi = 300)
cat("  → figures/08_fusions_driver_dotplot.png\n\n")

# Remove old figure (previous name, replaced by lineplot+barplot)
old_fig <- file.path(OUT_DIR, "figures", "08_fusions_driver_trajectory.png")
if (file.exists(old_fig)) invisible(file.remove(old_fig))

# ------------------------------------------------------------------
# 5. FUSION HEATMAP PER MODEL (top N by FFPM, stage order)
#    + combined heatmap with the union of the top N of each model
# ------------------------------------------------------------------
cat("3. Fusion heatmaps per model (top", TOP_N_HEATMAP, "FFPM)...\n")

# The same gene pair can be reported with several breakpoints in the same
# sample -> sum FFPM by Sample+FusionName before building any matrix
fusions_agg <- fusions_all %>%
  group_by(Sample, `#FusionName`) %>%
  summarise(FFPM = sum(FFPM, na.rm = TRUE), .groups = "drop")

get_top_fusions <- function(cl, n = TOP_N_HEATMAP) {
  fusions_agg %>%
    filter(Sample %in% stage_orders[[cl]]) %>%
    group_by(`#FusionName`) %>%
    summarise(max_FFPM = max(FFPM, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(max_FFPM)) %>%
    slice_head(n = n) %>%
    pull(`#FusionName`)
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

build_ann_colors <- function(ann_col) {
  ann_colors <- filter_annotation_colors(ann_col)
  ann_colors <- setNames(lapply(names(ann_colors), function(nm) {
    vals <- ann_colors[[nm]]
    vals[names(vals) %in% unique(ann_col[[nm]])]
  }), names(ann_colors))

  # RET-specific colors (same as the RET WES heatmap and oncoprint scripts,
  # COLOR_DRUG/COLOR_MODEL, for visual consistency).
  if (!is.null(ann_colors$Drug)) {
    ann_colors$Drug <- c(
      "BLU667"  = "#6C5CE7",
      "LOXO292" = "#00B894",
      "No"      = "#95A5A6"
    )[names(ann_colors$Drug)]
  }
  if (!is.null(ann_colors$Cell_Line)) {
    ann_colors$Cell_Line <- c(
      "RET1" = "#7570B3",
      "RET2"  = "#E7298A"
    )[names(ann_colors$Cell_Line)]
  }
  ann_colors
}

# Legend labels (same criterion as the other RNA heatmaps): only the legend
# TITLE is renamed (Condition->Condition, Drug->Drug, Cell_Line->Model,
# same criterion as RET_08/RET_09) and the Condition values
# (Resistance->Resistant). Drug (BLU667/LOXO292/No) and Cell_Line
# (RET1/RET2) are NOT translated: they are technical codes/model names, not
# descriptive text. Applied AFTER build_ann_colors() (which matches by English
# column name against config/colors.R), not before: if ann_col were renamed
# first, filter_annotation_colors() would no longer find a match.
TRADUCCION_TITULO_ES <- c(Condition = "Condition", Drug = "Drug", Cell_Line = "Model")

traducir_leyenda_es <- function(ann_col, ann_colors) {
  if ("Condition" %in% colnames(ann_col)) {
    ann_col$Condition <- ifelse(ann_col$Condition == "Resistance", "Resistant", ann_col$Condition)
    if (!is.null(ann_colors$Condition)) {
      names(ann_colors$Condition) <- ifelse(names(ann_colors$Condition) == "Resistance", "Resistant", names(ann_colors$Condition))
    }
  }
  renombrar <- function(nombres) {
    nuevos <- TRADUCCION_TITULO_ES[nombres]
    ifelse(is.na(nuevos), nombres, nuevos)
  }
  colnames(ann_col) <- renombrar(colnames(ann_col))
  names(ann_colors) <- renombrar(names(ann_colors))
  list(ann_col = ann_col, ann_colors = ann_colors)
}

plot_fusion_heatmap <- function(cl, n = TOP_N_HEATMAP) {
  samples     <- stage_orders[[cl]]
  top_fusions <- get_top_fusions(cl, n)

  if (length(top_fusions) < 2) {
    cat("  ⚠", cl, "— fewer than 2 fusions, heatmap skipped\n")
    return(invisible(NULL))
  }

  mat_log <- build_fusion_matrix(samples, top_fusions)
  ann_col <- meta_all[samples, c("Condition", "Drug"), drop = FALSE]
  trad <- traducir_leyenda_es(ann_col, build_ann_colors(ann_col))

  pheatmap(mat_log,
           annotation_col = trad$ann_col,
           annotation_colors = trad$ann_colors,
           cluster_rows = TRUE,
           cluster_cols = FALSE,
           show_rownames = TRUE,
           show_colnames = TRUE,
           fontsize_row = 9,
           main = paste0("RNA fusions (top ", n, " FFPM) — ", cl),
           filename = file.path(OUT_DIR, "figures", paste0("08_fusions_heatmap_", cl, ".png")),
           width = 10, height = max(4, length(top_fusions) * 0.5))
  cat("  →", paste0("figures/08_fusions_heatmap_", cl, ".png"), "\n")
}
invisible(lapply(CELL_LINES, plot_fusion_heatmap))

# Combined heatmap: union of the top N of RET1 and the top N of RET2,
# all samples of both models, stage order per model
top_fusions_combined <- Reduce(union, lapply(CELL_LINES, get_top_fusions))
combined_samples <- unlist(stage_orders[CELL_LINES], use.names = FALSE)

if (length(top_fusions_combined) >= 2) {
  mat_log_combined <- build_fusion_matrix(combined_samples, top_fusions_combined)
  ann_col_combined <- meta_all[combined_samples, c("Cell_Line", "Condition", "Drug"), drop = FALSE]
  trad_combined <- traducir_leyenda_es(ann_col_combined, build_ann_colors(ann_col_combined))

  combined_fig <- paste0("08_fusions_heatmap_combined_top", TOP_N_HEATMAP, ".png")
  pheatmap(mat_log_combined,
           annotation_col = trad_combined$ann_col,
           annotation_colors = trad_combined$ann_colors,
           cluster_rows = TRUE,
           cluster_cols = FALSE,
           gaps_col = length(stage_orders[["RET1"]]),
           show_rownames = TRUE,
           show_colnames = TRUE,
           fontsize_row = 9,
           fontsize_col = 7,
           main = paste0("RNA fusions — union top ", TOP_N_HEATMAP, " per model (RET1 + RET2)"),
           filename = file.path(OUT_DIR, "figures", combined_fig),
           width = 16, height = max(4, length(top_fusions_combined) * 0.5))
  cat("  →", file.path("figures", combined_fig), "(", length(top_fusions_combined), "fusions )\n")
}

# Remove the combined figure from the previous top5 version if it exists
old_combined <- file.path(OUT_DIR, "figures", "08_fusions_heatmap_combined_top5.png")
if (TOP_N_HEATMAP != 5 && file.exists(old_combined)) invisible(file.remove(old_combined))
cat("\n")

# ------------------------------------------------------------------
# 6. MOST STRIKING FUSIONS PER MODEL
#    Criterion: robust support (LargeAnchorSupport=YES), high FFPM,
#    evidence in known databases (ChimerKB/COSMIC/Mitelman in 'annots')
# ------------------------------------------------------------------
cat("4. Most striking fusions per model...\n")

top_striking <- fusions_all %>%
  mutate(known_db = grepl("ChimerKB|Cosmic|COSMIC|Mitelman|ChimerSeq|FusionGDB", annots, ignore.case = TRUE)) %>%
  group_by(Cell_Line, `#FusionName`) %>%
  summarise(
    n_samples          = n_distinct(Sample),
    samples            = paste(sort(unique(Sample)), collapse = ", "),
    mean_FFPM          = round(mean(FFPM, na.rm = TRUE), 2),
    max_FFPM           = round(max(FFPM, na.rm = TRUE), 2),
    max_JunctionReads  = max(JunctionReadCount, na.rm = TRUE),
    any_LargeAnchor    = any(LargeAnchorSupport == "YES", na.rm = TRUE),
    prot_effect        = paste(unique(na.omit(PROT_FUSION_TYPE)), collapse = "/"),
    known_db           = any(known_db),
    .groups = "drop"
  ) %>%
  mutate(is_driver = `#FusionName` == driver_fusion[Cell_Line]) %>%
  filter(any_LargeAnchor | is_driver) %>%
  arrange(Cell_Line, desc(is_driver), desc(mean_FFPM)) %>%
  group_by(Cell_Line) %>%
  slice_head(n = TOP_N) %>%
  ungroup()

write.csv(top_striking, file.path(OUT_DIR, "tables", "08_fusions_top_striking.csv"), row.names = FALSE)
cat("  → tables/08_fusions_top_striking.csv (top", TOP_N, "per model)\n\n")

# ------------------------------------------------------------------
# 7. OVERLAPS BETWEEN MODELS
# ------------------------------------------------------------------
cat("5. Fusions shared between RET1 and RET2...\n")

by_model <- fusions_all %>%
  group_by(`#FusionName`, Cell_Line) %>%
  summarise(n_samples = n_distinct(Sample), mean_FFPM = round(mean(FFPM, na.rm = TRUE), 2), .groups = "drop")

shared <- by_model %>%
  group_by(`#FusionName`) %>%
  filter(n_distinct(Cell_Line) == length(CELL_LINES)) %>%
  ungroup() %>%
  pivot_wider(names_from = Cell_Line, values_from = c(n_samples, mean_FFPM)) %>%
  arrange(desc(pmax(!!!rlang::syms(paste0("mean_FFPM_", CELL_LINES)))))

write.csv(shared, file.path(OUT_DIR, "tables", "08_fusions_shared_models.csv"), row.names = FALSE)
cat("  → tables/08_fusions_shared_models.csv (", nrow(shared), "shared fusions)\n\n")

# ------------------------------------------------------------------
# 8. REPORT
# ------------------------------------------------------------------
driver_check <- driver_traj %>%
  group_by(Cell_Line) %>%
  summarise(min_FFPM = round(min(FFPM), 2), max_FFPM = round(max(FFPM), 2), .groups = "drop")

report_text <- paste0(
  "# RNA fusions (STAR-Fusion) — RET: RET1 + RET2\n\n",
  "## Source\n",
  "- ", length(fusion_files), " xlsx files in `data/fusions/` (STAR-Fusion + FusionInspector)\n",
  "- ", nrow(fusions_all), " total fusion events, ", length(unique(fusions_all$`#FusionName`)), " unique fusions\n\n",
  "## Driver fusion — confirmation at RNA level\n",
  paste(capture.output(print(as.data.frame(driver_check), row.names = FALSE)), collapse = "\n"), "\n",
  "- RET1: ", driver_fusion["RET1"], "\n",
  "- RET2: ", driver_fusion["RET2"], "\n\n",
  "## Most striking fusions per model (top ", TOP_N, ")\n",
  "See `tables/08_fusions_top_striking.csv`. Criterion: LargeAnchorSupport=YES and/or driver fusion, ",
  "ordered by mean FFPM; `known_db` is flagged if present in ChimerKB/COSMIC/Mitelman.\n\n",
  "## Overlaps between models\n",
  "See `tables/08_fusions_shared_models.csv` — ", nrow(shared), " fusions detected in both models (RET1 and RET2).\n\n",
  "## Interpretation\n",
  "[INTERPRET]\n",
  "- Is the driver fusion maintained in all resistant lines or lost in any (driver switch)?\n",
  "- Are the fusions shared between models technical noise (read-through) or biological candidates?\n",
  "- Does any new fusion appear only in resistant lines (acquired) vs. already present in the parental?\n\n",
  "## Figures\n",
  "- 08_fusions_driver_lineplot.png\n",
  "- 08_fusions_driver_barplot.png\n",
  "- 08_fusions_driver_dotplot.png (size=JunctionReadCount, color=FFPM)\n",
  "- 08_fusions_heatmap_RET1.png (top ", TOP_N_HEATMAP, ")\n",
  "- 08_fusions_heatmap_RET2.png (top ", TOP_N_HEATMAP, ")\n",
  "- 08_fusions_heatmap_combined_top", TOP_N_HEATMAP, ".png (union RET1 + RET2)\n"
)
writeLines(report_text, file.path(OUT_DIR, "report", "08_fusions_report.md"))
cat("Report: RET/combined/report/08_fusions_report.md\n\n")

cat("================================================================\n")
cat("  RET FUSIONS finished —", OUT_PREFIX, "\n")
cat("================================================================\n\n")
