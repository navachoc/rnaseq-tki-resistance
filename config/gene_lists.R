# gene_lists.R
# Reference gene lists for TKI resistance in NSCLC (ALK, MET, EGFR, RET),
# shared by the RNAseq_bulk scripts.
# Usage: source(here::here("config", "gene_lists.R")). Exports all_gene_lists.

# ==============================================================================
# ONCOGENES AND DRIVERS
# ==============================================================================

oncogenes <- c(
  "ALK", "EML4", "RET", "ROS1", "MET", "AXL",
  "EGFR", "KRAS", "NRAS", "HRAS", "BRAF", "TP53", "MYC"
)

# ALK bypass (alternative receptors activated in resistance)
alk_bypass <- c(
  "EGFR", "MET", "ERBB2", "ERBB3", "IGF1R",
  "KIT", "SRC", "YES1", "AXL"
)

# ==============================================================================
# SIGNALING PATHWAYS
# ==============================================================================

pi3k_akt_pathway <- c(
  "PIK3CA", "PIK3CB", "AKT1", "AKT2", "AKT3",
  "MTOR", "PTEN", "PDPK1", "IRS1", "IGF1R",
  "TSC1", "TSC2", "FOXO1", "FOXO3"
)

mapk_erk_pathway <- c(
  "KRAS", "HRAS", "NRAS", "BRAF", "RAF1",
  "MAP2K1", "MAP2K2", "MAPK1", "MAPK3",
  "DUSP6", "SPRY2", "FOS", "JUN"
)

jak_stat_pathway <- c(
  "JAK1", "JAK2", "JAK3",
  "STAT1", "STAT3", "STAT5A", "STAT5B",
  "IL6", "IL6R"
)

# ==============================================================================
# EMT (Epithelial-mesenchymal transition)
# ==============================================================================

emt_tfs <- c("SNAI1", "SNAI2", "TWIST1", "TWIST2", "ZEB1", "ZEB2")

mesenchymal_markers <- c(
  "VIM", "FN1", "CDH2", "CDH11",
  "MMP2", "MMP9", "MMP14", "SPARC",
  "ITGB1", "ITGA5", "TNC"
)

epithelial_markers <- c(
  "CDH1", "OCLN", "CLDN1", "CLDN4", "EPCAM", "KRT8", "KRT18"
)

emt_signature <- unique(c(emt_tfs, mesenchymal_markers, epithelial_markers))

# ==============================================================================
# APOPTOSIS AND CELL CYCLE
# ==============================================================================

apoptosis_genes <- c(
  "BCL2", "BCL2L1", "MCL1", "BCL2L11", "BAX", "BAK1",
  "CASP3", "CASP8", "CASP9", "TP53", "MDM2", "CDKN1A", "CDKN2A"
)

cell_cycle <- c(
  "CCND1", "CCNE1", "CCNA2", "CCNB1",
  "CDK2", "CDK4", "CDK6",
  "CDKN1A", "CDKN1B", "CDKN2A",
  "RB1", "E2F1", "MKI67"
)

# ==============================================================================
# STEMNESS
# ==============================================================================

stemness_genes <- c(
  "NANOG", "SOX2", "KLF4", "MYC",
  "ALDH1A1", "CD44", "PROM1", "BMI1", "LIN28A"
)

# ==============================================================================
# ALK RESISTANCE CORE MASTER LIST
# Priority genes for visualization and WB validation
# ==============================================================================

alk_resistance_core <- unique(c(
  alk_bypass,
  "ALK", "EML4", "AXL", "MET", "EGFR", "ERBB2", "ERBB3", "IGF1R",
  "YAP1", "TEAD1", "KRAS", "BRAF", "RAF1", "MAP2K1", "MAPK1",
  "PIK3CA", "AKT1", "AKT2", "AKT3", "MTOR", "PTEN",
  "STAT3", "JAK1", "JAK2", "TP53", "RB1", "CDKN2A",
  "BCL2", "MCL1", "BCL2L1", "VIM", "CDH1", "CDH2", "ZEB1", "ZEB2",
  "SNAI1", "SNAI2", "FN1", "CCND1", "CDK4", "CDK6"
))

# ==============================================================================
# MET RESISTANCE — Specific signatures (Giménez-Xavier et al. 2017, Mol Cancer Ther)
# Original paper of the MET1 models (PMID 28396363)
# ==============================================================================

# MET genes and their signaling context
met_core <- c(
  "MET", "EGFR", "KRAS", "BRAF", "AXL", "ERBB2", "ERBB3", "IGF1R",
  "MAPK1", "AKT1", "MTOR", "PIK3CA", "TP53", "CDK4", "CDK6", "CCND1", "MYC"
)

# AKT/ERK/mTOR reactivation — common to ALL resistant MET1 models (Giménez-Xavier 2017)
gimenez_common <- c("AKT1", "AKT2", "AKT3", "MAPK1", "MAPK3", "MTOR",
                    "RPS6KB1", "EIF4EBP1", "MAP2K1", "MAP2K2")

# NF2 inactivation + YAP overexpression — group 2 of resistant MET1 models (Giménez-Xavier 2017)
gimenez_group2 <- c("NF2", "WWTR1", "YAP1", "CTGF", "CYR61", "ANKRD1",
                    "TEAD1", "TEAD2", "TEAD3", "TEAD4")

# MET-HER crosstalk → erlotinib sensitivity (Giménez-Xavier 2017)
gimenez_egfr_bypass <- c("EGFR", "ERBB2", "ERBB3", "ERBB4",
                          "HBEGF", "AREG", "EREG", "NRG1")

# MET bypass (alternative receptors activated in resistance)
met_bypass <- c("EGFR", "KRAS", "BRAF", "AXL", "RET", "FGFR1", "FGFR2",
                "ERBB2", "IGF1R")

# MET resistance core master list (for visualization and WB validation)
met_resistance_core <- unique(c(
  met_core, met_bypass,
  "NF2", "YAP1", "WWTR1", "CTGF", "CYR61",
  "STAT3", "JAK1", "JAK2",
  "BCL2", "MCL1", "BCL2L1",
  "VIM", "CDH1", "CDH2", "ZEB1", "ZEB2", "SNAI1", "SNAI2",
  "CCND1", "CDK4", "CDK6", "RB1", "CDKN2A"
))

# ==============================================================================
# EGFR resistance — EGFR1 / EGFR1_R (Quintanal-Villalonga et al. 2020)
# Histological transformation EGFR-TKI → SCLC + LuSCC (in vivo resistance)
# Mechanisms: RB1 inactivation, FGFR3-TACC3 fusion, MYCL reamplification
egfr_resistance_core <- unique(c(
  # Inactivated driver
  "RB1",
  # EGFR1_R fusion
  "FGFR3", "TACC3",
  # EGFR1_R amplification
  "MYCL", "MYC",
  # SCLC markers (Quintanal-Villalonga 2020, Fig 4E)
  "ASCL2", "SYP", "NCAM1", "NTRK2", "CXCR4", "LCP1", "ADH7", "KRT75",
  # LuSCC markers
  "TP63", "SOX2", "NKX2-1", "ALDH3A1", "KRT4", "KRT40",
  # EGFR pathway
  "EGFR", "ERBB2", "ERBB3", "KRAS", "BRAF", "MEK1",
  # Plasticity / dedifferentiation
  "VIM", "CDH1", "CDH2", "ZEB1", "SNAI1",
  # Potential bypass
  "AXL", "MET", "IGF1R", "FGFR1"
))

# ==============================================================================
# RET resistance — genes with an established role in TKI resistance
# (bypass RTKs, RAS/MAPK, PI3K/AKT, EMT/ECM, innate immunity)
# Used to screen the full signature (not just the heatmap top 15) and
# prioritize genes with real literature support over uncharacterized genes.
# ==============================================================================

ret_bypass <- c(
  "EGFR", "MET", "ALK", "LTK", "ROS1", "AXL",
  "ERBB2", "ERBB3", "IGF1R", "KIT", "FGFR1", "FGFR2"
)

ecm_adhesion <- c(
  "FN1", "MMP2", "MMP9", "MMP14", "MMP19",
  "COL1A1", "COL1A2", "COL3A1", "COL5A1",
  "ITGA5", "ITGB1", "ITGA9", "SPARC", "TNC",
  "FAK", "PTK2", "YAP1", "WWTR1"
)

innate_immune_chemotaxis <- c(
  "CXCL8", "CXCL1", "CXCL2", "CXCL3", "CCL2", "CCL5",
  "ITGAM", "ITGB2", "SELL", "CD44", "IL6", "IL1B"
)

xenobiotic_metabolism <- c(
  "CYP1A1", "CYP1A2", "CYP1B1",
  "UGT1A1", "UGT1A6", "UGT1A9"
)

ret_resistance_core <- unique(c(
  oncogenes, ret_bypass, mapk_erk_pathway, pi3k_akt_pathway,
  jak_stat_pathway, emt_signature, ecm_adhesion,
  innate_immune_chemotaxis, xenobiotic_metabolism
))

# ==============================================================================
# Curated WES panel — ALK
# List derived from the WES analysis of ALK1 (same panel used in the three
# ALK WES models: ALK1/ALK2/ALK3). Used to cross-reference
# WES findings (SNV/CNV/LOH) with the RNA-seq expression of those same genes.
# ==============================================================================

wes_panel_alk <- c(
  # RTK/RAS/RAF
  "EGFR", "ALK", "ROS1", "RET", "MET", "ERBB2", "KRAS", "NRAS", "HRAS", "BRAF",
  "MAP2K1", "FGFR1", "FGFR2", "FGFR3",
  # PI3K/AKT/mTOR
  "PIK3CA", "PTEN", "AKT1", "AKT2", "AKT3", "MTOR", "TSC1", "TSC2",
  # Cell cycle
  "CDKN2A", "CDK4", "CDK6", "CCND1", "RB1", "TP53",
  # DNA repair
  "BRCA1", "BRCA2", "ATM", "CHEK2", "PALB2", "MLH1", "MSH2", "MSH6", "PMS2",
  # Epigenetic regulators
  "SETD2", "KMT2A", "KMT2C", "KMT2D", "ARID1A", "SMARCA4", "SMARCB1",
  # Suppressors / frequent co-mutations in NSCLC
  "STK11", "KEAP1", "NFE2L2", "NF1", "NOTCH1", "FBXW7", "CTNNB1", "APC"
)

# ==============================================================================
# EXPORTED LIST
# ==============================================================================

all_gene_lists <- list(
  oncogenes           = oncogenes,
  alk_bypass          = alk_bypass,
  pi3k_akt_pathway    = pi3k_akt_pathway,
  mapk_erk_pathway    = mapk_erk_pathway,
  jak_stat_pathway    = jak_stat_pathway,
  emt_tfs             = emt_tfs,
  emt_signature       = emt_signature,
  mesenchymal_markers = mesenchymal_markers,
  epithelial_markers  = epithelial_markers,
  apoptosis_genes     = apoptosis_genes,
  cell_cycle          = cell_cycle,
  stemness_genes      = stemness_genes,
  alk_resistance_core = alk_resistance_core,
  wes_panel_alk        = wes_panel_alk,
  # MET-specific
  met_core            = met_core,
  met_bypass          = met_bypass,
  met_resistance_core = met_resistance_core,
  gimenez_common      = gimenez_common,
  gimenez_group2      = gimenez_group2,
  gimenez_egfr_bypass = gimenez_egfr_bypass,
  # EGFR-specific
  egfr_resistance_core = egfr_resistance_core,
  # RET-specific
  ret_bypass               = ret_bypass,
  ecm_adhesion             = ecm_adhesion,
  innate_immune_chemotaxis = innate_immune_chemotaxis,
  xenobiotic_metabolism    = xenobiotic_metabolism,
  ret_resistance_core      = ret_resistance_core
)

cat("gene_lists.R loaded\n")
cat("   alk_resistance_core: ", length(alk_resistance_core),  "genes\n")
cat("   wes_panel_alk:       ", length(wes_panel_alk),        "genes\n")
cat("   met_resistance_core: ", length(met_resistance_core),  "genes\n")
cat("   egfr_resistance_core:", length(egfr_resistance_core), "genes\n")
cat("   ret_resistance_core: ", length(ret_resistance_core),  "genes\n\n")
