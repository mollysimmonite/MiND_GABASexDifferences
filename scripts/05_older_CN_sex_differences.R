# =============================================================================
# 5_older_CN_sex_differences.R
# Sex Differences in GABA and Glx — Older CN Adults
#
# Author:       Molly Simmonite
# Date:         April 2026
# Project:      MiND Study — MRS Analysis
#
# Description:  Examines sex differences in GABA and Glx in older CN adults
#               (cross-sectional, wave 1 only). Primary analyses use bilateral
#               average voxels (AUD, SM, VV). Exploratory analyses repeat for
#               MEM and LPV.
#               MANOVA used as omnibus test across voxels; univariate linear
#               models run as follow-up if MANOVA is significant.
#               Age included as covariate in all models.
#
# Input:        data/MRS_data_analysis.csv
# Output:       outputs/5_older_CN_sex_differences_summary.txt
#
# Dependencies: helpers.R, tidyverse
# =============================================================================

source("scripts/helpers.R")

dir.create("outputs", recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# 1. LOAD AND PREPARE DATA
# =============================================================================

dat <- read_csv("data/MRS_data_analysis.csv")

# Older CN adults, wave 1 only
older_CN <- dat %>%
  filter(AgeCategory == "Older",
         Subgroup_f  == "CN",
         Wave_Num    == 1)

cat("Older CN sample (wave 1): n =", nrow(older_CN), "\n")
cat("Sex breakdown:\n")
print(table(older_CN$Sex))
cat("\n")


# =============================================================================
# 2. DEFINE ANALYSIS PARAMETERS
# =============================================================================

primary_voxels     <- c("AUD", "SM", "VV")
exploratory_voxels <- c("MEM", "LPV")
metabolites        <- c("G", "GLX")
corrections        <- c("UNC", "ATC")


# =============================================================================
# 3. RUN PRIMARY ANALYSES
# =============================================================================

results <- list()

for (met in metabolites) {
  for (corr in corrections) {
    label <- paste0(met, " — ", corr)
    results[[label]] <- run_manova_univariate(older_CN, met, corr, primary_voxels, label)
  }
}


# =============================================================================
# 4. RUN EXPLORATORY ANALYSES (MEM, LPV)
# =============================================================================

for (met in metabolites) {
  for (corr in corrections) {
    for (vox in exploratory_voxels) {
      label <- paste0(met, " — ", corr, " [EXPLORATORY: ", vox, "]")
      results[[label]] <- run_univariate_single(older_CN, met, corr, vox, label)
    }
  }
}


# =============================================================================
# 5. PRINT AND SAVE SUMMARY
# =============================================================================

summary_lines <- format_summary(results, title = "Sex Differences in Older CN Adults (Wave 1)")
cat(paste(summary_lines, collapse = "\n"), "\n")

writeLines(summary_lines, "outputs/5_older_CN_sex_differences_summary.txt")
cat("\nSummary saved to outputs/5_older_CN_sex_differences_summary.txt\n")