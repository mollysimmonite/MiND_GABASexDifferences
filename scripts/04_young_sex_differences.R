# =============================================================================
# 4_young_sex_differences.R
# Sex Differences in GABA and Glx — Young Adults
#
# Author:       Molly Simmonite
# Date:         April 2026
# Project:      MiND Study — MRS Analysis
#
# Description:  Examines sex differences in GABA and Glx in young CN adults
#               (cross-sectional). Primary analyses use bilateral average voxels
#               (AUD, SM, VV). Exploratory analyses repeat for LPV.
#               MANOVA used as omnibus test across voxels; univariate linear
#               models run as follow-up if MANOVA is significant.
#               Age included as covariate in all models.
#
# Input:        data/MRS_data_analysis.csv
# Output:       outputs/4_young_sex_differences_summary.txt
#
# Dependencies: helpers.R, tidyverse
# =============================================================================

source("scripts/helpers.R")

dir.create("outputs", recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# 1. LOAD AND PREPARE DATA
# =============================================================================

dat <- read_csv("data/MRS_data_analysis.csv")

# Young CN adults only, wave 1 only
young <- dat %>%
  filter(AgeCategory == "Young",
         Subgroup_f  == "CN",
         Wave_Num    == 1)

cat("Young CN sample: n =", nrow(young), "\n")
cat("Sex breakdown:\n")
print(table(young$Sex))
cat("\n")


# =============================================================================
# 2. DEFINE ANALYSIS PARAMETERS
# =============================================================================

primary_voxels     <- c("AUD", "SM", "VV")
exploratory_voxels <- c("LPV")
metabolites        <- c("G", "GLX")
corrections        <- c("UNC", "ATC")


# =============================================================================
# 3. RUN PRIMARY ANALYSES
# =============================================================================

results <- list()

for (met in metabolites) {
  for (corr in corrections) {
    label <- paste0(met, " — ", corr)
    results[[label]] <- run_manova_univariate(young, met, corr, primary_voxels, label)
  }
}


# =============================================================================
# 4. RUN EXPLORATORY ANALYSES (LPV)
# =============================================================================

for (met in metabolites) {
  for (corr in corrections) {
    label <- paste0(met, " — ", corr, " [EXPLORATORY: LPV]")
    results[[label]] <- run_univariate_single(young, met, corr, "LPV", label)
  }
}


# =============================================================================
# 5. PRINT AND SAVE SUMMARY
# =============================================================================

summary_lines <- format_summary(results, title = "Sex Differences in Young CN Adults")
cat(paste(summary_lines, collapse = "\n"), "\n")

writeLines(summary_lines, "outputs/4_young_sex_differences_summary.txt")
cat("\nSummary saved to outputs/4_young_sex_differences_summary.txt\n")