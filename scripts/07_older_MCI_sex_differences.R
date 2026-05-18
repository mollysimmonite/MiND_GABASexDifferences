# =============================================================================
# 7_older_MCI_sex_differences.R
# Sex Differences in GABA and Glx — Older MCI Adults
#
# Author:       Molly Simmonite
# Date:         April 2026
# Project:      MiND Study — MRS Analysis
#
# Description:  Examines sex differences in GABA and Glx in older MCI adults
#               (cross-sectional, wave 1 only). Primary analyses use bilateral
#               average voxels (AUD, SM, VV). Exploratory analyses repeat for
#               MEM. Note smaller sample than CN analyses — interpret with
#               caution. Age included as covariate in all models.
#
# Input:        data/MRS_data_analysis.csv
# Output:       outputs/7_older_MCI_sex_differences_summary.txt
#
# Dependencies: helpers.R, tidyverse
# =============================================================================

source("scripts/helpers.R")

dir.create("outputs", recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# 1. LOAD AND PREPARE DATA
# =============================================================================

dat <- read_csv("data/MRS_data_analysis.csv")

# Older MCI adults, wave 1 only
older_MCI <- dat %>%
  filter(AgeCategory == "Older",
         Subgroup_f  == "MCI",
         Wave_Num    == 1)

cat("Older MCI sample (wave 1): n =", nrow(older_MCI), "\n")
cat("Sex breakdown:\n")
print(table(older_MCI$Sex))
cat("\n")


# =============================================================================
# 2. DEFINE ANALYSIS PARAMETERS
# =============================================================================

primary_voxels     <- c("AUD", "SM", "VV")
exploratory_voxels <- c("MEM")
metabolites        <- c("G", "GLX")
corrections        <- c("UNC", "ATC")


# =============================================================================
# 3. RUN PRIMARY ANALYSES
# =============================================================================

results <- list()

for (met in metabolites) {
  for (corr in corrections) {
    label <- paste0(met, " — ", corr)
    results[[label]] <- run_manova_univariate(older_MCI, met, corr, primary_voxels, label)
  }
}


# =============================================================================
# 4. RUN EXPLORATORY ANALYSES (MEM)
# =============================================================================

for (met in metabolites) {
  for (corr in corrections) {
    label <- paste0(met, " — ", corr, " [EXPLORATORY: MEM]")
    results[[label]] <- run_univariate_single(older_MCI, met, corr, "MEM", label)
  }
}


# =============================================================================
# 5. PRINT AND SAVE SUMMARY
# =============================================================================

summary_lines <- format_summary(results, title = "Sex Differences in Older MCI Adults (Wave 1)")
cat(paste(summary_lines, collapse = "\n"), "\n")

writeLines(summary_lines, "outputs/7_older_MCI_sex_differences_summary.txt")
cat("\nSummary saved to outputs/7_older_MCI_sex_differences_summary.txt\n")