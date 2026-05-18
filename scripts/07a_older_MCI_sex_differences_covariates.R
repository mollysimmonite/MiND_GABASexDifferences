# =============================================================================
# 07a_older_MCI_sex_differences_covariates.R
# Sex Differences in GABA and Glx - Older MCI Adults (with covariates)
#
# Author:       Molly Simmonite
# Date:         April 2026
# Project:      MiND Study - MRS Analysis
#
# Description:  Follow-up to 07_older_MCI_sex_differences.R. Repeats all
#               analyses adding MRSTime (time of scan) and
#               FMRI_ScannerUpgrade as covariates to check robustness of
#               primary findings.
#
# Input:        data/MRS_data_analysis.csv
# Output:       outputs/07a_older_MCI_sex_differences_covariates_summary.txt
#
# Dependencies: helpers.R, tidyverse
# =============================================================================

source("scripts/helpers.R")

dir.create("outputs", recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# 1. LOAD AND PREPARE DATA
# =============================================================================

dat <- read_csv("data/MRS_data_analysis.csv")

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
extra_covs         <- c("MRSTime", "FMRI_ScannerUpgrade")


# =============================================================================
# 3. RUN PRIMARY ANALYSES WITH COVARIATES
# =============================================================================

results <- list()

for (met in metabolites) {
  for (corr in corrections) {
    label <- paste0(met, " - ", corr)
    results[[label]] <- run_manova_univariate(older_MCI, met, corr, primary_voxels, label,
                                              extra_covariates = extra_covs)
  }
}


# =============================================================================
# 4. RUN EXPLORATORY ANALYSES WITH COVARIATES (MEM)
# =============================================================================

for (met in metabolites) {
  for (corr in corrections) {
    label <- paste0(met, " - ", corr, " [EXPLORATORY: MEM]")
    results[[label]] <- run_univariate_single(older_MCI, met, corr, "MEM", label,
                                              extra_covariates = extra_covs)
  }
}


# =============================================================================
# 5. PRINT AND SAVE SUMMARY
# =============================================================================

summary_lines <- format_summary(
  results,
  title     = "Sex Differences in Older MCI Adults (Wave 1) - with covariates",
  predictor = "Sex (ref = Female); Covariates: Age, MRSTime, FMRI_ScannerUpgrade"
)
cat(paste(summary_lines, collapse = "\n"), "\n")

writeLines(summary_lines, "outputs/07a_older_MCI_sex_differences_covariates_summary.txt")
cat("\nSummary saved to outputs/07a_older_MCI_sex_differences_covariates_summary.txt\n")