# =============================================================================
# master_run.R
# Master Run Script
#
# Author:       Molly Simmonite
# Date:         April 2026
# Project:      MiND Study — MRS Analysis
#
# Description:  Runs all analysis scripts in order. Comment out any scripts
#               you don't want to run. Scripts must be run from the project
#               root directory (where data/ and scripts/ folders live).
#
# Order:
#   Data preparation (must run in order 1 -> 2 -> 3):
#     1_create_dataset.R        — load, clean, filter raw data
#     2_MRS_QC_check.R          — flag implausible MRS values
#     3_MRS_flags.R             — apply flags, GM exclusions, save analysis dataset
#
#   Descriptives (run before analyses):
#     0_descriptives.R          — sample characteristics tables
#
#   Primary analyses:
#     4_young_sex_differences.R          — sex differences, young CN
#     5_older_CN_sex_differences.R       — sex differences, older CN
#     6_older_CN_longitudinal.R          — longitudinal change, older CN
#     7_older_MCI_sex_differences.R      — sex differences, older MCI
#     8_CN_vs_MCI_sex.R                  — sex x diagnosis interaction
#
#   Follow-up analyses (run after reviewing primary results):
#     04a_young_sex_differences_covariates.R     — young CN with covariates
#     05a_older_CN_sex_differences_covariates.R  — older CN with covariates
#     06a_older_CN_longitudinal_covariates.R     — longitudinal with covariates
#     07a_older_MCI_sex_differences_covariates.R — older MCI with covariates
#     08a_CN_vs_MCI_sex_covariates.R             — sex x diagnosis with covariates
#
#   Bilateral hemisphere analyses:
#     04b_young_sex_differences_bilateral.R      — young CN, Sex x Hemisphere
#     05b_older_CN_sex_differences_bilateral.R   — older CN, Sex x Hemisphere
#     06b_older_CN_longitudinal_bilateral.R      — older CN longitudinal, Sex x Hemisphere x Time
#     07b_older_MCI_sex_differences_bilateral.R  — older MCI, Sex x Hemisphere
#     08b_CN_vs_MCI_sex_bilateral.R              — Sex x Hemisphere x Diagnosis
#
# =============================================================================

cat("================================================================\n")
cat("MiND Study MRS Analysis — Master Run Script\n")
cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("================================================================\n\n")

run_script <- function(path) {
  cat("----------------------------------------------------------------\n")
  cat("Running:", path, "\n")
  cat("----------------------------------------------------------------\n")
  tryCatch(
    source(path),
    error = function(e) cat("ERROR in", path, ":\n", conditionMessage(e), "\n")
  )
  cat("\n")
}

# =============================================================================
# DATA PREPARATION
# =============================================================================

run_script("scripts/01_create_dataset.R")
run_script("scripts/02_MRS_QC_check.R")
run_script("scripts/03_MRS_flags.R")

# =============================================================================
# DESCRIPTIVES
# =============================================================================

run_script("scripts/00_descriptives.R")

# =============================================================================
# PRIMARY ANALYSES
# =============================================================================

run_script("scripts/04_young_sex_differences.R")
run_script("scripts/05_older_CN_sex_differences.R")
run_script("scripts/06_older_CN_longitudinal.R")
run_script("scripts/07_older_MCI_sex_differences.R")
run_script("scripts/08_CN_vs_MCI_sex.R")

# =============================================================================
# FOLLOW-UP ANALYSES
# =============================================================================

run_script("scripts/04a_young_sex_differences_covariates.R")
run_script("scripts/05a_older_CN_sex_differences_covariates.R")
run_script("scripts/06a_older_CN_longitudinal_covariates.R")
run_script("scripts/07a_older_MCI_sex_differences_covariates.R")
run_script("scripts/08a_CN_vs_MCI_sex_covariates.R")

# =============================================================================
# BILATERAL HEMISPHERE ANALYSES
# =============================================================================

run_script("scripts/04b_young_sex_differences_bilateral.R")
run_script("scripts/05b_older_CN_sex_differences_bilateral.R")
run_script("scripts/06b_older_CN_longitudinal_bilateral.R")
run_script("scripts/07b_older_MCI_sex_differences_bilateral.R")
run_script("scripts/08b_CN_vs_MCI_sex_bilateral.R")

cat("================================================================\n")
cat("All scripts complete.\n")
cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("================================================================\n")