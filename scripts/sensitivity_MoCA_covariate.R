# =============================================================================
# sensitivity_MoCA_covariate.R
# Sensitivity Check: Older CN Sex Differences with MoCA as Covariate
#
# Author:       Molly Simmonite
# Date:         April 2026
# Project:      MiND Study - MRS Analysis
#
# Description:  Runs older CN W1 sex difference analyses (script 05) with MoCA
#               added as an additional covariate. For comparison against primary
#               findings — not intended for reporting.
#
# Input:        data/MRS_data_analysis.csv
# Output:       outputs/sensitivity_MoCA_covariate_summary.txt
#
# Dependencies: helpers.R, tidyverse
# =============================================================================

source("scripts/helpers.R")

dir.create("outputs", recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# 1. LOAD AND PREPARE DATA
# =============================================================================

dat <- read_csv("data/MRS_data_analysis.csv")

older_CN <- dat %>%
  filter(AgeCategory == "Older",
         Subgroup_f  == "CN",
         Wave_Num    == 1)

n_with_moca    <- sum(!is.na(older_CN$MoCA))
n_without_moca <- sum(is.na(older_CN$MoCA))
cat("Older CN sample (wave 1): n =", nrow(older_CN), "\n")
cat("MoCA available:", n_with_moca, "| Missing:", n_without_moca, "\n\n")


# =============================================================================
# 2. RUN PRIMARY ANALYSES WITH MoCA AS COVARIATE
# =============================================================================

primary_voxels <- c("AUD", "SM", "VV")
metabolites    <- c("G", "GLX")
corrections    <- c("UNC", "ATC")

results_moca <- list()

for (met in metabolites) {
  for (corr in corrections) {
    label <- paste0(met, " - ", corr, " [+MoCA]")
    results_moca[[label]] <- run_manova_univariate(older_CN, met, corr, primary_voxels, label,
                                                   extra_covariates = "MoCA")
  }
}


# =============================================================================
# 3. ALSO RUN WITHOUT MoCA (PRIMARY) FOR DIRECT COMPARISON
# =============================================================================

results_primary <- list()

for (met in metabolites) {
  for (corr in corrections) {
    label <- paste0(met, " - ", corr, " [primary]")
    results_primary[[label]] <- run_manova_univariate(older_CN, met, corr, primary_voxels, label)
  }
}


# =============================================================================
# 4. PRINT SIDE-BY-SIDE COMPARISON
# =============================================================================

comparison_lines <- c(
  "================================================================",
  "SENSITIVITY CHECK: Older CN Sex Differences — MoCA as Covariate",
  "Comparing: primary model (Age only) vs +MoCA model",
  "* = p < .05",
  "================================================================",
  ""
)

met_corr_combos <- expand.grid(met = metabolites, corr = corrections, stringsAsFactors = FALSE)

for (i in seq_len(nrow(met_corr_combos))) {
  met  <- met_corr_combos$met[i]
  corr <- met_corr_combos$corr[i]
  
  key_primary <- paste0(met, " - ", corr, " [primary]")
  key_moca    <- paste0(met, " - ", corr, " [+MoCA]")
  
  r_p <- results_primary[[key_primary]]
  r_m <- results_moca[[key_moca]]
  
  comparison_lines <- c(comparison_lines, paste0("--- ", met, " - ", corr, " ---"))
  
  if (!is.null(r_p$error) || !is.null(r_m$error)) {
    comparison_lines <- c(comparison_lines, "  ERROR in one or both models", "")
    next
  }

  if (!is.null(r_p$manova$error_msg) || !is.null(r_m$manova$error_msg)) {
    comparison_lines <- c(comparison_lines, "  MANOVA extraction failed in one or both models", "")
    next
  }

  # MANOVA comparison
  p_manova_p <- r_p$manova$p
  p_manova_m <- r_m$manova$p
  
  comparison_lines <- c(comparison_lines,
                        paste0("  MANOVA (Pillai):"),
                        paste0("    Primary:  F(", r_p$manova$df1, ",", r_p$manova$df2, ") = ", r_p$manova$F,
                               ", p ", format_p(p_manova_p), ifelse(p_manova_p < .05, " *", "")),
                        paste0("    +MoCA:    F(", r_m$manova$df1, ",", r_m$manova$df2, ") = ", r_m$manova$F,
                               ", p ", format_p(p_manova_m), ifelse(p_manova_m < .05, " *", ""),
                               "  [n = ", r_m$n, "]")
  )
  
  # Univariate comparison
  comparison_lines <- c(comparison_lines, "  Univariates (Sex beta, p):")
  for (vox in primary_voxels) {
    u_p <- r_p$univariate[[vox]]
    u_m <- r_m$univariate[[vox]]
    if (is.null(u_p) || is.null(u_m)) next
    
    comparison_lines <- c(comparison_lines,
                          paste0("    ", vox, ":"),
                          paste0("      Primary:  beta = ", u_p$beta, ", p ", format_p(u_p$p), ifelse(u_p$p < .05, " *", "")),
                          paste0("      +MoCA:    beta = ", u_m$beta, ", p ", format_p(u_m$p), ifelse(u_m$p < .05, " *", ""))
    )
  }
  comparison_lines <- c(comparison_lines, "")
}

cat(paste(comparison_lines, collapse = "\n"), "\n")

writeLines(comparison_lines, "outputs/sensitivity_MoCA_covariate_summary.txt")
cat("\nSummary saved to outputs/sensitivity_MoCA_covariate_summary.txt\n")