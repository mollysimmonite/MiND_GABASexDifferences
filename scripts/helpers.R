# =============================================================================
# helpers.R
# Shared Helper Functions for MRS Sex Differences Analyses
#
# Author:       Molly Simmonite
# Date:         April 2026
# Project:      MiND Study — MRS Analysis
#
# Description:  Defines shared functions used across analysis scripts 04–08.
#               Source this script at the top of each analysis script.
#               Functions are added here as needed across analyses.
#
# Functions:
#   run_manova_univariate()  — MANOVA + univariate follow-ups for a set of voxels
#   run_univariate_single()  — univariate linear model for a single voxel
#   format_summary()         — formats results list into printable summary lines
#
# Dependencies: tidyverse
# =============================================================================

library(tidyverse)


# =============================================================================
# run_manova_univariate()
# =============================================================================
# Runs MANOVA across a set of voxels with Sex as predictor and Age as covariate,
# followed by univariate linear models per voxel.
#
# Arguments:
#   data             — dataframe (already filtered to relevant sample)
#   metabolite       — "G" or "GLX"
#   correction       — "UNC" or "ATC"
#   voxels           — character vector of voxel names e.g. c("AUD", "SM", "VV")
#   label            — string label for output
#   extra_covariates — character vector of additional covariate column names beyond
#                      Sex and Age (default: NULL). Used in 'a' scripts.

run_manova_univariate <- function(data, metabolite, correction, voxels, label,
                                  extra_covariates = NULL) {
  
  cols <- paste0(metabolite, "_", correction, "_", voxels)
  
  # Check columns exist
  missing <- cols[!cols %in% names(data)]
  if (length(missing) > 0) {
    return(list(label = label, error = paste("Missing columns:", paste(missing, collapse = ", "))))
  }
  
  # Base covariates always included
  covariate_cols <- c("Sex", "Age_DuringParticipation")
  if (!is.null(extra_covariates)) {
    covariate_cols <- c(covariate_cols, extra_covariates)
  }
  
  # Drop rows with any NA in outcomes or covariates
  analysis_dat <- data %>%
    select(all_of(c(covariate_cols, cols))) %>%
    drop_na()
  
  n_total  <- nrow(analysis_dat)
  n_male   <- sum(analysis_dat$Sex == "Male")
  n_female <- sum(analysis_dat$Sex == "Female")
  
  if (n_total < 10 || n_male < 3 || n_female < 3) {
    return(list(label = label, error = paste("Insufficient sample size (n =", n_total, ")")))
  }
  
  # Build formula
  covariate_str  <- paste(c("Sex", "Age_DuringParticipation", extra_covariates), collapse = " + ")
  outcome_matrix <- as.matrix(analysis_dat[, cols])
  
  # MANOVA using base R manova() — cleaner extractable structure than car::Manova()
  manova_fit <- tryCatch(
    manova(as.formula(paste0("outcome_matrix ~ ", covariate_str)), data = analysis_dat),
    error = function(e) NULL
  )
  
  if (is.null(manova_fit)) {
    return(list(label = label, error = "MANOVA failed"))
  }
  
  manova_summary <- tryCatch({
    s         <- summary(manova_fit, test = "Pillai")
    sex_row   <- s$stats["Sex", ]
    list(
      pillai = round(sex_row["Pillai"], 3),
      F      = round(sex_row["approx F"], 3),
      df1    = round(sex_row["num Df"], 0),
      df2    = round(sex_row["den Df"], 0),
      p      = round(sex_row["Pr(>F)"], 4),
      sig    = ifelse(sex_row["Pr(>F)"] < .05, "*", "")
    )
  }, error = function(e) {
    list(error_msg = paste("MANOVA extraction failed:", e$message))
  })
  
  # Univariate follow-ups
  uni_results <- list()
  for (col in cols) {
    vox <- str_extract(col, "[^_]+$")
    fit <- lm(as.formula(paste0("`", col, "` ~ ", covariate_str)), data = analysis_dat)
    coef_sex <- summary(fit)$coefficients["SexMale", ]
    uni_results[[vox]] <- list(
      beta = round(coef_sex["Estimate"], 4),
      t    = round(coef_sex["t value"], 3),
      p    = round(coef_sex["Pr(>|t|)"], 4),
      sig  = ifelse(coef_sex["Pr(>|t|)"] < .05, "*", "")
    )
  }
  
  list(
    label    = label,
    n        = n_total,
    n_male   = n_male,
    n_female = n_female,
    manova   = manova_summary,
    univariate = uni_results,
    error    = NULL
  )
}


# =============================================================================
# run_univariate_single()
# =============================================================================
# Runs a univariate linear model for a single voxel.
# Used for exploratory single-voxel analyses (LPV, MEM).
#
# Arguments: same as run_manova_univariate() but voxel is a single string.

run_univariate_single <- function(data, metabolite, correction, voxel, label,
                                  extra_covariates = NULL) {
  
  col <- paste0(metabolite, "_", correction, "_", voxel)
  
  if (!col %in% names(data)) {
    return(list(label = label, error = paste("Column not found:", col)))
  }
  
  covariate_cols <- c("Sex", "Age_DuringParticipation")
  if (!is.null(extra_covariates)) {
    covariate_cols <- c(covariate_cols, extra_covariates)
  }
  
  analysis_dat <- data %>%
    select(all_of(c(covariate_cols, col))) %>%
    drop_na()
  
  n_total  <- nrow(analysis_dat)
  n_male   <- sum(analysis_dat$Sex == "Male")
  n_female <- sum(analysis_dat$Sex == "Female")
  
  if (n_total < 10 || n_male < 3 || n_female < 3) {
    return(list(label = label, error = paste("Insufficient sample size (n =", n_total, ")")))
  }
  
  covariate_str <- paste(c("Sex", "Age_DuringParticipation", extra_covariates), collapse = " + ")
  fit           <- lm(as.formula(paste0("`", col, "` ~ ", covariate_str)), data = analysis_dat)
  coef_sex      <- summary(fit)$coefficients["SexMale", ]
  
  uni_result <- list(
    beta = round(coef_sex["Estimate"], 4),
    t    = round(coef_sex["t value"], 3),
    p    = round(coef_sex["Pr(>|t|)"], 4),
    sig  = ifelse(coef_sex["Pr(>|t|)"] < .05, "*", "")
  )
  uni_list <- list()
  uni_list[[voxel]] <- uni_result
  
  list(
    label      = label,
    n          = n_total,
    n_male     = n_male,
    n_female   = n_female,
    manova     = NULL,
    univariate = uni_list,
    error      = NULL
  )
}


# =============================================================================
# format_summary()
# =============================================================================
# Formats a results list into a character vector of printable summary lines.
#
# Arguments:
#   results    — named list of result objects from run_manova_univariate()
#                or run_univariate_single()
#   title      — string title for the summary block
#   predictor  — string describing the predictor (default: "Sex (ref = Female)")

format_p <- function(p) {
  if (is.na(p)) return("NA")
  if (p < .0001) return("< .0001")
  paste0("= ", formatC(p, format = "f", digits = 4))
}

format_summary <- function(results, title, predictor = "Sex (ref = Female); Covariate: Age") {
  lines <- c(
    "================================================================",
    paste0("SUMMARY: ", title),
    paste0("Predictor: ", predictor),
    "* = p < .05",
    "================================================================",
    ""
  )
  
  for (res in results) {
    lines <- c(lines, paste0("--- ", res$label, " ---"))
    
    if (!is.null(res$error)) {
      lines <- c(lines, paste0("  ERROR: ", res$error), "")
      next
    }
    
    lines <- c(lines, paste0("  n = ", res$n, " (", res$n_male, "M / ", res$n_female, "F)"))
    
    if (!is.null(res$manova)) {
      m <- res$manova
      if (!is.null(m$error_msg)) {
        lines <- c(lines, paste0("  MANOVA: ", m$error_msg))
      } else {
        lines <- c(lines, paste0(
          "  MANOVA: Pillai = ", m$pillai,
          ", F(", m$df1, ",", m$df2, ") = ", m$F,
          ", p ", format_p(m$p), " ", m$sig
        ))
      }
    }
    
    lines <- c(lines, "  Univariate (+ = higher in males):")
    for (vox in names(res$univariate)) {
      u <- res$univariate[[vox]]
      lines <- c(lines, paste0(
        "    ", vox, ": beta = ", u$beta,
        ", t = ", u$t,
        ", p ", format_p(u$p), " ", u$sig
      ))
    }
    lines <- c(lines, "")
  }
  
  lines
}