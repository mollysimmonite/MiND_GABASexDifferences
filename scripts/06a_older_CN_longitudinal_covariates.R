# =============================================================================
# 06a_older_CN_longitudinal_covariates.R
# Longitudinal Sex Differences - Older CN Adults (with covariates)
#
# Author:       Molly Simmonite
# Date:         April 2026
# Project:      MiND Study - MRS Analysis
#
# Description:  Follow-up to 06_older_CN_longitudinal.R. Repeats significant
#               findings (GLX UNC VV) adding MRSTime and FMRI_ScannerUpgrade
#               as covariates. All voxels repeated for completeness.
#
# Input:        data/MRS_data_analysis.csv
#               data/Participation_Dates.csv
# Output:       outputs/06a_older_CN_longitudinal_covariates_summary.txt
#
# Dependencies: helpers.R, tidyverse, lme4, lmerTest
# =============================================================================

source("scripts/helpers.R")
library(lme4)
library(lmerTest)

dir.create("outputs", recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# 1. LOAD AND PREPARE DATA
# =============================================================================

dat   <- read_csv("data/MRS_data_analysis.csv")
dates <- read_csv("data/Participation_Dates.csv")

# Compute time between waves (years from wave 1)
dates_clean <- dates %>%
  mutate(SubNum  = as.character(SubNum),
         MRSDate = as.Date(MRSDate, format = "%m/%d/%y")) %>%
  filter(!is.na(MRSDate), !is.na(SubNum)) %>%
  select(SubNum, Wave_Num, MRSDate)

wave1_dates <- dates_clean %>%
  filter(Wave_Num == 1) %>%
  select(SubNum, Wave1_Date = MRSDate)

dates_with_time <- dates_clean %>%
  left_join(wave1_dates, by = "SubNum") %>%
  mutate(Time_Years = as.numeric(difftime(MRSDate, Wave1_Date, units = "days")) / 365.25) %>%
  select(SubNum, Wave_Num, Time_Years)

# Older CN adults, waves 1 and 2 only
older_CN_long <- dat %>%
  filter(AgeCategory == "Older",
         Subgroup_f  == "CN",
         Wave_Num    %in% c(1, 2)) %>%
  mutate(SubNum = as.character(SubNum)) %>%
  left_join(dates_with_time, by = c("SubNum", "Wave_Num")) %>%
  mutate(Time_Years = ifelse(Wave_Num == 1, 0, Time_Years))

cat("Older CN longitudinal sample:\n")
cat("  Total rows:", nrow(older_CN_long), "\n")
cat("  Participants:", n_distinct(older_CN_long$SubNum), "\n\n")


# =============================================================================
# 2. LME HELPER WITH EXTRA COVARIATES
# =============================================================================

run_lme_voxel_covs <- function(data, col, label, extra_covariates = NULL) {
  
  cov_cols <- c("SubNum", "Sex", "Time_Years", "Age_DuringParticipation", extra_covariates)
  
  analysis_dat <- data %>%
    select(all_of(cov_cols), all_of(col)) %>%
    drop_na()
  
  n_total  <- nrow(analysis_dat)
  n_subs   <- n_distinct(analysis_dat$SubNum)
  n_male   <- n_distinct(analysis_dat$SubNum[analysis_dat$Sex == "Male"])
  n_female <- n_distinct(analysis_dat$SubNum[analysis_dat$Sex == "Female"])
  
  if (n_total < 10 || n_male < 3 || n_female < 3) {
    return(list(label = label, col = col, error = paste("Insufficient sample (n =", n_total, ")")))
  }
  
  all_terms <- c("Sex * Time_Years", "Age_DuringParticipation", extra_covariates)
  formula <- as.formula(paste0(
    "`", col, "` ~ ", paste(all_terms, collapse = " + "), " + (1 | SubNum)"
  ))
  
  fit <- tryCatch(
    lmer(formula, data = analysis_dat, REML = FALSE),
    error = function(e) NULL
  )
  
  if (is.null(fit)) {
    return(list(label = label, col = col, error = "Model failed to converge"))
  }
  
  coefs <- summary(fit)$coefficients
  
  extract_row <- function(coefs, term) {
    if (!term %in% rownames(coefs)) return(NULL)
    row <- coefs[term, ]
    list(
      beta = round(row["Estimate"], 4),
      t    = round(row["t value"], 3),
      p    = round(row["Pr(>|t|)"], 4),
      sig  = ifelse(row["Pr(>|t|)"] < .05, "*", "")
    )
  }
  
  list(
    label       = label,
    col         = col,
    n           = n_total,
    n_subs      = n_subs,
    n_male      = n_male,
    n_female    = n_female,
    sex         = extract_row(coefs, "SexMale"),
    time        = extract_row(coefs, "Time_Years"),
    interaction = extract_row(coefs, "SexMale:Time_Years"),
    error       = NULL
  )
}


# =============================================================================
# 3. RUN ANALYSES
# =============================================================================

primary_voxels <- c("AUD", "SM", "VV")
metabolites    <- c("G", "GLX")
corrections    <- c("UNC", "ATC")
extra_covs     <- c("MRSTime", "FMRI_ScannerUpgrade")

results <- list()

for (met in metabolites) {
  for (corr in corrections) {
    for (vox in primary_voxels) {
      col   <- paste0(met, "_", corr, "_", vox)
      label <- paste0(met, " - ", corr, " - ", vox)
      results[[label]] <- run_lme_voxel_covs(older_CN_long, col, label,
                                             extra_covariates = extra_covs)
    }
  }
}


# =============================================================================
# 4. FORMAT AND SAVE SUMMARY
# =============================================================================

format_lme_summary <- function(results, title) {
  lines <- c(
    "================================================================",
    paste0("SUMMARY: ", title),
    "Model: outcome ~ Sex * Time_Years + Age + MRSTime + ScannerUpgrade + (1 | SubNum)",
    "Primary interest: Sex x Time interaction",
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
    
    lines <- c(lines, paste0(
      "  n = ", res$n, " observations, ",
      res$n_subs, " participants (",
      res$n_male, "M / ", res$n_female, "F)"
    ))
    
    format_term <- function(term_result, name) {
      if (is.null(term_result)) return(paste0("  ", name, ": not estimated"))
      paste0("  ", name, ": beta = ", term_result$beta,
             ", t = ", term_result$t,
             ", p ", format_p(term_result$p), " ", term_result$sig)
    }
    
    lines <- c(lines,
               format_term(res$sex,         "Sex (male vs female)  "),
               format_term(res$time,        "Time                  "),
               format_term(res$interaction, "Sex x Time (key)      ")
    )
    lines <- c(lines, "")
  }
  
  lines
}

summary_lines <- format_lme_summary(results, "Longitudinal Change in Older CN Adults - with covariates")
cat(paste(summary_lines, collapse = "\n"), "\n")

writeLines(summary_lines, "outputs/06a_older_CN_longitudinal_covariates_summary.txt")
cat("\nSummary saved to outputs/06a_older_CN_longitudinal_covariates_summary.txt\n")