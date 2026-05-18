# =============================================================================
# 06b_older_CN_longitudinal_bilateral.R
# Longitudinal Bilateral Hemisphere Analysis - Older CN Adults
#
# Author:       Molly Simmonite
# Date:         April 2026
# Project:      MiND Study - MRS Analysis
#
# Description:  Tests for sex differences, lateralisation, and their change
#               over time in older CN adults using bilateral voxel pairs.
#               Model: value ~ Sex * Hemisphere * Time_Years + Age + (1|SubNum)
#               Key tests:
#                 Sex x Hemisphere            — lateralisation of sex difference
#                 Sex x Time                  — change over time by sex
#                 Sex x Hemisphere x Time     — does lateralisation change over time?
#
# Voxel pairs:  AUD (LAUD/RAUD), SM (LSM/RSM), VV (LVV/RVV)
# Input:        data/MRS_data_analysis.csv, data/Participation_Dates.csv
# Output:       outputs/06b_older_CN_longitudinal_bilateral_summary.txt
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

# Compute Time_Years from wave 1
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

# Older CN, waves 1 and 2
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

voxel_pairs <- list(
  AUD = c("LAUD", "RAUD"),
  SM  = c("LSM",  "RSM"),
  VV  = c("LVV",  "RVV")
)

metabolites <- c("G", "GLX")
corrections <- c("UNC", "ATC")


# =============================================================================
# 2. HELPER: RESHAPE AND RUN THREE-WAY MODEL
# =============================================================================

run_hemisphere_longitudinal <- function(data, metabolite, correction, voxel_name, voxel_pair, label) {
  
  col_L <- paste0(metabolite, "_", correction, "_", voxel_pair[1])
  col_R <- paste0(metabolite, "_", correction, "_", voxel_pair[2])
  
  missing <- c(col_L, col_R)[!c(col_L, col_R) %in% names(data)]
  if (length(missing) > 0) {
    return(list(label = label, error = paste("Missing columns:", paste(missing, collapse = ", "))))
  }
  
  long_dat <- data %>%
    select(SubNum, Sex, Age_DuringParticipation, Time_Years,
           L = all_of(col_L), R = all_of(col_R)) %>%
    pivot_longer(cols = c(L, R), names_to = "Hemisphere", values_to = "value") %>%
    mutate(Hemisphere = factor(Hemisphere, levels = c("L", "R"))) %>%
    drop_na()
  
  n_subs   <- n_distinct(long_dat$SubNum)
  n_male   <- n_distinct(long_dat$SubNum[long_dat$Sex == "Male"])
  n_female <- n_distinct(long_dat$SubNum[long_dat$Sex == "Female"])
  n_obs    <- nrow(long_dat)
  
  if (n_subs < 10 || n_male < 3 || n_female < 3) {
    return(list(label = label, error = paste("Insufficient sample (n =", n_subs, ")")))
  }
  
  fit <- tryCatch(
    lmer(value ~ Sex * Hemisphere * Time_Years + Age_DuringParticipation + (1 | SubNum),
         data = long_dat, REML = FALSE),
    error = function(e) NULL
  )
  
  if (is.null(fit)) {
    return(list(label = label, error = "Model failed"))
  }
  
  coefs <- summary(fit)$coefficients
  
  extract_term <- function(term) {
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
    label            = label,
    n_subs           = n_subs,
    n_obs            = n_obs,
    n_male           = n_male,
    n_female         = n_female,
    sex              = extract_term("SexMale"),
    hemisphere       = extract_term("HemisphereR"),
    time             = extract_term("Time_Years"),
    sex_hemi         = extract_term("SexMale:HemisphereR"),
    sex_time         = extract_term("SexMale:Time_Years"),
    hemi_time        = extract_term("HemisphereR:Time_Years"),
    three_way        = extract_term("SexMale:HemisphereR:Time_Years"),
    error            = NULL
  )
}


# =============================================================================
# 3. RUN ANALYSES
# =============================================================================

results <- list()

for (met in metabolites) {
  for (corr in corrections) {
    for (vox_name in names(voxel_pairs)) {
      label <- paste0(met, " - ", corr, " - ", vox_name)
      results[[label]] <- run_hemisphere_longitudinal(
        older_CN_long, met, corr, vox_name, voxel_pairs[[vox_name]], label
      )
    }
  }
}


# =============================================================================
# 4. FORMAT AND SAVE SUMMARY
# =============================================================================

format_bilateral_long_summary <- function(results, title) {
  lines <- c(
    "================================================================",
    paste0("SUMMARY: ", title),
    "Model: value ~ Sex * Hemisphere * Time_Years + Age + (1|SubNum)",
    "Hemisphere: L = left (ref), R = right",
    "Key tests: Sex x Hemisphere, Sex x Time, Sex x Hemisphere x Time",
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
      "  n = ", res$n_subs, " participants (", res$n_male, "M / ", res$n_female, "F), ",
      res$n_obs, " observations"
    ))
    
    fmt <- function(term_result, name) {
      if (is.null(term_result)) return(paste0("  ", name, ": not estimated"))
      paste0("  ", name, ": beta = ", term_result$beta,
             ", t = ", term_result$t,
             ", p ", format_p(term_result$p), " ", term_result$sig)
    }
    
    lines <- c(lines,
               fmt(res$sex,       "Sex                           "),
               fmt(res$hemisphere,"Hemisphere (R vs L)           "),
               fmt(res$time,      "Time                          "),
               fmt(res$sex_hemi,  "Sex x Hemisphere [KEY]        "),
               fmt(res$sex_time,  "Sex x Time [KEY]              "),
               fmt(res$hemi_time, "Hemisphere x Time             "),
               fmt(res$three_way, "Sex x Hemisphere x Time [KEY] ")
    )
    lines <- c(lines, "")
  }
  
  lines
}

summary_lines <- format_bilateral_long_summary(results, "Longitudinal Bilateral Analysis in Older CN Adults")
cat(paste(summary_lines, collapse = "\n"), "\n")

writeLines(summary_lines, "outputs/06b_older_CN_longitudinal_bilateral_summary.txt")
cat("\nSummary saved to outputs/06b_older_CN_longitudinal_bilateral_summary.txt\n")