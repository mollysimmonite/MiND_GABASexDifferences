# =============================================================================
# 6_older_CN_longitudinal.R
# Longitudinal Sex Differences in GABA and Glx — Older CN Adults
#
# Author:       Molly Simmonite
# Date:         April 2026
# Project:      MiND Study — MRS Analysis
#
# Description:  Examines whether GABA and Glx change differently by sex over
#               time (wave 1 to wave 2) in older CN adults. Uses linear mixed
#               models with Sex, Time (years since wave 1), and Sex x Time
#               interaction as fixed effects, Age and SubNum as random
#               intercept. Primary interest is the Sex x Time interaction.
#               Wave 3 excluded due to small n.
#               Primary analyses use bilateral average voxels (AUD, SM, VV).
#               MEM excluded from longitudinal analysis as only 2 participants
#               have repeated MEM data.
#
# Input:        data/MRS_data_analysis.csv
#               data/Participation_Dates.csv
# Output:       outputs/6_older_CN_longitudinal_summary.txt
#
# Dependencies: helpers.R, tidyverse, lme4, lmerTest
# =============================================================================

source("scripts/helpers.R")
library(lme4)
library(lmerTest)  # gives p-values for lmer

dir.create("outputs", recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# 1. LOAD DATA
# =============================================================================

dat   <- read_csv("data/MRS_data_analysis.csv")
dates <- read_csv("data/Participation_Dates.csv")

# =============================================================================
# 2. COMPUTE TIME BETWEEN WAVES
# =============================================================================
# For each participant, calculate time in years from their wave 1 MRS date.
# This gives a continuous time variable rather than just wave number.

dates_clean <- dates %>%
  mutate(
    SubNum  = as.character(SubNum),
    MRSDate = as.Date(MRSDate, format = "%m/%d/%y")  # adjust format if needed
  ) %>%
  filter(!is.na(MRSDate), !is.na(SubNum)) %>%
  select(SubNum, Wave_Num, MRSDate)

# Get wave 1 date per participant as anchor
wave1_dates <- dates_clean %>%
  filter(Wave_Num == 1) %>%
  select(SubNum, Wave1_Date = MRSDate)

# Calculate time in years from wave 1 for each wave
dates_with_time <- dates_clean %>%
  left_join(wave1_dates, by = "SubNum") %>%
  mutate(
    Time_Years = as.numeric(difftime(MRSDate, Wave1_Date, units = "days")) / 365.25
  ) %>%
  select(SubNum, Wave_Num, Time_Years)

cat("=== Time Between Waves (years) ===\n")
dates_with_time %>%
  filter(Wave_Num == 2) %>%
  summarise(
    n    = n(),
    mean = round(mean(Time_Years, na.rm = TRUE), 2),
    sd   = round(sd(Time_Years,   na.rm = TRUE), 2),
    min  = round(min(Time_Years,  na.rm = TRUE), 2),
    max  = round(max(Time_Years,  na.rm = TRUE), 2)
  ) %>%
  print()
cat("\n")


# =============================================================================
# 3. PREPARE ANALYSIS DATASET
# =============================================================================

# Older CN adults, waves 1 and 2 only (exclude wave 3)
older_CN_long <- dat %>%
  filter(AgeCategory == "Older",
         Subgroup_f  == "CN",
         Wave_Num    %in% c(1, 2)) %>%
  mutate(SubNum = as.character(SubNum)) %>%
  left_join(dates_with_time, by = c("SubNum", "Wave_Num"))

# Wave 1 time should be 0
older_CN_long <- older_CN_long %>%
  mutate(Time_Years = ifelse(Wave_Num == 1, 0, Time_Years))

cat("Older CN longitudinal sample:\n")
cat("  Total rows:", nrow(older_CN_long), "\n")
cat("  Participants:", n_distinct(older_CN_long$SubNum), "\n")
cat("  Wave breakdown:\n")
print(table(older_CN_long$Wave_Num, older_CN_long$Sex))
cat("\n")

# Flag participants with missing Time_Years at wave 2
missing_time <- older_CN_long %>%
  filter(Wave_Num == 2, is.na(Time_Years))
if (nrow(missing_time) > 0) {
  cat("Warning:", nrow(missing_time), "wave 2 rows missing Time_Years — will be excluded from models\n\n")
}


# =============================================================================
# 4. HELPER FUNCTION FOR LME
# =============================================================================
# Separate from the MANOVA helpers — runs one LME per voxel and extracts
# the Sex x Time interaction term.

run_lme_voxel <- function(data, col, label) {
  
  analysis_dat <- data %>%
    select(SubNum, Sex, Time_Years, Age_DuringParticipation, all_of(col)) %>%
    drop_na()
  
  n_total  <- nrow(analysis_dat)
  n_subs   <- n_distinct(analysis_dat$SubNum)
  n_male   <- n_distinct(analysis_dat$SubNum[analysis_dat$Sex == "Male"])
  n_female <- n_distinct(analysis_dat$SubNum[analysis_dat$Sex == "Female"])
  
  if (n_total < 10 || n_male < 3 || n_female < 3) {
    return(list(label = label, col = col, error = paste("Insufficient sample (n =", n_total, ")")))
  }
  
  formula <- as.formula(paste0(
    "`", col, "` ~ Sex * Time_Years + Age_DuringParticipation + (1 | SubNum)"
  ))
  
  fit <- tryCatch(
    lmer(formula, data = analysis_dat, REML = FALSE),
    error = function(e) NULL
  )
  
  if (is.null(fit)) {
    return(list(label = label, col = col, error = "Model failed to converge"))
  }
  
  coefs <- summary(fit)$coefficients
  
  # Extract Sex main effect, Time main effect, and Sex x Time interaction
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
    label      = label,
    col        = col,
    n          = n_total,
    n_subs     = n_subs,
    n_male     = n_male,
    n_female   = n_female,
    sex        = extract_row(coefs, "SexMale"),
    time       = extract_row(coefs, "Time_Years"),
    interaction = extract_row(coefs, "SexMale:Time_Years"),
    error      = NULL
  )
}


# =============================================================================
# 5. RUN ANALYSES
# =============================================================================

primary_voxels     <- c("AUD", "SM", "VV")
metabolites        <- c("G", "GLX")
corrections        <- c("UNC", "ATC")

results <- list()

for (met in metabolites) {
  for (corr in corrections) {
    for (vox in primary_voxels) {
      col   <- paste0(met, "_", corr, "_", vox)
      label <- paste0(met, " — ", corr, " — ", vox)
      results[[label]] <- run_lme_voxel(older_CN_long, col, label)
    }
  }
}


# =============================================================================
# 6. FORMAT AND PRINT SUMMARY
# =============================================================================

format_lme_summary <- function(results, title) {
  lines <- c(
    "================================================================",
    paste0("SUMMARY: ", title),
    "Model: outcome ~ Sex * Time_Years + Age + (1 | SubNum)",
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

summary_lines <- format_lme_summary(results, "Longitudinal Change in Older CN Adults (W1-W2)")
cat(paste(summary_lines, collapse = "\n"), "\n")

writeLines(summary_lines, "outputs/6_older_CN_longitudinal_summary.txt")
cat("\nSummary saved to outputs/6_older_CN_longitudinal_summary.txt\n")