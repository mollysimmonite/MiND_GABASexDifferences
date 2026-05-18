# =============================================================================
# 0_descriptives.R
# Sample Descriptive Statistics
#
# Author:       Molly Simmonite
# Date:         April 2026
# Project:      MiND Study — MRS Analysis
#
# Description:  Generates descriptive statistics tables for three groups:
#               Young CN, Older CN, and Older MCI (wave 1 only).
#               Each table summarises n, age, MoCA, education, and MothersEd
#               by sex, with t-tests comparing males vs females.
#
# Input:        data/MRS_data_analysis.csv
# Output:       outputs/0_descriptives_summary.txt
#               outputs/0_descriptives_young_CN.csv
#               outputs/0_descriptives_older_CN.csv
#               outputs/0_descriptives_older_MCI.csv
#
# Dependencies: tidyverse
# =============================================================================

library(tidyverse)

dir.create("outputs", recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# 1. LOAD DATA
# =============================================================================

dat <- read_csv("data/MRS_data_analysis.csv")


# =============================================================================
# 2. HELPER FUNCTIONS
# =============================================================================

# Summarise a continuous variable by sex and run t-test
summarise_by_sex <- function(data, var) {
  var_sym <- sym(var)
  
  # Summary stats by sex
  stats <- data %>%
    group_by(Sex) %>%
    summarise(
      n    = sum(!is.na(!!var_sym)),
      mean = round(mean(!!var_sym, na.rm = TRUE), 2),
      sd   = round(sd(!!var_sym,   na.rm = TRUE), 2),
      .groups = "drop"
    ) %>%
    mutate(mean_sd = paste0(mean, " (", sd, ")")) %>%
    select(Sex, n, mean_sd)
  
  # T-test
  male_vals   <- data[[var]][data$Sex == "Male"]
  female_vals <- data[[var]][data$Sex == "Female"]
  
  male_vals   <- male_vals[!is.na(male_vals)]
  female_vals <- female_vals[!is.na(female_vals)]
  
  if (length(male_vals) < 2 || length(female_vals) < 2) {
    t_result <- list(statistic = NA, parameter = NA, p.value = NA)
  } else {
    t_result <- t.test(male_vals, female_vals)
  }
  
  list(
    variable = var,
    stats    = stats,
    t        = round(t_result$statistic, 3),
    df       = round(t_result$parameter, 1),
    p        = round(t_result$p.value, 4),
    sig      = ifelse(!is.na(t_result$p.value) & t_result$p.value < .05, "*", "")
  )
}

# Build a formatted table for a group
build_table <- function(data, group_label, variables) {
  
  lines <- c(
    paste0("=== ", group_label, " ==="),
    paste0("n = ", nrow(data),
           " (", sum(data$Sex == "Male"), "M / ",
           sum(data$Sex == "Female"), "F)"),
    ""
  )
  
  rows <- list()
  
  for (var in variables) {
    result <- summarise_by_sex(data, var)
    stats  <- result$stats
    
    male_val   <- stats$mean_sd[stats$Sex == "Male"]
    female_val <- stats$mean_sd[stats$Sex == "Female"]
    
    if (length(male_val)   == 0) male_val   <- "NA"
    if (length(female_val) == 0) female_val <- "NA"
    
    lines <- c(lines, paste0(
      sprintf("%-20s", var),
      "  Male: ", sprintf("%-15s", male_val),
      "  Female: ", sprintf("%-15s", female_val),
      "  t(", result$df, ") = ", result$t,
      ", p = ", result$p, " ", result$sig
    ))
    
    rows[[var]] <- data.frame(
      Variable      = var,
      Male_mean_sd  = male_val,
      Female_mean_sd = female_val,
      t             = result$t,
      df            = result$df,
      p             = result$p,
      sig           = result$sig,
      stringsAsFactors = FALSE
    )
  }
  
  list(lines = lines, rows = bind_rows(rows))
}


# =============================================================================
# 3. DEFINE VARIABLES AND GROUPS
# =============================================================================

variables <- c("Age_DuringParticipation", "MoCA", "Education", "MothersEd")

young_CN   <- dat %>% filter(AgeCategory == "Young", Subgroup_f == "CN",  Wave_Num == 1)
older_CN   <- dat %>% filter(AgeCategory == "Older", Subgroup_f == "CN",  Wave_Num == 1)
older_MCI  <- dat %>% filter(AgeCategory == "Older", Subgroup_f == "MCI", Wave_Num == 1)


# =============================================================================
# 4. BUILD TABLES
# =============================================================================

young_CN_table  <- build_table(young_CN,  "Young CN (Wave 1)",  variables)
older_CN_table  <- build_table(older_CN,  "Older CN (Wave 1)",  variables)
older_MCI_table <- build_table(older_MCI, "Older MCI (Wave 1)", variables)


# =============================================================================
# 5. PRINT AND SAVE
# =============================================================================

all_lines <- c(
  "================================================================",
  "DESCRIPTIVE STATISTICS — MiND Study MRS Analysis",
  "Values: Mean (SD). T-tests comparing males vs females.",
  "* = p < .05",
  "================================================================",
  "",
  young_CN_table$lines,
  "",
  older_CN_table$lines,
  "",
  older_MCI_table$lines
)

cat(paste(all_lines, collapse = "\n"), "\n")

writeLines(all_lines, "outputs/0_descriptives_summary.txt")
cat("\nSummary saved to outputs/0_descriptives_summary.txt\n")

# Save CSVs
write_csv(young_CN_table$rows,  "outputs/0_descriptives_young_CN.csv")
write_csv(older_CN_table$rows,  "outputs/0_descriptives_older_CN.csv")
write_csv(older_MCI_table$rows, "outputs/0_descriptives_older_MCI.csv")
cat("Tables saved to outputs/\n")