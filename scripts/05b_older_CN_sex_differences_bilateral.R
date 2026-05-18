# =============================================================================
# 05b_older_CN_sex_differences_bilateral.R
# Sex Differences - Older CN, Bilateral Hemisphere Analysis
#
# Author:       Molly Simmonite
# Date:         April 2026
# Project:      MiND Study - MRS Analysis
#
# Description:  Tests for sex differences and sex x hemisphere interactions
#               in older CN adults (wave 1) using bilateral voxel pairs (L/R).
#               Mixed model: value ~ Sex * Hemisphere + Age + (1|SubNum)
#               Key test: Sex x Hemisphere interaction (lateralisation).
#
# Voxel pairs:  AUD (LAUD/RAUD), SM (LSM/RSM), VV (LVV/RVV)
# Input:        data/MRS_data_analysis.csv
# Output:       outputs/05b_older_CN_sex_differences_bilateral_summary.txt
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

dat <- read_csv("data/MRS_data_analysis.csv")

older_CN <- dat %>%
  filter(AgeCategory == "Older",
         Subgroup_f  == "CN",
         Wave_Num    == 1)

cat("Older CN sample (wave 1): n =", nrow(older_CN), "\n")
cat("Sex breakdown:\n")
print(table(older_CN$Sex))
cat("\n")

voxel_pairs <- list(
  AUD = c("LAUD", "RAUD"),
  SM  = c("LSM",  "RSM"),
  VV  = c("LVV",  "RVV")
)

metabolites <- c("G", "GLX")
corrections <- c("UNC", "ATC")


# =============================================================================
# 2. HELPER: RESHAPE TO LONG AND RUN HEMISPHERE MODEL
# =============================================================================

run_hemisphere_model <- function(data, metabolite, correction, voxel_name, voxel_pair, label) {
  
  col_L <- paste0(metabolite, "_", correction, "_", voxel_pair[1])
  col_R <- paste0(metabolite, "_", correction, "_", voxel_pair[2])
  
  missing <- c(col_L, col_R)[!c(col_L, col_R) %in% names(data)]
  if (length(missing) > 0) {
    return(list(label = label, error = paste("Missing columns:", paste(missing, collapse = ", "))))
  }
  
  long_dat <- data %>%
    select(SubNum, Sex, Age_DuringParticipation, L = all_of(col_L), R = all_of(col_R)) %>%
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
    lmer(value ~ Sex * Hemisphere + Age_DuringParticipation + (1 | SubNum),
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
    label       = label,
    n_subs      = n_subs,
    n_obs       = n_obs,
    n_male      = n_male,
    n_female    = n_female,
    sex         = extract_term("SexMale"),
    hemisphere  = extract_term("HemisphereR"),
    interaction = extract_term("SexMale:HemisphereR"),
    error       = NULL
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
      results[[label]] <- run_hemisphere_model(
        older_CN, met, corr, vox_name, voxel_pairs[[vox_name]], label
      )
    }
  }
}


# =============================================================================
# 4. FORMAT AND SAVE SUMMARY
# =============================================================================

format_bilateral_summary <- function(results, title) {
  lines <- c(
    "================================================================",
    paste0("SUMMARY: ", title),
    "Model: value ~ Sex * Hemisphere + Age + (1|SubNum)",
    "Hemisphere: L = left (ref), R = right",
    "Key test: Sex x Hemisphere interaction (lateralisation)",
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
               fmt(res$sex,         "Sex (male vs female)          "),
               fmt(res$hemisphere,  "Hemisphere (R vs L)           "),
               fmt(res$interaction, "Sex x Hemisphere [KEY]        ")
    )
    lines <- c(lines, "")
  }
  
  lines
}

summary_lines <- format_bilateral_summary(results, "Bilateral Sex Differences in Older CN Adults (Wave 1)")
cat(paste(summary_lines, collapse = "\n"), "\n")

writeLines(summary_lines, "outputs/05b_older_CN_sex_differences_bilateral_summary.txt")
cat("\nSummary saved to outputs/05b_older_CN_sex_differences_bilateral_summary.txt\n")