# =============================================================================
# 08b_CN_vs_MCI_sex_bilateral.R
# Sex x Diagnosis Interaction - Bilateral Hemisphere Analysis
#
# Author:       Molly Simmonite
# Date:         April 2026
# Project:      MiND Study - MRS Analysis
#
# Description:  Tests for sex x diagnosis interactions in bilateral voxels.
#               Model: value ~ Sex * Diagnosis * Hemisphere + Age + (1|SubNum)
#               Key tests:
#                 Sex x Hemisphere            — lateralisation of sex difference
#                 Sex x Diagnosis             — diagnosis moderates sex difference
#                 Sex x Hemisphere x Diagnosis — lateralisation differs by diagnosis
#
# Voxel pairs:  AUD (LAUD/RAUD), SM (LSM/RSM), VV (LVV/RVV)
# Input:        data/MRS_data_analysis.csv
# Output:       outputs/08b_CN_vs_MCI_sex_bilateral_summary.txt
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

older_all <- dat %>%
  filter(AgeCategory == "Older",
         Subgroup_f  %in% c("CN", "MCI"),
         Wave_Num    == 1) %>%
  mutate(Subgroup_f = factor(Subgroup_f, levels = c("CN", "MCI")))

cat("Older adults sample (wave 1): n =", nrow(older_all), "\n")
cat("Sex x Diagnosis breakdown:\n")
print(table(older_all$Subgroup_f, older_all$Sex))
cat("\n")

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

run_hemisphere_diagnosis <- function(data, metabolite, correction, voxel_name, voxel_pair, label) {
  
  col_L <- paste0(metabolite, "_", correction, "_", voxel_pair[1])
  col_R <- paste0(metabolite, "_", correction, "_", voxel_pair[2])
  
  missing <- c(col_L, col_R)[!c(col_L, col_R) %in% names(data)]
  if (length(missing) > 0) {
    return(list(label = label, error = paste("Missing columns:", paste(missing, collapse = ", "))))
  }
  
  long_dat <- data %>%
    select(SubNum, Sex, Subgroup_f, Age_DuringParticipation,
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
    lmer(value ~ Sex * Hemisphere * Subgroup_f + Age_DuringParticipation + (1 | SubNum),
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
    label         = label,
    n_subs        = n_subs,
    n_obs         = n_obs,
    n_male        = n_male,
    n_female      = n_female,
    sex           = extract_term("SexMale"),
    hemisphere    = extract_term("HemisphereR"),
    diagnosis     = extract_term("Subgroup_fMCI"),
    sex_hemi      = extract_term("SexMale:HemisphereR"),
    sex_dx        = extract_term("SexMale:Subgroup_fMCI"),
    hemi_dx       = extract_term("HemisphereR:Subgroup_fMCI"),
    three_way     = extract_term("SexMale:HemisphereR:Subgroup_fMCI"),
    error         = NULL
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
      results[[label]] <- run_hemisphere_diagnosis(
        older_all, met, corr, vox_name, voxel_pairs[[vox_name]], label
      )
    }
  }
}


# =============================================================================
# 4. FORMAT AND SAVE SUMMARY
# =============================================================================

format_bilateral_dx_summary <- function(results, title) {
  lines <- c(
    "================================================================",
    paste0("SUMMARY: ", title),
    "Model: value ~ Sex * Hemisphere * Diagnosis + Age + (1|SubNum)",
    "Reference: Female, CN, Left hemisphere",
    "Key tests: Sex x Hemisphere, Sex x Diagnosis, Sex x Hemisphere x Diagnosis",
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
               fmt(res$sex,       "Sex                                    "),
               fmt(res$hemisphere,"Hemisphere (R vs L)                    "),
               fmt(res$diagnosis, "Diagnosis (MCI vs CN)                  "),
               fmt(res$sex_hemi,  "Sex x Hemisphere [KEY]                 "),
               fmt(res$sex_dx,    "Sex x Diagnosis [KEY]                  "),
               fmt(res$hemi_dx,   "Hemisphere x Diagnosis                 "),
               fmt(res$three_way, "Sex x Hemisphere x Diagnosis [KEY]     ")
    )
    lines <- c(lines, "")
  }
  
  lines
}

summary_lines <- format_bilateral_dx_summary(results, "Bilateral Sex x Diagnosis Analysis in Older Adults (Wave 1)")
cat(paste(summary_lines, collapse = "\n"), "\n")

writeLines(summary_lines, "outputs/08b_CN_vs_MCI_sex_bilateral_summary.txt")
cat("\nSummary saved to outputs/08b_CN_vs_MCI_sex_bilateral_summary.txt\n")