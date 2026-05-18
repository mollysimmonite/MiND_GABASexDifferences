# =============================================================================
# 8_CN_vs_MCI_sex.R
# Sex x Diagnosis Interaction in GABA and Glx — Older Adults
#
# Author:       Molly Simmonite
# Date:         April 2026
# Project:      MiND Study — MRS Analysis
#
# Description:  Examines whether sex differences in GABA and Glx differ
#               between CN and MCI older adults (cross-sectional, wave 1 only).
#               Primary interest is the Sex x Diagnosis interaction.
#               Primary analyses use bilateral average voxels (AUD, SM, VV).
#               Exploratory analyses repeat for MEM.
#               Age included as covariate in all models.
#
# Input:        data/MRS_data_analysis.csv
# Output:       outputs/8_CN_vs_MCI_sex_summary.txt
#
# Dependencies: helpers.R, tidyverse
# =============================================================================

source("scripts/helpers.R")

dir.create("outputs", recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# 1. LOAD AND PREPARE DATA
# =============================================================================

dat <- read_csv("data/MRS_data_analysis.csv")

# All older adults (CN and MCI), wave 1 only
older_all <- dat %>%
  filter(AgeCategory == "Older",
         Subgroup_f  %in% c("CN", "MCI"),
         Wave_Num    == 1) %>%
  mutate(Subgroup_f = factor(Subgroup_f, levels = c("CN", "MCI")))  # CN as reference

cat("Older adults sample (wave 1): n =", nrow(older_all), "\n")
cat("Sex x Diagnosis breakdown:\n")
print(table(older_all$Subgroup_f, older_all$Sex))
cat("\n")


# =============================================================================
# 2. DEFINE ANALYSIS PARAMETERS
# =============================================================================

primary_voxels     <- c("AUD", "SM", "VV")
exploratory_voxels <- c("MEM")
metabolites        <- c("G", "GLX")
corrections        <- c("UNC", "ATC")


# =============================================================================
# 3. HELPER FUNCTION FOR SEX x DIAGNOSIS INTERACTION
# =============================================================================
# Extends the standard MANOVA approach to include Diagnosis and Sex x Diagnosis
# as predictors. Univariate follow-ups extract the interaction term.

run_interaction_analysis <- function(data, metabolite, correction, voxels, label) {
  
  cols <- paste0(metabolite, "_", correction, "_", voxels)
  
  missing <- cols[!cols %in% names(data)]
  if (length(missing) > 0) {
    return(list(label = label, error = paste("Missing columns:", paste(missing, collapse = ", "))))
  }
  
  analysis_dat <- data %>%
    select(Sex, Subgroup_f, Age_DuringParticipation, all_of(cols)) %>%
    drop_na()
  
  n_total  <- nrow(analysis_dat)
  n_male   <- sum(analysis_dat$Sex == "Male")
  n_female <- sum(analysis_dat$Sex == "Female")
  
  if (n_total < 10 || n_male < 3 || n_female < 3) {
    return(list(label = label, error = paste("Insufficient sample size (n =", n_total, ")")))
  }
  
  outcome_matrix <- as.matrix(analysis_dat[, cols])
  
  # MANOVA with Sex x Diagnosis interaction
  manova_fit <- tryCatch(
    manova(outcome_matrix ~ Sex * Subgroup_f + Age_DuringParticipation,
           data = analysis_dat),
    error = function(e) NULL
  )
  
  if (is.null(manova_fit)) {
    return(list(label = label, error = "MANOVA failed"))
  }
  
  # Extract interaction term from MANOVA
  manova_summary <- tryCatch({
    s            <- summary(manova_fit, test = "Pillai")
    interaction_row <- s$stats["Sex:Subgroup_f", ]
    list(
      pillai = round(interaction_row["Pillai"], 3),
      F      = round(interaction_row["approx F"], 3),
      df1    = round(interaction_row["num Df"], 0),
      df2    = round(interaction_row["den Df"], 0),
      p      = round(interaction_row["Pr(>F)"], 4),
      sig    = ifelse(interaction_row["Pr(>F)"] < .05, "*", "")
    )
  }, error = function(e) {
    list(error_msg = paste("MANOVA extraction failed:", e$message))
  })
  
  # Univariate follow-ups — extract Sex x Diagnosis interaction
  uni_results <- list()
  for (col in cols) {
    vox <- str_extract(col, "[^_]+$")
    fit <- lm(as.formula(paste0("`", col, "` ~ Sex * Subgroup_f + Age_DuringParticipation")),
              data = analysis_dat)
    coefs <- summary(fit)$coefficients
    int_term <- "SexMale:Subgroup_fMCI"
    if (int_term %in% rownames(coefs)) {
      coef_int <- coefs[int_term, ]
      uni_results[[vox]] <- list(
        beta = round(coef_int["Estimate"], 4),
        t    = round(coef_int["t value"], 3),
        p    = round(coef_int["Pr(>|t|)"], 4),
        sig  = ifelse(coef_int["Pr(>|t|)"] < .05, "*", "")
      )
    }
  }
  
  list(
    label      = label,
    n          = n_total,
    n_male     = n_male,
    n_female   = n_female,
    manova     = manova_summary,
    univariate = uni_results,
    error      = NULL
  )
}

# Single voxel version for exploratory analyses
run_interaction_single <- function(data, metabolite, correction, voxel, label) {
  
  col <- paste0(metabolite, "_", correction, "_", voxel)
  if (!col %in% names(data)) {
    return(list(label = label, error = paste("Column not found:", col)))
  }
  
  analysis_dat <- data %>%
    select(Sex, Subgroup_f, Age_DuringParticipation, all_of(col)) %>%
    drop_na()
  
  n_total  <- nrow(analysis_dat)
  n_male   <- sum(analysis_dat$Sex == "Male")
  n_female <- sum(analysis_dat$Sex == "Female")
  
  if (n_total < 10 || n_male < 3 || n_female < 3) {
    return(list(label = label, error = paste("Insufficient sample size (n =", n_total, ")")))
  }
  
  fit      <- lm(as.formula(paste0("`", col, "` ~ Sex * Subgroup_f + Age_DuringParticipation")),
                 data = analysis_dat)
  coefs    <- summary(fit)$coefficients
  int_term <- "SexMale:Subgroup_fMCI"
  
  if (!int_term %in% rownames(coefs)) {
    return(list(label = label, error = "Interaction term not estimable"))
  }
  
  coef_int <- coefs[int_term, ]
  uni_list <- list()
  uni_list[[voxel]] <- list(
    beta = round(coef_int["Estimate"], 4),
    t    = round(coef_int["t value"], 3),
    p    = round(coef_int["Pr(>|t|)"], 4),
    sig  = ifelse(coef_int["Pr(>|t|)"] < .05, "*", "")
  )
  
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
# 4. RUN PRIMARY ANALYSES
# =============================================================================

results <- list()

for (met in metabolites) {
  for (corr in corrections) {
    label <- paste0(met, " — ", corr)
    results[[label]] <- run_interaction_analysis(older_all, met, corr, primary_voxels, label)
  }
}


# =============================================================================
# 5. RUN EXPLORATORY ANALYSES (MEM)
# =============================================================================

for (met in metabolites) {
  for (corr in corrections) {
    label <- paste0(met, " — ", corr, " [EXPLORATORY: MEM]")
    results[[label]] <- run_interaction_single(older_all, met, corr, "MEM", label)
  }
}


# =============================================================================
# 6. PRINT AND SAVE SUMMARY
# =============================================================================

summary_lines <- format_summary(
  results,
  title     = "Sex x Diagnosis Interaction in Older Adults (Wave 1)",
  predictor = "Sex x Diagnosis interaction (ref: Female, CN); Covariate: Age"
)
cat(paste(summary_lines, collapse = "\n"), "\n")

writeLines(summary_lines, "outputs/8_CN_vs_MCI_sex_summary.txt")
cat("\nSummary saved to outputs/8_CN_vs_MCI_sex_summary.txt\n")