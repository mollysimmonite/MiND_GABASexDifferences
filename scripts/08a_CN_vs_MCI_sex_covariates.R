# =============================================================================
# 08a_CN_vs_MCI_sex_covariates.R
# Sex x Diagnosis Interaction - Older Adults (with covariates)
#
# Author:       Molly Simmonite
# Date:         April 2026
# Project:      MiND Study - MRS Analysis
#
# Description:  Follow-up to 08_CN_vs_MCI_sex.R. Repeats all analyses adding
#               MRSTime (time of scan) and FMRI_ScannerUpgrade as covariates
#               to check robustness of primary findings.
#
# Input:        data/MRS_data_analysis.csv
# Output:       outputs/08a_CN_vs_MCI_sex_covariates_summary.txt
#
# Dependencies: helpers.R, tidyverse
# =============================================================================

source("scripts/helpers.R")

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


# =============================================================================
# 2. DEFINE ANALYSIS PARAMETERS
# =============================================================================

primary_voxels     <- c("AUD", "SM", "VV")
exploratory_voxels <- c("MEM")
metabolites        <- c("G", "GLX")
corrections        <- c("UNC", "ATC")
extra_covs         <- c("MRSTime", "FMRI_ScannerUpgrade")


# =============================================================================
# 3. HELPER FUNCTIONS FOR SEX x DIAGNOSIS INTERACTION WITH COVARIATES
# =============================================================================

run_interaction_analysis_covs <- function(data, metabolite, correction, voxels, label,
                                          extra_covariates = NULL) {
  
  cols <- paste0(metabolite, "_", correction, "_", voxels)
  
  missing <- cols[!cols %in% names(data)]
  if (length(missing) > 0) {
    return(list(label = label, error = paste("Missing columns:", paste(missing, collapse = ", "))))
  }
  
  fixed_terms <- c("Sex * Subgroup_f", "Age_DuringParticipation", extra_covariates)
  formula_str <- paste(fixed_terms, collapse = " + ")
  select_cols <- c("Sex", "Subgroup_f", "Age_DuringParticipation", extra_covariates, cols)
  
  analysis_dat <- data %>%
    select(all_of(select_cols)) %>%
    drop_na()
  
  n_total  <- nrow(analysis_dat)
  n_male   <- sum(analysis_dat$Sex == "Male")
  n_female <- sum(analysis_dat$Sex == "Female")
  
  if (n_total < 10 || n_male < 3 || n_female < 3) {
    return(list(label = label, error = paste("Insufficient sample size (n =", n_total, ")")))
  }
  
  outcome_matrix <- as.matrix(analysis_dat[, cols])
  
  manova_fit <- tryCatch(
    manova(as.formula(paste0("outcome_matrix ~ ", formula_str)),
           data = analysis_dat),
    error = function(e) NULL
  )
  
  if (is.null(manova_fit)) {
    return(list(label = label, error = "MANOVA failed"))
  }
  
  manova_summary <- tryCatch({
    s               <- summary(manova_fit, test = "Pillai")
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
  
  uni_results <- list()
  for (col in cols) {
    vox <- str_extract(col, "[^_]+$")
    fit <- lm(as.formula(paste0("`", col, "` ~ ", formula_str)),
              data = analysis_dat)
    coefs    <- summary(fit)$coefficients
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

run_interaction_single_covs <- function(data, metabolite, correction, voxel, label,
                                        extra_covariates = NULL) {
  
  col <- paste0(metabolite, "_", correction, "_", voxel)
  if (!col %in% names(data)) {
    return(list(label = label, error = paste("Column not found:", col)))
  }
  
  extra_str   <- if (!is.null(extra_covariates)) paste("+", paste(extra_covariates, collapse = " + ")) else ""
  select_cols <- c("Sex", "Subgroup_f", "Age_DuringParticipation", extra_covariates, col)
  
  analysis_dat <- data %>%
    select(all_of(select_cols)) %>%
    drop_na()
  
  n_total  <- nrow(analysis_dat)
  n_male   <- sum(analysis_dat$Sex == "Male")
  n_female <- sum(analysis_dat$Sex == "Female")
  
  if (n_total < 10 || n_male < 3 || n_female < 3) {
    return(list(label = label, error = paste("Insufficient sample size (n =", n_total, ")")))
  }
  
  fit      <- lm(as.formula(paste0("`", col, "` ~ Sex * Subgroup_f + Age_DuringParticipation ", extra_str)),
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
    label <- paste0(met, " - ", corr)
    results[[label]] <- run_interaction_analysis_covs(older_all, met, corr, primary_voxels, label,
                                                      extra_covariates = extra_covs)
  }
}


# =============================================================================
# 5. RUN EXPLORATORY ANALYSES (MEM)
# =============================================================================

for (met in metabolites) {
  for (corr in corrections) {
    label <- paste0(met, " - ", corr, " [EXPLORATORY: MEM]")
    results[[label]] <- run_interaction_single_covs(older_all, met, corr, "MEM", label,
                                                    extra_covariates = extra_covs)
  }
}


# =============================================================================
# 6. PRINT AND SAVE SUMMARY
# =============================================================================

summary_lines <- format_summary(
  results,
  title     = "Sex x Diagnosis Interaction in Older Adults (Wave 1) - with covariates",
  predictor = "Sex x Diagnosis interaction (ref: Female, CN); Covariates: Age, MRSTime, FMRI_ScannerUpgrade"
)
cat(paste(summary_lines, collapse = "\n"), "\n")

writeLines(summary_lines, "outputs/08a_CN_vs_MCI_sex_covariates_summary.txt")
cat("\nSummary saved to outputs/08a_CN_vs_MCI_sex_covariates_summary.txt\n")