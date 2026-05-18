# =============================================================================
# 4a_young_sex_differences_covariates.R
# Sex Differences in GABA and Glx — Young Adults (with covariates)
#
# Author:       Molly Simmonite
# Date:         April 2026
# Project:      MiND Study — MRS Analysis
#
# Description:  Follow-up to 4_young_sex_differences.R. Repeats all analyses
#               with additional covariates:
#                 - MRSTime (time of scan)
#                 - FMRI_ScannerUpgrade (scanner change)
#               Plus an exploratory analysis replacing Sex with a 3-level
#               hormonal factor (Male, Female-Follicular, Female-Luteal) using
#               MenstrualStage. Note: 15/29 young females have NA for
#               MenstrualStage — this analysis is based on a reduced sample
#               and should be interpreted cautiously.
#
# Input:        data/MRS_data_analysis.csv
# Output:       outputs/4a_young_sex_differences_covariates_summary.txt
#
# Dependencies: helpers.R, tidyverse
# =============================================================================

source("scripts/helpers.R")

dir.create("outputs", recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# 1. LOAD AND PREPARE DATA
# =============================================================================

dat <- read_csv("data/MRS_data_analysis.csv")

young <- dat %>%
  filter(AgeCategory == "Young",
         Subgroup_f  == "CN",
         Wave_Num    == 1)

cat("Young CN sample: n =", nrow(young), "\n")
cat("Sex breakdown:\n")
print(table(young$Sex))
cat("\n")


# =============================================================================
# 2. DEFINE ANALYSIS PARAMETERS
# =============================================================================

primary_voxels     <- c("AUD", "SM", "VV")
exploratory_voxels <- c("LPV")
metabolites        <- c("G", "GLX")
corrections        <- c("UNC", "ATC")
extra_covs         <- c("MRSTime", "FMRI_ScannerUpgrade")


# =============================================================================
# 3. RUN PRIMARY ANALYSES WITH COVARIATES
# =============================================================================

results <- list()

for (met in metabolites) {
  for (corr in corrections) {
    label <- paste0(met, " — ", corr)
    results[[label]] <- run_manova_univariate(young, met, corr, primary_voxels, label,
                                              extra_covariates = extra_covs)
  }
}


# =============================================================================
# 4. RUN EXPLORATORY ANALYSES WITH COVARIATES (LPV)
# =============================================================================

for (met in metabolites) {
  for (corr in corrections) {
    label <- paste0(met, " — ", corr, " [EXPLORATORY: LPV]")
    results[[label]] <- run_univariate_single(young, met, corr, "LPV", label,
                                              extra_covariates = extra_covs)
  }
}


# =============================================================================
# 5. EXPLORATORY: HORMONAL FACTOR ANALYSIS
# =============================================================================
# Replace binary Sex with a 3-level hormonal factor:
#   0 = Male, 1 = Female-Follicular, 2 = Female-Luteal
# MenstrualStage NAs result in reduced female sample (n = 14 follicular/luteal
# vs 29 total). Interpret with caution.

young_hormonal <- young %>%
  mutate(
    HormonalGroup = case_when(
      Sex == "Male"                    ~ "Male",
      Sex == "Female" & MenstrualStage == 1 ~ "Female-Follicular",
      Sex == "Female" & MenstrualStage == 2 ~ "Female-Luteal",
      TRUE                             ~ NA_character_
    ),
    HormonalGroup = factor(HormonalGroup,
                           levels = c("Male", "Female-Follicular", "Female-Luteal"))
  ) %>%
  filter(!is.na(HormonalGroup))

cat("Hormonal factor sample (excluding NA MenstrualStage):\n")
print(table(young_hormonal$HormonalGroup))
cat("\n")

hormonal_results <- list()

for (met in metabolites) {
  for (corr in corrections) {
    for (vox in primary_voxels) {
      col   <- paste0(met, "_", corr, "_", vox)
      label <- paste0(met, " — ", corr, " — ", vox, " [HORMONAL FACTOR]")
      
      if (!col %in% names(young_hormonal)) next
      
      analysis_dat <- young_hormonal %>%
        select(HormonalGroup, Age_DuringParticipation, MRSTime, FMRI_ScannerUpgrade, all_of(col)) %>%
        drop_na()
      
      if (nrow(analysis_dat) < 10) next
      
      fit      <- lm(as.formula(paste0("`", col, "` ~ HormonalGroup + Age_DuringParticipation + MRSTime + FMRI_ScannerUpgrade")),
                     data = analysis_dat)
      coefs    <- summary(fit)$coefficients
      
      # Extract Male vs Female-Follicular and Male vs Female-Luteal contrasts
      extract_coef <- function(coefs, term) {
        if (!term %in% rownames(coefs)) return(NULL)
        row <- coefs[term, ]
        list(
          beta = round(row["Estimate"], 4),
          t    = round(row["t value"], 3),
          p    = round(row["Pr(>|t|)"], 4),
          sig  = ifelse(row["Pr(>|t|)"] < .05, "*", "")
        )
      }
      
      foll <- extract_coef(coefs, "HormonalGroupFemale-Follicular")
      lut  <- extract_coef(coefs, "HormonalGroupFemale-Luteal")
      
      hormonal_results[[label]] <- list(
        label      = label,
        n          = nrow(analysis_dat),
        follicular = foll,
        luteal     = lut
      )
    }
  }
}


# =============================================================================
# 6. FORMAT AND SAVE SUMMARY
# =============================================================================

summary_lines <- format_summary(
  results,
  title     = "Sex Differences in Young CN Adults (with covariates: MRSTime, ScannerUpgrade)",
  predictor = "Sex (ref = Female); Covariates: Age, MRSTime, FMRI_ScannerUpgrade"
)

# Add hormonal factor results
summary_lines <- c(summary_lines,
                   "================================================================",
                   "EXPLORATORY: Hormonal Factor Analysis (Male / Female-Follicular / Female-Luteal)",
                   "Reference = Male; Covariates: Age, MRSTime, FMRI_ScannerUpgrade",
                   paste0("Note: reduced sample due to MenstrualStage NAs (n females = ",
                          sum(young_hormonal$HormonalGroup != "Male"), ")"),
                   "* = p < .05",
                   "================================================================",
                   ""
)

for (res in hormonal_results) {
  summary_lines <- c(summary_lines, paste0("--- ", res$label, " ---"))
  summary_lines <- c(summary_lines, paste0("  n = ", res$n))
  
  if (!is.null(res$follicular)) {
    f <- res$follicular
    summary_lines <- c(summary_lines, paste0(
      "  Female-Follicular vs Male: beta = ", f$beta,
      ", t = ", f$t, ", p ", format_p(f$p), " ", f$sig
    ))
  }
  if (!is.null(res$luteal)) {
    l <- res$luteal
    summary_lines <- c(summary_lines, paste0(
      "  Female-Luteal vs Male:     beta = ", l$beta,
      ", t = ", l$t, ", p ", format_p(l$p), " ", l$sig
    ))
  }
  summary_lines <- c(summary_lines, "")
}

cat(paste(summary_lines, collapse = "\n"), "\n")

writeLines(summary_lines, "outputs/4a_young_sex_differences_covariates_summary.txt")
cat("\nSummary saved to outputs/4a_young_sex_differences_covariates_summary.txt\n")


# =============================================================================
# 7. PLOTS — 4 GROUP RAW VALUES
# =============================================================================
# Visualise raw GABA and Glx values across four groups:
#   Male, Female-NA, Female-Follicular, Female-Luteal
# One PDF page per metabolite x correction, primary voxels as facets.
# Including NA females allows visual inspection of whether they resemble
# one of the other groups (e.g. if similar to luteal, may suggest
# hormonal contraception use).

dir.create("plots", recursive = TRUE, showWarnings = FALSE)

# Build 4-group factor including NA females
young_4groups <- young %>%
  mutate(
    HormonalGroup4 = case_when(
      Sex == "Male"                         ~ "Male",
      Sex == "Female" & MenstrualStage == 1 ~ "Female-Follicular",
      Sex == "Female" & MenstrualStage == 2 ~ "Female-Luteal",
      Sex == "Female" & is.na(MenstrualStage) ~ "Female-NA",
      TRUE                                  ~ NA_character_
    ),
    HormonalGroup4 = factor(HormonalGroup4,
                            levels = c("Male", "Female-Follicular", "Female-Luteal", "Female-NA"))
  )

cat("\n4-group sample breakdown:\n")
print(table(young_4groups$HormonalGroup4))
cat("\n")

group_colours <- c(
  "Male"               = "#4C72B0",
  "Female-Follicular"  = "#DD8452",
  "Female-Luteal"      = "#55A868",
  "Female-NA"          = "#C44E52"
)

pdf("plots/4a_hormonal_groups_raw.pdf", width = 10, height = 5)

for (met in metabolites) {
  for (corr in corrections) {
    
    cols <- paste0(met, "_", corr, "_", primary_voxels)
    cols_present <- cols[cols %in% names(young_4groups)]
    if (length(cols_present) == 0) next
    
    plot_dat <- young_4groups %>%
      select(HormonalGroup4, all_of(cols_present)) %>%
      pivot_longer(all_of(cols_present), names_to = "Variable", values_to = "Value") %>%
      mutate(Voxel = str_extract(Variable, "[^_]+$"),
             Voxel = factor(Voxel, levels = primary_voxels)) %>%
      filter(!is.na(Value), !is.na(HormonalGroup4))
    
    p <- ggplot(plot_dat, aes(x = HormonalGroup4, y = Value, fill = HormonalGroup4)) +
      geom_boxplot(outlier.shape = NA, alpha = 0.5) +
      geom_jitter(aes(colour = HormonalGroup4), width = 0.15, alpha = 0.6, size = 1.5) +
      scale_fill_manual(values = group_colours)   +
      scale_colour_manual(values = group_colours) +
      facet_wrap(~ Voxel, scales = "free_y") +
      labs(
        title    = paste0(met, " — ", corr, " correction"),
        subtitle = "Raw values by hormonal group. Female-NA: MenstrualStage not recorded.",
        x        = NULL,
        y        = "Value"
      ) +
      theme_minimal(base_size = 12) +
      theme(
        axis.text.x  = element_text(angle = 30, hjust = 1),
        legend.position = "none"
      )
    
    print(p)
  }
}

dev.off()
cat("Plots saved to plots/4a_hormonal_groups_raw.pdf\n")