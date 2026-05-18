# =============================================================================
# 1_create_dataset.R
# MRS Data Cleaning & Preparation
#
# Author:       Molly Simmonite
# Date:         April 2026
# Project:      MiND Study — MRS Analysis
#
# Description:  Loads raw MRS data, applies inclusion/exclusion criteria,
#               recodes variables, drops unused corrections (ATCG), and
#               outputs a clean dataset ready for analysis. Also produces
#               a GM fraction QC plot to inform outlier exclusion decisions.
#
# Input:        data/MiNDMaster_March26.csv
# Output:       data/MRS_data_filtered.csv
#               QC/GM_fraction_QC.png
#
# Dependencies: tidyverse, ggplot2
# =============================================================================

library(tidyverse)
library(ggplot2)

# =============================================================================
# 1. LOAD DATA
# =============================================================================

dat <- read_csv("data/MiNDMaster_March26.csv")

cat("Raw data dimensions:", nrow(dat), "rows x", ncol(dat), "cols\n")


# =============================================================================
# 2. INCLUSION / EXCLUSION FILTERS
# =============================================================================

dat_clean <- dat %>%
  
  # Keep only completed QC and enrollment
  filter(QC_Complete == "Completed",
         EnrollmentStatus == "Completed") %>%
  
  # Keep only CN (1) and MCI (2) — excludes uncertain (0), Alzheimer's (3), and any unexpected values
  filter(Subgroup %in% c(1, 2)) %>%
  
  # Exclude waves where participant was on Estrogen or GABApentin
  # (wave-level exclusion only — prior waves retained)
  filter(is.na(Estrogen)  | Estrogen  != 1) %>%
  filter(is.na(GABApentin) | GABApentin != 1)

cat("After exclusions:", nrow(dat_clean), "rows\n")
cat("Participants removed:", n_distinct(dat$SubNum) - n_distinct(dat_clean$SubNum), "\n")
cat("Rows removed:", nrow(dat) - nrow(dat_clean), "\n\n")


# =============================================================================
# 3. RECODE VARIABLES
# =============================================================================

dat_clean <- dat_clean %>%
  mutate(
    # Factors with meaningful labels
    Sex          = factor(Sex_M1,        levels = c(1, 2), labels = c("Male", "Female")),
    AgeCategory  = factor(AgeCategory_Young1, levels = c(1, 2), labels = c("Young", "Older")),
    Subgroup_f   = factor(Subgroup,      levels = c(1, 2), labels = c("CN", "MCI")),
    Wave_f       = factor(Wave_Num),
    
    # Ensure SubNum is treated as an ID (character), not numeric
    SubNum       = as.character(SubNum),
    
    # Recode MothersEd 999 to NA (missing data code)
    MothersEd    = ifelse(MothersEd == 999, NA, MothersEd)
  )


# =============================================================================
# 4. DROP ATCG COLUMNS
# (group correction — not to be used in analyses)
# =============================================================================

atcg_cols <- names(dat_clean)[str_detect(names(dat_clean), "_ATCG_")]
cat("Dropping", length(atcg_cols), "ATCG columns\n\n")

dat_clean <- dat_clean %>%
  select(-all_of(atcg_cols))


# =============================================================================
# 5. GM FRACTION QC — INSPECT FOR OUTLIERS
# =============================================================================
# Not hard-excluding yet — visual inspection first.
# Primary voxels only (L/R pairs + averages).

primary_voxels <- c("LAUD", "RAUD", "AUD", "LSM", "RSM", "SM", "LVV", "RVV", "VV")

gm_cols <- names(dat_clean)[str_detect(names(dat_clean), "^GM_FRA_")]

# Subset to primary voxels only for the summary
gm_primary <- gm_cols[str_extract(gm_cols, "(?<=GM_FRA_).+") %in% primary_voxels]

# Summary table
cat("=== GM Fraction Summary (primary voxels) ===\n")
dat_clean %>%
  select(all_of(gm_primary)) %>%
  pivot_longer(everything(), names_to = "Voxel", values_to = "GM_FRA") %>%
  mutate(Voxel = str_remove(Voxel, "GM_FRA_")) %>%
  group_by(Voxel) %>%
  summarise(
    n       = sum(!is.na(GM_FRA)),
    mean    = round(mean(GM_FRA, na.rm = TRUE), 3),
    sd      = round(sd(GM_FRA,   na.rm = TRUE), 3),
    min     = round(min(GM_FRA,  na.rm = TRUE), 3),
    max     = round(max(GM_FRA,  na.rm = TRUE), 3),
    n_below_0.4 = sum(GM_FRA < 0.4, na.rm = TRUE)  # rough low-GM flag
  ) %>%
  print(n = Inf)

# Boxplots of GM fraction per voxel
p_gm <- dat_clean %>%
  select(SubNum, Wave_f, all_of(gm_primary)) %>%
  pivot_longer(all_of(gm_primary), names_to = "Voxel", values_to = "GM_FRA") %>%
  mutate(Voxel = str_remove(Voxel, "GM_FRA_"),
         Voxel = factor(Voxel, levels = primary_voxels)) %>%
  ggplot(aes(x = Voxel, y = GM_FRA)) +
  geom_boxplot(outlier.colour = "red", outlier.size = 2, fill = "steelblue", alpha = 0.5) +
  geom_jitter(width = 0.15, alpha = 0.3, size = 1) +
  labs(title = "GM Fraction by Voxel",
       subtitle = "Red points = potential outliers. Inspect before deciding on exclusion threshold.",
       x = "Voxel", y = "GM Fraction") +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

dir.create("QC", recursive = TRUE, showWarnings = FALSE)
ggsave("QC/GM_fraction_QC.png", p_gm, width = 10, height = 5, dpi = 150)
cat("\nGM fraction plot saved to QC/GM_fraction_QC.png\n")


# =============================================================================
# 6. QUICK SAMPLE SUMMARY
# =============================================================================

cat("\n=== Sample Summary After Cleaning ===\n")

dat_clean %>%
  distinct(SubNum, AgeCategory, Subgroup_f, Sex) %>%
  count(AgeCategory, Subgroup_f, Sex) %>%
  print()

cat("\nWaves per participant:\n")
dat_clean %>%
  count(SubNum) %>%
  count(n, name = "n_participants") %>%
  rename(n_waves = n) %>%
  print()


# =============================================================================
# 7. SAVE CLEAN DATA
# =============================================================================

write_csv(dat_clean, "data/MRS_data_filtered.csv")
cat("\nClean data saved to data/MRS_data_filtered.csv\n")
cat("Final dimensions:", nrow(dat_clean), "rows x", ncol(dat_clean), "cols\n")