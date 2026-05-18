# =============================================================================
# 3_MRS_flags.R
# Apply MRS Sanity Flags — Set Implausible Values to NA
#
# Author:       Molly Simmonite
# Date:         April 2026
# Project:      MiND Study — MRS Analysis
#
# Description:  Takes the flagged values identified in 02_MRS_sanity_check.R
#               and sets those specific values to NA in the clean dataset.
#               Also sets MRS values to NA for any voxel where GM_FRA < 0.2,
#               indicating likely poor voxel placement.
#               Exclusions are value-level and voxel-level only:
#                 - A bad Glx value does NOT affect GABA from the same scan
#                 - A bad value in one voxel does NOT affect other voxels
#                 - No whole participants are excluded based on MRS flags
#
# Input:        data/MRS_data_filtered.csv
#               QC/MRS_QC_flags_review.csv
# Output:       data/MRS_data_analysis.csv
#               QC/MRS_exclusions_applied.csv
#               QC/MRS_available_data_summary.csv
#
# Dependencies: tidyverse
# =============================================================================

library(tidyverse)

# =============================================================================
# 1. LOAD DATA
# =============================================================================

dat   <- read_csv("data/MRS_data_filtered.csv")
flags <- read_csv("QC/MRS_QC_flags_review.csv")

cat("Clean data:", nrow(dat), "rows\n")
cat("Flagged values to NA:", nrow(flags), "\n\n")


# =============================================================================
# 2. SET FLAGGED VALUES TO NA
# =============================================================================
# For each flagged row, set the specific Variable column to NA for that
# SubNum x Wave_Num combination. All other columns untouched.

dat_flagged <- dat

for (i in seq_len(nrow(flags))) {
  sub  <- flags$SubNum[i]
  wave <- flags$Wave_Num[i]
  var  <- flags$Variable[i]
  
  # Only proceed if the variable exists in the dataset
  if (var %in% names(dat_flagged)) {
    dat_flagged[dat_flagged$SubNum == sub & dat_flagged$Wave_Num == wave, var] <- NA
  } else {
    cat("Warning: variable", var, "not found in dataset — skipping\n")
  }
}


# =============================================================================
# 3. SET MRS VALUES TO NA WHERE GM FRACTION < 0.2 (POOR VOXEL PLACEMENT)
# =============================================================================
# A GM fraction below 0.2 suggests the voxel was poorly placed and contains
# very little grey matter. This is a data quality issue (voxel placement),
# not a biological implausibility. All MRS metabolite values for that specific
# voxel are set to NA; other voxels and all other columns are untouched.
#
# Average voxels (AUD, SM, VV) are derived from bilateral pairs, so if either
# constituent voxel has GM_FRA < 0.2, the average is also set to NA.
# Mapping: LAUD/RAUD -> AUD, LSM/RSM -> SM, LVV/RVV -> VV

gm_threshold <- 0.2

# Bilateral voxels to check, and their corresponding average voxel
bilateral_map <- list(
  LAUD = "AUD",
  RAUD = "AUD",
  LSM  = "SM",
  RSM  = "SM",
  LVV  = "VV",
  RVV  = "VV"
)

cat("=== GM Fraction Exclusions (GM_FRA <", gm_threshold, ") ===\n")

for (vox in names(bilateral_map)) {
  avg_vox <- bilateral_map[[vox]]
  gm_col  <- paste0("GM_FRA_", vox)
  
  if (!gm_col %in% names(dat_flagged)) {
    cat("Warning: ", gm_col, "not found — skipping\n")
    next
  }
  
  # MRS columns for the bilateral voxel itself
  mrs_cols_vox <- names(dat_flagged)[str_detect(names(dat_flagged), paste0("^(G|GLX)_.*_", vox, "$")) &
                                       !str_detect(names(dat_flagged), "^(GM|WM|CSF)_FRA")]
  
  # MRS columns for the corresponding average voxel
  mrs_cols_avg <- names(dat_flagged)[str_detect(names(dat_flagged), paste0("^(G|GLX)_.*_", avg_vox, "$")) &
                                       !str_detect(names(dat_flagged), "^(GM|WM|CSF)_FRA")]
  
  low_gm_rows <- !is.na(dat_flagged[[gm_col]]) & dat_flagged[[gm_col]] < gm_threshold
  n_affected  <- sum(low_gm_rows)
  
  if (n_affected > 0) {
    dat_flagged[low_gm_rows, mrs_cols_vox] <- NA
    dat_flagged[low_gm_rows, mrs_cols_avg] <- NA
    cat(vox, ":", n_affected, "rows with GM_FRA <", gm_threshold,
        "-> NAs set for", vox, "and average voxel", avg_vox, "\n")
  }
}
cat("\n")


# =============================================================================
# 4. SUMMARY OF CHANGES
# =============================================================================

cat("=== Total NAs introduced (spectral flags + GM fraction exclusions) ===\n")

# Count NAs before and after for MRS columns only
mrs_cols <- names(dat)[str_detect(names(dat), "^(G|GLX)_") &
                         !str_detect(names(dat), "^(GM|WM|CSF)_FRA")]

before <- dat         %>% select(all_of(mrs_cols)) %>% summarise(across(everything(), ~ sum(is.na(.x)))) %>% pivot_longer(everything(), values_to = "na_before")
after  <- dat_flagged %>% select(all_of(mrs_cols)) %>% summarise(across(everything(), ~ sum(is.na(.x)))) %>% pivot_longer(everything(), values_to = "na_after")

changes <- left_join(before, after, by = "name") %>%
  mutate(new_nas = na_after - na_before) %>%
  filter(new_nas > 0) %>%
  arrange(desc(new_nas))

print(changes, n = Inf)
cat("\nTotal new NAs introduced:", sum(changes$new_nas), "\n")

# Save updated summary table of exclusions we introduced
# Compare before (dat) and after (dat_flagged) to find only newly introduced NAs

dat_long_before <- dat %>%
  select(SubNum, Wave_Num, all_of(mrs_cols)) %>%
  pivot_longer(all_of(mrs_cols), names_to = "Variable", values_to = "Value_before")

dat_long_after <- dat_flagged %>%
  select(SubNum, Wave_Num, all_of(mrs_cols)) %>%
  pivot_longer(all_of(mrs_cols), names_to = "Variable", values_to = "Value_after")

mrs_na_summary <- left_join(dat_long_before, dat_long_after, by = c("SubNum", "Wave_Num", "Variable")) %>%
  # Keep only rows where value was NOT NA before but IS NA after
  filter(!is.na(Value_before) & is.na(Value_after)) %>%
  mutate(
    Metabolite = str_extract(Variable, "^[^_]+"),
    Correction = str_extract(Variable, "(?<=_)[^_]+(?=_)"),
    Voxel      = str_extract(Variable, "[^_]+$")
  ) %>%
  select(SubNum, Wave_Num, Variable, Metabolite, Correction, Voxel)

write_csv(mrs_na_summary, "QC/MRS_exclusions_applied.csv")
cat("\nExclusions table saved to QC/MRS_exclusions_applied.csv\n",
    "(", nrow(mrs_na_summary), "values newly set to NA across",
    n_distinct(mrs_na_summary$SubNum), "participants)\n")


# =============================================================================
# 5. SAVE
# =============================================================================

write_csv(dat_flagged, "data/MRS_data_analysis.csv")
cat("\nSaved to data/MRS_data_analysis.csv\n")


# =============================================================================
# 6. AVAILABLE DATA SUMMARY TABLE
# =============================================================================
# Count of non-NA observations per voxel x metabolite x correction,
# broken down by group (AgeCategory x Wave x Subgroup x Sex).
# Primary corrections only (UNC and ATC).
# Primary voxels + MEM and LPV.

summary_voxels <- c("LAUD", "RAUD", "AUD", "LSM", "RSM", "SM", "LVV", "RVV", "VV", "MEM", "LPV")
summary_corrections <- c("UNC", "ATC")

# Build group label: Young CN is wave-independent; older adults split by wave
dat_flagged <- dat_flagged %>%
  mutate(
    Group = case_when(
      AgeCategory == "Young"               ~ paste("Young CN", Sex),
      AgeCategory == "Older"               ~ paste("Older W", Wave_Num, Subgroup_f, Sex),
      TRUE                                 ~ NA_character_
    )
  )

# Select relevant MRS columns
summary_cols <- names(dat_flagged)[
  str_detect(names(dat_flagged), paste0("^(G|GLX)_(", paste(summary_corrections, collapse = "|"), ")_(", paste(summary_voxels, collapse = "|"), ")$"))
]

# Pivot to long and count non-NA observations per group x voxel x metabolite x correction
avail_summary <- dat_flagged %>%
  select(SubNum, Group, all_of(summary_cols)) %>%
  pivot_longer(all_of(summary_cols), names_to = "Variable", values_to = "Value") %>%
  mutate(
    Metabolite  = str_extract(Variable, "^[^_]+"),
    Correction  = str_extract(Variable, "(?<=_)[^_]+(?=_)"),
    Voxel       = str_extract(Variable, "[^_]+$"),
    Voxel       = factor(Voxel, levels = summary_voxels)
  ) %>%
  filter(!is.na(Group)) %>%
  group_by(Group, Metabolite, Correction, Voxel) %>%
  summarise(n_available = sum(!is.na(Value)), .groups = "drop") %>%
  pivot_wider(names_from = Group, values_from = n_available) %>%
  arrange(Metabolite, Correction, Voxel)

write_csv(avail_summary, "QC/MRS_available_data_summary.csv")
cat("\nAvailable data summary saved to QC/MRS_available_data_summary.csv\n")
print(avail_summary, n = Inf)