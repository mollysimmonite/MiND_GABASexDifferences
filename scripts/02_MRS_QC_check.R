# =============================================================================
# 2_MRS_QC_check.R
# MRS Data Sanity Check
#
# Author:       Molly Simmonite
# Date:         April 2026
# Project:      MiND Study — MRS Analysis
#
# Description:  Checks MRS metabolite values (GABA and Glx) for implausible
#               values, including negatives, zeros, and extreme values (> 5 MAD
#               from the voxel median). Uses MAD rather than SD to avoid
#               masking where large outliers inflate the mean/SD. Produces a
#               summary table of flagged values and boxplots per voxel.
#               Checks all correction types (UNC, CSF, TC, ATC, CR).
#
# Input:        data/MRS_data_filtered.csv
# Output:       QC/MRS_QC_flags_review.csv          — flagged values for inspection
#               QC/MRS_sanity_plots.pdf          — boxplots per metabolite x voxel
#
# Dependencies: tidyverse, ggplot2
# =============================================================================

library(tidyverse)
library(ggplot2)

# =============================================================================
# 1. LOAD CLEAN DATA
# =============================================================================

dat <- read_csv("data/MRS_data_filtered.csv")

cat("Loaded:", nrow(dat), "rows\n\n")


# =============================================================================
# 2. RESHAPE MRS COLUMNS TO LONG FORMAT
# =============================================================================
# Columns follow the pattern: METABOLITE_CORRECTION_VOXEL
# e.g. G_UNC_LAUD, GLX_ATC_RSM
# We pivot to long so we can check all values in one pass.

# Identify all MRS columns (G_ and GLX_ prefixes, excluding GM/WM/CSF fractions)
mrs_cols <- names(dat)[str_detect(names(dat), "^(G|GLX)_") &
                         !str_detect(names(dat), "^(GM|WM|CSF)_FRA")]

dat_long <- dat %>%
  select(SubNum, Wave_Num, AgeCategory, Subgroup_f, Sex, all_of(mrs_cols)) %>%
  pivot_longer(
    cols      = all_of(mrs_cols),
    names_to  = "Variable",
    values_to = "Value"
  ) %>%
  # Parse the variable name into components
  mutate(
    Metabolite  = str_extract(Variable, "^[^_]+"),           # G or GLX
    Correction  = str_extract(Variable, "(?<=_)[^_]+(?=_)"), # UNC, CSF, TC, ATC, CR
    Voxel       = str_extract(Variable, "[^_]+$")            # LAUD, RSM, etc.
  )


# =============================================================================
# 3. FLAG IMPLAUSIBLE VALUES
# =============================================================================
# Flag if:
#   (a) Value is negative
#   (b) Value is exactly zero
#   (c) Value is > 5 MAD from the voxel median (robust to extreme values;
#       using MAD rather than SD avoids masking where huge outliers inflate
#       the mean/SD and make other bad values look less extreme)

dat_flagged <- dat_long %>%
  group_by(Metabolite, Correction, Voxel) %>%
  mutate(
    voxel_median = median(Value, na.rm = TRUE),
    voxel_mad    = mad(Value,    na.rm = TRUE),
    flag_negative = !is.na(Value) & Value < 0,
    flag_zero     = !is.na(Value) & Value == 0,
    flag_extreme  = !is.na(Value) & abs(Value - voxel_median) > 5 * voxel_mad,
    any_flag      = flag_negative | flag_zero | flag_extreme
  ) %>%
  ungroup()

# Summary of flags
cat("=== Flag Summary ===\n")
cat("Negative values:     ", sum(dat_flagged$flag_negative, na.rm = TRUE), "\n")
cat("Zero values:         ", sum(dat_flagged$flag_zero,     na.rm = TRUE), "\n")
cat("Extreme (>5 MAD):    ", sum(dat_flagged$flag_extreme,  na.rm = TRUE), "\n")
cat("Total flagged:       ", sum(dat_flagged$any_flag,       na.rm = TRUE), "\n\n")


# =============================================================================
# 4. EXPORT FLAGGED VALUES TABLE
# =============================================================================

flags_out <- dat_flagged %>%
  filter(any_flag) %>%
  select(SubNum, Wave_Num, Variable, Metabolite, Correction, Voxel,
         Value, voxel_median, voxel_mad,
         flag_negative, flag_zero, flag_extreme) %>%
  mutate(across(c(voxel_median, voxel_mad), ~ round(.x, 4)),
         Value = round(Value, 4)) %>%
  arrange(Metabolite, Correction, Voxel, SubNum)

dir.create("QC", recursive = TRUE, showWarnings = FALSE)

write_csv(flags_out, "QC/MRS_QC_flags_review.csv")
cat("Flagged values saved to QC/MRS_QC_flags_review.csv\n")
cat("(", nrow(flags_out), "flagged observations across", n_distinct(flags_out$SubNum), "participants)\n\n")


# =============================================================================
# 5. BOXPLOTS PER METABOLITE x CORRECTION x VOXEL
# =============================================================================
# One PDF page per metabolite x correction combination.
# Flagged points highlighted in red.

# Only plot primary voxels to keep it readable
primary_voxels <- c("LAUD", "RAUD", "AUD", "LSM", "RSM", "SM", "LVV", "RVV", "VV")

plot_dat <- dat_flagged %>%
  filter(Voxel %in% primary_voxels) %>%
  mutate(Voxel = factor(Voxel, levels = primary_voxels))

# Get all metabolite x correction combos present in the data
combos <- plot_dat %>%
  distinct(Metabolite, Correction) %>%
  arrange(Metabolite, Correction)

pdf("QC/MRS_sanity_plots.pdf", width = 11, height = 6)

for (i in seq_len(nrow(combos))) {
  met  <- combos$Metabolite[i]
  corr <- combos$Correction[i]
  
  plot_data <- plot_dat %>%
    filter(Metabolite == met, Correction == corr)
  
  # Set y axis upper limit to median + 10 MAD across all voxels in this combo
  # so extreme flagged values don't compress the bulk of the data
  # Points above the limit are still shown as red dots at the ceiling
  y_max <- median(plot_data$Value, na.rm = TRUE) + 10 * mad(plot_data$Value, na.rm = TRUE)
  y_min <- min(0, min(plot_data$Value, na.rm = TRUE))  # include negatives if present
  
  p <- plot_data %>%
    ggplot(aes(x = Voxel, y = Value)) +
    geom_boxplot(outlier.shape = NA, fill = "steelblue", alpha = 0.4) +
    geom_jitter(aes(colour = any_flag), width = 0.2, alpha = 0.5, size = 1.5) +
    scale_colour_manual(values = c("FALSE" = "grey40", "TRUE" = "red"),
                        labels = c("FALSE" = "OK", "TRUE" = "Flagged"),
                        name   = NULL) +
    coord_cartesian(ylim = c(y_min, y_max)) +  # clip without dropping data points
    labs(
      title    = paste0(met, " — ", corr, " correction"),
      subtitle = "Red = negative, zero, or > 5 MAD from voxel median. Y-axis capped at median + 10 MAD.",
      x        = "Voxel",
      y        = "Value"
    ) +
    theme_minimal(base_size = 13) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  print(p)
}

dev.off()
cat("Plots saved to QC/MRS_sanity_plots.pdf\n")