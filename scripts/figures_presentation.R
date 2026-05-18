# =============================================================================
# figures_presentation.R
# Generate Figures for Collaborator Presentation
#
# Author:       Molly Simmonite
# Date:         April 2026
# Project:      MiND Study - MRS Analysis
#
# Description:  Generates all figures needed for collaborator presentation.
#               Uses estimated marginal means (emmeans) for cross-sectional
#               plots and model-predicted trajectories for longitudinal plots.
#               All figures saved to plots/presentation/
#
# Figures:
#   fig1_GABA_older_CN.pdf       — Slide 3: GABA emmeans, older CN, by sex
#   fig2_GLX_older_CN.pdf        — Slide 5: GLX emmeans, older CN, by sex
#   fig3_GABA_diagnosis.pdf      — Slide 6: GABA emmeans, older adults, by sex x diagnosis
#   fig4_longitudinal.pdf        — Slide 7: GABA trajectories over time by sex
#   fig5_GABA_young.pdf          — Slide 8: GABA emmeans, young CN, by sex
#
# Dependencies: tidyverse, lme4, lmerTest, emmeans
# =============================================================================

source("scripts/helpers.R")
library(lme4)
library(lmerTest)
library(emmeans)

dir.create("plots/presentation", recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# SHARED AESTHETICS
# =============================================================================

sex_colours <- c("Female" = "#E07B8A", "Male" = "#5B8DB8")
sex_fills   <- c("Female" = "#E07B8A44", "Male" = "#5B8DB844")

voxel_labels <- c("AUD" = "Auditory\n(AUD)", "SM" = "Sensorimotor\n(SM)", "VV" = "Visual\n(VV)")

theme_presentation <- function() {
  theme_bw(base_size = 14) +
    theme(
      panel.grid.minor   = element_blank(),
      panel.grid.major.x = element_blank(),
      strip.background   = element_rect(fill = "grey92", colour = "grey70"),
      strip.text         = element_text(face = "bold", size = 13),
      axis.title         = element_text(size = 13),
      axis.text          = element_text(size = 11),
      legend.title       = element_blank(),
      legend.text        = element_text(size = 12),
      legend.position    = "top",
      plot.title         = element_text(face = "bold", size = 15),
      plot.subtitle      = element_text(size = 12, colour = "grey40")
    )
}

# Helper to extract emmeans from lm model
get_emmeans <- function(data, col, group_var = "Sex", covariates = "Age_DuringParticipation") {
  formula_str <- paste0("`", col, "` ~ ", group_var, " + ", covariates)
  fit         <- lm(as.formula(formula_str), data = data)
  em          <- emmeans(fit, specs = group_var)
  as.data.frame(em)
}

# Helper to add significance annotation
sig_label <- function(p) {
  if (p < .001) return("***")
  if (p < .01)  return("**")
  if (p < .05)  return("*")
  return("ns")
}


# =============================================================================
# LOAD DATA
# =============================================================================

dat   <- read_csv("data/MRS_data_analysis.csv")
dates <- read_csv("data/Participation_Dates.csv")

older_CN  <- dat %>% filter(AgeCategory == "Older", Subgroup_f == "CN",  Wave_Num == 1)
older_MCI <- dat %>% filter(AgeCategory == "Older", Subgroup_f == "MCI", Wave_Num == 1)
older_all <- dat %>% filter(AgeCategory == "Older", Subgroup_f %in% c("CN", "MCI"), Wave_Num == 1) %>%
  mutate(Subgroup_f = factor(Subgroup_f, levels = c("CN", "MCI")))
young_CN  <- dat %>% filter(AgeCategory == "Young", Subgroup_f == "CN",  Wave_Num == 1)


# =============================================================================
# FIG 1: GABA OLDER CN — ESTIMATED MARGINAL MEANS BY SEX
# =============================================================================

voxels      <- c("AUD", "SM", "VV")
corrections <- c("UNC", "ATC")

emm_list <- list()
for (corr in corrections) {
  for (vox in voxels) {
    col    <- paste0("G_", corr, "_", vox)
    emm    <- get_emmeans(older_CN %>% select(Sex, Age_DuringParticipation, all_of(col)) %>% drop_na(),
                          col)
    emm$correction <- corr
    emm$voxel      <- vox
    emm_list[[paste0(corr, "_", vox)]] <- emm
  }
}

emm_GABA_CN <- bind_rows(emm_list) %>%
  mutate(
    voxel      = factor(voxel, levels = voxels, labels = voxel_labels),
    correction = factor(correction, levels = c("UNC", "ATC"),
                        labels = c("Uncorrected (UNC)", "Alpha Tissue Corrected (ATC)"))
  )

# p-values from primary script for annotation
gaba_cn_pvals <- tribble(
  ~correction, ~voxel, ~p,
  "Uncorrected (UNC)", "Auditory\n(AUD)",      0.0001,
  "Uncorrected (UNC)", "Sensorimotor\n(SM)",   0.0012,
  "Uncorrected (UNC)", "Visual\n(VV)",         0.0001,
  "Alpha Tissue Corrected (ATC)", "Auditory\n(AUD)",    0.0018,
  "Alpha Tissue Corrected (ATC)", "Sensorimotor\n(SM)", 0.0723,
  "Alpha Tissue Corrected (ATC)", "Visual\n(VV)",       0.0001
) %>%
  mutate(sig = map_chr(p, sig_label))

# y position for voxel asterisk annotations
y_max <- emm_GABA_CN %>%
  group_by(correction, voxel) %>%
  summarise(y_pos = max(emmean + SE) * 1.12, .groups = "drop")

annot_GABA_CN <- left_join(gaba_cn_pvals, y_max, by = c("correction", "voxel"))

# MANOVA annotation — placed just above y_min of data
y_min_fig1 <- emm_GABA_CN %>%
  group_by(correction) %>%
  summarise(y_ann = min(emmean - SE) * 0.97, .groups = "drop")

manova_annot_fig1 <- tibble(
  correction = c("Uncorrected (UNC)", "Alpha Tissue Corrected (ATC)"),
  label      = c("MANOVA p < .001, females > males", "MANOVA p < .001, females > males")
) %>% left_join(y_min_fig1, by = "correction")

fig1 <- ggplot(emm_GABA_CN, aes(x = voxel, y = emmean, colour = Sex, group = Sex)) +
  geom_point(size = 3.5, position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(ymin = emmean - SE, ymax = emmean + SE),
                width = 0.15, linewidth = 0.8, position = position_dodge(width = 0.4)) +
  geom_text(data = annot_GABA_CN, aes(x = voxel, y = y_pos, label = sig),
            inherit.aes = FALSE, size = 5, colour = "black") +
  facet_wrap(~ correction) +
  scale_colour_manual(values = sex_colours) +
  labs(
    title    = "GABA: Sex Differences in Older CN Adults",
    subtitle = "Estimated marginal means +/- SE, adjusted for age.\nMANOVA: UNC p < .001, ATC p < .001 (females > males, all voxels).",
    x        = NULL,
    y        = "GABA (estimated marginal mean)"
  ) +
  theme_presentation()

ggsave("plots/presentation/fig1_GABA_older_CN.pdf", fig1, width = 10, height = 5)
cat("Saved fig1_GABA_older_CN.pdf\n")


# =============================================================================
# FIG 2: GLX OLDER CN — ESTIMATED MARGINAL MEANS BY SEX
# =============================================================================

emm_list <- list()
for (corr in corrections) {
  for (vox in voxels) {
    col <- paste0("GLX_", corr, "_", vox)
    emm <- get_emmeans(older_CN %>% select(Sex, Age_DuringParticipation, all_of(col)) %>% drop_na(), col)
    emm$correction <- corr
    emm$voxel      <- vox
    emm_list[[paste0(corr, "_", vox)]] <- emm
  }
}

emm_GLX_CN <- bind_rows(emm_list) %>%
  mutate(
    voxel      = factor(voxel, levels = voxels, labels = voxel_labels),
    correction = factor(correction, levels = c("UNC", "ATC"),
                        labels = c("Uncorrected (UNC)", "Alpha Tissue Corrected (ATC)"))
  )

glx_cn_pvals <- tribble(
  ~correction, ~voxel, ~p,
  "Uncorrected (UNC)", "Auditory\n(AUD)",      0.0001,
  "Uncorrected (UNC)", "Sensorimotor\n(SM)",   0.0018,
  "Uncorrected (UNC)", "Visual\n(VV)",         0.0351,
  "Alpha Tissue Corrected (ATC)", "Auditory\n(AUD)",    0.0057,
  "Alpha Tissue Corrected (ATC)", "Sensorimotor\n(SM)", 0.0883,
  "Alpha Tissue Corrected (ATC)", "Visual\n(VV)",       0.2331
) %>%
  mutate(sig = map_chr(p, sig_label))

y_max <- emm_GLX_CN %>%
  group_by(correction, voxel) %>%
  summarise(y_pos = max(emmean + SE) * 1.12, .groups = "drop")

annot_GLX_CN <- left_join(glx_cn_pvals, y_max, by = c("correction", "voxel"))

# MANOVA annotation
y_min_fig2 <- emm_GLX_CN %>%
  group_by(correction) %>%
  summarise(y_ann = min(emmean - SE) * 0.97, .groups = "drop")

manova_annot_fig2 <- tibble(
  correction = c("Uncorrected (UNC)", "Alpha Tissue Corrected (ATC)"),
  label      = c("MANOVA p < .001, females > males", "MANOVA p = .039, females > males")
) %>% left_join(y_min_fig2, by = "correction")

fig2 <- ggplot(emm_GLX_CN, aes(x = voxel, y = emmean, colour = Sex, group = Sex)) +
  geom_point(size = 3.5, position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(ymin = emmean - SE, ymax = emmean + SE),
                width = 0.15, linewidth = 0.8, position = position_dodge(width = 0.4)) +
  geom_text(data = annot_GLX_CN, aes(x = voxel, y = y_pos, label = sig),
            inherit.aes = FALSE, size = 5, colour = "black") +
  facet_wrap(~ correction) +
  scale_colour_manual(values = sex_colours) +
  labs(
    title    = "Glx: Sex Differences in Older CN Adults",
    subtitle = "Estimated marginal means +/- SE, adjusted for age.\nMANOVA: UNC p < .001, ATC p = .039 (females > males). ATC effect distributed across voxels.",
    x        = NULL,
    y        = "Glx (estimated marginal mean)"
  ) +
  theme_presentation()

ggsave("plots/presentation/fig2_GLX_older_CN.pdf", fig2, width = 10, height = 5)
cat("Saved fig2_GLX_older_CN.pdf\n")


# =============================================================================
# FIG 3: GABA BY SEX x DIAGNOSIS — OLDER ADULTS
# =============================================================================

emm_list <- list()
for (corr in corrections) {
  for (vox in voxels) {
    col  <- paste0("G_", corr, "_", vox)
    dat_sub <- older_all %>%
      select(Sex, Subgroup_f, Age_DuringParticipation, all_of(col)) %>%
      drop_na()
    fit  <- lm(as.formula(paste0("`", col, "` ~ Sex * Subgroup_f + Age_DuringParticipation")),
               data = dat_sub)
    em   <- as.data.frame(emmeans(fit, specs = ~ Sex * Subgroup_f))
    em$correction <- corr
    em$voxel      <- vox
    emm_list[[paste0(corr, "_", vox)]] <- em
  }
}

emm_GABA_dx <- bind_rows(emm_list) %>%
  mutate(
    voxel      = factor(voxel, levels = voxels, labels = voxel_labels),
    correction = factor(correction, levels = c("UNC", "ATC"),
                        labels = c("Uncorrected (UNC)", "Alpha Tissue Corrected (ATC)")),
    Group      = paste0(Sex, "\n", Subgroup_f)
  )

# MANOVA annotation for Sex main effect
y_min_fig3 <- emm_GABA_dx %>%
  group_by(correction) %>%
  summarise(y_ann = min(emmean - SE) * 0.97, .groups = "drop")

manova_annot_fig3 <- tibble(
  correction = c("Uncorrected (UNC)", "Alpha Tissue Corrected (ATC)"),
  label      = c("Sex MANOVA p < .001, females > males", "Sex MANOVA p < .001, females > males")
) %>% left_join(y_min_fig3, by = "correction")

fig3 <- ggplot(emm_GABA_dx, aes(x = voxel, y = emmean, colour = Sex, shape = Subgroup_f, group = interaction(Sex, Subgroup_f))) +
  geom_point(size = 3.5, position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = emmean - SE, ymax = emmean + SE),
                width = 0.15, linewidth = 0.8, position = position_dodge(width = 0.5)) +
  facet_wrap(~ correction) +
  scale_colour_manual(values = sex_colours) +
  scale_shape_manual(values = c("CN" = 16, "MCI" = 17), name = "Diagnosis") +
  labs(
    title    = "GABA: Sex Differences by Diagnosis in Older Adults",
    subtitle = "Estimated marginal means +/- SE, adjusted for age.\nSex MANOVA: UNC p < .001, ATC p < .001 (females > males). No significant Sex x Diagnosis interaction.",
    x        = NULL,
    y        = "GABA (estimated marginal mean)"
  ) +
  theme_presentation() +
  theme(legend.position = "top")

ggsave("plots/presentation/fig3_GABA_diagnosis.pdf", fig3, width = 10, height = 5)
cat("Saved fig3_GABA_diagnosis.pdf\n")


# =============================================================================
# FIG 4: LONGITUDINAL TRAJECTORIES — GABA OLDER CN
# =============================================================================

# Compute Time_Years
dates_clean <- dates %>%
  mutate(SubNum  = as.character(SubNum),
         MRSDate = as.Date(MRSDate, format = "%m/%d/%y")) %>%
  filter(!is.na(MRSDate)) %>%
  select(SubNum, Wave_Num, MRSDate)

wave1_dates <- dates_clean %>%
  filter(Wave_Num == 1) %>%
  select(SubNum, Wave1_Date = MRSDate)

dates_with_time <- dates_clean %>%
  left_join(wave1_dates, by = "SubNum") %>%
  mutate(Time_Years = as.numeric(difftime(MRSDate, Wave1_Date, units = "days")) / 365.25) %>%
  select(SubNum, Wave_Num, Time_Years)

older_CN_long <- dat %>%
  filter(AgeCategory == "Older", Subgroup_f == "CN", Wave_Num %in% c(1, 2)) %>%
  mutate(SubNum = as.character(SubNum)) %>%
  left_join(dates_with_time, by = c("SubNum", "Wave_Num")) %>%
  mutate(Time_Years = ifelse(Wave_Num == 1, 0, Time_Years))

# Fit LME and get predicted trajectories for each voxel x correction
traj_list <- list()
for (corr in c("UNC", "ATC")) {
  for (vox in voxels) {
    col     <- paste0("G_", corr, "_", vox)
    fit_dat <- older_CN_long %>%
      select(SubNum, Sex, Age_DuringParticipation, Time_Years, all_of(col)) %>%
      drop_na()
    
    fit <- lmer(as.formula(paste0("`", col, "` ~ Sex * Time_Years + Age_DuringParticipation + (1|SubNum)")),
                data = fit_dat, REML = FALSE)
    
    # Predicted values at W1 (0) and mean W2 (~4 years)
    new_dat <- expand.grid(
      Sex                    = c("Female", "Male"),
      Time_Years             = seq(0, 5, by = 0.1),
      Age_DuringParticipation = mean(fit_dat$Age_DuringParticipation, na.rm = TRUE)
    )
    new_dat$predicted <- predict(fit, newdata = new_dat, re.form = NA)
    new_dat$correction <- corr
    new_dat$voxel      <- vox
    traj_list[[paste0(corr, "_", vox)]] <- new_dat
  }
}

traj_GABA <- bind_rows(traj_list) %>%
  mutate(
    voxel      = factor(voxel, levels = voxels, labels = voxel_labels),
    correction = factor(correction, levels = c("UNC", "ATC"),
                        labels = c("Uncorrected (UNC)", "Alpha Tissue Corrected (ATC)"))
  )

fig4 <- ggplot(traj_GABA, aes(x = Time_Years, y = predicted, colour = Sex)) +
  geom_line(linewidth = 1.1) +
  facet_grid(voxel ~ correction, scales = "free_y") +
  scale_colour_manual(values = sex_colours) +
  labs(
    title    = "GABA: Longitudinal Trajectories by Sex — Older CN Adults",
    subtitle = "Model-predicted values (adjusted for age). Sex difference stable over time.",
    x        = "Time since Wave 1 (years)",
    y        = "GABA (model-predicted)"
  ) +
  theme_presentation()

ggsave("plots/presentation/fig4_longitudinal.pdf", fig4, width = 10, height = 8)
cat("Saved fig4_longitudinal.pdf\n")


# =============================================================================
# FIG 5: GABA YOUNG CN — FOR COMPARISON WITH OLDER CN
# =============================================================================

emm_list <- list()
for (corr in corrections) {
  for (vox in voxels) {
    col <- paste0("G_", corr, "_", vox)
    emm <- get_emmeans(young_CN %>% select(Sex, Age_DuringParticipation, all_of(col)) %>% drop_na(), col)
    emm$correction <- corr
    emm$voxel      <- vox
    emm_list[[paste0(corr, "_", vox)]] <- emm
  }
}

emm_GABA_young <- bind_rows(emm_list) %>%
  mutate(
    voxel      = factor(voxel, levels = voxels, labels = voxel_labels),
    correction = factor(correction, levels = c("UNC", "ATC"),
                        labels = c("Uncorrected (UNC)", "Alpha Tissue Corrected (ATC)"))
  )

gaba_young_pvals <- tribble(
  ~correction, ~voxel, ~p,
  "Uncorrected (UNC)", "Auditory\n(AUD)",      0.2464,
  "Uncorrected (UNC)", "Sensorimotor\n(SM)",   0.0904,
  "Uncorrected (UNC)", "Visual\n(VV)",         0.1103,
  "Alpha Tissue Corrected (ATC)", "Auditory\n(AUD)",    0.3381,
  "Alpha Tissue Corrected (ATC)", "Sensorimotor\n(SM)", 0.1291,
  "Alpha Tissue Corrected (ATC)", "Visual\n(VV)",       0.2634
) %>%
  mutate(sig = map_chr(p, sig_label))

y_max <- emm_GABA_young %>%
  group_by(correction, voxel) %>%
  summarise(y_pos = max(emmean + SE) * 1.12, .groups = "drop")

annot_GABA_young <- left_join(gaba_young_pvals, y_max, by = c("correction", "voxel"))

# Compute mean beta across voxels per correction to determine overall direction
# Betas from primary script: G UNC: AUD=0.0546, SM=-0.0878, VV=-0.0687
#                            G ATC: AUD=0.0700, SM=-0.1254, VV=-0.0745
mean_beta_unc <- mean(c(0.0546, -0.0878, -0.0687))  # -0.034 -> females higher overall
mean_beta_atc <- mean(c(0.0700, -0.1254, -0.0745))  # -0.043 -> females higher overall

y_min_fig5 <- emm_GABA_young %>%
  group_by(correction) %>%
  summarise(y_ann = min(emmean - SE) * 0.97, .groups = "drop")

manova_annot <- tibble(
  correction = c("Uncorrected (UNC)", "Alpha Tissue Corrected (ATC)"),
  label      = c(
    paste0("MANOVA p = .015, females > males\n(mean beta = ", round(mean_beta_unc, 3), "; AUD diverges)"),
    paste0("MANOVA p = .040, females > males\n(mean beta = ", round(mean_beta_atc, 3), "; AUD diverges)")
  )
) %>% left_join(y_min_fig5, by = "correction")

fig5 <- ggplot(emm_GABA_young, aes(x = voxel, y = emmean, colour = Sex, group = Sex)) +
  geom_point(size = 3.5, position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(ymin = emmean - SE, ymax = emmean + SE),
                width = 0.15, linewidth = 0.8, position = position_dodge(width = 0.4)) +
  geom_text(data = annot_GABA_young, aes(x = voxel, y = y_pos, label = sig),
            inherit.aes = FALSE, size = 5, colour = "black") +
  facet_wrap(~ correction) +
  scale_colour_manual(values = sex_colours) +
  labs(
    title    = "GABA: Sex Differences in Young CN Adults",
    subtitle = paste0("Estimated marginal means +/- SE, adjusted for age.\n",
                      "MANOVA: UNC p = .015, ATC p = .040 (females > males overall; mean beta = ",
                      round(mean_beta_unc, 3), " UNC, ", round(mean_beta_atc, 3), " ATC).\n",
                      "No individual voxel survives — AUD trends opposite to SM and VV."),
    x        = NULL,
    y        = "GABA (estimated marginal mean)"
  ) +
  theme_presentation()

ggsave("plots/presentation/fig5_GABA_young.pdf", fig5, width = 10, height = 5)
cat("Saved fig5_GABA_young.pdf\n")

cat("\nAll figures saved to plots/presentation/\n")