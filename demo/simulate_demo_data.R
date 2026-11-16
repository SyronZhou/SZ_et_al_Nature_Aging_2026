# ---------------------------------------------------------------------------
# simulate_demo_data.R
#
# Generates a small SIMULATED dataset with exactly the schema that the
# analysis notebooks (code/HRS_analysis_MSM.Rmd, code/HRS_analysis_Cox.Rmd,
# code/ELSA_analysis_MSM.Rmd, code/ELSA_analysis_Cox.Rmd) expect, so that they
# can be run end to end without access to the restricted HRS and ELSA source
# data.
#
# The values are entirely synthetic. They are generated to be structurally
# realistic (biennial waves, monotone dementia, frailty that accrues with age,
# sensory impairment that can onset and remit, deaths during follow-up)
# so that every model in the notebooks is exercised and returns a finite
# estimate. They are NOT calibrated to reproduce the published effect sizes,
# and no value derives from any real participant.
#
# Usage:  Rscript demo/simulate_demo_data.R
# Writes: data/hrs_long.csv, data/death_track.rds, data/hrs_age_months.csv,
#         data/hrs_fi_reduced.csv, data/hrs_fi_item_counts.csv,
#         data/elsa_long.csv, data/elsa_fi_reduced.csv,
#         data/elsa_fi_item_counts.csv
# ---------------------------------------------------------------------------

suppressPackageStartupMessages(library(dplyr))

set.seed(2026)

N_HRS   <- 600L
N_ELSA  <- 600L
HRS_WAVES  <- seq(2004L, 2020L, by = 2L)
ELSA_WAVES <- seq(2002L, 2018L, by = 2L)

# Paths are relative to the REPOSITORY ROOT, matching the analysis notebooks,
# so the script writes where the notebooks look no matter where it is run from.
repo_root <- if (basename(getwd()) %in% c("code", "demo")) normalizePath("..") else normalizePath(".")
data_dir  <- file.path(repo_root, "data")
dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)
cat("Repository root:", repo_root, "\n")

# ---------------------------------------------------------------------------
# Helper: build one cohort's person-wave panel
# ---------------------------------------------------------------------------
simulate_cohort <- function(n, waves, id_name,
                            fi_col, frail_col, fi_noiadl_col, frail_noiadl_col,
                            baseline_fi_mean) {

  people <- tibble(
    .id               = seq_len(n),
    birth_year        = sample(1915:1958, n, replace = TRUE),
    Gender            = sample(c(1L, 2L), n, replace = TRUE, prob = c(0.44, 0.56)),
    Race              = sample(c(0L, 1L), n, replace = TRUE, prob = c(0.75, 0.25)),
    Education         = sample(0:2, n, replace = TRUE, prob = c(0.25, 0.35, 0.40)),
    Marital_status    = sample(c(0L, 1L), n, replace = TRUE, prob = c(0.32, 0.68)),
    # person-level latent frailty level and rate of accrual
    fi_intercept      = pmax(0, rnorm(n, baseline_fi_mean, 0.10)),
    fi_slope          = pmax(0, rnorm(n, 0.012, 0.006)),
    hear_onset_wave   = sample(c(NA_integer_, seq_along(waves)), n, replace = TRUE,
                               prob = c(0.55, rep(0.45 / length(waves), length(waves)))),
    vis_onset_wave    = sample(c(NA_integer_, seq_along(waves)), n, replace = TRUE,
                               prob = c(0.60, rep(0.40 / length(waves), length(waves))))
  )

  panel <- tidyr::expand_grid(.id = people$.id, wave_idx = seq_along(waves)) %>%
    left_join(people, by = ".id") %>%
    mutate(
      year = waves[wave_idx],
      age  = year - birth_year,
      # frailty index accrues with time, bounded to (0, 1)
      fi_raw = fi_intercept + fi_slope * (wave_idx - 1) + rnorm(n(), 0, 0.03),
      fi     = pmin(pmax(fi_raw, 0), 0.95),
      frail  = as.integer(fi >= 0.25)
    )

  # Force a minority of people to acquire sensory impairment in the SAME wave as
  # frailty onset. Without this, same-interval onset (multi-state state 6) arises
  # only by chance, leaving that state too sparse for msm::qratio.msm() to form
  # a rate ratio out of it -- which fails with "object 'L' not found" whenever a
  # transition rate is exactly zero.
  first_frail <- panel %>%
    filter(frail == 1L) %>%
    group_by(.id) %>%
    summarise(first_frail_wave = min(wave_idx), .groups = "drop")

  people <- people %>%
    left_join(first_frail, by = ".id") %>%
    mutate(
      sim_hear = !is.na(first_frail_wave) & runif(n()) < 0.12,
      sim_vis  = !is.na(first_frail_wave) & runif(n()) < 0.12,
      hear_onset_wave = if_else(sim_hear, first_frail_wave, hear_onset_wave),
      vis_onset_wave  = if_else(sim_vis,  first_frail_wave, vis_onset_wave)
    ) %>%
    dplyr::select(-sim_hear, -sim_vis, -first_frail_wave)

  panel <- panel %>%
    dplyr::select(-hear_onset_wave, -vis_onset_wave) %>%
    left_join(people %>% dplyr::select(.id, hear_onset_wave, vis_onset_wave),
              by = ".id") %>%
    mutate(
      # sensory impairment: present from onset, with occasional remission so the
      # time-varying exposures in the Cox models are genuinely reversible
      hearing_loss = if_else(!is.na(hear_onset_wave) & wave_idx >= hear_onset_wave,
                             rbinom(n(), 1L, 0.85), 0L),
      vision_overall = if_else(!is.na(vis_onset_wave) & wave_idx >= vis_onset_wave,
                               rbinom(n(), 1L, 0.85), 0L),
      # time-varying lifestyle covariates
      Employment        = rbinom(n(), 1L, 0.45),
      Smoking           = rbinom(n(), 1L, 0.17),
      Alcohol           = rbinom(n(), 1L, 0.40),
      Physical_activity = rbinom(n(), 1L, 0.55),
      Social_activity   = rbinom(n(), 1L, 0.60)
    ) %>%
    group_by(.id) %>%
    arrange(wave_idx, .by_group = TRUE) %>%
    mutate(
      # a slightly leaner index for the (I)ADL-removed sensitivity analysis
      fi_noiadl    = pmin(pmax(fi - runif(n(), 0.00, 0.04), 0), 0.95),
      frail_noiadl = as.integer(fi_noiadl >= 0.25),
      # dementia hazard rises with age, frailty and sensory impairment
      lp = -5.6 + 0.075 * (age - 65) + 0.75 * frail +
             0.30 * hearing_loss + 0.35 * vision_overall,
      dem_draw = rbinom(n(), 1L, plogis(lp)),
      dementia = cummax(dem_draw)                      # absorbing, as in the real data
    ) %>%
    ungroup()

  # introduce a realistic sprinkle of item-level missingness
  blank <- function(x, p) { x[runif(length(x)) < p] <- NA; x }
  panel <- panel %>%
    mutate(
      Race              = blank(Race, 0.004),
      Education         = blank(Education, 0.010),
      Marital_status    = blank(Marital_status, 0.008),
      Smoking           = blank(Smoking, 0.012),
      Alcohol           = blank(Alcohol, 0.012),
      Employment        = blank(Employment, 0.010),
      hearing_loss      = blank(hearing_loss, 0.010),
      vision_overall    = blank(vision_overall, 0.010),
      Physical_activity = blank(Physical_activity, 0.015),
      Social_activity   = blank(Social_activity, 0.015)
    )

  out <- panel %>%
    transmute(
      !!id_name := sprintf("%06d", .id),
      source = paste0(tolower(substr(id_name, 1, 1)), year),
      year, age, birth_year,
      Gender, Race, Education, Marital_status,
      Employment, Smoking, Alcohol, Physical_activity, Social_activity,
      hearing_loss, vision_overall,
      !!fi_col           := round(fi, 4),
      !!frail_col        := frail,
      !!fi_noiadl_col    := round(fi_noiadl, 4),
      !!frail_noiadl_col := frail_noiadl,
      dementia
    ) %>%
    filter(age >= 45)   # keep the age filter in the notebooks meaningful

  out
}

# ---------------------------------------------------------------------------
# HRS: uppercase FI column names, plus a death tracker
# ---------------------------------------------------------------------------
hrs_long <- simulate_cohort(
  n = N_HRS, waves = HRS_WAVES, id_name = "hhidpn",
  fi_col = "FI_score", frail_col = "Frailty_status",
  fi_noiadl_col = "FI_score_noIADL", frail_noiadl_col = "Frailty_status_noIADL",
  baseline_fi_mean = 0.18
)

# Death tracker: ~28% die. A decedent is last interviewed at a wave between
# 2008 and 2020 and dies 1-2 years later, so no interview follows a death.
# Deaths after the 2020 wave are censored by the analysis notebooks, as in the
# real data.
death_draw <- hrs_long %>%
  distinct(hhidpn) %>%
  mutate(
    died      = rbinom(n(), 1L, 0.28),
    last_wave = if_else(died == 1L,
                        sample(seq(2008L, 2020L, by = 2L), n(), replace = TRUE),
                        NA_integer_)
  )

hrs_long <- hrs_long %>%
  left_join(death_draw, by = "hhidpn") %>%
  filter(is.na(last_wave) | year <= last_wave) %>%
  dplyr::select(-died, -last_wave)

readr::write_csv(hrs_long, file.path(data_dir, "hrs_long.csv"))

last_year <- hrs_long %>% group_by(hhidpn) %>% summarise(last = max(year), .groups = "drop")

death_track <- last_year %>%
  left_join(death_draw %>% dplyr::select(hhidpn, died), by = "hhidpn") %>%
  mutate(
    death_year  = if_else(died == 1L, last + sample(1:2, n(), replace = TRUE), NA_real_),
    death_month = if_else(died == 1L, as.numeric(sample(1:12, n(), replace = TRUE)), NA_real_),
    # exercise the "missing month -> mid-year" branch in the analysis notebook
    death_month = if_else(died == 1L & runif(n()) < 0.05, NA_real_, death_month),
    lastalive_year  = as.numeric(last),
    lastalive_month = 6
  ) %>%
  transmute(hhidpn, death_year, death_month,
            lastalive_year, lastalive_month, died)

saveRDS(death_track, file.path(data_dir, "death_track.rds"))

# Age in months at each interview (attained-age time scale): interview month
# June, birth month drawn at random.
birth_month <- hrs_long %>% distinct(hhidpn) %>%
  mutate(birth_month = sample(1:12, n(), replace = TRUE))
hrs_age_months <- hrs_long %>%
  left_join(birth_month, by = "hhidpn") %>%
  transmute(hhidpn, year,
            age_m = (year * 12 + 6) - (birth_year * 12 + birth_month),
            src   = "interview date - birth date")
readr::write_csv(hrs_age_months, file.path(data_dir, "hrs_age_months.csv"))

# Reduced frailty index ((I)ADL and affective CES-D items removed): a slightly
# leaner score than the main index, scored on 22 items.
simulate_reduced_fi <- function(long, id, fi_col, prefix, n_items, n_full, n_noiadl) {
  fi <- long[[fi_col]]
  red <- tibble(
    !!id := long[[id]],
    year = long$year,
    score = round(pmin(pmax(fi - runif(length(fi), 0.00, 0.05), 0), 0.95), 4)
  )
  red$frail <- as.integer(red$score >= 0.25)
  red$n_ans <- ifelse(is.na(red$score), NA_integer_, n_items)
  names(red)[3:5] <- paste0(prefix, c("score_reduced", "_reduced", "n_answered_reduced"))
  counts <- tibble(year = sort(unique(long$year)), items_full = n_full,
                   items_noIADL = n_noiadl, items_reduced = n_items,
                   cutoff_reduced = ceiling(0.80 * n_items))
  list(scores = red, counts = counts)
}
hrs_red <- simulate_reduced_fi(hrs_long, "hhidpn", "FI_score", "x", 22L, 37L, 26L)
names(hrs_red$scores)[3:5] <- c("FI_score_reduced", "Frailty_status_reduced", "fi_n_answered_reduced")
readr::write_csv(hrs_red$scores, file.path(data_dir, "hrs_fi_reduced.csv"))
readr::write_csv(hrs_red$counts, file.path(data_dir, "hrs_fi_item_counts.csv"))

# ---------------------------------------------------------------------------
# ELSA: lowercase FI column names, no mortality file
# ---------------------------------------------------------------------------
elsa_long <- simulate_cohort(
  n = N_ELSA, waves = ELSA_WAVES, id_name = "idauniq",
  fi_col = "fi_score", frail_col = "frailty",
  fi_noiadl_col = "fi_score_noIADL", frail_noiadl_col = "frailty_noIADL",
  baseline_fi_mean = 0.12
)

readr::write_csv(elsa_long, file.path(data_dir, "elsa_long.csv"))

elsa_red <- simulate_reduced_fi(elsa_long, "idauniq", "fi_score", "x", 31L, 49L, 36L)
names(elsa_red$scores)[3:5] <- c("fi_score_reduced", "frailty_reduced", "fi_n_answered_reduced")
readr::write_csv(elsa_red$scores, file.path(data_dir, "elsa_fi_reduced.csv"))
readr::write_csv(elsa_red$counts, file.path(data_dir, "elsa_fi_item_counts.csv"))

# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------
cat("\nSimulated demo data written to data/\n")
cat(sprintf("  data/hrs_long.csv     %6d rows, %3d people, %2d waves (%d-%d)\n",
            nrow(hrs_long), dplyr::n_distinct(hrs_long$hhidpn),
            length(HRS_WAVES), min(HRS_WAVES), max(HRS_WAVES)))
cat(sprintf("  data/death_track.rds  %6d people, %3d deaths (%d by the 2020 wave)\n",
            nrow(death_track), sum(death_track$died),
            sum(death_track$died == 1 & death_track$death_year <= 2020)))
cat("  data/hrs_age_months.csv, data/hrs_fi_reduced.csv, data/hrs_fi_item_counts.csv\n")
cat(sprintf("  data/elsa_long.csv    %6d rows, %3d people, %2d waves (%d-%d)\n",
            nrow(elsa_long), dplyr::n_distinct(elsa_long$idauniq),
            length(ELSA_WAVES), min(ELSA_WAVES), max(ELSA_WAVES)))
cat("  data/elsa_fi_reduced.csv, data/elsa_fi_item_counts.csv\n")
cat(sprintf("\n  HRS  person-waves with dementia: %d (%.1f%% of people ever affected)\n",
            sum(hrs_long$dementia == 1, na.rm = TRUE),
            100 * mean(tapply(hrs_long$dementia, hrs_long$hhidpn, max, na.rm = TRUE) == 1)))
cat(sprintf("  ELSA person-waves with dementia: %d (%.1f%% of people ever affected)\n",
            sum(elsa_long$dementia == 1, na.rm = TRUE),
            100 * mean(tapply(elsa_long$dementia, elsa_long$idauniq, max, na.rm = TRUE) == 1)))
cat("\nThese files are synthetic. Do not interpret any estimate derived from them.\n")
