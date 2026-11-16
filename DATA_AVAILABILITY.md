# Data availability

**No participant-level data are included in this repository.** Both source
cohorts are restricted-access, and their data-use agreements prohibit
redistribution. To reproduce the analyses you must obtain the data yourself
(both are free for approved researchers) and place the files under a local,
git-ignored `data/` directory so that the relative paths in the notebooks
resolve.

## Health and Retirement Study (HRS, United States)

- Provider: HRS, University of Michigan — https://hrsdata.isr.umich.edu
- Register for a free account and accept the conditions of use.
- Files used:
  - **RAND HRS Longitudinal File** (`randhrs1992_2022v1.dta`)
  - **RAND HRS Fat files**, one per wave (`h04f1c`, `h06f4b`, `h08f3b`,
    `hd10f6b`, `h12f3a`, `h14f2b`, `h16f2c`, `h18f2c`, `h20f1b`)
  - **Langa–Weir Classification of Cognitive Function**
    (`cogfinalimp_9522wide.dta`)
  - **Cross-Wave Tracker** (`trk2022tr_r.csv`)
- Expected location: `data/HRS/…` (see paths at the top of
  `code/HRS_dataextract.Rmd`).

## English Longitudinal Study of Ageing (ELSA)

- Provider: UK Data Service — https://ukdataservice.ac.uk (study **SN 5050**)
- Register, apply for the End User Licence, and download the Stata edition.
- Files used: core and derived data files for waves 1–9 (`wave_1_core_data_v3`,
  `wave_2_core_data_v4`, `wave_2_derived_variables`, … through wave 9).
- Expected location: `data/ELSA/UKDA-5050-stata/…` (see the path at the top of
  `code/ELSA_extractdata.Rmd`).

## Derived, harmonised datasets

The two extraction notebooks write intermediate files to `data/` that the
analysis notebooks read back in:

| File | Written by | Contents |
|------|------------|----------|
| `hrs_long.csv` | `HRS_dataextract.Rmd` | Harmonised HRS person-wave data, frailty index and dementia status |
| `death_track.rds` | `HRS_dataextract.Rmd` | One row per HRS participant: year and month last known alive, and of death |
| `hrs_age_months.csv` | `HRS_dataextract.Rmd` | Age in months at each HRS interview (attained-age analysis) |
| `hrs_fi_reduced.csv`, `hrs_fi_item_counts.csv` | `HRS_dataextract.Rmd` | Reduced frailty index without (I)ADL and affective items, and its item counts |
| `elsa_long.csv` | `ELSA_extractdata.Rmd` | Harmonised ELSA person-wave data, frailty index and dementia status |
| `elsa_fi_reduced.csv`, `elsa_fi_item_counts.csv` | `ELSA_extractdata.Rmd` | Reduced frailty index for ELSA, and its item counts |

These are derived from restricted data and are therefore **also git-ignored**;
regenerate them locally by running the extraction notebooks first. The files of
the same names under `demo/` are simulated and contain no participant data.
