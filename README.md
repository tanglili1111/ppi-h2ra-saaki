# PPI vs. H2RA in Sepsis-Associated Acute Kidney Injury (SA-AKI): A Target-Trial Emulation with External Validation

Reproducible analysis code for a target-trial emulation comparing **proton-pump
inhibitors (PPIs)** with **histamine-2 receptor antagonists (H2RAs)** as the first
in-ICU acid suppressant in patients with **sepsis-associated acute kidney injury
(SA-AKI)**. The primary analysis uses **MIMIC-IV**; the findings are externally
validated in the **eICU Collaborative Research Database**, with an additional
cohort built from the **Salzburg Intensive Care database (SICdb)**.

> **No patient-level data are contained in this repository.** MIMIC-IV, eICU, and
> SICdb are credentialed-access resources. Only analysis code and aggregate,
> shareable outputs (effect tables, balance diagnostics, figures) are distributed
> here.

---

## Table of contents

- [Study design](#study-design)
- [What this repository contains](#what-this-repository-contains)
- [Analysis pipeline](#analysis-pipeline)
  - [MIMIC-IV (primary)](#mimic-iv-primary)
  - [eICU (external validation)](#eicu-external-validation)
  - [SQL cohort-extraction pipelines](#sql-cohort-extraction-pipelines)
- [Requirements](#requirements)
- [Data access and the upstream SQL pipelines](#data-access-and-the-upstream-sql-pipelines)
- [Configuration: database credentials](#configuration-database-credentials)
- [How to run](#how-to-run)
- [Outputs](#outputs)
- [Reproducibility](#reproducibility)
- [Notes and limitations](#notes-and-limitations)
- [Disclaimer](#disclaimer)
- [Citation](#citation)
- [License](#license)

---

## Study design

| Element | Specification |
| --- | --- |
| **Framework** | Target-trial emulation; active-comparator, new-user design |
| **Exposure** | First in-ICU acid suppressant (PPI vs. H2RA) within the window `[T0, T0 + 24h]` (main) or `[T0, T0 + 48h]` (MIMIC-IV sensitivity) |
| **Time zero (`T0`)** | Onset of sepsis-associated AKI (KDIGO by creatinine and urine output) |
| **Primary estimand** | Per-protocol effect via **clone-censor-weight (CCW)**, grace period of 1 day (main) / 2 days (sensitivity) |
| **Co-primary / sensitivity estimand** | Head-to-head **IPTW** (inverse-probability-of-treatment weighting), an ITT-analog |
| **Primary outcomes** | **MAKE30** and **MAKE90** (Major Adverse Kidney Events) in MIMIC-IV; in-hospital **MAKE** and death in eICU |
| **Secondary outcomes** | All-cause mortality; renal replacement therapy (RRT) initiation; renal recovery; AKI progression to KDIGO stage 3 (competing-risk endpoints) |
| **Confounders** | 56 pre-specified baseline covariates in MIMIC-IV; the 24-covariate subset available in eICU, all measured at or before `T0` |
| **Missing data** | Multiple imputation by chained equations (MICE), `m = 20` |
| **Inference** | Estimates pooled across imputations by **Rubin's rules**; CCW confidence intervals via a cluster bootstrap that re-estimates the weights inside every resample; Benjamini-Hochberg FDR control across pre-specified secondary endpoints |
| **Robustness** | Negative-control outcomes, E-values, a prevalent-user-bias gradient (prevalent-only vs. all vs. new-user), and external replication in eICU and SICdb |

---

## What this repository contains

```
ppi-h2ra-saaki/
├── README.md                 # this file
├── LICENSE                   # MIT (review/replace as appropriate)
├── .gitignore                # excludes caches, generated outputs, credentials
├── .gitattributes            # enforces LF line endings
├── CITATION.cff              # machine-readable citation metadata ("Cite this repository")
├── mimic/                    # MIMIC-IV (primary cohort)
│   ├── ppi_h2ra_saaki_pipeline.sql   # cohort-extraction pipeline (builds analytic_cohort and _48h)
│   ├── 01_MIMIC_24h.R        # primary analysis (24h exposure window)
│   └── 01_MIMIC_48h.R        # sensitivity analysis (48h exposure window)
├── eicu/                     # eICU (external validation)
│   ├── eICU_pipeline.sql     # cohort-extraction pipeline (Stages 0-4)
│   └── eICU.R                # external-validation analysis
└── sicdb/                    # SICdb (external validation)
    └── SICdb_pipeline.sql    # SA-AKI cohort construction + exposure grouping
```

Not included (must be supplied by the user — see below):

- The databases themselves: MIMIC-IV, eICU, and SICdb (credentialed access via
  PhysioNet).
- The `public.elixhauser_scores` build script for MIMIC-IV (see "Prerequisites
  for the MIMIC-IV SQL" below); add it under `mimic/` for a fully self-contained
  pipeline.

---

## Analysis pipeline

### MIMIC-IV (primary)

`mimic/01_MIMIC_24h.R` is a single top-to-bottom pipeline organised into modules
**R0–R16**: environment and helpers; load and range-check; MICE (`m = 20`); IPTW
weighting and balance; Table 1; the primary IPTW effects for MAKE30/MAKE90; the
competing-risk secondary outcomes; negative controls and E-values; forest and
publication figures; pre-specified subgroups; a new-user sensitivity analysis;
the full clone-censor-weight (CCW) per-protocol analysis with a weight-re-estimating
cluster bootstrap; supplementary analyses that characterise prevalent-user bias;
and Benjamini-Hochberg FDR over the secondary endpoints.

`mimic/01_MIMIC_48h.R` repeats the corrected main analysis with a 48-hour exposure
window / CCW grace period, confirming that both the primary null and the
prevalent-user-bias signature are robust to the window width.

**Prerequisites for the MIMIC-IV SQL.** `mimic/ppi_h2ra_saaki_pipeline.sql` expects
the following to already exist in the database: the MIMIC-IV core schemas
(`mimiciv_hosp`, `mimiciv_icu`); the official mimic-code derived tables (schema
`mimiciv_derived`) — `sepsis3`, `age`, `icustay_detail`, `kdigo_stages`,
`charlson`, `rrt`, `sofa`, `ventilation`, `vasoactive_agent`, `first_day_weight`,
`urine_output`; and a `public.elixhauser_scores(hadm_id, elixhauser_score)` table
holding the van Walraven-weighted Elixhauser score per admission. The script then
builds `analytic_cohort` and `analytic_cohort_48h`, which the two R scripts read.

### eICU (external validation)

`eicu/eICU.R` reproduces the analysis in the eICU Collaborative Research Database
across modules **R0–R11**: environment, the 24 eICU propensity-score covariates,
and shared helpers; load/clean/EDA; MICE (auto-cached); IPTW weighting and
balance; Table 1; the primary in-hospital MAKE and death effects (doubly-robust
adjusting `gi_bleed`, plus a weight-only sensitivity); competing-risk secondary
outcomes with a hospital-stratified site correction for RRT; an RRT sensitivity
analysis; E-values; figures (propensity-score overlap, cumulative-incidence
curves, and an eICU forest plot with the MIMIC-IV estimates as reference); a
master results table; and a new-user (prevalent-user-bias) probe.

Two eICU-specific points: (1) eICU follow-up ends at hospital discharge (about
98% of patients are discharged within 30 days, so the in-hospital horizon
approximates 30 days), so the eICU script uses binary risk differences / risk
ratios and competing-risk models rather than the 30/90-day landmark MAKE of the
MIMIC-IV analysis; (2) RRT initiation shows a spurious decrease in eICU driven by
multi-center between-hospital confounding, which a hospital-stratified Cox model
returns to the null.

### SQL cohort-extraction pipelines

The per-cohort folders contain the database-side pipelines that build the analytic
cohorts from the raw critical-care tables. Each runs top to bottom against a
PostgreSQL instance holding the corresponding database.

**`mimic/ppi_h2ra_saaki_pipeline.sql`** builds the primary-analysis tables
`analytic_cohort` and `analytic_cohort_48h` from the MIMIC-IV core and
`mimiciv_derived` tables (see Prerequisites above). It defines the SA-AKI cohort
(time zero = first KDIGO >= 1 within 7 days of sepsis), applies the exclusions,
assigns PPI-vs-H2RA exposure and the new-initiator flag, and extracts the MAKE
components and the competing-risk secondary outcomes; each module ends with an
aggregate validation count.

**`eicu/eICU_pipeline.sql`** builds `eicu_analytic_cohort`:

| Stage | Purpose |
| --- | --- |
| **Stage 0** | (Optional, recommended) index the large tables (`lab`, `intakeoutput`, `treatment`, `diagnosis`, ...) to speed up joins |
| **Stage 1** | Cohort: sepsis -> AKI (KDIGO by creatinine and urine output) -> exclusions (ESRD, chronic-dialysis dependence, kidney transplant, pregnancy) -> baseline renal function (eGFR by CKD-EPI 2021) |
| **Stage 2** | Exposure: first acid suppressant (PPI vs. H2RA) within `[T0, T0+24h]`, plus new-user determination using eICU's home-medication history (`admissiondrug`) |
| **Stage 3** | In-hospital outcomes mirroring the MIMIC-IV definitions (MAKE, death, RRT, progression, recovery), with follow-up `= [T0, discharge]` |
| **Stage 4** | Propensity-score covariates (APACHE scores, comorbidities, day-1 physiology, ventilation, vasoactives, GI bleeding, lactate) and assembly of the final `eicu_analytic_cohort` |

**`sicdb/SICdb_pipeline.sql`** builds the SA-AKI cohort and exposure groups in the
Salzburg Intensive Care database (SICdb). It comprises two parts:

- **Part 1 — SA-AKI cohort construction** (output: `sic_cohort`, with an
  attrition table). Steps: base layer (first ICU, adult, ICU stay >= 24h, weight
  cleaning) -> sepsis by a sustained-antibiotics Sepsis-3 proxy (first dose <= 48h
  and administration span >= 2 days; an `AdmissionFormHasSepsis` flag is retained
  for sensitivity) -> baseline creatinine (three-tier hierarchy: pre-admission,
  ICU early 24h, then MDRD back-calculation) -> T0 derivation (earliest qualifying
  offset from creatinine or urine output) -> the official **KDIGO_AKI_168**
  published algorithm as the AKI gold standard (rather than unreliable ICD codes)
  plus exclusion flags -> merge.
- **Part 2 — exposure grouping** (output: `sic_exposure`). The first acid
  suppressant within `[T0, T0+24h]` defines the group (PPI = esomeprazole /
  pantoprazole, H2RA = ranitidine, otherwise none/both), with a new-user flag.

The cohort criteria deliberately follow transportable definitions from Luo et al.
(CKJ 2026). One practical limitation is visible in the by-year exposure query:
H2RA use is very rare in SICdb (ranitidine having been withdrawn from the market),
which constrains a head-to-head PPI-vs-H2RA contrast in this database. The
provided SICdb SQL covers cohort construction and exposure assignment only; an
outcomes/covariate stage and a SICdb R analysis are not part of this file.

Inline comments in all three SQL files include recorded cohort counts used as
reference checks during development; these are aggregate and contain no
patient-level data.

---

## Requirements

- **R** ≥ 4.3
- **PostgreSQL** (a running instance holding the MIMIC-IV, eICU, and/or SICdb tables)
- R packages (all on CRAN):

  ```r
  install.packages(c(
    "DBI", "RPostgres",                 # database access
    "dplyr", "tidyr", "tibble",         # data wrangling
    "mice",                             # multiple imputation
    "WeightIt", "cobalt",               # propensity weighting + balance
    "tableone", "survey", "survival",   # Table 1, weighted models, survival / Fine-Gray
    "naniar", "ggplot2", "scales",      # missingness viz + plotting
    "EValue"                            # sensitivity (E-values)
  ))
  ```

  `parallel` and `stats` are part of base R and require no installation.

For the SQL pipelines, the standard database tables must be present:

- **eICU:** `patient, lab, diagnosis, treatment, intakeoutput, medication,
  infusiondrug, admissiondrug, pasthistory, apachepredvar, apacheapsvar,
  apachepatientresult`.
- **SICdb:** `cases, medication, laboratory, data_float_h, data_ref,
  d_references` (with the database-specific constants documented in the file
  header, e.g. creatinine `laboratoryid IN (367,368)`, urine output `dataid=725`).

---

## Data access and the upstream SQL pipelines

This code consumes PostgreSQL tables that are **produced upstream** and are
**not** part of this repository:

- **MIMIC-IV:** `analytic_cohort` (24h window) and `analytic_cohort_48h` (48h),
  built by the companion MIMIC-IV SQL pipeline.
- **eICU:** `eicu_analytic_cohort`, built by `eicu/eICU_pipeline.sql`.
- **SICdb:** `sic_cohort` and `sic_exposure`, built by `sicdb/SICdb_pipeline.sql`.

MIMIC-IV, the eICU Collaborative Research Database, and SICdb are all
credentialed-access resources. To obtain them you must:

1. Complete the required human-subjects research training (CITI),
2. Sign the relevant PhysioNet credentialed-access Data Use Agreement, and
3. Be approved for the corresponding projects on
   [PhysioNet](https://physionet.org/).

No raw or patient-level data may be redistributed. Keep the databases local and
share only aggregate outputs.

---

## Configuration: database credentials

Credentials are read from **environment variables** so that they are never
written into the code or committed to version control. Set them in your shell or
in a local (git-ignored) `.Renviron`:

```r
Sys.setenv(
  PGHOST     = "localhost",
  PGPORT     = "5432",
  PGDATABASE = "mimiciv",          # or "eicu" / "sicdb"
  PGUSER     = "postgres",
  PGPASSWORD = "your_password_here" # never commit this value
)
```

The R scripts read `Sys.getenv("PGPASSWORD")` with **no default**, so a connection
fails fast if the password is not set — the intended behaviour for a public
repository.

---

## How to run

1. **Build the cohort tables in PostgreSQL.** Run the MIMIC-IV SQL pipeline to
   create `analytic_cohort` (and `analytic_cohort_48h`); run
   `eicu/eICU_pipeline.sql` to create `eicu_analytic_cohort`; and run
   `sicdb/SICdb_pipeline.sql` to create `sic_cohort` / `sic_exposure`.

2. **Run the MIMIC-IV analysis** (primary first; the 48h script reads the
   corrected 24h CCW output for its side-by-side comparison):

   ```r
   source("mimic/01_MIMIC_24h.R")   # primary (24h)
   source("mimic/01_MIMIC_48h.R")   # sensitivity (48h)
   ```

3. **Run the eICU validation:**

   ```r
   source("eicu/eICU.R")
   ```

Runtime notes:

- **MICE** (`m = 20`) is slow on first run. Each script caches the imputation to
  an `.rds` file and reuses it automatically; delete that `.rds` to force
  recomputation if the input table changes.
- The **CCW cluster bootstrap** in the MIMIC-IV scripts is the other slow step.
  On Unix/macOS, set `N_CORES > 1` near the top of those scripts to parallelize
  it.

---

## Outputs

All outputs are **aggregate and safe to share**. Running the pipelines produces:

**Effect tables and diagnostics (CSV)** — the MIMIC-IV primary/secondary effect
tables, Table 1 (unweighted and weighted), negative-control and E-value tables,
the CCW point-estimate and bootstrap-CI tables, the bias-gradient table, the 48h
sensitivity tables, and the eICU master results table.

**Figures (PDF)** — covariate-balance Love plots, MICE convergence/density checks,
study flow chart, propensity-score overlap, cumulative-incidence and
Kaplan-Meier curves, and the primary / subgroup / eICU forest plots.

**Cached intermediates (RDS)** — imputed long tables and weighted datasets; these
are large, machine-generated, and excluded from version control by default.

The SQL pipelines additionally create their cohort tables in the database and
print aggregate validation counts to the console.

---

## Reproducibility

- **Seeds are fixed**, so results are deterministic given the same input tables
  and package versions.
- After a run, record the exact environment with `sessionInfo()` (or an
  [`renv`](https://rstudio.github.io/renv/) lockfile) and include it with any
  release.
- Because the input data are credentialed, exact reproduction requires your own
  approved MIMIC-IV, eICU, and SICdb installations and the upstream SQL pipelines.

---

## Notes and limitations

- In MIMIC-IV the exposure timestamp derives from prescription start times, which
  are coarse; the 48h-window sensitivity analysis is partly intended to be more
  forgiving of this timing granularity.
- eICU endpoints are in-hospital (follow-up ends at discharge), so the eICU and
  MIMIC-IV outcome definitions are aligned but not identical; the eICU horizon
  approximates 30 days.
- In eICU, RRT initiation is subject to strong between-hospital ecology; the
  hospital-stratified site correction and the RRT sensitivity analysis address it.
- In SICdb, H2RA exposure is very rare (ranitidine having been withdrawn), which
  limits a head-to-head PPI-vs-H2RA comparison in that database; SICdb also lacks
  a pre-admission medication history, so the new-user flag is an in-ICU
  approximation.
- The prevalent-user-bias modules (the new-user arms in every analysis) are
  central to interpretation: a signal concentrated in the prevalent-user
  population and absent under a new-user design points to prevalent-user bias
  rather than a causal effect.

---

## Disclaimer

This repository is provided for research and educational purposes only. It is
**not** medical advice and must not be used to guide clinical care. The code is
released "as is," without warranty of any kind.

---

## Citation

If you use this code, please cite the associated manuscript. This repository also
includes a `CITATION.cff`, so GitHub's **"Cite this repository"** button will
generate a formatted citation. A suggested entry (update once published):

```
[Authors]. Proton-pump inhibitors versus H2-receptor antagonists in
sepsis-associated acute kidney injury: a target-trial emulation of MIMIC-IV with
external validation in eICU and SICdb. [Journal], [Year]. [DOI]
```

Please also cite the MIMIC-IV, eICU, and SICdb databases and PhysioNet per their
requirements.

---

## License

Released under the [MIT License](LICENSE). Review and replace this with the
license appropriate for your project and institution before publishing, and fill
in the copyright holder in the `LICENSE` file.
