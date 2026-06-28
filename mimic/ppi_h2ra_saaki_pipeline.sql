-- ============================================================================
-- PPI vs H2RA in Sepsis-Associated Acute Kidney Injury (SA-AKI)
-- MIMIC-IV cohort-building and outcome-extraction pipeline
-- ----------------------------------------------------------------------------
-- Design : target-trial emulation, active-comparator, new-initiator (sensitivity)
-- Estimands : per-protocol via clone-censor-weight (CCW) + IPTW head-to-head
-- Primary outcome : MAKE30   Co-primary : MAKE90
-- Time zero (T0)  : onset of SA-AKI = first KDIGO>=1 within 7 d of sepsis (aki_time)
-- Engine : PostgreSQL (developed on pgAdmin 4), database `mimiciv`
--
-- HOW TO RUN
--   Execute this file top-to-bottom. Tables are created in dependency order and
--   keyed on `stay_id` throughout. Each module ends with a validation SELECT;
--   the "[ref]" values in those comments are from the authors' run and reproduce
--   exactly for the cohort-defining counts. A few outcome counts (flagged [d])
--   can differ by a handful of rows between MIMIC-IV point releases.
--
-- PREREQUISITES (must already exist in the database)
--   * MIMIC-IV core schemas: mimiciv_hosp, mimiciv_icu
--   * MIMIC-IV derived tables built with the official mimic-code repository
--     (schema mimiciv_derived): sepsis3, age, icustay_detail, kdigo_stages,
--     charlson, rrt, sofa, ventilation, vasoactive_agent, first_day_weight,
--     urine_output
--   * public.elixhauser_scores(hadm_id, elixhauser_score): van Walraven-weighted
--     Elixhauser score per admission (see README for the build script)
--
-- DATA-USE NOTE
--   MIMIC-IV is credentialed PhysioNet data governed by a Data Use Agreement.
--   Every table created here is patient-level and must remain on your machine.
--   Only aggregate output (counts, estimates, confidence intervals) may leave it.
--   Do NOT commit any data extract to this repository.
--
-- KEY IMPLEMENTATION CHOICES
--   * No ICU length-of-stay >=24 h filter (would select against early deaths).
--   * Death date = COALESCE(patients.dod, admissions.deathtime::date), so in-
--     hospital deaths are captured robustly.
--   * Dialysis dependence judged in the 48 h before min(discharge, horizon)
--     -> rrt_dep_30 / rrt_dep_90.
--   * MAKE time-to-event pair make_tte_time/_event (earlier of death or RRT
--     start) is provided for RMST analyses.
--   * Creatinine itemid = 50912; physiologic upper bound valuenum < 50 for
--     follow-up creatinine; non-recovery and recovery windows both start at T0.
--   * Censored survival times surv_28 / surv_90 are precomputed for R.
--
-- OUTPUT TABLES OF INTEREST
--   analytic_cohort       main analytic dataset (~21,731 rows x 121 cols)
--   analytic_cohort_48h   48 h exposure-window sensitivity dataset
-- ============================================================================


-- ############################################################################
-- PART A — COHORT CONSTRUCTION (steps A1-A6)
-- ############################################################################

-- ============================================================================
-- A1. Inclusion cohort: cohort_inclusion
--     age >= 18  AND  Sepsis-3  AND  first KDIGO>=1 within 7 d of sepsis (sets T0)
--     AND first ICU stay of the admission.
--     No ICU-LOS >= 24 h filter (avoids selection bias against early deaths).
-- ============================================================================
DROP TABLE IF EXISTS cohort_inclusion;
CREATE TABLE cohort_inclusion AS
WITH sepsis_raw AS (
    SELECT s3.stay_id, ie.subject_id, ie.hadm_id,
           s3.suspected_infection_time AS sepsis_time,
           ie.intime, ie.outtime,
           EXTRACT(EPOCH FROM (ie.outtime - ie.intime))/3600/24 AS icu_los_days
    FROM mimiciv_derived.sepsis3 s3
    INNER JOIN mimiciv_icu.icustays ie ON s3.stay_id = ie.stay_id
    WHERE s3.sepsis3 = TRUE AND s3.suspected_infection_time IS NOT NULL
),
with_age AS (
    SELECT sr.*, age.age AS age_at_icu
    FROM sepsis_raw sr
    INNER JOIN mimiciv_derived.age age ON sr.hadm_id = age.hadm_id
    WHERE age.age >= 18
),
first_stay AS (
    SELECT wa.*
    FROM with_age wa
    INNER JOIN mimiciv_derived.icustay_detail icd ON wa.stay_id = icd.stay_id
    WHERE icd.first_icu_stay = TRUE
),
aki_first AS (   -- first KDIGO>=1 over the whole stay (aki_stage = max of creat/UO staging, earliest time)
    SELECT stay_id, aki_time, first_aki_stage
    FROM (SELECT stay_id, charttime AS aki_time, aki_stage AS first_aki_stage,
                 ROW_NUMBER() OVER (PARTITION BY stay_id ORDER BY charttime) AS rn
          FROM mimiciv_derived.kdigo_stages WHERE aki_stage >= 1) k
    WHERE rn = 1
)
SELECT fs.subject_id, fs.hadm_id, fs.stay_id, fs.intime, fs.outtime, fs.icu_los_days,
       fs.age_at_icu, fs.sepsis_time, af.aki_time, af.first_aki_stage
FROM first_stay fs
INNER JOIN aki_first af ON fs.stay_id = af.stay_id
WHERE af.aki_time >= fs.sepsis_time
  AND af.aki_time <= fs.sepsis_time + INTERVAL '7 days';
SELECT 'A1_inclusion' AS step, COUNT(*) AS n FROM cohort_inclusion;   -- [ref] 23,745

-- ============================================================================
-- A2. Exclusions + final cohort: cohort_final
--     (1) ESRD / chronic dialysis   (2) kidney transplant   (3) pregnancy/puerperium
--     ESRD uses only "ESRD + dialysis-dependence" codes (not dialysis-encounter
--     codes V56/Z49), so patients on acute RRT this admission are not excluded.
-- ============================================================================
DROP TABLE IF EXISTS cohort_exclusion_flags;
CREATE TABLE cohort_exclusion_flags AS
WITH dx AS (
    SELECT ci.stay_id,
      MAX(CASE WHEN (d.icd_version=9  AND d.icd_code IN ('5856','V4511'))
                OR (d.icd_version=10 AND d.icd_code IN ('N186','Z992'))
           THEN 1 ELSE 0 END) AS excl_esrd_dialysis,
      MAX(CASE WHEN (d.icd_version=9  AND d.icd_code IN ('V420','99681'))
                OR (d.icd_version=10 AND (d.icd_code='Z940' OR d.icd_code LIKE 'T861%'))
           THEN 1 ELSE 0 END) AS excl_transplant,
      MAX(CASE WHEN (d.icd_version=9  AND (SUBSTR(d.icd_code,1,3) BETWEEN '630' AND '679'
                          OR d.icd_code LIKE 'V22%' OR d.icd_code LIKE 'V23%'
                          OR d.icd_code LIKE 'V24%' OR d.icd_code LIKE 'V27%' OR d.icd_code LIKE 'V28%'))
                OR (d.icd_version=10 AND (d.icd_code LIKE 'O%' OR d.icd_code LIKE 'Z33%'
                          OR d.icd_code LIKE 'Z34%' OR d.icd_code LIKE 'Z3A%'
                          OR d.icd_code LIKE 'Z37%' OR d.icd_code LIKE 'Z39%'))
           THEN 1 ELSE 0 END) AS excl_pregnancy
    FROM cohort_inclusion ci
    LEFT JOIN mimiciv_hosp.diagnoses_icd d ON ci.hadm_id = d.hadm_id
    GROUP BY ci.stay_id
)
SELECT ci.*, COALESCE(dx.excl_esrd_dialysis,0) AS excl_esrd_dialysis,
       COALESCE(dx.excl_transplant,0) AS excl_transplant,
       COALESCE(dx.excl_pregnancy,0)  AS excl_pregnancy
FROM cohort_inclusion ci LEFT JOIN dx ON ci.stay_id = dx.stay_id;
DROP TABLE IF EXISTS cohort_final;
CREATE TABLE cohort_final AS
SELECT subject_id, hadm_id, stay_id, intime, outtime, icu_los_days,
       age_at_icu, sepsis_time, aki_time, first_aki_stage
FROM cohort_exclusion_flags
WHERE excl_esrd_dialysis=0 AND excl_transplant=0 AND excl_pregnancy=0;
SELECT COUNT(*) AS n_inclusion,                                    -- [ref] 23,745
       SUM(excl_esrd_dialysis) AS n_esrd,                         -- [ref]  1,542
       SUM(excl_transplant)    AS n_transplant,                   -- [ref]    475
       SUM(excl_pregnancy)     AS n_pregnancy,                    -- [ref]     55
       COUNT(*) FILTER (WHERE excl_esrd_dialysis=0 AND excl_transplant=0
                          AND excl_pregnancy=0) AS n_final         -- [ref] 21,883
FROM cohort_exclusion_flags;

-- ============================================================================
-- A3. Baseline creatinine: cohort_baseline_creat
--     Three-tier priority (creatinine itemid=50912, floor 0.2 mg/dL):
--       (1) lowest pre-admission value 7-365 d before admittime
--       (2) lowest value in the first 24 h of ICU
--       (3) back-calculated from MDRD assuming eGFR = 75
-- ============================================================================
DROP TABLE IF EXISTS cohort_baseline_creat;
CREATE TABLE cohort_baseline_creat AS
WITH base AS (
    SELECT cf.*, adm.admittime, pat.gender
    FROM cohort_final cf
    INNER JOIN mimiciv_hosp.admissions adm ON cf.hadm_id = adm.hadm_id
    INNER JOIN mimiciv_hosp.patients   pat ON cf.subject_id = pat.subject_id
),
creat_pre AS (
    SELECT b.stay_id, ROUND(MIN(le.valuenum)::numeric,2) AS creat_pre
    FROM base b JOIN mimiciv_hosp.labevents le
      ON le.subject_id=b.subject_id AND le.itemid=50912
     AND le.valuenum IS NOT NULL AND le.valuenum>=0.2 AND le.valuenum<20
     AND le.charttime BETWEEN b.admittime - INTERVAL '365 days' AND b.admittime - INTERVAL '7 days'
    GROUP BY b.stay_id
),
creat_24h AS (
    SELECT b.stay_id, ROUND(MIN(le.valuenum)::numeric,2) AS creat_24h
    FROM base b JOIN mimiciv_hosp.labevents le
      ON le.subject_id=b.subject_id AND le.itemid=50912
     AND le.valuenum IS NOT NULL AND le.valuenum>=0.2 AND le.valuenum<20
     AND le.charttime BETWEEN b.intime AND b.intime + INTERVAL '24 hours'
    GROUP BY b.stay_id
)
SELECT b.subject_id, b.hadm_id, b.stay_id, b.intime, b.outtime, b.icu_los_days,
       b.age_at_icu, b.gender, b.sepsis_time, b.aki_time, b.first_aki_stage,
       cp.creat_pre, c24.creat_24h,
       ROUND((POWER(175.0*POWER(b.age_at_icu,-0.203)
              *(CASE WHEN b.gender='F' THEN 0.742 ELSE 1.0 END)/75.0, 1.0/1.154))::numeric,2) AS creat_mdrd,
       COALESCE(cp.creat_pre, c24.creat_24h,
                ROUND((POWER(175.0*POWER(b.age_at_icu,-0.203)
                       *(CASE WHEN b.gender='F' THEN 0.742 ELSE 1.0 END)/75.0,1.0/1.154))::numeric,2)) AS baseline_creat,
       CASE WHEN cp.creat_pre IS NOT NULL THEN 'preadmission'
            WHEN c24.creat_24h IS NOT NULL THEN 'first24h_icu'
            ELSE 'mdrd_egfr75' END AS baseline_creat_source
FROM base b
LEFT JOIN creat_pre cp ON cp.stay_id=b.stay_id
LEFT JOIN creat_24h c24 ON c24.stay_id=b.stay_id;
SELECT baseline_creat_source, COUNT(*) AS n,
       ROUND(AVG(baseline_creat)::numeric,2) AS mean_base,
       ROUND(MIN(baseline_creat)::numeric,2) AS min_base
FROM cohort_baseline_creat GROUP BY baseline_creat_source ORDER BY n DESC;
-- [ref] first24h_icu ~11,083 (1.29) | preadmission ~10,739 (0.90) | mdrd_egfr 61 (0.88)

-- ============================================================================
-- A4. Baseline eGFR + CKD flag: cohort_baseline_egfr
--     eGFR via CKD-EPI 2021 (race-free). CKD = baseline eGFR < 60 OR a CKD
--     diagnosis code (ICD-9 585x / ICD-10 N18x).
-- ============================================================================
DROP TABLE IF EXISTS cohort_baseline_egfr;
CREATE TABLE cohort_baseline_egfr AS
WITH ckd_dx AS (
    SELECT bc.stay_id,
      MAX(CASE WHEN (d.icd_version=9 AND d.icd_code LIKE '585%')
                OR (d.icd_version=10 AND d.icd_code LIKE 'N18%') THEN 1 ELSE 0 END) AS ckd_dx
    FROM cohort_baseline_creat bc
    LEFT JOIN mimiciv_hosp.diagnoses_icd d ON bc.hadm_id = d.hadm_id
    GROUP BY bc.stay_id
),
egfr AS (
    SELECT bc.*,
      ROUND((142.0
        * POWER(LEAST   (bc.baseline_creat/(CASE WHEN bc.gender='F' THEN 0.7 ELSE 0.9 END),1.0),
                (CASE WHEN bc.gender='F' THEN -0.241 ELSE -0.302 END))
        * POWER(GREATEST(bc.baseline_creat/(CASE WHEN bc.gender='F' THEN 0.7 ELSE 0.9 END),1.0), -1.200)
        * POWER(0.9938, bc.age_at_icu)
        * (CASE WHEN bc.gender='F' THEN 1.012 ELSE 1.0 END))::numeric,1) AS baseline_egfr
    FROM cohort_baseline_creat bc
)
SELECT e.*, COALESCE(cd.ckd_dx,0) AS ckd_dx,
       CASE WHEN e.baseline_egfr<60 OR COALESCE(cd.ckd_dx,0)=1 THEN 1 ELSE 0 END AS ckd
FROM egfr e LEFT JOIN ckd_dx cd ON cd.stay_id=e.stay_id;
SELECT COUNT(*) AS n, ROUND(AVG(baseline_egfr)::numeric,1) AS mean_egfr,   -- [ref] 21,883 / ~78.6
       SUM(ckd) AS n_ckd, ROUND(100.0*SUM(ckd)/COUNT(*),1) AS pct_ckd      -- [ref] ~7,653 / 35%
FROM cohort_baseline_egfr;

-- ============================================================================
-- A5. Exposure PPI vs H2RA: cohort_exposure  (source = prescriptions = treatment decision)
--     exposure = class of the FIRST acid suppressant in [T0, T0+24h];
--     newuser  = no acid-suppressant order before T0 this admission.
--     Drugs matched on full generic names (so -prazole/-tidine do not catch
--     aripiprazole / azacitidine). newuser is operationalized as "no prior acid-
--     suppressant order this admission" (~ SAP "none in 7 d before T0"; MIMIC has
--     no pre-admission meds and prescriptions.starttime is the initiation time).
-- ============================================================================
DROP TABLE IF EXISTS cohort_exposure;
CREATE TABLE cohort_exposure AS
WITH acid AS (
    SELECT e.stay_id, e.aki_time, rx.starttime, rx.drug,
      CASE WHEN rx.drug ILIKE '%pantoprazole%' OR rx.drug ILIKE '%omeprazole%'
             OR rx.drug ILIKE '%esomeprazole%' OR rx.drug ILIKE '%lansoprazole%'
             OR rx.drug ILIKE '%rabeprazole%'  OR rx.drug ILIKE '%dexlansoprazole%'
           THEN 'PPI' ELSE 'H2RA' END AS acid_class
    FROM cohort_baseline_egfr e
    JOIN mimiciv_hosp.prescriptions rx
      ON rx.hadm_id=e.hadm_id AND rx.starttime IS NOT NULL
     AND (rx.drug ILIKE '%pantoprazole%' OR rx.drug ILIKE '%omeprazole%'
       OR rx.drug ILIKE '%esomeprazole%' OR rx.drug ILIKE '%lansoprazole%'
       OR rx.drug ILIKE '%rabeprazole%'  OR rx.drug ILIKE '%dexlansoprazole%'
       OR rx.drug ILIKE '%famotidine%'   OR rx.drug ILIKE '%ranitidine%'
       OR rx.drug ILIKE '%cimetidine%'   OR rx.drug ILIKE '%nizatidine%')
),
first_in_window AS (
    SELECT DISTINCT ON (stay_id)
        stay_id, acid_class AS first_acid_class, drug AS first_acid_drug, starttime AS first_acid_time
    FROM acid
    WHERE starttime >= aki_time AND starttime <= aki_time + INTERVAL '24 hours'
    ORDER BY stay_id, starttime ASC, acid_class
),
both_flag AS (
    SELECT stay_id, CASE WHEN COUNT(DISTINCT acid_class)>1 THEN 1 ELSE 0 END AS both_in_window
    FROM acid WHERE starttime>=aki_time AND starttime<=aki_time+INTERVAL '24 hours' GROUP BY stay_id
),
exact_tie_flag AS (
    SELECT stay_id, CASE WHEN COUNT(DISTINCT acid_class)>1 THEN 1 ELSE 0 END AS exact_tie
    FROM (SELECT stay_id, acid_class, starttime, MIN(starttime) OVER (PARTITION BY stay_id) AS min_st
          FROM acid WHERE starttime>=aki_time AND starttime<=aki_time+INTERVAL '24 hours') w
    WHERE starttime=min_st GROUP BY stay_id
),
prior_flag AS (
    SELECT DISTINCT stay_id, 1 AS prior_acid FROM acid WHERE starttime < aki_time
)
SELECT e.*,
       COALESCE(fw.first_acid_class,'neither') AS exposure_group,
       fw.first_acid_drug, fw.first_acid_time,
       CASE WHEN fw.first_acid_time IS NOT NULL
            THEN ROUND((EXTRACT(EPOCH FROM (fw.first_acid_time-e.aki_time))/3600.0)::numeric,1) END AS hours_to_exposure,
       COALESCE(bf.both_in_window,0) AS both_in_window,
       COALESCE(et.exact_tie,0)      AS exact_tie,
       COALESCE(pf.prior_acid,0)     AS prior_acid,
       CASE WHEN COALESCE(pf.prior_acid,0)=0 THEN 1 ELSE 0 END AS newuser
FROM cohort_baseline_egfr e
LEFT JOIN first_in_window fw ON fw.stay_id=e.stay_id
LEFT JOIN both_flag bf ON bf.stay_id=e.stay_id
LEFT JOIN exact_tie_flag et ON et.stay_id=e.stay_id
LEFT JOIN prior_flag pf ON pf.stay_id=e.stay_id;
SELECT exposure_group, COUNT(*) AS n, SUM(newuser) AS n_newuser,
       SUM(both_in_window) AS n_both, SUM(exact_tie) AS n_tie
FROM cohort_exposure GROUP BY exposure_group ORDER BY n DESC;
-- [ref] neither 15,351 (nu 4,050) | PPI 3,449 (nu 1,535) | H2RA 3,083 (nu 1,363)

-- ============================================================================
-- A6. New-initiator base cohort: cohort_newuser (for IPTW head-to-head / CCW)
-- ============================================================================
DROP TABLE IF EXISTS cohort_newuser;
CREATE TABLE cohort_newuser AS SELECT * FROM cohort_exposure WHERE newuser=1;
SELECT COUNT(*) AS n_newuser,                                              -- [ref] 6,948
       SUM((exposure_group='PPI')::int)  AS n_ppi,                        -- [ref] 1,535
       SUM((exposure_group='H2RA')::int) AS n_h2ra,                       -- [ref] 1,363
       SUM((exposure_group='neither')::int) AS n_neither                  -- [ref] 4,050
FROM cohort_newuser;


-- ############################################################################
-- PART B — OUTCOME EXTRACTION (steps B1-B8; full cohort_exposure -> cohort_outcome)
--   The clock starts at T0 (not the dosing time) to avoid immortal-time bias;
--   CCW handles the T0 -> first-dose grace period downstream.
--   Kidney outcomes are censored at hospital discharge; death is a competing event.
-- ############################################################################

-- ============================================================================
-- B1. cohort_outcome_base: death + follow-up / censoring windows
--     death_date = COALESCE(dod, deathtime::date); death_icu flag added.
-- ============================================================================
DROP TABLE IF EXISTS cohort_outcome_base;
CREATE TABLE cohort_outcome_base AS
WITH dd AS (
    SELECT ce.subject_id, ce.hadm_id, ce.stay_id, ce.aki_time, ce.outtime AS icu_outtime,
           adm.admittime, adm.dischtime, adm.deathtime, pat.dod,
           COALESCE(pat.dod, adm.deathtime::date) AS death_date
    FROM cohort_exposure ce
    INNER JOIN mimiciv_hosp.admissions adm ON ce.hadm_id    = adm.hadm_id
    INNER JOIN mimiciv_hosp.patients   pat ON ce.subject_id = pat.subject_id
)
SELECT
    subject_id, hadm_id, stay_id, aki_time AS t0,
    admittime, dischtime, deathtime, dod, icu_outtime, death_date,
    (death_date - aki_time::date) AS days_to_death,
    CASE WHEN death_date IS NOT NULL AND (death_date-aki_time::date)<=28 THEN 1 ELSE 0 END AS death_28,
    CASE WHEN death_date IS NOT NULL AND (death_date-aki_time::date)<=30 THEN 1 ELSE 0 END AS death_30,
    CASE WHEN death_date IS NOT NULL AND (death_date-aki_time::date)<=90 THEN 1 ELSE 0 END AS death_90,
    CASE WHEN deathtime IS NOT NULL THEN 1 ELSE 0 END AS death_inhosp,
    CASE WHEN deathtime IS NOT NULL AND deathtime <= icu_outtime THEN 1 ELSE 0 END AS death_icu,
    ROUND((EXTRACT(EPOCH FROM (dischtime-aki_time))/86400.0)::numeric,2) AS los_after_t0_days,
    CASE WHEN dischtime < aki_time + INTERVAL '30 days' THEN 1 ELSE 0 END AS disch_before_d30,
    CASE WHEN dischtime < aki_time + INTERVAL '90 days' THEN 1 ELSE 0 END AS disch_before_d90
FROM dd;
-- clean impossible timestamps (death/discharge before T0; null discharge)
DELETE FROM cohort_outcome_base
WHERE days_to_death < 0 OR dischtime < t0 OR dischtime IS NULL;
SELECT COUNT(*) AS n_clean FROM cohort_outcome_base;                       -- [ref] ~21,828
SELECT SUM(death_inhosp) AS n_inhosp, SUM(death_icu) AS n_icu,             -- [ref] ~3,996 / ~2,808
       SUM(death_28) AS n_d28, SUM(death_90) AS n_d90,                     -- [ref] ~4,909 / ~6,421
       ROUND(100.0*SUM(death_28)/COUNT(*),1) AS pct_d28                    -- [ref] ~22.5
FROM cohort_outcome_base;

-- ============================================================================
-- B2. cohort_rrt: RRT initiation + dialysis dependence
--     source mimiciv_derived.rrt with dialysis_active = 1.
--     Dependence window = [min(discharge, horizon) - 48h, min(discharge, horizon)]
--     -> rrt_dep_30 / rrt_dep_90.
-- ============================================================================
DROP TABLE IF EXISTS cohort_rrt;
CREATE TABLE cohort_rrt AS
WITH base AS (SELECT stay_id, t0, dischtime FROM cohort_outcome_base),
rrt_active AS (
    SELECT r.stay_id, r.charttime
    FROM mimiciv_derived.rrt r INNER JOIN base b ON r.stay_id=b.stay_id
    WHERE r.dialysis_active = 1
),
agg AS (
    SELECT b.stay_id, b.t0, b.dischtime,
      MAX((ra.charttime IS NOT NULL)::int) AS rrt_ever,
      MAX((ra.charttime <  b.t0)::int)     AS rrt_pre_t0,
      MIN(CASE WHEN ra.charttime > b.t0 THEN ra.charttime END) AS rrt_init_time,
      MAX(CASE WHEN ra.charttime >= LEAST(b.dischtime, b.t0+INTERVAL '30 days') - INTERVAL '48 hours'
                AND ra.charttime <= LEAST(b.dischtime, b.t0+INTERVAL '30 days') THEN 1 ELSE 0 END) AS rrt_dep_30,
      MAX(CASE WHEN ra.charttime >= LEAST(b.dischtime, b.t0+INTERVAL '90 days') - INTERVAL '48 hours'
                AND ra.charttime <= LEAST(b.dischtime, b.t0+INTERVAL '90 days') THEN 1 ELSE 0 END) AS rrt_dep_90
    FROM base b LEFT JOIN rrt_active ra ON ra.stay_id=b.stay_id
    GROUP BY b.stay_id, b.t0, b.dischtime
)
SELECT stay_id,
       COALESCE(rrt_ever,0)   AS rrt_ever,
       COALESCE(rrt_pre_t0,0) AS rrt_pre_t0,
       rrt_init_time,
       CASE WHEN rrt_init_time IS NOT NULL
            THEN ROUND((EXTRACT(EPOCH FROM (rrt_init_time-t0))/86400.0)::numeric,2) END AS days_to_rrt,
       CASE WHEN rrt_init_time IS NOT NULL THEN 1 ELSE 0 END AS rrt_init,
       COALESCE(rrt_dep_30,0) AS rrt_dep_30,
       COALESCE(rrt_dep_90,0) AS rrt_dep_90
FROM agg;
SELECT SUM(rrt_ever) AS n_ever, SUM(rrt_pre_t0) AS n_pre, SUM(rrt_init) AS n_init,   -- [ref] ~1,242 / 103 / 1,227
       SUM(rrt_dep_30) AS n_dep30, SUM(rrt_dep_90) AS n_dep90,                       -- [ref] ~621 / ~600
       PERCENTILE_CONT(0.5) WITHIN GROUP (ORDER BY days_to_rrt) AS med_days_to_rrt   -- [ref] ~1.49
FROM cohort_rrt;

-- ============================================================================
-- B3. cohort_nonrecovery: creatinine non-recovery (MAKE component c); source labevents 50912
--     Creatinine window >= T0 and valuenum < 50; last-creatinine eGFR uses the
--     same CKD-EPI 2021. Non-recovery = last creatinine > 1.5x baseline OR
--     eGFR < 0.75 x baseline eGFR; no creatinine -> NULL.
-- ============================================================================
DROP TABLE IF EXISTS cohort_nonrecovery;
CREATE TABLE cohort_nonrecovery AS
WITH base AS (
    SELECT b.subject_id, b.stay_id, b.t0, b.dischtime,
           ce.baseline_creat, ce.baseline_egfr, ce.age_at_icu, ce.gender
    FROM cohort_outcome_base b
    INNER JOIN cohort_exposure ce ON b.stay_id = ce.stay_id
),
creat AS (
    SELECT b.stay_id, b.t0, le.charttime, le.valuenum
    FROM base b INNER JOIN mimiciv_hosp.labevents le
      ON le.subject_id=b.subject_id AND le.itemid=50912
     AND le.valuenum>0 AND le.valuenum<50
     AND le.charttime>=b.t0 AND le.charttime<=b.dischtime
),
last30 AS (SELECT DISTINCT ON (stay_id) stay_id, valuenum AS last_creat_30
           FROM creat WHERE charttime<=t0+INTERVAL '30 days' ORDER BY stay_id, charttime DESC),
last90 AS (SELECT DISTINCT ON (stay_id) stay_id, valuenum AS last_creat_90
           FROM creat WHERE charttime<=t0+INTERVAL '90 days' ORDER BY stay_id, charttime DESC),
joined AS (
    SELECT b.stay_id, b.baseline_creat, b.baseline_egfr, b.age_at_icu, b.gender,
           l30.last_creat_30, l90.last_creat_90
    FROM base b LEFT JOIN last30 l30 ON l30.stay_id=b.stay_id
                LEFT JOIN last90 l90 ON l90.stay_id=b.stay_id
),
with_egfr AS (
    SELECT j.*,
      CASE WHEN last_creat_30 IS NOT NULL THEN ROUND((142.0
        * POWER(LEAST   (last_creat_30/(CASE WHEN gender='F' THEN 0.7 ELSE 0.9 END),1.0), CASE WHEN gender='F' THEN -0.241 ELSE -0.302 END)
        * POWER(GREATEST(last_creat_30/(CASE WHEN gender='F' THEN 0.7 ELSE 0.9 END),1.0), -1.200)
        * POWER(0.9938, age_at_icu) * (CASE WHEN gender='F' THEN 1.012 ELSE 1.0 END))::numeric,1) END AS egfr_30,
      CASE WHEN last_creat_90 IS NOT NULL THEN ROUND((142.0
        * POWER(LEAST   (last_creat_90/(CASE WHEN gender='F' THEN 0.7 ELSE 0.9 END),1.0), CASE WHEN gender='F' THEN -0.241 ELSE -0.302 END)
        * POWER(GREATEST(last_creat_90/(CASE WHEN gender='F' THEN 0.7 ELSE 0.9 END),1.0), -1.200)
        * POWER(0.9938, age_at_icu) * (CASE WHEN gender='F' THEN 1.012 ELSE 1.0 END))::numeric,1) END AS egfr_90
    FROM joined j
)
SELECT stay_id, baseline_creat, baseline_egfr, last_creat_30, last_creat_90, egfr_30, egfr_90,
       CASE WHEN last_creat_30 IS NULL THEN NULL
            WHEN last_creat_30>1.5*baseline_creat OR egfr_30<0.75*baseline_egfr THEN 1 ELSE 0 END AS non_recovery_30,
       CASE WHEN last_creat_90 IS NULL THEN NULL
            WHEN last_creat_90>1.5*baseline_creat OR egfr_90<0.75*baseline_egfr THEN 1 ELSE 0 END AS non_recovery_90
FROM with_egfr;
SELECT SUM((non_recovery_30=1)::int) AS n_nonrec30,                        -- [ref] ~5,501
       SUM((non_recovery_30=0)::int) AS n_rec30,                           -- [ref] ~16,024
       SUM((non_recovery_30 IS NULL)::int) AS n_missing30                  -- [ref] ~303
FROM cohort_nonrecovery;

-- ============================================================================
-- B4. cohort_make: assemble MAKE30 / MAKE90 (uses rrt_dep_30/90)
--     MAKE = death OR dialysis dependence OR creatinine non-recovery, by horizon.
--     NULL only when the non-recovery component is unknown (no follow-up creat).
-- ============================================================================
DROP TABLE IF EXISTS cohort_make;
CREATE TABLE cohort_make AS
SELECT
    b.stay_id, b.death_30, b.death_90,
    r.rrt_dep_30, r.rrt_dep_90,
    nr.non_recovery_30, nr.non_recovery_90,
    CASE WHEN b.death_30=1 OR r.rrt_dep_30=1 OR nr.non_recovery_30=1 THEN 1
         WHEN nr.non_recovery_30 IS NULL THEN NULL ELSE 0 END AS make30,
    CASE WHEN b.death_90=1 OR r.rrt_dep_90=1 OR nr.non_recovery_90=1 THEN 1
         WHEN nr.non_recovery_90 IS NULL THEN NULL ELSE 0 END AS make90
FROM cohort_outcome_base b
INNER JOIN cohort_rrt         r  ON r.stay_id  = b.stay_id
INNER JOIN cohort_nonrecovery nr ON nr.stay_id = b.stay_id;
SELECT SUM((make30=1)::int) AS n_make30,                                   -- [ref] ~8,187
       ROUND(100.0*SUM((make30=1)::int)/NULLIF(SUM((make30 IS NOT NULL)::int),0),1) AS pct_make30,  -- [ref] ~37.7
       SUM((make90=1)::int) AS n_make90,                                   -- [ref] ~9,178
       ROUND(100.0*SUM((make90=1)::int)/NULLIF(SUM((make90 IS NOT NULL)::int),0),1) AS pct_make90    -- [ref] ~42.2
FROM cohort_make;

-- ============================================================================
-- B5. cohort_progression: AKI progression to KDIGO 3 (only those < stage 3 at T0;
--     death is a competing event -> analyze with cause-specific HR).
-- ============================================================================
DROP TABLE IF EXISTS cohort_progression;
CREATE TABLE cohort_progression AS
WITH base AS (
    SELECT b.stay_id, b.t0, b.dischtime, ce.first_aki_stage
    FROM cohort_outcome_base b INNER JOIN cohort_exposure ce ON b.stay_id=ce.stay_id
),
stage3 AS (
    SELECT b.stay_id, MIN(k.charttime) AS prog_time
    FROM base b INNER JOIN mimiciv_derived.kdigo_stages k
      ON k.stay_id=b.stay_id AND k.aki_stage=3
     AND k.charttime>b.t0 AND k.charttime<=b.dischtime
    WHERE b.first_aki_stage<3
    GROUP BY b.stay_id
)
SELECT b.stay_id, b.first_aki_stage,
       CASE WHEN b.first_aki_stage<3 THEN 1 ELSE 0 END AS prog_eligible,
       s.prog_time,
       CASE WHEN s.prog_time IS NOT NULL
            THEN ROUND((EXTRACT(EPOCH FROM (s.prog_time-b.t0))/86400.0)::numeric,2) END AS days_to_prog,
       CASE WHEN s.prog_time IS NOT NULL THEN 1 ELSE 0 END AS prog_event
FROM base b LEFT JOIN stage3 s ON s.stay_id=b.stay_id;
SELECT SUM(prog_eligible) AS n_elig, SUM(prog_event) AS n_prog,            -- [ref] ~21,361 / ~4,990
       ROUND(100.0*SUM(prog_event)/NULLIF(SUM(prog_eligible),0),1) AS pct  -- [ref] ~23.4
FROM cohort_progression;

-- ============================================================================
-- B6. cohort_recovery: renal recovery (secondary; death competing); only patients
--     with creatinine-based AKI. Definition: first post-T0 creatinine <= 1.5x
--     baseline that is sustained (no later rebound) with no active RRT in the
--     prior 48 h. recovery_eligible = ever > 1.5x baseline.
-- ============================================================================
DROP TABLE IF EXISTS cohort_recovery;
CREATE TABLE cohort_recovery AS
WITH base AS (
    SELECT b.subject_id, b.stay_id, b.t0, b.dischtime, ce.baseline_creat
    FROM cohort_outcome_base b INNER JOIN cohort_exposure ce ON b.stay_id=ce.stay_id
),
creat AS (
    SELECT b.stay_id, le.charttime, le.valuenum,
           CASE WHEN le.valuenum<=1.5*b.baseline_creat THEN 1 ELSE 0 END AS at_goal
    FROM base b INNER JOIN mimiciv_hosp.labevents le
      ON le.subject_id=b.subject_id AND le.itemid=50912
     AND le.valuenum>0 AND le.valuenum<50
     AND le.charttime>=b.t0 AND le.charttime<=b.dischtime
),
last_above AS (SELECT stay_id, MAX(charttime) AS last_above_time FROM creat WHERE at_goal=0 GROUP BY stay_id),
rrt_active AS (
    SELECT r.stay_id, r.charttime AS rrt_time
    FROM mimiciv_derived.rrt r INNER JOIN base b ON r.stay_id=b.stay_id WHERE r.dialysis_active=1
),
cand AS (
    SELECT DISTINCT ON (c.stay_id) c.stay_id, c.charttime AS rec_time
    FROM creat c INNER JOIN last_above la ON la.stay_id=c.stay_id
    WHERE c.at_goal=1 AND c.charttime>la.last_above_time
      AND NOT EXISTS (SELECT 1 FROM rrt_active ra
                      WHERE ra.stay_id=c.stay_id
                        AND ra.rrt_time>=c.charttime-INTERVAL '48 hours' AND ra.rrt_time<=c.charttime)
    ORDER BY c.stay_id, c.charttime ASC
)
SELECT b.stay_id,
       CASE WHEN la.last_above_time IS NOT NULL THEN 1 ELSE 0 END AS recovery_eligible,
       cd.rec_time,
       CASE WHEN cd.rec_time IS NOT NULL
            THEN ROUND((EXTRACT(EPOCH FROM (cd.rec_time-b.t0))/86400.0)::numeric,2) END AS days_to_recovery,
       CASE WHEN cd.rec_time IS NOT NULL THEN 1 ELSE 0 END AS recovery_event
FROM base b
LEFT JOIN last_above la ON la.stay_id=b.stay_id
LEFT JOIN cand cd       ON cd.stay_id=b.stay_id;
SELECT SUM(recovery_eligible) AS n_elig,                                   -- [ref] ~8,282
       SUM(CASE WHEN recovery_eligible=1 THEN recovery_event ELSE 0 END) AS n_recov,   -- [ref] ~4,101
       PERCENTILE_CONT(0.5) WITHIN GROUP (ORDER BY days_to_recovery) AS med_days,      -- [ref] ~4.1
       SUM((recovery_eligible=0)::int) AS n_never_elevated                 -- [ref] ~13,546
FROM cohort_recovery;

-- ============================================================================
-- B7. cohort_nco: negative-control outcomes — incident VTE (primary NCO) +
--     stage>=2 pressure injury (secondary NCO); source diagnoses_icd.
--     Limitation: diagnosis codes carry no in-hospital timestamp, so "post-T0"
--     VTE is approximated as "occurring this admission with no prior history".
-- ============================================================================
DROP TABLE IF EXISTS cohort_nco;
CREATE TABLE cohort_nco AS
WITH base AS (SELECT subject_id, hadm_id, stay_id FROM cohort_outcome_base),
adm AS (SELECT b.stay_id, a.admittime FROM base b INNER JOIN mimiciv_hosp.admissions a ON a.hadm_id=b.hadm_id),
vte_index AS (
    SELECT DISTINCT b.stay_id
    FROM base b INNER JOIN mimiciv_hosp.diagnoses_icd d ON d.hadm_id=b.hadm_id
    WHERE (d.icd_version=9  AND (d.icd_code LIKE '4151%' OR d.icd_code LIKE '4534%'
                               OR d.icd_code LIKE '4538%' OR d.icd_code IN ('4532','4539')))
       OR (d.icd_version=10 AND (d.icd_code LIKE 'I26%'  OR d.icd_code LIKE 'I824%'
                               OR d.icd_code LIKE 'I826%' OR d.icd_code LIKE 'I822%'))
),
vte_prior AS (
    SELECT DISTINCT b.stay_id
    FROM base b INNER JOIN adm ON adm.stay_id=b.stay_id
    INNER JOIN mimiciv_hosp.admissions a2 ON a2.subject_id=b.subject_id AND a2.admittime<adm.admittime
    INNER JOIN mimiciv_hosp.diagnoses_icd d ON d.hadm_id=a2.hadm_id
    WHERE (d.icd_version=9  AND (d.icd_code LIKE '4151%' OR d.icd_code LIKE '4534%'
                               OR d.icd_code LIKE '4538%' OR d.icd_code IN ('4532','4539')))
       OR (d.icd_version=10 AND (d.icd_code LIKE 'I26%'  OR d.icd_code LIKE 'I824%'
                               OR d.icd_code LIKE 'I826%' OR d.icd_code LIKE 'I822%'))
),
pi_index AS (
    SELECT DISTINCT b.stay_id
    FROM base b INNER JOIN mimiciv_hosp.diagnoses_icd d ON d.hadm_id=b.hadm_id
    WHERE (d.icd_version=9  AND d.icd_code IN ('70722','70723','70724','70725'))
       OR (d.icd_version=10 AND d.icd_code LIKE 'L89%' AND RIGHT(TRIM(d.icd_code),1) IN ('2','3','4'))
)
SELECT b.stay_id,
       CASE WHEN vi.stay_id IS NOT NULL AND vp.stay_id IS NULL THEN 1 ELSE 0 END AS nco_vte,
       CASE WHEN pi.stay_id IS NOT NULL THEN 1 ELSE 0 END AS nco_pressure
FROM base b
LEFT JOIN vte_index vi ON vi.stay_id=b.stay_id
LEFT JOIN vte_prior vp ON vp.stay_id=b.stay_id
LEFT JOIN pi_index  pi ON pi.stay_id=b.stay_id;
SELECT SUM(nco_vte) AS n_vte, SUM(nco_pressure) AS n_pi FROM cohort_nco;   -- [ref] ~1,403 / ~916

-- ============================================================================
-- B8. cohort_outcome: assemble all outcomes + drop rows with MAKE fully missing.
--     make_tte_time/_event = earlier of death or RRT start (for RMST).
--     surv_28 / surv_90 = survival time censored at tau (event = death_28/90).
-- ============================================================================
DROP TABLE IF EXISTS cohort_outcome;
CREATE TABLE cohort_outcome AS
SELECT
    b.subject_id, b.hadm_id, b.stay_id, b.t0,
    b.dischtime, b.days_to_death, b.los_after_t0_days, b.disch_before_d30, b.disch_before_d90,
    -- death
    b.death_28, b.death_30, b.death_90, b.death_inhosp, b.death_icu,
    LEAST(COALESCE(b.days_to_death,28),28) AS surv_28,    -- event = death_28
    LEAST(COALESCE(b.days_to_death,90),90) AS surv_90,    -- event = death_90
    -- MAKE (binary + time-to-event)
    m.make30, m.make90,
    CASE WHEN b.death_date IS NOT NULL OR r.rrt_init=1 THEN 1 ELSE 0 END AS make_tte_event,
    CASE WHEN b.death_date IS NOT NULL OR r.rrt_init=1
         THEN LEAST(COALESCE(b.days_to_death,1e9), COALESCE(r.days_to_rrt,1e9)) END AS make_tte_time,
    -- RRT
    r.rrt_init, r.days_to_rrt, r.rrt_dep_30, r.rrt_dep_90, r.rrt_pre_t0,
    -- non-recovery
    nr.non_recovery_30, nr.non_recovery_90,
    -- progression / recovery
    pr.prog_eligible, pr.prog_event, pr.days_to_prog,
    rc.recovery_eligible, rc.recovery_event, rc.days_to_recovery,
    -- negative controls
    nc.nco_vte, nc.nco_pressure
FROM cohort_outcome_base b
INNER JOIN cohort_make        m  ON m.stay_id  = b.stay_id
INNER JOIN cohort_rrt         r  ON r.stay_id  = b.stay_id
INNER JOIN cohort_nonrecovery nr ON nr.stay_id = b.stay_id
INNER JOIN cohort_progression pr ON pr.stay_id = b.stay_id
INNER JOIN cohort_recovery    rc ON rc.stay_id = b.stay_id
INNER JOIN cohort_nco         nc ON nc.stay_id = b.stay_id
WHERE NOT (m.make30 IS NULL AND m.make90 IS NULL);
-- final N + outcome overview
SELECT COUNT(*) AS n_final,                                                -- [ref] ~21,731
       ROUND(100.0*AVG(make30::numeric),1) AS pct_make30,                 -- [ref] ~37.7
       ROUND(100.0*AVG(make90::numeric),1) AS pct_make90,                 -- [ref] ~42.2
       ROUND(100.0*AVG(death_28::numeric),1) AS pct_d28,                  -- [ref] ~22.6
       ROUND(100.0*AVG(death_90::numeric),1) AS pct_d90,                  -- [ref] ~29.5
       ROUND(100.0*AVG(rrt_init::numeric),1) AS pct_rrtinit,              -- [ref] ~5.6
       SUM(prog_event) AS n_prog, SUM(recovery_event) AS n_recov,         -- [ref] ~4,987 / ~4,101
       ROUND(100.0*AVG(nco_vte::numeric),2) AS pct_vte,                   -- [ref] ~6.44
       ROUND(100.0*AVG(nco_pressure::numeric),2) AS pct_pi                -- [ref] ~4.19
FROM cohort_outcome;


-- ############################################################################
-- PART C — COVARIATE EXTRACTION (for the propensity-score model)
--   Base table: cohort_exposure (N = 21,883), keyed on stay_id.
--   All covariates measured "at or before T0" (SAP section 9).
--
--   Notes:
--     * cov_baseline carries optional SAPS-II / OASIS joins, commented out by
--       default: those are first-day scores that can include post-T0 data for
--       early-T0 patients, so they are for Table 1 description only and never
--       enter the PS model.
--     * cov_labs caps extreme values (INR>10->10, bilirubin>30->30, WBC>100->100)
--       AFTER removing physiologically impossible values, so true missingness is
--       preserved for multiple imputation.
--     * SOFA-at-T0 is the single severity score in the PS model (cleanly pre-T0).
--     * KDIGO stage at T0, age, sex, baseline_egfr and exposure already live in
--       cohort_exposure and are pulled in when the analytic table is assembled;
--       they are not re-extracted here.
-- ############################################################################

-- ============================================================================
-- C1. cov_baseline: demographics + comorbidities (key = stay_id)
-- ============================================================================
DROP TABLE IF EXISTS cov_baseline;
CREATE TABLE cov_baseline AS
WITH htn AS (   -- hypertension (not in Charlson; coded separately. ICD codes are space-padded, so LIKE 'X%' is safe)
    SELECT DISTINCT ce.stay_id
    FROM cohort_exposure ce
    INNER JOIN mimiciv_hosp.diagnoses_icd d ON d.hadm_id = ce.hadm_id
    WHERE (d.icd_version = 9  AND (d.icd_code LIKE '401%' OR d.icd_code LIKE '402%'
            OR d.icd_code LIKE '403%' OR d.icd_code LIKE '404%' OR d.icd_code LIKE '405%'))
       OR (d.icd_version = 10 AND (d.icd_code LIKE 'I10%' OR d.icd_code LIKE 'I11%'
            OR d.icd_code LIKE 'I12%' OR d.icd_code LIKE 'I13%'
            OR d.icd_code LIKE 'I15%' OR d.icd_code LIKE 'I16%'))
)
SELECT
    ce.stay_id, ce.hadm_id, ce.subject_id,
    -- demographics / context
    adm.race,
    adm.admission_type,                        -- admission type
    adm.admission_location,                    -- admission source
    pat.anchor_year_group,                     -- admission era (time confounder: ranitidine 2019-20 recall)
    icu.first_careunit,                        -- ICU type (proxy for medical/surgical)
    wt.weight AS weight_kg,
    -- comorbidities (Charlson 0/1 flags taken as-is, not hand-coded)
    ch.myocardial_infarct, ch.congestive_heart_failure, ch.peripheral_vascular_disease,
    ch.cerebrovascular_disease, ch.dementia, ch.chronic_pulmonary_disease, ch.rheumatic_disease,
    ch.peptic_ulcer_disease, ch.mild_liver_disease, ch.severe_liver_disease,
    ch.diabetes_without_cc, ch.diabetes_with_cc, ch.paraplegia, ch.renal_disease,
    ch.malignant_cancer, ch.metastatic_solid_tumor, ch.aids,
    ch.charlson_comorbidity_index,
    -- hypertension (coded separately above)
    CASE WHEN h.stay_id IS NOT NULL THEN 1 ELSE 0 END AS hypertension,
    -- Elixhauser summary score (your table, joined on hadm_id)
    elix.elixhauser_score
    -- Optional SAPS-II / OASIS (first day; Table 1 only, NOT in the PS model).
    -- If your derived schema lacks sapsii/oasis, leave commented; to enable,
    -- build them with mimic-code first, then uncomment the 2 columns + 2 joins.
    -- , sps.sapsii AS sapsii_firstday
    -- , oas.oasis  AS oasis_firstday
FROM cohort_exposure ce
LEFT JOIN mimiciv_hosp.admissions          adm  ON adm.hadm_id    = ce.hadm_id
LEFT JOIN mimiciv_hosp.patients            pat  ON pat.subject_id = ce.subject_id
LEFT JOIN mimiciv_icu.icustays             icu  ON icu.stay_id    = ce.stay_id
LEFT JOIN mimiciv_derived.first_day_weight wt   ON wt.stay_id     = ce.stay_id
LEFT JOIN mimiciv_derived.charlson         ch   ON ch.hadm_id     = ce.hadm_id
LEFT JOIN public.elixhauser_scores         elix ON elix.hadm_id   = ce.hadm_id
-- LEFT JOIN mimiciv_derived.sapsii         sps  ON sps.stay_id    = ce.stay_id   -- uncomment to enable
-- LEFT JOIN mimiciv_derived.oasis          oas  ON oas.stay_id    = ce.stay_id   -- uncomment to enable
LEFT JOIN htn h ON h.stay_id = ce.stay_id;
SELECT COUNT(*) AS n_rows, COUNT(DISTINCT stay_id) AS n_stays FROM cov_baseline;   -- [ref] 21,883
SELECT
    ROUND(AVG(weight_kg)::numeric,1)               AS mean_weight,        -- [ref] ~84.9
    SUM((weight_kg IS NULL)::int)                  AS n_weight_missing,   -- [ref] ~59
    SUM(congestive_heart_failure)                  AS n_chf,              -- [ref] ~33%
    SUM(diabetes_without_cc)+SUM(diabetes_with_cc) AS n_diabetes,        -- [ref] ~36%
    SUM(hypertension)                              AS n_htn,             -- [ref] ~66%
    SUM((mild_liver_disease=1 OR severe_liver_disease=1)::int) AS n_liver, -- [ref] ~17%
    SUM(renal_disease)                             AS n_ckd_charlson,    -- [ref] ~21%
    SUM(malignant_cancer)                          AS n_cancer,          -- [ref] ~14%
    ROUND(AVG(charlson_comorbidity_index)::numeric,1) AS mean_charlson,  -- [ref] ~5.4
    SUM((charlson_comorbidity_index IS NULL)::int)    AS n_charlson_miss, -- [ref] 0
    ROUND(AVG(elixhauser_score)::numeric,1)        AS mean_elix,         -- [ref] ~15.0
    SUM((elixhauser_score IS NULL)::int)           AS n_elix_miss        -- [ref] ~5
FROM cov_baseline;
-- category distributions (confirm populated; see which levels to collapse for modeling)
SELECT first_careunit,      COUNT(*) AS n FROM cov_baseline GROUP BY first_careunit      ORDER BY n DESC;
SELECT admission_type,      COUNT(*) AS n FROM cov_baseline GROUP BY admission_type      ORDER BY n DESC;
SELECT anchor_year_group,   COUNT(*) AS n FROM cov_baseline GROUP BY anchor_year_group   ORDER BY anchor_year_group;

-- ============================================================================
-- C2. cov_labs: closest-to-T0 labs (key = stay_id)
--     Window [T0-48h, min(first acid-suppressant time, T0+24h)]; per patient per
--     analyte take the value nearest T0 (ties prefer pre-T0). Impossible values
--     removed first; then bilirubin/INR/WBC capped while KEEPING NULLs so true
--     missingness is imputable. (Uses CASE caps, not LEAST, since LEAST ignores NULL.)
-- ============================================================================
DROP TABLE IF EXISTS cov_labs;
CREATE TABLE cov_labs AS
WITH base AS (
    SELECT ce.stay_id, ce.subject_id, ce.aki_time AS t0,
           COALESCE(ce.first_acid_time, ce.aki_time + INTERVAL '24 hours') AS cap
    FROM cohort_exposure ce
),
labs AS (
    SELECT b.stay_id, b.t0, le.charttime, le.valuenum,
      CASE le.itemid
        WHEN 50813 THEN 'lactate' WHEN 50862 THEN 'albumin' WHEN 50885 THEN 'bili'
        WHEN 51237 THEN 'inr'     WHEN 51265 THEN 'plt'
        WHEN 51301 THEN 'wbc'     WHEN 51300 THEN 'wbc'
        WHEN 50983 THEN 'na'      WHEN 50824 THEN 'na'
        WHEN 50971 THEN 'k'       WHEN 50822 THEN 'k'
        WHEN 50882 THEN 'hco3'    WHEN 50868 THEN 'ag'  WHEN 50960 THEN 'mg'
        WHEN 51222 THEN 'hgb'     WHEN 50811 THEN 'hgb'
      END AS lab
    FROM base b
    INNER JOIN mimiciv_hosp.labevents le
      ON le.subject_id = b.subject_id
     AND le.itemid IN (50813,50862,50885,51237,51265,51301,51300,
                       50983,50824,50971,50822,50882,50868,50960,51222,50811)
     AND le.valuenum IS NOT NULL
     AND le.charttime >= b.t0 - INTERVAL '48 hours'
     AND le.charttime <= b.cap
),
labs_clean AS (   -- drop physiologically impossible values
    SELECT stay_id, t0, charttime, lab, valuenum
    FROM labs
    WHERE (lab='lactate' AND valuenum BETWEEN 0   AND 50)
       OR (lab='albumin' AND valuenum BETWEEN 0.5 AND 7)
       OR (lab='bili'    AND valuenum BETWEEN 0   AND 100)
       OR (lab='inr'     AND valuenum BETWEEN 0   AND 30)
       OR (lab='plt'     AND valuenum BETWEEN 0   AND 2000)
       OR (lab='wbc'     AND valuenum BETWEEN 0   AND 300)
       OR (lab='na'      AND valuenum BETWEEN 100 AND 180)
       OR (lab='k'       AND valuenum BETWEEN 1   AND 10)
       OR (lab='hco3'    AND valuenum BETWEEN 2   AND 60)
       OR (lab='ag'      AND valuenum BETWEEN 0   AND 60)
       OR (lab='mg'      AND valuenum BETWEEN 0   AND 10)
       OR (lab='hgb'     AND valuenum BETWEEN 2   AND 25)
),
pick AS (   -- per patient per analyte: value nearest T0 (ties prefer pre-T0)
    SELECT DISTINCT ON (stay_id, lab) stay_id, lab, valuenum
    FROM labs_clean
    ORDER BY stay_id, lab,
             ABS(EXTRACT(EPOCH FROM (charttime - t0))) ASC,
             (charttime <= t0) DESC
),
agg AS (   -- aggregate raw values (uncapped); NULL stays NULL
    SELECT
        b.stay_id,
        MAX(CASE WHEN lab='lactate' THEN valuenum END) AS lactate,
        MAX(CASE WHEN lab='albumin' THEN valuenum END) AS albumin,
        MAX(CASE WHEN lab='bili'    THEN valuenum END) AS bili_raw,
        MAX(CASE WHEN lab='inr'     THEN valuenum END) AS inr_raw,
        MAX(CASE WHEN lab='plt'     THEN valuenum END) AS platelet,
        MAX(CASE WHEN lab='wbc'     THEN valuenum END) AS wbc_raw,
        MAX(CASE WHEN lab='na'      THEN valuenum END) AS sodium,
        MAX(CASE WHEN lab='k'       THEN valuenum END) AS potassium,
        MAX(CASE WHEN lab='hco3'    THEN valuenum END) AS bicarbonate,
        MAX(CASE WHEN lab='ag'      THEN valuenum END) AS anion_gap,
        MAX(CASE WHEN lab='mg'      THEN valuenum END) AS magnesium,
        MAX(CASE WHEN lab='hgb'     THEN valuenum END) AS hemoglobin
    FROM base b
    LEFT JOIN pick p ON p.stay_id = b.stay_id
    GROUP BY b.stay_id
)
SELECT
    stay_id, lactate, albumin,
    CASE WHEN bili_raw > 30  THEN 30  ELSE bili_raw END AS bilirubin_total,  -- cap, keep NULL
    CASE WHEN inr_raw  > 10  THEN 10  ELSE inr_raw  END AS inr,              -- cap, keep NULL
    platelet,
    CASE WHEN wbc_raw  > 100 THEN 100 ELSE wbc_raw  END AS wbc,              -- cap, keep NULL
    sodium, potassium, bicarbonate, anion_gap, magnesium, hemoglobin
FROM agg;
SELECT COUNT(*) AS n_rows FROM cov_labs;   -- [ref] 21,883
SELECT
    ROUND(PERCENTILE_CONT(0.5)  WITHIN GROUP (ORDER BY bilirubin_total)::numeric,2) AS bili_med,  -- [ref] ~0.8-1.2
    ROUND(PERCENTILE_CONT(0.75) WITHIN GROUP (ORDER BY bilirubin_total)::numeric,2) AS bili_p75,  -- [ref] ~1.5-3
    ROUND(100.0*AVG((bilirubin_total IS NULL)::int),1) AS miss_bili,
    ROUND(100.0*AVG((inr IS NULL)::int),1)             AS miss_inr,
    ROUND(100.0*AVG((wbc IS NULL)::int),1)             AS miss_wbc,
    MAX(bilirubin_total) AS bili_max,  -- <= 30
    MAX(inr) AS inr_max,               -- <= 10
    MAX(wbc) AS wbc_max                -- <= 100
FROM cov_labs;

-- ============================================================================
-- C3. cov_support: organ support + fluids/urine output (key = stay_id)
--     Fluid/UO window [max(intime, T0-24h), T0]; support flags started before T0.
--     Negative fluid balance is expected: MIMIC inputevents are ICU-only and miss
--     ED resuscitation, so fluid balance is a relative covariate (note in paper).
-- ============================================================================
DROP TABLE IF EXISTS cov_support;
CREATE TABLE cov_support AS
WITH base AS (
    SELECT ce.stay_id, ce.hadm_id, ce.subject_id, ce.intime, ce.aki_time AS t0,
           GREATEST(ce.intime, ce.aki_time - INTERVAL '24 hours') AS win_start,
           cb.weight_kg
    FROM cohort_exposure ce
    LEFT JOIN cov_baseline cb ON cb.stay_id = ce.stay_id
),
vent AS (   -- invasive mechanical ventilation started before T0
    SELECT DISTINCT b.stay_id
    FROM base b INNER JOIN mimiciv_derived.ventilation v ON v.stay_id = b.stay_id
    WHERE v.ventilation_status IN ('InvasiveVent','Trach') AND v.starttime <= b.t0
),
vaso AS (   -- vasopressors before T0
    SELECT DISTINCT b.stay_id
    FROM base b INNER JOIN mimiciv_derived.vasoactive_agent va ON va.stay_id = b.stay_id
    WHERE va.starttime <= b.t0
      AND (va.norepinephrine>0 OR va.epinephrine>0 OR va.dopamine>0
        OR va.phenylephrine>0 OR va.vasopressin>0)
),
steroid AS (   -- systemic corticosteroids before T0
    SELECT DISTINCT b.stay_id
    FROM base b INNER JOIN mimiciv_hosp.prescriptions rx ON rx.hadm_id = b.hadm_id
    WHERE rx.starttime < b.t0
      AND (rx.drug ILIKE '%hydrocortisone%' OR rx.drug ILIKE '%methylprednisolone%'
        OR rx.drug ILIKE '%prednisone%' OR rx.drug ILIKE '%prednisolone%'
        OR rx.drug ILIKE '%dexamethasone%')
),
fluid_in AS (
    SELECT b.stay_id, SUM(ie.amount) AS in_ml
    FROM base b INNER JOIN mimiciv_icu.inputevents ie ON ie.stay_id = b.stay_id
    WHERE ie.amountuom = 'ml' AND ie.starttime >= b.win_start AND ie.starttime <= b.t0
    GROUP BY b.stay_id
),
fluid_out AS (
    SELECT b.stay_id, SUM(oe.value) AS out_ml
    FROM base b INNER JOIN mimiciv_icu.outputevents oe ON oe.stay_id = b.stay_id
    WHERE oe.charttime >= b.win_start AND oe.charttime <= b.t0
    GROUP BY b.stay_id
),
uo AS (
    SELECT b.stay_id, SUM(u.urineoutput) AS uo_ml
    FROM base b INNER JOIN mimiciv_derived.urine_output u ON u.stay_id = b.stay_id
    WHERE u.charttime >= b.win_start AND u.charttime <= b.t0
    GROUP BY b.stay_id
)
SELECT
    b.stay_id,
    CASE WHEN v.stay_id  IS NOT NULL THEN 1 ELSE 0 END AS mech_vent,
    CASE WHEN va.stay_id IS NOT NULL THEN 1 ELSE 0 END AS vasopressor,
    CASE WHEN s.stay_id  IS NOT NULL THEN 1 ELSE 0 END AS corticosteroid,
    CASE WHEN fi.in_ml IS NULL AND fo.out_ml IS NULL THEN NULL
         ELSE COALESCE(fi.in_ml,0) - COALESCE(fo.out_ml,0) END AS fluid_balance_ml,
    CASE WHEN b.weight_kg IS NOT NULL
              AND EXTRACT(EPOCH FROM (b.t0 - b.win_start))/3600.0 > 0
         THEN ROUND((u.uo_ml / b.weight_kg
              / (EXTRACT(EPOCH FROM (b.t0 - b.win_start))/3600.0))::numeric, 2) END AS uo_rate_ml_kg_h
FROM base b
LEFT JOIN vent v    ON v.stay_id  = b.stay_id
LEFT JOIN vaso va   ON va.stay_id = b.stay_id
LEFT JOIN steroid s ON s.stay_id  = b.stay_id
LEFT JOIN fluid_in fi  ON fi.stay_id = b.stay_id
LEFT JOIN fluid_out fo ON fo.stay_id = b.stay_id
LEFT JOIN uo u      ON u.stay_id  = b.stay_id;
SELECT COUNT(*) AS n_rows FROM cov_support;
SELECT
    ROUND(100.0*AVG(mech_vent::numeric),1)      AS pct_vent,      -- [ref] ~53.7%
    ROUND(100.0*AVG(vasopressor::numeric),1)    AS pct_vaso,      -- [ref] ~40%
    ROUND(100.0*AVG(corticosteroid::numeric),1) AS pct_steroid,   -- [ref] ~16.1%
    PERCENTILE_CONT(0.5) WITHIN GROUP (ORDER BY fluid_balance_ml) AS med_fluid_bal,  -- [ref] ~-1220
    ROUND(100.0*AVG((fluid_balance_ml IS NULL)::int),1) AS miss_fluid,               -- [ref] ~7.1%
    PERCENTILE_CONT(0.5) WITHIN GROUP (ORDER BY uo_rate_ml_kg_h) AS med_uo_rate,     -- [ref] ~0.69
    ROUND(100.0*AVG((uo_rate_ml_kg_h IS NULL)::int),1) AS miss_uo                    -- [ref] ~8.1%
FROM cov_support;

-- ============================================================================
-- C4. cov_sofa: SOFA total + 6 subscores (key = stay_id)
--     The only severity score in the PS model. Take the *_24hours columns from
--     the most recent SOFA row at or before T0 = "worst SOFA in 24 h before T0".
-- ============================================================================
DROP TABLE IF EXISTS cov_sofa;
CREATE TABLE cov_sofa AS
WITH base AS (
    SELECT ce.stay_id, ce.aki_time AS t0 FROM cohort_exposure ce
),
pick AS (   -- the latest hourly row at or before T0
    SELECT DISTINCT ON (b.stay_id)
        b.stay_id,
        s.sofa_24hours,
        s.respiration_24hours, s.coagulation_24hours, s.liver_24hours,
        s.cardiovascular_24hours, s.cns_24hours, s.renal_24hours
    FROM base b
    INNER JOIN mimiciv_derived.sofa s ON s.stay_id = b.stay_id
    WHERE s.starttime <= b.t0
    ORDER BY b.stay_id, s.starttime DESC
)
SELECT
    b.stay_id,
    p.sofa_24hours           AS sofa_total,
    p.respiration_24hours    AS sofa_resp,
    p.coagulation_24hours    AS sofa_coag,
    p.liver_24hours          AS sofa_liver,
    p.cardiovascular_24hours AS sofa_cardio,
    p.cns_24hours            AS sofa_cns,
    p.renal_24hours          AS sofa_renal
FROM base b
LEFT JOIN pick p ON p.stay_id = b.stay_id;
SELECT COUNT(*) AS n_rows FROM cov_sofa;
SELECT
    PERCENTILE_CONT(0.25) WITHIN GROUP (ORDER BY sofa_total) AS q1_sofa,   -- [ref] ~3
    PERCENTILE_CONT(0.5)  WITHIN GROUP (ORDER BY sofa_total) AS med_sofa,  -- [ref] ~4
    PERCENTILE_CONT(0.75) WITHIN GROUP (ORDER BY sofa_total) AS q3_sofa,   -- [ref] ~7
    ROUND(AVG(sofa_total)::numeric,1)             AS mean_sofa,            -- [ref] ~5.0
    ROUND(100.0*AVG((sofa_total IS NULL)::int),1) AS miss_sofa,            -- [ref] ~4.1%
    ROUND(AVG(sofa_renal)::numeric,2)  AS mean_renal,                      -- [ref] ~0.70
    ROUND(AVG(sofa_cardio)::numeric,2) AS mean_cardio,                     -- [ref] ~1.51
    ROUND(AVG(sofa_resp)::numeric,2)   AS mean_resp                        -- [ref] ~1.07
FROM cov_sofa;

-- ============================================================================
-- C5. cov_confounders: key confounders (key = stay_id)
-- ============================================================================
DROP TABLE IF EXISTS cov_confounders;
CREATE TABLE cov_confounders AS
WITH base AS (
    SELECT ce.stay_id, ce.hadm_id, ce.subject_id, ce.aki_time AS t0 FROM cohort_exposure ce
),
gibleed AS (   -- GI bleeding (this admission, billing code; no in-hospital timestamp -> imprecise timing)
    SELECT DISTINCT b.stay_id
    FROM base b INNER JOIN mimiciv_hosp.diagnoses_icd d ON d.hadm_id = b.hadm_id
    WHERE (d.icd_version=9 AND (d.icd_code LIKE '578%' OR d.icd_code IN ('4560','45620','53082')
            OR d.icd_code LIKE '5310%' OR d.icd_code LIKE '5312%' OR d.icd_code LIKE '5314%' OR d.icd_code LIKE '5316%'
            OR d.icd_code LIKE '5320%' OR d.icd_code LIKE '5322%' OR d.icd_code LIKE '5324%' OR d.icd_code LIKE '5326%'
            OR d.icd_code LIKE '5330%' OR d.icd_code LIKE '5332%' OR d.icd_code LIKE '5334%' OR d.icd_code LIKE '5336%'
            OR d.icd_code LIKE '5340%' OR d.icd_code LIKE '5342%' OR d.icd_code LIKE '5344%' OR d.icd_code LIKE '5346%'))
       OR (d.icd_version=10 AND (d.icd_code LIKE 'K920%' OR d.icd_code LIKE 'K921%' OR d.icd_code LIKE 'K922%'
            OR d.icd_code IN ('I8501','I8511','K2211')
            OR d.icd_code LIKE 'K250%' OR d.icd_code LIKE 'K252%' OR d.icd_code LIKE 'K254%' OR d.icd_code LIKE 'K256%'
            OR d.icd_code LIKE 'K260%' OR d.icd_code LIKE 'K262%' OR d.icd_code LIKE 'K264%' OR d.icd_code LIKE 'K266%'
            OR d.icd_code LIKE 'K270%' OR d.icd_code LIKE 'K272%' OR d.icd_code LIKE 'K274%' OR d.icd_code LIKE 'K276%'
            OR d.icd_code LIKE 'K280%' OR d.icd_code LIKE 'K282%' OR d.icd_code LIKE 'K284%' OR d.icd_code LIKE 'K286%'))
),
rbc AS (   -- pre-T0 red-cell transfusion (timestamped generic bleeding/anemia signal, correctly pre-T0)
    SELECT DISTINCT b.stay_id
    FROM base b INNER JOIN mimiciv_icu.inputevents ie
      ON ie.stay_id = b.stay_id AND ie.starttime < b.t0
    WHERE ie.itemid IN (SELECT itemid FROM mimiciv_icu.d_items
                        WHERE label ILIKE '%packed red%' OR label ILIKE '%prbc%'
                           OR label ILIKE '%red blood cell%')
),
dnr AS (   -- DNR / comfort care (code status not full code before T0)
    SELECT DISTINCT b.stay_id
    FROM base b INNER JOIN mimiciv_icu.chartevents c ON c.stay_id = b.stay_id
    WHERE c.itemid = 223758 AND c.charttime <= b.t0
      AND (c.value ILIKE '%DNR%' OR c.value ILIKE '%DNI%' OR c.value ILIKE '%comfort%')
),
nephro AS (   -- nephrotoxins before T0
    SELECT b.stay_id,
      MAX((rx.drug ILIKE '%vancomycin%' AND rx.route ILIKE '%IV%')::int) AS vanco,
      MAX((rx.drug ILIKE '%gentamicin%' OR rx.drug ILIKE '%tobramycin%' OR rx.drug ILIKE '%amikacin%')::int) AS aminoglyc,
      MAX((rx.drug ILIKE '%ibuprofen%' OR rx.drug ILIKE '%ketorolac%' OR rx.drug ILIKE '%naproxen%'
        OR rx.drug ILIKE '%indomethacin%' OR rx.drug ILIKE '%diclofenac%')::int) AS nsaid,
      MAX((rx.drug ILIKE '%lisinopril%' OR rx.drug ILIKE '%enalapril%' OR rx.drug ILIKE '%captopril%'
        OR rx.drug ILIKE '%ramipril%' OR rx.drug ILIKE '%benazepril%' OR rx.drug ILIKE '%losartan%'
        OR rx.drug ILIKE '%valsartan%' OR rx.drug ILIKE '%candesartan%' OR rx.drug ILIKE '%irbesartan%'
        OR rx.drug ILIKE '%olmesartan%' OR rx.drug ILIKE '%telmisartan%')::int) AS raasi,
      MAX((rx.drug ILIKE '%amphotericin%')::int) AS ampho,
      -- IV iodinated contrast (generic + common brand names; under-recorded in MIMIC = lower bound, note in paper)
      MAX((rx.drug ILIKE '%iohexol%'    OR rx.drug ILIKE '%iopamidol%'  OR rx.drug ILIKE '%iodixanol%'
        OR rx.drug ILIKE '%ioversol%'   OR rx.drug ILIKE '%iopromide%'  OR rx.drug ILIKE '%diatrizoate%'
        OR rx.drug ILIKE '%ioxaglate%'  OR rx.drug ILIKE '%omnipaque%'  OR rx.drug ILIKE '%visipaque%'
        OR rx.drug ILIKE '%isovue%'     OR rx.drug ILIKE '%optiray%'    OR rx.drug ILIKE '%ultravist%')::int) AS iv_contrast
    FROM base b INNER JOIN mimiciv_hosp.prescriptions rx ON rx.hadm_id = b.hadm_id
    WHERE rx.starttime < b.t0
    GROUP BY b.stay_id
),
en AS (   -- enteral nutrition before T0
    SELECT DISTINCT b.stay_id
    FROM base b INNER JOIN mimiciv_icu.inputevents ie ON ie.stay_id = b.stay_id
    WHERE ie.starttime < b.t0
      AND ie.itemid IN (SELECT itemid FROM mimiciv_icu.d_items WHERE category ILIKE '%enteral%')
),
antiplt AS (   -- antiplatelets (P2Y12 inhibitors) before T0 (PPI-clopidogrel interaction confounder)
    SELECT DISTINCT b.stay_id
    FROM base b INNER JOIN mimiciv_hosp.prescriptions rx ON rx.hadm_id = b.hadm_id
    WHERE rx.starttime < b.t0
      AND (rx.drug ILIKE '%clopidogrel%' OR rx.drug ILIKE '%prasugrel%' OR rx.drug ILIKE '%ticagrelor%')
)
SELECT
    b.stay_id,
    CASE WHEN g.stay_id  IS NOT NULL THEN 1 ELSE 0 END AS gi_bleed_dx,            -- billing code, imprecise timing
    CASE WHEN r.stay_id  IS NOT NULL THEN 1 ELSE 0 END AS rbc_transfusion_pre_t0, -- timestamped, pre-T0
    CASE WHEN dn.stay_id IS NOT NULL THEN 1 ELSE 0 END AS dnr_comfort,
    COALESCE(n.vanco,0)       AS nephro_vanco,
    COALESCE(n.aminoglyc,0)   AS nephro_aminoglyc,
    COALESCE(n.nsaid,0)       AS nephro_nsaid,
    COALESCE(n.raasi,0)       AS nephro_acei_arb,
    COALESCE(n.ampho,0)       AS nephro_ampho,
    COALESCE(n.iv_contrast,0) AS nephro_iv_contrast,
    CASE WHEN COALESCE(n.vanco,0)+COALESCE(n.aminoglyc,0)+COALESCE(n.nsaid,0)
              +COALESCE(n.raasi,0)+COALESCE(n.ampho,0)+COALESCE(n.iv_contrast,0) > 0
         THEN 1 ELSE 0 END AS nephrotoxin_any,
    CASE WHEN e.stay_id  IS NOT NULL THEN 1 ELSE 0 END AS enteral_nutrition,
    CASE WHEN ap.stay_id IS NOT NULL THEN 1 ELSE 0 END AS antiplatelet
FROM base b
LEFT JOIN gibleed g  ON g.stay_id  = b.stay_id
LEFT JOIN rbc r      ON r.stay_id  = b.stay_id
LEFT JOIN dnr dn     ON dn.stay_id = b.stay_id
LEFT JOIN nephro n   ON n.stay_id  = b.stay_id
LEFT JOIN en e       ON e.stay_id  = b.stay_id
LEFT JOIN antiplt ap ON ap.stay_id = b.stay_id;
SELECT COUNT(*) AS n_rows FROM cov_confounders;
SELECT
    ROUND(100.0*AVG(gi_bleed_dx::numeric),1)            AS pct_gibleed_dx,   -- [ref] ~8.7%
    ROUND(100.0*AVG(rbc_transfusion_pre_t0::numeric),1) AS pct_rbc_pre_t0,   -- [ref] ~12.7%
    ROUND(100.0*AVG(dnr_comfort::numeric),1)            AS pct_dnr,          -- [ref] ~3.2%
    ROUND(100.0*AVG(nephrotoxin_any::numeric),1)        AS pct_nephro_any,   -- [ref] ~51%
    ROUND(100.0*AVG(nephro_vanco::numeric),1)           AS pct_vanco,        -- [ref] ~43.1%
    ROUND(100.0*AVG(nephro_iv_contrast::numeric),1)     AS pct_contrast,     -- [ref] low (<5%, under-recorded)
    ROUND(100.0*AVG(nephro_nsaid::numeric),1)           AS pct_nsaid,        -- [ref] ~6.0%
    ROUND(100.0*AVG(nephro_acei_arb::numeric),1)        AS pct_acei_arb,     -- [ref] ~8.5%
    ROUND(100.0*AVG(enteral_nutrition::numeric),1)      AS pct_en,           -- [ref] ~6.1%
    ROUND(100.0*AVG(antiplatelet::numeric),1)           AS pct_antiplt       -- [ref] ~5.7%
FROM cov_confounders;


-- ############################################################################
-- PART D — ANALYTIC TABLE
--   analytic_cohort = cohort_outcome (anchor) LEFT JOIN the 6 covariate tables,
--   key = stay_id. o.* carries all cohort_outcome columns (IDs + outcomes); each
--   other table contributes only its non-ID columns, so there are no duplicates.
--   Final N ~ 21,731 (= cohort_outcome rows), 121 columns.
--   The table keeps all exposure groups (PPI / H2RA / neither); the head-to-head
--   analysis subsets to PPI vs H2RA in R, with the neither group retained for CCW.
--   ** Patient-level table -> stays local, goes into R; only aggregate output
--      may leave the machine (Data Use Agreement). **
-- ############################################################################
DROP TABLE IF EXISTS analytic_cohort;
CREATE TABLE analytic_cohort AS
SELECT
    o.*,                                                   -- cohort_outcome (IDs + outcomes)
    -- ---- from cohort_exposure (excl. subject_id/hadm_id/stay_id) ----
    e.intime, e.outtime, e.icu_los_days, e.age_at_icu, e.gender,
    e.sepsis_time, e.aki_time, e.first_aki_stage,
    e.creat_pre, e.creat_24h, e.creat_mdrd,
    e.baseline_creat, e.baseline_creat_source, e.baseline_egfr, e.ckd_dx, e.ckd,
    e.exposure_group, e.first_acid_drug, e.first_acid_time, e.hours_to_exposure,
    e.both_in_window, e.exact_tie, e.prior_acid, e.newuser,
    -- ---- from cov_baseline (excl. stay_id/hadm_id/subject_id) ----
    b.race, b.admission_type, b.admission_location, b.anchor_year_group,
    b.first_careunit, b.weight_kg,
    b.myocardial_infarct, b.congestive_heart_failure, b.peripheral_vascular_disease,
    b.cerebrovascular_disease, b.dementia, b.chronic_pulmonary_disease, b.rheumatic_disease,
    b.peptic_ulcer_disease, b.mild_liver_disease, b.severe_liver_disease,
    b.diabetes_without_cc, b.diabetes_with_cc, b.paraplegia, b.renal_disease,
    b.malignant_cancer, b.metastatic_solid_tumor, b.aids,
    b.charlson_comorbidity_index, b.hypertension, b.elixhauser_score,
    -- ---- from cov_labs (excl. stay_id) ----
    l.lactate, l.albumin, l.bilirubin_total, l.inr, l.platelet, l.wbc,
    l.sodium, l.potassium, l.bicarbonate, l.anion_gap, l.magnesium, l.hemoglobin,
    -- ---- from cov_support (excl. stay_id) ----
    s.mech_vent, s.vasopressor, s.corticosteroid, s.fluid_balance_ml, s.uo_rate_ml_kg_h,
    -- ---- from cov_sofa (excl. stay_id) ----
    f.sofa_total, f.sofa_resp, f.sofa_coag, f.sofa_liver, f.sofa_cardio, f.sofa_cns, f.sofa_renal,
    -- ---- from cov_confounders (excl. stay_id) ----
    c.gi_bleed_dx, c.rbc_transfusion_pre_t0, c.dnr_comfort,
    c.nephro_vanco, c.nephro_aminoglyc, c.nephro_nsaid, c.nephro_acei_arb,
    c.nephro_ampho, c.nephro_iv_contrast, c.nephrotoxin_any,
    c.enteral_nutrition, c.antiplatelet
FROM cohort_outcome o
LEFT JOIN cohort_exposure  e ON e.stay_id = o.stay_id
LEFT JOIN cov_baseline     b ON b.stay_id = o.stay_id
LEFT JOIN cov_labs         l ON l.stay_id = o.stay_id
LEFT JOIN cov_support      s ON s.stay_id = o.stay_id
LEFT JOIN cov_sofa         f ON f.stay_id = o.stay_id
LEFT JOIN cov_confounders  c ON c.stay_id = o.stay_id;
-- ---- checks ----
-- (1) row + column count (expect rows ~ 21,731, columns = 121)
SELECT
    (SELECT COUNT(*) FROM analytic_cohort)                        AS n_rows,
    (SELECT COUNT(*) FROM cohort_outcome)                         AS n_rows_outcome_should_match,
    (SELECT COUNT(*) FROM information_schema.columns
      WHERE table_schema='public' AND table_name='analytic_cohort') AS n_cols;
-- (2) join integrity: all should be 0 (every row matched each covariate table)
SELECT
    SUM((exposure_group IS NULL)::int) AS miss_join_exposure,   -- expect 0 (exposure_group is never NULL in source)
    SUM((weight_kg      IS NULL)::int) AS miss_join_baseline,   -- ~59: weight's own missingness, acceptable
    SUM((sofa_total     IS NULL)::int) AS sofa_null,            -- ~4%: SOFA's own missingness (-> MICE)
    SUM((lactate        IS NULL)::int) AS lactate_null          -- ~11%: lactate's own missingness (-> MICE)
FROM analytic_cohort;
-- (3) exposure-group distribution (head-to-head analysis uses PPI + H2RA)
SELECT exposure_group, COUNT(*) AS n,
       ROUND(100.0*COUNT(*)/SUM(COUNT(*)) OVER (),1) AS pct
FROM analytic_cohort
GROUP BY exposure_group ORDER BY n DESC;
-- [ref] neither ~15,200 / PPI ~3,430 / H2RA ~3,050 (slightly fewer than the 21,883 pre-clean)
-- (4) new-initiator subset size (new-user sensitivity uses newuser=1)
SELECT newuser, COUNT(*) AS n FROM analytic_cohort GROUP BY newuser ORDER BY newuser;


-- ############################################################################
-- PART E — SENSITIVITY: 48 h EXPOSURE WINDOW
--   Re-classify exposure over [T0, T0+48h] instead of [T0, T0+24h], then rebuild
--   the analytic table. Only the exposure-window bound changes; the new-user /
--   prior-acid definitions (which look BEFORE T0) are unchanged.
-- ############################################################################

-- ============================================================================
-- E1. cohort_exposure_48h: 48 h version of A5 (window bound only; others identical)
-- ============================================================================
DROP TABLE IF EXISTS cohort_exposure_48h;
CREATE TABLE cohort_exposure_48h AS
WITH acid AS (
    SELECT e.stay_id, e.aki_time, rx.starttime, rx.drug,
      CASE WHEN rx.drug ILIKE '%pantoprazole%' OR rx.drug ILIKE '%omeprazole%'
             OR rx.drug ILIKE '%esomeprazole%' OR rx.drug ILIKE '%lansoprazole%'
             OR rx.drug ILIKE '%rabeprazole%'  OR rx.drug ILIKE '%dexlansoprazole%'
           THEN 'PPI' ELSE 'H2RA' END AS acid_class
    FROM cohort_baseline_egfr e
    JOIN mimiciv_hosp.prescriptions rx
      ON rx.hadm_id=e.hadm_id AND rx.starttime IS NOT NULL
     AND (rx.drug ILIKE '%pantoprazole%' OR rx.drug ILIKE '%omeprazole%'
       OR rx.drug ILIKE '%esomeprazole%' OR rx.drug ILIKE '%lansoprazole%'
       OR rx.drug ILIKE '%rabeprazole%'  OR rx.drug ILIKE '%dexlansoprazole%'
       OR rx.drug ILIKE '%famotidine%'   OR rx.drug ILIKE '%ranitidine%'
       OR rx.drug ILIKE '%cimetidine%'   OR rx.drug ILIKE '%nizatidine%')
),
first_in_window AS (
    SELECT DISTINCT ON (stay_id)
        stay_id, acid_class AS first_acid_class, drug AS first_acid_drug, starttime AS first_acid_time
    FROM acid
    WHERE starttime >= aki_time AND starttime <= aki_time + INTERVAL '48 hours'   -- 48h
    ORDER BY stay_id, starttime ASC, acid_class
),
both_flag AS (
    SELECT stay_id, CASE WHEN COUNT(DISTINCT acid_class)>1 THEN 1 ELSE 0 END AS both_in_window
    FROM acid WHERE starttime>=aki_time AND starttime<=aki_time+INTERVAL '48 hours' GROUP BY stay_id  -- 48h
),
exact_tie_flag AS (
    SELECT stay_id, CASE WHEN COUNT(DISTINCT acid_class)>1 THEN 1 ELSE 0 END AS exact_tie
    FROM (SELECT stay_id, acid_class, starttime, MIN(starttime) OVER (PARTITION BY stay_id) AS min_st
          FROM acid WHERE starttime>=aki_time AND starttime<=aki_time+INTERVAL '48 hours') w           -- 48h
    WHERE starttime=min_st GROUP BY stay_id
),
prior_flag AS (
    SELECT DISTINCT stay_id, 1 AS prior_acid FROM acid WHERE starttime < aki_time   -- unchanged (before T0)
)
SELECT e.stay_id,
       COALESCE(fw.first_acid_class,'neither') AS exposure_group_48h,
       fw.first_acid_drug AS first_acid_drug_48h, fw.first_acid_time AS first_acid_time_48h,
       CASE WHEN fw.first_acid_time IS NOT NULL
            THEN ROUND((EXTRACT(EPOCH FROM (fw.first_acid_time-e.aki_time))/3600.0)::numeric,1) END AS hours_to_exposure_48h,
       COALESCE(bf.both_in_window,0) AS both_in_window_48h,
       COALESCE(et.exact_tie,0)      AS exact_tie_48h,
       COALESCE(pf.prior_acid,0)     AS prior_acid_48h,
       CASE WHEN COALESCE(pf.prior_acid,0)=0 THEN 1 ELSE 0 END AS newuser_48h
FROM cohort_baseline_egfr e
LEFT JOIN first_in_window fw ON fw.stay_id=e.stay_id
LEFT JOIN both_flag bf ON bf.stay_id=e.stay_id
LEFT JOIN exact_tie_flag et ON et.stay_id=e.stay_id
LEFT JOIN prior_flag pf ON pf.stay_id=e.stay_id;
SELECT exposure_group_48h, COUNT(*) AS n, SUM(newuser_48h) AS n_newuser,
       SUM(both_in_window_48h) AS n_both, SUM(exact_tie_48h) AS n_tie
FROM cohort_exposure_48h GROUP BY exposure_group_48h ORDER BY n DESC;
-- [ref] neither 13,010 (nu 3,538) | PPI 4,714 (nu 1,813) | H2RA 4,159 (nu 1,597)
-- migration 24h -> 48h: confirms 24h initiators keep their class; only neithers
-- gain a classification when their first acid suppressant falls in 24-48h.
SELECT c24.exposure_group AS grp_24h, c48.exposure_group_48h AS grp_48h, COUNT(*) AS n
FROM cohort_exposure c24 JOIN cohort_exposure_48h c48 USING (stay_id)
GROUP BY c24.exposure_group, c48.exposure_group_48h
ORDER BY grp_24h, n DESC;

-- ============================================================================
-- E2. analytic_cohort_48h: clone of analytic_cohort with exposure columns sourced
--     from cohort_exposure_48h (e2). Baseline/ID/creatinine/eGFR columns still come
--     from cohort_exposure (e) since they are window-independent. 121 columns, so it
--     drops straight into the same R imputation pipeline.
-- ============================================================================
DROP TABLE IF EXISTS analytic_cohort_48h;
CREATE TABLE analytic_cohort_48h AS
SELECT
    o.*,
    -- baseline columns: still from cohort_exposure e (window-independent)
    e.intime, e.outtime, e.icu_los_days, e.age_at_icu, e.gender,
    e.sepsis_time, e.aki_time, e.first_aki_stage,
    e.creat_pre, e.creat_24h, e.creat_mdrd,
    e.baseline_creat, e.baseline_creat_source, e.baseline_egfr, e.ckd_dx, e.ckd,
    -- exposure-classification columns: from cohort_exposure_48h e2, aliased to standard names
    e2.exposure_group_48h    AS exposure_group,
    e2.first_acid_drug_48h   AS first_acid_drug,
    e2.first_acid_time_48h   AS first_acid_time,
    e2.hours_to_exposure_48h AS hours_to_exposure,
    e2.both_in_window_48h    AS both_in_window,
    e2.exact_tie_48h         AS exact_tie,
    e2.prior_acid_48h        AS prior_acid,
    e2.newuser_48h           AS newuser,
    -- cov_baseline
    b.race, b.admission_type, b.admission_location, b.anchor_year_group,
    b.first_careunit, b.weight_kg,
    b.myocardial_infarct, b.congestive_heart_failure, b.peripheral_vascular_disease,
    b.cerebrovascular_disease, b.dementia, b.chronic_pulmonary_disease, b.rheumatic_disease,
    b.peptic_ulcer_disease, b.mild_liver_disease, b.severe_liver_disease,
    b.diabetes_without_cc, b.diabetes_with_cc, b.paraplegia, b.renal_disease,
    b.malignant_cancer, b.metastatic_solid_tumor, b.aids,
    b.charlson_comorbidity_index, b.hypertension, b.elixhauser_score,
    -- cov_labs
    l.lactate, l.albumin, l.bilirubin_total, l.inr, l.platelet, l.wbc,
    l.sodium, l.potassium, l.bicarbonate, l.anion_gap, l.magnesium, l.hemoglobin,
    -- cov_support
    s.mech_vent, s.vasopressor, s.corticosteroid, s.fluid_balance_ml, s.uo_rate_ml_kg_h,
    -- cov_sofa
    f.sofa_total, f.sofa_resp, f.sofa_coag, f.sofa_liver, f.sofa_cardio, f.sofa_cns, f.sofa_renal,
    -- cov_confounders
    c.gi_bleed_dx, c.rbc_transfusion_pre_t0, c.dnr_comfort,
    c.nephro_vanco, c.nephro_aminoglyc, c.nephro_nsaid, c.nephro_acei_arb,
    c.nephro_ampho, c.nephro_iv_contrast, c.nephrotoxin_any,
    c.enteral_nutrition, c.antiplatelet
FROM cohort_outcome o
LEFT JOIN cohort_exposure      e  ON e.stay_id  = o.stay_id   -- baseline / ID columns
LEFT JOIN cohort_exposure_48h  e2 ON e2.stay_id = o.stay_id   -- 48h exposure classification
LEFT JOIN cov_baseline     b ON b.stay_id = o.stay_id
LEFT JOIN cov_labs         l ON l.stay_id = o.stay_id
LEFT JOIN cov_support      s ON s.stay_id = o.stay_id
LEFT JOIN cov_sofa         f ON f.stay_id = o.stay_id
LEFT JOIN cov_confounders  c ON c.stay_id = o.stay_id;
-- checks: column count = 121, exposure_group never NULL, same row count as 24h
SELECT (SELECT COUNT(*) FROM analytic_cohort_48h) AS n_rows,
       (SELECT COUNT(*) FROM information_schema.columns
         WHERE table_schema='public' AND table_name='analytic_cohort_48h') AS n_cols;
SELECT SUM((exposure_group IS NULL)::int) AS miss_exposure,    -- expect 0
       SUM((aki_time       IS NULL)::int) AS miss_akitime      -- expect 0
FROM analytic_cohort_48h;
SELECT exposure_group, COUNT(*) AS n FROM analytic_cohort_48h GROUP BY exposure_group ORDER BY n DESC;

-- ============================================================================
-- Pipeline complete. Export analytic_cohort (and analytic_cohort_48h) to R
-- locally, e.g.:
--   COPY (SELECT * FROM analytic_cohort) TO '/local/path/analytic_cohort.csv'
--        WITH (FORMAT csv, HEADER true);
--   -- or, from R:  DBI::dbReadTable(con, "analytic_cohort")
-- Everything stays local; only aggregate results leave the machine.
-- ============================================================================
