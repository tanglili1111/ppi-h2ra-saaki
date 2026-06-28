-- ############################################################################
-- ############################################################################
--  PPI vs. H2RA in sepsis-associated AKI (SA-AKI) | eICU external validation | full data pipeline
--  ===========================================================================
--  Integrates Stage 1-4; run top to bottom. Aligned with the MIMIC-IV primary cohort.
--  Final output table: eicu_analytic_cohort (cohort + exposure + in-hospital outcomes + propensity-score covariates)
--
--  Stage 1: sepsis -> AKI [creatinine + urine output, KDIGO] -> exclusions (ESRD/transplant/pregnancy/chronic dialysis) -> baseline renal function
--  Stage 2: exposure PPI vs. H2RA (first acid suppressant in [T0, T0+24h]) + new-user determination
--  Stage 3: in-hospital outcomes (MAKE/death/RRT/progression/recovery; follow-up = [T0, discharge]; mirrors MIMIC B1-B6)
--  Stage 4: covariates (APACHE scores/comorbidities/day-1 physiology/ventilation/vasoactives/lactate) + assemble analytic_cohort
--
--  Time unit: offset (minutes since ICU admission). Creatinine in mg/dL. eGFR by CKD-EPI 2021.
--  Prerequisite: the standard eICU tables (patient, lab, diagnosis, treatment, intakeoutput, medication,
--        infusiondrug, admissiondrug, pasthistory, apachepredvar, apacheapsvar,
--        apachepatientresult) are already deployed in the current database.
-- ############################################################################
-- ############################################################################


-- ============================================================================
-- STAGE 0 (optional, strongly recommended): index the large tables to speed up JOINs several-fold.
--   lab has millions of rows; intakeoutput/treatment/diagnosis are also large. Run once before building the tables.
-- ============================================================================
CREATE INDEX IF NOT EXISTS idx_lab_pid   ON lab(patientunitstayid);
CREATE INDEX IF NOT EXISTS idx_io_pid    ON intakeoutput(patientunitstayid);
CREATE INDEX IF NOT EXISTS idx_treat_pid ON treatment(patientunitstayid);
CREATE INDEX IF NOT EXISTS idx_dx_pid    ON diagnosis(patientunitstayid);
CREATE INDEX IF NOT EXISTS idx_med_pid   ON medication(patientunitstayid);
CREATE INDEX IF NOT EXISTS idx_inf_pid   ON infusiondrug(patientunitstayid);
CREATE INDEX IF NOT EXISTS idx_ph_pid    ON pasthistory(patientunitstayid);
CREATE INDEX IF NOT EXISTS idx_admdrug_pid ON admissiondrug(patientunitstayid);



-- ============================================================================
--  eICU extraction pipeline | STAGE 1 (v7): cohort (sepsis -> AKI [creatinine + urine output] -> exclusions -> baseline)
--  ---------------------------------------------------------------------------
--  Aligned with MIMIC (whose kdigo_stages is also creatinine + urine output). Time in offset (minutes since ICU admission).
--  AKI = any one KDIGO criterion:
--    - Creatinine: >= 1.5x nadir baseline, or a rise >= 0.3 within 48h (measured at off >= 0; includes AKI present on admission)
--    - Urine output: 6h rolling < 0.5 mL/kg/h (window span >= 5h; requires admissionweight)
--        v6: urine-output labels keep true volumes only (exclude count labels such as Urine Count/Occurrence/Incontinence) + cap at 5000 mL
--    T0 = the earlier of the two criteria; record the by_creat / by_uo source.
--    [v7] first_aki_stage = creatinine stage at T0 (true baseline severity), no longer the peak (the peak carries post-treatment information,
--         which would corrupt the progression endpoint and introduce collider bias in the propensity score). Aligned with MIMIC's "KDIGO stage at T0".
--  Baseline creatinine = lowest creatinine (nadir) within the measurement window; if none, fall back to MDRD-75.
--  At the end: a urine-output structure-calibration query (recommended to run separately first, to confirm values are per-interval volumes).
-- ============================================================================


-- E1. Sepsis + adult (>=18) + first ICU stay -> eicu_cohort_base
DROP TABLE IF EXISTS eicu_cohort_base;
CREATE TABLE eicu_cohort_base AS
WITH sepsis AS (
  SELECT DISTINCT patientunitstayid FROM diagnosis WHERE LOWER(diagnosisstring) ~ 'sepsis|septic'
)
SELECT
  p.patientunitstayid AS stay_id, p.patienthealthsystemstayid AS hadm_id, p.gender,
  CASE WHEN p.age = '> 89' THEN 90 WHEN p.age ~ '^[0-9]+$' THEN p.age::int END AS age_at_icu,
  p.ethnicity, p.unittype, p.unitadmitsource, p.hospitaladmitsource, p.hospitalid,
  p.admissionweight AS weight_kg,
  p.hospitaldischargeoffset, p.hospitaldischargestatus,
  p.unitdischargeoffset, p.unitdischargestatus, p.unitvisitnumber
FROM patient p
JOIN sepsis s ON s.patientunitstayid = p.patientunitstayid
WHERE p.unitvisitnumber = 1
  AND (p.age = '> 89' OR (p.age ~ '^[0-9]+$' AND p.age::int >= 18))
  AND p.hospitaldischargeoffset IS NOT NULL AND p.hospitaldischargeoffset > 0;
SELECT 'E1_sepsis_adult_firstICU' AS step, COUNT(*) AS n FROM eicu_cohort_base;   -- [known: 22,828]


-- E2. AKI/KDIGO (creatinine + urine output) -> T0 + source + stage + baseline + ckd_hist
DROP TABLE IF EXISTS eicu_cohort_aki;
CREATE TABLE eicu_cohort_aki AS
WITH creat AS (
  SELECT b.stay_id, l.labresultoffset AS off, l.labresult AS creat
  FROM eicu_cohort_base b
  JOIN lab l ON l.patientunitstayid = b.stay_id
   AND LOWER(l.labname) = 'creatinine'
   AND l.labresult IS NOT NULL AND l.labresult > 0.1 AND l.labresult < 50
   AND l.labresultoffset BETWEEN -1440 AND 7*1440
),
nadir AS (
  SELECT stay_id, MIN(creat) AS base_creat FROM creat GROUP BY stay_id
),
ckd_hist AS (
  SELECT b.stay_id,
    MAX(CASE WHEN LOWER(h.pasthistoryvalue) ~ 'chronic kidney|ckd|renal insufficiency|chronic renal|renal failure|dialysis|esrd|end-stage renal'
              AND LOWER(h.pasthistoryvalue) !~ 'kidney stones' THEN 1 ELSE 0 END) AS ckd_hist
  FROM eicu_cohort_base b LEFT JOIN pasthistory h ON h.patientunitstayid = b.stay_id
  GROUP BY b.stay_id
),
staged AS (
  SELECT cr.stay_id, cr.off, cr.creat, n.base_creat,
         cr.creat - MIN(cr.creat) OVER (
           PARTITION BY cr.stay_id ORDER BY cr.off
           RANGE BETWEEN 2880 PRECEDING AND CURRENT ROW) AS rise_48h
  FROM creat cr JOIN nadir n ON n.stay_id = cr.stay_id
),
creat_aki AS (   -- first AKI by the creatinine criterion (off >= 0)
  SELECT stay_id, MIN(off) AS aki_off
  FROM staged
  WHERE (creat >= 1.5*base_creat OR rise_48h >= 0.3) AND off >= 0
  GROUP BY stay_id
),
-- --- Urine-output KDIGO ---
uo AS (
  SELECT b.stay_id, io.intakeoutputoffset AS off, io.cellvaluenumeric AS vol, b.weight_kg
  FROM eicu_cohort_base b
  JOIN intakeoutput io ON io.patientunitstayid = b.stay_id
  WHERE LOWER(io.celllabel) ~ 'urine|foley|void'
    AND LOWER(io.celllabel) !~ 'count|occurrence|incontinen|number of|unmeasured|stool'  -- exclude non-volume labels such as counts/incontinence/stool
    AND io.cellvaluenumeric IS NOT NULL
    AND io.cellvaluenumeric >= 0 AND io.cellvaluenumeric < 5000   -- cap: exclude data-entry errors (e.g., 1.8M mL)
    AND io.intakeoutputoffset BETWEEN 0 AND 7*1440
    AND b.weight_kg IS NOT NULL AND b.weight_kg > 0
),
uo_roll AS (
  SELECT stay_id, off, weight_kg,
    SUM(vol) OVER w AS uo_6h,
    off - MIN(off) OVER w AS span_min
  FROM uo
  WINDOW w AS (PARTITION BY stay_id ORDER BY off RANGE BETWEEN 360 PRECEDING AND CURRENT ROW)
),
uo_aki AS (   -- 6h rolling < 0.5 mL/kg/h (window span >= 5h)
  SELECT stay_id, MIN(off) AS aki_off
  FROM uo_roll
  WHERE span_min >= 300 AND uo_6h < 0.5 * weight_kg * (span_min/60.0)
  GROUP BY stay_id
),
-- --- Combine the two criteria ---
aki_all AS (
  SELECT stay_id, MIN(aki_off) AS aki_offset,
    BOOL_OR(src = 'creat') AS by_creat,
    BOOL_OR(src = 'uo')    AS by_uo
  FROM (SELECT stay_id, aki_off, 'creat' AS src FROM creat_aki
        UNION ALL
        SELECT stay_id, aki_off, 'uo'    AS src FROM uo_aki) z
  GROUP BY stay_id
),
stage_at_t0 AS (   -- [v7] creatinine stage at T0 (true baseline-severity covariate; UO-triggered with no creatinine rise -> 1)
  SELECT a.stay_id,
    MAX(CASE WHEN s.creat >= 3.0*s.base_creat OR s.creat >= 4.0 THEN 3
             WHEN s.creat >= 2.0*s.base_creat THEN 2
             WHEN s.creat >= 1.5*s.base_creat THEN 1 ELSE 0 END) AS creat_stage_t0
  FROM aki_all a JOIN staged s ON s.stay_id = a.stay_id AND s.off = a.aki_offset
  GROUP BY a.stay_id
)
SELECT b.*, a.aki_offset AS t0_offset, a.by_creat, a.by_uo,
       GREATEST(COALESCE(st.creat_stage_t0,0), 1) AS first_aki_stage,
       COALESCE(n.base_creat,
         ROUND((POWER(175.0*POWER(b.age_at_icu,-0.203)*(CASE WHEN b.gender='Female' THEN 0.742 ELSE 1.0 END)/75.0, 1.0/1.154))::numeric,2)
       ) AS baseline_creat,
       COALESCE(ck.ckd_hist,0) AS ckd_hist
FROM eicu_cohort_base b
JOIN aki_all a ON a.stay_id = b.stay_id
LEFT JOIN nadir n ON n.stay_id = b.stay_id
LEFT JOIN stage_at_t0 st ON st.stay_id = b.stay_id
LEFT JOIN ckd_hist ck ON ck.stay_id = b.stay_id
WHERE a.aki_offset BETWEEN 0 AND 7*1440;
SELECT 'E2_with_AKI_T0' AS step, COUNT(*) AS n,
       SUM(by_creat::int) AS n_by_creat, SUM(by_uo::int) AS n_by_uo,
       SUM((by_uo AND NOT by_creat)::int) AS n_uo_only,
       ROUND(AVG(t0_offset)/1440.0,2) AS mean_t0_day,
       ROUND(AVG(baseline_creat)::numeric,2) AS mean_base,
       ROUND(100.0*AVG((first_aki_stage=1)::int),1) AS pct_s1,
       ROUND(100.0*AVG((first_aki_stage=3)::int),1) AS pct_s3
FROM eicu_cohort_aki;   -- ["E2_with_AKI_T0"	13383	12325	3680	1058	0.74	1.39	63.6	19.7] 


-- E3. Exclude ESRD/chronic-dialysis dependence, kidney transplant, pregnancy -> eicu_cohort_final
DROP TABLE IF EXISTS eicu_cohort_final;
CREATE TABLE eicu_cohort_final AS
WITH ph AS (
  SELECT b.stay_id,
    MAX(CASE WHEN LOWER(h.pasthistoryvalue) ~ 'dialysis|hemodialysis|esrd|end-stage renal|end stage renal'
              AND LOWER(h.pasthistoryvalue) !~ 'not on dialysis|not currently dialyzed|no dialysis'
             THEN 1 ELSE 0 END) AS esrd_ph,
    MAX(CASE WHEN LOWER(h.pasthistoryvalue) ~ 'renal transplant|kidney transplant' THEN 1 ELSE 0 END) AS tx_ph
  FROM eicu_cohort_aki b LEFT JOIN pasthistory h ON h.patientunitstayid = b.stay_id
  GROUP BY b.stay_id
),
dx AS (
  SELECT b.stay_id,
    MAX(CASE WHEN LOWER(d.diagnosisstring) ~ 'end-stage renal|esrd|dialysis dependent|chronic.*dialysis' THEN 1 ELSE 0 END) AS esrd_dx,
    MAX(CASE WHEN LOWER(d.diagnosisstring) ~ 'renal transplant|kidney transplant' THEN 1 ELSE 0 END) AS tx_dx,
    MAX(CASE WHEN LOWER(d.diagnosisstring) ~ 'pregnan|obstetric|postpartum|peripartum' THEN 1 ELSE 0 END) AS preg_dx
  FROM eicu_cohort_aki b LEFT JOIN diagnosis d ON d.patientunitstayid = b.stay_id
  GROUP BY b.stay_id
),
chronic_dial AS (   -- [v7] chronic-dialysis evidence (treatment table): dialysis explicitly labeled "chronic renal failure", or already on dialysis before T0 (before AKI onset)
  SELECT b.stay_id,
    MAX(CASE WHEN LOWER(t.treatmentstring) ~ 'dialysis|hemodialysis|hemofiltration|cvvh|cavhd|crrt|renal replacement|continuous renal'
              AND (LOWER(t.treatmentstring) ~ 'chronic renal failure' OR t.treatmentoffset < b.t0_offset)
             THEN 1 ELSE 0 END) AS chronic_dial
  FROM eicu_cohort_aki b LEFT JOIN treatment t ON t.patientunitstayid = b.stay_id
  GROUP BY b.stay_id
),
flagged AS (
  SELECT b.*,
    GREATEST(COALESCE(ph.esrd_ph,0), COALESCE(dx.esrd_dx,0)) AS excl_esrd,
    GREATEST(COALESCE(ph.tx_ph,0),   COALESCE(dx.tx_dx,0))   AS excl_transplant,
    COALESCE(dx.preg_dx,0)                                   AS excl_pregnancy,
    COALESCE(cd.chronic_dial,0)                              AS excl_chronic_dial
  FROM eicu_cohort_aki b
  LEFT JOIN ph ON ph.stay_id = b.stay_id
  LEFT JOIN dx ON dx.stay_id = b.stay_id
  LEFT JOIN chronic_dial cd ON cd.stay_id = b.stay_id
)
SELECT * FROM flagged
WHERE excl_esrd = 0 AND excl_transplant = 0 AND excl_pregnancy = 0 AND excl_chronic_dial = 0;
SELECT 'E3_after_exclusions' AS step, COUNT(*) AS n FROM eicu_cohort_final;   -- [ref] ?
-- number excluded for chronic dialysis (treatment hit, within the AKI cohort)
SELECT COUNT(DISTINCT b.stay_id) AS n_chronic_dial_caught
FROM eicu_cohort_aki b JOIN treatment t ON t.patientunitstayid = b.stay_id
WHERE LOWER(t.treatmentstring) ~ 'dialysis|hemodialysis|hemofiltration|cvvh|cavhd|crrt|renal replacement|continuous renal'
  AND (LOWER(t.treatmentstring) ~ 'chronic renal failure' OR t.treatmentoffset < b.t0_offset);   -- [824] ?


-- E4. eGFR (CKD-EPI 2021) + CKD (eGFR < 60 or ckd_hist)
DROP TABLE IF EXISTS eicu_cohort_baseline;
CREATE TABLE eicu_cohort_baseline AS
SELECT f.*,
  ROUND((142.0
    * POWER(LEAST   (f.baseline_creat/(CASE WHEN f.gender='Female' THEN 0.7 ELSE 0.9 END),1.0),
            (CASE WHEN f.gender='Female' THEN -0.241 ELSE -0.302 END))
    * POWER(GREATEST(f.baseline_creat/(CASE WHEN f.gender='Female' THEN 0.7 ELSE 0.9 END),1.0), -1.200)
    * POWER(0.9938, f.age_at_icu)
    * (CASE WHEN f.gender='Female' THEN 1.012 ELSE 1.0 END))::numeric,1) AS baseline_egfr,
  CASE WHEN (142.0
    * POWER(LEAST   (f.baseline_creat/(CASE WHEN f.gender='Female' THEN 0.7 ELSE 0.9 END),1.0),
            (CASE WHEN f.gender='Female' THEN -0.241 ELSE -0.302 END))
    * POWER(GREATEST(f.baseline_creat/(CASE WHEN f.gender='Female' THEN 0.7 ELSE 0.9 END),1.0), -1.200)
    * POWER(0.9938, f.age_at_icu)
    * (CASE WHEN f.gender='Female' THEN 1.012 ELSE 1.0 END)) < 60
       OR f.ckd_hist = 1 THEN 1 ELSE 0 END AS ckd
FROM eicu_cohort_final f;
SELECT 'E4_baseline' AS step, COUNT(*) AS n,
       ROUND(AVG(baseline_creat)::numeric,2) AS mean_base,
       ROUND(AVG(baseline_egfr)::numeric,1) AS mean_egfr,
       ROUND(100.0*AVG(ckd)::numeric,1) AS pct_ckd
FROM eicu_cohort_baseline;   -- ["E4_baseline"	11855	1.13	77.9	37.4] ?


-- ----------------------------------------------------------------------------
-- Calibration (recommended to run separately first): urine-output record structure -- celllabel + value range
--   If mean_vol is a few hundred (per-interval volume) -> our SUM approach is correct; if it is in the thousands and rising over time -> it may be cumulative and must be changed
-- ----------------------------------------------------------------------------
SELECT celllabel, COUNT(*) AS n, COUNT(DISTINCT patientunitstayid) AS n_stays,
       ROUND(AVG(cellvaluenumeric)::numeric,0) AS mean_vol,
       ROUND(PERCENTILE_CONT(0.5) WITHIN GROUP (ORDER BY cellvaluenumeric)::numeric,0) AS median_vol,
       ROUND(MAX(cellvaluenumeric)::numeric,0) AS max_vol
FROM intakeoutput
WHERE LOWER(celllabel) ~ 'urine|foley|void'
GROUP BY celllabel ORDER BY n DESC LIMIT 30;



-- ============================================================================
--  eICU extraction pipeline | STAGE 2: exposure PPI vs. H2RA + new-user determination
--  ---------------------------------------------------------------------------
--  Aligned with MIMIC's cohort_exposure:
--    - exposure_group = class of the first acid suppressant within [T0, T0+24h] (PPI/H2RA/neither)
--    - drug sources = medication (oral/IV push) + infusiondrug (IV infusion pantoprazole/esomeprazole -> PPI)
--    - both_in_window = both classes appear within the window
--    - [v7] newuser = no in-ICU acid suppressant before T0 AND no home acid suppressant in admissiondrug
--      (eICU has a home-medication history, so it identifies pre-admission chronic users better than MIMIC -> a stronger bias probe)
--  Time in offset (minutes since ICU admission); prerequisite: eicu_cohort_baseline already built.
-- ============================================================================

DROP TABLE IF EXISTS eicu_cohort_exposure;
CREATE TABLE eicu_cohort_exposure AS
WITH acid_med AS (   -- medication table (drop cancelled orders)
  SELECT b.stay_id, m.drugstartoffset AS off,
    CASE
      WHEN LOWER(m.drugname) ~ 'pantoprazole|omeprazole|esomeprazole|lansoprazole|rabeprazole|dexlansoprazole|protonix|prilosec|nexium|prevacid|aciphex|dexilant' THEN 'PPI'
      WHEN LOWER(m.drugname) ~ 'famotidine|ranitidine|cimetidine|nizatidine|pepcid|zantac|tagamet|axid' THEN 'H2RA'
    END AS acid_class
  FROM eicu_cohort_baseline b
  JOIN medication m ON m.patientunitstayid = b.stay_id
   AND m.drugordercancelled = 'No' AND m.drugstartoffset IS NOT NULL
  WHERE LOWER(m.drugname) ~ 'prazole|tidine|protonix|prilosec|nexium|prevacid|aciphex|dexilant|pepcid|zantac|tagamet|axid'
),
acid_inf AS (   -- infusiondrug IV infusions (all are PPIs)
  SELECT b.stay_id, i.infusionoffset AS off, 'PPI'::text AS acid_class
  FROM eicu_cohort_baseline b
  JOIN infusiondrug i ON i.patientunitstayid = b.stay_id AND i.infusionoffset IS NOT NULL
  WHERE LOWER(i.drugname) ~ 'pantoprazole|esomeprazole|protonix'
),
acid AS (
  SELECT stay_id, off, acid_class FROM acid_med WHERE acid_class IS NOT NULL
  UNION ALL
  SELECT stay_id, off, acid_class FROM acid_inf
),
acid_t0 AS (   -- attach each patient's T0
  SELECT a.stay_id, a.off, a.acid_class, b.t0_offset
  FROM acid a JOIN eicu_cohort_baseline b ON b.stay_id = a.stay_id
),
first_in_window AS (   -- first acid suppressant within the window (break ties alphabetically)
  SELECT DISTINCT ON (stay_id)
    stay_id, acid_class AS first_acid_class, off AS first_acid_offset
  FROM acid_t0
  WHERE off >= t0_offset AND off <= t0_offset + 1440
  ORDER BY stay_id, off ASC, acid_class
),
both_flag AS (
  SELECT stay_id, CASE WHEN COUNT(DISTINCT acid_class) > 1 THEN 1 ELSE 0 END AS both_in_window
  FROM acid_t0 WHERE off >= t0_offset AND off <= t0_offset + 1440
  GROUP BY stay_id
),
prior_icu AS (   -- in-ICU acid suppressant before T0
  SELECT DISTINCT stay_id FROM acid_t0 WHERE off < t0_offset
),
prior_home AS (   -- home acid suppressant from admissiondrug (pre-admission medication history)
  SELECT DISTINCT b.stay_id
  FROM eicu_cohort_baseline b
  JOIN admissiondrug ad ON ad.patientunitstayid = b.stay_id
  WHERE LOWER(ad.drugname) ~ 'pantoprazole|omeprazole|esomeprazole|lansoprazole|rabeprazole|dexlansoprazole|protonix|prilosec|nexium|prevacid|aciphex|dexilant|famotidine|ranitidine|cimetidine|nizatidine|pepcid|zantac|tagamet|axid'
)
SELECT b.*,
  COALESCE(fw.first_acid_class, 'neither') AS exposure_group,
  fw.first_acid_offset,
  CASE WHEN fw.first_acid_offset IS NOT NULL
       THEN ROUND(((fw.first_acid_offset - b.t0_offset)/60.0)::numeric, 1) END AS hours_to_exposure,
  COALESCE(bf.both_in_window, 0) AS both_in_window,
  CASE WHEN pi.stay_id IS NULL AND ph.stay_id IS NULL THEN 1 ELSE 0 END AS newuser
FROM eicu_cohort_baseline b
LEFT JOIN first_in_window fw ON fw.stay_id = b.stay_id
LEFT JOIN both_flag       bf ON bf.stay_id = b.stay_id
LEFT JOIN prior_icu       pi ON pi.stay_id = b.stay_id
LEFT JOIN prior_home      ph ON ph.stay_id = b.stay_id;

-- Check (1): distribution of the three groups + new-user/both counts (vs. MIMIC: neither 15,230 / PPI 3,439 / H2RA 3,075)
SELECT exposure_group, COUNT(*) AS n,
       SUM(newuser) AS n_newuser,
       SUM(both_in_window) AS n_both,
       ROUND(AVG(hours_to_exposure)::numeric,1) AS mean_hrs_to_exp
FROM eicu_cohort_exposure
GROUP BY exposure_group ORDER BY n DESC;
-- "neither"	8923	4874	0	
-- "PPI"	2055	1204	29	6.5
-- "H2RA"	877	605	70	7.1



-- Check (2): head-to-head analyzable counts (PPI/H2RA after excluding both_in_window)
SELECT
  COUNT(*) FILTER (WHERE exposure_group='PPI'  AND both_in_window=0) AS n_ppi,
  COUNT(*) FILTER (WHERE exposure_group='H2RA' AND both_in_window=0) AS n_h2ra,
  COUNT(*) FILTER (WHERE exposure_group='neither')                  AS n_neither,
  COUNT(*) FILTER (WHERE exposure_group IN ('PPI','H2RA') AND both_in_window=0 AND newuser=1) AS n_initiator_newuser
FROM eicu_cohort_exposure;
-- 2026	807	8923	1747


-- ============================================================================
--  eICU extraction pipeline | STAGE 3: in-hospital outcomes (exact mirror of MIMIC B1-B6, follow-up switched to in-hospital)
--  ---------------------------------------------------------------------------
--  Follow-up window = [T0, discharge] (hospitaldischargeoffset, minutes); 98% discharged within 30d -> in-hospital ~ 30d.
--  Outcomes (all in-hospital), each defined to match MIMIC:
--    - death_inhosp = hospitaldischargestatus='Expired'
--    - rrt_init = first dialysis after T0 (treatment table); rrt_dep = dialysis within 48h before discharge (-> a MAKE component)
--    - non_recovery = last creatinine > 1.5x baseline OR last eGFR < 0.75x baseline (no creatinine -> NULL)
--    - make_inhosp = death OR rrt_dep OR non_recovery (non_recovery unknown -> NULL)  [primary outcome]
--    - progression = reaching KDIGO 3 after T0 (only if T0 stage < 3; death competes)
--    - recovery = creatinine falling to <= 1.5x baseline after T0, sustained with no rebound and no dialysis in the prior 48h (only if ever > 1.5x; death competes)
--  CKD-EPI 2021; prerequisite: eicu_cohort_exposure already built.
-- ============================================================================

DROP TABLE IF EXISTS eicu_cohort_outcomes;
CREATE TABLE eicu_cohort_outcomes AS
WITH post_creat AS (   -- creatinine within [T0, discharge]
  SELECT e.stay_id, l.labresultoffset AS off, l.labresult AS creat,
         e.t0_offset, e.baseline_creat, e.baseline_egfr, e.gender, e.age_at_icu
  FROM eicu_cohort_exposure e
  JOIN lab l ON l.patientunitstayid = e.stay_id
   AND LOWER(l.labname) = 'creatinine'
   AND l.labresult IS NOT NULL AND l.labresult > 0.1 AND l.labresult < 50
   AND l.labresultoffset >= e.t0_offset
   AND l.labresultoffset <= e.hospitaldischargeoffset
),
last_creat AS (   -- last creatinine
  SELECT DISTINCT ON (stay_id) stay_id, creat AS last_creat,
         baseline_creat, baseline_egfr, gender, age_at_icu
  FROM post_creat ORDER BY stay_id, off DESC
),
nonrec_flag AS (   -- non-recovery (last eGFR by CKD-EPI 2021)
  SELECT stay_id,
    CASE WHEN last_creat > 1.5*baseline_creat
           OR (142.0
               * POWER(LEAST   (last_creat/(CASE WHEN gender='Female' THEN 0.7 ELSE 0.9 END),1.0),
                       (CASE WHEN gender='Female' THEN -0.241 ELSE -0.302 END))
               * POWER(GREATEST(last_creat/(CASE WHEN gender='Female' THEN 0.7 ELSE 0.9 END),1.0), -1.200)
               * POWER(0.9938, age_at_icu)
               * (CASE WHEN gender='Female' THEN 1.012 ELSE 1.0 END)) < 0.75*baseline_egfr
         THEN 1 ELSE 0 END AS non_recovery
  FROM last_creat
),
-- dialysis events (treatment table; chronic-dialysis patients already excluded -> this is acute RRT)
dial AS (
  SELECT e.stay_id, t.treatmentoffset AS dial_off, e.t0_offset, e.hospitaldischargeoffset AS dc_off
  FROM eicu_cohort_exposure e
  JOIN treatment t ON t.patientunitstayid = e.stay_id
  WHERE LOWER(t.treatmentstring) ~ 'dialysis|hemodialysis|hemofiltration|cvvh|cavhd|crrt|renal replacement|continuous renal'
),
rrt AS (
  SELECT stay_id,
    MIN(CASE WHEN dial_off > t0_offset THEN dial_off END) AS rrt_init_offset,
    MAX((dial_off >= dc_off - 2880 AND dial_off <= dc_off AND dial_off > t0_offset)::int) AS rrt_dep_inhosp
  FROM dial GROUP BY stay_id
),
-- progression to KDIGO 3
prog AS (
  SELECT stay_id, MIN(off) AS prog_offset
  FROM post_creat
  WHERE (creat >= 3.0*baseline_creat OR creat >= 4.0) AND off > t0_offset
  GROUP BY stay_id
),
-- recovery: was once > 1.5x baseline, then fell to <= 1.5x with no rebound and no dialysis in the prior 48h
above AS (
  SELECT stay_id, MAX(off) AS last_above_off
  FROM post_creat WHERE creat > 1.5*baseline_creat GROUP BY stay_id
),
recov AS (
  SELECT DISTINCT ON (pc.stay_id) pc.stay_id, pc.off AS rec_offset
  FROM post_creat pc JOIN above a ON a.stay_id = pc.stay_id
  WHERE pc.creat <= 1.5*pc.baseline_creat AND pc.off > a.last_above_off
    AND NOT EXISTS (SELECT 1 FROM dial d
                    WHERE d.stay_id = pc.stay_id
                      AND d.dial_off >= pc.off - 2880 AND d.dial_off <= pc.off)
  ORDER BY pc.stay_id, pc.off ASC
)
SELECT e.*,
  -- death
  CASE WHEN e.hospitaldischargestatus = 'Expired' THEN 1 ELSE 0 END AS death_inhosp,
  ROUND((e.hospitaldischargeoffset/1440.0)::numeric,2) AS dc_day,
  -- RRT
  CASE WHEN r.rrt_init_offset IS NOT NULL THEN 1 ELSE 0 END AS rrt_init,
  CASE WHEN r.rrt_init_offset IS NOT NULL
       THEN ROUND(((r.rrt_init_offset - e.t0_offset)/1440.0)::numeric,2) END AS days_to_rrt,
  COALESCE(r.rrt_dep_inhosp,0) AS rrt_dep_inhosp,
  -- non-recovery + MAKE (three-valued logic: if death/RRT is true then MAKE=1; if non-recovery is unknown and the rest are 0 then NULL)
  nf.non_recovery,
  CASE WHEN (CASE WHEN e.hospitaldischargestatus='Expired' THEN 1 ELSE 0 END)=1
            OR COALESCE(r.rrt_dep_inhosp,0)=1 OR nf.non_recovery=1 THEN 1
       WHEN nf.non_recovery IS NULL THEN NULL ELSE 0 END AS make_inhosp,
  -- progression (only T0 stage < 3)
  CASE WHEN e.first_aki_stage < 3 THEN 1 ELSE 0 END AS prog_eligible,
  CASE WHEN e.first_aki_stage < 3 AND p.prog_offset IS NOT NULL THEN 1 ELSE 0 END AS prog_event,
  CASE WHEN e.first_aki_stage < 3 AND p.prog_offset IS NOT NULL
       THEN ROUND(((p.prog_offset - e.t0_offset)/1440.0)::numeric,2) END AS days_to_prog,
  -- recovery (only if ever > 1.5x baseline)
  CASE WHEN ab.last_above_off IS NOT NULL THEN 1 ELSE 0 END AS recovery_eligible,
  CASE WHEN ab.last_above_off IS NOT NULL AND rc.rec_offset IS NOT NULL THEN 1 ELSE 0 END AS recovery_event,
  CASE WHEN rc.rec_offset IS NOT NULL
       THEN ROUND(((rc.rec_offset - e.t0_offset)/1440.0)::numeric,2) END AS days_to_recovery
FROM eicu_cohort_exposure e
LEFT JOIN nonrec_flag nf ON nf.stay_id = e.stay_id
LEFT JOIN rrt   r  ON r.stay_id  = e.stay_id
LEFT JOIN prog  p  ON p.stay_id  = e.stay_id
LEFT JOIN above ab ON ab.stay_id = e.stay_id
LEFT JOIN recov rc ON rc.stay_id = e.stay_id;


-- Check (1): primary outcome + components (vs. MIMIC in-hospital: death ~18%, MAKE30 ~37.7%, RRT init ~5.6%)
SELECT COUNT(*) AS n,
  ROUND(100.0*AVG(death_inhosp),1) AS pct_death,
  ROUND(100.0*AVG(rrt_init),1) AS pct_rrt_init,
  ROUND(100.0*AVG(rrt_dep_inhosp),1) AS pct_rrt_dep,
  ROUND(100.0*AVG(non_recovery::numeric),1) AS pct_nonrec,
  SUM((make_inhosp IS NULL)::int) AS n_make_null,
  ROUND(100.0*SUM((make_inhosp=1)::int)/NULLIF(SUM((make_inhosp IS NOT NULL)::int),0),1) AS pct_make
FROM eicu_cohort_outcomes;
-- 11855	20.8	6.4	2.5	27.4	72	37.2

-- Check (2): MAKE + death by exposure group (head-to-head preview, excluding both_in_window)
SELECT exposure_group, COUNT(*) AS n,
  ROUND(100.0*AVG(death_inhosp),1) AS pct_death,
  ROUND(100.0*SUM((make_inhosp=1)::int)/NULLIF(SUM((make_inhosp IS NOT NULL)::int),0),1) AS pct_make
FROM eicu_cohort_outcomes
WHERE both_in_window = 0
GROUP BY exposure_group ORDER BY n DESC;
-- "neither"	8923	21.1	38.2
-- "PPI"	2026	20.2	33.7
-- "H2RA"	807	19.8	34.7

-- Check (3): progression / recovery (competing-risk secondary outcomes; vs. MIMIC progression ~23.4%, recovery ~49.5%)
SELECT
  SUM(prog_eligible) AS n_prog_elig, SUM(prog_event) AS n_prog,
  ROUND(100.0*SUM(prog_event)/NULLIF(SUM(prog_eligible),0),1) AS pct_prog,
  SUM(recovery_eligible) AS n_rec_elig, SUM(recovery_event) AS n_rec,
  ROUND(100.0*SUM(recovery_event)/NULLIF(SUM(recovery_eligible),0),1) AS pct_rec
FROM eicu_cohort_outcomes;
-- 10207	1098	10.8	9288	6826	73.5

-- Calibration (recommended to run separately first): dialysis treatmentstring hits -- confirm RRT capture is reasonable
SELECT t.treatmentstring, COUNT(*) AS n, COUNT(DISTINCT t.patientunitstayid) AS n_stays
FROM treatment t
WHERE LOWER(t.treatmentstring) ~ 'dialysis|hemodialysis|hemofiltration|cvvh|cavhd|crrt|renal replacement|continuous renal'
GROUP BY t.treatmentstring ORDER BY n DESC LIMIT 25;



-- ============================================================================
--  eICU extraction pipeline | STAGE 4: covariates (for the propensity score) + assemble the final analytic_cohort
--  ---------------------------------------------------------------------------
--  Aligned with the MIMIC propensity-score variables, taking the subset available in eICU. Already present (in the cohort): age/sex/race/weight/
--  unittype/baseline_egfr/ckd/first_aki_stage. This stage adds:
--    - Severity: apachescore, acutephysiologyscore, predictedhospitalmortality (apachepatientresult, deduplicated)
--    - Comorbidities: cirrhosis/hepaticfailure/diabetes/immunosupp/aids/leukemia/lymphoma/metastatic (apachepredvar)
--    - Day-1 physiology: HR/MAP/RR/Temp/GCS/WBC/albumin/bilirubin/BUN/glucose/sodium/hematocrit (apacheapsvar)
--    - Interventions: mechanical ventilation ventday1 (apachepredvar), vasoactive drugs (infusiondrug, around T0)
--    - Confounding-by-indication: GI bleeding (diagnosis), liver disease (cirrhosis or hepaticfailure)
--    - Lactate: lab (worst value around T0; not in the APACHE tables)
--  Missing values left as NULL -> MICE imputation in R (consistent with MIMIC). Prerequisite: eicu_cohort_outcomes already built.
-- ============================================================================

-- C1. APACHE: scores + comorbidities + day-1 physiology ----------------------------------------
DROP TABLE IF EXISTS eicu_cov_apache;
CREATE TABLE eicu_cov_apache AS
WITH apr AS (   -- deduplicate scores: one row per patient, prefer a valid score / IVa
  SELECT DISTINCT ON (patientunitstayid) patientunitstayid,
    NULLIF(apachescore, -1) AS apachescore,
    NULLIF(acutephysiologyscore, -1) AS acutephysiologyscore,
    CASE WHEN predictedhospitalmortality ~ '^[0-9.]+$'
         THEN predictedhospitalmortality::numeric END AS pred_hosp_mort
  FROM apachepatientresult
  ORDER BY patientunitstayid, (apachescore > 0) DESC, apacheversion DESC
)
SELECT b.stay_id,
  apr.apachescore, apr.acutephysiologyscore, apr.pred_hosp_mort,
  -- comorbidities (1 = present; -1/NULL -> 0)
  CASE WHEN pv.cirrhosis=1        THEN 1 ELSE 0 END AS cm_cirrhosis,
  CASE WHEN pv.hepaticfailure=1   THEN 1 ELSE 0 END AS cm_hepaticfailure,
  CASE WHEN pv.diabetes=1         THEN 1 ELSE 0 END AS cm_diabetes,
  CASE WHEN pv.immunosuppression=1 THEN 1 ELSE 0 END AS cm_immunosupp,
  CASE WHEN pv.aids=1             THEN 1 ELSE 0 END AS cm_aids,
  CASE WHEN pv.leukemia=1         THEN 1 ELSE 0 END AS cm_leukemia,
  CASE WHEN pv.lymphoma=1         THEN 1 ELSE 0 END AS cm_lymphoma,
  CASE WHEN pv.metastaticcancer=1 THEN 1 ELSE 0 END AS cm_metastatic,
  CASE WHEN pv.cirrhosis=1 OR pv.hepaticfailure=1 THEN 1 ELSE 0 END AS liver_disease,
  CASE WHEN pv.ventday1=1         THEN 1 ELSE 0 END AS vent_day1,
  -- day-1 physiology (-1 -> NULL)
  NULLIF(av.heartrate,-1)       AS hr,
  NULLIF(av.meanbp,-1)          AS map,
  NULLIF(av.respiratoryrate,-1) AS rr,
  NULLIF(av.temperature,-1)     AS temp,
  (NULLIF(av.eyes,-1) + NULLIF(av.motor,-1) + NULLIF(av.verbal,-1)) AS gcs,
  NULLIF(av.wbc,-1)        AS wbc,
  NULLIF(av.albumin,-1)    AS albumin,
  NULLIF(av.bilirubin,-1)  AS bilirubin,
  NULLIF(av.bun,-1)        AS bun,
  NULLIF(av.glucose,-1)    AS glucose,
  NULLIF(av.sodium,-1)     AS sodium,
  NULLIF(av.hematocrit,-1) AS hematocrit
FROM eicu_cohort_outcomes b
LEFT JOIN apr           ON apr.patientunitstayid = b.stay_id
LEFT JOIN apachepredvar pv ON pv.patientunitstayid = b.stay_id
LEFT JOIN apacheapsvar  av ON av.patientunitstayid = b.stay_id;

-- C2. Vasoactive drugs (infusiondrug, [T0-24h, T0+6h]) -----------------------------
DROP TABLE IF EXISTS eicu_cov_vaso;
CREATE TABLE eicu_cov_vaso AS
SELECT e.stay_id,
  MAX(CASE WHEN i.infusionoffset BETWEEN e.t0_offset - 1440 AND e.t0_offset + 360 THEN 1 ELSE 0 END) AS vasoactive
FROM eicu_cohort_outcomes e
JOIN infusiondrug i ON i.patientunitstayid = e.stay_id
WHERE LOWER(i.drugname) ~ 'norepinephrine|epinephrine|dopamine|dobutamine|vasopressin|phenylephrine|levophed|neo-?synephrine|milrinone|vasopressor'
GROUP BY e.stay_id;

-- C3. GI bleeding (diagnosis; the core confounding-by-indication variable for PPI) ----------------------
DROP TABLE IF EXISTS eicu_cov_gibleed;
CREATE TABLE eicu_cov_gibleed AS
SELECT e.stay_id,
  MAX(CASE WHEN LOWER(d.diagnosisstring) ~ 'gi bleed|gastrointestinal bleed|gastrointestinal hemorrhage|gi hemorrhage|peptic ulcer|duodenal ulcer|gastric ulcer|melena|hematemesis|variceal|esophageal varices'
           THEN 1 ELSE 0 END) AS gi_bleed
FROM eicu_cohort_outcomes e
JOIN diagnosis d ON d.patientunitstayid = e.stay_id
GROUP BY e.stay_id;

-- C4. Lactate (lab, worst/maximum value in [T0-24h, T0+6h]; not in the APACHE tables) -------------------
DROP TABLE IF EXISTS eicu_cov_lactate;
CREATE TABLE eicu_cov_lactate AS
SELECT e.stay_id, MAX(l.labresult) AS lactate
FROM eicu_cohort_outcomes e
JOIN lab l ON l.patientunitstayid = e.stay_id
 AND LOWER(l.labname) ~ 'lactate'
 AND l.labresult IS NOT NULL AND l.labresult > 0 AND l.labresult < 50
 AND l.labresultoffset BETWEEN e.t0_offset - 1440 AND e.t0_offset + 360
GROUP BY e.stay_id;

-- ============================================================================
-- Assemble the final analytic_cohort (cohort + exposure + in-hospital outcomes + covariates)
-- ============================================================================
DROP TABLE IF EXISTS eicu_analytic_cohort;
CREATE TABLE eicu_analytic_cohort AS
SELECT o.*,
  ap.apachescore, ap.acutephysiologyscore, ap.pred_hosp_mort,
  ap.cm_cirrhosis, ap.cm_hepaticfailure, ap.cm_diabetes, ap.cm_immunosupp,
  ap.cm_aids, ap.cm_leukemia, ap.cm_lymphoma, ap.cm_metastatic,
  ap.liver_disease, ap.vent_day1,
  ap.hr, ap.map, ap.rr, ap.temp, ap.gcs, ap.wbc, ap.albumin,
  ap.bilirubin, ap.bun, ap.glucose, ap.sodium, ap.hematocrit,
  COALESCE(v.vasoactive,0) AS vasoactive,
  COALESCE(g.gi_bleed,0)   AS gi_bleed,
  lc.lactate
FROM eicu_cohort_outcomes o
LEFT JOIN eicu_cov_apache  ap ON ap.stay_id = o.stay_id
LEFT JOIN eicu_cov_vaso    v  ON v.stay_id  = o.stay_id
LEFT JOIN eicu_cov_gibleed g  ON g.stay_id  = o.stay_id
LEFT JOIN eicu_cov_lactate lc ON lc.stay_id = o.stay_id;


-- Check (1): covariate means + missingness (to plan MICE in R; and confirm physiology values are real values, not APACHE points)
SELECT COUNT(*) AS n,
  ROUND(AVG(apachescore)::numeric,1) AS apache, ROUND(100.0*AVG((apachescore IS NULL)::int),1) AS miss_apache,
  ROUND(AVG(hr)::numeric,0) AS hr, ROUND(AVG(map)::numeric,0) AS map, ROUND(AVG(gcs)::numeric,1) AS gcs,
  ROUND(AVG(albumin)::numeric,2) AS alb, ROUND(100.0*AVG((albumin IS NULL)::int),1) AS miss_alb,
  ROUND(AVG(bilirubin)::numeric,2) AS bili,
  ROUND(AVG(lactate)::numeric,2) AS lact, ROUND(100.0*AVG((lactate IS NULL)::int),1) AS miss_lact,
  ROUND(100.0*AVG(vent_day1),1) AS pct_vent, ROUND(100.0*AVG(vasoactive),1) AS pct_vaso,
  ROUND(100.0*AVG(liver_disease),1) AS pct_liver, ROUND(100.0*AVG(gi_bleed),1) AS pct_gib,
  ROUND(100.0*AVG(cm_diabetes),1) AS pct_dm
FROM eicu_analytic_cohort;
-- 11855	76.9	9.9	114	77	12.0	2.43	38.2	1.66	3.95	29.9	34.8	30.7	3.4	5.5	23.5

-- Check (2): confounding-indication preview (PPI vs. H2RA, excluding both) -- whether the direction of confounding-by-indication matches MIMIC
SELECT exposure_group, COUNT(*) AS n,
  ROUND(100.0*AVG(gi_bleed),1) AS pct_gibleed,
  ROUND(100.0*AVG(liver_disease),1) AS pct_liver,
  ROUND(AVG(apachescore),1) AS apache,
  ROUND(100.0*AVG(vasoactive),1) AS pct_vaso,
  ROUND(100.0*AVG(vent_day1),1) AS pct_vent
FROM eicu_analytic_cohort
WHERE both_in_window=0 AND exposure_group IN ('PPI','H2RA')
GROUP BY exposure_group;
-- "H2RA"	807	2.6	2.6	81.3	44.4	44.7
--"PPI"	2026	12.1	5.4	82.3	36.1	34.2

-- Check (3): final-table dimensions (column count should be ~ outcomes table + 28 covariate columns)
SELECT COUNT(*) AS n_rows,
  (SELECT COUNT(*) FROM information_schema.columns WHERE table_name='eicu_analytic_cohort') AS n_cols
FROM eicu_analytic_cohort;
-- 11855	73