-- ############################################################################
-- SICdb | Part 1 FINAL: SA-AKI cohort construction + attrition table
-- ----------------------------------------------------------------------------
-- Design rationale:
--   - Target-trial emulation (PPI vs. H2RA), not a prediction model -> preserve causal temporality (sepsis first, AKI <=7d, T0 = AKI onset)
--   - Adopts transportable criteria from Luo et al. CKJ 2026 (SICdb SA-AKI n=1893):
--       [literature] first admission + first ICU; ICU stay >=24h; physiologic-boundary outlier cleaning
--   - AKI definition improves on that paper: does not use ICD codes (SICdb has only a primary diagnosis; N17 is unreliable),
--       instead uses the official KDIGO_AKI_168 (creatinine + urine output, a published algorithm) as the AKI gold standard,
--       and derives the T0 time point from laboratory creatinine / data_float_h urine output (needed for causal inference).
--
-- Confirmed constants (measured in this database):
--   creatinine laboratoryid IN (367,368) mg/dl | urine output dataid=725 ml/h |
--   weight weightonadmission /1000 = kg | Yes/No/Unk = 740/738/739 | female sex=736 |
--   KDIGO_AKI_168 = data_ref.fieldid 3113, refid 3116/3117/3118 = stage1/2/3 (3115=0) |
--   "Offset" columns are capitalized | antibiotics/drugs d_references.referencename='Drug'
--
-- AKI determination logic (dual):
--   Official determination = KDIGO_AKI_168 stage >= 1 (this decides whether it is AKI; comparable to the literature)
--   T0 time = the earliest qualifying Offset derived from the creatinine or urine-output method (for causal temporality and follow-up anchoring)
--   -> entry requirement: official KDIGO positive AND a derivable T0 AND T0 within 7d after sepsis
--
-- Output: sic_cohort (one row per SA-AKI case; includes t_sepsis / t0 / baseline_creat / criterion flags)
-- ############################################################################

-- ============================================================================
-- Step 0: base layer -- first ICU + adult + ICU stay >=24h, weight cleaning
--   [literature] first admission/first ICU + LOS >=24h; ICU-stay duration uses timeofstay (seconds)
-- ============================================================================
DROP TABLE IF EXISTS sic_base;
CREATE TABLE sic_base AS
WITH ranked AS (
    SELECT caseid, patientid, ageonadmission, sex, admissionyear,
           admissionformhassepsis, icd10main, hoursofcrrt, offsetofdeath,
           timeofstay, icuoffset, weightonadmission, heightonadmission,
           ROW_NUMBER() OVER (PARTITION BY patientid
                              ORDER BY offsetafterfirstadmission ASC, caseid ASC) AS rn
    FROM cases
)
SELECT caseid AS stay_id, patientid, ageonadmission AS age, sex, admissionyear,
       admissionformhassepsis, icd10main, hoursofcrrt, offsetofdeath, timeofstay, icuoffset,
       CASE WHEN weightonadmission/1000.0 BETWEEN 30 AND 250
            THEN ROUND((weightonadmission/1000.0)::numeric,1) END AS weight_kg,   -- physiologic-boundary cleaning
       CASE WHEN heightonadmission BETWEEN 130 AND 220
            THEN heightonadmission END AS height_cm
FROM ranked
WHERE rn = 1
  AND ageonadmission >= 18
  AND timeofstay >= 24*3600;          -- [literature] ICU stay >=24h
-- Note: if timeofstay actually denotes total hospital length of stay rather than ICU duration, an ICU-stay metric can be used instead; here timeofstay is hospital seconds.
--     [VERIFY] if a strict "ICU >=24h" is required, confirm the semantics of timeofstay (Appendix 0 includes a check).

-- Appendix 0: confirm whether the semantics of timeofstay (days) are reasonable, to help check the ICU >=24h filter
SELECT ROUND(PERCENTILE_CONT(0.5) WITHIN GROUP (ORDER BY timeofstay/86400.0)::numeric,2) AS median_los_days,
       ROUND(MIN(timeofstay/86400.0)::numeric,2) AS min_d, ROUND(MAX(timeofstay/86400.0)::numeric,1) AS max_d,
       COUNT(*) AS n_after_24h_filter
FROM sic_base;   --2.22	1.00	126.0	14803

-- ============================================================================
-- Step 1: sepsis -- sustained antibiotics (Sepsis-3 suspected-infection proxy) + time anchor t_sepsis
--   [improvement] no longer "any single antibiotic dose": require early initiation (<=48h) and sustained delivery (coverage >=? or multiple doses),
--          filtering out SICdb's many single-dose perioperative prophylactic antibiotics, improving specificity and approximating Sepsis-3.
--   Criterion A (native sensitivity): AdmissionFormHasSepsis=740; Criterion C (main analysis): sustained antibiotics
-- ============================================================================
DROP TABLE IF EXISTS sic_sepsis;
CREATE TABLE sic_sepsis AS
WITH abx AS (   -- antibiotic administration events (ILIKE per column; within 7d)
    SELECT m.caseid, m."Offset" AS off, m.offsetdrugend AS off_end
    FROM medication m
    JOIN d_references d ON d.referenceglobalid = m.drugid
    WHERE d.referencename='Drug'
      AND m."Offset" >= 0 AND m."Offset" <= 7*86400
      AND (
        d.referencevalue ILIKE '%piperacillin%' OR d.referencevalue ILIKE '%tazobactam%' OR
        d.referencevalue ILIKE '%meropenem%'    OR d.referencevalue ILIKE '%imipenem%'  OR
        d.referencevalue ILIKE '%ceftriaxon%'   OR d.referencevalue ILIKE '%cefuroxim%' OR
        d.referencevalue ILIKE '%cefepim%'      OR d.referencevalue ILIKE '%cefotaxim%' OR
        d.referencevalue ILIKE '%ciprofloxacin%' OR d.referencevalue ILIKE '%levofloxacin%' OR
        d.referencevalue ILIKE '%moxifloxacin%' OR d.referencevalue ILIKE '%metronidazol%' OR
        d.referencevalue ILIKE '%vancomycin%'   OR d.referencevalue ILIKE '%gentamic%'  OR
        d.referencevalue ILIKE '%tobramycin%'   OR d.referencevalue ILIKE '%amikacin%'  OR
        d.referencevalue ILIKE '%linezolid%'    OR d.referencevalue ILIKE '%clindamycin%' OR
        d.referencevalue ILIKE '%ampicillin%'   OR d.referencevalue ILIKE '%amoxicillin%' OR
        d.referencevalue ILIKE '%sultamicillin%' OR d.referencevalue ILIKE '%fosfomycin%' OR
        d.referencevalue ILIKE '%tigecyclin%'   OR d.referencevalue ILIKE '%colistin%'
      )
),
-- [!] SICdb administration form: medication expands one continuous infusion into minute-level records [one row per 60 seconds].
--   So (a) coverage duration must use the MIN(off)->MAX(off) span, not a single offsetdrugend (only 60 larger than off -> collapses);
--      (b) COUNT(*) is "number of minutes", not "number of doses"; do not use it.
--   abx already takes all administrations within 7d (not limited to 48h) to compute the true span; first-dose timing is judged separately via abx_first48.
abx_span AS (
    SELECT caseid,
           MIN(off) AS first_off,
           (MAX(off) - MIN(off))/86400.0 AS span_days         -- first->last span (days) = sustained duration
    FROM abx
    GROUP BY caseid
),
abx_first48 AS (   -- whether the first dose is within 48h (septic, not late nosocomial)
    SELECT caseid, MIN(off) AS first48_off
    FROM abx WHERE off <= 2*86400
    GROUP BY caseid
)
SELECT b.stay_id,
       CASE WHEN b.admissionformhassepsis = 740 THEN 1 ELSE 0 END AS sepsis_A,
       0::int AS t_sepsis_A,
       -- Criterion C (final): sustained antibiotics = first dose initiated <=48h AND administration span >= 2 days (the classic Sepsis-3 "suspected infection":
       --   antimicrobial therapy sustained >=48h). Uses only the MIN/MAX(off) span, avoiding the collapse caused by minute-level expansion.
       CASE WHEN f.caseid IS NOT NULL AND sp.span_days >= 2 THEN 1 ELSE 0 END AS sepsis_C,
       sp.first_off AS t_sepsis_C,
       sp.span_days
FROM sic_base b
LEFT JOIN abx_span    sp ON sp.caseid = b.stay_id
LEFT JOIN abx_first48 f  ON f.caseid  = b.stay_id;

-- ============================================================================
-- Step 2: baseline creatinine -- SAP Section 8 three-tier hierarchy (pre-admission -> ICU early 24h -> MDRD back-calculation)
-- ============================================================================
DROP TABLE IF EXISTS sic_baseline;
CREATE TABLE sic_baseline AS
WITH cr AS (
    SELECT b.stay_id, l."Offset" AS off, l.laboratoryvalue AS v, b.age, b.sex
    FROM sic_base b
    JOIN laboratory l ON l.caseid = b.stay_id
    WHERE l.laboratoryid IN (367,368)
      AND l.laboratoryvalue IS NOT NULL
      AND l.laboratoryvalue > 0.1 AND l.laboratoryvalue < 20     -- physiologic-boundary cleaning
),
t1 AS (SELECT stay_id, MIN(v) AS v FROM cr WHERE off < 0                  GROUP BY stay_id),
t2 AS (SELECT stay_id, MIN(v) AS v FROM cr WHERE off BETWEEN 0 AND 86400  GROUP BY stay_id)
SELECT b.stay_id,
   COALESCE(t1.v, t2.v,
     ROUND((POWER(175.0*POWER(b.age,-0.203)*(CASE WHEN b.sex=736 THEN 0.742 ELSE 1.0 END)/75.0,
                  1.0/1.154))::numeric,2)
   ) AS baseline_creat,
   CASE WHEN t1.v IS NOT NULL THEN 'pre_icu'
        WHEN t2.v IS NOT NULL THEN 'icu_early24h'
        ELSE 'mdrd75' END AS baseline_src
FROM sic_base b
LEFT JOIN t1 ON t1.stay_id=b.stay_id
LEFT JOIN t2 ON t2.stay_id=b.stay_id;

-- ============================================================================
-- Step 3: T0 derivation -- the earliest qualifying Offset from the creatinine or urine-output method (for temporality / follow-up anchoring)
--   Note: T0 only provides the "time point"; "whether it is AKI" is decided by the official KDIGO_AKI_168 (Step 4).
-- ============================================================================
DROP TABLE IF EXISTS sic_t0;
CREATE TABLE sic_t0 AS
WITH
cr AS (
    SELECT b.stay_id, l."Offset" AS off, l.laboratoryvalue AS v, bl.baseline_creat
    FROM sic_base b
    JOIN sic_baseline bl ON bl.stay_id=b.stay_id
    JOIN laboratory l ON l.caseid=b.stay_id
    WHERE l.laboratoryid IN (367,368)
      AND l.laboratoryvalue IS NOT NULL
      AND l.laboratoryvalue > 0.1 AND l.laboratoryvalue < 20
      AND l."Offset" BETWEEN 0 AND 14*86400
),
cr_roll AS (
    SELECT stay_id, off, v, baseline_creat,
           v - MIN(v) OVER (PARTITION BY stay_id ORDER BY off
                  RANGE BETWEEN 172800 PRECEDING AND CURRENT ROW) AS rise48
    FROM cr
),
aki_cr AS (
    SELECT stay_id, MIN(off) AS off FROM cr_roll
    WHERE v >= 1.5*baseline_creat OR rise48 >= 0.3
    GROUP BY stay_id
),
uo AS (
    SELECT b.stay_id, d."Offset" AS off, d.val AS ml, b.weight_kg
    FROM sic_base b
    JOIN data_float_h d ON d.caseid=b.stay_id
    WHERE d.dataid=725 AND d.val IS NOT NULL AND d.val >= 0 AND d.val < 4000
      AND d."Offset" BETWEEN 0 AND 14*86400 AND b.weight_kg IS NOT NULL
),
uo_roll AS (
    SELECT stay_id, off, weight_kg,
           SUM(ml) OVER w AS ml_6h, (off - MIN(off) OVER w) AS span, COUNT(*) OVER w AS nh
    FROM uo
    WINDOW w AS (PARTITION BY stay_id ORDER BY off RANGE BETWEEN 21600 PRECEDING AND CURRENT ROW)
),
aki_uo AS (
    SELECT stay_id, MIN(off) AS off FROM uo_roll
    WHERE span >= 18000 AND nh >= 5 AND ml_6h < 0.5*weight_kg*(span/3600.0)
    GROUP BY stay_id
)
SELECT stay_id, MIN(off) AS t0_aki,
       BOOL_OR(src='cr') AS t0_by_creat, BOOL_OR(src='uo') AS t0_by_uo
FROM (
    SELECT stay_id, off, 'cr' AS src FROM aki_cr
    UNION ALL
    SELECT stay_id, off, 'uo' AS src FROM aki_uo
) z GROUP BY stay_id;  --8252

-- ============================================================================
-- Step 4: official KDIGO_AKI_168 (AKI gold standard) + exclusion flags
-- ============================================================================
DROP TABLE IF EXISTS sic_kref;
CREATE TABLE sic_kref AS
SELECT b.stay_id,
   CASE WHEN k.caseid IS NOT NULL THEN 1 ELSE 0 END AS kdigo168_pos,
   CASE WHEN b.icd10main ~* '^(N185|N186|Z992|Z49)'    THEN 1 ELSE 0 END AS excl_esrd,
   CASE WHEN b.icd10main ~* '^(Z940|T861)'             THEN 1 ELSE 0 END AS excl_transplant,
   CASE WHEN b.icd10main ~* '^(O|Z33|Z34|Z3A|Z37|Z39)' THEN 1 ELSE 0 END AS excl_pregnancy
FROM sic_base b
LEFT JOIN (SELECT DISTINCT caseid FROM data_ref
           WHERE fieldid=3113 AND refid IN (3116,3117,3118)) k
       ON k.caseid=b.stay_id;

-- ============================================================================
-- Step 5: merge -> sic_cohort
--   SA-AKI entry = official KDIGO positive AND has T0 AND T0 within [t_sepsis, t_sepsis+7d]
--   Main criterion = C (sustained antibiotics); criterion A (flag) retained for sensitivity analysis
-- ============================================================================
DROP TABLE IF EXISTS sic_cohort;
CREATE TABLE sic_cohort AS
SELECT
   b.stay_id, b.patientid, b.age, b.sex, b.admissionyear, b.weight_kg, b.height_cm,
   b.hoursofcrrt, b.offsetofdeath, b.timeofstay, b.icuoffset,
   s.sepsis_A, s.sepsis_C, s.t_sepsis_A, s.t_sepsis_C, s.span_days,
   bl.baseline_creat, bl.baseline_src,
   t.t0_aki AS t0, t.t0_by_creat, t.t0_by_uo,
   k.kdigo168_pos,
   (k.excl_esrd + k.excl_transplant + k.excl_pregnancy) AS excl_any,
   k.excl_esrd, k.excl_transplant, k.excl_pregnancy,
   -- SA-AKI temporality determination (official KDIGO positive + has T0 + T0 within 7d after sepsis)
   CASE WHEN s.sepsis_A=1 AND k.kdigo168_pos=1 AND t.t0_aki IS NOT NULL
         AND t.t0_aki BETWEEN s.t_sepsis_A AND s.t_sepsis_A + 7*86400 THEN 1 ELSE 0 END AS saaki_A,
   CASE WHEN s.sepsis_C=1 AND k.kdigo168_pos=1 AND t.t0_aki IS NOT NULL AND s.t_sepsis_C IS NOT NULL
         AND t.t0_aki BETWEEN s.t_sepsis_C AND s.t_sepsis_C + 7*86400 THEN 1 ELSE 0 END AS saaki_C
FROM sic_base b
JOIN sic_sepsis  s  ON s.stay_id=b.stay_id
JOIN sic_baseline bl ON bl.stay_id=b.stay_id
LEFT JOIN sic_t0   t ON t.stay_id=b.stay_id
JOIN sic_kref     k  ON k.stay_id=b.stay_id;

-- ============================================================================
-- [*] Attrition table (A vs. C side by side) -- with literature comparison anchors
-- ============================================================================
SELECT step, label, n_A, n_C FROM (
  VALUES
   (1,'All cases',                          (SELECT COUNT(*) FROM cases), (SELECT COUNT(*) FROM cases)),
   (2,'First ICU + adult >=18 + ICU stay >=24h',          (SELECT COUNT(*) FROM sic_base), (SELECT COUNT(*) FROM sic_base)),
   (3,'+ sepsis (A=flag / C=sustained antibiotics)',
       (SELECT COUNT(*) FROM sic_cohort WHERE sepsis_A=1),
       (SELECT COUNT(*) FROM sic_cohort WHERE sepsis_C=1)),
   (4,'+ official KDIGO positive',
       (SELECT COUNT(*) FROM sic_cohort WHERE sepsis_A=1 AND kdigo168_pos=1),
       (SELECT COUNT(*) FROM sic_cohort WHERE sepsis_C=1 AND kdigo168_pos=1)),
   (5,'+ has T0 and temporal SA-AKI (sepsis first, <=7d)',
       (SELECT COUNT(*) FROM sic_cohort WHERE saaki_A=1),
       (SELECT COUNT(*) FROM sic_cohort WHERE saaki_C=1)),
   (6,'- exclude ESRD/transplant/pregnancy = final cohort',
       (SELECT COUNT(*) FROM sic_cohort WHERE saaki_A=1 AND excl_any=0),
       (SELECT COUNT(*) FROM sic_cohort WHERE saaki_C=1 AND excl_any=0))
) AS t(step,label,n_A,n_C) ORDER BY step;
-- Literature anchor: SICdb SA-AKI in Luo et al. CKJ 2026 = 1893 (Sepsis-3 + AKI + first + ICU >=24h)
-- 1	"All cases"	27350	27350
-- 2	"First ICU + adult >=18 + ICU stay >=24h"	14803	14803
-- 3	"+ sepsis (A=flag / C=sustained antibiotics)"	579	4528
-- 4	"+ official KDIGO positive"	443	3342
-- 5	"+ has T0 and temporal SA-AKI (sepsis first, <=7d)"	404	2741
-- 6	"- exclude ESRD/transplant/pregnancy = final cohort"	404	2737



-- Appendix 1: T0 source (how many cases the urine-output method added) -- C final cohort
SELECT
  COUNT(*) FILTER (WHERE saaki_C=1 AND excl_any=0) AS final_C,
  COUNT(*) FILTER (WHERE saaki_C=1 AND excl_any=0 AND t0_by_creat AND NOT t0_by_uo) AS t0_only_creat,
  COUNT(*) FILTER (WHERE saaki_C=1 AND excl_any=0 AND t0_by_uo AND NOT t0_by_creat) AS t0_only_uo,
  COUNT(*) FILTER (WHERE saaki_C=1 AND excl_any=0 AND t0_by_creat AND t0_by_uo)     AS t0_both
FROM sic_cohort;  --2737	139	1338	1260

-- Appendix 2: effect of tightening the sepsis criterion (share with sustained-antibiotic span >= 2 days; how many short-course / perioperative-prophylaxis cases are filtered out)
SELECT
  COUNT(*) FILTER (WHERE span_days IS NOT NULL) AS any_abx,
  COUNT(*) FILTER (WHERE sepsis_C=1)            AS sustained_abx_2d,
  ROUND(PERCENTILE_CONT(0.5) WITHIN GROUP (ORDER BY span_days) FILTER (WHERE span_days IS NOT NULL)::numeric,2) AS median_span_days
FROM sic_cohort;  --10748	4528	1.92

-- Appendix 3: baseline-source distribution (C final cohort)
SELECT baseline_src, COUNT(*) n,
       ROUND(100.0*COUNT(*)/SUM(COUNT(*)) OVER (),1) AS pct
FROM sic_cohort WHERE saaki_C=1 AND excl_any=0 GROUP BY baseline_src ORDER BY n DESC;
-- "icu_early24h"	1723	63.0
-- "pre_icu"	952	34.8
-- "mdrd75"	62	2.3






-- ############################################################################
-- SICdb | Part 2: exposure grouping (PPI vs. H2RA vs. none)
-- ----------------------------------------------------------------------------
-- Input: sic_cohort (finalized in Part 1; C-path main analysis saaki_C=1 AND excl_any=0; about 2737 cases)
-- Exposure definition (aligned with MIMIC/eICU):
--   group = class of the FIRST acid suppressant within the [T0, T0+24h] window (grace period 24h)
--   PPI : Esomeprazol(1426) / Pantoprazol(1427)
--   H2RA: Ranitidin(1599)
--   none: neither PPI nor H2RA within the window
--   both: both PPI and H2RA appear within the window (dropped / sensitivity)
--   newuser: no acid suppressant before T0 (off < T0) (SICdb has no pre-admission medication history -> in-ICU approximation)
--
-- [*] Administration form verified: acid suppressants are [discrete single-dose events] (bid/qd, each lasting 60s, 12/24h apart),
--   not the minute-level continuous expansion of antibiotics -> taking the first dose MIN(off) is correct, and COUNT is the true number of doses.
-- [*] Time unit: all Offset values are in seconds; window [t0, t0+86400].
--
-- Output: sic_exposure (sic_cohort + exposure_group / newuser / both, etc.)
-- ############################################################################

-- ============================================================================
-- Step 1: acid-suppressant administration events (restricted to the main-analysis cohort, with T0)
-- ============================================================================
DROP TABLE IF EXISTS sic_acid_events;
CREATE TABLE sic_acid_events AS
SELECT c.stay_id, c.t0,
       m."Offset" AS off,
       CASE WHEN m.drugid IN (1426,1427) THEN 'PPI'
            WHEN m.drugid = 1599         THEN 'H2RA' END AS acid_class
FROM sic_cohort c
JOIN medication m ON m.caseid = c.stay_id
WHERE c.saaki_C = 1 AND c.excl_any = 0           -- main-analysis cohort
  AND m.drugid IN (1426,1427,1599)
  AND m."Offset" >= 0;  --26597

-- ============================================================================
-- Step 2: within-window (=[T0, T0+24h]) exposure determination
--   first_in_window: class of the earliest acid suppressant within the window (ties broken by PPI/H2RA alphabetical order; affects only the very few at the same instant)
--   both_in_window : whether both PPI and H2RA occur within the window
-- ============================================================================
DROP TABLE IF EXISTS sic_exp_window;
CREATE TABLE sic_exp_window AS
WITH in_win AS (
    SELECT stay_id, off, acid_class
    FROM sic_acid_events
    WHERE off >= t0 AND off <= t0 + 86400
),
first_in_window AS (
    SELECT DISTINCT ON (stay_id) stay_id, acid_class AS first_acid_class, off AS first_acid_off
    FROM in_win
    ORDER BY stay_id, off ASC, acid_class ASC
),
both_flag AS (
    SELECT stay_id, CASE WHEN COUNT(DISTINCT acid_class) > 1 THEN 1 ELSE 0 END AS both_in_window
    FROM in_win GROUP BY stay_id
)
SELECT fw.stay_id, fw.first_acid_class, fw.first_acid_off, COALESCE(bf.both_in_window,0) AS both_in_window
FROM first_in_window fw
LEFT JOIN both_flag bf ON bf.stay_id = fw.stay_id;  --2422

-- ============================================================================
-- Step 3: pre-T0 acid suppressant (newuser determination)
-- ============================================================================
DROP TABLE IF EXISTS sic_prior_acid;
CREATE TABLE sic_prior_acid AS
SELECT DISTINCT stay_id, 1 AS prior_acid
FROM sic_acid_events
WHERE off < t0;  --1832

-- ============================================================================
-- Step 4: assemble sic_exposure
-- ============================================================================
DROP TABLE IF EXISTS sic_exposure;
CREATE TABLE sic_exposure AS
SELECT c.*,
   COALESCE(w.first_acid_class, 'none') AS exposure_group,
   w.first_acid_off,
   CASE WHEN w.first_acid_off IS NOT NULL
        THEN ROUND(((w.first_acid_off - c.t0)/3600.0)::numeric,1) END AS hours_to_exposure,
   COALESCE(w.both_in_window, 0) AS both_in_window,
   COALESCE(p.prior_acid, 0)     AS prior_acid,
   CASE WHEN COALESCE(p.prior_acid,0)=0 THEN 1 ELSE 0 END AS newuser
FROM sic_cohort c
LEFT JOIN sic_exp_window w ON w.stay_id = c.stay_id
LEFT JOIN sic_prior_acid p ON p.stay_id = c.stay_id
WHERE c.saaki_C = 1 AND c.excl_any = 0;          -- main-analysis cohort  --2737

-- ============================================================================
-- [*] Exposure grouping results
-- ============================================================================
-- (1) overall distribution of the three groups + new-user/both split
SELECT exposure_group,
       COUNT(*) AS n,
       ROUND(100.0*COUNT(*)/SUM(COUNT(*)) OVER (),1) AS pct,
       SUM(newuser)         AS n_newuser,
       SUM(both_in_window)  AS n_both,
       ROUND(AVG(hours_to_exposure)::numeric,1) AS mean_hrs_to_exposure
FROM sic_exposure
GROUP BY exposure_group ORDER BY n DESC;
--"PPI"	2417	88.3	671	4	10.2
--"none"	315	11.5	232	0	
--"H2RA"	5	0.2	2	0	7.3


-- (2) [*] head-to-head analyzable counts (excluding both): how many PPI vs. H2RA + the new-user subset
SELECT
  COUNT(*) FILTER (WHERE exposure_group='PPI'  AND both_in_window=0) AS n_ppi,
  COUNT(*) FILTER (WHERE exposure_group='H2RA' AND both_in_window=0) AS n_h2ra,
  COUNT(*) FILTER (WHERE exposure_group='none')                      AS n_none,
  COUNT(*) FILTER (WHERE exposure_group='PPI'  AND both_in_window=0 AND newuser=1) AS n_ppi_newuser,
  COUNT(*) FILTER (WHERE exposure_group='H2RA' AND both_in_window=0 AND newuser=1) AS n_h2ra_newuser
FROM sic_exposure;
--2413	5	315	670	2

-- (3) H2RA distribution by year (to see how much of the market-withdrawal effect falls within the cohort)
SELECT admissionyear,
       COUNT(*) FILTER (WHERE exposure_group='PPI'  AND both_in_window=0) AS ppi,
       COUNT(*) FILTER (WHERE exposure_group='H2RA' AND both_in_window=0) AS h2ra
FROM sic_exposure
GROUP BY admissionyear ORDER BY admissionyear;
-- 2013	268	0
-- 2014	261	3
-- 2015	262	0
-- 2016	296	2
-- 2017	289	0
-- 2018	286	0
-- 2019	269	0
-- 2020	270	0
-- 2021	212	0

