# #############################################################################
#  Proton-pump inhibitors vs. H2-receptor antagonists in sepsis-associated
#  acute kidney injury (SA-AKI): a target-trial emulation
#  ===========================================================================
#  MIMIC-IV cohort -- main analysis pipeline (R)
#
#  Design     : target-trial emulation; active-comparator, new-user design.
#  Exposure   : first in-ICU acid suppressant (PPI vs. H2RA) within [T0, T0+24h];
#               T0 = onset of SA-AKI.
#  Estimands  : (i) per-protocol = clone-censor-weight (CCW), grace period 1 day;
#               (ii) head-to-head IPTW = ITT-analog (co-primary / sensitivity).
#  Primary    : MAKE30 (major adverse kidney events at 30 d) and MAKE90.
#  Inference  : 20 multiply-imputed datasets, pooled by Rubin's rules.
#
#  Reproducibility
#    * Requires a PostgreSQL table `analytic_cohort` produced by the companion
#      SQL pipeline (see ppi_h2ra_saaki_pipeline.sql). MIMIC-IV requires
#      credentialed access (PhysioNet); NO patient-level data are in this repo.
#    * DB credentials are read from environment variables (see R0). Never commit
#      credentials to a public repository.
#    * Set seed is fixed; results are deterministic given the same input table
#      and package versions. Session info should be recorded with the release.
#
#  Outputs (all aggregate; safe to share): *.csv effect tables, Table 1,
#    and publication figures (forest / flow diagram / CIF / KM).
#
# #############################################################################

# =============================================================================
# PART 0 -- Environment, dependencies, global helpers (defined ONCE; reused)
# =============================================================================
Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
rm(list = ls())

library(DBI); library(RPostgres); library(dplyr); library(tidyr); library(tibble)
library(mice); library(WeightIt); library(cobalt)
library(tableone); library(survey); library(survival)
library(naniar); library(ggplot2)

## Propensity-score covariates (measured at or before T0; shared by all
## weighting and CCW modules). 56 baseline confounders.
ps_covars <- c(
  "age_at_icu","gender","race_cat","admission_type","first_careunit_cat","anchor_year_group",
  "weight_kg","baseline_egfr","first_aki_stage","ckd",
  "sofa_total","lactate","mech_vent","vasopressor","corticosteroid",
  "albumin","bilirubin_total","inr","platelet","wbc","sodium","potassium",
  "bicarbonate","anion_gap","magnesium","hemoglobin","fluid_balance_ml","uo_rate_ml_kg_h",
  "myocardial_infarct","congestive_heart_failure","peripheral_vascular_disease",
  "cerebrovascular_disease","dementia","chronic_pulmonary_disease","rheumatic_disease",
  "peptic_ulcer_disease","mild_liver_disease","severe_liver_disease",
  "diabetes_without_cc","diabetes_with_cc","paraplegia","renal_disease",
  "malignant_cancer","metastatic_solid_tumor","aids","elixhauser_score","hypertension",
  "gi_bleed_dx","rbc_transfusion_pre_t0","dnr_comfort",
  "nephrotoxin_any","nephro_vanco","nephro_acei_arb","nephro_nsaid",
  "enteral_nutrition","antiplatelet")

## Rubin's rules for a scalar estimate (pool on the given scale; exponentiate
## log-scale quantities afterwards).
pool_rubin <- function(est, v) {
  m <- length(est); qbar <- mean(est); ubar <- mean(v); b <- var(est)
  tvar <- ubar + (1 + 1/m) * b; se <- sqrt(tvar)
  df <- if (b > 0) (m - 1) * (1 + ubar / ((1 + 1/m) * b))^2 else 1e6
  c(est = qbar, se = se,
    lo = qbar - qt(.975, df) * se, hi = qbar + qt(.975, df) * se,
    p  = 2 * pt(-abs(qbar / se), df))
}
expci <- function(p) c(est = exp(unname(p["est"])), lo = exp(unname(p["lo"])),
                       hi = exp(unname(p["hi"])), p = unname(p["p"]))

## E-value (VanderWeele & Ding closed form; HR / sHR approximated as RR).
evalue_rr <- function(est, lo, hi) {
  f <- function(x) { x <- ifelse(x < 1, 1/x, x); x + sqrt(x*(x-1)) }
  ev_pt <- f(est)
  ev_ci <- if (est < 1) { if (hi >= 1) 1 else f(hi) } else { if (lo <= 1) 1 else f(lo) }
  c(round(ev_pt, 2), round(ev_ci, 2))
}

## Small utilities: effective sample size; factor-to-numeric.
ess  <- function(w) sum(w)^2 / sum(w^2)
numv <- function(x) if (is.numeric(x)) as.numeric(x) else as.numeric(as.integer(factor(x)) - 1)


# #############################################################################
#  R0 -- Load + type + range checks + missingness (no imputation yet)
# #############################################################################

## (1) Connect and read. Credentials come from environment variables so they
##     are never hard-coded in a public repo. Set e.g.:
##       Sys.setenv(PGHOST="localhost", PGDATABASE="mimiciv",
##                  PGUSER="postgres", PGPASSWORD="********")
con <- dbConnect(RPostgres::Postgres(),
                 host     = Sys.getenv("PGHOST",     "localhost"),
                 port     = as.integer(Sys.getenv("PGPORT", "5432")),
                 dbname   = Sys.getenv("PGDATABASE", "mimiciv"),
                 user     = Sys.getenv("PGUSER",     "postgres"),
                 password = Sys.getenv("PGPASSWORD"))
dat <- dbGetQuery(con, "SELECT * FROM analytic_cohort;")
dbDisconnect(con)
cat("Dimensions:", nrow(dat), "rows x", ncol(dat), "cols\n")   # expect 21744 x 121
cat("stay_id unique?", length(unique(dat$stay_id)) == nrow(dat), "\n")

## (2) Declare variable types -------------------------------------------------
cont_vars <- c(
  "age_at_icu","weight_kg",
  "creat_pre","creat_24h","creat_mdrd","baseline_creat","baseline_egfr",
  "lactate","albumin","bilirubin_total","inr","platelet","wbc",
  "sodium","potassium","bicarbonate","anion_gap","magnesium","hemoglobin",
  "sofa_total","sofa_resp","sofa_coag","sofa_liver","sofa_cardio","sofa_cns","sofa_renal",
  "fluid_balance_ml","uo_rate_ml_kg_h","elixhauser_score","charlson_comorbidity_index",
  "hours_to_exposure","icu_los_days","los_after_t0_days",
  "surv_28","surv_90","make_tte_time",
  "days_to_death","days_to_rrt","days_to_prog","days_to_recovery")
factor_multi <- c("gender","race","admission_type","admission_location","first_careunit",
                  "anchor_year_group","exposure_group","baseline_creat_source",
                  "first_acid_drug","first_aki_stage")
factor_bin <- c(
  "ckd","ckd_dx","newuser","both_in_window","exact_tie","prior_acid",
  "myocardial_infarct","congestive_heart_failure","peripheral_vascular_disease",
  "cerebrovascular_disease","dementia","chronic_pulmonary_disease","rheumatic_disease",
  "peptic_ulcer_disease","mild_liver_disease","severe_liver_disease",
  "diabetes_without_cc","diabetes_with_cc","paraplegia","renal_disease",
  "malignant_cancer","metastatic_solid_tumor","aids","hypertension",
  "gi_bleed_dx","rbc_transfusion_pre_t0","dnr_comfort",
  "nephro_vanco","nephro_aminoglyc","nephro_nsaid","nephro_acei_arb",
  "nephro_ampho","nephro_iv_contrast","nephrotoxin_any",
  "enteral_nutrition","antiplatelet","mech_vent","vasopressor","corticosteroid",
  "make30","make90","make_tte_event","death_28","death_30","death_90",
  "death_inhosp","death_icu","disch_before_d30","disch_before_d90",
  "rrt_init","rrt_dep_30","rrt_dep_90","rrt_pre_t0",
  "non_recovery_30","non_recovery_90","prog_eligible","prog_event",
  "recovery_eligible","recovery_event","nco_vte","nco_pressure")
dat <- dat %>% mutate(across(any_of(cont_vars), as.numeric),
                      across(any_of(c(factor_multi, factor_bin)), as.factor))
cat("\nColumns still numeric after typing (should be cont_vars + IDs):\n")
print(names(dat)[sapply(dat, is.numeric)])

## (2b) Timestamp consistency (each should be 0) ------------------------------
cat("\naki_time < sepsis_time:", sum(dat$aki_time  < dat$sepsis_time, na.rm = TRUE), "\n")
cat("dischtime  < t0        :", sum(dat$dischtime < dat$t0,          na.rm = TRUE), "\n")
cat("days_to_death < 0      :", sum(dat$days_to_death < 0,           na.rm = TRUE), "\n")

## (3) Physiologic range checks ----------------------------------------------
print(summary(dat[, intersect(cont_vars, names(dat))]))
range_check <- tibble::tribble(
  ~var,               ~lo,   ~hi,
  "lactate",           0,     50,
  "albumin",           0.5,   7,
  "sodium",            100,   180,
  "potassium",         1,     10,
  "bicarbonate",       2,     60,
  "hemoglobin",        2,     25,
  "platelet",          0,     2000,
  "uo_rate_ml_kg_h",   0,     30,
  "weight_kg",         20,    400)
oob <- range_check %>% rowwise() %>%
  mutate(n_out = sum(dat[[var]] < lo | dat[[var]] > hi, na.rm = TRUE)) %>% ungroup()
cat("\nOut-of-range counts (should be ~0; a few set to NA below):\n"); print(oob)

pdf("eda_continuous_distributions.pdf", width = 11, height = 8)
op <- par(mfrow = c(3, 3))
for (v in intersect(cont_vars, names(dat))) {
  x <- dat[[v]]; if (all(is.na(x))) next
  hist(x, main = v, xlab = "", breaks = 40, col = "grey80", border = "white")
}
par(op); dev.off()

## (4) Missingness ------------------------------------------------------------
miss_tbl <- data.frame(variable = names(dat), n_missing = colSums(is.na(dat)),
                       pct_missing = round(100 * colSums(is.na(dat)) / nrow(dat), 1)) %>% arrange(desc(pct_missing))
cat("\nVariables with missing values:\n");  print(miss_tbl[miss_tbl$n_missing > 0, ], row.names = FALSE)
cat("\nVariables with >40% missing:\n");     print(miss_tbl[miss_tbl$pct_missing > 40, ], row.names = FALSE)
miss_cols <- names(dat)[colSums(is.na(dat)) > 0]
pdf("eda_missingness.pdf", width = 11, height = 8)
print(gg_miss_var(dat[, miss_cols]))
print(vis_miss(dat[, miss_cols], warn_large_data = FALSE))
dev.off()

cat("\nExposure groups:\n"); print(table(dat$exposure_group))   # H2RA 3075 / neither 15230 / PPI 3439


# #############################################################################
#  R1 -- Multiple imputation (MICE, m = 20)
#   Only genuine-missing covariates are imputed; structural / intermediate
#   variables are not. SOFA total is a passive sum of its components.
#   Exposure + make90/death are predictors (for congeniality with the analysis).
# #############################################################################

## (1) Out-of-range -> NA (from `oob`: weight, urine-output rate) -------------
dat$weight_kg[dat$weight_kg < 20 | dat$weight_kg > 400]                 <- NA
dat$uo_rate_ml_kg_h[dat$uo_rate_ml_kg_h < 0 | dat$uo_rate_ml_kg_h > 30] <- NA
cat("After set-NA: weight missing", sum(is.na(dat$weight_kg)),
    "; uo_rate missing", sum(is.na(dat$uo_rate_ml_kg_h)), "\n")

## (2) Collapse high-cardinality categoricals --------------------------------
dat$race_cat <- factor(case_when(
  grepl("WHITE", dat$race)           ~ "White",
  grepl("BLACK", dat$race)           ~ "Black",
  grepl("HISPANIC|LATINO", dat$race) ~ "Hispanic",
  grepl("ASIAN", dat$race)           ~ "Asian",
  TRUE                               ~ "Other/Unknown"))
fc <- as.character(dat$first_careunit); rare <- names(which(table(fc) < 200))
dat$first_careunit_cat <- factor(ifelse(fc %in% rare, "Other", fc))

## (3) Imputation targets / predictors ---------------------------------------
vars_impute <- c(
  "weight_kg","albumin","lactate","bilirubin_total","inr","wbc",
  "magnesium","anion_gap","platelet",
  "sodium","potassium","bicarbonate","hemoglobin",
  "uo_rate_ml_kg_h","fluid_balance_ml","elixhauser_score",
  "sofa_resp","sofa_coag","sofa_liver","sofa_cardio","sofa_cns","sofa_renal","sofa_total")
vars_predict <- c(
  "age_at_icu","gender","race_cat","admission_type","first_careunit_cat","anchor_year_group",
  "baseline_creat","baseline_egfr","first_aki_stage","ckd","ckd_dx","prior_acid","newuser",
  "charlson_comorbidity_index","hypertension",
  "myocardial_infarct","congestive_heart_failure","peripheral_vascular_disease",
  "cerebrovascular_disease","dementia","chronic_pulmonary_disease","rheumatic_disease",
  "peptic_ulcer_disease","mild_liver_disease","severe_liver_disease",
  "diabetes_without_cc","diabetes_with_cc","paraplegia","renal_disease",
  "malignant_cancer","metastatic_solid_tumor","aids",
  "nephro_vanco","nephro_aminoglyc","nephro_nsaid","nephro_acei_arb",
  "nephro_ampho","nephro_iv_contrast","nephrotoxin_any",
  "mech_vent","vasopressor","corticosteroid",
  "gi_bleed_dx","rbc_transfusion_pre_t0","dnr_comfort","enteral_nutrition","antiplatelet",
  "exposure_group","make90","death_28","death_90")
imp_dat  <- dat[, c(vars_impute, vars_predict)]
stay_key <- dat$stay_id
cat("\nVariables being imputed:\n"); print(sort(names(which(colSums(is.na(imp_dat)) > 0))))

## (4) MICE method matrix: PMM for continuous; passive sum for SOFA total -----
ini  <- mice(imp_dat, maxit = 0, printFlag = FALSE)
meth <- ini$method; pred <- ini$predictorMatrix
meth[vars_predict] <- ""
meth[setdiff(vars_impute, "sofa_total")] <- "pmm"
meth["sofa_total"] <- "~ I(sofa_resp + sofa_coag + sofa_liver + sofa_cardio + sofa_cns + sofa_renal)"
pred["sofa_total", ] <- 0; pred[, "sofa_total"] <- 0
cat("\nNumber of variables per imputation method:\n"); print(table(meth[meth != ""]))

## (5) Run MICE (m = 20, maxit = 20). Auto-cached: delete imp_mice_m20.rds to
##     force recomputation if the INPUT TABLE changed. First run ~15-40 min.
if (file.exists("imp_mice_m20.rds")) {
  imp <- readRDS("imp_mice_m20.rds"); cat("Loaded cached MICE object.\n")
} else {
  set.seed(20240601)
  imp <- mice(imp_dat, m = 20, maxit = 20,
              method = meth, predictorMatrix = pred, printFlag = TRUE)
  saveRDS(imp, "imp_mice_m20.rds")
}

## (6) Assemble analysis-ready long table (stay_id + all outcome/other cols) --
long  <- complete(imp, "long", include = TRUE)
long$stay_id <- rep(stay_key, imp$m + 1)
long  <- left_join(long, dat[, c("stay_id", setdiff(names(dat), names(imp_dat)))], by = "stay_id")
saveRDS(long, "imp_long_analysis.rds")
cat("\nAnalysis-ready long table:", nrow(long), "x", ncol(long), "\n")

## (7) Imputation-quality diagnostics ----------------------------------------
pdf("mice_convergence.pdf", width = 11, height = 8); plot(imp); dev.off()
pdf("mice_density.pdf", width = 11, height = 7)
densityplot(imp, ~ albumin + lactate + bilirubin_total + inr + uo_rate_ml_kg_h + sofa_renal)
dev.off()
diag_one <- function(imp, var) {
  x <- imp$data[[var]]; mi <- is.na(x); if (sum(mi) == 0) return(NULL)
  iv <- unlist(lapply(1:imp$m, function(i) complete(imp, i)[[var]][mi]))
  data.frame(variable = var, n_missing = sum(mi),
             obs_median = round(median(x[!mi]), 2), imp_median = round(median(iv), 2),
             obs_mean   = round(mean(x[!mi]),   2), imp_mean   = round(mean(iv),   2))
}
imp_diag <- do.call(rbind, lapply(setdiff(vars_impute, "sofa_total"), function(v) diag_one(imp, v)))
cat("\nImputation diagnostics (observed vs imputed; should be close):\n")
print(imp_diag, row.names = FALSE)


# #############################################################################
#  R2 -- IPTW (propensity-score weighting) + balance, PPI vs H2RA head-to-head
# #############################################################################

long <- readRDS("imp_long_analysis.rds")
ps_formula <- reformulate(ps_covars, "treat")

## Sanity check: PS covariates must have no NA after imputation --------------
na_chk <- colSums(is.na(long[long$.imp >= 1, ps_covars, drop = FALSE]))
if (any(na_chk > 0)) { cat("WARNING: PS covariates still NA after imputation:\n"); print(na_chk[na_chk > 0]) } else
  cat("PS covariates complete after imputation.\n")

## Per imputation: fit PS -> stabilized weights -> trim at 1st/99th pct -------
m <- 20
hh_list <- vector("list", m); bal_list <- vector("list", m); wt_diag <- data.frame()
for (i in 1:m) {
  di <- long %>% filter(.imp == i, exposure_group %in% c("PPI","H2RA")) %>%
    mutate(treat = as.integer(exposure_group == "PPI"))
  W  <- weightit(ps_formula, data = di, method = "glm",   # older WeightIt: method = "ps"
                 estimand = "ATE", stabilize = TRUE)
  Wt <- trim(W, at = 0.99, lower = TRUE); di$ipw <- Wt$weights
  hh_list[[i]]  <- di
  bal_list[[i]] <- bal.tab(Wt, un = TRUE, stats = "mean.diffs", binary = "std", continuous = "std")
  w <- di$ipw
  wt_diag <- rbind(wt_diag, data.frame(
    imp = i, w_max = round(max(w), 2), w_p99 = round(quantile(w, .99), 2),
    n_ppi = sum(di$treat == 1), ess_ppi = round(ess(w[di$treat == 1])),
    n_h2  = sum(di$treat == 0), ess_h2  = round(ess(w[di$treat == 0]))))
}

## Diagnostic 1: weights / effective sample size (mean over 20) --------------
cat("\nWeights / effective sample size (mean of 20):\n"); print(round(colMeans(wt_diag[, -1]), 1))

## Diagnostic 2: balance table (|SMD| before/after, mean over 20) ------------
covs    <- rownames(bal_list[[1]]$Balance)
un_mat  <- sapply(bal_list, function(b) abs(b$Balance$Diff.Un))
adj_mat <- sapply(bal_list, function(b) abs(b$Balance$Diff.Adj))
bal_summary <- data.frame(covariate = covs,
                          smd_unweighted = round(rowMeans(un_mat), 3),
                          smd_weighted   = round(rowMeans(adj_mat), 3)) %>%
  filter(covariate != "prop.score") %>% arrange(desc(smd_weighted))
cat("\n# covariates with weighted SMD > 0.1:", sum(bal_summary$smd_weighted > 0.1), "(ideal 0)\n")
cat("Max weighted SMD:", max(bal_summary$smd_weighted), "\n\n")
print(bal_summary, row.names = FALSE)

## Diagnostic 3: Love plot ---------------------------------------------------
lp <- pivot_longer(bal_summary, c(smd_unweighted, smd_weighted), names_to = "type", values_to = "smd") %>%
  mutate(type = recode(type, smd_unweighted = "Before weighting", smd_weighted = "After weighting"))
ggplot(lp, aes(smd, reorder(covariate, smd), color = type)) +
  geom_point(size = 2) + geom_vline(xintercept = 0.1, linetype = "dashed") +
  labs(x = "|Standardized mean difference|", y = NULL, color = NULL,
       title = "Covariate balance before/after IPTW (mean of 20 imputations)") +
  theme_minimal(base_size = 11)
ggsave("iptw_loveplot.pdf", width = 8, height = 11)

hh_all <- bind_rows(hh_list)
saveRDS(hh_all, "iptw_weighted_long.rds")
cat("\nWeighted head-to-head long table:", nrow(hh_all), "rows\n")


# #############################################################################
#  R3 -- Table 1 (baseline characteristics; unweighted + IPTW-weighted, SMD)
# #############################################################################

hh_all <- readRDS("iptw_weighted_long.rds")
d1 <- hh_all %>% filter(.imp == 1) %>% mutate(exposure_group = droplevels(factor(exposure_group)))

vars <- c(
  "age_at_icu","gender","race_cat","admission_type","first_careunit_cat","anchor_year_group",
  "weight_kg","baseline_creat","baseline_egfr","first_aki_stage","ckd",
  "sofa_total","lactate","mech_vent","vasopressor","corticosteroid",
  "albumin","bilirubin_total","inr","platelet","wbc","sodium","potassium",
  "bicarbonate","anion_gap","magnesium","hemoglobin","fluid_balance_ml","uo_rate_ml_kg_h",
  "charlson_comorbidity_index","elixhauser_score","hypertension",
  "myocardial_infarct","congestive_heart_failure","peripheral_vascular_disease",
  "cerebrovascular_disease","chronic_pulmonary_disease","mild_liver_disease",
  "severe_liver_disease","diabetes_without_cc","diabetes_with_cc","renal_disease",
  "malignant_cancer","metastatic_solid_tumor",
  "gi_bleed_dx","rbc_transfusion_pre_t0","dnr_comfort",
  "nephrotoxin_any","nephro_vanco","nephro_acei_arb","nephro_nsaid",
  "enteral_nutrition","antiplatelet")
nonnormal <- c("lactate","bilirubin_total","inr","fluid_balance_ml","uo_rate_ml_kg_h",
               "baseline_creat","sofa_total","charlson_comorbidity_index","elixhauser_score",
               "platelet","wbc")
factor_vars <- c("gender","race_cat","admission_type","first_careunit_cat","anchor_year_group",
                 "first_aki_stage","ckd","mech_vent","vasopressor","corticosteroid","hypertension",
                 "myocardial_infarct","congestive_heart_failure","peripheral_vascular_disease",
                 "cerebrovascular_disease","chronic_pulmonary_disease","mild_liver_disease",
                 "severe_liver_disease","diabetes_without_cc","diabetes_with_cc","renal_disease",
                 "malignant_cancer","metastatic_solid_tumor",
                 "gi_bleed_dx","rbc_transfusion_pre_t0","dnr_comfort","nephrotoxin_any",
                 "nephro_vanco","nephro_acei_arb","nephro_nsaid","enteral_nutrition","antiplatelet")

## Unweighted ----------------------------------------------------------------
t1_unw <- CreateTableOne(vars = vars, strata = "exposure_group", data = d1,
                         factorVars = factor_vars, addOverall = TRUE, test = FALSE)
cat("\n========== Table 1 (unweighted) ==========\n")
print(t1_unw, nonnormal = nonnormal, smd = TRUE, showAllLevels = TRUE)

## IPTW-weighted -------------------------------------------------------------
dsgn <- svydesign(ids = ~1, weights = ~ipw, data = d1)
t1_w <- svyCreateTableOne(vars = vars, strata = "exposure_group", data = dsgn,
                          factorVars = factor_vars, test = FALSE)
cat("\n========== Table 1 (IPTW-weighted) ==========\n")
print(t1_w, nonnormal = nonnormal, smd = TRUE, showAllLevels = TRUE)

write.csv(print(t1_unw, nonnormal = nonnormal, smd = TRUE, showAllLevels = TRUE, printToggle = FALSE),
          "table1_unweighted.csv")
write.csv(print(t1_w,  nonnormal = nonnormal, smd = TRUE, showAllLevels = TRUE, printToggle = FALSE),
          "table1_weighted.csv")
cat("\nTable 1 written (unweighted + weighted).\n")



# #############################################################################
#  R4 -- Primary outcome causal effects (IPTW / ITT-analog)
#    Outcomes : MAKE30 (primary), MAKE90 (co-primary)
#    Estimands: risk difference (RD) + risk ratio (RR) [binary] + restricted
#               mean survival time (RMST) difference (event-free days)
#    Per imputation: svyglm with robust SE -> 20 datasets pooled by Rubin's rules
#  Output is an AGGREGATE effect table -> safe to share (primary_effects_iptw.csv)
# #############################################################################

hh_all <- readRDS("iptw_weighted_long.rds")

## Single outcome x single horizon tau: one estimate per imputation -> pool.
##  With a fixed horizon, fully observed death, and in-hospital RRT, there is no
##  loss to follow-up before tau, so RMST(tau) = E[min(T, tau)] and a weighted
##  linear model estimates the between-group difference (robust SE). Out-of-
##  hospital RRT is unobserved (same limitation as the binary MAKE endpoint).
analyze_make <- function(dat, make_var, tau) {
  m <- max(dat$.imp); R <- data.frame()
  for (i in 1:m) {
    di <- dat %>% filter(.imp == i, exposure_group %in% c("PPI","H2RA")) %>%
      mutate(treat = as.integer(exposure_group == "PPI"),
             y  = as.integer(as.character(.data[[make_var]])),           # 0/1
             ev = as.integer(as.character(make_tte_event)),              # death or RRT
             tt = pmin(ifelse(ev == 1, make_tte_time, tau), tau)) %>%    # min(event time, tau)
      filter(!is.na(y))                                                  # drop outcome NA (9 for make30)
    dz <- svydesign(ids = ~1, weights = ~ipw, data = di)
    f_rd <- svyglm(y  ~ treat, design = dz, family = gaussian())                 # linear prob -> RD
    f_rr <- svyglm(y  ~ treat, design = dz, family = quasipoisson(link = "log")) # robust Poisson -> logRR
    f_rm <- svyglm(tt ~ treat, design = dz, family = gaussian())                 # RMST difference
    w0 <- di$treat == 0; w1 <- di$treat == 1
    R <- rbind(R, data.frame(
      rd  = coef(f_rd)[["treat"]],  v_rd  = vcov(f_rd)["treat","treat"],
      lrr = coef(f_rr)[["treat"]],  v_lrr = vcov(f_rr)["treat","treat"],
      rm  = coef(f_rm)[["treat"]],  v_rm  = vcov(f_rm)["treat","treat"],
      risk_h2  = weighted.mean(di$y[w0],  di$ipw[w0]),
      risk_ppi = weighted.mean(di$y[w1],  di$ipw[w1]),
      rmst_h2  = weighted.mean(di$tt[w0], di$ipw[w0]),
      rmst_ppi = weighted.mean(di$tt[w1], di$ipw[w1])))
  }
  rd <- pool_rubin(R$rd, R$v_rd)
  rr <- pool_rubin(R$lrr, R$v_lrr)        # pool on the log scale, then exponentiate
  rm <- pool_rubin(R$rm, R$v_rm)
  list(
    n_ppi = sum(dat$.imp == 1 & dat$exposure_group == "PPI"),
    n_h2  = sum(dat$.imp == 1 & dat$exposure_group == "H2RA"),
    risk_ppi = mean(R$risk_ppi), risk_h2 = mean(R$risk_h2),
    rmst_ppi = mean(R$rmst_ppi), rmst_h2 = mean(R$rmst_h2),
    RD = rd,
    RR = c(est = exp(unname(rr["est"])), lo = exp(unname(rr["lo"])),
           hi = exp(unname(rr["hi"])), p = unname(rr["p"])),
    RMST = rm)
}

res30 <- analyze_make(hh_all, "make30", 30)
res90 <- analyze_make(hh_all, "make90", 90)

## Print --------------------------------------------------------------------
show_res <- function(r, lab) {
  cat("\n==========", lab, "==========\n")
  cat(sprintf("N (head-to-head): PPI %d  vs  H2RA %d\n", r$n_ppi, r$n_h2))
  cat(sprintf("Weighted risk : PPI %.1f%%   H2RA %.1f%%\n", 100*r$risk_ppi, 100*r$risk_h2))
  cat(sprintf("Risk diff RD  : %+.1f%%  (95%%CI %+.1f%% , %+.1f%%)   p=%.3f\n",
              100*unname(r$RD["est"]), 100*unname(r$RD["lo"]), 100*unname(r$RD["hi"]), unname(r$RD["p"])))
  cat(sprintf("Risk ratio RR : %.2f   (95%%CI %.2f , %.2f)\n",
              unname(r$RR["est"]), unname(r$RR["lo"]), unname(r$RR["hi"])))
  cat(sprintf("RMST (days)   : PPI %.2f   H2RA %.2f\n", r$rmst_ppi, r$rmst_h2))
  cat(sprintf("RMST diff     : %+.2f days (95%%CI %+.2f , %+.2f)   p=%.3f\n",
              unname(r$RMST["est"]), unname(r$RMST["lo"]), unname(r$RMST["hi"]), unname(r$RMST["p"])))
  cat("(Interpretation: negative RD/RMST diff = better with PPI; RR<1 = lower risk with PPI)\n")
}
show_res(res30, "MAKE30 (primary)")
show_res(res90, "MAKE90 (co-primary)")

## Export aggregate table ---------------------------------------------------
to_row <- function(r, lab) data.frame(
  outcome = lab, n_ppi = r$n_ppi, n_h2 = r$n_h2,
  risk_ppi_pct = round(100*r$risk_ppi, 1), risk_h2_pct = round(100*r$risk_h2, 1),
  RD_pct = round(100*unname(r$RD["est"]), 2),
  RD_lo  = round(100*unname(r$RD["lo"]), 2), RD_hi = round(100*unname(r$RD["hi"]), 2),
  RD_p   = signif(unname(r$RD["p"]), 3),
  RR     = round(unname(r$RR["est"]), 3),
  RR_lo  = round(unname(r$RR["lo"]), 3), RR_hi = round(unname(r$RR["hi"]), 3),
  RMSTd_day = round(unname(r$RMST["est"]), 3),
  RMSTd_lo  = round(unname(r$RMST["lo"]), 3), RMSTd_hi = round(unname(r$RMST["hi"]), 3),
  RMSTd_p   = signif(unname(r$RMST["p"]), 3),
  row.names = NULL)
res_tbl <- rbind(to_row(res30, "MAKE30"), to_row(res90, "MAKE90"))
write.csv(res_tbl, "primary_effects_iptw.csv", row.names = FALSE)
cat("\n---- Summary (written to primary_effects_iptw.csv) ----\n"); print(res_tbl, row.names = FALSE)
# Checks: does RD 95%CI include 0, RR include 1, RMST diff include 0?
# MAKE30 and MAKE90 are co-primary (pre-specified); no multiplicity adjustment.



# #############################################################################
#  R5 -- Secondary outcomes
#    Batch 1: (a) 28/90-day mortality -- RD + RR + RMST (survival days)
#             (b) RRT initiation (death as competing risk) -- cause-specific HR
#                 + Fine-Gray sHR + weighted cumulative incidence
#    Batch 2: (c) renal recovery and (d) AKI progression to KDIGO 3, both as
#             in-hospital competing-risk endpoints (death competes), restricted
#             to their eligible subsets.
#    IPTW-weighted; 20 imputations pooled by Rubin's rules.
#  Outputs are AGGREGATE effect tables -> safe to share.
# #############################################################################

hh_all <- readRDS("iptw_weighted_long.rds")

# =============================================================================
#  Batch 1 (a) -- Mortality (28 / 90 days)
# =============================================================================
analyze_mortality <- function(dat, death_var, surv_var) {
  m <- max(dat$.imp); R <- data.frame()
  for (i in 1:m) {
    di <- dat %>% filter(.imp == i, exposure_group %in% c("PPI","H2RA")) %>%
      mutate(treat = as.integer(exposure_group == "PPI"),
             y = as.integer(as.character(.data[[death_var]])),     # 0/1 death
             s = as.numeric(.data[[surv_var]]))                    # min(time to death, tau)
    dz <- svydesign(ids = ~1, weights = ~ipw, data = di)
    f_rd <- svyglm(y ~ treat, design = dz, family = gaussian())
    f_rr <- svyglm(y ~ treat, design = dz, family = quasipoisson(link = "log"))
    f_rm <- svyglm(s ~ treat, design = dz, family = gaussian())
    w0 <- di$treat == 0; w1 <- di$treat == 1
    R <- rbind(R, data.frame(
      rd = coef(f_rd)[["treat"]], v_rd = vcov(f_rd)["treat","treat"],
      lrr = coef(f_rr)[["treat"]], v_lrr = vcov(f_rr)["treat","treat"],
      rm = coef(f_rm)[["treat"]], v_rm = vcov(f_rm)["treat","treat"],
      risk_h2 = weighted.mean(di$y[w0], di$ipw[w0]), risk_ppi = weighted.mean(di$y[w1], di$ipw[w1]),
      rmst_h2 = weighted.mean(di$s[w0], di$ipw[w0]), rmst_ppi = weighted.mean(di$s[w1], di$ipw[w1])))
  }
  rd <- pool_rubin(R$rd, R$v_rd); rr <- pool_rubin(R$lrr, R$v_lrr); rm <- pool_rubin(R$rm, R$v_rm)
  list(n_ppi = sum(dat$.imp==1 & dat$exposure_group=="PPI"),
       n_h2  = sum(dat$.imp==1 & dat$exposure_group=="H2RA"),
       risk_ppi = mean(R$risk_ppi), risk_h2 = mean(R$risk_h2),
       rmst_ppi = mean(R$rmst_ppi), rmst_h2 = mean(R$rmst_h2),
       RD = rd, RR = expci(rr), RMST = rm)
}

res_d28 <- analyze_mortality(hh_all, "death_28", "surv_28")
res_d90 <- analyze_mortality(hh_all, "death_90", "surv_90")

show_mort <- function(r, lab) {
  cat("\n==========", lab, "==========\n")
  cat(sprintf("N: PPI %d vs H2RA %d\n", r$n_ppi, r$n_h2))
  cat(sprintf("Weighted mortality: PPI %.1f%%   H2RA %.1f%%\n", 100*r$risk_ppi, 100*r$risk_h2))
  cat(sprintf("Risk diff RD : %+.1f%% (95%%CI %+.1f%% , %+.1f%%)  p=%.3f   [>0 = more deaths with PPI]\n",
              100*unname(r$RD["est"]), 100*unname(r$RD["lo"]), 100*unname(r$RD["hi"]), unname(r$RD["p"])))
  cat(sprintf("Risk ratio RR : %.2f (95%%CI %.2f , %.2f)            [>1 = worse with PPI]\n",
              unname(r$RR["est"]), unname(r$RR["lo"]), unname(r$RR["hi"])))
  cat(sprintf("Survival days : PPI %.2f   H2RA %.2f\n", r$rmst_ppi, r$rmst_h2))
  cat(sprintf("RMST diff    : %+.2f days (95%%CI %+.2f , %+.2f)  p=%.3f   [>0 = PPI survives longer = better]\n",
              unname(r$RMST["est"]), unname(r$RMST["lo"]), unname(r$RMST["hi"]), unname(r$RMST["p"])))
}
show_mort(res_d28, "28-day mortality")
show_mort(res_d90, "90-day mortality")

# =============================================================================
#  Batch 1 (b) -- RRT initiation (death as competing risk)
# =============================================================================
analyze_rrt <- function(dat, tau) {
  m <- max(dat$.imp); csh <- data.frame(); shr <- data.frame(); cif <- data.frame()
  cif_arm <- function(sub) {
    sf <- survfit(Surv(ftime, ef) ~ 1, data = sub, weights = sub$ipw)
    ss <- summary(sf, times = tau, extend = TRUE)
    j <- which(colnames(ss$pstate) == "RRT"); as.numeric(ss$pstate[1, j])
  }
  for (i in 1:m) {
    di <- dat %>% filter(.imp == i, exposure_group %in% c("PPI","H2RA"),
                         as.integer(as.character(rrt_pre_t0)) == 0) %>%       # exclude pre-T0 dialysis
      mutate(treat = as.integer(exposure_group == "PPI"),
             rinit = as.integer(as.character(rrt_init)),
             t_rrt = ifelse(rinit == 1, days_to_rrt, Inf),
             t_dth = ifelse(!is.na(days_to_death), days_to_death, Inf),
             fe    = pmin(t_rrt, t_dth),
             status = ifelse(fe > tau, 0L, ifelse(t_rrt <= t_dth, 1L, 2L)),    # 0 cens / 1 RRT / 2 death
             ftime  = pmin(fe, tau),
             ef = factor(status, levels = c(0,1,2), labels = c("cens","RRT","death")),
             rid = row_number())                                          # patient id (cluster-robust SE)
    # cause-specific HR (death censored)
    cox <- coxph(Surv(ftime, status == 1) ~ treat, data = di, weights = ipw, robust = TRUE)
    csh <- rbind(csh, data.frame(b = coef(cox)[["treat"]], v = vcov(cox)["treat","treat"]))
    # Fine-Gray sHR (death kept in risk set; weight = FG weight x IPTW; id = rid for robust SE)
    fg <- finegray(Surv(ftime, ef) ~ treat + ipw + rid, data = di, etype = "RRT")
    fgc <- coxph(Surv(fgstart, fgstop, fgstatus) ~ treat, data = fg,
                 weights = fgwt * ipw, id = rid, robust = TRUE)
    shr <- rbind(shr, data.frame(b = coef(fgc)[["treat"]], v = vcov(fgc)["treat","treat"]))
    # weighted cumulative incidence @ tau
    cif <- rbind(cif, data.frame(cif_ppi = cif_arm(di[di$treat==1,]),
                                 cif_h2  = cif_arm(di[di$treat==0,])))
  }
  cs <- pool_rubin(csh$b, csh$v); sh <- pool_rubin(shr$b, shr$v)
  keep <- dat$.imp==1 & as.integer(as.character(dat$rrt_pre_t0))==0
  list(n_ppi = sum(keep & dat$exposure_group=="PPI"),
       n_h2  = sum(keep & dat$exposure_group=="H2RA"),
       cif_ppi = mean(cif$cif_ppi), cif_h2 = mean(cif$cif_h2),
       csHR = expci(cs), sHR = expci(sh))
}

res_rrt30 <- analyze_rrt(hh_all, 30)
res_rrt90 <- analyze_rrt(hh_all, 90)

show_rrt <- function(r, lab) {
  cat("\n==========", lab, "==========\n")
  cat(sprintf("N (excl. pre-T0 dialysis): PPI %d vs H2RA %d\n", r$n_ppi, r$n_h2))
  cat(sprintf("Weighted cumulative incidence @tau: PPI %.1f%%   H2RA %.1f%%\n", 100*r$cif_ppi, 100*r$cif_h2))
  cat(sprintf("Cause-specific HR : %.2f (95%%CI %.2f , %.2f)  p=%.3f   [>1 = higher RRT rate with PPI]\n",
              unname(r$csHR["est"]), unname(r$csHR["lo"]), unname(r$csHR["hi"]), unname(r$csHR["p"])))
  cat(sprintf("Fine-Gray sHR     : %.2f (95%%CI %.2f , %.2f)  p=%.3f   [>1 = more cumulative RRT with PPI]\n",
              unname(r$sHR["est"]), unname(r$sHR["lo"]), unname(r$sHR["hi"]), unname(r$sHR["p"])))
}
show_rrt(res_rrt30, "RRT initiation @30 days")
show_rrt(res_rrt90, "RRT initiation @90 days")

## Export Batch 1 -----------------------------------------------------------
mort_row <- function(r, lab) data.frame(
  outcome = lab, n_ppi = r$n_ppi, n_h2 = r$n_h2,
  ppi_pct = round(100*r$risk_ppi,1), h2_pct = round(100*r$risk_h2,1),
  RD_pct = round(100*unname(r$RD["est"]),2), RD_lo = round(100*unname(r$RD["lo"]),2),
  RD_hi = round(100*unname(r$RD["hi"]),2), RD_p = signif(unname(r$RD["p"]),3),
  RR = round(unname(r$RR["est"]),3), RR_lo = round(unname(r$RR["lo"]),3), RR_hi = round(unname(r$RR["hi"]),3),
  RMSTd = round(unname(r$RMST["est"]),3), RMSTd_lo = round(unname(r$RMST["lo"]),3),
  RMSTd_hi = round(unname(r$RMST["hi"]),3), RMSTd_p = signif(unname(r$RMST["p"]),3), row.names = NULL)
rrt_row <- function(r, lab) data.frame(
  outcome = lab, n_ppi = r$n_ppi, n_h2 = r$n_h2,
  ppi_cif_pct = round(100*r$cif_ppi,1), h2_cif_pct = round(100*r$cif_h2,1),
  csHR = round(unname(r$csHR["est"]),3), csHR_lo = round(unname(r$csHR["lo"]),3),
  csHR_hi = round(unname(r$csHR["hi"]),3), csHR_p = signif(unname(r$csHR["p"]),3),
  sHR = round(unname(r$sHR["est"]),3), sHR_lo = round(unname(r$sHR["lo"]),3),
  sHR_hi = round(unname(r$sHR["hi"]),3), sHR_p = signif(unname(r$sHR["p"]),3), row.names = NULL)

cat("\n---- Mortality summary ----\n")
print(rbind(mort_row(res_d28,"death_28"), mort_row(res_d90,"death_90")), row.names = FALSE)
cat("\n---- RRT summary ----\n")
print(rbind(rrt_row(res_rrt30,"RRT_30"), rrt_row(res_rrt90,"RRT_90")), row.names = FALSE)
write.csv(rbind(mort_row(res_d28,"death_28"), mort_row(res_d90,"death_90")),
          "secondary_mortality.csv", row.names = FALSE)
write.csv(rbind(rrt_row(res_rrt30,"RRT_30"), rrt_row(res_rrt90,"RRT_90")),
          "secondary_rrt.csv", row.names = FALSE)

# =============================================================================
#  Batch 2 -- Generic competing-risk analysis (in-hospital endpoint; death
#  competes; restricted to an eligible subset). elig / evvar / tvar are the
#  column names (as strings) for eligibility / event / event-time.
# =============================================================================
analyze_cr <- function(dat, elig, evvar, tvar, tau) {
  m <- max(dat$.imp); csh <- data.frame(); shr <- data.frame(); cif <- data.frame()
  cif_arm <- function(sub) {
    sf <- survfit(Surv(ftime, ef) ~ 1, data = sub, weights = sub$ipw)
    ss <- summary(sf, times = tau, extend = TRUE)
    j <- which(colnames(ss$pstate) == "event"); as.numeric(ss$pstate[1, j])
  }
  for (i in 1:m) {
    di <- dat %>% filter(.imp == i, exposure_group %in% c("PPI","H2RA"),
                         as.integer(as.character(.data[[elig]])) == 1) %>%   # eligible only
      mutate(treat = as.integer(exposure_group == "PPI"),
             ev1   = as.integer(as.character(.data[[evvar]])),               # target event 0/1
             t_ev  = ifelse(ev1 == 1, .data[[tvar]], Inf),
             t_dth = ifelse(as.integer(as.character(death_inhosp)) == 1, days_to_death, Inf), # in-hosp death competes
             ctime = pmin(los_after_t0_days, tau),                           # censor at discharge or tau
             fe    = pmin(t_ev, t_dth),
             status = ifelse(fe > ctime, 0L, ifelse(t_ev <= t_dth, 1L, 2L)),  # 0 cens / 1 event / 2 death
             ftime  = pmin(fe, ctime),
             ef = factor(status, levels = c(0,1,2), labels = c("cens","event","death")),
             rid = row_number())
    # cause-specific HR (death censored)
    cox <- coxph(Surv(ftime, status == 1) ~ treat, data = di, weights = ipw, robust = TRUE)
    csh <- rbind(csh, data.frame(b = coef(cox)[["treat"]], v = vcov(cox)["treat","treat"]))
    # Fine-Gray sHR (death kept in risk set; weight = FG x IPTW; id for cluster-robust SE)
    fg <- finegray(Surv(ftime, ef) ~ treat + ipw + rid, data = di, etype = "event")
    fgc <- coxph(Surv(fgstart, fgstop, fgstatus) ~ treat, data = fg,
                 weights = fgwt * ipw, id = rid, robust = TRUE)
    shr <- rbind(shr, data.frame(b = coef(fgc)[["treat"]], v = vcov(fgc)["treat","treat"]))
    # weighted cumulative incidence @ tau
    cif <- rbind(cif, data.frame(cif_ppi = cif_arm(di[di$treat==1,]),
                                 cif_h2  = cif_arm(di[di$treat==0,])))
  }
  cs <- pool_rubin(csh$b, csh$v); sh <- pool_rubin(shr$b, shr$v)
  keep <- dat$.imp==1 & as.integer(as.character(dat[[elig]]))==1
  list(n_ppi = sum(keep & dat$exposure_group=="PPI"),
       n_h2  = sum(keep & dat$exposure_group=="H2RA"),
       cif_ppi = mean(cif$cif_ppi), cif_h2 = mean(cif$cif_h2),
       csHR = expci(cs), sHR = expci(sh))
}

res_rec30  <- analyze_cr(hh_all, "recovery_eligible", "recovery_event", "days_to_recovery", 30)
res_rec90  <- analyze_cr(hh_all, "recovery_eligible", "recovery_event", "days_to_recovery", 90)
res_prog30 <- analyze_cr(hh_all, "prog_eligible",     "prog_event",     "days_to_prog",     30)
res_prog90 <- analyze_cr(hh_all, "prog_eligible",     "prog_event",     "days_to_prog",     90)

## Print (good = TRUE -> HR > 1 means "better") ------------------------------
show_cr <- function(r, lab, good) {
  dir_cs <- if (good) "[>1 = more/faster recovery with PPI = better]" else "[>1 = more progression with PPI = worse]"
  cat("\n==========", lab, "==========\n")
  cat(sprintf("N (eligible): PPI %d vs H2RA %d\n", r$n_ppi, r$n_h2))
  cat(sprintf("Weighted cumulative incidence @tau: PPI %.1f%%   H2RA %.1f%%\n", 100*r$cif_ppi, 100*r$cif_h2))
  cat(sprintf("Cause-specific HR : %.2f (95%%CI %.2f , %.2f)  p=%.3f   %s\n",
              unname(r$csHR["est"]), unname(r$csHR["lo"]), unname(r$csHR["hi"]), unname(r$csHR["p"]), dir_cs))
  cat(sprintf("Fine-Gray sHR     : %.2f (95%%CI %.2f , %.2f)  p=%.3f\n",
              unname(r$sHR["est"]), unname(r$sHR["lo"]), unname(r$sHR["hi"]), unname(r$sHR["p"])))
}
show_cr(res_rec30,  "Renal recovery @30 days",         TRUE)
show_cr(res_rec90,  "Renal recovery @90 days",         TRUE)
show_cr(res_prog30, "AKI progression to KDIGO3 @30 days", FALSE)
show_cr(res_prog90, "AKI progression to KDIGO3 @90 days", FALSE)

## Export Batch 2 -----------------------------------------------------------
cr_row <- function(r, lab) data.frame(
  outcome = lab, n_ppi = r$n_ppi, n_h2 = r$n_h2,
  ppi_cif_pct = round(100*r$cif_ppi,1), h2_cif_pct = round(100*r$cif_h2,1),
  csHR = round(unname(r$csHR["est"]),3), csHR_lo = round(unname(r$csHR["lo"]),3),
  csHR_hi = round(unname(r$csHR["hi"]),3), csHR_p = signif(unname(r$csHR["p"]),3),
  sHR = round(unname(r$sHR["est"]),3), sHR_lo = round(unname(r$sHR["lo"]),3),
  sHR_hi = round(unname(r$sHR["hi"]),3), sHR_p = signif(unname(r$sHR["p"]),3), row.names = NULL)

cat("\n---- Renal recovery summary ----\n")
print(rbind(cr_row(res_rec30,"recovery_30"), cr_row(res_rec90,"recovery_90")), row.names = FALSE)
cat("\n---- AKI progression summary ----\n")
print(rbind(cr_row(res_prog30,"progression_30"), cr_row(res_prog90,"progression_90")), row.names = FALSE)
write.csv(rbind(cr_row(res_rec30,"recovery_30"), cr_row(res_rec90,"recovery_90")),
          "secondary_recovery.csv", row.names = FALSE)
write.csv(rbind(cr_row(res_prog30,"progression_30"), cr_row(res_prog90,"progression_90")),
          "secondary_progression.csv", row.names = FALSE)
# Checks: do recovery/progression HRs depart from 1? Are csHR and sHR consistent?
# Together with mortality and RRT, all three MAKE components (death / RRT /
# non-recovery) have now been examined separately.




# #############################################################################
#  R6 -- Robustness: negative-control outcomes (NCO) + E-values
#    (1) NCO: incident VTE (primary) + stage>=2 pressure injury (secondary),
#        same IPTW analysis -> should be null if confounding is well controlled.
#    (2) E-value: minimum strength of unmeasured confounding that could explain
#        away each estimate (VanderWeele & Ding closed form).
#  Outputs are AGGREGATE -> safe to share (nco_results.csv / evalue_table.csv).
# #############################################################################

hh_all <- readRDS("iptw_weighted_long.rds")

# =============================================================================
#  (1) Negative-control outcomes (binary, IPTW-weighted, Rubin-pooled)
# =============================================================================
analyze_binary <- function(dat, yvar) {
  m <- max(dat$.imp); R <- data.frame()
  for (i in 1:m) {
    di <- dat %>% filter(.imp == i, exposure_group %in% c("PPI","H2RA")) %>%
      mutate(treat = as.integer(exposure_group == "PPI"),
             y = as.integer(as.character(.data[[yvar]])))
    dz <- svydesign(ids = ~1, weights = ~ipw, data = di)
    f_rd <- svyglm(y ~ treat, design = dz, family = gaussian())
    f_rr <- svyglm(y ~ treat, design = dz, family = quasipoisson(link = "log"))
    w0 <- di$treat == 0; w1 <- di$treat == 1
    R <- rbind(R, data.frame(
      rd = coef(f_rd)[["treat"]], v_rd = vcov(f_rd)["treat","treat"],
      lrr = coef(f_rr)[["treat"]], v_lrr = vcov(f_rr)["treat","treat"],
      risk_h2 = weighted.mean(di$y[w0], di$ipw[w0]), risk_ppi = weighted.mean(di$y[w1], di$ipw[w1])))
  }
  rd <- pool_rubin(R$rd, R$v_rd); rr <- pool_rubin(R$lrr, R$v_lrr)
  list(n_ppi = sum(dat$.imp==1 & dat$exposure_group=="PPI"),
       n_h2  = sum(dat$.imp==1 & dat$exposure_group=="H2RA"),
       risk_ppi = mean(R$risk_ppi), risk_h2 = mean(R$risk_h2), RD = rd, RR = expci(rr))
}

res_vte <- analyze_binary(hh_all, "nco_vte")
res_pi  <- analyze_binary(hh_all, "nco_pressure")

show_nco <- function(r, lab) {
  cat("\n==========", lab, "(negative control, should be ~null) ==========\n")
  cat(sprintf("Weighted incidence: PPI %.1f%%   H2RA %.1f%%\n", 100*r$risk_ppi, 100*r$risk_h2))
  cat(sprintf("Risk diff RD : %+.1f%% (95%%CI %+.1f%% , %+.1f%%)  p=%.3f\n",
              100*unname(r$RD["est"]), 100*unname(r$RD["lo"]), 100*unname(r$RD["hi"]), unname(r$RD["p"])))
  cat(sprintf("Risk ratio RR : %.2f (95%%CI %.2f , %.2f)\n",
              unname(r$RR["est"]), unname(r$RR["lo"]), unname(r$RR["hi"])))
  cat("(RR~1 with CI including 1 -> no clear residual-confounding signal)\n")
}
show_nco(res_vte, "incident VTE")
show_nco(res_pi,  "stage>=2 pressure injury")

nco_row <- function(r, lab) data.frame(
  nco = lab, ppi_pct = round(100*r$risk_ppi,2), h2_pct = round(100*r$risk_h2,2),
  RD_pct = round(100*unname(r$RD["est"]),2), RD_lo = round(100*unname(r$RD["lo"]),2),
  RD_hi = round(100*unname(r$RD["hi"]),2), RD_p = signif(unname(r$RD["p"]),3),
  RR = round(unname(r$RR["est"]),3), RR_lo = round(unname(r$RR["lo"]),3),
  RR_hi = round(unname(r$RR["hi"]),3), row.names = NULL)
nco_tbl <- rbind(nco_row(res_vte,"VTE"), nco_row(res_pi,"pressure_injury"))
write.csv(nco_tbl, "nco_results.csv", row.names = FALSE)
cat("\n---- NCO summary (written to nco_results.csv) ----\n"); print(nco_tbl, row.names = FALSE)

# =============================================================================
#  (2) E-value (VanderWeele & Ding closed form; = EValue::evalues.RR)
#      Computed for the point estimate and for the CI bound closest to the null.
#      HR / sHR are treated approximately as RR. Reads the CSVs from R4-R6.
# =============================================================================
add_ev <- function(label, est, lo, hi) {
  e <- evalue_rr(est, lo, hi)
  data.frame(outcome = label, ratio = round(est,3), lo = round(lo,3), hi = round(hi,3),
             evalue_point = e[1], evalue_ci = e[2], row.names = NULL)
}

ev_tbl <- data.frame()
grab <- function(path, fun) tryCatch(fun(read.csv(path)), error = function(e) {
  cat("(skipped, file not found:", path, ")\n"); NULL })

ev_tbl <- rbind(ev_tbl,
                grab("primary_effects_iptw.csv", function(d) rbind(
                  add_ev("MAKE30 (RR)", d$RR[1], d$RR_lo[1], d$RR_hi[1]),
                  add_ev("MAKE90 (RR)", d$RR[2], d$RR_lo[2], d$RR_hi[2]))),
                grab("secondary_mortality.csv", function(d) rbind(
                  add_ev("Death 28d (RR)", d$RR[1], d$RR_lo[1], d$RR_hi[1]),
                  add_ev("Death 90d (RR)", d$RR[2], d$RR_lo[2], d$RR_hi[2]))),
                grab("secondary_rrt.csv", function(d)
                  add_ev("RRT (sHR)", d$sHR[1], d$sHR_lo[1], d$sHR_hi[1])),
                grab("secondary_recovery.csv", function(d)
                  add_ev("Recovery 90d (sHR)", d$sHR[2], d$sHR_lo[2], d$sHR_hi[2])),
                grab("secondary_progression.csv", function(d)
                  add_ev("Progression (sHR)", d$sHR[1], d$sHR_lo[1], d$sHR_hi[1])),
                grab("nco_results.csv", function(d) rbind(
                  add_ev("NCO: VTE (RR)", d$RR[1], d$RR_lo[1], d$RR_hi[1]),
                  add_ev("NCO: pressure (RR)", d$RR[2], d$RR_lo[2], d$RR_hi[2]))))

cat("\n---- E-value table (written to evalue_table.csv) ----\n")
print(ev_tbl, row.names = FALSE)
write.csv(ev_tbl, "evalue_table.csv", row.names = FALSE)
# Interpretation:
#  - If the NCO RRs are ~1 with CI including 1, there is no clear residual-
#    confounding signal; a systematic NCO shift toward PPI would suggest the
#    main estimates' small departures are residual confounding.
#  - E-values are naturally small for near-null results (point ~1.1-1.3, CI ~1);
#    this is expected. The safety conclusion rests on CI upper bounds excluding
#    harm together with clean negative controls.




# #############################################################################
#  R7 -- Forest plot (primary + secondary + negative controls)
#    Reads the aggregate effect CSVs from R4-R6 and renders a reproducible PDF
#    (main manuscript figure). Group-level aggregates only -> safe to share.
#    NOTE: the original builds the plot object for interactive export; here we
#    add a ggsave() so the figure is actually written. The plot spec is unchanged.
# #############################################################################

read_safe <- function(p) tryCatch(read.csv(p), error = function(e) {
  stop(paste("Missing", p, "-- run R4-R6 first in the same working directory")) })
prim <- read_safe("primary_effects_iptw.csv")
mort <- read_safe("secondary_mortality.csv")
rrt  <- read_safe("secondary_rrt.csv")
rec  <- read_safe("secondary_recovery.csv")
prog <- read_safe("secondary_progression.csv")
nco  <- read_safe("nco_results.csv")

## 1. Assemble one long table -----------------------------------------------
fp <- bind_rows(
  data.frame(label="MAKE30",           est=prim$RR[1],  lo=prim$RR_lo[1],  hi=prim$RR_hi[1],  measure="RR",  group="Primary outcomes"),
  data.frame(label="MAKE90",           est=prim$RR[2],  lo=prim$RR_lo[2],  hi=prim$RR_hi[2],  measure="RR",  group="Primary outcomes"),
  data.frame(label="28-day mortality", est=mort$RR[1],  lo=mort$RR_lo[1],  hi=mort$RR_hi[1],  measure="RR",  group="Secondary outcomes"),
  data.frame(label="90-day mortality", est=mort$RR[2],  lo=mort$RR_lo[2],  hi=mort$RR_hi[2],  measure="RR",  group="Secondary outcomes"),
  data.frame(label="RRT initiation",   est=rrt$sHR[1],  lo=rrt$sHR_lo[1],  hi=rrt$sHR_hi[1],  measure="sHR", group="Secondary outcomes"),
  data.frame(label="Renal recovery",   est=rec$sHR[2],  lo=rec$sHR_lo[2],  hi=rec$sHR_hi[2],  measure="sHR", group="Secondary outcomes"),
  data.frame(label="AKI progression",  est=prog$sHR[1], lo=prog$sHR_lo[1], hi=prog$sHR_hi[1], measure="sHR", group="Secondary outcomes"),
  data.frame(label="VTE",              est=nco$RR[1],   lo=nco$RR_lo[1],   hi=nco$RR_hi[1],   measure="RR",  group="Negative controls"),
  data.frame(label="Pressure injury",  est=nco$RR[2],   lo=nco$RR_lo[2],   hi=nco$RR_hi[2],   measure="RR",  group="Negative controls"))

fp$group <- factor(fp$group, levels = c("Primary outcomes","Secondary outcomes","Negative controls"))
fp$label <- factor(fp$label, levels = rev(fp$label))           # reverse -> first row on top
fp$txt   <- sprintf("%s %.2f (%.2f-%.2f)", fp$measure, fp$est, fp$lo, fp$hi)
fp$col   <- ifelse(fp$group == "Negative controls", "nco", "real")

## 2. Plot ------------------------------------------------------------------
xmax_txt <- 2.8   # right-hand range reserved for the numeric labels (log axis)
p <- ggplot(fp, aes(x = est, y = label, color = col)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey55", linewidth = 0.4) +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.20, linewidth = 0.6) +
  geom_point(shape = 15, size = 2.8) +
  geom_text(aes(x = 1.42, label = txt), hjust = 0, size = 2.9, color = "grey20") +
  facet_grid(rows = vars(group), scales = "free_y", space = "free_y", switch = "y") +
  scale_x_log10(breaks = c(0.8, 0.9, 1.0, 1.1, 1.25), limits = c(0.74, xmax_txt)) +
  scale_color_manual(values = c(real = "#185FA5", nco = "#5F5E5A"), guide = "none") +
  coord_cartesian(clip = "off") +
  labs(x = "Effect size RR / sHR (95% CI; dashed line 1.0 = no difference)", y = NULL,
       subtitle = "PPI vs H2RA - IPTW-weighted - 20 imputations pooled (Rubin's rules)",
       caption  = "RR = risk ratio; sHR = Fine-Gray subdistribution hazard ratio. Except renal recovery, estimates <1 favor PPI.") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor   = element_blank(),
        panel.grid.major.x = element_line(linewidth = 0.25, color = "grey90"),
        strip.placement = "outside",
        strip.text.y.left = element_text(angle = 0, face = "bold", hjust = 1, size = 10),
        axis.text.y = element_text(size = 10, color = "black"),
        plot.subtitle = element_text(size = 9, color = "grey40"),
        plot.caption  = element_text(size = 8, color = "grey45", hjust = 0),
        plot.margin = margin(10, 16, 10, 10))

## 3. Save (cairo_pdf if available, else default device) --------------------
ok <- tryCatch({ ggsave("figure_forest_main.pdf", p, width = 8, height = 5, device = cairo_pdf); TRUE },
               error = function(e) FALSE)
if (!ok) ggsave("figure_forest_main.pdf", p, width = 8, height = 5)
cat("\nForest plot written to figure_forest_main.pdf\n")






# #############################################################################
#  R8 -- Publication tables and figures (5 parts, runnable standalone)
#    (1) Table 2 outcomes table (reads the R4-R6 CSVs)
#    (2) Figure 1 study flow diagram (STROBE)
#    (3) Propensity-score overlap (positivity)
#    (4) Weighted cumulative-incidence curves (RRT / recovery / progression)
#    (5) Weighted survival curve (mortality)
#  All group-level aggregates / curves -> safe to share.
#  NOTE: figures are saved with ggsave() here so the script produces the files;
#        the plot specifications are unchanged from the locked engine.
# #############################################################################

# =============================================================================
#  (1) Table 2 -- primary + secondary + negative-control outcomes
# =============================================================================
rd <- function(p) tryCatch(read.csv(p), error = function(e) NULL)
prim<-rd("primary_effects_iptw.csv"); mort<-rd("secondary_mortality.csv")
rrt<-rd("secondary_rrt.csv"); rec<-rd("secondary_recovery.csv")
prog<-rd("secondary_progression.csv"); nco<-rd("nco_results.csv"); ev<-rd("evalue_table.csv")

ci <- function(e,l,h,d=2) sprintf(paste0("%.",d,"f (%.",d,"f to %.",d,"f)"), e,l,h)
evget <- function(name){ if(is.null(ev)) return(NA); r<-ev[ev$outcome==name,]; if(nrow(r)) r$evalue_point[1] else NA }
row2 <- function(grp,out,ppi,h2,absd,eff,p,evn) data.frame(
  Group=grp, Outcome=out, PPI=ppi, H2RA=h2, `Abs.diff`=absd,
  `Effect (95% CI)`=eff, P=p, `E-value`=evget(evn), check.names=FALSE)

t2 <- bind_rows(
  row2("Primary","MAKE30", sprintf("%.1f%%",prim$risk_ppi_pct[1]), sprintf("%.1f%%",prim$risk_h2_pct[1]),
       paste0(ci(prim$RD_pct[1],prim$RD_lo[1],prim$RD_hi[1]),"%"),
       paste0("RR ",ci(prim$RR[1],prim$RR_lo[1],prim$RR_hi[1])), prim$RD_p[1], "MAKE30 (RR)"),
  row2("Primary","MAKE90", sprintf("%.1f%%",prim$risk_ppi_pct[2]), sprintf("%.1f%%",prim$risk_h2_pct[2]),
       paste0(ci(prim$RD_pct[2],prim$RD_lo[2],prim$RD_hi[2]),"%"),
       paste0("RR ",ci(prim$RR[2],prim$RR_lo[2],prim$RR_hi[2])), prim$RD_p[2], "MAKE90 (RR)"),
  row2("Secondary","28-day mortality", sprintf("%.1f%%",mort$ppi_pct[1]), sprintf("%.1f%%",mort$h2_pct[1]),
       paste0(ci(mort$RD_pct[1],mort$RD_lo[1],mort$RD_hi[1]),"%"),
       paste0("RR ",ci(mort$RR[1],mort$RR_lo[1],mort$RR_hi[1])), mort$RD_p[1], "Death 28d (RR)"),
  row2("Secondary","90-day mortality", sprintf("%.1f%%",mort$ppi_pct[2]), sprintf("%.1f%%",mort$h2_pct[2]),
       paste0(ci(mort$RD_pct[2],mort$RD_lo[2],mort$RD_hi[2]),"%"),
       paste0("RR ",ci(mort$RR[2],mort$RR_lo[2],mort$RR_hi[2])), mort$RD_p[2], "Death 90d (RR)"),
  row2("Secondary","RRT initiation", sprintf("%.1f%%",rrt$ppi_cif_pct[1]), sprintf("%.1f%%",rrt$h2_cif_pct[1]),
       sprintf("%+.1f%%",rrt$ppi_cif_pct[1]-rrt$h2_cif_pct[1]),
       paste0("sHR ",ci(rrt$sHR[1],rrt$sHR_lo[1],rrt$sHR_hi[1])), rrt$sHR_p[1], "RRT (sHR)"),
  row2("Secondary","Renal recovery", sprintf("%.1f%%",rec$ppi_cif_pct[2]), sprintf("%.1f%%",rec$h2_cif_pct[2]),
       sprintf("%+.1f%%",rec$ppi_cif_pct[2]-rec$h2_cif_pct[2]),
       paste0("sHR ",ci(rec$sHR[2],rec$sHR_lo[2],rec$sHR_hi[2])), rec$sHR_p[2], "Recovery 90d (sHR)"),
  row2("Secondary","AKI progression", sprintf("%.1f%%",prog$ppi_cif_pct[1]), sprintf("%.1f%%",prog$h2_cif_pct[1]),
       sprintf("%+.1f%%",prog$ppi_cif_pct[1]-prog$h2_cif_pct[1]),
       paste0("sHR ",ci(prog$sHR[1],prog$sHR_lo[1],prog$sHR_hi[1])), prog$sHR_p[1], "Progression (sHR)"),
  row2("Negative control","VTE", sprintf("%.1f%%",nco$ppi_pct[1]), sprintf("%.1f%%",nco$h2_pct[1]),
       paste0(ci(nco$RD_pct[1],nco$RD_lo[1],nco$RD_hi[1]),"%"),
       paste0("RR ",ci(nco$RR[1],nco$RR_lo[1],nco$RR_hi[1])), nco$RD_p[1], "NCO: VTE (RR)"),
  row2("Negative control","Pressure injury", sprintf("%.1f%%",nco$ppi_pct[2]), sprintf("%.1f%%",nco$h2_pct[2]),
       paste0(ci(nco$RD_pct[2],nco$RD_lo[2],nco$RD_hi[2]),"%"),
       paste0("RR ",ci(nco$RR[2],nco$RR_lo[2],nco$RR_hi[2])), nco$RD_p[2], "NCO: pressure (RR)"))

# RMST sub-table (MAKE / mortality only; difference in days)
t2_rmst <- bind_rows(
  data.frame(Outcome="MAKE30",  `RMST diff (days, 95% CI)`=ci(prim$RMSTd_day[1],prim$RMSTd_lo[1],prim$RMSTd_hi[1]), P=prim$RMSTd_p[1], check.names=FALSE),
  data.frame(Outcome="MAKE90",  `RMST diff (days, 95% CI)`=ci(prim$RMSTd_day[2],prim$RMSTd_lo[2],prim$RMSTd_hi[2]), P=prim$RMSTd_p[2], check.names=FALSE),
  data.frame(Outcome="28-day mortality", `RMST diff (days, 95% CI)`=ci(mort$RMSTd[1],mort$RMSTd_lo[1],mort$RMSTd_hi[1]), P=mort$RMSTd_p[1], check.names=FALSE),
  data.frame(Outcome="90-day mortality", `RMST diff (days, 95% CI)`=ci(mort$RMSTd[2],mort$RMSTd_lo[2],mort$RMSTd_hi[2]), P=mort$RMSTd_p[2], check.names=FALSE))

write.csv(t2, "table2_outcomes.csv", row.names = FALSE)
write.csv(t2_rmst, "table2_rmst.csv", row.names = FALSE)
cat("===== Table 2 =====\n"); print(t2, row.names = FALSE)
cat("\n===== Table 2 (RMST) =====\n"); print(t2_rmst, row.names = FALSE)

# =============================================================================
#  (2) Figure 1 -- study flow diagram (STROBE)
# =============================================================================
box <- function(x,y,lab,type) data.frame(x=x,y=y,lab=lab,type=type)
boxes <- bind_rows(
  box(3.5,11.0,"Inclusion criteria\nSepsis-3 | age >=18 | first ICU stay\nFirst AKI within 7 d of sepsis (KDIGO >=1)\nN = 23,745","main"),
  box(8.0, 9.8,"Excluded:\nESRD / chronic dialysis  1,542\nKidney transplant  475\nPregnancy / puerperium  55","excl"),
  box(3.5, 8.6,"Eligible cohort\nN = 21,883","main"),
  box(8.0, 7.7,"Excluded: timestamp error /\nmissing MAKE  139","excl"),
  box(3.5, 6.5,"Analytic cohort\nN = 21,744","main"),
  box(1.3, 4.6,"PPI\n3,439","grp"),
  box(3.5, 4.6,"H2RA\n3,075","grp"),
  box(5.7, 4.6,"No acid suppression\n15,230","grp"),
  box(2.4, 2.4,"Head-to-head IPTW analysis\nPPI 3,439  vs  H2RA 3,075\nN = 6,514","hl"))

seg <- function(x,xe,y,ye) data.frame(x=x,xend=xe,y=y,yend=ye)
arrows <- bind_rows(
  seg(3.5,3.5,10.5,9.1), seg(3.5,6.55,9.8,9.8),     # A->B, ->excl1
  seg(3.5,3.5, 8.1,6.9), seg(3.5,6.55,7.7,7.7),     # B->C, ->excl2
  seg(3.5,3.5, 6.1,5.25),                            # C->allocation bar
  seg(1.3,5.7, 5.25,5.25),                          # horizontal bar
  seg(1.3,1.3, 5.25,5.0), seg(3.5,3.5,5.25,5.0), seg(5.7,5.7,5.25,5.0),  # bar -> three groups
  seg(1.3,1.3, 4.1,3.0),  seg(3.5,3.5,4.1,3.0),     # PPI / H2RA downward
  seg(1.3,3.5, 3.0,3.0),  seg(2.4,2.4,3.0,2.75))    # merge -> head-to-head
no_arrow <- 6:7  # horizontal-bar segments drawn without arrowheads

f1 <- ggplot() +
  geom_segment(data=arrows[-no_arrow,], aes(x=x,xend=xend,y=y,yend=yend),
               arrow=arrow(length=unit(0.16,"cm"),type="closed"), linewidth=0.4, color="grey40") +
  geom_segment(data=arrows[no_arrow,], aes(x=x,xend=xend,y=y,yend=yend), linewidth=0.4, color="grey40") +
  geom_label(data=boxes, aes(x=x,y=y,label=lab,fill=type),
             color="grey15", size=2.9, lineheight=0.95, label.size=0.3,
             label.r=unit(0.12,"lines"), label.padding=unit(0.4,"lines")) +
  scale_fill_manual(values=c(main="#E6F1FB", excl="#F1EFE8", grp="#E1F5EE", hl="#FAEEDA"), guide="none") +
  scale_x_continuous(limits=c(-0.3,10)) + scale_y_continuous(limits=c(1.6,11.8)) +
  theme_void()

# =============================================================================
#  (3) Propensity-score overlap (PS refit on imp #1)
# =============================================================================
hh <- readRDS("iptw_weighted_long.rds")
d1 <- hh %>% filter(.imp==1, exposure_group %in% c("PPI","H2RA")) %>%
  mutate(treat = as.integer(exposure_group=="PPI"))
psfit <- glm(reformulate(ps_covars, "treat"), data=d1, family=binomial)
d1$ps <- fitted(psfit)

f_ps <- ggplot(d1, aes(x=ps, fill=factor(treat, labels=c("H2RA","PPI")))) +
  geom_density(alpha=0.45, color=NA) +
  scale_fill_manual(values=c(H2RA="#888780", PPI="#185FA5"), name=NULL) +
  labs(x="Propensity score (P[receive PPI])", y="Density",
       subtitle="Overlapping PS distributions -> positivity assumption met") +
  theme_minimal(base_size=11) +
  theme(legend.position=c(0.85,0.85), panel.grid.minor=element_blank())

# =============================================================================
#  (4) Weighted cumulative-incidence curves (RRT / recovery / progression; imp #1)
#      Competing event = death; curves are illustrative, formal inference in Table 2
# =============================================================================
build_cr <- function(dat, elig, evvar, tvar, tau, drop_pre_rrt=FALSE) {
  d <- dat %>% filter(.imp==1, exposure_group %in% c("PPI","H2RA"))
  if (!is.null(elig)) d <- d %>% filter(as.integer(as.character(.data[[elig]]))==1)
  if (drop_pre_rrt)   d <- d %>% filter(as.integer(as.character(rrt_pre_t0))==0)
  d %>% mutate(
    treat = factor(ifelse(exposure_group=="PPI","PPI","H2RA"), levels=c("H2RA","PPI")),
    ev1   = as.integer(as.character(.data[[evvar]])),
    t_ev  = ifelse(ev1==1, .data[[tvar]], Inf),
    t_dth = ifelse(as.integer(as.character(death_inhosp))==1, days_to_death, Inf),
    ctime = pmin(los_after_t0_days, tau),
    fe    = pmin(t_ev, t_dth),
    status= ifelse(fe>ctime, 0L, ifelse(t_ev<=t_dth, 1L, 2L)),
    ftime = pmin(fe, ctime),
    ef    = factor(status, levels=c(0,1,2), labels=c("cens","event","death")))
}
cif_curve <- function(d, lab) {
  do.call(rbind, lapply(c("H2RA","PPI"), function(g){
    s <- d[d$treat==g,]
    sf <- survfit(Surv(ftime, ef) ~ 1, data=s, weights=s$ipw)
    j <- which(colnames(sf$pstate)=="event")
    rbind(data.frame(time=0,cif=0,treat=g,outcome=lab),
          data.frame(time=sf$time, cif=sf$pstate[,j], treat=g, outcome=lab))
  }))
}
cif_all <- bind_rows(
  cif_curve(build_cr(hh, NULL,                "rrt_init",      "days_to_rrt",      90, TRUE),  "RRT initiation"),
  cif_curve(build_cr(hh, "recovery_eligible", "recovery_event","days_to_recovery", 90),        "Renal recovery"),
  cif_curve(build_cr(hh, "prog_eligible",     "prog_event",    "days_to_prog",     90),        "AKI progression"))
cif_all$outcome <- factor(cif_all$outcome, levels=c("RRT initiation","Renal recovery","AKI progression"))

f_cif <- ggplot(cif_all, aes(time, cif, color=treat)) +
  geom_step(linewidth=0.7) +
  facet_wrap(~outcome, scales="free_y") +
  scale_color_manual(values=c(H2RA="#888780", PPI="#185FA5"), name=NULL) +
  scale_y_continuous(labels=scales::percent) +
  labs(x="Days since T0", y="Weighted cumulative incidence",
       subtitle="Competing event = death; curves illustrative (inference: Table 2)") +
  theme_minimal(base_size=11) +
  theme(legend.position="top", panel.grid.minor=element_blank())

# =============================================================================
#  (5) Weighted survival curve (mortality, 90 days; imp #1)
# =============================================================================
dm <- hh %>% filter(.imp==1, exposure_group %in% c("PPI","H2RA")) %>%
  mutate(treat = factor(ifelse(exposure_group=="PPI","PPI","H2RA"), levels=c("H2RA","PPI")),
         d90 = as.integer(as.character(death_90)), s90 = as.numeric(surv_90))
km <- bind_rows(lapply(c("H2RA","PPI"), function(g){
  s <- dm[dm$treat==g,]
  sf <- survfit(Surv(s90, d90) ~ 1, data=s, weights=s$ipw)
  rbind(data.frame(time=0, cm=0, treat=g),
        data.frame(time=sf$time, cm=1-sf$surv, treat=g))   # cumulative mortality = 1 - survival
}))
f_km <- ggplot(km, aes(time, cm, color=treat)) +
  geom_step(linewidth=0.7) +
  scale_color_manual(values=c(H2RA="#888780", PPI="#185FA5"), name=NULL) +
  scale_y_continuous(labels=scales::percent) +
  labs(x="Days since T0", y="Weighted cumulative mortality", subtitle="90-day all-cause mortality (IPTW-weighted)") +
  theme_minimal(base_size=11) +
  theme(legend.position=c(0.15,0.85), panel.grid.minor=element_blank())

## Save figures (cairo_pdf if available, else default device) ---------------
save_fig <- function(obj, file, w, h) {
  ok <- tryCatch({ ggsave(file, obj, width=w, height=h, device=cairo_pdf); TRUE },
                 error = function(e) FALSE)
  if (!ok) ggsave(file, obj, width=w, height=h)
}
save_fig(f1,    "figure1_flowchart.pdf",    8, 9)
save_fig(f_ps,  "figure_ps_overlap.pdf",    6, 4.5)
save_fig(f_cif, "figure_cif_curves.pdf",    9, 3.5)
save_fig(f_km,  "figure_km_mortality.pdf",  6, 4.5)
cat("\nFigures written: figure1_flowchart.pdf, figure_ps_overlap.pdf,",
    "figure_cif_curves.pdf, figure_km_mortality.pdf\n")
# Outputs: table2_outcomes.csv / table2_rmst.csv / the four figure PDFs above.






# #############################################################################
#  R9 -- Subgroup analysis (5 pre-specified subgroups) + forest plot
#    Subgroups (biological / clinical rationale):
#      (1) AKI stage (KDIGO 1 vs >=2)  (2) pre-existing CKD
#      (3) septic shock (vasopressor)  (4) age (<65 / >=65)  (5) sex
#    Method: within each subgroup cell, re-estimate PS weights (cell-specific
#      balance), compute weighted MAKE30/90 RR per imputation -> Rubin pooled.
#      Effect modification = between-subgroup heterogeneity test (binary Wald /
#      multi-level Cochran Q, same source as the within-cell estimates).
#      Framed as "showing the null is consistent across subgroups".
#    AKI is collapsed to 1 vs >=2: stage 3 alone is too small (n~197, H2RA ~67)
#      and causes PS separation (SMD ~0.27); collapsed >=2 (n~662) balances
#      better (SMD ~0.14). For 1/2/3 levels, change aki_grp.
#    Uses the helpers and ps_covars already defined in PART 0.
#  Outputs: mimic_subgroup_effects.csv + mimic_subgroup_forest_make30.pdf
#  Note: many cell re-weightings -> may take a few minutes; small cells may emit
#        fitted-probability 0/1 warnings (expected). Group-level aggregates only.
#  This module replaces the old (age/sex-only) subgroup block of the locked engine.
# #############################################################################

# =============================================================================
#  (1) Load + derive the 5 subgroup variables (AKI collapsed to 1 vs >=2)
# =============================================================================
hh <- readRDS("iptw_weighted_long.rds")
hh <- hh %>% mutate(
  aki_grp   = factor(ifelse(as.character(first_aki_stage) == "1", "Stage1", "Stage2plus"),
                     levels = c("Stage1","Stage2plus")),
  ckd_grp   = factor(ifelse(as.integer(as.character(ckd)) == 1, "CKD", "NoCKD"), levels = c("NoCKD","CKD")),
  shock_grp = factor(ifelse(as.integer(as.character(vasopressor)) == 1, "Yes", "No"),  levels = c("No","Yes")),
  age_group = factor(ifelse(age_at_icu < 65, "<65", ">=65"), levels = c("<65",">=65")),
  sex       = droplevels(factor(as.character(gender))))
cat("AKI collapsed (imp#1, head-to-head):\n")
print(table(hh$aki_grp[hh$.imp == 1 & hh$exposure_group %in% c("PPI","H2RA")]))

## Subgroup definitions (name / variable / levels / display label) ----
subgroups <- list(
  list(name = "AKI stage (KDIGO)", var = "aki_grp",   levs = c("Stage1","Stage2plus"),
       disp = c("KDIGO stage 1","KDIGO stage >=2")),
  list(name = "Pre-existing CKD",  var = "ckd_grp",   levs = c("NoCKD","CKD"), disp = c("No CKD","CKD")),
  list(name = "Septic shock", var = "shock_grp", levs = c("No","Yes"),
       disp = c("No (no vasopressor)","Yes (vasopressor)")),
  list(name = "Age",      var = "age_group", levs = c("<65",">=65"),  disp = c("<65 y",">=65 y")),
  list(name = "Sex",      var = "sex",       levs = levels(hh$sex),
       disp = ifelse(levels(hh$sex) == "F", "Female", "Male")))

# =============================================================================
#  (2) One cell: re-estimate weights within cell -> MAKE30 & MAKE90 weighted RR
#      (same weights for both endpoints, saving half the PS refits)
# =============================================================================
cell_rr <- function(dat, svar = NULL, level = NULL) {
  m <- max(dat$.imp); R <- data.frame(); smd <- c(); np <- NA; nh <- NA
  for (i in 1:m) {
    di <- dat %>% filter(.imp == i, exposure_group %in% c("PPI","H2RA"))
    if (!is.null(svar)) di <- di %>% filter(.data[[svar]] == level)
    di <- di %>% mutate(treat = as.integer(exposure_group == "PPI"))
    if (i == 1) { np <- sum(di$treat == 1); nh <- sum(di$treat == 0) }
    if (min(sum(di$treat == 1), sum(di$treat == 0)) < 10) return(NULL)   # cell too small -> drop
    kc <- ps_covars[sapply(ps_covars, function(vv) length(unique(na.omit(di[[vv]]))) > 1)]
    res <- tryCatch({
      W  <- weightit(reformulate(kc, "treat"), data = di, method = "glm",
                     estimand = "ATE", stabilize = TRUE)
      Wt <- trim(W, at = 0.99, lower = TRUE); di$w <- Wt$weights
      sm <- max(abs(bal.tab(Wt, un = FALSE, binary = "std", continuous = "std")$Balance$Diff.Adj), na.rm = TRUE)
      get1 <- function(yv) {
        di$y <- as.integer(as.character(di[[yv]])); d2 <- di[!is.na(di$y), ]
        dz <- svydesign(ids = ~1, weights = ~w, data = d2)
        fr <- svyglm(y ~ treat, design = dz, family = quasipoisson(link = "log"))
        c(lrr = unname(coef(fr)[["treat"]]), v = unname(vcov(fr)["treat","treat"]),
          rp = weighted.mean(d2$y[d2$treat == 1], d2$w[d2$treat == 1]),
          rh = weighted.mean(d2$y[d2$treat == 0], d2$w[d2$treat == 0]))
      }
      list(m30 = get1("make30"), m90 = get1("make90"), sm = sm)
    }, error = function(e) NULL)
    if (is.null(res)) next
    R <- rbind(R, data.frame(lrr30 = res$m30["lrr"], v30 = res$m30["v"],
                             rp30 = res$m30["rp"], rh30 = res$m30["rh"], lrr90 = res$m90["lrr"], v90 = res$m90["v"]))
    smd <- c(smd, res$sm)
  }
  if (nrow(R) < 2) return(NULL)
  rr30 <- pool_rubin(R$lrr30, R$v30); rr90 <- pool_rubin(R$lrr90, R$v90)
  list(n_ppi = np, n_h2 = nh, max_smd = round(mean(smd), 3),
       RR30 = expci(rr30), logrr30 = unname(rr30["est"]), v_lrr30 = unname(rr30["se"])^2,
       RR90 = expci(rr90), logrr90 = unname(rr90["est"]), v_lrr90 = unname(rr90["se"])^2,
       risk_ppi30 = mean(R$rp30), risk_h2_30 = mean(R$rh30))
}

# =============================================================================
#  (3) Overall + each subgroup level + between-subgroup heterogeneity P
# =============================================================================
het_p <- function(logrr, v) {                 # binary Wald / multi-level Cochran Q
  K <- length(logrr)
  if (K < 2 || any(is.na(v)) || any(v <= 0)) return(NA_real_)
  if (K == 2) { z <- (logrr[1] - logrr[2]) / sqrt(v[1] + v[2]); 2*(1 - pnorm(abs(z))) }
  else { wf <- 1/v; lf <- sum(wf*logrr)/sum(wf); Q <- sum(wf*(logrr - lf)^2); 1 - pchisq(Q, K - 1) }
}

cat("\n... computing overall effect (full head-to-head, PS refit)\n")
ov <- cell_rr(hh)

sg_rows <- list(); het_rows <- data.frame()
for (sg in subgroups) {
  cat(sprintf("... subgroup: %s (%d levels)\n", sg$name, length(sg$levs)))
  lr30 <- c(); vv30 <- c(); lr90 <- c(); vv90 <- c()
  for (k in seq_along(sg$levs)) {
    e <- cell_rr(hh, sg$var, sg$levs[k])
    if (is.null(e)) { cat(sprintf("    [skip] %s too small / not converged\n", sg$disp[k])); next }
    sg_rows[[length(sg_rows)+1]] <- data.frame(
      subgroup = sg$name, level = sg$disp[k], n_ppi = e$n_ppi, n_h2 = e$n_h2, max_smd = e$max_smd,
      RR30 = round(e$RR30["est"],3), lo30 = round(e$RR30["lo"],3), hi30 = round(e$RR30["hi"],3),
      RR90 = round(e$RR90["est"],3), lo90 = round(e$RR90["lo"],3), hi90 = round(e$RR90["hi"],3))
    lr30 <- c(lr30, e$logrr30); vv30 <- c(vv30, e$v_lrr30)
    lr90 <- c(lr90, e$logrr90); vv90 <- c(vv90, e$v_lrr90)
  }
  het_rows <- rbind(het_rows, data.frame(subgroup = sg$name,
                                         p_int_make30 = round(het_p(lr30, vv30), 3), p_int_make90 = round(het_p(lr90, vv90), 3)))
}
sg_tbl <- bind_rows(sg_rows) %>% left_join(het_rows, by = "subgroup")
write.csv(sg_tbl, "mimic_subgroup_effects.csv", row.names = FALSE)
cat("\n========== Subgroup effects (written to mimic_subgroup_effects.csv) ==========\n")
cat(sprintf("Overall MAKE30 RR: %.3f (%.3f-%.3f)  [PPI %d / H2RA %d]\n",
            ov$RR30["est"], ov$RR30["lo"], ov$RR30["hi"], ov$n_ppi, ov$n_h2))
print(sg_tbl, row.names = FALSE)
cat("\nCheck: are all P-interaction(make30) > 0.05/0.10 -> null consistent across subgroups, no effect modification.\n")

# =============================================================================
#  (4) Forest plot (MAKE30, primary; numbers placed inside the panel, no clip)
# =============================================================================
fr <- list(); kk <- 0
addrow <- function(label, est = NA, lo = NA, hi = NA, txt = "") { kk <<- kk + 1
fr[[kk]] <<- data.frame(ord = kk, label = label, est = est, lo = lo, hi = hi, txt = txt) }

addrow(sprintf("Overall (PPI %d / H2RA %d)", ov$n_ppi, ov$n_h2),
       ov$RR30["est"], ov$RR30["lo"], ov$RR30["hi"],
       sprintf("%.2f (%.2f-%.2f)", ov$RR30["est"], ov$RR30["lo"], ov$RR30["hi"]))
addrow("")
for (sg in subgroups) {
  pint <- het_rows$p_int_make30[het_rows$subgroup == sg$name]
  addrow(sprintf("%s   (P-interaction = %s)", sg$name, ifelse(length(pint)==0 || is.na(pint), "NA", sprintf("%.2f", pint))))
  sub <- sg_tbl[sg_tbl$subgroup == sg$name, ]
  for (j in seq_len(nrow(sub)))
    addrow(sprintf("    %s  (n=%d)", sub$level[j], sub$n_ppi[j] + sub$n_h2[j]),
           sub$RR30[j], sub$lo30[j], sub$hi30[j],
           sprintf("%.2f (%.2f-%.2f)", sub$RR30[j], sub$lo30[j], sub$hi30[j]))
  addrow("")
}
F <- bind_rows(fr); F$y <- max(F$ord) - F$ord + 1

xlo  <- min(0.55, min(F$lo, na.rm = TRUE) * 0.92)
xdh  <- max(1.55, max(F$hi, na.rm = TRUE) * 1.08)   # data / CI right edge
xnum <- xdh * 1.08                                  # numeric label start (inside panel)
xph  <- xnum * 1.70                                 # panel right edge (room for numbers)

p <- ggplot(F, aes(y = y)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey55") +
  geom_errorbarh(aes(xmin = lo, xmax = hi), na.rm = TRUE, height = .3, color = "#185FA5", linewidth = .6) +
  geom_point(aes(x = est), na.rm = TRUE, size = 2.5, color = "#185FA5") +
  geom_text(aes(x = xnum, label = txt), hjust = 0, size = 2.9, color = "grey20", na.rm = TRUE) +
  scale_x_log10(breaks = c(0.5, 0.7, 1, 1.5, 2), limits = c(xlo, xph)) +
  scale_y_continuous(breaks = F$y, labels = F$label, expand = expansion(add = 0.8)) +
  labs(x = "MAKE30 risk ratio RR (PPI vs H2RA, log scale; <1 favors PPI)", y = NULL,
       title = "PPI vs H2RA in SA-AKI - MAKE30 subgroup analysis (MIMIC)",
       subtitle = "IPTW re-weighted within subgroups, 20 imputations pooled (Rubin); P-interaction = between-subgroup heterogeneity test") +
  theme_minimal(base_size = 10) +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
        axis.text.y = element_text(hjust = 0), plot.title = element_text(face = "bold"))

## Prefer cairo_pdf; fall back to default device if Cairo unavailable --------
ok <- tryCatch({ ggsave("mimic_subgroup_forest_make30.pdf", p, width = 10.5, height = 7, device = cairo_pdf); TRUE },
               error = function(e) { ggsave("mimic_subgroup_forest_make30.pdf", p, width = 10.5, height = 7); FALSE })
cat(sprintf("\nSubgroup forest plot written to mimic_subgroup_forest_make30.pdf (%s)\n",
            if (ok) "cairo_pdf" else "default device"))
# For a MAKE90 forest, swap RR30/lo30/hi30 for RR90/lo90/hi90.





# #############################################################################
#  R10 -- New-user sensitivity (re-estimate weights within the new-user subset)
#    Restricts to first-ever users (newuser == 1), re-estimates PS / IPTW within
#    that subset (subset-specific balance), and re-runs the primary MAKE30/90
#    analysis (RD / RR / RMST). Compared with the full-cohort primary, a shift
#    toward the null indicates prevalent-user bias.
#    IPTW-weighted; 20 imputations pooled by Rubin's rules.
#  Output: newuser_results.csv (aggregate -> safe to share).
# #############################################################################

hh <- readRDS("iptw_weighted_long.rds")

newuser_all <- function(dat, tau30=30, tau90=90) {
  m <- max(dat$.imp); R30<-data.frame(); R90<-data.frame()
  maxsmd<-numeric(m); nppi<-numeric(m); nh2<-numeric(m)
  one <- function(di, yvar, tau) {
    di$y  <- as.integer(as.character(di[[yvar]]))
    di$tt <- pmin(ifelse(as.integer(as.character(di$make_tte_event))==1, di$make_tte_time, tau), tau)
    d2 <- di[!is.na(di$y),]
    dz2 <- svydesign(ids=~1, weights=~w, data=d2)
    dzf <- svydesign(ids=~1, weights=~w, data=di)
    frd <- svyglm(y~treat, design=dz2, family=gaussian())
    frr <- svyglm(y~treat, design=dz2, family=quasipoisson(link="log"))
    frm <- svyglm(tt~treat, design=dzf, family=gaussian())
    w0<-d2$treat==0; w1<-d2$treat==1
    data.frame(rd=coef(frd)[["treat"]],v_rd=vcov(frd)["treat","treat"],
               lrr=coef(frr)[["treat"]],v_lrr=vcov(frr)["treat","treat"],
               rm=coef(frm)[["treat"]],v_rm=vcov(frm)["treat","treat"],
               rppi=weighted.mean(d2$y[w1],d2$w[w1]), rh2=weighted.mean(d2$y[w0],d2$w[w0]))
  }
  for (i in 1:m) {
    di <- dat %>% filter(.imp==i, exposure_group %in% c("PPI","H2RA"),
                         as.integer(as.character(newuser))==1) %>%
      mutate(treat=as.integer(exposure_group=="PPI"))
    kc <- ps_covars[sapply(ps_covars, function(vv) length(unique(di[[vv]]))>1)]   # drop constant-in-subset columns
    W  <- weightit(reformulate(kc,"treat"), data=di, method="glm", estimand="ATE", stabilize=TRUE)
    Wt <- trim(W, at=0.99, lower=TRUE); di$w <- Wt$weights
    bt <- bal.tab(Wt, un=FALSE, binary="std", continuous="std")
    maxsmd[i] <- max(abs(bt$Balance$Diff.Adj), na.rm=TRUE)
    nppi[i] <- sum(di$treat==1); nh2[i] <- sum(di$treat==0)
    R30 <- rbind(R30, one(di,"make30",tau30)); R90 <- rbind(R90, one(di,"make90",tau90))
  }
  fin <- function(R) list(rppi=mean(R$rppi), rh2=mean(R$rh2),
                          RD=pool_rubin(R$rd,R$v_rd), RR=expci(pool_rubin(R$lrr,R$v_lrr)),
                          RMST=pool_rubin(R$rm,R$v_rm))
  list(n_ppi=round(mean(nppi)), n_h2=round(mean(nh2)), max_smd=round(mean(maxsmd),3),
       make30=fin(R30), make90=fin(R90))
}

nu <- newuser_all(hh)

cat("\n===== New-user sensitivity =====\n")
cat(sprintf("N (new-user head-to-head): PPI %d vs H2RA %d; after re-weighting max|SMD| = %.3f (<0.1 desirable)\n",
            nu$n_ppi, nu$n_h2, nu$max_smd))
show_nu <- function(r, lab) {
  cat(sprintf("\n[%s] weighted risk PPI %.1f%% / H2RA %.1f%%\n", lab, 100*r$rppi, 100*r$rh2))
  cat(sprintf("  RD  : %+.1f%% (95%%CI %+.1f%% , %+.1f%%)  p=%.3f\n",
              100*unname(r$RD["est"]),100*unname(r$RD["lo"]),100*unname(r$RD["hi"]),unname(r$RD["p"])))
  cat(sprintf("  RR  : %.2f (95%%CI %.2f , %.2f)\n", unname(r$RR["est"]),unname(r$RR["lo"]),unname(r$RR["hi"])))
  cat(sprintf("  RMST: %+.2f days (95%%CI %+.2f , %+.2f)  p=%.3f\n",
              unname(r$RMST["est"]),unname(r$RMST["lo"]),unname(r$RMST["hi"]),unname(r$RMST["p"])))
}
show_nu(nu$make30, "MAKE30"); show_nu(nu$make90, "MAKE90")

nu_tbl <- bind_rows(
  data.frame(outcome="MAKE30", n_ppi=nu$n_ppi, n_h2=nu$n_h2, max_smd=nu$max_smd,
             PPI_pct=round(100*nu$make30$rppi,1), H2RA_pct=round(100*nu$make30$rh2,1),
             RD_pct=round(100*unname(nu$make30$RD["est"]),2), RD_lo=round(100*unname(nu$make30$RD["lo"]),2),
             RD_hi=round(100*unname(nu$make30$RD["hi"]),2),
             RR=round(unname(nu$make30$RR["est"]),3), RR_lo=round(unname(nu$make30$RR["lo"]),3), RR_hi=round(unname(nu$make30$RR["hi"]),3),
             RMSTd=round(unname(nu$make30$RMST["est"]),3)),
  data.frame(outcome="MAKE90", n_ppi=nu$n_ppi, n_h2=nu$n_h2, max_smd=nu$max_smd,
             PPI_pct=round(100*nu$make90$rppi,1), H2RA_pct=round(100*nu$make90$rh2,1),
             RD_pct=round(100*unname(nu$make90$RD["est"]),2), RD_lo=round(100*unname(nu$make90$RD["lo"]),2),
             RD_hi=round(100*unname(nu$make90$RD["hi"]),2),
             RR=round(unname(nu$make90$RR["est"]),3), RR_lo=round(unname(nu$make90$RR["lo"]),3), RR_hi=round(unname(nu$make90$RR["hi"]),3),
             RMSTd=round(unname(nu$make90$RMST["est"]),3)))
write.csv(nu_tbl, "newuser_results.csv", row.names=FALSE)
cat("\n---- New-user summary (written to newuser_results.csv) ----\n"); print(nu_tbl, row.names = FALSE)
# Checks: max|SMD| < 0.1 (re-balancing succeeded) + RD/RR/RMST relative to the
# full-cohort primary. A shift toward the null = prevalent-user bias in the
# full-cohort estimate (the IPTW MAKE estimate is already near-null, so the
# shift is modest here; the clearest demonstration is in the CCW analysis).






# =============================================================================
#  R11 -- CCW step 1: weight construction + positivity diagnostics (no outcome yet)
#    IPTW  = PPI vs H2RA head-to-head PS weights (within initiators)
#    IPCW  = initiator vs neither PS weights = 1/P(initiate | L)
#    final = IPTW x IPCW = multinomial ATE weights (stabilized), trimmed at 1/99
#  Purpose: assess positivity (especially initiator vs neither overlap) to decide
#           whether to clone everyone or first apply a common-support trim.
#  Aggregate diagnostics only -> safe to share.
# =============================================================================

hh <- readRDS("imp_long_analysis.rds")   # full-cohort imputations (must include neither)

# ---- key check: data must contain neither + exposure-timing variables ----
cat("==== data self-check (imp 1) ====\n")
cat("exposure_group:\n"); print(table(hh$exposure_group[hh$.imp==1], useNA="ifany"))
need <- c("exposure_group","hours_to_exposure","both_in_window")
miss <- need[!need %in% names(hh)]
if (length(miss)) stop(paste("Missing CCW-required columns:", paste(miss, collapse=", "),
                             "\n-> merge these by id from analytic_cohort, or re-run R1 on the full cohort."))
if (!"neither" %in% hh$exposure_group) stop(
  "No 'neither' in the data! R1 may have imputed only the head-to-head set.\n-> re-run R1 (MICE) on the full analytic_cohort, then do CCW.")

# ---- exclude both-drugs-in-window (SAP section 14.9) ----
n_both <- sum(hh$.imp==1 & as.integer(as.character(hh$both_in_window))==1)
hh <- hh %>% filter(as.integer(as.character(both_in_window))==0)
hh <- hh %>% mutate(initiator = as.integer(exposure_group %in% c("PPI","H2RA")))
cat(sprintf("\nExcluded %d both-drug-in-window; the rest enter cloning.\n", n_both))

m <- max(hh$.imp)

D <- data.frame()
pinit_dec <- NULL
for (i in 1:m) {
  di <- hh %>% filter(.imp==i)
  
  ## IPCW: initiator vs neither (everyone)
  fc <- glm(reformulate(ps_covars,"initiator"), data=di, family=binomial())
  di$pinit <- pmin(pmax(fitted(fc), 1e-6), 1-1e-6)
  pb <- mean(di$initiator)
  di$ipcw <- ifelse(di$initiator==1, pb/di$pinit, (1-pb)/(1-di$pinit))
  
  ## IPTW: PPI vs H2RA (within initiators)
  init <- di %>% filter(initiator==1) %>% mutate(treat=as.integer(exposure_group=="PPI"))
  ft <- glm(reformulate(ps_covars,"treat"), data=init, family=binomial())
  init$pppi <- pmin(pmax(fitted(ft),1e-6),1-1e-6)
  pp <- mean(init$treat)
  init$iptw <- ifelse(init$treat==1, pp/init$pppi, (1-pp)/(1-init$pppi))
  
  ## combine + trim at 1/99
  init$wraw <- init$iptw * init$ipcw
  qlo <- quantile(init$wraw,.01); qhi <- quantile(init$wraw,.99)
  init$wfin <- pmin(pmax(init$wraw,qlo),qhi)
  
  ## balance: weighted PPI/H2RA arms vs full-cohort ATE target (numeric / binary covariates)
  num_covs <- ps_covars[sapply(ps_covars, function(v){ x<-di[[v]]; is.numeric(x)||length(unique(x))==2 })]
  smd_arm <- function(adf) sapply(num_covs, function(v){
    xf<-numv(di[[v]]); xa<-numv(adf[[v]]); s<-sd(xf)
    if(is.na(s)||s==0) NA_real_ else (weighted.mean(xa, adf$wfin)-mean(xf))/s })
  maxsmd <- max(abs(c(smd_arm(init[init$treat==1,]), smd_arm(init[init$treat==0,]))), na.rm=TRUE)
  
  ## positivity: initiator vs neither overlap in P(initiate | L)
  pI <- di$pinit[di$initiator==1]; pN <- di$pinit[di$initiator==0]
  thr <- quantile(pI,.01)                         # 1st percentile of initiator P
  frac_nei_below <- mean(pN < thr)                # fraction of neither below initiator support (extrapolation zone)
  if (i==1) pinit_dec <- data.frame(
    quantile=c("p01","p05","p10","p25","p50","p75","p90","p95","p99"),
    initiator=round(quantile(pI,c(.01,.05,.10,.25,.50,.75,.90,.95,.99)),3),
    neither  =round(quantile(pN,c(.01,.05,.10,.25,.50,.75,.90,.95,.99)),3))
  
  D <- rbind(D, data.frame(imp=i,
                           n_ppi=sum(init$treat==1), n_h2=sum(init$treat==0), n_neither=sum(di$initiator==0),
                           pinit_init_med=median(pI), pinit_nei_med=median(pN), pinit_init_min=min(pI),
                           frac_nei_below=frac_nei_below,
                           ipcw_max=max(di$ipcw[di$initiator==1]), ipcw_p99=as.numeric(quantile(di$ipcw[di$initiator==1],.99)),
                           wfin_max=max(init$wfin), wfin_p99=as.numeric(quantile(init$wfin,.99)),
                           ess_ppi=ess(init$wfin[init$treat==1]), ess_h2=ess(init$wfin[init$treat==0]),
                           eff_ppi=ess(init$wfin[init$treat==1])/sum(init$treat==1),
                           eff_h2 =ess(init$wfin[init$treat==0])/sum(init$treat==0),
                           maxsmd=maxsmd))
}

A <- colMeans(D[,-1])
cat("\n==================== CCW step 1 diagnostics (mean over 20 imputations) ====================\n")
cat(sprintf("Initiators PPI %.0f / H2RA %.0f; neither %.0f\n", A["n_ppi"], A["n_h2"], A["n_neither"]))
cat("\n--- positivity (initiator vs neither) ---\n")
cat(sprintf("Median P(initiate | L): initiators %.3f vs neither %.3f\n", A["pinit_init_med"], A["pinit_nei_med"]))
cat(sprintf("Minimum initiator P(initiate | L) = %.4f\n", A["pinit_init_min"]))
cat(sprintf("Fraction of neither below initiator support (extrapolation zone) = %.1f%%\n", 100*A["frac_nei_below"]))
cat("\nP(initiate | L) quantile comparison (imp 1):\n"); print(pinit_dec, row.names=FALSE)
cat("\n--- weights ---\n")
cat(sprintf("IPCW (initiators) max %.1f / p99 %.1f\n", A["ipcw_max"], A["ipcw_p99"]))
cat(sprintf("Final weight max %.1f / p99 %.1f (trimmed at 1/99)\n", A["wfin_max"], A["wfin_p99"]))
cat(sprintf("Effective sample size ESS: PPI %.0f (efficiency %.0f%%) / H2RA %.0f (efficiency %.0f%%)\n",
            A["ess_ppi"], 100*A["eff_ppi"], A["ess_h2"], 100*A["eff_h2"]))
cat("\n--- balance (combined weights vs full-cohort ATE target) ---\n")
cat(sprintf("Max |SMD| = %.3f (<0.1 desirable)\n", A["maxsmd"]))

write.csv(D, "ccw_step1_diagnostics.csv", row.names=FALSE)
write.csv(pinit_dec, "ccw_step1_pinit_overlap.csv", row.names=FALSE)
# Interpretation:
#  - ESS/n efficiency: >50% excellent; 20-50% usable; <10% caution.
#  - Extrapolation-zone neither fraction: <5% positivity broadly OK (trim a few);
#    >20% means many neither have no comparable initiator -> common-support trim.
#  - wfin p99/max: if still large after trimming (e.g. p99>10) -> heavy tail.
#  - maxsmd<0.1: combined weights achieve ATE balance.






# #############################################################################
#  R12 -- CCW step 2: clone-censor-weight + weighted KM -> MAKE30/90 point est.
#    - combined weight IPTW x IPCW (multinomial ATE, stabilized, 1/99 trimmed)
#    - light common-support trim: drop neither below initiator support
#    - clone everyone into both arms; day0 [0,24h] weight = 1 (shared early
#      experience); censor non-adherent clones at end of day0; from day1
#      adherent clones carry the combined weight
#    - MAKE time-to-event: death at its real time; RRT-dependence / non-recovery
#      placed at the landmark day h
#    - weighted KM -> risk@h (= MAKE_h risk) -> RD/RR; RMST = MAKE-free survival
#      restricted mean
#  This stage gives point estimates only (mean over 20 imputations); CIs (cluster
#  bootstrap) come in R14. Uses the helpers from R11. Aggregate -> safe to share.
# #############################################################################

hh <- readRDS("imp_long_analysis.rds")


mk_tte <- function(d, h, make_var) {
  died <- !is.na(d$days_to_death) & d$days_to_death <= h
  Tk <- ifelse(died, pmin(pmax(d$days_to_death,0)+0.5, h), h)
  mv <- as.integer(as.character(d[[make_var]]))
  ev <- ifelse(died, 1L, ifelse(!is.na(mv) & mv==1, 1L, ifelse(!is.na(mv) & mv==0, 0L, NA)))
  data.frame(Tk=Tk, ev=ev)
}
build_arm <- function(d, a) {
  d <- d %>% mutate(adher = as.integer(exposure_group==a))
  grace <- d %>% transmute(start=0, stop=pmin(Tk,1),
                           event=ifelse(Tk<1, ev, 0L), w=1, arm=a, pid=pid)
  post  <- d %>% filter(adher==1, Tk>=1) %>%
    transmute(start=1, stop=Tk, event=ev, w=wfin, arm=a, pid=pid)
  bind_rows(grace, post)
}
km_summ <- function(tm, sv, h) {
  tt <- c(0, tm); ss <- c(1, sv)
  k <- tt <= h; tt <- tt[k]; ss <- ss[k]
  if (tail(tt,1) < h) { tt <- c(tt,h); ss <- c(ss, tail(ss,1)) }
  K <- length(tt)
  list(risk = 1 - ss[K], rmst = sum(ss[-K]*diff(tt)))
}



hh <- hh %>% filter(as.integer(as.character(both_in_window))==0) %>%
  mutate(initiator = as.integer(exposure_group %in% c("PPI","H2RA")),
         pid = row_number())   # placeholder; rebuilt as a within-.imp row id below

m <- max(hh$.imp)
RES <- data.frame()              # per-imputation point estimates
worst_cov <- character(0)
km_store <- NULL                 # store imp1 KM curves for plotting

for (i in 1:m) {
  di <- hh %>% filter(.imp==i) %>% mutate(pid = row_number())
  
  ## ---- weights ----
  fc <- glm(reformulate(ps_covars,"initiator"), data=di, family=binomial())
  di$pinit <- pmin(pmax(fitted(fc),1e-6),1-1e-6)
  pb <- mean(di$initiator)
  di$ipcw <- ifelse(di$initiator==1, pb/di$pinit, (1-pb)/(1-di$pinit))
  
  ## ---- common-support trim: drop neither below initiator support ----
  thr <- quantile(di$pinit[di$initiator==1], .01)
  drop_id <- di$initiator==0 & di$pinit < thr
  n_drop <- sum(drop_id)
  di <- di[!drop_id,]
  
  init <- di %>% filter(initiator==1) %>% mutate(treat=as.integer(exposure_group=="PPI"))
  ft <- glm(reformulate(ps_covars,"treat"), data=init, family=binomial())
  init$pppi <- pmin(pmax(fitted(ft),1e-6),1-1e-6)
  pp <- mean(init$treat)
  init$iptw <- ifelse(init$treat==1, pp/init$pppi, (1-pp)/(1-init$pppi))
  init$wraw <- init$iptw * init$ipcw
  qlo<-quantile(init$wraw,.01); qhi<-quantile(init$wraw,.99)
  init$wfin <- pmin(pmax(init$wraw,qlo),qhi)
  
  # write wfin back to di (neither and non-adherent opposite-arm clones never use it; set 1 as placeholder)
  di$wfin <- 1
  di$wfin[match(init$pid, di$pid)] <- init$wfin
  
  # record worst-balanced covariates (imp1 only, diagnostic)
  if (i==1) {
    num_covs <- ps_covars[sapply(ps_covars, function(v){x<-di[[v]];is.numeric(x)||length(unique(x))==2})]
    smd1 <- sapply(num_covs, function(v){
      xf<-numv(di[[v]]); s<-sd(xf); if(is.na(s)||s==0) return(NA_real_)
      a1<-(weighted.mean(numv(init[[v]][init$treat==1]), init$wfin[init$treat==1])-mean(xf))/s
      a0<-(weighted.mean(numv(init[[v]][init$treat==0]), init$wfin[init$treat==0])-mean(xf))/s
      max(abs(c(a1,a0))) })
    worst_cov <- names(sort(smd1, decreasing=TRUE))[1:5]
    worst_val <- round(sort(smd1, decreasing=TRUE)[1:5],3)
  }
  
  ## ---- two landmarks h=30, 90 ----
  for (h in c(30,90)) {
    te <- mk_tte(di, h, if (h==30) "make30" else "make90")
    d2 <- di; d2$Tk <- te$Tk; d2$ev <- te$ev
    d2 <- d2[!is.na(d2$ev),]                       # drop non-recovery missing (~2%)
    cl <- bind_rows(build_arm(d2,"PPI"), build_arm(d2,"H2RA"))
    sf <- survfit(Surv(start,stop,event)~arm, data=cl, weights=w)
    # split by strata index
    str_lv <- names(sf$strata); idx <- rep(seq_along(sf$strata), sf$strata)
    ut <- sf$time; us <- sf$surv
    get <- function(armlab) {
      sel <- idx == which(grepl(armlab, str_lv)); km_summ(ut[sel], us[sel], h) }
    P <- get("PPI"); H <- get("H2RA")
    RES <- rbind(RES, data.frame(imp=i, horizon=h,
                                 risk_ppi=P$risk, risk_h2=H$risk, rd=P$risk-H$risk, rr=P$risk/H$risk,
                                 rmst_ppi=P$rmst, rmst_h2=H$rmst, rmst_diff=P$rmst-H$rmst,
                                 n_drop_neither=n_drop,
                                 ess_ppi=ess(init$wfin[init$treat==1]), ess_h2=ess(init$wfin[init$treat==0])))
    # store imp1 curves (h=30 and h=90)
    if (i==1) {
      for (armlab in c("PPI","H2RA")) {
        sel <- idx==which(grepl(armlab,str_lv))
        km_store <- rbind(km_store, data.frame(horizon=h, arm=armlab, time=ut[sel], surv=us[sel]))
      }
    }
  }
}

# ---- pool (mean over 20 imputations) ----
summ <- RES %>% group_by(horizon) %>% summarise(across(c(risk_ppi,risk_h2,rd,rr,
                                                         rmst_ppi,rmst_h2,rmst_diff,ess_ppi,ess_h2,n_drop_neither), mean), .groups="drop")

cat("\n========== CCW per-protocol -- MAKE30/90 point estimates (mean over 20 imputations) ==========\n")
cat(sprintf("Extrapolation neither trimmed: ~%.0f per imputation; effective sample ESS PPI %.0f / H2RA %.0f\n",
            summ$n_drop_neither[1], summ$ess_ppi[1], summ$ess_h2[1]))
cat(sprintf("\nWorst-balanced covariates (imp1): %s\n", paste(sprintf("%s=%.3f", worst_cov, worst_val), collapse=", ")))
for (r in 1:nrow(summ)) {
  s <- summ[r,]
  cat(sprintf("\n----- MAKE%d -----\n", s$horizon))
  cat(sprintf("Weighted cumulative incidence: PPI %.1f%%  H2RA %.1f%%\n", 100*s$risk_ppi, 100*s$risk_h2))
  cat(sprintf("RD = %+.2f%%   RR = %.3f   [>0/>1 = worse with PPI]\n", 100*s$rd, s$rr))
  cat(sprintf("MAKE-free survival RMST: PPI %.2f  H2RA %.2f days   RMST diff = %+.3f days [>0 = better with PPI]\n",
              s$rmst_ppi, s$rmst_h2, s$rmst_diff))
}

write.csv(RES,  "ccw_make_byimp.csv", row.names=FALSE)
write.csv(summ, "ccw_make_pooled_point.csv", row.names=FALSE)
write.csv(km_store, "ccw_km_curve_imp1.csv", row.names=FALSE)
# Checks: are RD/RR still near the null (consistent with the ITT-analog null)?
# risk@30 magnitude should be close to the ITT-analog MAKE30 (~39%) but differs
# because neither / early deaths contribute person-time. The point estimate is
# the registered CCW primary; its CI (cluster bootstrap) comes in R14.







# #############################################################################
#  R13 -- CCW new-user point estimate (the decisive test: is the main-analysis
#         benefit real or prevalent-user bias?)
#    = R12 restricted to newuser == 1, with weights re-estimated within the
#      subset (constant-in-subset covariates dropped). MAKE30/90 point estimates
#      (mean over 20 imputations). Uses the helpers from R11.
#  Aggregate -> safe to share.
# #############################################################################

hh <- readRDS("imp_long_analysis.rds")

# ** new-user restriction **
hh <- hh %>% filter(as.integer(as.character(both_in_window))==0,
                    as.integer(as.character(newuser))==1) %>%
  mutate(initiator = as.integer(exposure_group %in% c("PPI","H2RA")))
cat(sprintf("New-user cohort (imp1): initiators %d, neither %d\n",
            sum(hh$.imp==1 & hh$initiator==1), sum(hh$.imp==1 & hh$initiator==0)))

m <- max(hh$.imp); RES <- data.frame()

for (i in 1:m) {
  di <- hh %>% filter(.imp==i) %>% mutate(pid=row_number())
  kc <- ps_covars[sapply(ps_covars, function(v) length(unique(di[[v]]))>1)]   # drop constant-in-subset columns
  
  fc <- glm(reformulate(kc,"initiator"), data=di, family=binomial())
  di$pinit <- pmin(pmax(fitted(fc),1e-6),1-1e-6)
  pb <- mean(di$initiator)
  di$ipcw <- ifelse(di$initiator==1, pb/di$pinit, (1-pb)/(1-di$pinit))
  
  thr <- quantile(di$pinit[di$initiator==1],.01)
  n_drop <- sum(di$initiator==0 & di$pinit<thr); di <- di[!(di$initiator==0 & di$pinit<thr),]
  
  init <- di %>% filter(initiator==1) %>% mutate(treat=as.integer(exposure_group=="PPI"))
  ft <- glm(reformulate(kc,"treat"), data=init, family=binomial())
  init$pppi <- pmin(pmax(fitted(ft),1e-6),1-1e-6)
  pp <- mean(init$treat)
  init$iptw <- ifelse(init$treat==1, pp/init$pppi, (1-pp)/(1-init$pppi))
  init$wraw <- init$iptw*init$ipcw
  qlo<-quantile(init$wraw,.01); qhi<-quantile(init$wraw,.99); init$wfin<-pmin(pmax(init$wraw,qlo),qhi)
  di$wfin <- 1; di$wfin[match(init$pid, di$pid)] <- init$wfin
  
  for (h in c(30,90)) {
    te <- mk_tte(di,h, if(h==30)"make30" else "make90"); d2<-di; d2$Tk<-te$Tk; d2$ev<-te$ev
    d2 <- d2[!is.na(d2$ev),]
    cl <- bind_rows(build_arm(d2,"PPI"), build_arm(d2,"H2RA"))
    sf <- survfit(Surv(start,stop,event)~arm, data=cl, weights=w)
    str_lv<-names(sf$strata); idx<-rep(seq_along(sf$strata), sf$strata)
    get <- function(a){ sel<-idx==which(grepl(a,str_lv)); km_summ(sf$time[sel], sf$surv[sel], h) }
    P<-get("PPI"); H<-get("H2RA")
    RES <- rbind(RES, data.frame(imp=i, horizon=h,
                                 risk_ppi=P$risk, risk_h2=H$risk, rd=P$risk-H$risk, rr=P$risk/H$risk,
                                 rmst_diff=P$rmst-H$rmst, n_drop=n_drop,
                                 ess_ppi=ess(init$wfin[init$treat==1]), ess_h2=ess(init$wfin[init$treat==0])))
  }
}

summ <- RES %>% group_by(horizon) %>% summarise(across(everything(), mean), .groups="drop") %>% select(-imp)
cat("\n========== New-user CCW -- MAKE30/90 point estimates (mean over 20 imputations) ==========\n")
cat(sprintf("Effective sample ESS PPI %.0f / H2RA %.0f; extrapolation neither trimmed ~%.0f/imputation\n",
            summ$ess_ppi[1], summ$ess_h2[1], summ$n_drop[1]))
for (r in 1:nrow(summ)) { s<-summ[r,]
cat(sprintf("\n----- MAKE%d -----\n", s$horizon))
cat(sprintf("Weighted cumulative incidence: PPI %.1f%%  H2RA %.1f%%\n", 100*s$risk_ppi, 100*s$risk_h2))
cat(sprintf("RD = %+.2f%%   RR = %.3f   RMST diff = %+.3f days\n", 100*s$rd, s$rr, s$rmst_diff))
}
write.csv(RES, "ccw_newuser_byimp.csv", row.names=FALSE)
# Contrast: main CCW MAKE30 RR 0.921 / MAKE90 RR 0.945.
#  - If new-user CCW RR returns to ~1.0 -> the main-analysis benefit was
#    prevalent-user bias; the robust conclusion is null.
#  - If it stays ~0.92 -> a real positive finding.






# #############################################################################
#  R14 (REVISED) -- CCW cluster bootstrap with WEIGHT RE-ESTIMATION per replicate
#  ---------------------------------------------------------------------------
#  WHY THIS REPLACES THE OLD R14
#    The previous R14 fit the PS / IPCW / IPTW ONCE per imputation and held the
#    weights fixed while resampling (multiplier method). That omits the
#    uncertainty of estimating the weights -> the CI is anti-conservative
#    (too narrow), which matters because the whole "primary positive finding"
#    rests on whether the Main-CCW MAKE30 CI excludes 1.
#
#    Here, INSIDE every bootstrap replicate we:
#      (1) resample PATIENTS with replacement (true resample, row duplication;
#          the cluster unit is the original patient),
#      (2) RE-FIT the initiator-vs-neither PS (IPCW) and the PPI-vs-H2RA PS
#          (IPTW) on the resampled data,
#      (3) re-do the common-support trim + 1/99 weight winsorizing,
#      (4) rebuild clones and recompute RD / log-RR / RMST-difference.
#    Bootstrap variance is taken within each imputation; point estimates and
#    within-imputation variances are then pooled by Rubin's rules.
#
#  ESTIMAND / WEIGHTS unchanged from R11-R12:
#    final weight = stabilized IPTW(PPI|init,L) x IPCW(initiate|L), trimmed 1/99;
#    grace block [0,1] carries weight 1 (shared early experience, immortal-time
#    safe); from day 1 only adherent clones continue with the combined weight.
#
#  OUTPUTS (all aggregate -> safe to share):
#    ccw_bootstrap_ci_refit.csv   full results (Main + New-user, MAKE30/90)
#    ccw_ci_sidebyside.csv        prevalent vs new-user, RR point/CI/width/excl-1
#    ccw_ci_width_inflation.csv   refit CI width vs old fixed-weight CI width
#                                 (only if the old ccw_bootstrap_ci.csv exists)
#
#  REQUIRES (already in the session if PART 0 + R1 were run):
#    - ps_covars (PART 0)
#    - imp_long_analysis.rds (full-cohort imputations, INCLUDING 'neither')
#    - columns: both_in_window, newuser, exposure_group, days_to_death,
#               make30, make90, initiator-derivable from exposure_group
#    - packages: survival, dplyr (and parallel if N_CORES > 1)
#
#  RUNTIME: heavy. Start with B = 200 and N_CORES > 1 (Unix/macOS) for a first
#    look; use B = 500-1000 for the final number. With weights re-fit twice per
#    replicate this is the decisive, correct run for the primary inference.
# #############################################################################

suppressPackageStartupMessages({ library(survival); library(dplyr) })

# ---- tunables --------------------------------------------------------------
set.seed(20260610)
B       <- 500   # bootstrap replicates (200 for a quick look; 500-1000 final)
N_CORES <- 1     # >1 parallelizes the bootstrap via mclapply (Unix/macOS only)

stopifnot(exists("ps_covars"))   # defined once in PART 0

# ---- self-contained helpers (so this block runs even if R12 was not sourced)-
mk_tte <- function(d, h, mvn) {
  died <- !is.na(d$days_to_death) & d$days_to_death <= h
  Tk <- ifelse(died, pmin(pmax(d$days_to_death, 0) + 0.5, h), h)
  mv <- as.integer(as.character(d[[mvn]]))
  ev <- ifelse(died, 1L, ifelse(!is.na(mv) & mv == 1, 1L,
                                ifelse(!is.na(mv) & mv == 0, 0L, NA)))
  data.frame(Tk = Tk, ev = ev)
}
km_summ <- function(tm, sv, h) {
  tt <- c(0, tm); ss <- c(1, sv); k <- tt <= h; tt <- tt[k]; ss <- ss[k]
  if (tail(tt, 1) < h) { tt <- c(tt, h); ss <- c(ss, tail(ss, 1)) }
  K <- length(tt); list(risk = 1 - ss[K], rmst = sum(ss[-K] * diff(tt)))
}

# ---- single landmark estimate from a weighted, cloned data set -------------
.est_one <- function(di, h) {
  te <- mk_tte(di, h, if (h == 30) "make30" else "make90")
  d <- di; d$Tk <- te$Tk; d$ev <- te$ev; d <- d[!is.na(d$ev), ]
  mk <- function(a) {
    adher <- as.integer(d$exposure_group == a)
    grace <- data.frame(start = 0, stop = pmin(d$Tk, 1),
                        event = ifelse(d$Tk < 1, d$ev, 0L), w = 1, arm = a)
    pi <- adher == 1 & d$Tk >= 1
    post <- data.frame(start = 1, stop = d$Tk[pi], event = d$ev[pi],
                       w = d$wfin[pi], arm = a)
    rbind(grace, post)
  }
  cl <- rbind(mk("PPI"), mk("H2RA"))
  sf <- survfit(Surv(start, stop, event) ~ arm, data = cl, weights = cl$w)
  lv <- names(sf$strata); idx <- rep(seq_along(sf$strata), sf$strata)
  g <- function(a) { sel <- idx == which(grepl(a, lv)); km_summ(sf$time[sel], sf$surv[sel], h) }
  P <- g("PPI"); H <- g("H2RA")
  if (!is.finite(H$risk) || !is.finite(P$risk) || H$risk <= 0 || P$risk <= 0)
    return(c(NA_real_, NA_real_, NA_real_))
  c(P$risk - H$risk, log(P$risk / H$risk), P$rmst - H$rmst)
}

# ---- fit weights ONCE on `di`, return both landmarks (h = 30 and 90) -------
#      di must contain `initiator` and the PS covariates + outcome columns.
.fit_and_estimate <- function(di) {
  res <- c(rd30 = NA_real_, lrr30 = NA_real_, rmstd30 = NA_real_,
           rd90 = NA_real_, lrr90 = NA_real_, rmstd90 = NA_real_)
  try({
    di$pid <- seq_len(nrow(di))                 # unique clone-unit id per (resampled) row
    kc <- ps_covars[vapply(ps_covars, function(v) length(unique(di[[v]])) > 1, logical(1))]
    ## IPCW: P(initiate | L)
    fc <- glm(reformulate(kc, "initiator"), data = di, family = binomial())
    di$pinit <- pmin(pmax(fitted(fc), 1e-6), 1 - 1e-6)
    pb <- mean(di$initiator)
    di$ipcw <- ifelse(di$initiator == 1, pb / di$pinit, (1 - pb) / (1 - di$pinit))
    ## common-support trim: drop neither below initiator support
    thr <- quantile(di$pinit[di$initiator == 1], 0.01)
    di <- di[!(di$initiator == 0 & di$pinit < thr), ]
    ## IPTW: P(PPI | initiate, L)
    init <- di[di$initiator == 1, ]
    init$treat <- as.integer(init$exposure_group == "PPI")
    if (min(sum(init$treat == 1), sum(init$treat == 0)) < 5) stop("treatment cell < 5")
    ft <- glm(reformulate(kc, "treat"), data = init, family = binomial())
    init$pppi <- pmin(pmax(fitted(ft), 1e-6), 1 - 1e-6)
    pp <- mean(init$treat)
    init$iptw <- ifelse(init$treat == 1, pp / init$pppi, (1 - pp) / (1 - init$pppi))
    ## combined stabilized weight = IPTW x IPCW, winsorized 1/99
    init$wraw <- init$iptw * init$ipcw
    ql <- quantile(init$wraw, 0.01); qh <- quantile(init$wraw, 0.99)
    init$wfin <- pmin(pmax(init$wraw, ql), qh)
    di$wfin <- 1
    di$wfin[match(init$pid, di$pid)] <- init$wfin
    ## both landmarks share the same weights
    e30 <- .est_one(di, 30); e90 <- .est_one(di, 90)
    res[1:3] <- e30; res[4:6] <- e90
  }, silent = TRUE)
  res
}

# ---- one cohort (prevalent-inclusive or new-user) with full refit bootstrap -
run_cohort_refit <- function(newuser_only, B, n_cores = 1) {
  hh <- hh0[as.integer(as.character(hh0$both_in_window)) == 0, ]
  if (newuser_only) hh <- hh[as.integer(as.character(hh$newuser)) == 1, ]
  hh$initiator <- as.integer(hh$exposure_group %in% c("PPI", "H2RA"))
  M <- max(hh$.imp)
  cols <- c("rd30", "lrr30", "rmstd30", "rd90", "lrr90", "rmstd90")
  PT <- matrix(NA_real_, M, 6, dimnames = list(NULL, cols))   # per-imp point estimate
  VR <- matrix(NA_real_, M, 6, dimnames = list(NULL, cols))   # per-imp bootstrap variance
  NF <- numeric(M)
  tag <- if (newuser_only) "new-user" else "main"
  for (i in 1:M) {
    di0 <- hh[hh$.imp == i, ]; n <- nrow(di0)
    PT[i, ] <- .fit_and_estimate(di0)                          # observed point estimate
    one_boot <- function(b) .fit_and_estimate(di0[sample.int(n, n, replace = TRUE), , drop = FALSE])
    bl <- if (n_cores > 1 && .Platform$OS.type == "unix")
      parallel::mclapply(1:B, one_boot, mc.cores = n_cores)
    else lapply(1:B, one_boot)
    bm <- do.call(rbind, bl)
    ok <- stats::complete.cases(bm); NF[i] <- sum(!ok)
    if (sum(ok) >= 2) VR[i, ] <- apply(bm[ok, , drop = FALSE], 2, var)
    cat(sprintf("  [%s] imp %d/%d  (boot fails %d/%d)\n", tag, i, M, NF[i], B))
  }
  ## Rubin pooling (t reference df), per estimand x horizon
  poolt <- function(pt, vr) {
    keep <- is.finite(pt) & is.finite(vr); pt <- pt[keep]; vr <- vr[keep]; m <- length(pt)
    qb <- mean(pt); Ub <- mean(vr); Bb <- var(pt)
    Tt <- Ub + (1 + 1 / m) * Bb; se <- sqrt(Tt)
    df <- if (Bb > 0) (m - 1) * (1 + Ub / ((1 + 1 / m) * Bb))^2 else 1e6
    tc <- qt(0.975, df); c(est = qb, lo = qb - tc * se, hi = qb + tc * se, se = se)
  }
  mk_row <- function(h) {
    s <- as.character(h)
    rd <- poolt(PT[, paste0("rd", s)],    VR[, paste0("rd", s)])
    rr <- poolt(PT[, paste0("lrr", s)],   VR[, paste0("lrr", s)])
    rm <- poolt(PT[, paste0("rmstd", s)], VR[, paste0("rmstd", s)])
    data.frame(
      cohort = if (newuser_only) "New-user CCW" else "Main CCW (prevalent-incl.)",
      horizon = h,
      RD_pct = 100 * rd["est"], RD_lo = 100 * rd["lo"], RD_hi = 100 * rd["hi"],
      RR = exp(rr["est"]), RR_lo = exp(rr["lo"]), RR_hi = exp(rr["hi"]),
      RMSTd = rm["est"], RMSTd_lo = rm["lo"], RMSTd_hi = rm["hi"],
      mean_boot_fail = round(mean(NF), 1), row.names = NULL)
  }
  rbind(mk_row(30), mk_row(90))
}

# ===========================================================================
#  RUN
# ===========================================================================
hh0 <- as.data.frame(readRDS("imp_long_analysis.rds"))
if (!"neither" %in% hh0$exposure_group)
  stop("No 'neither' in imp_long_analysis.rds -> re-run R1 (MICE) on the FULL cohort.")

t0 <- Sys.time()
cat("Running MAIN CCW (prevalent-inclusive), refit bootstrap ...\n")
main <- run_cohort_refit(FALSE, B, N_CORES)
cat("Running NEW-USER CCW, refit bootstrap ...\n")
nu   <- run_cohort_refit(TRUE,  B, N_CORES)
res  <- rbind(main, nu)
cat(sprintf("\nElapsed: %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))

# ---- add CI widths + null-exclusion flags ---------------------------------
res2 <- res %>% mutate(
  RD_width   = round(RD_hi   - RD_lo, 2),
  RR_width   = round(RR_hi   - RR_lo, 3),
  RMST_width = round(RMSTd_hi - RMSTd_lo, 3),
  RR_excl1   = (RR_lo > 1)   | (RR_hi < 1),
  RD_excl0   = (RD_lo > 0)   | (RD_hi < 0),
  RMST_excl0 = (RMSTd_lo > 0) | (RMSTd_hi < 0)
) %>%
  mutate(across(c(RD_pct, RD_lo, RD_hi), ~round(., 2)),
         across(c(RR, RR_lo, RR_hi), ~round(., 3)),
         across(c(RMSTd, RMSTd_lo, RMSTd_hi), ~round(., 3)))

write.csv(res2, "ccw_bootstrap_ci_refit.csv", row.names = FALSE)
cat("\n================ CCW refit-bootstrap results (Main vs New-user) ================\n")
print(res2, row.names = FALSE)

# ---- prevalent vs new-user, side by side (headline RR) --------------------
mrow <- res2[res2$cohort == "Main CCW (prevalent-incl.)", ]
nrow_ <- res2[res2$cohort == "New-user CCW", ]
cmp <- data.frame(
  horizon          = mrow$horizon,
  Main_RR          = sprintf("%.3f (%.3f-%.3f)", mrow$RR, mrow$RR_lo, mrow$RR_hi),
  Main_RRwidth     = mrow$RR_width,
  Main_excl1       = mrow$RR_excl1,
  NewUser_RR       = sprintf("%.3f (%.3f-%.3f)", nrow_$RR, nrow_$RR_lo, nrow_$RR_hi),
  NewUser_RRwidth  = nrow_$RR_width,
  NewUser_excl1    = nrow_$RR_excl1
)
write.csv(cmp, "ccw_ci_sidebyside.csv", row.names = FALSE)
cat("\n---- Prevalent-inclusive vs New-user (MAKE RR, point / CI / width / excludes 1) ----\n")
print(cmp, row.names = FALSE)

# ---- CI-width inflation vs the OLD fixed-weight bootstrap (if available) ---
old <- tryCatch(read.csv("ccw_bootstrap_ci.csv"), error = function(e) NULL)
if (!is.null(old) && all(c("RR_lo", "RR_hi", "cohort", "horizon") %in% names(old))) {
  old$key <- ifelse(grepl("New", old$cohort), "newuser", "main")
  old$RR_width_fixed <- round(old$RR_hi - old$RR_lo, 3)
  res2$key <- ifelse(grepl("New", res2$cohort), "newuser", "main")
  infl <- merge(res2[, c("key", "horizon", "RR_width")],
                old[, c("key", "horizon", "RR_width_fixed")],
                by = c("key", "horizon"))
  names(infl)[names(infl) == "RR_width"] <- "RR_width_refit"
  infl$width_ratio_refit_vs_fixed <- round(infl$RR_width_refit / infl$RR_width_fixed, 2)
  write.csv(infl, "ccw_ci_width_inflation.csv", row.names = FALSE)
  cat("\n---- CI-width inflation: refit (correct) vs old fixed-weight bootstrap ----\n")
  print(infl[order(infl$key, infl$horizon), ], row.names = FALSE)
} else {
  cat("\n(No old ccw_bootstrap_ci.csv found -> skipping width-inflation comparison.)\n")
}

# ===========================================================================
#  DECISION RULE (prints the branch you are on)
# ===========================================================================
main30 <- res2[res2$cohort == "Main CCW (prevalent-incl.)" & res2$horizon == 30, ]
cat("\n================ INTERPRETATION ================\n")
if (!isTRUE(main30$RR_excl1)) {
  cat("BRANCH 1: with correctly propagated weight uncertainty, the Main-CCW MAKE30\n")
  cat(sprintf("  RR = %.3f (95%% CI %.3f-%.3f) NO LONGER excludes 1.\n",
              main30$RR, main30$RR_lo, main30$RR_hi))
  cat("  -> There is no 'positive' primary finding to explain. Report the registered\n")
  cat("     primary (prevalent-inclusive CCW) as null; the ITT-analog and new-user\n")
  cat("     analyses are concordant. The dilemma dissolves.\n")
} else {
  cat(sprintf("BRANCH 2: Main-CCW MAKE30 RR = %.3f (95%% CI %.3f-%.3f) still excludes 1.\n",
              main30$RR, main30$RR_lo, main30$RR_hi))
  cat("  -> Report the registered primary as found, then use the new-user attenuation\n")
  cat("     + clean negative controls + ITT-analog null as triangulation: the lone\n")
  cat("     signal is the analysis most exposed to prevalent-user bias and vanishes\n")
  cat("     once that bias is removed. Compare Main vs New-user POINT estimates (bias)\n")
  cat("     and CI WIDTHS (power) in ccw_ci_sidebyside.csv before concluding.\n")
}
cat("\nFiles written: ccw_bootstrap_ci_refit.csv, ccw_ci_sidebyside.csv",
    if (!is.null(old)) ", ccw_ci_width_inflation.csv" else "", "\n", sep = "")






# #############################################################################
#  R15 -- SUPPLEMENTARY ANALYSES TO CHARACTERIZE THE PREVALENT-USER BIAS
#  ---------------------------------------------------------------------------
#  Append this AFTER R14 (revised). It reuses R14's functions:
#     .fit_and_estimate(), .est_one(), mk_tte(), km_summ()
#  and the same weight construction (IPTW x IPCW, common-support trim, 1/99
#  winsorize), so every number here is produced by the SAME machinery that
#  produced the registered primary RR = 0.921.
#
#  Three pieces, each answering one reviewer question:
#   (15A) E-VALUES on the CCW primary RRs
#         -> "How fragile is the protective signal to residual confounding/
#            selection?"  (small E-value = easily explained by prevalent-user
#            differences.)
#   (15B) NEGATIVE CONTROLS (incident VTE primary; pressure injury secondary)
#         run through the MAIN CCW pipeline, with new-user CCW for contrast
#         -> "Does the prevalent-inclusive pipeline manufacture a spurious
#            association on an outcome with no causal link to PPI vs H2RA?"
#   (15C) PREVALENT-ONLY vs ALL vs NEW-USER CCW (bias gradient)
#         -> "Is the signal concentrated in prevalent users and absent in
#            incident users (dose-response of the bias)?"
#
#  DEPENDENCIES (in session after R14): .fit_and_estimate, .est_one, mk_tte,
#     km_summ, ps_covars, and imp_long_analysis.rds (full cohort incl neither).
#  OUTPUTS (all aggregate -> shareable):
#     ccw_primary_evalues.csv, ccw_negative_controls.csv, ccw_bias_gradient.csv
# #############################################################################

suppressPackageStartupMessages(library(survival))
stopifnot(exists(".fit_and_estimate"), exists("ps_covars"))

set.seed(20260611)
B_GRAD <- 300    # bootstrap reps for the prevalent-only MAKE run (survfit; heavier)
B_NCO  <- 500    # bootstrap reps for the binary negative controls (cheap)
N_CORES <- 1     # >1 parallelizes (mclapply; Unix/macOS only)

hh0 <- as.data.frame(readRDS("imp_long_analysis.rds"))
if (!"neither" %in% hh0$exposure_group)
  stop("imp_long_analysis.rds has no 'neither' -> re-run R1 on the FULL cohort.")

## shared Rubin pooler (t reference df) --------------------------------------
.poolt <- function(pt, vr) {
  k <- is.finite(pt) & is.finite(vr); pt <- pt[k]; vr <- vr[k]; m <- length(pt)
  qb <- mean(pt); Ub <- mean(vr); Bb <- var(pt)
  Tt <- Ub + (1 + 1/m) * Bb; se <- sqrt(Tt)
  df <- if (Bb > 0) (m - 1) * (1 + Ub / ((1 + 1/m) * Bb))^2 else 1e6
  tc <- qt(.975, df); c(est = qb, lo = qb - tc * se, hi = qb + tc * se)
}

# ===========================================================================
#  15A.  E-VALUES on the CCW primary risk ratios
#        (VanderWeele & Ding 2017). For RR<1 use RR*=1/RR.
# ===========================================================================
ev_one <- function(x) { x <- ifelse(x < 1, 1 / x, x); x + sqrt(x * (x - 1)) }
ev_rr <- function(est, lo, hi) {
  Ep <- ev_one(est)
  if (est < 1)      Eci <- if (hi >= 1) 1 else ev_one(hi)   # bound nearer null = upper
  else if (est > 1) Eci <- if (lo <= 1) 1 else ev_one(lo)   # bound nearer null = lower
  else              Eci <- 1
  c(E_point = round(Ep, 3), E_ci = round(Eci, 3))
}

ci <- tryCatch(read.csv("ccw_bootstrap_ci_refit.csv"), error = function(e) NULL)
if (is.null(ci) && exists("res2")) ci <- res2
if (is.null(ci)) stop("Run R14 first (need ccw_bootstrap_ci_refit.csv or res2).")

ev <- do.call(rbind, lapply(seq_len(nrow(ci)), function(j) {
  e <- ev_rr(ci$RR[j], ci$RR_lo[j], ci$RR_hi[j])
  data.frame(cohort = ci$cohort[j], horizon = ci$horizon[j],
             RR = ci$RR[j], RR_lo = ci$RR_lo[j], RR_hi = ci$RR_hi[j],
             E_point = e["E_point"], E_ci = e["E_ci"], row.names = NULL)
}))
write.csv(ev, "ccw_primary_evalues.csv", row.names = FALSE)
cat("\n================ 15A. E-values on CCW primary RRs ================\n")
print(ev, row.names = FALSE)
cat("\nReading: an unmeasured confounder/selection associated with BOTH exposure and\n",
    "outcome by RR >= E_point would explain the point estimate away; >= E_ci nullifies\n",
    "the CI. Small E-values => the signal is fragile and easily produced by the way\n",
    "prevalent users differ from incident users.\n", sep = "")

# ===========================================================================
#  15B.  NEGATIVE CONTROLS through the MAIN CCW pipeline
# ===========================================================================
## verbatim copy of R14's weight construction (so weights are identical) -----
.ccw_weights <- function(di) {
  di$pid <- seq_len(nrow(di))
  kc <- ps_covars[vapply(ps_covars, function(v) length(unique(di[[v]])) > 1, logical(1))]
  fc <- glm(reformulate(kc, "initiator"), data = di, family = binomial())
  di$pinit <- pmin(pmax(fitted(fc), 1e-6), 1 - 1e-6)
  pb <- mean(di$initiator)
  di$ipcw <- ifelse(di$initiator == 1, pb / di$pinit, (1 - pb) / (1 - di$pinit))
  thr <- quantile(di$pinit[di$initiator == 1], 0.01)
  di <- di[!(di$initiator == 0 & di$pinit < thr), ]
  init <- di[di$initiator == 1, ]; init$treat <- as.integer(init$exposure_group == "PPI")
  if (min(sum(init$treat == 1), sum(init$treat == 0)) < 5) stop("treatment cell < 5")
  ft <- glm(reformulate(kc, "treat"), data = init, family = binomial())
  init$pppi <- pmin(pmax(fitted(ft), 1e-6), 1 - 1e-6)
  pp <- mean(init$treat)
  init$iptw <- ifelse(init$treat == 1, pp / init$pppi, (1 - pp) / (1 - init$pppi))
  init$wraw <- init$iptw * init$ipcw
  ql <- quantile(init$wraw, 0.01); qh <- quantile(init$wraw, 0.99)
  init$wfin <- pmin(pmax(init$wraw, ql), qh)
  di$wfin <- 1; di$wfin[match(init$pid, di$pid)] <- init$wfin
  di
}
## weighted binary RR/RD of a (timeless) NCO among the per-protocol initiators-
.nco_estimate <- function(di, nco) {
  out <- c(rd = NA_real_, lrr = NA_real_)
  try({
    di <- .ccw_weights(di)
    y  <- as.integer(as.character(di[[nco]]))
    ii <- di$initiator == 1 & !is.na(y)
    init <- di[ii, ]; init$y <- y[ii]; init$treat <- as.integer(init$exposure_group == "PPI")
    if (min(sum(init$treat == 1), sum(init$treat == 0)) < 5) stop("cell < 5")
    wr <- function(a) { s <- init[init$treat == a, ]; sum(s$wfin * s$y) / sum(s$wfin) }
    rP <- wr(1); rH <- wr(0)
    if (rP <= 0 || rH <= 0) stop("zero risk")
    out <- c(rd = rP - rH, lrr = log(rP / rH))
  }, silent = TRUE)
  out
}
## one NCO x one population, with refit cluster bootstrap --------------------
ccw_run_nco <- function(nco, pop, B, n_cores = 1) {
  hh <- hh0[as.integer(as.character(hh0$both_in_window)) == 0, ]
  if (pop == "newuser") hh <- hh[as.integer(as.character(hh$newuser)) == 1, ]
  hh$initiator <- as.integer(hh$exposure_group %in% c("PPI", "H2RA"))
  M <- max(hh$.imp)
  PT <- matrix(NA_real_, M, 2, dimnames = list(NULL, c("rd", "lrr")))
  VR <- matrix(NA_real_, M, 2, dimnames = list(NULL, c("rd", "lrr"))); NF <- numeric(M)
  for (i in 1:M) {
    di0 <- hh[hh$.imp == i, ]; n <- nrow(di0)
    PT[i, ] <- .nco_estimate(di0, nco)
    ob <- function(b) .nco_estimate(di0[sample.int(n, n, replace = TRUE), , drop = FALSE], nco)
    bl <- if (n_cores > 1 && .Platform$OS.type == "unix")
      parallel::mclapply(1:B, ob, mc.cores = n_cores) else lapply(1:B, ob)
    bm <- do.call(rbind, bl); ok <- stats::complete.cases(bm); NF[i] <- sum(!ok)
    if (sum(ok) >= 2) VR[i, ] <- apply(bm[ok, , drop = FALSE], 2, var)
    cat(sprintf("  [NCO %s/%s] imp %d/%d (fails %d/%d)\n", nco, pop, i, M, NF[i], B))
  }
  rr <- .poolt(PT[, "lrr"], VR[, "lrr"]); rd <- .poolt(PT[, "rd"], VR[, "rd"])
  data.frame(nco = nco, pop = pop,
             RR = round(exp(rr["est"]), 3), RR_lo = round(exp(rr["lo"]), 3), RR_hi = round(exp(rr["hi"]), 3),
             RD_pct = round(100 * rd["est"], 2), RD_lo = round(100 * rd["lo"], 2), RD_hi = round(100 * rd["hi"], 2),
             excl1 = (exp(rr["lo"]) > 1) | (exp(rr["hi"]) < 1), row.names = NULL)
}

## resolve NCO column names (adjust here if your analytic table differs) -----
pick <- function(cands) { hit <- cands[cands %in% names(hh0)]; if (length(hit) == 0) NA_character_ else hit[1] }
vte_col  <- pick(c("vte_incident", "nco_vte", "incident_vte", "vte"))
pinj_col <- pick(c("pressure_injury", "nco_pressure_injury", "pressure_injury_incident", "pinj"))
cat(sprintf("\nNCO columns resolved -> VTE: %s | pressure injury: %s\n",
            ifelse(is.na(vte_col), "NOT FOUND", vte_col),
            ifelse(is.na(pinj_col), "NOT FOUND", pinj_col)))
if (is.na(vte_col))
  cat("  (no VTE column matched; columns containing 'vte':",
      paste(grep("vte", names(hh0), value = TRUE, ignore.case = TRUE), collapse = ", "),
      ")\n  -> set vte_col manually and re-run 15B.\n")

cat("\n================ 15B. Negative controls through MAIN CCW pipeline ================\n")
nco_rows <- list()
for (nm in c(vte_col, pinj_col)) {
  if (is.na(nm)) next
  nco_rows[[paste0(nm, "_all")]]     <- ccw_run_nco(nm, "all",     B_NCO, N_CORES)
  nco_rows[[paste0(nm, "_newuser")]] <- ccw_run_nco(nm, "newuser", B_NCO, N_CORES)
}
if (length(nco_rows)) {
  nco_tab <- do.call(rbind, nco_rows)
  write.csv(nco_tab, "ccw_negative_controls.csv", row.names = FALSE)
  print(nco_tab, row.names = FALSE)
  cat("\nReading: under a valid negative control the truth is RR = 1. If the MAIN (all)\n",
      "CCW pushes the NCO away from 1 -- especially in the SAME protective direction as\n",
      "MAKE -- that is DIRECT evidence the prevalent-inclusive pipeline carries residual\n",
      "bias. A null NCO under new-user that mirrors the MAKE pattern strengthens the\n",
      "prevalent-user-bias interpretation. (Use these for negative-control calibration.)\n",
      sep = "")
} else {
  cat("No NCO columns found -> skipped. Set vte_col / pinj_col manually.\n")
}

# ===========================================================================
#  15C.  PREVALENT-ONLY vs ALL vs NEW-USER CCW (bias gradient)
# ===========================================================================
ccw_run_pop <- function(pop, B, n_cores = 1) {
  hh <- hh0[as.integer(as.character(hh0$both_in_window)) == 0, ]
  if (pop == "newuser")   hh <- hh[as.integer(as.character(hh$newuser)) == 1, ]
  if (pop == "prevalent") hh <- hh[as.integer(as.character(hh$newuser)) == 0, ]
  hh$initiator <- as.integer(hh$exposure_group %in% c("PPI", "H2RA"))
  M <- max(hh$.imp)
  cols <- c("rd30", "lrr30", "rmstd30", "rd90", "lrr90", "rmstd90")
  PT <- matrix(NA_real_, M, 6, dimnames = list(NULL, cols))
  VR <- matrix(NA_real_, M, 6, dimnames = list(NULL, cols)); NF <- numeric(M)
  for (i in 1:M) {
    di0 <- hh[hh$.imp == i, ]; n <- nrow(di0)
    PT[i, ] <- .fit_and_estimate(di0)
    ob <- function(b) .fit_and_estimate(di0[sample.int(n, n, replace = TRUE), , drop = FALSE])
    bl <- if (n_cores > 1 && .Platform$OS.type == "unix")
      parallel::mclapply(1:B, ob, mc.cores = n_cores) else lapply(1:B, ob)
    bm <- do.call(rbind, bl); ok <- stats::complete.cases(bm); NF[i] <- sum(!ok)
    if (sum(ok) >= 2) VR[i, ] <- apply(bm[ok, , drop = FALSE], 2, var)
    cat(sprintf("  [%s] imp %d/%d (fails %d/%d)\n", pop, i, M, NF[i], B))
  }
  row_h <- function(h) {
    s <- as.character(h)
    rr <- .poolt(PT[, paste0("lrr", s)], VR[, paste0("lrr", s)])
    rd <- .poolt(PT[, paste0("rd", s)],  VR[, paste0("rd", s)])
    data.frame(population = pop, horizon = h,
               RR = round(exp(rr["est"]), 3), RR_lo = round(exp(rr["lo"]), 3), RR_hi = round(exp(rr["hi"]), 3),
               RD_pct = round(100 * rd["est"], 2), RD_lo = round(100 * rd["lo"], 2), RD_hi = round(100 * rd["hi"], 2),
               RR_excl1 = (exp(rr["lo"]) > 1) | (exp(rr["hi"]) < 1), row.names = NULL)
  }
  rbind(row_h(30), row_h(90))
}

cat("\n================ 15C. Bias gradient: prevalent-only vs all vs new-user ================\n")
prev_tab <- ccw_run_pop("prevalent", B_GRAD, N_CORES)

## pull ALL and NEW-USER rows from the R14 output for a single gradient table
ci2 <- tryCatch(read.csv("ccw_bootstrap_ci_refit.csv"), error = function(e) NULL)
to_pop <- function(df, lab) data.frame(
  population = lab, horizon = df$horizon,
  RR = df$RR, RR_lo = df$RR_lo, RR_hi = df$RR_hi,
  RD_pct = df$RD_pct, RD_lo = df$RD_lo, RD_hi = df$RD_hi,
  RR_excl1 = (df$RR_lo > 1) | (df$RR_hi < 1))
grad <- prev_tab
if (!is.null(ci2)) {
  all_rows <- to_pop(ci2[grepl("Main", ci2$cohort), ], "all (prevalent+incident)")
  nu_rows  <- to_pop(ci2[grepl("New",  ci2$cohort), ], "new-user (incident only)")
  grad <- rbind(prev_tab, all_rows, nu_rows)
}
grad$population <- factor(grad$population,
                          levels = c("prevalent", "all (prevalent+incident)", "new-user (incident only)"))
grad <- grad[order(grad$horizon, grad$population), ]
write.csv(grad, "ccw_bias_gradient.csv", row.names = FALSE)
print(grad, row.names = FALSE)
cat("\nReading: the expected bias dose-response is\n",
    "  prevalent-only RR < all RR (0.921) < new-user RR (~1.00).\n",
    "If prevalent-only is the most protective and new-user sits at the null, the MAKE30\n",
    "signal is concentrated exactly where prevalent-user bias lives and absent where the\n",
    "design matches the T0 decision -- the cleanest single piece of evidence for your\n",
    "'prevalent-user bias' interpretation.\n", sep = "")

cat("\nFiles written: ccw_primary_evalues.csv, ccw_negative_controls.csv, ccw_bias_gradient.csv\n")






# ===========================================================================
#   R16.  BH-FDR over the pre-specified SECONDARY endpoints (SAP 13.4)
#  ---------------------------------------------------------------------------
#  SET (locked BEFORE computing q; state this in the manuscript):
#    5 pre-specified secondary endpoints, ONE p each, from the PRE-SPECIFIED
#    primary estimator and at its primary horizon:
#       kidney recovery        -> Fine-Gray sHR p, 30-day   [recovery_30 / sHR_p]
#       progression to KDIGO3  -> cause-specific HR p, 30-day [progression_30 / csHR_p]
#       RRT initiation         -> cause-specific HR p, 30-day [RRT_30 / csHR_p]
#       28-day all-cause death -> RMST-difference p           [death_28 / RMSTd_p]
#       90-day all-cause death -> RMST-difference p           [death_90 / RMSTd_p]
#    Recovery/progression/RRT are in-hospital endpoints, more complete at day 30;
#    the 90-day versions are SUPPORTIVE and are NOT in the BH set. Both mortality
#    horizons are included because SAP 13.4 lists 28- and 90-day mortality as
#    distinct secondary endpoints.
#
#    EXCLUDED from BH:
#       MAKE30 / MAKE90        (primary & co-primary; nested, no adjustment)
#       negative controls       (diagnostic, not hypothesis tests)
#       subgroup interactions   (hypothesis-generating, SAP 13.5)
#       sensitivity analyses    (interpreted qualitatively)
#       alternate-horizon recovery/progression/RRT (supportive)
# ===========================================================================

# Pull one Rubin-pooled p-value from an R5 output file (errors loudly if absent)
grab <- function(file, row_outcome, p_col) {
  d <- tryCatch(read.csv(file), error = function(e)
    stop(sprintf("Cannot open %s (run R5 first).", file)))
  if (!"outcome" %in% names(d)) stop(sprintf("No 'outcome' column in %s.", file))
  if (!p_col %in% names(d))     stop(sprintf("No '%s' column in %s. Columns: %s",
                                             p_col, file, paste(names(d), collapse = ", ")))
  v <- d[d$outcome == row_outcome, p_col]
  if (length(v) != 1 || is.na(v))
    stop(sprintf("p not found for outcome '%s' in %s [%s].", row_outcome, file, p_col))
  as.numeric(v)
}

sec <- data.frame(
  endpoint  = c("Kidney recovery (Fine-Gray sHR, 30d)",
                "Progression to KDIGO 3 (cause-specific HR, 30d)",
                "RRT initiation (cause-specific HR, 30d)",
                "28-day all-cause mortality (RMST diff)",
                "90-day all-cause mortality (RMST diff)"),
  estimator = c("sHR", "csHR", "csHR", "RMSTd", "RMSTd"),
  horizon   = c(30, 30, 30, 28, 90),
  p = c(
    grab("secondary_recovery.csv",    "recovery_30",    "sHR_p"),
    grab("secondary_progression.csv", "progression_30", "csHR_p"),
    grab("secondary_rrt.csv",         "RRT_30",         "csHR_p"),
    grab("secondary_mortality.csv",   "death_28",       "RMSTd_p"),
    grab("secondary_mortality.csv",   "death_90",       "RMSTd_p")
  ),
  stringsAsFactors = FALSE
)

sec$q_BH <- p.adjust(sec$p, method = "BH")
sec$significant_q0.10 <- sec$q_BH < 0.10
sec <- sec[order(sec$p), ]

write.csv(sec, "secondary_BH_FDR.csv", row.names = FALSE)
cat("\n===== 16A. Secondary endpoints -- Benjamini-Hochberg FDR (q<0.10), SAP 13.4 =====\n")
print(sec, row.names = FALSE, digits = 4)
cat("\nSet = 5 pre-specified secondary endpoints, one p each at its primary horizon\n",
    "(recovery/progression/RRT at 30d; 28- and 90-day mortality). EXCLUDED: MAKE30/\n",
    "MAKE90 (primary/co-primary), negative controls, subgroup interactions, sensitivity\n",
    "analyses, and the supportive 90-day recovery/progression/RRT results.\n", sep = "")

