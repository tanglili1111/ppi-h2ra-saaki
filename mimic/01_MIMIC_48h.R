# #############################################################################
#  PPI vs H2RA in SA-AKI -- SENSITIVITY ANALYSIS: exposure window / CCW grace = 48h
#  ---------------------------------------------------------------------------
#  Pre-specified supplementary analysis (SAP 14, item 6: 48-72h window).
#  Role: show that BOTH the primary null AND the prevalent-user-bias signature
#  (main CCW protective -> new-user CCW null) are robust to widening the grace
#  period from 24h to 48h. A wider window is also more forgiving of the coarse
#  prescriptions.starttime timing, partly addressing that data-quality limitation.
#
#  *** METHOD MATCHES THE CORRECTED 24h MAIN ANALYSIS (R14 REVISED) ***
#    The CCW cluster bootstrap RE-ESTIMATES the PS / IPCW / IPTW inside every
#    resample (patient-level resampling with row duplication), so the weight-
#    estimation uncertainty is propagated -- identical convention to 24h R14.
#
#  Structure
#    PART 1 (R0+R1): data prep + MICE m=20 on analytic_cohort_48h (skipped if RDS exists)
#    PART 2 (R2/R4/R10): IPTW ITT-analog MAKE30/90, full + new-user (svyglm analytic SE)
#    PART 3 (R11-R14): CCW point estimate + refit cluster-bootstrap CI, main + new-user
#  Run: source top to bottom. MICE and the refit bootstrap are both slow.
#    Set N_CORES > 1 on Unix/macOS to parallelize the bootstrap.
#  Prereq: build analytic_cohort_48h in PostgreSQL (same SQL pipeline as 24h with
#    the exposure window = 48h); Sys.setenv(PGPASSWORD="..."); R >= 4.3.
#  DUA: runs locally; only aggregate CSVs leave the environment.
# #############################################################################

# =============================================================================
# PART 0 -- environment, libraries, helpers (CCW refit bootstrap, grace GD=2)
# =============================================================================
Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
remove(list = ls())

library(DBI); library(RPostgres); library(dplyr); library(mice); library(survival)
library(WeightIt); library(cobalt); library(survey)   # for the IPTW (ITT-analog) part

path_imp_48 <- "imp_long_analysis_48h.rds"   # single constant for save/read
GD          <- 2                              # grace period = 2 days (48h window)
B_BOOT      <- 500                            # bootstrap replicates (200 for a quick look)
N_CORES     <- 1                              # >1 parallelizes bootstrap (mclapply; Unix/macOS)
set.seed(20260610)

ps_covars <- c(
  "age_at_icu","gender","race_cat","admission_type","first_careunit_cat","anchor_year_group",
  "weight_kg","baseline_egfr","first_aki_stage","ckd","sofa_total","lactate","mech_vent",
  "vasopressor","corticosteroid","albumin","bilirubin_total","inr","platelet","wbc","sodium",
  "potassium","bicarbonate","anion_gap","magnesium","hemoglobin","fluid_balance_ml","uo_rate_ml_kg_h",
  "myocardial_infarct","congestive_heart_failure","peripheral_vascular_disease","cerebrovascular_disease",
  "dementia","chronic_pulmonary_disease","rheumatic_disease","peptic_ulcer_disease","mild_liver_disease",
  "severe_liver_disease","diabetes_without_cc","diabetes_with_cc","paraplegia","renal_disease",
  "malignant_cancer","metastatic_solid_tumor","aids","elixhauser_score","hypertension","gi_bleed_dx",
  "rbc_transfusion_pre_t0","dnr_comfort","nephrotoxin_any","nephro_vanco","nephro_acei_arb",
  "nephro_nsaid","enteral_nutrition","antiplatelet")

## -- MAKE time-to-event (death at real time; non-fatal components at landmark h) --
mk_tte <- function(d, h, mvn) {
  died <- !is.na(d$days_to_death) & d$days_to_death <= h
  Tk <- ifelse(died, pmin(pmax(d$days_to_death, 0) + 0.5, h), h)
  mv <- as.integer(as.character(d[[mvn]]))
  ev <- ifelse(died, 1L, ifelse(!is.na(mv) & mv == 1, 1L, ifelse(!is.na(mv) & mv == 0, 0L, NA)))
  data.frame(Tk = Tk, ev = ev)
}
## -- risk@h and RMST(0,h) from a weighted KM --
km_summ <- function(tm, sv, h) {
  tt <- c(0, tm); ss <- c(1, sv); k <- tt <= h; tt <- tt[k]; ss <- ss[k]
  if (tail(tt, 1) < h) { tt <- c(tt, h); ss <- c(ss, tail(ss, 1)) }
  K <- length(tt); list(risk = 1 - ss[K], rmst = sum(ss[-K] * diff(tt)))
}
## -- Rubin pooling with t reference df (matches corrected 24h R14) --
.poolt <- function(pt, vr) {
  k <- is.finite(pt) & is.finite(vr); pt <- pt[k]; vr <- vr[k]; m <- length(pt)
  qb <- mean(pt); Ub <- mean(vr); Bb <- var(pt)
  Tt <- Ub + (1 + 1/m) * Bb; se <- sqrt(Tt)
  df <- if (Bb > 0) (m - 1) * (1 + Ub / ((1 + 1/m) * Bb))^2 else 1e6
  tc <- qt(.975, df); c(est = qb, lo = qb - tc * se, hi = qb + tc * se)
}
## -- Rubin pooling for IPTW effects (t-df) + RR exp helper (matches main R4/R10) --
pool_rubin <- function(est, v) {
  m <- length(est); qbar <- mean(est); ubar <- mean(v); b <- var(est)
  tvar <- ubar + (1 + 1/m) * b; se <- sqrt(tvar)
  df <- if (b > 0) (m - 1) * (1 + ubar / ((1 + 1/m) * b))^2 else 1e6
  c(est = qbar, se = se, lo = qbar - qt(.975, df) * se, hi = qbar + qt(.975, df) * se,
    p = 2 * pt(-abs(qbar / se), df))
}
expci <- function(p) c(est = exp(unname(p["est"])), lo = exp(unname(p["lo"])),
                       hi = exp(unname(p["hi"])), p = unname(p["p"]))

## -- one landmark estimate from a weighted, cloned data set (grace = gd) --
.est_one48 <- function(di, h, gd) {
  te <- mk_tte(di, h, if (h == 30) "make30" else "make90")
  d <- di; d$Tk <- te$Tk; d$ev <- te$ev; d <- d[!is.na(d$ev), ]
  mk <- function(a) {
    adher <- as.integer(d$exposure_group == a)
    grace <- data.frame(start = 0, stop = pmin(d$Tk, gd),
                        event = ifelse(d$Tk < gd, d$ev, 0L), w = 1, arm = a)
    pi <- adher == 1 & d$Tk >= gd
    post <- data.frame(start = gd, stop = d$Tk[pi], event = d$ev[pi], w = d$wfin[pi], arm = a)
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
## -- fit weights ONCE on di, estimate both landmarks (h=30 and 90), grace = gd --
.fit_and_estimate48 <- function(di, gd) {
  res <- c(rd30 = NA_real_, lrr30 = NA_real_, rmstd30 = NA_real_,
           rd90 = NA_real_, lrr90 = NA_real_, rmstd90 = NA_real_)
  try({
    di$pid <- seq_len(nrow(di))
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
    init <- di[di$initiator == 1, ]; init$treat <- as.integer(init$exposure_group == "PPI")
    if (min(sum(init$treat == 1), sum(init$treat == 0)) < 5) stop("treatment cell < 5")
    ft <- glm(reformulate(kc, "treat"), data = init, family = binomial())
    init$pppi <- pmin(pmax(fitted(ft), 1e-6), 1 - 1e-6)
    pp <- mean(init$treat)
    init$iptw <- ifelse(init$treat == 1, pp / init$pppi, (1 - pp) / (1 - init$pppi))
    ## combined stabilized weight = IPTW x IPCW, winsorized 1/99
    init$wraw <- init$iptw * init$ipcw
    ql <- quantile(init$wraw, 0.01); qh <- quantile(init$wraw, 0.99)
    init$wfin <- pmin(pmax(init$wraw, ql), qh)
    di$wfin <- 1; di$wfin[match(init$pid, di$pid)] <- init$wfin
    res[1:3] <- .est_one48(di, 30, gd); res[4:6] <- .est_one48(di, 90, gd)
  }, silent = TRUE)
  res
}


# #############################################################################
#  PART 1 (R0+R1) -- data prep + MICE (input: analytic_cohort_48h)
#    Skipped if the 48h imputation RDS already exists.
# #############################################################################
if (!file.exists(path_imp_48)) {
  
  con <- dbConnect(RPostgres::Postgres(),
                   host     = Sys.getenv("PGHOST", "localhost"),
                   port     = as.integer(Sys.getenv("PGPORT", "5432")),
                   dbname   = Sys.getenv("PGDATABASE", "mimiciv"),
                   user     = Sys.getenv("PGUSER", "postgres"),
                   password = Sys.getenv("PGPASSWORD"))
  dat <- dbGetQuery(con, "SELECT * FROM analytic_cohort_48h;")
  dbDisconnect(con)
  cat("Dimensions:", nrow(dat), "rows x", ncol(dat), "cols\n")
  
  cont_vars <- c(
    "age_at_icu","weight_kg","creat_pre","creat_24h","creat_mdrd","baseline_creat","baseline_egfr",
    "lactate","albumin","bilirubin_total","inr","platelet","wbc",
    "sodium","potassium","bicarbonate","anion_gap","magnesium","hemoglobin",
    "sofa_total","sofa_resp","sofa_coag","sofa_liver","sofa_cardio","sofa_cns","sofa_renal",
    "fluid_balance_ml","uo_rate_ml_kg_h","elixhauser_score","charlson_comorbidity_index",
    "hours_to_exposure","icu_los_days","los_after_t0_days","surv_28","surv_90","make_tte_time",
    "days_to_death","days_to_rrt","days_to_prog","days_to_recovery")
  factor_multi <- c("gender","race","admission_type","admission_location","first_careunit",
                    "anchor_year_group","exposure_group","baseline_creat_source","first_acid_drug","first_aki_stage")
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
  
  dat$weight_kg[dat$weight_kg < 20 | dat$weight_kg > 400] <- NA
  dat$uo_rate_ml_kg_h[dat$uo_rate_ml_kg_h < 0 | dat$uo_rate_ml_kg_h > 30] <- NA
  
  dat$race_cat <- factor(case_when(
    grepl("WHITE", dat$race) ~ "White", grepl("BLACK", dat$race) ~ "Black",
    grepl("HISPANIC|LATINO", dat$race) ~ "Hispanic", grepl("ASIAN", dat$race) ~ "Asian",
    TRUE ~ "Other/Unknown"))
  fc <- as.character(dat$first_careunit); rare <- names(which(table(fc) < 200))
  dat$first_careunit_cat <- factor(ifelse(fc %in% rare, "Other", fc))
  
  vars_impute <- c("weight_kg","albumin","lactate","bilirubin_total","inr","wbc",
                   "magnesium","anion_gap","platelet","sodium","potassium","bicarbonate","hemoglobin",
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
    "exposure_group","make30","make90","death_28","death_90")   # make30 added for congeniality
  imp_dat <- dat[, c(vars_impute, vars_predict)]; stay_key <- dat$stay_id
  
  ini <- mice(imp_dat, maxit = 0, printFlag = FALSE)
  meth <- ini$method; pred <- ini$predictorMatrix
  meth[vars_predict] <- ""; meth[setdiff(vars_impute, "sofa_total")] <- "pmm"
  meth["sofa_total"] <- "~ I(sofa_resp + sofa_coag + sofa_liver + sofa_cardio + sofa_cns + sofa_renal)"
  pred["sofa_total", ] <- 0; pred[, "sofa_total"] <- 0
  
  set.seed(20240601)
  imp <- mice(imp_dat, m = 20, maxit = 20, method = meth, predictorMatrix = pred, printFlag = TRUE)
  saveRDS(imp, "imp_mice_m20_48h.rds")
  long <- complete(imp, "long", include = TRUE)
  long$stay_id <- rep(stay_key, imp$m + 1)
  extra <- dat[, c("stay_id", setdiff(names(dat), names(imp_dat)))]
  long <- left_join(long, extra, by = "stay_id")
  saveRDS(long, path_imp_48)
  cat("\nSaved 48h imputed long table ->", path_imp_48, "(", nrow(long), "rows)\n")
  
} else cat("Detected", path_imp_48, "; skipping PART 1 (MICE).\n")


# #############################################################################
#  PART 2 (R2 + R4 + R10, 48h) -- IPTW (ITT-analog) primary + new-user
#    svyglm analytic SE (Rubin pooled) -- identical convention to main R4/R10.
# #############################################################################
hh0 <- as.data.frame(readRDS(path_imp_48))

## R2: head-to-head IPTW weights (PPI vs H2RA, ATE, stabilized, 1/99 trimmed)
ww <- list(); smd_full <- NA
for (i in 1:max(hh0$.imp)) {
  di <- hh0 %>% filter(.imp == i, exposure_group %in% c("PPI","H2RA")) %>%
    mutate(treat = as.integer(exposure_group == "PPI"))
  W  <- weightit(reformulate(ps_covars, "treat"), data = di, method = "glm", estimand = "ATE", stabilize = TRUE)
  Wt <- trim(W, at = 0.99, lower = TRUE); di$ipw <- Wt$weights
  if (i == 1) smd_full <- max(abs(bal.tab(Wt, un = FALSE, binary = "std", continuous = "std")$Balance$Diff.Adj), na.rm = TRUE)
  ww[[i]] <- di
}
hh_iptw <- bind_rows(ww)
cat(sprintf("\n48h head-to-head IPTW: weighted max|SMD| (imp1) = %.3f (<0.1 desirable)\n", smd_full))

## R4: IPTW primary MAKE30/90 (full cohort; RD + RR + RMST)
analyze_make <- function(dat, make_var, tau) {
  m <- max(dat$.imp); R <- data.frame()
  for (i in 1:m) {
    di <- dat %>% filter(.imp == i) %>%
      mutate(treat = as.integer(exposure_group == "PPI"),
             y = as.integer(as.character(.data[[make_var]])),
             ev = as.integer(as.character(make_tte_event)),
             tt = pmin(ifelse(ev == 1, make_tte_time, tau), tau)) %>% filter(!is.na(y))
    dz <- svydesign(ids = ~1, weights = ~ipw, data = di)
    frd <- svyglm(y ~ treat, design = dz, family = gaussian())
    frr <- svyglm(y ~ treat, design = dz, family = quasipoisson(link = "log"))
    frm <- svyglm(tt ~ treat, design = dz, family = gaussian())
    w0 <- di$treat == 0; w1 <- di$treat == 1
    R <- rbind(R, data.frame(rd = coef(frd)[["treat"]], v_rd = vcov(frd)["treat","treat"],
                             lrr = coef(frr)[["treat"]], v_lrr = vcov(frr)["treat","treat"],
                             rm = coef(frm)[["treat"]], v_rm = vcov(frm)["treat","treat"],
                             rh = weighted.mean(di$y[w0], di$ipw[w0]), rp = weighted.mean(di$y[w1], di$ipw[w1])))
  }
  list(n_ppi = sum(dat$.imp == 1 & dat$exposure_group == "PPI"), n_h2 = sum(dat$.imp == 1 & dat$exposure_group == "H2RA"),
       risk_ppi = mean(R$rp), risk_h2 = mean(R$rh),
       RD = pool_rubin(R$rd, R$v_rd), RR = expci(pool_rubin(R$lrr, R$v_lrr)), RMST = pool_rubin(R$rm, R$v_rm))
}
ip30 <- analyze_make(hh_iptw, "make30", 30)
ip90 <- analyze_make(hh_iptw, "make90", 90)

## R10: new-user IPTW (re-estimate weights within the 48h new-user subset)
newuser_make <- function(dat, make_var, tau) {
  m <- max(dat$.imp); R <- data.frame(); nppi <- 0; nh2 <- 0; sm <- c()
  for (i in 1:m) {
    di <- dat %>% filter(.imp == i, exposure_group %in% c("PPI","H2RA"),
                         as.integer(as.character(newuser)) == 1) %>%
      mutate(treat = as.integer(exposure_group == "PPI"))
    kc <- ps_covars[sapply(ps_covars, function(v) length(unique(di[[v]])) > 1)]
    W  <- weightit(reformulate(kc, "treat"), data = di, method = "glm", estimand = "ATE", stabilize = TRUE)
    Wt <- trim(W, at = 0.99, lower = TRUE); di$w <- Wt$weights
    if (i == 1) { nppi <- sum(di$treat == 1); nh2 <- sum(di$treat == 0) }
    sm <- c(sm, max(abs(bal.tab(Wt, un = FALSE, binary = "std", continuous = "std")$Balance$Diff.Adj), na.rm = TRUE))
    di$y <- as.integer(as.character(di[[make_var]])); d2 <- di[!is.na(di$y), ]
    dz <- svydesign(ids = ~1, weights = ~w, data = d2)
    frd <- svyglm(y ~ treat, design = dz, family = gaussian())
    frr <- svyglm(y ~ treat, design = dz, family = quasipoisson(link = "log"))
    w0 <- d2$treat == 0; w1 <- d2$treat == 1
    R <- rbind(R, data.frame(rd = coef(frd)[["treat"]], v_rd = vcov(frd)["treat","treat"],
                             lrr = coef(frr)[["treat"]], v_lrr = vcov(frr)["treat","treat"],
                             rh = weighted.mean(d2$y[w0], d2$w[w0]), rp = weighted.mean(d2$y[w1], d2$w[w1])))
  }
  list(n_ppi = nppi, n_h2 = nh2, max_smd = round(mean(sm), 3), risk_ppi = mean(R$rp), risk_h2 = mean(R$rh),
       RD = pool_rubin(R$rd, R$v_rd), RR = expci(pool_rubin(R$lrr, R$v_lrr)))
}
nu30 <- newuser_make(hh_iptw, "make30", 30)
nu90 <- newuser_make(hh_iptw, "make90", 90)

cat("\n========== 48h IPTW (ITT-analog) MAKE -- full cohort + new-user ==========\n")
fmt <- function(lab, r, extra = "") cat(sprintf("%-18s RR %.3f (%.3f-%.3f)  RD %+.2f%%  risk PPI %.1f%% / H2RA %.1f%%%s\n",
                                                lab, r$RR["est"], r$RR["lo"], r$RR["hi"], 100*r$RD["est"], 100*r$risk_ppi, 100*r$risk_h2, extra))
fmt("Full MAKE30",     ip30, sprintf("  [n PPI %d / H2RA %d]", ip30$n_ppi, ip30$n_h2))
fmt("Full MAKE90",     ip90)
fmt("New-user MAKE30", nu30, sprintf("  [n PPI %d / H2RA %d, max|SMD| %.3f]", nu30$n_ppi, nu30$n_h2, nu30$max_smd))
fmt("New-user MAKE90", nu90)

iptw48 <- bind_rows(
  data.frame(analysis="IPTW full",     outcome="MAKE30", n_ppi=ip30$n_ppi, n_h2=ip30$n_h2,
             RR=round(ip30$RR["est"],3), RR_lo=round(ip30$RR["lo"],3), RR_hi=round(ip30$RR["hi"],3),
             RD_pct=round(100*ip30$RD["est"],2), RD_lo=round(100*ip30$RD["lo"],2), RD_hi=round(100*ip30$RD["hi"],2)),
  data.frame(analysis="IPTW full",     outcome="MAKE90", n_ppi=ip90$n_ppi, n_h2=ip90$n_h2,
             RR=round(ip90$RR["est"],3), RR_lo=round(ip90$RR["lo"],3), RR_hi=round(ip90$RR["hi"],3),
             RD_pct=round(100*ip90$RD["est"],2), RD_lo=round(100*ip90$RD["lo"],2), RD_hi=round(100*ip90$RD["hi"],2)),
  data.frame(analysis="IPTW new-user", outcome="MAKE30", n_ppi=nu30$n_ppi, n_h2=nu30$n_h2,
             RR=round(nu30$RR["est"],3), RR_lo=round(nu30$RR["lo"],3), RR_hi=round(nu30$RR["hi"],3),
             RD_pct=round(100*nu30$RD["est"],2), RD_lo=round(100*nu30$RD["lo"],2), RD_hi=round(100*nu30$RD["hi"],2)),
  data.frame(analysis="IPTW new-user", outcome="MAKE90", n_ppi=nu90$n_ppi, n_h2=nu90$n_h2,
             RR=round(nu90$RR["est"],3), RR_lo=round(nu90$RR["lo"],3), RR_hi=round(nu90$RR["hi"],3),
             RD_pct=round(100*nu90$RD["est"],2), RD_lo=round(100*nu90$RD["lo"],2), RD_hi=round(100*nu90$RD["hi"],2)))
write.csv(iptw48, "iptw_48h_primary_newuser.csv", row.names = FALSE)
cat("\nWritten to iptw_48h_primary_newuser.csv\n")
cat("Reference (24h): IPTW full MAKE30 0.947 / MAKE90 0.967; IPTW new-user MAKE30 1.045.\n")


# #############################################################################
#  PART 3 (R11-R14, 48h) -- CCW point estimate + REFIT cluster-bootstrap CI
#    Grace GD=2. Weights re-estimated inside every resample (matches 24h R14).
# #############################################################################
hh0 <- as.data.frame(readRDS(path_imp_48))
if (!"neither" %in% hh0$exposure_group)
  stop("imp_long_analysis_48h.rds has no 'neither' -> re-run PART 1 on the FULL cohort.")
cat(sprintf("\nLoaded 48h imputations: %d rows, imp 1..%d\n", nrow(hh0), max(hh0$.imp)))

run_cohort_refit48 <- function(newuser_only, gd, B, n_cores = 1) {
  hh <- hh0[as.integer(as.character(hh0$both_in_window)) == 0, ]
  if (newuser_only) hh <- hh[as.integer(as.character(hh$newuser)) == 1, ]
  hh$initiator <- as.integer(hh$exposure_group %in% c("PPI","H2RA"))
  M <- max(hh$.imp); cols <- c("rd30","lrr30","rmstd30","rd90","lrr90","rmstd90")
  PT <- matrix(NA_real_, M, 6, dimnames = list(NULL, cols))
  VR <- matrix(NA_real_, M, 6, dimnames = list(NULL, cols)); NF <- numeric(M)
  tag <- if (newuser_only) "new-user" else "main"
  for (i in 1:M) {
    di0 <- hh[hh$.imp == i, ]; n <- nrow(di0)
    PT[i, ] <- .fit_and_estimate48(di0, gd)
    ob <- function(b) .fit_and_estimate48(di0[sample.int(n, n, replace = TRUE), , drop = FALSE], gd)
    bl <- if (n_cores > 1 && .Platform$OS.type == "unix")
      parallel::mclapply(1:B, ob, mc.cores = n_cores) else lapply(1:B, ob)
    bm <- do.call(rbind, bl); ok <- stats::complete.cases(bm); NF[i] <- sum(!ok)
    if (sum(ok) >= 2) VR[i, ] <- apply(bm[ok, , drop = FALSE], 2, var)
    cat(sprintf("  [%s 48h] imp %d/%d (boot fails %d/%d)\n", tag, i, M, NF[i], B))
  }
  row_h <- function(h) {
    s <- as.character(h)
    rd <- .poolt(PT[, paste0("rd", s)],    VR[, paste0("rd", s)])
    rr <- .poolt(PT[, paste0("lrr", s)],   VR[, paste0("lrr", s)])
    rm <- .poolt(PT[, paste0("rmstd", s)], VR[, paste0("rmstd", s)])
    data.frame(cohort = if (newuser_only) "New-user CCW (48h)" else "Main CCW (48h)", horizon = h,
               RD_pct = round(100*rd["est"],2), RD_lo = round(100*rd["lo"],2), RD_hi = round(100*rd["hi"],2),
               RR = round(exp(rr["est"]),3), RR_lo = round(exp(rr["lo"]),3), RR_hi = round(exp(rr["hi"]),3),
               RMSTd = round(rm["est"],3), RMSTd_lo = round(rm["lo"],3), RMSTd_hi = round(rm["hi"],3),
               RR_excl1 = (exp(rr["lo"]) > 1) | (exp(rr["hi"]) < 1), mean_boot_fail = round(mean(NF),1), row.names = NULL)
  }
  rbind(row_h(30), row_h(90))
}

cat("\nRunning main CCW (48h) refit bootstrap ... (slow)\n")
main <- run_cohort_refit48(FALSE, GD, B_BOOT, N_CORES)
cat("Running new-user CCW (48h) refit bootstrap ... (slow)\n")
nu   <- run_cohort_refit48(TRUE,  GD, B_BOOT, N_CORES)
res48 <- rbind(main, nu)
write.csv(res48, "ccw_48h_bootstrap_ci_refit.csv", row.names = FALSE)

cat("\n========== 48h sensitivity -- CCW (refit cluster bootstrap + Rubin, t-df) ==========\n")
print(res48, row.names = FALSE)

## Side by side with the corrected 24h result (ccw_bootstrap_ci_refit.csv from R14)
cat("\n---------- comparison with 24h main analysis (both refit) ----------\n")
c24 <- tryCatch(read.csv("ccw_bootstrap_ci_refit.csv"), error = function(e) NULL)
if (!is.null(c24)) {
  show_cmp <- function(coh24, coh48, lab) {
    a <- c24[grepl(coh24, c24$cohort), ]; b <- res48[res48$cohort == coh48, ]
    for (h in c(30, 90)) {
      x <- a[a$horizon == h, ]; y <- b[b$horizon == h, ]
      cat(sprintf("%-13s MAKE%d | 24h RR %.3f (%.3f-%.3f)   48h RR %.3f (%.3f-%.3f)\n",
                  lab, h, x$RR, x$RR_lo, x$RR_hi, y$RR, y$RR_lo, y$RR_hi))
    }
  }
  show_cmp("Main", "Main CCW (48h)", "Main CCW")
  show_cmp("New",  "New-user CCW (48h)", "New-user CCW")
} else {
  cat("(24h ccw_bootstrap_ci_refit.csv not found; run the corrected R14 first.)\n")
}

cat("\nInterpretation (window-robustness of the bias signature):\n")
cat("  - 48h Main CCW < 1 AND New-user CCW CI including 1 -> prevalent-user-bias\n")
cat("    conclusion holds in both windows (not a 24h artifact).\n")
cat("  - If New-user CCW (48h) excludes 1, soften wording to 'largely, not entirely,\n")
cat("    explained by prevalent-user bias'.\n")
cat("\nFiles written: iptw_48h_primary_newuser.csv, ccw_48h_bootstrap_ci_refit.csv\n")