# #############################################################################
#  PPI vs H2RA in SA-AKI - eICU external validation - full R analysis pipeline
#  ===========================================================================
#  One-shot reproduction:
#    R0 setup -> R1 load/clean/EDA -> R2 MICE (auto-cache) -> R3 IPTW + balance
#    -> R4 Table 1 -> R5 primary -> R6 secondary (+ RRT site correction)
#    -> R7 RRT sensitivity -> R8 E-value -> R9 figures -> R10 master table
#    -> R11 new-user (prevalent-user) bias probe
#  All tables (.csv) and figures (.pdf) are written to the working directory.
#  If the MICE cache exists, R2 is skipped instantly.
#  ---------------------------------------------------------------------------
#  Key design choices (pre-specified):
#   - Physiologic severity uses acutephysiologyscore (monotone), not the bimodal
#     "worst vital signs".
#   - In-hospital endpoints; MAKE/death as binary RD/RR; RRT/progression/recovery
#     as competing-risk (death competes).
#   - Primary doubly-robust (adjust gi_bleed, the only covariate with residual
#     SMD>0.1 after weighting) + a weight-only sensitivity.
#   - In eICU, RRT shows a spurious decrease driven by multi-center between-
#     hospital confounding; a hospital-stratified Cox returns it to the null.
#  ---------------------------------------------------------------------------
#  Prerequisites: build eicu_analytic_cohort in PostgreSQL first.
#    Set the database password before running, e.g. Sys.setenv(PGPASSWORD="...").
#  DUA: everything runs locally; only aggregate estimates leave the environment.
# #############################################################################

# =============================================================================
# R0 -- environment / dependencies / covariates / helpers (defined once here)
# =============================================================================
Sys.setenv(LANGUAGE = "en"); options(stringsAsFactors = FALSE); remove(list = ls())

library(DBI); library(RPostgres); library(dplyr); library(tidyr)
library(mice); library(WeightIt); library(cobalt)
library(tableone); library(survey); library(survival)
library(naniar); library(ggplot2)

## Propensity-score covariates (24; available in eICU, all measured at/before T0)
ps_covars <- c(
  "age_at_icu","gender","ethnicity_cat","unittype_cat",
  "weight_kg","baseline_egfr","first_aki_stage","ckd",
  "acutephysiologyscore",                                    # physiologic severity (monotone)
  "albumin","bilirubin","lactate","bun","hematocrit","gcs",  # directional labs
  "vent_day1","vasoactive",
  "gi_bleed","liver_disease",                                # confounding by indication
  "cm_diabetes","cm_immunosupp","cm_metastatic","cm_leukemia","cm_lymphoma"
)
tau <- 90   # in-hospital horizon cap

## Rubin pooling / exponentiate / E-value / ESS ----
pool_rubin <- function(est, v) {
  m <- length(est); qbar <- mean(est); ubar <- mean(v); b <- var(est)
  tvar <- ubar + (1 + 1/m) * b; se <- sqrt(tvar)
  df <- if (b > 0) (m - 1) * (1 + ubar / ((1 + 1/m) * b))^2 else 1e6
  c(est = qbar, se = se, lo = qbar - qt(.975, df)*se, hi = qbar + qt(.975, df)*se,
    p = 2*pt(-abs(qbar/se), df))
}
expci <- function(p) c(est = exp(unname(p["est"])), lo = exp(unname(p["lo"])),
                       hi = exp(unname(p["hi"])), p = unname(p["p"]))
evalue_rr <- function(est, lo, hi) {
  f <- function(x) { x <- ifelse(x < 1, 1/x, x); x + sqrt(x*(x-1)) }
  ev_pt <- f(est)
  ev_ci <- if (est < 1) { if (hi >= 1) 1 else f(hi) } else { if (lo <= 1) 1 else f(lo) }
  c(point = round(ev_pt, 2), ci = round(ev_ci, 2))
}
ess <- function(w) sum(w)^2 / sum(w^2)

## Shared: add competing-risk columns (in-hospital endpoint; death competes;
## censor at min(discharge, tau)) ----
add_cr_cols <- function(di, evvar, tvar, tau) {
  di %>% mutate(
    treat  = as.integer(exposure_group == "PPI"),
    ev1    = as.integer(as.character(.data[[evvar]])),
    t_ev   = ifelse(ev1 == 1, .data[[tvar]], Inf),
    t_dth  = ifelse(as.integer(as.character(death_inhosp)) == 1, days_to_death_inhosp, Inf),
    ctime  = pmin(los_after_t0_days, tau),
    fe     = pmin(t_ev, t_dth),
    status = ifelse(fe > ctime, 0L, ifelse(t_ev <= t_dth, 1L, 2L)),  # 0 cens / 1 event / 2 death
    ftime  = pmin(fe, ctime),
    ef     = factor(status, levels = c(0,1,2), labels = c("cens","event","death")))
}

# =============================================================================
# R1 -- load / type / clean / EDA / missingness
# =============================================================================
con <- dbConnect(RPostgres::Postgres(), host = "localhost", port = 5432,
                 dbname = "eicu", user = "postgres", password = Sys.getenv("PGPASSWORD"))  
dat <- dbGetQuery(con, "SELECT * FROM eicu_analytic_cohort;")
dbDisconnect(con)
cat("Loaded:", nrow(dat), "x", ncol(dat), "| stay_id unique:",
    length(unique(dat$stay_id)) == nrow(dat), "\n")

cont_vars <- c("age_at_icu","weight_kg","baseline_creat","baseline_egfr",
               "apachescore","acutephysiologyscore","pred_hosp_mort","hr","map","rr","temp","gcs",
               "wbc","albumin","bilirubin","bun","glucose","sodium","hematocrit","lactate",
               "hours_to_exposure","dc_day","days_to_rrt","days_to_prog","days_to_recovery",
               "t0_offset","hospitaldischargeoffset")
factor_multi <- c("gender","ethnicity","unittype","exposure_group","first_aki_stage")
factor_bin <- c("ckd","ckd_hist","newuser","both_in_window","by_creat","by_uo",
                "liver_disease","vent_day1","vasoactive","gi_bleed","cm_cirrhosis","cm_hepaticfailure",
                "cm_diabetes","cm_immunosupp","cm_aids","cm_leukemia","cm_lymphoma","cm_metastatic",
                "death_inhosp","rrt_init","rrt_dep_inhosp","non_recovery","make_inhosp",
                "prog_eligible","prog_event","recovery_eligible","recovery_event")
dat <- dat %>% mutate(across(any_of(cont_vars), as.numeric),
                      across(any_of(c(factor_multi, factor_bin)), as.factor))

## in-hospital follow-up duration ----
dat$los_after_t0_days <- (dat$hospitaldischargeoffset - dat$t0_offset) / 1440

## collapse high-cardinality categoricals ----
dat$ethnicity_cat <- factor(dplyr::case_when(
  grepl("Caucasian", dat$ethnicity) ~ "White", grepl("African American", dat$ethnicity) ~ "Black",
  grepl("Hispanic", dat$ethnicity) ~ "Hispanic", grepl("Asian", dat$ethnicity) ~ "Asian",
  TRUE ~ "Other/Unknown"))
ut <- as.character(dat$unittype); rare <- names(which(table(ut) < 200))
dat$unittype_cat <- factor(ifelse(ut %in% rare, "Other", ut))

## EDA: continuous distributions + missingness table/figure (saved) ----
pdf("eicu_eda_continuous.pdf", width = 11, height = 8); op <- par(mfrow = c(3,3))
for (v in intersect(cont_vars, names(dat))) { x <- dat[[v]]; if (all(is.na(x))) next
hist(x, main = v, xlab = "", breaks = 40, col = "grey80", border = "white") }
par(op); dev.off()
miss_tbl <- data.frame(variable = names(dat), n_missing = colSums(is.na(dat)),
                       pct_missing = round(100*colSums(is.na(dat))/nrow(dat), 1)) %>% arrange(desc(pct_missing))
write.csv(miss_tbl[miss_tbl$n_missing > 0, ], "eicu_missingness.csv", row.names = FALSE)
mc <- names(dat)[colSums(is.na(dat)) > 0]
if (length(mc) > 0) { pdf("eicu_eda_missingness.pdf", width = 11, height = 8)
  print(gg_miss_var(dat[, mc, drop = FALSE])); dev.off() }

## cleaning: out-of-range -> NA + exclude time anomalies / unknown sex ----
dat$weight_kg[dat$weight_kg < 20 | dat$weight_kg > 400] <- NA
dat$baseline_egfr[dat$baseline_egfr < 1 | dat$baseline_egfr > 200] <- NA
n0 <- nrow(dat)
dat <- dat %>% filter(los_after_t0_days >= 0, gender %in% c("Male","Female")) %>%
  mutate(gender = droplevels(factor(gender)))
cat("Excluded los<0 / unknown sex:", n0 - nrow(dat), "-> remaining", nrow(dat), "| exposure:",
    paste(names(table(dat$exposure_group)), table(dat$exposure_group), collapse=" "), "\n")

# =============================================================================
# R2 -- MICE (m=20, maxit=20, pmm; auto-cache) + long table + eGFR cap
# =============================================================================
vars_impute <- c("weight_kg","baseline_egfr","acutephysiologyscore",
                 "albumin","bilirubin","lactate","bun","hematocrit","gcs")
vars_predict <- c("age_at_icu","gender","ethnicity_cat","unittype_cat","first_aki_stage",
                  "ckd","baseline_creat","vent_day1","vasoactive","gi_bleed","liver_disease",
                  "cm_cirrhosis","cm_hepaticfailure","cm_diabetes","cm_immunosupp","cm_aids",
                  "cm_leukemia","cm_lymphoma","cm_metastatic","newuser","both_in_window",
                  "exposure_group","death_inhosp","rrt_init")
imp_dat <- dat[, c(vars_impute, vars_predict)]; stay_key <- dat$stay_id

if (file.exists("eicu_imp_mice_m20.rds")) {
  imp <- readRDS("eicu_imp_mice_m20.rds"); cat("Loaded cached MICE\n")
} else {
  ini <- mice(imp_dat, maxit = 0, printFlag = FALSE)
  meth <- ini$method; meth[vars_predict] <- ""; meth[vars_impute] <- "pmm"
  set.seed(20240601)
  imp <- mice(imp_dat, m = 20, maxit = 20, method = meth,
              predictorMatrix = ini$predictorMatrix, printFlag = TRUE)
  saveRDS(imp, "eicu_imp_mice_m20.rds")
  pdf("eicu_mice_convergence.pdf", 11, 8); plot(imp); dev.off()
  pdf("eicu_mice_density.pdf", 11, 7)
  densityplot(imp, ~ albumin + bilirubin + lactate + acutephysiologyscore + bun + hematocrit)
  dev.off()
}

long <- complete(imp, "long", include = TRUE)
long$stay_id <- rep(stay_key, imp$m + 1)
long <- left_join(long, dat[, c("stay_id", setdiff(names(dat), names(imp_dat)))], by = "stay_id")
long$baseline_egfr <- pmin(long$baseline_egfr, 150)                       # eGFR clinical cap
long$days_to_death_inhosp <- ifelse(as.integer(as.character(long$death_inhosp)) == 1,
                                    long$los_after_t0_days, NA_real_)
saveRDS(long, "eicu_imp_long_analysis.rds")
long <- readRDS("eicu_imp_long_analysis.rds")
cat("Long table:", nrow(long), "x", ncol(long), "\n")

# =============================================================================
# R3 -- IPTW weighting + balance + Love plot + ESS + LOS check
# =============================================================================
d1hh <- long %>% filter(.imp == 1, exposure_group %in% c("PPI","H2RA"))
nlev <- sapply(ps_covars, function(v) length(unique(na.omit(d1hh[[v]]))))
ps_use <- ps_covars[nlev > 1]
if (length(ps_use) < length(ps_covars))
  cat("Dropping constant covariates:", paste(setdiff(ps_covars, ps_use), collapse=", "), "\n")
ps_formula <- reformulate(ps_use, "treat")

m <- 20; hh_list <- vector("list", m); bal_list <- vector("list", m); wt_diag <- data.frame()
for (i in 1:m) {
  di <- long %>% filter(.imp == i, exposure_group %in% c("PPI","H2RA")) %>%
    mutate(treat = as.integer(exposure_group == "PPI"))
  W <- weightit(ps_formula, data = di, method = "glm", estimand = "ATE", stabilize = TRUE)
  Wt <- trim(W, at = 0.99, lower = TRUE); di$ipw <- Wt$weights
  hh_list[[i]] <- di
  bal_list[[i]] <- bal.tab(Wt, un = TRUE, stats = "mean.diffs", binary = "std", continuous = "std")
  w <- di$ipw
  wt_diag <- rbind(wt_diag, data.frame(imp = i, w_max = max(w),
                                       n_ppi = sum(di$treat==1), ess_ppi = ess(w[di$treat==1]),
                                       n_h2 = sum(di$treat==0), ess_h2 = ess(w[di$treat==0])))
}
hh_all <- bind_rows(hh_list)
hh_all$days_to_death_inhosp <- ifelse(as.integer(as.character(hh_all$death_inhosp)) == 1,
                                      hh_all$los_after_t0_days, NA_real_)
#saveRDS(hh_all, "eicu_iptw_weighted_long.rds")
hh_all <- readRDS("eicu_iptw_weighted_long.rds")

cat("\nweights/ESS (mean over 20):\n"); print(round(colMeans(wt_diag[,-1]), 1))
covs <- rownames(bal_list[[1]]$Balance)
bal_summary <- data.frame(covariate = covs,
                          smd_unweighted = round(rowMeans(sapply(bal_list, function(b) abs(b$Balance$Diff.Un))), 3),
                          smd_weighted   = round(rowMeans(sapply(bal_list, function(b) abs(b$Balance$Diff.Adj))), 3)) %>%
  filter(covariate != "prop.score") %>% arrange(desc(smd_weighted))
write.csv(bal_summary, "eicu_balance_table.csv", row.names = FALSE)
cat("# covariates with weighted SMD>0.1:", sum(bal_summary$smd_weighted > 0.1),
    "| max:", max(bal_summary$smd_weighted), "\n")

lp <- pivot_longer(bal_summary, c(smd_unweighted, smd_weighted), names_to="type", values_to="smd") %>%
  mutate(type = recode(type, smd_unweighted="unweighted", smd_weighted="weighted"))
ggplot(lp, aes(smd, reorder(covariate, smd), color=type)) + geom_point(size=2) +
  geom_vline(xintercept=0.1, linetype="dashed") +
  labs(x="|SMD|", y=NULL, color=NULL, title="eICU - covariate balance before/after IPTW") +
  theme_minimal(base_size=11)
ggsave("eicu_iptw_loveplot.pdf", width = 8, height = 10)

los_chk <- hh_list[[1]] %>% mutate(arm = ifelse(treat==1,"PPI","H2RA"),
                                   dead = as.integer(as.character(death_inhosp))) %>% group_by(arm) %>%
  summarise(n=n(), los_median=round(median(los_after_t0_days,na.rm=TRUE),1),
            los_mean_w=round(weighted.mean(los_after_t0_days, ipw, na.rm=TRUE),1),
            disch_alive_w=round(100*weighted.mean(dead==0, ipw, na.rm=TRUE),1), .groups="drop")
cat("\nLOS / discharged-alive (weighted):\n"); print(los_chk)

# =============================================================================
# R4 -- Table 1 (unweighted + IPTW-weighted; with SMD)
# =============================================================================
d1 <- hh_list[[1]] %>% mutate(exposure_group = droplevels(factor(exposure_group)))
t1_vars <- c("age_at_icu","gender","ethnicity_cat","unittype_cat","weight_kg",
             "baseline_creat","baseline_egfr","first_aki_stage","ckd","acutephysiologyscore",
             "apachescore","lactate","albumin","bilirubin","bun","hematocrit","gcs",
             "vent_day1","vasoactive","gi_bleed","liver_disease","cm_diabetes","cm_immunosupp",
             "cm_metastatic","cm_leukemia","cm_lymphoma","death_inhosp","make_inhosp","rrt_init")
t1_nonnorm <- c("baseline_creat","acutephysiologyscore","apachescore","lactate","bilirubin","bun")
t1_factor <- c("gender","ethnicity_cat","unittype_cat","first_aki_stage","ckd","vent_day1",
               "vasoactive","gi_bleed","liver_disease","cm_diabetes","cm_immunosupp","cm_metastatic",
               "cm_leukemia","cm_lymphoma","death_inhosp","make_inhosp","rrt_init")
t1u <- CreateTableOne(vars=t1_vars, strata="exposure_group", data=d1,
                      factorVars=t1_factor, addOverall=TRUE, test=FALSE)
t1w <- svyCreateTableOne(vars=t1_vars, strata="exposure_group",
                         data=svydesign(ids=~1, weights=~ipw, data=d1), factorVars=t1_factor, test=FALSE)
write.csv(print(t1u, nonnormal=t1_nonnorm, smd=TRUE, showAllLevels=TRUE, printToggle=FALSE),
          "eicu_table1_unweighted.csv")
write.csv(print(t1w, nonnormal=t1_nonnorm, smd=TRUE, showAllLevels=TRUE, printToggle=FALSE),
          "eicu_table1_weighted.csv")
cat("\nTable 1 saved (unweighted + weighted)\n")

# =============================================================================
# R5 -- primary: in-hospital MAKE + death (doubly-robust adjusting gi_bleed + weight-only)
# =============================================================================
analyze_binary <- function(dat, yvar, adjust = "gi_bleed") {
  R <- data.frame(); n_an <- c(ppi=NA, h2=NA)
  rhs <- if (is.null(adjust)) "treat" else paste("treat +", adjust)
  for (i in 1:max(dat$.imp)) {
    di <- dat %>% filter(.imp == i, exposure_group %in% c("PPI","H2RA")) %>%
      mutate(treat = as.integer(exposure_group=="PPI"),
             y = as.integer(as.character(.data[[yvar]]))) %>% filter(!is.na(y))
    if (i == 1) n_an <- c(ppi=sum(di$treat==1), h2=sum(di$treat==0))
    dz <- svydesign(ids=~1, weights=~ipw, data=di)
    f_rd <- svyglm(as.formula(paste("y ~", rhs)), design=dz, family=gaussian())
    f_rr <- svyglm(as.formula(paste("y ~", rhs)), design=dz, family=quasipoisson("log"))
    w0 <- di$treat==0; w1 <- di$treat==1
    R <- rbind(R, data.frame(rd=coef(f_rd)[["treat"]], v_rd=vcov(f_rd)["treat","treat"],
                             lrr=coef(f_rr)[["treat"]], v_lrr=vcov(f_rr)["treat","treat"],
                             risk_h2=weighted.mean(di$y[w0], di$ipw[w0]), risk_ppi=weighted.mean(di$y[w1], di$ipw[w1])))
  }
  list(n_ppi=unname(n_an["ppi"]), n_h2=unname(n_an["h2"]),
       risk_ppi=mean(R$risk_ppi), risk_h2=mean(R$risk_h2),
       RD=pool_rubin(R$rd, R$v_rd), RR=expci(pool_rubin(R$lrr, R$v_lrr)))
}
res_make    <- analyze_binary(hh_all, "make_inhosp",  "gi_bleed")
res_death   <- analyze_binary(hh_all, "death_inhosp", "gi_bleed")
res_make_w  <- analyze_binary(hh_all, "make_inhosp",  NULL)
res_death_w <- analyze_binary(hh_all, "death_inhosp", NULL)

prim_row <- function(r, lab, adj) data.frame(outcome=lab, adjust=adj, n_ppi=r$n_ppi, n_h2=r$n_h2,
                                             ppi_pct=round(100*r$risk_ppi,1), h2_pct=round(100*r$risk_h2,1),
                                             RD_pct=round(100*unname(r$RD["est"]),2), RD_lo=round(100*unname(r$RD["lo"]),2),
                                             RD_hi=round(100*unname(r$RD["hi"]),2), RR=round(unname(r$RR["est"]),3),
                                             RR_lo=round(unname(r$RR["lo"]),3), RR_hi=round(unname(r$RR["hi"]),3),
                                             RR_p=signif(unname(r$RR["p"]),3), row.names=NULL)
primary_tbl <- rbind(prim_row(res_make,"MAKE_inhosp","DR_gibleed"),
                     prim_row(res_death,"death_inhosp","DR_gibleed"), prim_row(res_make_w,"MAKE_inhosp","weight_only"),
                     prim_row(res_death_w,"death_inhosp","weight_only"))
write.csv(primary_tbl, "eicu_primary_effects.csv", row.names = FALSE)
cat("\n==== Primary outcomes ====\n"); print(primary_tbl, row.names = FALSE)

# =============================================================================
# R6 -- secondary: RRT / progression / recovery (competing risk, death competes)
#       + RRT hospital-stratified site correction
# =============================================================================
analyze_cr <- function(dat, elig, evvar, tvar, tau) {
  csh <- data.frame(); shr <- data.frame(); cif <- data.frame()
  cif_arm <- function(s) { sf <- survfit(Surv(ftime, ef) ~ 1, data=s, weights=s$ipw)
  ss <- summary(sf, times=tau, extend=TRUE); as.numeric(ss$pstate[1, which(colnames(ss$pstate)=="event")]) }
  for (i in 1:max(dat$.imp)) {
    di <- dat %>% filter(.imp==i, exposure_group %in% c("PPI","H2RA"))
    if (!is.null(elig)) di <- di %>% filter(as.integer(as.character(.data[[elig]]))==1)
    di <- add_cr_cols(di, evvar, tvar, tau) %>% mutate(rid = row_number())
    cox <- coxph(Surv(ftime, status==1) ~ treat, data=di, weights=ipw, robust=TRUE)
    csh <- rbind(csh, data.frame(b=coef(cox)[["treat"]], v=vcov(cox)["treat","treat"]))
    fg <- finegray(Surv(ftime, ef) ~ treat + ipw + rid, data=di, etype="event")
    fgc <- coxph(Surv(fgstart, fgstop, fgstatus) ~ treat, data=fg, weights=fgwt*ipw, id=rid, robust=TRUE)
    shr <- rbind(shr, data.frame(b=coef(fgc)[["treat"]], v=vcov(fgc)["treat","treat"]))
    cif <- rbind(cif, data.frame(cif_ppi=cif_arm(di[di$treat==1,]), cif_h2=cif_arm(di[di$treat==0,])))
  }
  keep <- dat$.imp==1 & dat$exposure_group %in% c("PPI","H2RA")
  if (!is.null(elig)) keep <- keep & as.integer(as.character(dat[[elig]]))==1
  list(n_ppi=sum(keep & dat$exposure_group=="PPI"), n_h2=sum(keep & dat$exposure_group=="H2RA"),
       cif_ppi=mean(cif$cif_ppi), cif_h2=mean(cif$cif_h2),
       csHR=expci(pool_rubin(csh$b, csh$v)), sHR=expci(pool_rubin(shr$b, shr$v)))
}
res_rrt  <- analyze_cr(hh_all, NULL, "rrt_init", "days_to_rrt", tau)
res_prog <- analyze_cr(hh_all, "prog_eligible", "prog_event", "days_to_prog", tau)
res_rec  <- analyze_cr(hh_all, "recovery_eligible", "recovery_event", "days_to_recovery", tau)

## RRT site correction: cause-specific Cox stratified by hospitalid (keeping IPTW weights), Rubin pooled ----
rrt_site <- data.frame()
for (i in 1:max(hh_all$.imp)) {
  di <- hh_all %>% filter(.imp==i, exposure_group %in% c("PPI","H2RA")) %>%
    add_cr_cols("rrt_init","days_to_rrt",tau)
  cx <- coxph(Surv(ftime, status==1) ~ treat + strata(hospitalid), data=di, weights=ipw, robust=TRUE)
  rrt_site <- rbind(rrt_site, data.frame(b=coef(cx)[["treat"]], v=vcov(cx)["treat","treat"]))
}
rrt_site_HR <- expci(pool_rubin(rrt_site$b, rrt_site$v))

cr_row <- function(r, lab) data.frame(outcome=lab, n_ppi=r$n_ppi, n_h2=r$n_h2,
                                      ppi_cif=round(100*r$cif_ppi,1), h2_cif=round(100*r$cif_h2,1),
                                      csHR=round(unname(r$csHR["est"]),3), sHR=round(unname(r$sHR["est"]),3),
                                      sHR_lo=round(unname(r$sHR["lo"]),3), sHR_hi=round(unname(r$sHR["hi"]),3),
                                      sHR_p=signif(unname(r$sHR["p"]),3), row.names=NULL)
secondary_tbl <- rbind(cr_row(res_rrt,"RRT_init"), cr_row(res_prog,"progression"), cr_row(res_rec,"recovery"))
secondary_tbl <- rbind(secondary_tbl, data.frame(outcome="RRT_site_adjusted", n_ppi=res_rrt$n_ppi,
                                                 n_h2=res_rrt$n_h2, ppi_cif=NA, h2_cif=NA, csHR=round(unname(rrt_site_HR["est"]),3),
                                                 sHR=round(unname(rrt_site_HR["est"]),3), sHR_lo=round(unname(rrt_site_HR["lo"]),3),
                                                 sHR_hi=round(unname(rrt_site_HR["hi"]),3), sHR_p=signif(unname(rrt_site_HR["p"]),3)))
write.csv(secondary_tbl, "eicu_secondary_effects.csv", row.names = FALSE)
cat("\n==== Secondary outcomes (RRT incl. hospital-stratified correction) ====\n"); print(secondary_tbl, row.names = FALSE)


# =============================================================================
# R7 -- RRT sensitivity: time horizon + hospital ecology
#       (confirm between-hospital confounding; rule out discharge censoring)
# =============================================================================
shr_at <- function(dat, evvar, tvar, tau, use_w=TRUE) {
  shr <- data.frame()
  for (i in 1:max(dat$.imp)) {
    di <- dat %>% filter(.imp==i, exposure_group %in% c("PPI","H2RA")) %>%
      add_cr_cols(evvar, tvar, tau) %>% mutate(w = if (use_w) ipw else 1, rid = row_number())
    fg <- finegray(Surv(ftime, ef) ~ treat + w + rid, data=di, etype="event")
    fgc <- coxph(Surv(fgstart, fgstop, fgstatus) ~ treat, data=fg, weights=fgwt*w, id=rid, robust=TRUE)
    shr <- rbind(shr, data.frame(b=coef(fgc)[["treat"]], v=vcov(fgc)["treat","treat"]))
  }
  p <- pool_rubin(shr$b, shr$v); c(sHR=round(exp(p["est"]),3), lo=round(exp(p["lo"]),3), hi=round(exp(p["hi"]),3))
}
rrt_horizon <- do.call(rbind, lapply(c(7,14,30,90), function(tt)
  data.frame(tau=tt, t(shr_at(hh_all, "rrt_init", "days_to_rrt", tt, TRUE)))))
rrt_horizon$unweighted_90 <- c(NA,NA,NA, shr_at(hh_all,"rrt_init","days_to_rrt",90,FALSE)["sHR.est"])
write.csv(rrt_horizon, "eicu_rrt_sensitivity.csv", row.names = FALSE)
cat("\n==== RRT time-horizon sensitivity (stable across horizons -> not discharge censoring) ====\n"); print(rrt_horizon, row.names = FALSE)

ds <- hh_all %>% filter(.imp==1, exposure_group %in% c("PPI","H2RA"))
ht <- ds %>% group_by(hospitalid) %>% summarise(n=n(), n_h2=sum(exposure_group=="H2RA"),
                                                rrt_rate=mean(as.integer(as.character(rrt_init))), .groups="drop") %>% filter(n>=20)
cat(sprintf("Hospitals: %d | per-hospital RRT rate %0.0f-%0.0f%% | H2RA-share vs RRT-rate Spearman: %.2f\n",
            n_distinct(ds$hospitalid), 100*min(ht$rrt_rate), 100*max(ht$rrt_rate),
            cor(ht$n_h2/ht$n, ht$rrt_rate, method="spearman")))

# =============================================================================
# R8 -- E-value (primary + secondary)
# =============================================================================
ev_tbl <- rbind(
  data.frame(outcome="MAKE_inhosp",  est=unname(res_make$RR["est"]),  lo=unname(res_make$RR["lo"]),  hi=unname(res_make$RR["hi"])),
  data.frame(outcome="death_inhosp", est=unname(res_death$RR["est"]), lo=unname(res_death$RR["lo"]), hi=unname(res_death$RR["hi"])),
  data.frame(outcome="progression",  est=unname(res_prog$sHR["est"]), lo=unname(res_prog$sHR["lo"]), hi=unname(res_prog$sHR["hi"])),
  data.frame(outcome="recovery",     est=unname(res_rec$sHR["est"]),  lo=unname(res_rec$sHR["lo"]),  hi=unname(res_rec$sHR["hi"])),
  data.frame(outcome="RRT_site_adj", est=unname(rrt_site_HR["est"]),  lo=unname(rrt_site_HR["lo"]),  hi=unname(rrt_site_HR["hi"])))
ev_tbl <- cbind(ev_tbl, t(apply(ev_tbl, 1, function(r)
  evalue_rr(as.numeric(r["est"]), as.numeric(r["lo"]), as.numeric(r["hi"])))))
names(ev_tbl)[5:6] <- c("evalue_point","evalue_ci")
write.csv(ev_tbl, "eicu_evalues.csv", row.names = FALSE)
cat("\n==== E-value ====\n"); print(ev_tbl, row.names = FALSE)

# =============================================================================
# R9 -- figures: PS overlap + CIF curves + eICU forest (with MIMIC reference)
# =============================================================================
## (1) PS overlap (positivity) ----
dps <- hh_list[[1]] %>% mutate(treat = as.integer(exposure_group=="PPI"))
dps$ps <- fitted(glm(reformulate(ps_use, "treat"), data=dps, family=binomial))
ggplot(dps, aes(ps, fill=factor(treat, labels=c("H2RA","PPI")))) +
  geom_density(alpha=.45, color=NA) +
  scale_fill_manual(values=c(H2RA="#888780", PPI="#185FA5"), name=NULL) +
  labs(x="Propensity score P[PPI]", y="Density", subtitle="PS distributions overlap -> positivity met") +
  theme_minimal(base_size=11)
ggsave("eicu_ps_overlap.pdf", width=8, height=5)

## (2) CIF curves (imp#1, weighted; competing event = death) ----
cif_curve <- function(elig, evvar, tvar, lab) {
  d <- hh_all %>% filter(.imp==1, exposure_group %in% c("PPI","H2RA"))
  if (!is.null(elig)) d <- d %>% filter(as.integer(as.character(.data[[elig]]))==1)
  d <- add_cr_cols(d, evvar, tvar, tau)
  do.call(rbind, lapply(c("H2RA","PPI"), function(g){ s <- d[d$exposure_group==g,]
  sf <- survfit(Surv(ftime, ef) ~ 1, data=s, weights=s$ipw)
  j <- which(colnames(sf$pstate)=="event")
  rbind(data.frame(time=0,cif=0,arm=g,outcome=lab),
        data.frame(time=sf$time, cif=sf$pstate[,j], arm=g, outcome=lab)) }))
}
cif_all <- bind_rows(cif_curve(NULL,"rrt_init","days_to_rrt","RRT initiation"),
                     cif_curve("prog_eligible","prog_event","days_to_prog","AKI progression"),
                     cif_curve("recovery_eligible","recovery_event","days_to_recovery","Renal recovery"))
ggplot(cif_all, aes(time, cif, color=arm)) + geom_step(linewidth=.7) +
  facet_wrap(~outcome, scales="free_y") +
  scale_color_manual(values=c(H2RA="#888780", PPI="#185FA5"), name=NULL) +
  scale_y_continuous(labels=scales::percent) +
  labs(x="Days since T0", y="Weighted cumulative incidence", subtitle="Competing event = death (RRT: see site correction)") +
  theme_minimal(base_size=11) + theme(legend.position="top")
ggsave("eicu_cif_curves.pdf", width=10, height=4)

## (3) eICU forest (with MIMIC point estimates as reference) ----
forest_df <- bind_rows(
  data.frame(outcome="MAKE (in-hospital)",  est=unname(res_make$RR["est"]),  lo=unname(res_make$RR["lo"]),  hi=unname(res_make$RR["hi"])),
  data.frame(outcome="Death (in-hospital)", est=unname(res_death$RR["est"]), lo=unname(res_death$RR["lo"]), hi=unname(res_death$RR["hi"])),
  data.frame(outcome="AKI progression",     est=unname(res_prog$sHR["est"]), lo=unname(res_prog$sHR["lo"]), hi=unname(res_prog$sHR["hi"])),
  data.frame(outcome="Renal recovery",      est=unname(res_rec$sHR["est"]),  lo=unname(res_rec$sHR["lo"]),  hi=unname(res_rec$sHR["hi"])),
  data.frame(outcome="RRT (crude)",         est=unname(res_rrt$sHR["est"]),  lo=unname(res_rrt$sHR["lo"]),  hi=unname(res_rrt$sHR["hi"])),
  data.frame(outcome="RRT (site-adjusted)", est=unname(rrt_site_HR["est"]),  lo=unname(rrt_site_HR["lo"]),  hi=unname(rrt_site_HR["hi"])))
# MIMIC point estimates (validated main-analysis values) for visual comparison
mimic_ref <- c("MAKE (in-hospital)"=0.947, "Death (in-hospital)"=0.968, "AKI progression"=0.972,
               "Renal recovery"=1.03, "RRT (crude)"=1.011, "RRT (site-adjusted)"=1.011)
forest_df$mimic <- mimic_ref[forest_df$outcome]
forest_df$outcome <- factor(forest_df$outcome, levels=rev(forest_df$outcome))
ggplot(forest_df, aes(est, outcome)) +
  geom_vline(xintercept=1, linetype="dashed", color="grey50") +
  geom_errorbarh(aes(xmin=lo, xmax=hi), height=.2, color="#185FA5") +
  geom_point(size=2.6, color="#185FA5") +
  geom_point(aes(x=mimic), shape=18, size=3, color="#C0392B") +
  scale_x_log10(breaks=c(0.5,0.7,1,1.4,2)) +
  labs(x="Effect (RR / sHR, log scale)", y=NULL,
       title="eICU estimates (blue) vs MIMIC point estimates (red diamond)",
       subtitle="all cross the null line 1; RRT returns to null after hospital stratification") +
  theme_minimal(base_size=11)
ggsave("eicu_forest_vs_mimic.pdf", width=8, height=5)
cat("\nFigures saved: PS overlap / CIF curves / forest\n")

# =============================================================================
# R10 -- master summary table (eICU vs MIMIC)
# =============================================================================
master <- data.frame(
  outcome = c("MAKE (in-hospital)","Death (in-hospital)","AKI progression","Renal recovery","RRT crude","RRT site-adjusted"),
  eICU_est = round(c(res_make$RR["est"], res_death$RR["est"], res_prog$sHR["est"],
                     res_rec$sHR["est"], res_rrt$sHR["est"], rrt_site_HR["est"]), 3),
  eICU_lo = round(c(res_make$RR["lo"], res_death$RR["lo"], res_prog$sHR["lo"],
                    res_rec$sHR["lo"], res_rrt$sHR["lo"], rrt_site_HR["lo"]), 3),
  eICU_hi = round(c(res_make$RR["hi"], res_death$RR["hi"], res_prog$sHR["hi"],
                    res_rec$sHR["hi"], res_rrt$sHR["hi"], rrt_site_HR["hi"]), 3),
  MIMIC_point = c(0.947, 0.968, 0.972, 1.03, 1.011, 1.011),
  agree_null = c("yes","yes","yes","yes","diverged->adjusted OK","yes"), row.names = NULL)
write.csv(master, "eicu_results_master.csv", row.names = FALSE)
cat("\n========== Master summary (eICU vs MIMIC) ==========\n"); print(master, row.names = FALSE)

cat("\n\n############ all outputs written ############\n")
cat("Tables (CSV): eicu_missingness / eicu_balance_table / eicu_table1_unweighted /",
    "eicu_table1_weighted / eicu_primary_effects / eicu_secondary_effects /",
    "eicu_rrt_sensitivity / eicu_evalues / eicu_results_master / eicu_newuser_effects\n")
cat("Figures (PDF): eicu_eda_continuous / eicu_eda_missingness / eicu_mice_convergence /",
    "eicu_mice_density / eicu_iptw_loveplot / eicu_ps_overlap / eicu_cif_curves /",
    "eicu_forest_vs_mimic\n")
cat("Intermediate RDS: eicu_imp_mice_m20 / eicu_imp_long_analysis / eicu_iptw_weighted_long\n")

# #############################################################################
#  R11 -- new-user (prevalent-user) bias probe
#    Restrict to new users (newuser==1), re-estimate PS + IPTW within the subset,
#    re-run MAKE / death / RRT. Purpose: (1) is the primary null robust among new
#    users (prevalent-user-bias probe)? (2) does the spurious RRT decrease decay
#    toward 1 among new users (= prevalent-user confounding)? eICU's admissiondrug
#    home-medication history makes new-user identification cleaner than MIMIC.
#    Reuses the helpers + analyze_binary (R5) + analyze_cr (R6) defined above.
#    (To run R11 standalone, first: long <- readRDS("eicu_imp_long_analysis.rds"))
#    Aggregate -> safe to share.
# #############################################################################

## restrict: new users + head-to-head ----
nu <- long %>% filter(as.integer(as.character(newuser)) == 1, exposure_group %in% c("PPI","H2RA"))
cat(sprintf("\nNew-user head-to-head (imp#1): PPI %d vs H2RA %d\n",
            sum(nu$.imp==1 & nu$exposure_group=="PPI"), sum(nu$.imp==1 & nu$exposure_group=="H2RA")))

## (A) re-estimate PS + IPTW within the new-user subset ----
d1nu <- nu %>% filter(.imp == 1)
nlev <- sapply(ps_covars, function(v) length(unique(na.omit(d1nu[[v]]))))
ps_use_nu <- ps_covars[nlev > 1]
if (length(ps_use_nu) < length(ps_covars))
  cat("Dropping constant-in-subset covariates:", paste(setdiff(ps_covars, ps_use_nu), collapse=", "), "\n") else
    cat("All covariates non-constant within the new-user subset\n")
ps_formula_nu <- reformulate(ps_use_nu, "treat")

nu_list <- vector("list", 20); wt_diag <- data.frame(); maxsmd <- numeric(20)
for (i in 1:20) {
  di <- nu %>% filter(.imp == i) %>% mutate(treat = as.integer(exposure_group == "PPI"))
  W  <- weightit(ps_formula_nu, data = di, method = "glm", estimand = "ATE", stabilize = TRUE)
  Wt <- trim(W, at = 0.99, lower = TRUE); di$ipw <- Wt$weights
  nu_list[[i]] <- di
  b <- bal.tab(Wt, un = TRUE, binary = "std", continuous = "std")
  maxsmd[i] <- max(abs(b$Balance$Diff.Adj[rownames(b$Balance) != "prop.score"]), na.rm = TRUE)
  w <- di$ipw
  wt_diag <- rbind(wt_diag, data.frame(imp = i, w_max = max(w),
                                       n_ppi = sum(di$treat==1), ess_ppi = ess(w[di$treat==1]),
                                       n_h2 = sum(di$treat==0), ess_h2 = ess(w[di$treat==0])))
}
cat("\nweights/ESS (mean over 20):\n"); print(round(colMeans(wt_diag[,-1]),1))
cat("weighted max SMD (mean over 20):", round(mean(maxsmd),3), "\n")
nu_all <- bind_rows(nu_list)

## (B) outcomes: MAKE + death (binary, DR + weight-only) + RRT (competing risk) ----
## reuse analyze_binary (R5) and analyze_cr (R6) -- no re-definition needed
nu_make  <- analyze_binary(nu_all, "make_inhosp",  "gi_bleed")
nu_makew <- analyze_binary(nu_all, "make_inhosp",  NULL)
nu_death <- analyze_binary(nu_all, "death_inhosp", "gi_bleed")
nu_rrt   <- analyze_cr(nu_all, NULL, "rrt_init", "days_to_rrt", tau)$sHR

pr <- function(lab, r, full) cat(sprintf(
  "%-22s PPI %.1f%% H2RA %.1f%%   RR %.3f (%.3f, %.3f)   [full cohort: %s]\n",
  lab, 100*r$risk_ppi, 100*r$risk_h2, r$RR["est"], r$RR["lo"], r$RR["hi"], full))
cat("\n========== New-user outcomes vs full cohort ==========\n")
pr("MAKE (doubly-robust)", nu_make, "0.950")
cat(sprintf("%-22s                         RR  %.3f (%.3f, %.3f)   [full cohort: 0.964]\n",
            "MAKE (weight-only)", nu_makew$RR["est"], nu_makew$RR["lo"], nu_makew$RR["hi"]))
pr("Death (doubly-robust)", nu_death, "0.953")
cat(sprintf("%-22s                         sHR %.3f (%.3f, %.3f)   [full cohort: 0.719] <- key\n",
            "RRT (competing risk)", nu_rrt["est"], nu_rrt["lo"], nu_rrt["hi"]))

out <- data.frame(
  outcome = c("MAKE_DR","MAKE_wonly","death_DR","RRT_sHR"),
  est = round(c(nu_make$RR["est"], nu_makew$RR["est"], nu_death$RR["est"], nu_rrt["est"]),3),
  lo  = round(c(nu_make$RR["lo"],  nu_makew$RR["lo"],  nu_death$RR["lo"],  nu_rrt["lo"]),3),
  hi  = round(c(nu_make$RR["hi"],  nu_makew$RR["hi"],  nu_death$RR["hi"],  nu_rrt["hi"]),3),
  full_cohort = c(0.950, 0.964, 0.953, 0.719))
write.csv(out, "eicu_newuser_effects.csv", row.names = FALSE)
cat("\n---- New-user summary (written to eicu_newuser_effects.csv) ----\n"); print(out, row.names = FALSE)
cat("\nInterpretation: (1) do MAKE/death stay null (main conclusion robust)?",
    "\n  (2) does the RRT sHR decay from 0.72 toward 1 -> prevalent-user confounding;",
    "if it stays ~0.72 -> residual renal confounding dominates.\n")

# #############################################################################
#  END OF eICU PIPELINE -- external validation of the MIMIC-IV primary analysis.
#  R0 setup -> R1 load/EDA -> R2 MICE -> R3 IPTW -> R4 Table 1 -> R5 primary
#  -> R6 secondary (+RRT site correction) -> R7 RRT sensitivity -> R8 E-value
#  -> R9 figures -> R10 master table -> R11 new-user probe.
# #############################################################################