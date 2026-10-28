# Model qualification on held-out data (never used in calibration).
# evERA experimental arm (giredestrant + everolimus), BOLERO-4, SERENA-6 (SERD switch proxy),
# postMONARCH, MAINTAIN, TRINITI-1, evERA 18-month OS and evERA experimental-arm toxicity.
# Trial emulations use the published assessment schedules, eligibility and dosing (R/arms.R).
# Writes results/validation.csv, results/validation_curves.rds, results/validation_summary.json

if (!exists("ROOT")) ROOT <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), ".."), winslash = "/")
source(file.path(ROOT, "R", "setup.R"))
NVAL <- as.integer(Sys.getenv("NVAL", "3000"))
p <- load_params(PFILES)
v <- make_vpop(p, NVAL, 31)
grid <- make_grid(p)
arms <- c(validation_arms(schedules = TRUE), pada1_arms())
res <- lapply(arms, function(a) simulate(p, v, grid, a, 2190L))
VT <- VAL_TARGETS
rows <- list()
add <- function(trial, endpoint, observed, obs_CI, model, criterion, passed)
  rows[[length(rows) + 1]] <<- data.frame(trial = trial, endpoint = endpoint, observed = as.character(observed),
                                          obs_CI = obs_CI, model = as.character(model), criterion = criterion,
                                          passed = passed, stringsAsFactors = FALSE)
within <- function(x, lo, hi) is.finite(x) && x >= lo && x <= hi
f2 <- function(x, d = 2) formatC(x, format = "f", digits = d)

ev <- evera(res$EV_exp, res$EV_ctrl)
for (lab in c("esr1", "all", "noesr1")) {
  o <- VT$EV[[lab]]
  m_ <- ev[[paste0("mPFS_exp_", lab)]]
  if (is.finite(o[2])) add("evERA", sprintf("median PFS, gir+eve (%s)", lab), o[1], sprintf("%g-%g", o[2], o[3]), f2(m_),
                           "within observed 95% CI", within(m_, o[2], o[3]))
  else add("evERA", sprintf("median PFS, gir+eve (%s)", lab), o[1], "", f2(m_), "descriptive", NA)
  add("evERA", sprintf("median PFS, ET+eve (%s)", lab), 5.5, "", f2(ev[[paste0("mPFS_ctl_", lab)]]),
      "calibrated (control arm)", NA)
  h <- ev[[paste0("HR_", lab)]][1]
  add("evERA", sprintf("HR PFS (%s)", lab), o[4], sprintf("%g-%g", o[5], o[6]), f2(h, 3), "within observed 95% CI",
      within(h, o[5], o[6]))
}
add("evERA", "HR weaker without ESR1 mutation", "0.84 vs 0.38", "",
    sprintf("%.2f vs %.2f", ev$HR_noesr1[1], ev$HR_esr1[1]), "HR(no ESR1) > HR(ESR1)", ev$HR_noesr1[1] > ev$HR_esr1[1])
for (lm in VT$EV_landmarks) {
  lab <- lm[1]; mo <- as.integer(lm[2])
  for (k in 1:2) {
    arm_ <- c("exp", "ctl")[k]; o <- as.numeric(lm[2 + k])
    mval <- ev[[sprintf("PFS%d_%s_%s", mo, arm_, lab)]]
    add("evERA", sprintf("%d-mo PFS %% %s (%s)", mo, arm_, lab), o, "", f2(mval, 1), "|diff| <= 10 points",
        abs(mval - o) <= 10)
  }
}
for (lab in names(VT$EV_ORR)) {
  o <- VT$EV_ORR[[lab]]
  add("evERA", sprintf("ORR gir+eve (%s)", lab), o[1], "", f2(ev[[paste0("ORR_exp_", lab)]], 3), "|diff| <= 0.10",
      abs(ev[[paste0("ORR_exp_", lab)]] - o[1]) <= 0.10)
  add("evERA", sprintf("ORR ET+eve (%s)", lab), o[2], "", f2(ev[[paste0("ORR_ctl_", lab)]], 3), "|diff| <= 0.10",
      abs(ev[[paste0("ORR_ctl_", lab)]] - o[2]) <= 0.10)
}
for (lab in names(VT$EV_OS18)) {
  o <- VT$EV_OS18[[lab]]
  me <- ev[[paste0("OS18_exp_", lab)]]; mc <- ev[[paste0("OS18_ctl_", lab)]]
  add("evERA", sprintf("18-mo OS %% gir+eve (%s)", lab), o[1], "", f2(me, 1), "|diff| <= 10 points", abs(me - o[1]) <= 10)
  if (lab == "all") add("evERA", sprintf("18-mo OS %% ET+eve (%s)", lab), o[2], "", f2(mc, 1), "calibrated (stage 4)", NA)
  else add("evERA", sprintf("18-mo OS %% ET+eve (%s)", lab), o[2], "", f2(mc, 1), "|diff| <= 10 points", abs(mc - o[2]) <= 10)
}
add("evERA", "ESR1 prevalence at entry (natural, before enrichment)", 0.40, "", f2(ev$esr1_frac_natural, 3),
    "descriptive", NA)
# BOLERO-4
b4 <- months(km_median(res$BOLERO4$pfs1, res$BOLERO4$pfs1_ev))
add("BOLERO-4", "median PFS, everolimus + letrozole (1L)", VT$BOLERO4[1], sprintf("%g-%g", VT$BOLERO4[2], VT$BOLERO4[3]),
    f2(b4), "within +/-20%", abs(b4 / VT$BOLERO4[1] - 1) <= 0.2)
# SERENA-6 (giredestrant as SERD proxy for camizestrant): extended analysis (Turner 2026),
# the trial's ctDNA schedule, effective detection limit, randomization delay and imaging check (R/arms.R)
s6 <- serena6(res$SER6_switch, res$SER6_cont)
S6 <- VT$SERENA6x
add("SERENA-6 (SERD proxy)", "median PFS after switch", S6$switch, "", f2(s6$mPFS_switch), "within +/-25%",
    abs(s6$mPFS_switch / S6$switch - 1) <= 0.25)
add("SERENA-6 (SERD proxy)", "median PFS, continue AI", S6$cont, "", f2(s6$mPFS_cont), "within +/-25%",
    abs(s6$mPFS_cont / S6$cont - 1) <= 0.25)
add("SERENA-6 (SERD proxy)", "HR", S6$HR[1], sprintf("%g-%g", S6$HR[2], S6$HR[3]),
    f2(s6$HR[1], 3), "within observed 95% CI", within(s6$HR[1], S6$HR[2], S6$HR[3]))
add("SERENA-6 (SERD proxy)", "24-month PFS % after switch", S6$PFS24[1], "", f2(s6$PFS24_switch, 1), "|diff| <= 10 points",
    abs(s6$PFS24_switch - S6$PFS24[1]) <= 10)
add("SERENA-6 (SERD proxy)", "24-month PFS %, continue AI", S6$PFS24[2], "", f2(s6$PFS24_cont, 1), "|diff| <= 10 points",
    abs(s6$PFS24_cont - S6$PFS24[2]) <= 10)
add("SERENA-6 (SERD proxy)", "OS HR (immature)", S6$OS_HR[1], sprintf("%g-%g", S6$OS_HR[2], S6$OS_HR[3]), f2(s6$OS_HR[1], 3),
    "within observed 95% CI", within(s6$OS_HR[1], S6$OS_HR[2], S6$OS_HR[3]))
det6 <- res$SER6_cont$esr1_detect >= 183 & res$SER6_cont$esr1_detect <= res$SER6_cont$pfs1
conc6 <- mean(res$SER6_cont$d3_switch[det6] < 0)
add("SERENA-6 (SERD proxy)", "ESR1 detected with concurrent progression (not randomized)", S6$conc_pd, "", f2(conc6, 3),
    "descriptive", NA)
add("SERENA-6 (SERD proxy)", "median time from start of line 1 to detection, randomized (mo)", S6$t_det, "", f2(s6$t_det, 1),
    "descriptive", NA)
add("SERENA-6 (SERD proxy)", "median allele fraction at detection, randomized (%)", 100 * S6$vaf_det, "", f2(100 * s6$vaf_det, 2),
    "descriptive (published: subsequent-test detections, read from figure)", NA)
# PADA-1 step 2: switch to fulvestrant + palbociclib vs continue AI + palbociclib on a rising ESR1 mutation (held out)
pd2 <- serena6(res$PD_switch, res$PD_cont, min_days = 0, excl_base = attr(res$PD_cont, "ct_lod"))
P1 <- VT$PADA1
add("PADA-1 step 2", "median PFS after switch to fulvestrant + palbociclib", P1$switch[1], sprintf("%g-%g", P1$switch[2], P1$switch[3]),
    f2(pd2$mPFS_switch), "within observed 95% CI", within(pd2$mPFS_switch, P1$switch[2], P1$switch[3]))
add("PADA-1 step 2", "median PFS, continue AI + palbociclib", P1$cont[1], sprintf("%g-%g", P1$cont[2], P1$cont[3]),
    f2(pd2$mPFS_cont), "within observed 95% CI", within(pd2$mPFS_cont, P1$cont[2], P1$cont[3]))
add("PADA-1 step 2", "HR", P1$HR[1], sprintf("%g-%g", P1$HR[2], P1$HR[3]), f2(pd2$HR[1], 3), "within observed 95% CI",
    within(pd2$HR[1], P1$HR[2], P1$HR[3]))
add("PADA-1 step 2", "update (Cabel 2026): medians switch / continue, HR", sprintf("%g / %g, %g", P1$upd[1], P1$upd[2], P1$upd[3]),
    sprintf("%g-%g", P1$upd[4], P1$upd[5]), sprintf("%.1f / %.1f, %.2f", pd2$mPFS_switch, pd2$mPFS_cont, pd2$HR[1]), "descriptive", NA)
add("PADA-1 step 2", "median allele fraction at the rising detection (%)", 100 * P1$vaf_det, "", f2(100 * pd2$vaf_det, 2),
    "descriptive (deterministic detection in the model)", NA)
# postMONARCH, MAINTAIN
pm <- two_arm_line2(res$pM_abema, res$pM_ful)
add("postMONARCH", "HR abemaciclib + fulvestrant vs fulvestrant", VT$postMONARCH$HR[1],
    sprintf("%g-%g", VT$postMONARCH$HR[2], VT$postMONARCH$HR[3]), f2(pm$HR[1], 3), "within observed 95% CI",
    within(pm$HR[1], VT$postMONARCH$HR[2], VT$postMONARCH$HR[3]))
add("postMONARCH", "median PFS abema + fulv / fulv", "6.0 / 5.3", "", sprintf("%.1f / %.1f", pm$mPFS_a, pm$mPFS_b),
    "descriptive", NA)
ma <- two_arm_line2(res$MA_ribo, res$MA_ful)
add("MAINTAIN", "HR ribociclib + fulvestrant vs fulvestrant", VT$MAINTAIN[1], sprintf("%g-%g", VT$MAINTAIN[2], VT$MAINTAIN[3]),
    f2(ma$HR[1], 3), "within observed 95% CI", within(ma$HR[1], VT$MAINTAIN[2], VT$MAINTAIN[3]))
# TRINITI-1
s1 <- second_line(res$TRINITI1)
tr <- months(km_median(s1$t, s1$e))
add("TRINITI-1", "median PFS ribociclib 300 (continuous) + everolimus 2.5 + exemestane", VT$TRINITI1[1],
    sprintf("%g-%g", VT$TRINITI1[2], VT$TRINITI1[3]), f2(tr), "within +/-30%", abs(tr / VT$TRINITI1[1] - 1) <= 0.3)
# ---- MIRACLE (first-line everolimus + letrozole), PARSIFAL and FALCON (first-line fulvestrant)
rfl <- lapply(c(miracle_arms(), fulvestrant_1l_arms()), function(a) simulate(p, v, grid, a, 2190L))
fl <- function(R) first_line(R)
mi_e <- fl(rfl$MI_eve); mi_l <- fl(rfl$MI_let)
hm <- cox_hr(mi_e$t, mi_e$e, mi_l$t, mi_l$e)
add("MIRACLE", "HR everolimus + letrozole vs letrozole (1L)", VT$MIRACLE$HR[1], sprintf("%g-%g", VT$MIRACLE$HR[2], VT$MIRACLE$HR[3]),
    f2(hm[1], 3), "within observed 95% CI", within(hm[1], VT$MIRACLE$HR[2], VT$MIRACLE$HR[3]))
add("MIRACLE", "median PFS everolimus + letrozole / letrozole", "19.4 / 12.9", "", sprintf("%.1f / %.1f", mi_e$mPFS, mi_l$mPFS),
    "descriptive (supplement: 16.5 / 10.0)", NA)
add("MIRACLE", "response everolimus + letrozole / letrozole", "0.500 / 0.393", "", sprintf("%.3f / %.3f", mi_e$ORR, mi_l$ORR),
    "descriptive", NA)
mde <- mean(rfl$MI_eve$ae_disc_eve | rfl$MI_eve$ae_disc_et)
add("MIRACLE", "stopped treatment for adverse events (everolimus arm)", VT$MIRACLE$ae_disc, "", f2(mde, 3), "|diff| <= 0.05",
    abs(mde - VT$MIRACLE$ae_disc) <= 0.05)
add("MIRACLE", "stomatitis any grade (everolimus arm, no mouthwash)", VT$MIRACLE$stom[["any"]], "", f2(mean(rfl$MI_eve$stom_any), 3),
    "|diff| <= 0.10", abs(mean(rfl$MI_eve$stom_any) - VT$MIRACLE$stom[["any"]]) <= 0.10)
add("MIRACLE", "stomatitis grade >= 3 (everolimus arm, no mouthwash)", VT$MIRACLE$stom[["g3"]], "", f2(mean(rfl$MI_eve$stom_g3), 3),
    "|diff| <= 0.05", abs(mean(rfl$MI_eve$stom_g3) - VT$MIRACLE$stom[["g3"]]) <= 0.05)
pf_f <- fl(rfl$PF_ful); pf_l <- fl(rfl$PF_let)
hp <- cox_hr(pf_f$t, pf_f$e, pf_l$t, pf_l$e)
add("PARSIFAL", "HR fulvestrant + palbociclib vs letrozole + palbociclib (1L)", VT$PARSIFAL$HR[1],
    sprintf("%g-%g", VT$PARSIFAL$HR[2], VT$PARSIFAL$HR[3]), f2(hp[1], 3), "within observed 95% CI",
    within(hp[1], VT$PARSIFAL$HR[2], VT$PARSIFAL$HR[3]))
add("PARSIFAL", "median PFS fulvestrant + palbo / letrozole + palbo", "27.9 / 32.8", "", sprintf("%.1f / %.1f", pf_f$mPFS, pf_l$mPFS),
    "descriptive", NA)
fa_f <- fl(rfl$FA_ful); fa_a <- fl(rfl$FA_ai)
hf <- cox_hr(fa_f$t, fa_f$e, fa_a$t, fa_a$e)
add("FALCON", "HR fulvestrant vs AI (1L, endocrine-naive)", VT$FALCON$HR[1], sprintf("%g-%g", VT$FALCON$HR[2], VT$FALCON$HR[3]),
    f2(hf[1], 3), "within observed 95% CI", within(hf[1], VT$FALCON$HR[2], VT$FALCON$HR[3]))
add("FALCON", "median PFS fulvestrant / AI", "16.6 / 13.8", "", sprintf("%.1f / %.1f", fa_f$mPFS, fa_a$mPFS), "descriptive", NA)

# evERA experimental-arm toxicity (exposure-driven extrapolation): toxicity-only runs with the
# simulated line-2 treatment durations
te <- second_line(res$EV_exp)$t
tc <- second_line(res$EV_ctrl)$t
vt <- make_vpop(p, NVAL, 41)
set.seed(5)
for (cfg in list(list(nm = "gir+eve", reg = regimen(NULL, "giredestrant", 10), dur = te, keys = VT$EV_tox_exp, calib = FALSE),
                 list(nm = "ET+eve", reg = regimen(NULL, "exemestane", 10), dur = tc,
                      keys = list(stom_any = 0.489, pneu_any = 0.124, ae_disc_eve = 0.118), calib = TRUE))) {
  d <- sample(pmax(cfg$dur, 1), NVAL, replace = TRUE)
  R <- simulate(p, vt, NULL, list(lines = list(cfg$reg), fixed_duration = as.integer(d)), 1900L)
  for (k in names(cfg$keys)) {
    mval <- mean(R[[k]]); o <- cfg$keys[[k]]
    add("evERA safety", sprintf("%s (%s)", k, cfg$nm), o, "", f2(mval, 3),
        if (cfg$calib) "calibrated" else "within +/-0.05 absolute", if (cfg$calib) NA else abs(mval - o) <= 0.05)
  }
}
df <- do.call(rbind, rows)
utils::write.csv(df, file.path(RES, "validation.csv"), row.names = FALSE)
print(df[, c("trial", "endpoint", "observed", "model", "passed")], right = FALSE)
saveRDS(c(ev$curves, list(s6_n = s6$n, b4_t = res$BOLERO4$pfs1, b4_e = res$BOLERO4$pfs1_ev, pd2_n = pd2$n, s6c = s6$curves,
                           pd2c = pd2$curves)),
        file.path(RES, "validation_curves.rds"))
np <- df$passed[!is.na(df$passed)]
cat(sprintf("PASSED %d of %d validation criteria\n", sum(np), length(np)))
save_json(list(n_pass = sum(np), n_crit = length(np)), "validation_summary.json")
progress_note("04_validate", sprintf("%d of %d validation criteria passed", sum(np), length(np)))
