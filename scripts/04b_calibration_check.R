# Re-simulate all calibration targets at the final point estimate with a larger, independent virtual
# population (out-of-sample check of the calibration), including PALOMA-2 extended follow-up, persevERA landmark PFS,
# PALOMA-3 and PADA-1 step 1 (ESR1 emergence).
# Also re-evaluates the joint stage 1 + 3 loss on the calibration population at the final parameters and checks the
# sampling noise of the PADA-1 step-1 targets on further 800-patient populations.
# Writes results/calibration_table.csv, calibration_curves.rds, calibration_extra.json and calib_stage13_final.json

if (!exists("ROOT")) ROOT <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), ".."), winslash = "/")
source(file.path(ROOT, "R", "setup.R"))
p <- load_params(PFILES)
v <- make_vpop(p, as.integer(Sys.getenv("N04B", "3000")), 51)
grid <- make_grid(p)
arms <- c(calibration_arms(schedules = TRUE), persevera_arms(), paloma3_arms(), pada1_arms()["PD_cont"])
res <- lapply(arms, function(a) simulate(p, v, grid, a, 2920L))
pg <- first_line(res$PV_gir); pl <- first_line(res$PV_let)
pv_hr <- cox_hr(pg$t, pg$e, pl$t, pl$e)
a <- first_line(res$P2_cdk); b <- first_line(res$P2_ai)
hr <- cox_hr(a$t, a$e, b$t, b$e)
b2 <- two_arm_line2(res$B2_eve, res$B2_exe)
s <- second_line(res$EV_ctrl)
es <- res$EV_ctrl$esr1_at_p1[s$m]
w <- evera_weights(es)
set.seed(99); pc <- stats::runif(v$n) < 0.42
ac <- acelera(res, pc)
mos <- months(km_median(res$ML2_os$os, res$ML2_os$os_ev))
ose <- (res$EV_ctrl$os - res$EV_ctrl$line2_start)[s$m]
os18 <- 100 * surv_at(ose, res$EV_ctrl$os_ev[s$m], 18 * MONTH, w)
p3 <- two_arm_line2(res$P3_pf, res$P3_f)
m3 <- res$P3_pf$in_line2 & res$P3_f$in_line2
pd <- pada1_step1(res$PD_cont, lag = arms$PD_cont$ct_lag)
T <- CAL_TARGETS
LM <- T$PV_landmarks
lmk <- function(r, m) 100 * surv_at(r$t, r$e, m * MONTH)
# type: "rel" (relative error, criterion |error| <= 15%) for medians and hazard ratios, or "pts" (absolute difference
# in percentage points, criterion |difference| <= 5 points) for every proportion; proportions are entered as fractions
row <- function(trial, target, published, model, type = "rel") data.frame(trial = trial, target = target, published = published,
                                                                         model = model, type = type, stringsAsFactors = FALSE)
tab <- rbind(
  row("PALOMA-2", "median PFS palbociclib + letrozole (mo; extended follow-up)", T$P2_cdk_mPFS, a$mPFS),
  row("PALOMA-2", "median PFS letrozole (mo)", T$P2_ai_mPFS, b$mPFS),
  row("PALOMA-2", "HR (extended follow-up)", T$P2_HR, hr[1]),
  row("PALOMA-2", "ORR palbociclib + letrozole", T$P2_ORR_cdk, a$ORR, "pts"),
  row("PALOMA-2", "ORR letrozole", T$P2_ORR_ai, b$ORR, "pts"),
  row("BOLERO-2", "median PFS exemestane (mo)", T$B2_exe_mPFS, b2$mPFS_b),
  row("BOLERO-2", "median PFS everolimus + exemestane (mo)", T$B2_eve_mPFS, b2$mPFS_a),
  row("BOLERO-2", "HR", T$B2_HR, b2$HR[1]),
  row("evERA control", "median PFS ET + everolimus, overall (mo)", T$EV_ctrl_mPFS_all, months(km_median(s$t, s$e, w))),
  row("evERA control", "median PFS ET + everolimus, ESR1-mutated (mo)", T$EV_ctrl_mPFS_esr1, months(km_median(s$t[es], s$e[es]))),
  row("acelERA", "median PFS giredestrant (mo)", T$AC_gir_mPFS, ac$mPFS_gir),
  row("acelERA", "median PFS physician choice ET (mo)", T$AC_ctl_mPFS, ac$mPFS_ctl),
  row("acelERA", "HR", T$AC_HR, ac$HR[1]),
  row("acelERA", "HR ESR1-mutated", T$AC_HR_esr1, ac$HR_esr1[1]),
  row("MONALEESA-2", "median OS ribociclib + letrozole (mo)", T$ML2_OS_cdk, mos),
  row("evERA control", "18-month OS ET + everolimus, ESR1-enriched", T$EV_ctrl_OS18 / 100, os18 / 100, "pts"),
  row("persevERA", "median PFS giredestrant + palbociclib (mo)", T$PV_gir_mPFS, pg$mPFS),
  row("persevERA", "median PFS letrozole + palbociclib (mo)", T$PV_let_mPFS, pl$mPFS),
  row("persevERA", "HR", T$PV_HR, pv_hr[1]),
  row("persevERA", "ORR difference, giredestrant - letrozole", T$PV_ORR_diff, pg$ORR - pl$ORR, "pts"),
  do.call(rbind, lapply(seq_along(LM$m), function(j) rbind(
    row("persevERA", sprintf("PFS at %d months, giredestrant + palbociclib", LM$m[j]), LM$gir[j] / 100, lmk(pg, LM$m[j]) / 100, "pts"),
    row("persevERA", sprintf("PFS at %d months, letrozole + palbociclib", LM$m[j]), LM$let[j] / 100, lmk(pl, LM$m[j]) / 100, "pts")))),
  row("PALOMA-3", "median PFS fulvestrant + palbociclib (mo)", T$P3_pf_mPFS, p3$mPFS_a),
  row("PALOMA-3", "median PFS fulvestrant (mo)", T$P3_f_mPFS, p3$mPFS_b),
  row("PALOMA-3", "HR", T$P3_HR, p3$HR[1]),
  row("PALOMA-3", "ORR fulvestrant + palbociclib (measurable disease)", T$P3_ORR[["pf"]], mean(res$P3_pf$resp2[m3]), "pts"),
  row("PALOMA-3", "ORR fulvestrant (measurable disease)", T$P3_ORR[["f"]], mean(res$P3_f$resp2[m3]), "pts"),
  row("PADA-1", "ESR1 mutation detectable at baseline", T$PD1$base_prev, pd$base_prev, "pts"),
  row("PADA-1", "rising ESR1 mutation, share of patients with progression or rising mutation", T$PD1$rise_share, pd$rise_share, "pts"),
  row("PADA-1", "synchronous progression, share of rising mutations", T$PD1$sync_share, pd$sync_share, "pts"),
  row("PADA-1", "median time on AI + palbociclib before randomization (mo)", T$PD1$t_rand, pd$t_rand),
  row("PADA-1", "median PFS, AI + palbociclib (mo)", T$PD1$mPFS, pd$mPFS))
tab$model <- round(tab$model, 4)
tab$rel_error_pct <- ifelse(tab$type == "rel", round(100 * (tab$model - tab$published) / tab$published, 1), NA)
tab$abs_error_points <- ifelse(tab$type == "pts", round(100 * (tab$model - tab$published), 1), NA)
utils::write.csv(tab[, c("trial", "target", "published", "model", "rel_error_pct", "abs_error_points")],
                 file.path(RES, "calibration_table.csv"), row.names = FALSE)
print(tab)
# sampling noise of the PADA-1 step-1 targets: the same metrics on the calibration population (800 patients, seed 21)
# and on five further populations of 800 patients
pd_seeds <- do.call(rbind, lapply(c(21, 52:56), function(sd) {
  vs <- make_vpop(p, 800L, sd)
  q <- pada1_step1(simulate(p, vs, grid, arms$PD_cont, 2920L), lag = arms$PD_cont$ct_lag)
  data.frame(seed = sd, base_prev = q$base_prev, rise_share = q$rise_share, sync_share = q$sync_share, t_rand = q$t_rand,
             mPFS = q$mPFS, n_rise = round(q$rise_all * vs$n))
}))
print(pd_seeds)
# descriptive: ESR1 detectable at progression on palbociclib + letrozole,
# median allele fraction at a rising detection (PADA-1: 0.8%)
cat(sprintf("ESR1 at progression (P2_cdk, default limit): %.3f; PADA-1 emulation: median VAF at rising detection %.4f, effective limit %.4f\n",
            a$esr1_prog, pd$vaf_det, attr(res$PD_cont, "ct_lod")))
save_json(list(esr1_prog_P2 = a$esr1_prog, pada1_vaf_det = pd$vaf_det, pada1_lod = attr(res$PD_cont, "ct_lod"),
               pada1_rise_all = pd$rise_all, pada1_net_inc = pd$net_inc, pada1_seeds = as.list(pd_seeds)), "calibration_extra.json")
t2b <- second_line(res$B2_eve); t2x <- second_line(res$B2_exe)
saveRDS(list(p2c_t = a$t, p2c_e = a$e, p2a_t = b$t, p2a_e = b$e, b2e_t = t2b$t, b2e_e = t2b$e, b2x_t = t2x$t,
             b2x_e = t2x$e, evc_t = s$t, evc_e = s$e, evc_w = w, ml2_t = res$ML2_os$os, ml2_e = res$ML2_os$os_ev,
             pvg_t = pg$t, pvg_e = pg$e, pvl_t = pl$t, pvl_e = pl$e,
             p3pf_t = second_line(res$P3_pf, m3)$t, p3pf_e = second_line(res$P3_pf, m3)$e,
             p3f_t = second_line(res$P3_f, m3)$t, p3f_e = second_line(res$P3_f, m3)$e,
             pd_net_t = pd$curves$t_net, pd_net_e = pd$curves$e_net, pd_pfs_t = res$PD_cont$pfs1, pd_pfs_e = res$PD_cont$pfs1_ev),
        file.path(RES, "calibration_curves.rds"))
ok <- ifelse(tab$type == "rel", abs(tab$rel_error_pct) <= 15, abs(tab$abs_error_points) <= 5)

# joint stage 1 + 3 loss of the calibration (03_calibrate_efficacy.R: 800 patients, seed 21) at the final parameters,
# i.e. after the stage-2 and stage-4 refits that follow the joint fit (stage 2's eps_g, the dependence of growth on S6K
# activity, also acts without everolimus)
CAL_AS_LIBRARY <- TRUE
source(file.path(ROOT, "scripts", "03_calibrate_efficacy.R"))
cal_loss <- function(q) {
  vq <- make_vpop(q, NCAL, SEED); gq <- make_grid(q)
  loss_stage("stage13", lapply(stats::setNames(STAGES$stage13$arms, STAGES$stage13$arms),
                               function(a) simulate(q, vq, gq, ARMS[[a]], STAGES$stage13$horizon)), vq)
}
Lfin <- cal_loss(p)
cat(sprintf("joint stage 1 + 3 loss on the calibration population: final parameters %.2f (end of the joint fit %.2f)\n",
            Lfin$loss, load_json("calib_stage13.json")$loss))
save_json(list(final_loss = Lfin$loss, final_info = Lfin$info), "calib_stage13_final.json")
progress_note("04b_calibration_check", sprintf("%d of %d targets met (relative error <= 15%%, or <= 5 points for proportions); joint loss at the final parameters %.2f",
                                               sum(ok), nrow(tab), Lfin$loss))
