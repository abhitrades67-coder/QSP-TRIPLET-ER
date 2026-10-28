# Tables and headline numbers computed from saved results.
# Writes results/report_numbers.json and results/tables/*.md

if (!exists("ROOT")) ROOT <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), ".."), winslash = "/")
source(file.path(ROOT, "R", "setup.R"))

fmt_num <- function(v, f = "%.2f") ifelse(is.na(v), "NA", sprintf(f, v))
md <- function(df, fmt = list(), default = "%.2f") {
  cols <- names(df)
  out <- c(paste0("| ", paste(cols, collapse = " | "), " |"), paste0("|", paste(rep("---", length(cols)), collapse = "|"), "|"))
  for (i in seq_len(nrow(df))) {
    cells <- sapply(cols, function(c) {
      v <- df[[c]][i]
      if (is.numeric(v)) { if (!is.finite(v)) "NA" else sprintf(if (!is.null(fmt[[c]])) fmt[[c]] else default, v) } else as.character(v)
    })
    out <- c(out, paste0("| ", paste(gsub("|", "\\|", cells, fixed = TRUE), collapse = " | "), " |"))
  }
  paste(out, collapse = "\n")
}
wr <- function(name, text) writeLines(text, file.path(TAB, paste0(name, ".md")), useBytes = TRUE)
q3 <- function(x) { x <- x[is.finite(x)]; if (!length(x)) return(c(NA, NA, NA)); c(stats::median(x), stats::quantile(x, .05, names = FALSE), stats::quantile(x, .95, names = FALSE)) }
fmt_ci <- function(x, d = 1) { m <- q3(x); if (!is.finite(m[1])) "NA" else sprintf(paste0("%.", d, "f (%.", d, "f to %.", d, "f)"), m[1], m[2], m[3]) }

N <- list()
p <- load_params(PFILES)
# ---------------- verification
ver <- utils::read.csv(file.path(RES, "verification.csv"), stringsAsFactors = FALSE)
ver$pass <- ifelse(is.na(ver$pass), "descriptive", ifelse(ver$pass, "yes", "NO"))
wr("verification", md(ver))
N$verification_pass <- sum(ver$pass == "yes"); N$verification_n <- sum(ver$pass != "descriptive")
# ---------------- toxicity calibration (event toxicities)
tc <- utils::read.csv(file.path(RES, "tox_calibration.csv"), stringsAsFactors = FALSE)
tc <- tc[!startsWith(tc$target, "info"), ]
wr("tox_calibration", md(data.frame(target = tc$target, regimen = tc$regimen, `published %` = 100 * tc$observed_target,
                                    `model %` = 100 * tc$model, parameter = tc$parameter, check.names = FALSE), default = "%.1f"))
N$tox_max_abs_err <- max(abs(tc$model - tc$observed_target)); N$tox_n <- nrow(tc)
# ---------------- AE-discontinuation sub-model
dc <- utils::read.csv(file.path(RES, "disc_calibration.csv"), stringsAsFactors = FALSE)
wr("disc_calibration", md(data.frame(`published rate` = dc$label, `published %` = 100 * dc$target, n = dc$n,
                                     `model %` = 100 * dc$model, `deviation (SE of the published rate)` = dc$z,
                                     source = dc$source, check.names = FALSE),
                          list(n = "%.0f", `deviation (SE of the published rate)` = "%+.1f"), default = "%.1f"))
N$disc_max_z <- max(abs(dc$z)); N$disc_n <- nrow(dc); N$disc_n_within2 <- sum(abs(dc$z) <= 2)
N$disc_worst <- dc$label[which.max(abs(dc$z))]; N$disc_worst_model <- dc$model[which.max(abs(dc$z))]
N$disc_worst_target <- dc$target[which.max(abs(dc$z))]
N$disc_tau <- p$tox$disc_tau; N$disc_floor <- p$tox$disc_floor
N$disc_pars <- p$tox[c("disc_cdk", "disc_eve", "disc_et", "disc_gir", "disc_et_eve", "disc_abema_mult", "disc_tau",
                       "disc_floor", "disc_t0", "disc_eve_gamma")]
if (file.exists(file.path(RES, "disc_shape_grid.csv"))) {
  sgd <- utils::read.csv(file.path(RES, "disc_shape_grid.csv"))
  N$disc_loss_const <- sgd$loss[sgd$floor >= 1][1]
  N$disc_loss_grid_best <- min(sgd$loss)
  N$disc_loss_final <- sum(dc$z^2)
}
# ---------------- efficacy calibration
cal <- utils::read.csv(file.path(RES, "calibration_table.csv"), stringsAsFactors = FALSE)
iod <- which(!is.na(cal$abs_error_points))
cal_ok <- ifelse(is.na(cal$rel_error_pct), abs(cal$abs_error_points) <= 5, abs(cal$rel_error_pct) <= 15)
wr("calibration", md(data.frame(trial = cal$trial, target = cal$target, published = cal$published, model = cal$model,
                                error = ifelse(is.na(cal$rel_error_pct), sprintf("%+.1f points", cal$abs_error_points),
                                                      sprintf("%+.1f%%", cal$rel_error_pct)),
                                check.names = FALSE), list(published = "%.3g", model = "%.3g")))
N$calib_within15 <- sum(cal_ok); N$calib_n <- nrow(cal)
N$calib_within25 <- sum(ifelse(is.na(cal$rel_error_pct), abs(cal$abs_error_points) <= 10, abs(cal$rel_error_pct) <= 25))
N$calib_miss <- paste0(cal$trial[!cal_ok], ": ", cal$target[!cal_ok])
pvr <- cal[cal$trial == "persevERA", ]
N$pv_cal <- stats::setNames(as.list(pvr$model), pvr$target)
N$cal_rows <- cal[, c("trial", "target", "published", "model", "rel_error_pct", "abs_error_points")]
N$cal_extra <- if (file.exists(file.path(RES, "calibration_extra.json"))) load_json("calibration_extra.json") else list()
# ---------------- validation (held out)
val <- utils::read.csv(file.path(RES, "validation.csv"), stringsAsFactors = FALSE)
val$passed <- ifelse(is.na(val$passed), "", ifelse(val$passed, "yes", "NO"))
wr("validation", md(data.frame(trial = val$trial, endpoint = val$endpoint, published = val$observed,
                               `published 95% CI` = ifelse(is.na(val$obs_CI), "", val$obs_CI), model = val$model,
                               criterion = val$criterion, passed = val$passed, check.names = FALSE)))
vv <- val[val$passed != "", ]
N$val_pass <- sum(vv$passed == "yes"); N$val_n <- nrow(vv)
N$val_fail <- paste(vv$trial[vv$passed == "NO"], vv$endpoint[vv$passed == "NO"], sep = ": ")
gng <- val[val$trial == "evERA" & grepl("^(median PFS, gir|HR PFS \\((esr1|all)\\)|HR weaker)", val$endpoint) & val$passed != "", ]
N$evera_gonogo_pass <- sum(gng$passed == "yes"); N$evera_gonogo_n <- nrow(gng)
# ---------------- external validation (SONIA; persevERA overall survival) and discontinuation in trial emulations
ext <- utils::read.csv(file.path(RES, "external_validation.csv"), stringsAsFactors = FALSE, check.names = FALSE)
sc <- function(m) c(sum(m == "yes", na.rm = TRUE), sum(m %in% c("yes", "NO")))
N$ext_met <- sc(ext$met); N$ext_sonia_met <- sc(ext$met[startsWith(ext$trial, "SONIA")])
dk <- utils::read.csv(file.path(RES, "disc_check.csv"), stringsAsFactors = FALSE, check.names = FALSE)
N$disc_check <- c(sum(dk$met == "yes"), nrow(dk))
# ---------------- ensemble (ABC)
ens <- load_json("ensemble.json")
N$abc_candidates <- ens$n_candidates; N$abc_accepted <- length(ens$members); N$abc_best_loss <- ens$best_loss; N$abc_threshold <- ens$threshold
cand <- utils::read.csv(file.path(RES, "abc_candidates.csv"), check.names = FALSE, stringsAsFactors = FALSE)
acc <- cand[cand$id %in% sapply(ens$members, `[[`, "id"), ]
pcols <- grep(".", names(cand), fixed = TRUE, value = TRUE)
wr("abc_params", md(data.frame(parameter = pcols, best_fit = sapply(pcols, function(c) get_param(p, c)),
                               ens_median = sapply(pcols, function(c) stats::median(acc[[c]])),
                               ens_5 = sapply(pcols, function(c) stats::quantile(acc[[c]], .05)),
                               ens_95 = sapply(pcols, function(c) stats::quantile(acc[[c]], .95))), default = "%.4g"))
# ---------------- scenario 1
es <- utils::read.csv(file.path(RES, "ensemble_summary.csv"), check.names = FALSE, stringsAsFactors = FALSE)
ref <- utils::read.csv(file.path(RES, "ref_summary.csv"), check.names = FALSE, stringsAsFactors = FALSE)
rownames(ref) <- ref$arm
arms6 <- c("A", "B", "C", "D1", "D2", "D3")
mets <- list(c("EDCT_RMST48", "EDCT RMST 48 mo", 1), c("EDCT_median", "EDCT median, mo", 1), c("OS_RMST60", "OS RMST 60 mo", 1),
             c("OS_median", "OS median, mo", 1), c("OS_60", "OS at 5 years", 2), c("PFS1_median", "PFS line 1 median, mo", 1),
             c("ORR1", "ORR line 1", 2), c("frac_line2", "started line 2", 2), c("PFS_line2_median", "PFS line 2 median, mo", 1),
             c("ae_disc_any", "AE-related discontinuation (any drug)", 2), c("sympt_g3", "symptomatic grade >=3 AE", 2),
             c("neut_g3", "grade >=3 neutropenia", 2), c("stom_g3", "grade 3 stomatitis", 3), c("pneu_g3", "grade >=3 pneumonitis", 3),
             c("QTWiST", "Q-TWiST 60 mo, mo", 1), c("ESR1_detect_24m", "ctDNA ESR1 detected by 24 mo", 2),
             c("frac_E_at_EDCT", "ESR1-mutant cell fraction at strategy failure", 2),
             c("frac_R_at_EDCT", "CDK4/6i-escape cell fraction at strategy failure", 2))
tab1 <- do.call(rbind, lapply(mets, function(m) { r <- data.frame(endpoint = m[2])
  for (a in arms6) r[[a]] <- fmt_ci(es[es$arm == a, m[1]], as.integer(m[3])); r }))
wr("scenario1_main", md(tab1))
tab2 <- do.call(rbind, lapply(mets, function(m) { r <- data.frame(endpoint = m[2])
  for (a in arms6) r[[a]] <- fmt_num(ref[a, m[1]], paste0("%.", m[3], "f")); r }))
wr("scenario1_ref", md(tab2))
piv <- function(metric) { w <- stats::reshape(es[, c("member", "arm", metric)], idvar = "member", timevar = "arm", direction = "wide")
  names(w) <- sub(paste0("^", metric, "\\."), "", names(w)); w }
PV <- lapply(stats::setNames(c("EDCT_RMST48", "OS_RMST60", "QTWiST", "ae_disc_any", "sympt_g3", "PFS1_median", "EDCT_median"),
                             c("EDCT_RMST48", "OS_RMST60", "QTWiST", "ae_disc_any", "sympt_g3", "PFS1_median", "EDCT_median")), piv)
crow <- list()
CONTRASTS <- list(c("A", "B"), c("A", "C"), c("B", "C"), c("D1", "B"), c("D2", "B"), c("D3", "B"), c("A", "D1"), c("A", "D3"),
                  c("A", "B_trial"), c("A_let", "B"))
clab <- function(a, b) if (b == "B_trial") "A vs B, trial-like attrition (15%)" else if (a == "A_let") "A_let (triplet with letrozole instead of giredestrant) vs B" else paste(a, "vs", b)
for (cc in CONTRASTS) {
  a <- cc[1]; b <- cc[2]
  if (!(b %in% names(PV$EDCT_RMST48))) next
  dE <- PV$EDCT_RMST48[[a]] - PV$EDCT_RMST48[[b]]; dO <- PV$OS_RMST60[[a]] - PV$OS_RMST60[[b]]
  dQ <- PV$QTWiST[[a]] - PV$QTWiST[[b]]; dT <- PV$ae_disc_any[[a]] - PV$ae_disc_any[[b]]
  crow[[length(crow) + 1]] <- data.frame(contrast = clab(a, b), `dEDCT48 (mo)` = fmt_ci(dE), `P(dEDCT48>=3)` = sprintf("%.2f", mean(dE >= 3)),
                                         `P(dEDCT48>0)` = sprintf("%.2f", mean(dE > 0)), `dOS60 (mo)` = fmt_ci(dO),
                                         `P(dOS60>0)` = sprintf("%.2f", mean(dO > 0)), `dQ-TWiST (mo)` = fmt_ci(dQ),
                                         `d AE-disc (pp)` = fmt_ci(100 * dT), check.names = FALSE)
  k <- paste0(a, b)
  N[[paste0("dE_", k)]] <- q3(dE); N[[paste0("dO_", k)]] <- q3(dO); N[[paste0("dQ_", k)]] <- q3(dQ); N[[paste0("dT_", k)]] <- q3(100 * dT)
  N[[paste0("PdE3_", k)]] <- mean(dE >= 3); N[[paste0("PdE0_", k)]] <- mean(dE > 0)
  N[[paste0("PdO0_", k)]] <- mean(dO > 0); N[[paste0("PdQ0_", k)]] <- mean(dQ > 0)
}
wr("contrasts", md(do.call(rbind, crow)))
# AE-related discontinuation of any drug by arm (ensemble medians): context for the tolerability criterion
tol <- do.call(rbind, lapply(intersect(c("A", "A_let", "B", "C", "D1", "D2", "D3", "A_eve5", "A_eve2.5", "B_trial"), names(PV$ae_disc_any)),
                             function(a) data.frame(arm = a, `AE-related discontinuation of any drug, % (ensemble median, 90% interval)` =
                                                      fmt_ci(100 * PV$ae_disc_any[[a]]), `symptomatic grade >=3 AE, %` =
                                                      fmt_ci(100 * PV$sympt_g3[[a]]), check.names = FALSE)))
wr("tolerability", md(tol))
for (a in c("B", "C", "D2", "A_eve5")) if (a %in% names(PV$ae_disc_any)) N[[paste0("disc_", a, "_median")]] <- stats::median(PV$ae_disc_any[[a]])
# ---------------- decision rule (A vs B, pre-specified in the protocol)
RA <- readRDS(file.path(RES, "ref_arms.rds"))
A <- RA$A; B <- RA$B; C <- RA$C
tau <- 60 * MONTH
qt <- function(R, ut, ur) { E <- pmin(R$edct, tau); O <- pmin(R$os, tau); T <- pmin(R$tox_days, E); mean(ut * T + (E - T) + ur * (O - E)) / MONTH }
us <- seq(0, 1, 0.1)
decide <- function(bname) {
  dE <- PV$EDCT_RMST48$A - PV$EDCT_RMST48[[bname]]; dO <- PV$OS_RMST60$A - PV$OS_RMST60[[bname]]; disc <- PV$ae_disc_any$A
  Bx <- RA[[bname]]
  qpos <- unlist(lapply(us, function(ut) sapply(us, function(ur) qt(A, ut, ur) - qt(Bx, ut, ur) > 0)))
  crit <- data.frame(criterion = c("P(dEDCT48 A-B >= 3 months) >= 0.80", "P(dOS60 A-B > 0) >= 0.70",
                                   "AE-related discontinuation of any drug in A <= 25% (ensemble median)",
                                   "Q-TWiST A-B > 0 across utility thresholds (fraction of 121 utility pairs)"),
                     value = c(sprintf("%.2f", mean(dE >= 3)), sprintf("%.2f", mean(dO > 0)), sprintf("%.1f%%", 100 * stats::median(disc)),
                               sprintf("%.2f", mean(qpos))),
                     met = c(mean(dE >= 3) >= 0.80, mean(dO > 0) >= 0.70, stats::median(disc) <= 0.25, mean(qpos) >= 0.999))
  fav <- all(crit$met); notfav <- (mean(dE <= 0) >= 0.5) || (stats::median(disc) > 0.25)
  list(crit = crit, decision = if (fav) "Triplet favored" else if (notfav) "Triplet not favored" else "Uncertain",
       P_dE_le0 = mean(dE <= 0), disc = stats::median(disc), qpos = mean(qpos))
}
D0 <- decide("B")
N$decision <- D0$decision; N$decision_criteria <- D0$crit; N$qtwist_frac_positive <- D0$qpos
N$P_dE_le0_AB <- D0$P_dE_le0; N$disc_A_median <- D0$disc
cr <- D0$crit; cr$met <- ifelse(cr$met, "yes", "no"); wr("decision", md(cr))
if ("B_trial" %in% names(PV$EDCT_RMST48) && !is.null(RA$B_trial)) {   # secondary, not pre-specified
  D1 <- decide("B_trial")
  N$decision_trial <- D1$decision; N$decision_criteria_trial <- D1$crit
  cr <- D1$crit; cr$met <- ifelse(cr$met, "yes", "no"); wr("decision_trial", md(cr))
}
# ---------------- threshold sensitivity (post hoc, defined after the results were known): the thresholds of the
# rule were set by judgement in the study protocol; tabulated so that readers can apply their own
BARS <- c(1, 1.5, 2, 2.5, 3, 4)
LIMITS <- c(0.25, 0.30, 0.40, 0.50)
bar_max <- function(dE) { s <- sort(dE); s[length(s) - ceiling(0.80 * length(s)) + 1] }   # largest bar with P(gain >= bar) >= 0.80
ATT <- list(rw = list(lab = "real-world (pre-specified analysis)", b = "B", D = D0))
if (exists("D1")) ATT$tr <- list(lab = "trial-like, 15% (secondary analysis)", b = "B_trial", D = D1)
decide_at <- function(dE, crit, bar, limit, disc) {
  fav <- mean(dE >= bar) >= 0.80 && crit$met[2] && disc <= limit && crit$met[4]
  if (fav) "favored" else if (mean(dE <= 0) >= 0.5 || disc > limit) "not favored" else "uncertain"
}
te <- list(); tdec <- list()
for (k in names(ATT)) {
  dE <- PV$EDCT_RMST48$A - PV$EDCT_RMST48[[ATT[[k]]$b]]
  r <- data.frame(`attrition in Arm B` = ATT[[k]]$lab, `gain A-B, months (median, 90% interval)` = fmt_ci(dE), check.names = FALSE)
  for (b in BARS) r[[sprintf("P(gain >= %s mo)%s", format(b), if (b == 3) ", pre-specified bar" else "")]] <- sprintf("%.2f", mean(dE >= b))
  r[["largest bar met with P >= 0.80 (mo)"]] <- sprintf("%.1f", floor(10 * bar_max(dE) + 1e-9) / 10)   # rounded down
  te[[k]] <- r
  N[[paste0("thr_P_", k)]] <- stats::setNames(as.list(sapply(BARS, function(b) mean(dE >= b))), as.character(BARS))
  N[[paste0("thr_barmax_", k)]] <- bar_max(dE)
  for (lim in LIMITS) {
    r <- data.frame(`attrition in Arm B` = ATT[[k]]$lab,
                    `tolerability limit, any drug` = sprintf("%.0f%%%s", 100 * lim, if (lim == 0.25) " (pre-specified)" else ""), check.names = FALSE)
    for (b in BARS) r[[sprintf("bar %s mo%s", format(b), if (b == 3) " (pre-specified)" else "")]] <- decide_at(dE, ATT[[k]]$D$crit, b, lim, ATT[[k]]$D$disc)
    tdec[[length(tdec) + 1]] <- r
  }
}
wr("threshold_efficacy", md(do.call(rbind, te)))
wr("threshold_decision", md(do.call(rbind, tdec)))
N$thr_disc_A <- D0$disc
N$crit_os_qt_rw <- D0$crit$met[2] && D0$crit$met[4]
if (!is.null(ATT$tr)) N$crit_os_qt_tr <- ATT$tr$D$crit$met[2] && ATT$tr$D$crit$met[4]
N$thr_disc_variants <- lapply(stats::setNames(nm = intersect(c("A_let", "A_eve5", "A_eve2.5", "D2"), names(PV$ae_disc_any))),
                              function(a) stats::median(PV$ae_disc_any[[a]]))
# ---------------- secondary tolerability criterion (post hoc, defined after the results were known): the backbone
# (CDK4/6 inhibitor or endocrine agent) stopped for adverse events; stopping everolimus alone counts as a step-down to the
# standard doublet. Same 25% limit and ensemble median as the pre-specified criterion (06c_disc_components.R)
if (file.exists(file.path(RES, "ensemble_disc_components.csv"))) {
  dcc <- utils::read.csv(file.path(RES, "ensemble_disc_components.csv"), stringsAsFactors = FALSE)
  carms <- intersect(c("A", "A_let", "A_eve5", "A_eve2.5", "A_tox1.5", "B", "B_trial", "C", "D1", "D2", "D3"), unique(dcc$arm))
  wr("disc_components", md(do.call(rbind, lapply(carms, function(a) { g <- dcc[dcc$arm == a, ]
    data.frame(arm = a, `any drug (pre-specified criterion)` = fmt_ci(100 * g$any), `CDK4/6 inhibitor` = fmt_ci(100 * g$cdk),
               everolimus = fmt_ci(100 * g$eve), `endocrine agent` = fmt_ci(100 * g$et),
               `backbone: CDK4/6 inhibitor or endocrine agent (secondary criterion)` = fmt_ci(100 * g$backbone),
               `everolimus only` = fmt_ci(100 * g$eve_only), check.names = FALSE) }))))
  bb <- function(a, col = "backbone") { g <- dcc[dcc$arm == a, ]; stats::setNames(g[[col]], g$member) }
  N$bb <- lapply(stats::setNames(nm = carms), function(a) q3(100 * bb(a)))
  N$eve_only <- lapply(stats::setNames(nm = carms), function(a) q3(100 * bb(a, "eve_only")))
  N$comp_cdk <- lapply(stats::setNames(nm = carms), function(a) q3(100 * bb(a, "cdk")))
  N$comp_eve <- lapply(stats::setNames(nm = carms), function(a) q3(100 * bb(a, "eve")))
  N$comp_et <- lapply(stats::setNames(nm = carms), function(a) q3(100 * bb(a, "et")))
  bA <- bb("A"); bB <- bb("B")
  N$dBB_AB <- q3(100 * (bA - bB[names(bA)]))
  sec <- list()
  for (k in names(ATT)) {
    D_ <- ATT[[k]]$D; dE <- PV$EDCT_RMST48$A - PV$EDCT_RMST48[[ATT[[k]]$b]]; disc_bb <- stats::median(bA)
    met <- c(D_$crit$met[1:2], disc_bb <= 0.25, D_$crit$met[4])
    dec <- if (all(met)) "Triplet favored" else if (mean(dE <= 0) >= 0.5 || disc_bb > 0.25) "Triplet not favored" else "Uncertain"
    N[[paste0("decision_bb_", k)]] <- dec
    sec[[k]] <- data.frame(criterion = c(D_$crit$criterion[1:2],
                                         "AE-related discontinuation of the backbone (CDK4/6 inhibitor or endocrine agent) in A <= 25% (ensemble median)",
                                         D_$crit$criterion[4], "Decision"),
                           value = c(D_$crit$value[1:2], sprintf("%.1f%%", 100 * disc_bb), D_$crit$value[4], dec),
                           met = c(ifelse(met, "yes", "no"), ""), check.names = FALSE)
  }
  st_ <- sec$rw; names(st_) <- c("criterion", "value, real-world attrition", "met")
  if (!is.null(sec$tr)) { st_[["value, trial-like attrition"]] <- sec$tr$value; st_[["met (trial-like)"]] <- sec$tr$met }
  wr("decision_secondary", md(st_))
  N$bb_A_median <- stats::median(bA)
  N$bb_A_range <- range(100 * bA); N$bb_A_share_le25 <- mean(bA <= 0.25)
}
# ---------------- variants (reference population)
va <- intersect(c("A", "A_let", "A_eve5", "A_eve2.5", "A_ribo_eve2.5", "A_ribo_eve10", "A_ribo300_eve2.5", "A_abema_eve10", "A_abema_eve5",
                  "A_tox1.5", "A_tox2.0", "A_eve5_tox2.0", "A_nodex", "A_dexstrict", "A_eve5_dexstrict", "B", "B_ribo", "B_abema",
                  "B_nodex", "B_dexstrict", "B_trial", "B_trial_10", "B_trial_20", "C", "C_ribo", "C_abema", "C_ful", "C_trial",
                  "B_ful1L", "C_ful1L"), ref$arm)
vt <- data.frame(arm = va, `EDCT RMST48 (mo)` = ref[va, "EDCT_RMST48"], `OS RMST60 (mo)` = ref[va, "OS_RMST60"],
                 `PFS1 median (mo)` = ref[va, "PFS1_median"], `AE disc any %` = 100 * ref[va, "ae_disc_any"],
                 `sympt. G>=3 %` = 100 * ref[va, "sympt_g3"], `neutrop. G>=3 %` = 100 * ref[va, "neut_g3"],
                 `stomat. G3 %` = 100 * ref[va, "stom_g3"], `pneum. G>=3 %` = 100 * ref[va, "pneu_g3"],
                 `diarr. G3 %` = 100 * ref[va, "diar_g3"], `RDI eve %` = 100 * ref[va, "RDI_eve"], `RDI CDK4/6i %` = 100 * ref[va, "RDI_cdk"],
                 `Q-TWiST (mo)` = ref[va, "QTWiST"], check.names = FALSE)
wr("variants", md(vt, default = "%.1f"))
# ---------------- scenario 2
s2 <- c("S2_A'(abema)", "S2_A'(palbo cont.)", "S2_B'", "S2_C'")
m2 <- RA[["S2_C'"]]$in_line2
l2 <- function(R) list(t = (R$pfs2 - R$line2_start)[m2], e = R$pfs2_ev[m2])
cR <- l2(RA[["S2_C'"]]); bR <- l2(RA[["S2_B'"]])
s2rows <- do.call(rbind, lapply(s2, function(k) { R <- RA[[k]]; x <- l2(R); osr <- (R$os - R$line2_start)[m2]
  data.frame(arm = sub("S2_", "", k), `median PFS (mo)` = months(km_median(x$t, x$e)), `HR vs C'` = cox_hr(x$t, x$e, cR$t, cR$e)[1],
             `HR vs B'` = cox_hr(x$t, x$e, bR$t, bR$e)[1], `RMST PFS 24 mo` = mean(pmin(x$t, 24 * MONTH)) / MONTH,
             `OS RMST 36 mo` = mean(pmin(osr, 36 * MONTH)) / MONTH, check.names = FALSE) }))
wr("scenario2", md(s2rows, list(`HR vs C'` = "%.2f", `HR vs B'` = "%.2f"), default = "%.1f"))
pv2 <- piv("PFS_line2_median")
N$S2_ApvsBp <- q3(pv2[["S2_A'(abema)"]] - pv2[["S2_B'"]]); N$S2_BpvsCp <- q3(pv2[["S2_B'"]] - pv2[["S2_C'"]])
N$S2_rows <- s2rows
# ---------------- subgroups
sg <- utils::read.csv(file.path(RES, "ensemble_subgroups.csv"), stringsAsFactors = FALSE)
sg_table <- function(contrast) {
  d <- sg[sg$contrast == contrast, ]
  do.call(rbind, lapply(unique(d$subgroup), function(nm) { g <- d[d$subgroup == nm, ]
    data.frame(subgroup = nm, `n (per member)` = as.integer(stats::median(g$n)), `dEDCT48 A-B (mo)` = fmt_ci(g$dEDCT48),
               `P(>=3 mo)` = sprintf("%.2f", mean(g$dEDCT48 >= 3)), `dOS60 A-B (mo)` = fmt_ci(g$dOS60), check.names = FALSE) }))
}
sgr <- sg_table("A-B"); wr("subgroups", md(sgr)); N$subgroups <- sgr
if (any(sg$contrast == "A-B_trial")) { sgt <- sg_table("A-B_trial"); wr("subgroups_trial", md(sgt)); N$subgroups_trial <- sgt }
# ---------------- GSA / EVPPI
if (file.exists(file.path(RES, "gsa_prcc.csv"))) {
  pr <- utils::read.csv(file.path(RES, "gsa_prcc.csv"), stringsAsFactors = FALSE)
  outs <- c("dEDCT48_AB", "dOS60_AB", "dQTW_AB", "AEdisc_A")
  pars <- unique(pr$parameter)
  W <- sapply(outs, function(o) sapply(pars, function(pp) { r <- pr[pr$output == o & pr$parameter == pp, ]
    sprintf("%+.2f%s", r$PRCC, if (isTRUE(r$p_value < 0.05)) "*" else "") }))
  ordv <- sapply(pars, function(pp) abs(pr$PRCC[pr$output == "dEDCT48_AB" & pr$parameter == pp]))
  tt <- data.frame(parameter = pars, W, check.names = FALSE)[order(-ordv), ]
  wr("prcc", md(tt))
  ev <- utils::read.csv(file.path(RES, "gsa_evppi.csv"), stringsAsFactors = FALSE)
  et <- stats::reshape(ev, idvar = "parameter", timevar = "output", direction = "wide")
  names(et) <- sub("^value\\.", "", names(et))
  et <- et[order(-et$dEDCT48_AB), ]
  wr("evppi", md(et, default = "%.3f"))
  gs <- utils::read.csv(file.path(RES, "gsa_samples.csv"), check.names = FALSE)
  N$gsa_n <- nrow(gs); N$gsa_dE_AB <- q3(gs$dEDCT48_AB); N$gsa_P3 <- mean(gs$dEDCT48_AB >= 3); N$gsa_P0 <- mean(gs$dEDCT48_AB > 0)
  N$gsa_dO_AB <- q3(gs$dOS60_AB); N$gsa_discA <- q3(100 * gs$AEdisc_A); N$gsa_discA_le25 <- mean(gs$AEdisc_A <= 0.25)
  top <- pr[pr$output == "dEDCT48_AB", ]; top <- top[order(-abs(top$PRCC)), ][1:5, ]
  N$prcc_top <- lapply(seq_len(nrow(top)), function(i) list(top$parameter[i], top$PRCC[i]))
  topd <- pr[pr$output == "AEdisc_A", ]; topd <- topd[order(-abs(topd$PRCC)), ][1:3, ]
  N$prcc_top_disc <- lapply(seq_len(nrow(topd)), function(i) list(topd$parameter[i], topd$PRCC[i]))
  N$evpi_E <- et$dEDCT48_AB[et$parameter == "EVPI (all)"]
  fit <- stats::lm(dEDCT48_AB ~ tum.p_attr + tum.p_attr_mid, data = gs)
  N$attr_slope <- unname(stats::coef(fit)[2])     # months of A-B gain per unit attrition (> 12 months)
  # value of information for the decision rule's efficacy threshold: net benefit = gain - 3 months
  evppi_bins <- function(x, nb, nbins = 8) { o <- order(x); bins <- split(o, cut(seq_along(o), nbins, labels = FALSE))
    cond <- numeric(length(nb)); for (b in bins) cond[b] <- mean(nb[b]); mean(pmax(cond, 0)) - max(mean(nb), 0) }
  pnames <- grep(".", names(gs), fixed = TRUE, value = TRUE)
  nb3 <- gs$dEDCT48_AB - 3
  e3 <- data.frame(parameter = c("EVPI (all)", pnames),
                   value = c(mean(pmax(nb3, 0)) - max(mean(nb3), 0), sapply(pnames, function(pp) evppi_bins(gs[[pp]], nb3))))
  e3 <- e3[order(-e3$value), ]
  names(e3) <- c("parameter", "EVPPI (months), gain of A over B vs 3 months")
  wr("evppi3", md(e3, default = "%.3f"))
  N$evpi3 <- e3[e3$parameter == "EVPI (all)", 2]
  e3p <- e3[e3$parameter != "EVPI (all)", ]
  N$evppi3_top <- lapply(1:3, function(i) list(e3p$parameter[i], e3p[i, 2]))
}
# ---------------- power
if (file.exists(file.path(RES, "power.csv"))) {
  pw <- utils::read.csv(file.path(RES, "power.csv"), stringsAsFactors = FALSE)
  pt <- stats::reshape(pw[, c("contrast", "endpoint", "n_per_arm", "power")], idvar = c("contrast", "endpoint"),
                       timevar = "n_per_arm", direction = "wide")
  names(pt) <- sub("^power\\.", "n = ", names(pt))
  for (c in grep("^n = ", names(pt), value = TRUE)) pt[[c]] <- 100 * pt[[c]]
  gains <- if ("gain_EDCT48" %in% names(pw)) tapply(pw$gain_EDCT48, pw$contrast, function(x) x[1]) else NULL
  if (!is.null(gains)) pt <- cbind(pt[, 1:2], `EDCT gain (mo)` = as.numeric(gains[pt$contrast]), pt[, -(1:2)])
  wr("power", md(pt, list(`EDCT gain (mo)` = "%.1f"), default = "%.0f"))
  N$power <- stats::setNames(as.list(pw$power), paste(pw$contrast, pw$endpoint, pw$n_per_arm, sep = "|"))
  # smallest tested size with >= 80% power (EDCT log-rank), per design
  n80 <- tapply(seq_len(nrow(pw)), pw$contrast, function(ix) { x <- pw[ix, ]; x <- x[x$endpoint == "EDCT_logrank" & x$power >= 0.8, ]
    if (nrow(x)) min(x$n_per_arm) else NA })
  N$n80 <- as.list(n80)
  if (!is.null(gains)) N$power_gain <- as.list(gains)
}
# ---------------- structural
if (file.exists(file.path(RES, "structural.csv"))) {
  st <- utils::read.csv(file.path(RES, "structural.csv"), stringsAsFactors = FALSE)
  stage <- c("no S6K-ER cross-talk" = "stage 2", "no adaptive escape" = "stages 1 + 3 with persevERA (SERD potency fixed)",
             "reversible escape" = "stages 1 + 3 with persevERA (SERD potency fixed)")
  # base model: joint loss at the final parameters (04b; the structural refits start from them), else at the end of the joint fit
  c1 <- if (file.exists(file.path(RES, "calib_stage13_final.json"))) list(loss = load_json("calib_stage13_final.json")$final_loss)
        else if (file.exists(file.path(RES, "calib_stage13.json"))) load_json("calib_stage13.json") else load_json("calib_stage1.json")
  c2 <- load_json("calib_stage2.json")
  st$recal <- sapply(seq_len(nrow(st)), function(i) {
    s_ <- st$structure[i]
    if (s_ == "base") sprintf("stages 1 + 3 with persevERA: %.2f; stage 2: %.2f", c1$loss, c2$loss)
    else if (s_ %in% names(stage) && is.finite(st$calib_loss[i])) sprintf("%s: %.2f", stage[[s_]], st$calib_loss[i])
    else "not recalibrated" })
  wr("structural", md(data.frame(structure = st$structure, contrast = st$contrast, `dEDCT RMST48 (mo)` = st$dEDCT48,
                                 `dOS RMST60 (mo)` = st$dOS60, `recalibration loss` = st$recal, check.names = FALSE), default = "%+.2f"))
  ab <- st[st$contrast == "A-B", ]
  N$struct_range_AB <- range(ab$dEDCT48); N$struct_n <- nrow(ab); N$struct_n3 <- sum(ab$dEDCT48 >= 3)
  N$struct_AB <- stats::setNames(as.list(ab$dEDCT48), ab$structure)
  N$stage13_loss <- c1$loss; N$stage13_info <- c1$info
}
# ---------------- alternative AE-discontinuation models (06b)
if (file.exists(file.path(RES, "disc_sensitivity.csv"))) {
  ds <- utils::read.csv(file.path(RES, "disc_sensitivity.csv"), check.names = FALSE, stringsAsFactors = FALSE); rownames(ds) <- ds$arm
  cols <- c("EDCT_RMST48", "OS_RMST60", "ae_disc_cdk", "ae_disc_eve", "ae_disc_et", "ae_disc_any", "RDI_eve", "QTWiST")
  VLAB <- c(base = "falling hazard; everolimus scaled by exposure", const = "constant hazards refitted to the same rates",
            gamma0 = "everolimus discontinuation independent of exposure")
  rws <- list()
  for (b in c("A", "A_eve5", "B", "C", "D2")) for (vn in names(VLAB)) {
    key <- paste0(b, "_", vn)
    if (!(key %in% rownames(ds))) next
    r <- data.frame(arm = b, `discontinuation model` = VLAB[[vn]], check.names = FALSE)
    for (c in cols) r[[c]] <- if (startsWith(c, "ae_") || startsWith(c, "RDI")) 100 * ds[key, c] else ds[key, c]
    rws[[length(rws) + 1]] <- r
    N[[paste0("disc_", vn, "_", b)]] <- ds[key, "ae_disc_any"]
    N[[paste0("edct_", vn, "_", b)]] <- ds[key, "EDCT_RMST48"]
  }
  dt <- do.call(rbind, rws)
  names(dt) <- c("arm", "discontinuation model", "EDCT RMST48", "OS RMST60", "disc CDK4/6i %", "disc eve %", "disc ET/SERD %",
                 "disc any %", "RDI eve %", "QTWiST")
  wr("disc_sensitivity", md(dt, default = "%.1f"))
}
# ---------------- reference headline numbers
for (k in c(arms6, "A_let", "A_eve5", "A_eve2.5", "A_dexstrict", "A_eve5_dexstrict", "B_dexstrict", "A_nodex", "B_trial", "C_trial", "A_abema_eve5",
            "A_abema_eve10"))
  if (k %in% rownames(ref)) for (m_ in c("EDCT_RMST48", "EDCT_median", "OS_RMST60", "OS_median", "PFS1_median", "ae_disc_any", "QTWiST",
                                         "sympt_g3", "frac_line2", "ae_disc_eve", "ae_disc_cdk", "ae_disc_et"))
    N[[paste0("ref_", k, "_", m_)]] <- ref[k, m_]
dEi <- (pmin(A$edct, 48 * MONTH) - pmin(B$edct, 48 * MONTH)) / MONTH
N$ind_benefit <- mean(dEi >= 3); N$ind_harm <- mean(dEi <= -3)
dOi <- (pmin(A$os, 60 * MONTH) - pmin(B$os, 60 * MONTH)) / MONTH
N$ind_os_benefit <- mean(dOi >= 3); N$ind_os_harm <- mean(dOi <= -3)
# ---------------- decision tree: who benefits (rpart, depth 3, at least 200 patients per leaf)
f <- readRDS(file.path(RES, "ref_features.rds"))
X <- data.frame(b0 = f$b0, PI3K_altered = as.numeric(f$Ptrunc), log10_ESR1_fraction = log10(f$fE0), log10_escape_fraction = log10(f$fR0),
                kmax = f$kmax, log10_tumor_cells = log10(f$N0), everolimus_exposure = f$expo_eve, AE_susceptibility = f$frailty, age = f$age)
X$y <- factor(ifelse(dEi >= 3, "benefit", ifelse(dEi <= -3, "harm", "neutral")))
tr <- rpart::rpart(y ~ ., data = X, method = "class", control = rpart::rpart.control(maxdepth = 3, minbucket = 200, cp = 0, xval = 0))
txt <- utils::capture.output(print(tr, digits = 3))
wr("tree", paste(c("```", txt, "```"), collapse = "\n"))
imp <- tr$variable.importance
N$tree_importance <- if (length(imp)) lapply(names(imp)[seq_len(min(4, length(imp)))], function(nm) list(nm, unname(imp[nm] / sum(imp)))) else list()
if (!is.null(f$TF)) N$tf10_frac <- mean(f$TF >= 0.10)
save_json(N, "report_numbers.json")
cat("decision:", N$decision, "\n"); print(D0$crit)
progress_note("10_tables", paste("tables written; decision:", N$decision))
