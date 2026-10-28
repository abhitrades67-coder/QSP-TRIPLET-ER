# Tables of the article, from saved results only:
#   table1.md, table2.md (main tables) and si_tables.md (Supplementary Tables 1-21, captions included) in QSP_MS_TAB
#   (default results/article_tables)
#   5_Supplementary_Data_1.xlsx (inputs with sources and all simulation outputs) in QSP_SUB (default results/)
# Usage: Rscript scripts/12_article_tables.R

if (!exists("ROOT")) ROOT <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), ".."), winslash = "/")
source(file.path(ROOT, "R", "setup.R"))
MS <- Sys.getenv("QSP_MS_TAB", file.path(RES, "article_tables")); dir.create(MS, showWarnings = FALSE, recursive = TRUE)
SUB <- Sys.getenv("QSP_SUB", RES)
N <- load_json("report_numbers.json")
ND <- "–"; GE <- "≥"; LE <- "≤"; MINUS <- "−"

num <- function(x, d = 1) ifelse(is.na(x), "NA", gsub("-", MINUS, formatC(x, format = "f", digits = d), fixed = TRUE))
q3 <- function(x) { x <- x[is.finite(x)]; c(stats::median(x), stats::quantile(x, .05, names = FALSE), stats::quantile(x, .95, names = FALSE)) }
ci <- function(x, d = 1, mult = 1) { q <- q3(mult * x)   # en dash between non-negative bounds, "to" when a bound is negative
  sprintf("%s (%s%s%s)", num(q[1], d), num(q[2], d), if (isTRUE(q[2] < 0)) " to " else ND, num(q[3], d)) }
mdt <- function(df) {
  esc <- function(v) gsub("|", "\\|", as.character(v), fixed = TRUE)
  out <- c(paste0("| ", paste(esc(names(df)), collapse = " | "), " |"), paste0("|", paste(rep("---", ncol(df)), collapse = "|"), "|"))
  for (i in seq_len(nrow(df))) out <- c(out, paste0("| ", paste(esc(unlist(lapply(df[i, ], as.character))), collapse = " | "), " |"))
  paste(out, collapse = "\n")
}
rd <- function(f) utils::read.csv(file.path(RES, f), check.names = FALSE, stringsAsFactors = FALSE)

ARMN <- c(A = "A, upfront triplet", A_let = "A with letrozole instead of giredestrant", A_eve5 = "A with everolimus 5 mg",
          A_eve2.5 = "A with everolimus 2.5 mg", A_ribo_eve2.5 = "A with ribociclib 600 mg + everolimus 2.5 mg",
          A_ribo_eve10 = "A with ribociclib 600 mg + everolimus 10 mg", A_ribo300_eve2.5 = "A with ribociclib 300 mg + everolimus 2.5 mg",
          A_abema_eve10 = "A with abemaciclib + everolimus 10 mg", A_abema_eve5 = "A with abemaciclib + everolimus 5 mg",
          A_tox1.5 = "A, toxicity interaction factor 1.5", A_tox2.0 = "A, toxicity interaction factor 2.0",
          A_eve5_tox2.0 = "A with everolimus 5 mg, interaction factor 2.0", A_nodex = "A without stomatitis prophylaxis",
          A_dexstrict = "A with strict stomatitis prophylaxis", A_eve5_dexstrict = "A with everolimus 5 mg and strict prophylaxis",
          B = "B, sequential", B_ribo = "B with ribociclib", B_abema = "B with abemaciclib", B_nodex = "B without stomatitis prophylaxis",
          B_dexstrict = "B with strict stomatitis prophylaxis", B_trial = "B, trial-like attrition (15%)", B_trial_10 = "B, trial-like attrition (10%)",
          B_trial_20 = "B, trial-like attrition (20%)", C = "C, standard", C_ribo = "C with ribociclib", C_abema = "C with abemaciclib",
          C_ful = "C with fulvestrant in line 2", C_trial = "C, trial-like attrition (15%)", B_ful1L = "B with fulvestrant in line 1",
          C_ful1L = "C with fulvestrant in line 1", D1 = "D1, SERD first", D2 = "D2, everolimus first", D3 = "D3, ctDNA-guided switch")
an <- function(a) ifelse(a %in% names(ARMN), ARMN[a], a)
PN <- c(tum.kmax_med = "Maximal proliferation rate, median (per day)", tum.rd = "Basal death fraction",
        tum.b0_mean = "Endocrine independence b~0~, mean", tum.b0_conc = "Endocrine independence b~0~, concentration",
        tum.log10_fE0 = "Initial ESR1-mutant fraction, median (log10)", tum.log10_fR0 = "Initial CDK4/6i-escape fraction (log10)",
        tum.sigma_R = "Adaptive escape rate (per day)", ph.IC50_c = "CDK4/6 inhibitor potency, IC~50~ (relative exposure)",
        tum.kmax_sd = "Maximal proliferation rate, spread (log scale)", tum.sd_log10_fE0 = "Initial ESR1-mutant fraction, spread (log10)",
        tum.pE_base = "Share of patients with a baseline ESR1 subclone", tum.LOD_pd = "Effective ctDNA detection limit, PADA-1 (allele fraction)",
        tum.lam_NL = "New-lesion hazard at nadir (per day)", ph.uG_ref = "Giredestrant potency index", sig.s_mut = "SERD potency shift against mutant ER",
        ph.uF_ref = "Fulvestrant potency index", ph.I10_e = "mTORC1 inhibition at everolimus 10 mg", sig.eps_g = "Dependence of growth on S6K activity (ε~g~)",
        tum.h_t = "Tumour-burden death hazard coefficient", tum.h_post = "Death hazard after endocrine therapy (per day)",
        tum.kappa = "mTORC1 dependence of adaptive escape (κ)", ph.IC50_e = "Everolimus dose–response IC~50~", ph.gir_ddi = "Giredestrant effect on everolimus exposure (fold)",
        tum.p_attr = "Attrition, progression after more than 12 months", tum.p_attr_mid = "Attrition, progression after 6 to 12 months",
        tum.dc0_med = "Post-endocrine therapy effect (per day)", tox.tox_int = "Three-drug toxicity interaction factor",
        tox.disc_eve_gamma = "Exposure dependence of everolimus discontinuation", tox.disc_tau = "Time constant of the discontinuation hazard (days)",
        tox.disc_floor = "Late discontinuation hazard (share of initial)", tum.early_prog_days = "Early-progression cut-off on line 1 (days)")
pn <- function(p) ifelse(p %in% names(PN), PN[p], p)

es <- rd("ensemble_summary.csv")
piv <- function(metric) { w <- stats::reshape(es[, c("member", "arm", metric)], idvar = "member", timevar = "arm", direction = "wide")
  names(w) <- sub(paste0("^", metric, "\\."), "", names(w)); w[order(w$member), ] }
dcc <- rd("ensemble_disc_components.csv")
bbw <- stats::reshape(dcc[, c("member", "arm", "backbone")], idvar = "member", timevar = "arm", direction = "wide")
names(bbw) <- sub("^backbone\\.", "", names(bbw)); bbw <- bbw[order(bbw$member), ]
ARMS6 <- c("A", "B", "C", "D1", "D2", "D3")

# ---------------------------------------------------------------------------- Table 1
E <- piv("EDCT_RMST48"); O <- piv("OS_RMST60"); Q <- piv("QTWiST")
row <- function(lab, f) { r <- data.frame(Outcome = lab, check.names = FALSE); for (a in ARMS6) r[[paste("Arm", a)]] <- f(a); r }
t1 <- rbind(
  row("EDCT, restricted mean over 48 months (months)", function(a) ci(E[[a]])),
  row("EDCT gain over Arm B (months)", function(a) if (a == "B") "reference" else ci(E[[a]] - E$B)),
  row(paste0("Probability of an EDCT gain ", GE, " 3 months"), function(a) if (a == "B") "reference" else num(mean(E[[a]] - E$B >= 3), 2)),
  row("First-line PFS, median (months)", function(a) ci(es$PFS1_median[es$arm == a])),
  row("Started a second endocrine-based line (%)", function(a) if (a == "A") "none" else ci(es$frac_line2[es$arm == a], 0, 100)),
  row("Overall survival, restricted mean over 60 months (months)", function(a) ci(O[[a]])),
  row("Overall-survival gain over Arm B (months)", function(a) if (a == "B") "reference" else ci(O[[a]] - O$B)),
  row("Q-TWiST over 60 months (months)", function(a) ci(Q[[a]])),
  row("Stopped at least one drug for AEs (%)", function(a) ci(es$ae_disc_any[es$arm == a], 1, 100)),
  row("Stopped the CDK4/6 inhibitor or endocrine agent for AEs (%)", function(a) ci(bbw[[a]], 1, 100)),
  row(paste0("Selected grade ", GE, "3 AE (%)"), function(a) ci(es$sympt_g3[es$arm == a], 1, 100)),
  row(paste0("Grade ", GE, "3 pneumonitis (%)"), function(a) ci(es$pneu_g3[es$arm == a], 1, 100)))
writeLines(mdt(t1), file.path(MS, "table1.md"), useBytes = TRUE)

# ---------------------------------------------------------------------------- Table 2
DI <- piv("ae_disc_any")
dE <- list(rw = E$A - E$B, tr = E$A - E$B_trial); dO <- list(rw = O$A - O$B, tr = O$A - O$B_trial)
discA <- stats::median(DI$A); bbA <- stats::median(bbw$A)
# Q-TWiST criterion in the reference population (as in 10_tables.R), against Arm B with each attrition assumption
RA2 <- readRDS(file.path(RES, "ref_arms.rds")); us2 <- seq(0, 1, 0.1)
qt2 <- function(R, ut, ur) { E_ <- pmin(R$edct, 60 * MONTH); O_ <- pmin(R$os, 60 * MONTH); T_ <- pmin(R$tox_days, E_); mean(ut * T_ + (E_ - T_) + ur * (O_ - E_)) / MONTH }
qshare <- function(b) mean(unlist(lapply(us2, function(ut) sapply(us2, function(ur) qt2(RA2$A, ut, ur) - qt2(RA2[[b]], ut, ur) > 0))))
qpos <- list(rw = qshare("B"), tr = qshare("B_trial"))
stopifnot(abs(qpos$rw - N$qtwist_frac_positive) < 1e-9)
yn <- function(x) if (isTRUE(x)) "yes" else "no"
dec <- function(k, disc) { fav <- mean(dE[[k]] >= 3) >= 0.8 && mean(dO[[k]] > 0) >= 0.7 && disc <= 0.25 && qpos[[k]] >= 0.999
  if (fav) "favoured" else if (mean(dE[[k]] <= 0) >= 0.5 || disc > 0.25) "not favoured" else "uncertain" }
t2 <- data.frame(
  Criterion = c(paste0("Probability of an EDCT gain ", GE, " 3 months (A vs B)"), "Probability of any overall-survival gain (A vs B)",
                "Stopped at least one drug for AEs, Arm A (ensemble median)", "Q-TWiST gain > 0 (share of 121 utility pairs, reference population)",
                "**Decision, pre-specified rule**",
                "Post hoc: stopped the CDK4/6 inhibitor or endocrine agent for AEs, Arm A (ensemble median)",
                "**Decision with the post hoc criterion in place of the third criterion**"),
  Threshold = c(paste(GE, "0.80"), paste(GE, "0.70"), paste(LE, "25%"), "all pairs", "", paste(LE, "25%"), ""),
  `Real-world attrition (primary analysis)` = c(num(mean(dE$rw >= 3), 2), num(mean(dO$rw > 0), 2), paste0(num(100 * discA), "%"), num(qpos$rw, 2),
                                                paste0("**", dec("rw", discA), "**"), paste0(num(100 * bbA), "%"), paste0("**", dec("rw", bbA), "**")),
  `Met` = c(yn(mean(dE$rw >= 3) >= 0.8), yn(mean(dO$rw > 0) >= 0.7), yn(discA <= 0.25), yn(qpos$rw >= 0.999), "", yn(bbA <= 0.25), ""),
  `Trial-like attrition, 15% (secondary analysis)` = c(num(mean(dE$tr >= 3), 2), num(mean(dO$tr > 0), 2), paste0(num(100 * discA), "%"), num(qpos$tr, 2),
                                                       paste0("**", dec("tr", discA), "**"), paste0(num(100 * bbA), "%"), paste0("**", dec("tr", bbA), "**")),
  `Met ` = c(yn(mean(dE$tr >= 3) >= 0.8), yn(mean(dO$tr > 0) >= 0.7), yn(discA <= 0.25), yn(qpos$tr >= 0.999), "", yn(bbA <= 0.25), ""),
  check.names = FALSE)
writeLines(mdt(t2), file.path(MS, "table2.md"), useBytes = TRUE)
stopifnot(sub("favoured", "favored", dec("rw", discA)) == tolower(sub("Triplet ", "", N$decision)))   # same verdict as 10_tables.R

# ---------------------------------------------------------------------------- Supplementary Tables
SI <- character(0)
add <- function(n, title, legend, tab) SI <<- c(SI, sprintf("**Supplementary Table %d | %s** %s", n, title, legend), "", tab, "")
pct1 <- function(x) num(100 * x, 1)
# 1 calibration targets
cal <- rd("calibration_table.csv")
okc <- ifelse(is.na(cal$rel_error_pct), abs(cal$abs_error_points) <= 5, abs(cal$rel_error_pct) <= 15)
isp <- !is.na(cal$abs_error_points)
pu <- ifelse(grepl("difference", cal$target), " points", "%")   # a difference of two proportions is in percentage points
add(1, "Efficacy calibration targets.", "Published values and model values at the final parameters on an independent population of 3,000 virtual patients. Proportions are shown as percentages (the persevERA response difference in percentage points); error is relative for times and hazard ratios and in percentage points for proportions. Tolerance: 15% or 5 percentage points.",
    mdt(data.frame(Trial = cal$trial, Target = cal$target,
                   Published = ifelse(isp, paste0(num(100 * cal$published), pu), trimws(formatC(cal$published, format = "fg", digits = 3))),
                   Model = ifelse(isp, paste0(num(100 * cal$model), pu), formatC(cal$model, format = "fg", digits = 3, flag = "#")),
                   Error = ifelse(isp, paste0(num(cal$abs_error_points), " points"), paste0(num(cal$rel_error_pct), "%")),
                   `Within tolerance` = ifelse(okc, "yes", "no"), check.names = FALSE)))
# 2 toxicity calibration
tc <- rd("tox_calibration.csv"); tc <- tc[!startsWith(tc$target, "info"), ]
TN <- c(neut_g3_palbo = "Grade 3-4 neutropenia, palbociclib + letrozole (PALOMA-2)", neut_g3_ribo = "Grade 3-4 neutropenia, ribociclib + letrozole (MONALEESA-2)",
        neut_g3_abema = "Grade 3-4 neutropenia, abemaciclib + AI (MONARCH 3)", stom_any_eve = "Stomatitis, any grade, ET + everolimus (evERA control arm)",
        stom_g3_eve = "Stomatitis, grade 3 (evERA control arm)", pneu_any_eve = "Pneumonitis, any grade (evERA control arm)",
        hyperg_g3_eve = "Hyperglycaemia, grade 3-4 (evERA control arm)", diar_g3_eve = "Diarrhoea, grade 3 (evERA control arm)",
        anemia_g3_eve = "Anaemia, grade 3-4 (evERA control arm)", diar_g2_abema = "Diarrhoea, grade 2 or higher, abemaciclib (MONARCH 3)",
        diar_g3_abema = "Diarrhoea, grade 3, abemaciclib (MONARCH 3)", alt_g3_ribo = "ALT elevation, grade 3-4, ribociclib (MONALEESA-2)",
        qtc_ribo = "QTcF above 480 ms, ribociclib (MONALEESA-2)", brady_gir = "Bradycardia, giredestrant + everolimus (evERA experimental arm)")
add(2, "Event-toxicity calibration.", "Published incidence and model incidence in toxicity-only runs with each regimen's published treatment duration.",
    mdt(data.frame(Target = ifelse(tc$target %in% names(TN), TN[tc$target], tc$target), `Published (%)` = pct1(tc$observed_target),
                   `Model (%)` = pct1(tc$model), check.names = FALSE)))
# 3 discontinuation calibration
dc <- rd("disc_calibration.csv")
add(3, "Calibration of discontinuation for adverse events.", "Published rates, number of patients, model rates and the deviation in standard errors of the published rate.",
    mdt(data.frame(`Published rate` = dc$label, `Published (%)` = pct1(dc$target), n = dc$n, `Model (%)` = pct1(dc$model),
                   `Deviation (s.e.)` = num(dc$z), Source = dc$source, check.names = FALSE)))
# 4 validation
val <- rd("validation.csv")
ep <- val$endpoint
for (pr in list(c("gir\\+eve", "giredestrant + everolimus"), c("ET\\+eve", "ET + everolimus"), c("\\(esr1\\)", "(ESR1-mutated)"),
                c("\\(all\\)", "(overall)"), c("\\(noesr1\\)", "(no ESR1 mutation)"), c("% exp", "%, experimental arm"), c("% ctl", "%, control arm"),
                c("^stom_any", "stomatitis, any grade"), c("^stom_g3", "stomatitis, grade 3"), c("^pneu_any", "pneumonitis, any grade"),
                c("^pneu_g3", "pneumonitis, grade 3 or higher"), c("^ae_disc_eve", "everolimus stopped for AEs")))
  ep <- gsub(pr[1], pr[2], ep)
vo <- suppressWarnings(as.numeric(val$observed)); vm <- suppressWarnings(as.numeric(val$model))
isprop <- (grepl("<= 0\\.|\\+/-0\\.", val$criterion) | grepl("^(stom|pneu|ae_disc)_|^ORR", val$endpoint)) & is.finite(vo) & vo <= 1   # shown as %
crit <- val$criterion
tol <- suppressWarnings(as.numeric(sub(".*(<= |\\+/-)(0\\.[0-9]+).*", "\\2", crit)))
crit[isprop & is.finite(tol)] <- sprintf("|diff| <= %g points", 100 * tol[isprop & is.finite(tol)])
# descriptive proportions as percentages too (single values or "a / b" pairs); medians in months to one decimal
isfrac <- grepl("prevalence at entry|concurrent progression|^response ", val$endpoint)
pctv <- function(x) vapply(x, function(s) { p <- strsplit(s, " / ", fixed = TRUE)[[1]]; v <- suppressWarnings(as.numeric(p))
  if (all(is.finite(v)) && all(v <= 1)) paste(paste0(num(100 * v, 1), "%"), collapse = " / ") else s }, "", USE.NAMES = FALSE)
ismed <- grepl("^median (PFS|time)", val$endpoint)
fmtv <- function(x, v, pr) ifelse(pr & is.finite(v), paste0(num(100 * v, 1), "%"), ifelse(isfrac, pctv(x), ifelse(ismed & is.finite(v), num(v, 1), x)))
vci <- ifelse(is.na(val$obs_CI), "", sub("^0\\.57-1\\.3$", "0.57-1.30", val$obs_CI))   # SERENA-6 overall survival: interval as published
add(4, "Held-out validation.", "Results not used in any calibration step, with the acceptance criterion fixed in code before the validation simulations. Proportions are shown as percentages. Rows without a pass or fail are calibrated or descriptive values shown for context; calibrated values come from these full-model emulations and can differ from the fitting runs in Supplementary Tables 2 and 3.",
    mdt(data.frame(Trial = val$trial, Endpoint = ep, Published = fmtv(val$observed, vo, isprop),
                   `Published 95% CI` = vci,
                   Model = fmtv(val$model, vm, isprop), Criterion = crit,
                   Met = ifelse(is.na(val$passed), "", ifelse(val$passed, "yes", "no")), check.names = FALSE)))
# 5 external validation
exv <- rd("external_validation.csv")
calx <- exv$criterion == "calibrated"   # calibration targets, not counted
add(5, "External validation against SONIA and persevERA.", "Criteria fixed in code before the validation simulations; both trials had been published before the model was built. persevERA PFS and discontinuation are calibration targets and are not counted. SONIA reports that 9% of patients in the CDK4/6 inhibitor-first group stopped the CDK4/6 inhibitor while continuing endocrine therapy; the tolerance of 6% to 15% for the model's rate of stopping it for AEs was set by the authors. Model values come from a separate population of 6,000 virtual patients and can differ from the calibration check in Supplementary Table 1 (3,000 patients); medians can differ by one assessment interval, because progression is detected at scheduled assessments. TF, ctDNA tumour fraction; the two definitions of high tumour DNA give hazard ratios that agree to two decimals because the model's effect hardly varies with tumour DNA (for example 0.904 and 0.902 in the high groups).",
    mdt(data.frame(Trial = exv$trial, Item = exv$item, Published = exv$published, `95% CI` = ifelse(is.na(exv$published_ci), "", exv$published_ci),
                   Model = exv$model, Criterion = ifelse(calx, "", gsub("+/-", "±", exv$criterion, fixed = TRUE)), Met = ifelse(calx, "calibration target", ifelse(exv$met == "yes", "yes", ifelse(exv$met == "NO", "no", ""))),
                   check.names = FALSE)))
# 6 outcomes by strategy (ensemble)
mets <- list(c("EDCT_RMST48", "EDCT, restricted mean over 48 months (months)", 1, 1), c("EDCT_median", "EDCT, median (months)", 1, 1),
             c("OS_RMST60", "Overall survival, restricted mean over 60 months (months)", 1, 1), c("OS_median", "Overall survival, median (months)", 1, 1),
             c("OS_60", "Alive at 5 years (%)", 0, 100), c("PFS1_median", "First-line PFS, median (months)", 1, 1), c("ORR1", "First-line objective response (%)", 0, 100),
             c("frac_line2", "Started a second endocrine-based line (%)", 0, 100), c("PFS_line2_median", "Second-line PFS, median (months)", 1, 1),
             c("ae_disc_any", "Stopped at least one drug for AEs (%)", 1, 100), c("sympt_g3", "Selected grade 3 or higher AE (%)", 1, 100),
             c("neut_g3", "Grade 3 or higher neutropenia (%)", 0, 100), c("stom_g3", "Grade 3 stomatitis (%)", 1, 100), c("pneu_g3", "Grade 3 or higher pneumonitis (%)", 1, 100),
             c("QTWiST", "Q-TWiST over 60 months (months)", 1, 1), c("ESR1_detect_24m", "ESR1 mutation detected in ctDNA by 24 months (%)", 0, 100),
             c("frac_E_at_EDCT", "ESR1-mutant share of the tumour at strategy failure (%)", 0, 100),
             c("frac_R_at_EDCT", "CDK4/6i-escape share of the tumour at strategy failure (%)", 0, 100))
t6 <- do.call(rbind, lapply(mets, function(m) { r <- data.frame(Outcome = m[2], check.names = FALSE)
  for (a in ARMS6) { x <- es[[m[1]]][es$arm == a]; r[[paste("Arm", a)]] <- if (m[1] == "frac_line2" && a == "A") "none" else if (all(!is.finite(x))) "NA" else ci(x, as.integer(m[3]), as.numeric(m[4])) }; r }))
add(6, "Outcomes by strategy across the parameter ensemble.", "Median (90% interval) over the 60 accepted parameter sets, each with 6,000 virtual patients per arm. None and NA: no second endocrine-based line in Arm A. Selected grade 3 or higher AE: stomatitis, pneumonitis, hyperglycaemia, diarrhoea or anaemia of grade 3 or higher.", mdt(t6))
# 7 paired contrasts
ctr <- list(c("A", "B"), c("A", "C"), c("B", "C"), c("D1", "B"), c("D2", "B"), c("D3", "B"), c("A", "D1"), c("A", "D3"), c("A", "B_trial"), c("A_let", "B"),
            c("A_eve5", "B"), c("A_eve2.5", "B"), c("A_tox1.5", "B"))
t7 <- do.call(rbind, lapply(ctr, function(cc) { a <- cc[1]; b <- cc[2]; g <- E[[a]] - E[[b]]; o <- O[[a]] - O[[b]]; qq <- Q[[a]] - Q[[b]]; dd <- DI[[a]] - DI[[b]]
  data.frame(Contrast = paste(an(a), "vs", an(b)), `EDCT gain (months)` = ci(g), `P(gain ${GE} 3 months)` = num(mean(g >= 3), 2), `P(gain > 0)` = num(mean(g > 0), 2),
             `OS gain (months)` = ci(o), `P(OS gain > 0)` = num(mean(o > 0), 2), `Q-TWiST gain (months)` = ci(qq), `Difference in stopping for AEs (points)` = ci(dd, 1, 100),
             check.names = FALSE) }))
names(t7) <- sub("${GE}", GE, names(t7), fixed = TRUE)
add(7, "Paired strategy contrasts across the parameter ensemble.", "Median (90% interval) of the paired differences over the 60 parameter sets; probabilities are shares of parameter sets. Trial-like attrition: 15% in Arm B.", mdt(t7))
# 8 discontinuation by drug
carms <- intersect(c("A", "A_let", "A_eve5", "A_eve2.5", "A_tox1.5", "B", "B_trial", "C", "D1", "D2", "D3"), unique(dcc$arm))
t8 <- do.call(rbind, lapply(carms, function(a) { g <- dcc[dcc$arm == a, ]
  data.frame(Strategy = an(a), `Any drug` = ci(g$any, 1, 100), `CDK4/6 inhibitor` = ci(g$cdk, 1, 100), Everolimus = ci(g$eve, 1, 100),
             `Endocrine agent` = ci(g$et, 1, 100), `CDK4/6 inhibitor or endocrine agent` = ci(g$backbone, 1, 100), `Everolimus only` = ci(g$eve_only, 1, 100), check.names = FALSE) }))
gdi <- sapply(load_json("ensemble.json")$members, function(m) m$params$ph$gir_ddi)   # giredestrant effect on everolimus exposure per set
r8 <- stats::cor(gdi, DI$A)
add(8, "Discontinuation for adverse events by drug.", sprintf("Percentage of patients who stopped each drug over the whole treatment period; median (90%% interval) over the 60 parameter sets. A stopped endocrine agent is replaced by another one in the model, so a switch counts as a stop. In Arm A, the spread mainly reflects the assumed giredestrant effect on everolimus exposure, sampled from 1.0- to 1.4-fold (correlation with stopping any drug %.2f); without this effect, %.0f%% of patients stopped any drug in the reference population.",
                                                             r8, 100 * N$ref_A_ae_disc_any), mdt(t8))
# 9 variants (reference population)
ref <- rd("ref_summary.csv"); rownames(ref) <- ref$arm
va <- intersect(c("A", "A_let", "A_eve5", "A_eve2.5", "A_ribo_eve2.5", "A_ribo_eve10", "A_ribo300_eve2.5", "A_abema_eve10", "A_abema_eve5", "A_tox1.5", "A_tox2.0",
                  "A_eve5_tox2.0", "A_nodex", "A_dexstrict", "A_eve5_dexstrict", "B", "B_ribo", "B_abema", "B_nodex", "B_dexstrict", "B_trial", "B_trial_10",
                  "B_trial_20", "C", "C_ribo", "C_abema", "C_ful", "C_trial", "B_ful1L", "C_ful1L", "D1", "D2", "D3"), ref$arm)
t9 <- data.frame(Strategy = an(va), `EDCT, RMST 48 months` = num(ref[va, "EDCT_RMST48"]), `OS, RMST 60 months` = num(ref[va, "OS_RMST60"]),
                 `First-line PFS, median` = num(ref[va, "PFS1_median"]), `Stopped any drug for AEs (%)` = pct1(ref[va, "ae_disc_any"]),
                 `Selected grade ≥3 AE (%)` = pct1(ref[va, "sympt_g3"]), `Grade ≥3 neutropenia (%)` = pct1(ref[va, "neut_g3"]),
                 `Grade 3 stomatitis (%)` = pct1(ref[va, "stom_g3"]), `Grade ≥3 pneumonitis (%)` = pct1(ref[va, "pneu_g3"]),
                 `Everolimus dose intensity (%)` = ifelse(is.na(ref[va, "RDI_eve"]), "", pct1(ref[va, "RDI_eve"])), `Q-TWiST` = num(ref[va, "QTWiST"]), check.names = FALSE)
add(9, "Triplet variants and comparators in the reference population.", "Calibrated parameters, 6,000 virtual patients per arm; times in months; percentages are cumulative over all endocrine-based lines. Selected grade ≥3 AE: stomatitis, pneumonitis, hyperglycaemia, diarrhoea or anaemia of grade 3 or higher.", mdt(t9))
# 10 threshold sensitivity
te <- paste(readLines(file.path(TAB, "threshold_efficacy.md"), encoding = "UTF-8", warn = FALSE), collapse = "\n")
td <- paste(readLines(file.path(TAB, "threshold_decision.md"), encoding = "UTF-8", warn = FALSE), collapse = "\n")
hz <- do.call(rbind, lapply(c(36, 48, 60), function(h) { w <- piv(paste0("EDCT_RMST", h)); g <- w$A - w$B; gt <- w$A - w$B_trial
  data.frame(`Horizon (months)` = if (h == 48) "48 (pre-specified)" else as.character(h), `Gain, real-world attrition (months)` = ci(g),
             `P(gain ${GE} 3 months)` = num(mean(g >= 3), 2), `Gain, trial-like attrition (months)` = ci(gt),
             `P(gain ${GE} 3 months), trial-like` = num(mean(gt >= 3), 2), check.names = FALSE) }))
names(hz) <- gsub("${GE}", GE, names(hz), fixed = TRUE)
te <- gsub("\\(([0-9][0-9.]*) to ([0-9][0-9.]*)\\)", paste0("(\\1", ND, "\\2)"), te)   # en dash when no bound is negative
SI <- c(SI, "**Supplementary Table 10 | Threshold and horizon sensitivity (post hoc).** a, probability (share of the 60 parameter sets) that the EDCT gain of Arm A over Arm B reaches each efficacy bar, and the largest bar reached with probability 0.80 or more. b, decision for each combination of efficacy bar and tolerability limit (stopping of any drug, ensemble median), with the other criteria as pre-specified; selected limits are shown, and Fig. 5c shows the full grid. c, gain of Arm A over Arm B in the restricted mean of EDCT over 36, 48 and 60 months (median and 90% interval over the 60 parameter sets) and the probability of a gain of at least 3 months.", "",
        "a", "", gsub(">=", GE, te, fixed = TRUE), "", "b", "", gsub("favored", "favoured", td, fixed = TRUE), "", "c", "", mdt(hz), "")
# 11 secondary criterion
tsec <- paste(readLines(file.path(TAB, "decision_secondary.md"), encoding = "UTF-8", warn = FALSE), collapse = "\n")
tsec <- gsub("Triplet not favored", "Not favoured", gsub("Triplet favored", "Favoured", tsec, fixed = TRUE), fixed = TRUE)
for (pr in list(c("| criterion | value, real-world attrition | met | value, trial-like attrition | met (trial-like) |",
                  "| Criterion | Real-world attrition | Met | Trial-like attrition, 15% | Met |"),
                c("P(dEDCT48 A-B >= 3 months) >= 0.80", "Probability of an EDCT gain >= 3 months (A vs B) >= 0.80"),
                c("P(dOS60 A-B > 0) >= 0.70", "Probability of any overall-survival gain (A vs B) >= 0.70"),
                c("AE-related discontinuation of the backbone (CDK4/6 inhibitor or endocrine agent) in A <= 25% (ensemble median)",
                  "Stopped the CDK4/6 inhibitor or endocrine agent for AEs, Arm A (ensemble median) <= 25%"),
                c("Q-TWiST A-B > 0 across utility thresholds (fraction of 121 utility pairs)", "Q-TWiST gain > 0 (share of 121 utility pairs, reference population)"),
                c("| Uncertain |", "| uncertain |")))
  tsec <- gsub(pr[1], pr[2], tsec, fixed = TRUE)
add(11, "Decision rule with the post hoc secondary tolerability criterion.", "Discontinuation of the backbone (the CDK4/6 inhibitor or the endocrine agent) replaces discontinuation of any drug; the other criteria are as pre-specified. Defined after the results were known.",
    gsub(">=", GE, gsub("<=", LE, tsec, fixed = TRUE), fixed = TRUE))
# 12 subgroups
sg <- rd("ensemble_subgroups.csv")
sgn <- unique(sg$subgroup)
SGN <- c("Tumor" = "Tumour", "low tercile" = "lowest third", "high tercile" = "highest third", "dependence mid" = "dependence, middle third")
t12 <- do.call(rbind, lapply(sgn, function(s) { g <- sg[sg$subgroup == s & sg$contrast == "A-B", ]; h <- sg[sg$subgroup == s & sg$contrast == "A-B_trial", ]
  nm <- s; for (k in names(SGN)) nm <- gsub(k, SGN[[k]], nm, fixed = TRUE)
  data.frame(Subgroup = gsub(">=", GE, nm, fixed = TRUE), `Patients per parameter set (median)` = format(round(stats::median(g$n)), big.mark = ","), `EDCT gain, real-world attrition` = ci(g$dEDCT48),
             `P(gain ${GE} 3 months)` = num(mean(g$dEDCT48 >= 3), 2), `EDCT gain, trial-like attrition` = ci(h$dEDCT48), `P(gain ${GE} 3 months), trial-like` = num(mean(h$dEDCT48 >= 3), 2),
             `OS gain, real-world attrition` = ci(g$dOS60), check.names = FALSE) }))
names(t12) <- gsub("${GE}", GE, names(t12), fixed = TRUE)
# reference population (as in 06_insilico_trial.R): how much of the slowest-proliferating third the 48-month horizon truncates
ra12 <- readRDS(file.path(RES, "ref_arms.rds")); vp12 <- make_vpop(load_params(PFILES), length(ra12$A$edct), 1001)
qk <- stats::quantile(vp12$kmax, c(1, 2) / 3, names = FALSE); slow <- vp12$kmax <= qk[1]; fast <- vp12$kmax > qk[2]
g72 <- function(i) mean(pmin(ra12$A$edct[i], 72 * MONTH) - pmin(ra12$B$edct[i], 72 * MONTH)) / MONTH
HZP <- list(slowA = mean(ra12$A$edct[slow] > 48 * MONTH), slowB = mean(ra12$B$edct[slow] > 48 * MONTH), slow72 = g72(slow), fast72 = g72(fast))
save_json(HZP, "horizon_proliferation.json")
hzp <- sprintf(" In the reference population, %.0f%% (Arm A) and %.0f%% (Arm B) of the slowest-proliferating third were still in endocrine-based disease control at 48 months; restricted to 72 months, the gain was %.1f months in the slowest third and %.1f months in the fastest third (post hoc).",
               100 * HZP$slowA, 100 * HZP$slowB, HZP$slow72, HZP$fast72)
ti <- N$tree_importance; ti1 <- if (is.matrix(ti)) ti[1, ] else unlist(ti[[1]]); stopifnot(ti1[1] == "kmax")
add(12, "Subgroup contrasts, Arm A versus Arm B.", paste0("Median (90% interval) over the 60 parameter sets; gains in months. Thirds are terciles of the virtual population. ",
    sprintf("In a classification tree of individual benefit (paired gain of 3 months or more), harm (loss of 3 months or more) or neither in the reference population (depth 3, at least 200 patients per leaf), the maximal proliferation rate carried %.0f%% of the variable importance.", 100 * as.numeric(ti1[2])), hzp), mdt(t12))
# 13 power
pw <- rd("power.csv")
pt <- stats::reshape(pw[, c("contrast", "endpoint", "n_per_arm", "power")], idvar = c("contrast", "endpoint"), timevar = "n_per_arm", direction = "wide")
names(pt) <- sub("^power\\.", "n = ", names(pt))
nn <- grepl("^n = ", names(pt)); names(pt)[nn] <- paste0("n = ", format(as.integer(sub("^n = ", "", names(pt)[nn])), big.mark = ",", trim = TRUE))
for (c in grep("^n = ", names(pt), value = TRUE)) pt[[c]] <- ifelse(is.na(pt[[c]]), "NA", num(100 * pt[[c]], 0))
gains <- tapply(pw$gain_EDCT48, pw$contrast, function(x) x[1])
EPN <- c(EDCT_logrank = "EDCT, log-rank", EDCT_RMST36 = "EDCT, RMST at 36 months", OS_logrank = "OS, log-rank")
pt <- data.frame(Design = gsub("tumor", "tumour", gsub(">=", GE, pt$contrast, fixed = TRUE)), Endpoint = EPN[pt$endpoint], `EDCT gain (months)` = num(gains[pt$contrast]),
                 pt[, grep("^n = ", names(pt))], check.names = FALSE)
add(13, "Simulated power of a two-arm randomised trial (%).", "1,000 virtual trials per design drawn from the reference population; 1:1 randomization, 24-month accrual, 36-month minimum follow-up, two-sided test at 0.05. n: patients per arm. NA: more patients than the eligible virtual population.", mdt(pt))
# 14 PRCC, 15 EVPPI
prc <- rd("gsa_prcc.csv")
outs <- c(dEDCT48_AB = "EDCT gain, A vs B", dOS60_AB = "OS gain, A vs B", dQTW_AB = "Q-TWiST gain, A vs B", AEdisc_A = "Stopping for AEs, Arm A")
pars <- unique(prc$parameter)
W <- sapply(names(outs), function(o) sapply(pars, function(pp) { r <- prc[prc$output == o & prc$parameter == pp, ]
  paste0(gsub("-", MINUS, sprintf("%+.2f", r$PRCC), fixed = TRUE), if (isTRUE(r$p_value < 0.05)) "*" else "") }))
ordv <- sapply(pars, function(pp) abs(prc$PRCC[prc$output == "dEDCT48_AB" & prc$parameter == pp]))
gs <- rd("gsa_samples.csv"); pb <- load_params(PFILES); g3 <- function(x) formatC(x, format = "g", digits = 3)
cv <- sapply(pars, function(pp) { v <- get_param(pb, pp); if (is.null(v)) "" else g3(v) }); sr <- sapply(pars, function(pp) paste(g3(min(gs[[pp]])), "to", g3(max(gs[[pp]]))))
t14 <- data.frame(Parameter = pn(pars), cv, sr, W, check.names = FALSE)[order(-ordv), ]; names(t14) <- c("Parameter", "Base-case value", "Sampled range", unname(outs))
add(14, "Global sensitivity analysis: partial rank correlation coefficients.", "160 Latin hypercube samples of 22 parameters over ranges wider than those of the parameter ensemble, 1,200 virtual patients each; the samples were not filtered by fit to the calibration targets, and the toxicity settings were mostly less favourable than the calibrated values. *P < 0.05.", mdt(t14))
ev <- rd("gsa_evppi.csv")
e3 <- utils::read.csv(text = paste(readLines(file.path(TAB, "evppi3.md"))[-2], collapse = "\n"), sep = "|", strip.white = TRUE, check.names = FALSE)
e3 <- e3[, !grepl("^V|^$", names(e3)) & nzchar(names(e3))]
names(e3) <- c("Parameter", "EVPPI (months)"); e3[[2]] <- sprintf("%.3f", as.numeric(e3[[2]]))
e3$Parameter <- ifelse(e3$Parameter == "EVPI (all)", "All parameters (EVPI)", pn(e3$Parameter))
add(15, "Expected value of partial perfect information.", "Value of information, in months of EDCT, for the decision whether the gain of Arm A over Arm B reaches 3 months (net benefit = gain minus 3 months), estimated by regression over the global sensitivity samples.", mdt(e3))
# 16 structural
st <- rd("structural.csv")
SN <- c(base = "Base model", `kappa = 0` = "κ fixed at 0", `kappa = 1` = "κ fixed at 1", `burden-only death hazard` = "Tumour-burden-only death hazard",
        `constant attrition 10%` = "Constant 10% attrition", `no S6K-ER cross-talk` = "No S6K–ER cross-talk (stage 2 recalibrated)",
        `no adaptive escape` = "No adaptive escape (stages 1 and 3 recalibrated)", `reversible escape` = "Reversible escape (stages 1 and 3 recalibrated)")
t16 <- data.frame(Structure = ifelse(st$structure %in% names(SN), SN[st$structure], st$structure), Contrast = gsub("-", " vs ", st$contrast, fixed = TRUE),
                  `EDCT gain (months)` = gsub("-", MINUS, sprintf("%+.2f", st$dEDCT48), fixed = TRUE), `OS gain (months)` = gsub("-", MINUS, sprintf("%+.2f", st$dOS60), fixed = TRUE),
                  `Recalibration loss` = ifelse(is.finite(st$calib_loss), num(st$calib_loss), ""), check.names = FALSE)
add(16, "Strategy contrasts under alternative model structures.", "Calibrated parameters, 3,000 virtual patients per arm (fewer than in the reference population, so base values differ slightly from Supplementary Table 9). κ fixed at 0 or 1 and no adaptive escape reproduce the base model because the calibrated rate of drug-induced escape is negligible; the tumour-burden-only death hazard (scale re-matched to MONALEESA-2 overall survival) changes overall survival but not EDCT. Recalibration loss: loss of the refitted stage (base model: stage 2 loss 0.3; joint loss of stages 1 and 3 at the final parameters, 55.9).", mdt(t16))
# 17 discontinuation in full-model emulations
dk <- rd("disc_check.csv")
add(17, "Discontinuation for adverse events in full-model trial emulations.", "Treatment until progression, follow-up equal to each trial's median follow-up. Values can differ slightly from Supplementary Table 4, which uses separate validation emulations. Fitted: the rate is a calibration target of the discontinuation sub-model. Criterion: within the larger of 3 percentage points and 20%.",
    mdt(data.frame(Trial = dk$trial, `Stopped for AEs` = dk$item, `Follow-up (months)` = num(dk$follow_up_mo), `Published (%)` = num(dk$published_pct),
                   `Model (%)` = num(dk$model_pct), Fitted = dk$fitted, Met = ifelse(dk$met == "yes", "yes", "no"), check.names = FALSE)))
# 18 parameters
ens <- load_json("ensemble.json"); cand <- rd("abc_candidates.csv"); acc <- cand[cand$id %in% sapply(ens$members, `[[`, "id"), ]
p <- load_params(PFILES)
pcols <- grep(".", names(cand), fixed = TRUE, value = TRUE)
sig <- function(x) formatC(x, format = "g", digits = 3)
n_elig <- sum(cand$loss <= ens$threshold)
add(18, "Calibrated parameters and their spread across the parameter ensemble.", sprintf("Best fit (calibrated point) and the median and 5th to 95th percentiles over the %d accepted parameter sets, those with the lowest loss among the %d of %d candidates within the chi-square bound. Parameters that the calibration targets cannot inform (κ, everolimus IC~50~, giredestrant interaction, attrition and post-endocrine effect) were sampled from prior ranges; the best-fit attrition values are the ESME-MBC shares. Toxicity and discontinuation parameters were not varied in the ensemble.",
                                                                                         nrow(acc), n_elig, nrow(cand)),
    mdt(data.frame(Parameter = pn(pcols), `Best fit` = sapply(pcols, function(c) sig(get_param(p, c))), `Ensemble median` = sapply(pcols, function(c) sig(stats::median(acc[[c]]))),
                   `5th percentile` = sapply(pcols, function(c) sig(stats::quantile(acc[[c]], .05))), `95th percentile` = sapply(pcols, function(c) sig(stats::quantile(acc[[c]], .95))),
                   check.names = FALSE)))
# 19 scenario 2
s2 <- paste(readLines(file.path(TAB, "scenario2.md"), encoding = "UTF-8", warn = FALSE), collapse = "\n")
s2 <- sub("| arm |", "| Second-line regimen |", s2, fixed = TRUE)
for (pr in list(c("| A'(abema) |", "| Abemaciclib + giredestrant + everolimus |"), c("| A'(palbo cont.) |", "| Palbociclib continued + giredestrant + everolimus |"),
                c("| B' |", "| Giredestrant + everolimus |"), c("| C' |", "| Exemestane + everolimus |"),
                c("HR vs C'", "HR vs exemestane + everolimus"), c("HR vs B'", "HR vs giredestrant + everolimus")))
  s2 <- gsub(pr[1], pr[2], s2, fixed = TRUE)
add(19, "Second-line scenario in the evERA population.", "Patients who progressed after at least 6 months of palbociclib + letrozole. PFS and OS from the start of the second line; reference population.", s2)
# 20 alternative discontinuation models
ds <- rd("disc_sensitivity.csv")
VL <- c(base = "Base (falling hazard; everolimus scaled by exposure)", const = "Constant hazards refitted to the same rates",
        gamma0 = "Everolimus discontinuation independent of exposure")
ds$model <- VL[sub("^.*_", "", ds$arm)]; ds$strategy <- sub("_[^_]+$", "", ds$arm)
add(20, "Outcomes under alternative models of discontinuation for adverse events.", "Reference population; times in months.",
    mdt(data.frame(Strategy = an(ds$strategy), Model = ds$model, `EDCT, RMST 48` = num(ds$EDCT_RMST48), `OS, RMST 60` = num(ds$OS_RMST60),
                   `Stopped CDK4/6i (%)` = pct1(ds$ae_disc_cdk), `Stopped everolimus (%)` = pct1(ds$ae_disc_eve), `Stopped endocrine agent (%)` = pct1(ds$ae_disc_et),
                   `Stopped any drug (%)` = pct1(ds$ae_disc_any), check.names = FALSE)))
# 21 verification
vf <- readLines(file.path(TAB, "verification.md"), encoding = "UTF-8", warn = FALSE)
if (file.exists(file.path(RES, "power_check.json")))   # log-rank test of the power analysis (08_power.R)
  vf <- c(vf, sprintf("| V8 log-rank test of the power analysis vs R 'survival' package (survdiff) | max abs difference in P value | %.2g | yes |",
                      load_json("power_check.json")$logrank_p_difference))
add(21, "Code and numerical verification.", "Checks of the implementation; descriptive rows report predictions without a pass criterion.", gsub("->", "to", paste(vf, collapse = "\n"), fixed = TRUE))
# house style for labels that come verbatim from results files: en dash in numeric ranges, minus sign for a difference,
# >= and <= as symbols, British spelling as in the text
SI <- gsub("(\\d)-(\\d)", paste0("\\1", ND, "\\2"), SI, perl = TRUE)
SI <- gsub(" - ", paste0(" ", MINUS, " "), SI, fixed = TRUE)
SI <- gsub(">=", GE, gsub("<=", LE, SI, fixed = TRUE), fixed = TRUE)
SI <- gsub("\\btumor\\b", "tumour", SI, perl = TRUE)
SI <- gsub("randomiz(ation|ed)", "randomis\\1", SI, perl = TRUE)
SI <- gsub("\\bb0\\b", "b~0~", SI, perl = TRUE)   # b0 as b with subscript 0
SI <- gsub("([Gg]rade) ≥ ([0-9])", "\\1 ≥\\2", SI, perl = TRUE)
# readable sources, words for abbreviations used only in labels, and typography
src <- c(NEJMoa1607303 = "Finn et al. 2016", persevERA_Turner2026_LBA1006 = "Turner et al. 2026 (abstract)",
         jamaoncology_royce_2018 = "Royce et al. 2018", NEJMoa2602457 = "Mayer et al. 2026", goetz2017 = "Goetz et al. 2017")
for (k in names(src)) SI <- gsub(k, src[[k]], SI, fixed = TRUE)
SI <- gsub("(et al\\. 20\\d\\d(?: \\(abstract\\))?) p(\\d+), p(\\d+) \\|", "\\1, pp. \\2, \\3 |", SI, perl = TRUE)
SI <- gsub("(et al\\. 20\\d\\d(?: \\(abstract\\))?) p(\\d+) \\|", "\\1, p. \\2 |", SI, perl = TRUE)
for (pr in list(c("Cabel 2026", "Cabel et al. 2026"), c("+/-", "±"), c("abema + fulv / fulv", "abemaciclib + fulvestrant / fulvestrant"),
                c("gir vs let", "giredestrant vs letrozole"), c("gir − let", "giredestrant − letrozole"), c("gir / let", "giredestrant / letrozole")))
  SI <- gsub(pr[1], pr[2], SI, fixed = TRUE)
for (arm in c("giredestrant + everolimus", "ET + everolimus"))
  SI <- gsub(paste0("18-mo OS % ", arm, " (overall)"), paste0("18-mo OS % ", arm, " (overall, ESR1-enriched)"), SI, fixed = TRUE)
SI <- gsub("\\bpalbo\\b", "palbociclib", SI, perl = TRUE)
SI <- gsub("\\b1L\\b", "first line", SI, perl = TRUE)
SI <- gsub(paste0("(?<![0-9.])", MINUS, "0(\\.0+)(?![0-9])"), "0\\1", SI, perl = TRUE)   # no negative zero
SI <- gsub("(?<=[0-9])e-0*([0-9]+)\\b", paste0(" × 10^", MINUS, "\\1^"), SI, perl = TRUE)
SI <- gsub("(?<=[0-9])e\\+?0*([0-9]+)\\b", " × 10^\\1^", SI, perl = TRUE)
SI <- gsub("(?<=[|( ])-(?=[0-9])", MINUS, SI, perl = TRUE)
SI <- gsub("(?<=[0-9])x(?=[ ,)|])", "×", SI, perl = TRUE); SI <- gsub(" x (?=[0-9])", " × ", SI, perl = TRUE)
SI <- gsub("(?<=\\| )TRUE(?= \\|)", "yes", SI, perl = TRUE); SI <- gsub("(?<=\\| )FALSE(?= \\|)", "no", SI, perl = TRUE)
SI <- gsub("(?<![*A-Za-z_])ESR1(?![*A-Za-z_0-9])", "*ESR1*", SI, perl = TRUE)
SI <- gsub("\\bAI\\b", "aromatase inhibitor", SI, perl = TRUE)   # aromatase inhibitor written out
writeLines(enc2utf8(SI), file.path(MS, "si_tables.md"), useBytes = TRUE)
cat("tables written:", file.path(MS, "table1.md"), file.path(MS, "table2.md"), file.path(MS, "si_tables.md"), "\n")

# ---------------------------------------------------------------------------- gain decomposition (reference population)
# paired EDCT gain of Arm A over Arm B (restricted at 48 months) split by what happened to each patient in Arm B
ra <- readRDS(file.path(RES, "ref_arms.rds"))
gd <- do.call(rbind, lapply(c("B", "B_trial"), function(cmp) {
  b <- ra[[cmp]]; d <- (pmin(ra$A$edct, 48 * MONTH) - pmin(b$edct, 48 * MONTH)) / MONTH
  g <- ifelse(b$in_line2, "started a second endocrine-based line",
       ifelse(b$prog1 & b$pfs1 < 6 * MONTH, "progressed within 6 months; no second endocrine-based line",
       ifelse(b$prog1, "progressed after 6 months or more; no second endocrine-based line (attrition)",
              "no first-line progression (death, or still on the first line at 6 years)")))
  do.call(rbind, lapply(unique(g), function(k) data.frame(comparator = cmp, pathway_in_arm_B = k, share = mean(g == k),
    mean_gain_months = mean(d[g == k]), contribution_months = sum(d[g == k]) / length(d), share_of_gain = sum(d[g == k]) / sum(d),
    stringsAsFactors = FALSE)))
}))
utils::write.csv(gd, file.path(RES, "gain_decomposition.csv"), row.names = FALSE)
for (cmp in c("B", "B_trial")) {
  x <- gd[gd$comparator == cmp & grepl("^progressed", gd$pathway_in_arm_B), ]; y <- gd[gd$comparator == cmp & grepl("^started", gd$pathway_in_arm_B), ]
  cat(sprintf("%s: no second line after progression %.1f%%, mean gain %.2f, share of gain %.1f%%; started %.1f%%, mean gain %.2f, share %.1f%%\n", cmp,
      100 * sum(x$share), sum(x$contribution_months) / sum(x$share), 100 * sum(x$share_of_gain), 100 * y$share, y$mean_gain_months, 100 * y$share_of_gain))
}

# ---------------------------------------------------------------------------- Supplementary Data 1 (Excel)
wb <- openxlsx::createWorkbook(creator = "QSP-TRIPLET-ER")
hs <- openxlsx::createStyle(textDecoration = "bold", fgFill = "#E8F1FA")
sheet <- function(name, df, note = NULL) {
  df[] <- lapply(df, function(v) if (is.character(v)) gsub("\\bAI\\b", "aromatase inhibitor", v, perl = TRUE) else v)
  for (pr in list(c("\\b([Tt])umor", "\\1umour"), c("([Rr])andomiz", "\\1andomis"), c("hyperglycemia", "hyperglycaemia"),
                  c("diarrhea", "diarrhoea"), c("\\b([Aa])nemia", "\\1naemia")))   # British spelling in labels and values
    df[] <- lapply(df, function(v) if (is.character(v)) gsub(pr[1], pr[2], v, perl = TRUE) else v)
  openxlsx::addWorksheet(wb, name)
  r0 <- 1
  if (!is.null(note)) { openxlsx::writeData(wb, name, note, startRow = 1, startCol = 1); r0 <- 3 }
  openxlsx::writeData(wb, name, df, startRow = r0, headerStyle = hs)
  openxlsx::setColWidths(wb, name, cols = seq_len(ncol(df)), widths = "auto")
}
readme <- data.frame(Sheet = c("Inputs_and_sources", "Calibration_targets", "Toxicity_calibration", "Discontinuation_calib", "Validation_held_out",
                               "External_validation", "Ensemble_parameters", "Ensemble_outcomes", "Ensemble_discontinuation", "Ensemble_subgroups",
                               "Reference_outcomes", "Gain_decomposition", "Global_sensitivity", "Power", "Structural", "Discontinuation_models"),
                     Content = c("Every published value used by the model, with its source, page and the value used in the model",
                                 "36 efficacy calibration targets: published and model values", "14 event-toxicity targets: published and model values",
                                 "10 discontinuation rates: published and model values", "Held-out validation criteria and results",
                                 "External validation (SONIA and persevERA overall survival): published and model values; persevERA PFS and discontinuation are calibration targets, shown for context",
                                 "The 60 accepted parameter sets (approximate Bayesian computation): those with the lowest loss among the candidates within the chi-square bound", "Outcomes of every strategy for each parameter set (6,000 virtual patients each); member = parameter set",
                                 "Discontinuation for adverse events by drug, per strategy and parameter set", "Subgroup contrasts per parameter set",
                                 "Outcomes of all strategies and variants in the reference population (calibrated parameters)", "Paired EDCT gain of Arm A over Arm B in the reference population, split by what happened to each patient in Arm B (real-world and trial-like attrition)", "Global sensitivity samples (parameters and outputs)",
                                 "Simulated trial power by design, endpoint and sample size", "Strategy contrasts under alternative model structures",
                                 "Outcomes under alternative discontinuation models"))
openxlsx::addWorksheet(wb, "README")
openxlsx::writeData(wb, "README", "Supplementary Data 1. Model inputs with sources and simulation outputs. Times in months unless stated; proportions as fractions.", startRow = 1)
openxlsx::writeData(wb, "README", readme, startRow = 3, headerStyle = hs); openxlsx::setColWidths(wb, "README", 1:2, widths = c(28, 120))
# inputs with sources: every Markdown table of verification/SOURCE_VERIFICATION.md, under its section heading
sv <- readLines(file.path(ROOT, "verification", "SOURCE_VERIFICATION.md"), encoding = "UTF-8", warn = FALSE)
openxlsx::addWorksheet(wb, "Inputs_and_sources")
r <- 1; head <- ""; i <- 1
while (i <= length(sv)) {
  if (grepl("^## ", sv[i])) head <- sub("^## ", "", sv[i])
  if (grepl("^\\|", sv[i])) {
    j <- i; while (j + 1 <= length(sv) && grepl("^\\|", sv[j + 1])) j <- j + 1
    rows <- sv[i:j]; rows <- rows[!grepl("^\\|[-| :]+\\|$", rows)]
    cells <- lapply(rows, function(x) trimws(strsplit(sub("^\\|", "", sub("\\|$", "", x)), "|", fixed = TRUE)[[1]]))
    nc <- max(lengths(cells)); m <- do.call(rbind, lapply(cells, function(x) c(x, rep("", nc - length(x)))))
    openxlsx::writeData(wb, "Inputs_and_sources", head, startRow = r); r <- r + 1
    df <- as.data.frame(m[-1, , drop = FALSE], stringsAsFactors = FALSE); names(df) <- m[1, ]
    openxlsx::writeData(wb, "Inputs_and_sources", df, startRow = r, headerStyle = hs); r <- r + nrow(df) + 2
    i <- j
  }
  i <- i + 1
}
openxlsx::setColWidths(wb, "Inputs_and_sources", 1:6, widths = c(55, 30, 40, 45, 30, 30))
sheet("Calibration_targets", rd("calibration_table.csv")); sheet("Toxicity_calibration", rd("tox_calibration.csv"))
vh <- val; vh$endpoint <- ep; vh$obs_CI <- sub("^0\\.57-1\\.3$", "0.57-1.30", vh$obs_CI)   # endpoint names as in Supplementary Table 4
sheet("Discontinuation_calib", rd("disc_calibration.csv")); sheet("Validation_held_out", vh)
ext <- exv; ext$criterion[which(calx)] <- "calibration target"; ext$met[which(calx)] <- ""; ext$met[which(ext$met == "NO")] <- "no"   # as in Supplementary Table 5
for (pr in list(c("gir vs let", "giredestrant vs letrozole"), c("gir - let", "giredestrant - letrozole"), c("gir / let", "giredestrant / letrozole")))
  ext$item <- gsub(pr[1], pr[2], ext$item, fixed = TRUE)
sheet("External_validation", ext)
sheet("Ensemble_parameters", acc); sheet("Ensemble_outcomes", es); sheet("Ensemble_discontinuation", dcc); sheet("Ensemble_subgroups", sg)
sheet("Reference_outcomes", ref); sheet("Gain_decomposition", gd); sheet("Global_sensitivity", rd("gsa_samples.csv")); sheet("Power", pw); sheet("Structural", st); sheet("Discontinuation_models", ds[, setdiff(names(ds), c("model", "strategy"))])
dir.create(SUB, showWarnings = FALSE, recursive = TRUE)
openxlsx::saveWorkbook(wb, file.path(SUB, "5_Supplementary_Data_1.xlsx"), overwrite = TRUE)
# openxlsx links every sheet to drawing parts that it does not write. Excel ignores the dangling links, but openpyxl (and
# so pandas.read_excel) refuses the file; drop links and content-type entries whose part is missing.
tidy_xlsx <- function(path) {
  td <- tempfile("xlsx_"); utils::unzip(path, exdir = td)
  rd_ <- function(f) paste(readLines(f, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  wr_ <- function(s, f) writeLines(enc2utf8(s), f, useBytes = TRUE, sep = "")
  for (rel in list.files(file.path(td, "xl", "worksheets", "_rels"), full.names = TRUE)) {
    x <- rd_(rel)
    for (r in regmatches(x, gregexpr("<Relationship [^>]*/>", x))[[1]])
      if (!file.exists(file.path(td, "xl", "worksheets", sub('.*Target="([^"]+)".*', "\\1", r)))) x <- sub(r, "", x, fixed = TRUE)
    wr_(x, rel)
  }
  ct <- file.path(td, "[Content_Types].xml"); x <- rd_(ct)
  for (o in regmatches(x, gregexpr("<Override [^>]*/>", x))[[1]])
    if (!file.exists(file.path(td, sub('.*PartName="/([^"]+)".*', "\\1", o)))) x <- sub(o, "", x, fixed = TRUE)
  wr_(x, ct)
  unlink(path); zip::zip(path, files = list.files(td, recursive = TRUE, all.files = TRUE), root = td, include_directories = FALSE)
  unlink(td, recursive = TRUE)
}
tidy_xlsx(file.path(SUB, "5_Supplementary_Data_1.xlsx"))
cat("written", file.path(SUB, "5_Supplementary_Data_1.xlsx"), "\n")
