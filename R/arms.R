# Treatment strategies (in silico arms) and published-trial emulations.
#
# With `schedules = TRUE` each trial emulation uses the tumor-assessment schedule and dosing reported in its
# publication (see verification/SOURCE_VERIFICATION.md):
#   PALOMA-2 every 12 weeks; BOLERO-2 every 6 weeks; BOLERO-4 every 8 weeks; evERA, acelERA,
#   MONALEESA-2, SERENA-6 every 8 weeks for 18 months then every 12 weeks; postMONARCH and
#   TRINITI-1 every 8 weeks for 12 months then every 12 weeks; MAINTAIN every 12 weeks;
#   TRINITI-1 gave ribociclib continuously (28 of 28 days); postMONARCH had no minimum duration
#   of prior CDK4/6 inhibitor therapy.

R_ <- regimen
NEVER <- 100000L   # 'switch day' that is never reached

SCHED <- list(
  default = c(56L, 56L, 540L, 84L),      # every 8 weeks for 18 months, then every 12 weeks
  q12w = c(84L, 84L, NEVER, 84L),
  q8w = c(56L, 56L, NEVER, 56L),
  q6w = c(42L, 42L, NEVER, 42L),
  q8w_12m = c(56L, 56L, 365L, 84L)
)

# trial emulation settings: every eligible patient enters line 2 (no attrition); eligibility =
# progression after >= 6 months on line 1 (as in evERA)
TRIAL <- list(p_attr = 0.0, p_attr_mid = 0.0, early_prog_days = 183)
trial_arm <- function(lines, ...) c(list(lines = lines), TRIAL, list(...))

strategy_arms <- function(cdk = "palbo", eve_A = 10, cdk_doses = NULL) {
  list(
    A = list(lines = list(R_(cdk, "giredestrant", eve_A, cdk_doses))),
    B = list(lines = list(R_(cdk, "letrozole", NULL, cdk_doses), R_(NULL, "giredestrant", 10))),
    C = list(lines = list(R_(cdk, "letrozole", NULL, cdk_doses), R_(NULL, "exemestane", 10))),
    D1 = list(lines = list(R_(cdk, "giredestrant", NULL, cdk_doses), R_(NULL, "giredestrant", 10))),
    D2 = list(lines = list(R_(cdk, "letrozole", 10, cdk_doses), R_(NULL, "giredestrant", 10))),
    # ctDNA monitoring as in SERENA-6 (every 10 weeks from 6 months, effective detection limit, switch about 5 weeks
    # after detection if imaging shows no progression)
    D3 = c(list(lines = list(R_(cdk, "letrozole", NULL, cdk_doses), R_(NULL, "giredestrant", 10)), d3_switch = TRUE),
           SER6_CT[setdiff(names(SER6_CT), "ass1")])
  )
}

# Scenario 2 (evERA population): line 1 = palbociclib + letrozole, compare line 2
scenario2_arms <- function() {
  l1 <- R_("palbo", "letrozole", NULL)
  list(
    "A'(abema)" = trial_arm(list(l1, R_("abema", "giredestrant", 10))),
    "A'(palbo cont.)" = trial_arm(list(l1, R_("palbo", "giredestrant", 10))),
    "B'" = trial_arm(list(l1, R_(NULL, "giredestrant", 10))),
    "C'" = trial_arm(list(l1, R_(NULL, "exemestane", 10)))
  )
}

calibration_arms <- function(schedules = FALSE) {
  l1c <- R_("palbo", "letrozole", NULL)
  l1a <- R_(NULL, "letrozole", NULL)
  a <- list(
    P2_cdk = list(lines = list(l1c)),
    P2_ai = list(lines = list(l1a)),
    B2_eve = trial_arm(list(l1a, R_(NULL, "exemestane", 10))),
    B2_exe = trial_arm(list(l1a, R_(NULL, "exemestane", NULL))),
    EV_ctrl = trial_arm(list(l1c, R_(NULL, "exemestane", 10))),
    AC_gir_c = trial_arm(list(l1c, R_(NULL, "giredestrant", NULL))),
    AC_ful_c = trial_arm(list(l1c, R_(NULL, "fulvestrant", NULL))),
    AC_gir_a = trial_arm(list(l1a, R_(NULL, "giredestrant", NULL))),
    AC_ful_a = trial_arm(list(l1a, R_(NULL, "fulvestrant", NULL))),
    ML2_os = list(lines = list(R_("ribo", "letrozole", NULL), R_(NULL, "exemestane", 10)))
  )
  if (schedules) {
    a$P2_cdk$ass1 <- SCHED$q12w
    a$P2_ai$ass1 <- SCHED$q12w
    for (k in c("B2_eve", "B2_exe")) a[[k]]$ass2 <- SCHED$q6w
    # evERA, acelERA (line 2) and MONALEESA-2 (line 1): every 8 weeks for 18 months, then every
    # 12 weeks = default schedule
  }
  a
}

validation_arms <- function(schedules = FALSE) {
  l1c <- R_("palbo", "letrozole", NULL)
  a <- list(
    EV_exp = trial_arm(list(l1c, R_(NULL, "giredestrant", 10))),
    EV_ctrl = trial_arm(list(l1c, R_(NULL, "exemestane", 10))),
    BOLERO4 = list(lines = list(R_(NULL, "letrozole", 10))),
    SER6_switch = list(lines = list(l1c), d3_switch = TRUE),
    SER6_cont = list(lines = list(l1c), d3_monitor = TRUE),
    pM_abema = trial_arm(list(l1c, R_("abema", "fulvestrant", NULL))),
    pM_ful = trial_arm(list(l1c, R_(NULL, "fulvestrant", NULL))),
    MA_ribo = trial_arm(list(l1c, R_("ribo", "fulvestrant", NULL))),
    TRINITI1 = trial_arm(list(l1c, R_("ribo", "exemestane", 2.5, c(300, 200))))
  )
  if (schedules) {
    a$BOLERO4$ass1 <- SCHED$q8w
    for (k in c("pM_abema", "pM_ful")) {
      a[[k]]$ass2 <- SCHED$q8w_12m
      a[[k]]$early_prog_days <- 0   # no minimum duration of prior CDK4/6i in postMONARCH
    }
    a$MA_ribo$ass2 <- SCHED$q12w
    a$MA_ful <- trial_arm(list(l1c, R_(NULL, "fulvestrant", NULL)), ass2 = SCHED$q12w)
    a$TRINITI1 <- trial_arm(list(l1c, R_("ribo", "exemestane", 2.5, c(300, 200), cdk_onoff = FALSE)),
                            ass2 = SCHED$q8w_12m)
    # SERENA-6 (extended analysis, Turner 2026): ctDNA every 2-3 months (here every 10 weeks) from 6
    # months on first-line AI + CDK4/6i; Guardant360 CDx (detection limit not reported): the PADA-1 effective limit
    # scaled by the ratio of the median allele fractions at a first detection after negative tests (about 0.26% vs
    # 0.8%); randomization about 5 weeks after the detection sample (results within 21 days, screening within 28
    # days) if imaging shows no progression; imaging every 8 weeks. Giredestrant is the SERD proxy for camizestrant.
    for (k in c("SER6_switch", "SER6_cont")) a[[k]] <- c(a[[k]], SER6_CT)
  }
  a
}
SER6_CT <- list(ct_first = 183L, ct_int = 70L, ct_lod_par = "LOD_pd", ct_lod_mult = 0.26 / 0.8, ct_lag = 35L, ct_confirm = TRUE,
                ass1 = c(56L, 56L, 100000L, 56L))

# persevERA (ASCO 2026, LBA1006): line 1 giredestrant or letrozole, each with palbociclib. The scan
# schedule is not in the abstract; the model default is used (every 8 weeks for 18 months, then 12).
persevera_arms <- function() list(
  PV_gir = list(lines = list(R_("palbo", "giredestrant", NULL))),
  PV_let = list(lines = list(R_("palbo", "letrozole", NULL)))
)

# PALOMA-3 (Cristofanilli, Lancet Oncol 2016; calibration): progression on or within 12 months of
# endocrine therapy (85% had received an AI) -> line 1 = letrozole alone; every patient who progresses on it starts
# line 2 = fulvestrant +/- palbociclib (no attrition, no minimum duration); imaging every 8 weeks for the first
# year, then every 12 weeks.
paloma3_arms <- function() list(
  P3_pf = list(lines = list(R_(NULL, "letrozole", NULL), R_("palbo", "fulvestrant", NULL)), p_attr = 0, p_attr_mid = 0,
               early_prog_days = 0, ass2 = SCHED$q8w_12m),
  P3_f = list(lines = list(R_(NULL, "letrozole", NULL), R_(NULL, "fulvestrant", NULL)), p_attr = 0, p_attr_mid = 0,
              early_prog_days = 0, ass2 = SCHED$q8w_12m)
)

# PADA-1 (Bidard, Lancet Oncol 2022; Cabel, Ann Oncol 2026): first-line AI + palbociclib (letrozole as the
# AI) with blood ESR1 testing (multiplex ddPCR) at inclusion, after one cycle and then every two cycles; imaging after
# three cycles and then every two cycles. A rising ESR1 mutation without progression on imaging led to randomization
# (about one month after the detection sample) to fulvestrant + palbociclib or to continuing AI + palbociclib. The
# assay's detection limit is not reported (lowest allele fraction at detection 0.1%, median 0.8%, and 47% of mutations
# not re-detected a month later): the model's deterministic threshold is an effective limit, tum.LOD_pd, calibrated
# on step 1 (ESR1 emergence). Step 2 (switch vs continue) is held out for validation.
PADA1_CT <- list(ct_first = 28L, ct_int = 56L, ct_lod_par = "LOD_pd", ct_lag = 30L, ct_confirm = TRUE,
                 ass1 = c(84L, 56L, NEVER, 56L))
pada1_arms <- function() list(
  PD_cont = c(list(lines = list(R_("palbo", "letrozole", NULL)), d3_monitor = TRUE), PADA1_CT),
  PD_switch = c(list(lines = list(R_("palbo", "letrozole", NULL)), d3_switch = TRUE, switch_et = "fulvestrant"), PADA1_CT)
)

# MIRACLE (Fan, JAMA Oncol 2021; validation): premenopausal patients after tamoxifen, first-line
# letrozole (with goserelin) +/- everolimus 10 mg; imaging every 8 weeks, every 12 weeks after the first year.
# The model has no menopausal status or tamoxifen history: letrozole with ovarian suppression is emulated as
# letrozole; no mouthwash prophylaxis (the trial predates its routine use).
miracle_arms <- function() list(
  MI_eve = list(lines = list(R_(NULL, "letrozole", 10)), ass1 = SCHED$q8w_12m, tox_over = NODEX),
  MI_let = list(lines = list(R_(NULL, "letrozole", NULL)), ass1 = SCHED$q8w_12m, tox_over = NODEX)
)

# First-line fulvestrant (validation; abstracts only, so the scan schedule is the model default):
# PARSIFAL (Llombart-Cussac, JAMA Oncol 2021): fulvestrant + palbociclib vs letrozole + palbociclib;
# FALCON (Robertson, Lancet 2016): fulvestrant vs anastrozole (letrozole as the non-steroidal AI proxy).
fulvestrant_1l_arms <- function() list(
  PF_ful = list(lines = list(R_("palbo", "fulvestrant", NULL))),
  PF_let = list(lines = list(R_("palbo", "letrozole", NULL))),
  FA_ful = list(lines = list(R_(NULL, "fulvestrant", NULL))),
  FA_ai = list(lines = list(R_(NULL, "letrozole", NULL)))
)

# SONIA (Sonke, Nature 2024): AI + CDK4/6i then fulvestrant (first), or AI then fulvestrant + CDK4/6i
# (second); imaging every 12 weeks; every progressing patient eligible for the protocol second line.
# No endocrine-based second line after stopping line 1: 60 of 315 (19%) and 51 of 404 (13%).
sonia_arms <- function(p_first = 0.19, p_second = 0.13) list(
  SO_first = list(lines = list(R_("palbo", "letrozole", NULL), R_(NULL, "fulvestrant", NULL)),
                  p_attr = p_first, p_attr_mid = p_first, early_prog_days = 0,
                  ass1 = SCHED$q12w, ass2 = SCHED$q12w),
  SO_second = list(lines = list(R_(NULL, "letrozole", NULL), R_("palbo", "fulvestrant", NULL)),
                   p_attr = p_second, p_attr_mid = p_second, early_prog_days = 0,
                   ass1 = SCHED$q12w, ass2 = SCHED$q12w)
)
