# Published targets for calibration and validation, toxicity-only emulations and scenario settings. Every value was
# checked against the primary publication or label; sources and page references: verification/SOURCE_VERIFICATION.md.

# ---------------------------------------------------------------------------- calibration targets
CAL_TARGETS <- list(
  # PALOMA-2 extended follow-up (Rugo, Breast Cancer Res Treat 2019;174:719-729; doi 10.1007/s10549-018-05125-4):
  # median follow-up about 38 months; median PFS 27.6 vs 14.5 months, HR 0.563. Response rates: primary analysis.
  P2_cdk_mPFS = 27.6, P2_ai_mPFS = 14.5, P2_HR = 0.563,
  P2_ORR_cdk = 0.553, P2_ORR_ai = 0.444,                  # measurable disease (NEJM 2016, Table 3)
  B2_eve_mPFS = 7.8, B2_exe_mPFS = 3.2, B2_HR = 0.45,     # BOLERO-2 final, local (Yardley 2013)
  EV_ctrl_mPFS_all = 5.5, EV_ctrl_mPFS_esr1 = 5.5,        # evERA control (NEJM 2026 p7)
  AC_gir_mPFS = 5.6, AC_ctl_mPFS = 5.4, AC_HR = 0.81, AC_HR_esr1 = 0.60,   # acelERA (JCO 2024 p6)
  ML2_OS_cdk = 63.9,                                      # MONALEESA-2 OS (NEJM 2022)
  EV_ctrl_OS18 = 62.0,                                    # evERA control 18-month OS (NEJM 2026 p8)
  # persevERA (ASCO 2026, LBA1006 p1): investigator PFS HR 0.89 (0.76-1.05); medians 33.1 and 28.2 months;
  # objective response 60.2% vs 58.8% (the difference is fitted, not the levels, which depend on the scan schedule)
  PV_HR = 0.89, PV_gir_mPFS = 33.1, PV_let_mPFS = 28.2, PV_ORR_diff = 0.014,
  # persevERA landmark PFS (giredestrant vs letrozole, each with palbociclib), percent at 12, 24 and 36 months:
  # reports of the ASCO presentation (Seymour 2026, OncLive; Targeted Oncology Staff 2026), identical values
  PV_landmarks = list(m = c(12, 24, 36), gir = c(77.8, 59.7, 45.8), let = c(77.3, 57.3, 41.9)),
  # PALOMA-3 final PFS analysis (Cristofanilli, Lancet Oncol 2016; pp. 425, 433; Table 2): investigator PFS 9.5
  # (9.2-11.0) vs 4.6 (3.5-5.6) months, HR 0.46 (0.36-0.59); objective response in measurable disease 24.6% vs 10.9%
  P3_pf_mPFS = 9.5, P3_f_mPFS = 4.6, P3_HR = 0.46, P3_ORR = c(pf = 0.246, f = 0.109),
  # PADA-1 step 1 (first-line AI + palbociclib, serial ctDNA ESR1 testing): baseline prevalence 33/1017 = 3.2%
  # (Cabel 2026 p3); rising mutation in 283 of the 678 patients with progression or a rising mutation before
  # monitoring stopped = 41.7% (Cabel p7); synchronous progression in 71/283 = 25% (Cabel p3 Fig 1; Bidard 2022:
  # 69/279); median 15.6 months on AI + palbociclib before randomization (Bidard p7); median PFS 26.8 months (Bidard p7)
  PD1 = list(base_prev = 0.032, rise_share = 0.417, sync_share = 0.25, t_rand = 15.6, mPFS = 26.8)
)

TOX_TARGETS <- list(
  neut_g3_palbo = 0.664,   # PALOMA-2 grade 3-4 neutropenia 66.4%
  neut_g3_ribo = 0.593,    # MONALEESA-2 first report 59.3%
  neut_g3_abema = 0.211,   # MONARCH 3 interim (primary) 21.1%
  disc_cdk = 0.074,        # PALOMA-2 palbociclib permanently stopped for AEs 7.4%
  stom_any_eve = 0.489, stom_g3_eve = 0.032, pneu_any_eve = 0.124, pneu_g3_eve = 0.027,
  hyperg_g3_eve = 0.022, diar_g3_eve = 0.011, anemia_g3_eve = 0.086,
  disc_eve = 0.118, disc_et = 0.065,                       # evERA control arm (Table 2, p10)
  diar_g2_abema = 0.367,   # MONARCH 3 grade >= 2 diarrhea 27.2% + 9.5%
  diar_g3_abema = 0.095,   # MONARCH 3 grade 3 diarrhea
  alt_g3_ribo = 0.093, qtc_ribo = 0.033,                   # MONALEESA-2 first report
  brady_gir = 0.038, disc_gir = 0.082                      # evERA experimental arm
)

# Median treatment durations (days) for the toxicity-only calibration runs
TOX_DURATIONS <- c(
  palbo = round(19.8 * MONTH),   # palbociclib label: PALOMA-2 palbociclib + letrozole 19.8 months
  ribo = round(13.0 * MONTH),    # MONALEESA-2 first report: median exposure 13.0 months
  abema = round(15.1 * MONTH),   # abemaciclib label: MONARCH 3 median treatment 15.1 months
  evectl = round(4.6 * MONTH),   # evERA control: standard endocrine therapy 4.6 months
  evexp = round(7.1 * MONTH)     # evERA experimental: giredestrant 7.1 months
)

# ---------------------------------------------------------------------------- validation targets (held out)
VAL_TARGETS <- list(
  EV = list(esr1 = c(10.0, 8.1, 12.9, 0.38, 0.27, 0.54), all = c(8.8, 6.6, 9.6, 0.56, 0.44, 0.71),
            noesr1 = c(5.7, NA, NA, 0.84, 0.59, 1.18)),
  EV_landmarks = list(c("esr1", 6, 66.1, 38.1), c("esr1", 12, 40.5, 15.2), c("all", 6, 57.4, 39.6),
                      c("all", 12, 34.1, 18.1)),
  EV_ORR = list(esr1 = c(0.266, 0.138), all = c(0.238, 0.117)),
  EV_OS18 = list(all = c(75.1, 62.0), esr1 = c(71.5, 50.9)),
  EV_tox_exp = list(stom_any = 0.473, stom_g3 = 0.027, pneu_any = 0.170, pneu_g3 = 0.027, ae_disc_eve = 0.170),
  BOLERO4 = c(22.0, 18.1, 25.1),
  postMONARCH = list(HR = c(0.73, 0.57, 0.95), mPFS = c(6.0, 5.3)),
  MAINTAIN = c(0.57, 0.39, 0.85),
  TRINITI1 = c(5.7, 3.6, 9.1),
  # SONIA (external validation) and its ctDNA analysis
  SONIA = list(PFS1_HR = c(0.59, 0.51, 0.69), PFS1 = c(24.7, 16.1), PFS2_HR = c(0.87, 0.74, 1.03),
               PFS2 = c(31.0, 26.8), OS_HR = c(0.98, 0.80, 1.20), CDK_months = c(24.6, 8.1),
               ct_high_HR = c(0.58, 0.38, 0.88), ct_low_HR = c(1.36, 0.95, 1.96), ct_high_frac = 0.34),
  PV = list(OS_HR = c(1.03, 0.83, 1.28), ET_disc = c(gir = 0.065, let = 0.055)),
  # First-line fulvestrant, abstracts (PubMed): PARSIFAL (Llombart-Cussac, JAMA Oncol 2021;7:1791-1799,
  # doi 10.1001/jamaoncol.2021.4301): fulvestrant + palbociclib vs letrozole + palbociclib, median PFS 27.9 vs 32.8
  # months, HR 1.13 (0.89-1.45), response 46.5% vs 50.2%; FALCON (Robertson, Lancet 2016;388:2997-3005,
  # doi 10.1016/S0140-6736(16)32389-3): fulvestrant vs anastrozole in endocrine-therapy-naive patients, 16.6 vs 13.8
  # months, HR 0.797 (0.637-0.999).
  PARSIFAL = list(HR = c(1.13, 0.89, 1.45), mPFS = c(27.9, 32.8), ORR = c(0.465, 0.502)),
  FALCON = list(HR = c(0.797, 0.637, 0.999), mPFS = c(16.6, 13.8)),
  # MIRACLE (Fan, JAMA Oncol 2021; pp. 1, 4, 5; Supplement 2 eFigure 2 and eTable): everolimus + letrozole vs
  # letrozole (premenopausal, with goserelin, after tamoxifen): PFS HR 0.64 (0.46-0.89); medians 19.4 vs 12.9 months
  # (main text; eFigure 2 gives 16.5 vs 10.0 with the same HR, so the HR is the validation criterion); response in
  # measurable disease 50.0% vs 39.3%; adverse events leading to discontinuation 8.9% (everolimus arm) vs 0%;
  # stomatitis any grade 79.2%, grade >= 3 9.9% (no mouthwash prophylaxis).
  MIRACLE = list(HR = c(0.64, 0.46, 0.89), mPFS = c(19.4, 12.9), ORR = c(0.500, 0.393), ae_disc = 0.089,
                 stom = c(any = 0.792, g3 = 0.099)),
  # PADA-1 step 2: switch to fulvestrant + palbociclib vs continue AI + palbociclib on a rising ESR1 mutation
  # without progression: PFS from randomization 11.9 (9.1-13.6) vs 5.7 (3.9-7.5) months, stratified HR 0.61 (0.43-0.86)
  # (Bidard p7); update 12.8 vs 5.9 months, HR 0.54 (0.38-0.75) (Cabel p7); median allele fraction at the rising
  # detection 0.8% (Cabel p8)
  PADA1 = list(switch = c(11.9, 9.1, 13.6), cont = c(5.7, 3.9, 7.5), HR = c(0.61, 0.43, 0.86),
               upd = c(12.8, 5.9, 0.54, 0.38, 0.75), vaf_det = 0.008),
  # SERENA-6 extended analysis (Turner 2026, pp. 8-10): PFS 16.8 (14.7-19.4) vs 9.2
  # (7.2-9.7) months, HR 0.45 (0.34-0.59); 24-month PFS 34.9% vs 14.2%; PFS2 HR 0.63 (0.46-0.86); OS HR 0.87
  # (0.57-1.30, immature); start of first-line therapy to detection, median 22 months (randomized); ESR1 mutation with
  # concurrent progression 117/548 (21%); median allele fraction at a subsequent-test detection about 0.26%
  # (Fig 1C, read from the figure)
  SERENA6x = list(switch = 16.8, cont = 9.2, HR = c(0.45, 0.34, 0.59), PFS24 = c(34.9, 14.2),
                  OS_HR = c(0.87, 0.57, 1.30), t_det = 22, conc_pd = 0.21, vaf_det = 0.0026)
)

# ---------------------------------------------------------------------------- scenarios
# Dexamethasone mouthwash. Base case = evERA practice (mouthwash strongly recommended).
NODEX <- list(dex_mult = 1.63, dex_sev = 1.86)       # BOLERO-2 (no prophylaxis) vs evERA control
DEX_STRICT <- list(dex_mult = 0.47, dex_sev = 0.62)  # SWISH (strict prophylaxis) vs evERA control

# Uncertainty ranges for the attrition parameters
ATTR_PRIOR <- list(p_attr = c(0.20, 0.40), p_attr_mid = c(0.45, 0.70))   # parameter ensemble
ATTR_GSA <- list(p_attr = c(0.10, 0.50), p_attr_mid = c(0.30, 0.80))     # global sensitivity

# ---------------------------------------------------------------------------- AE discontinuation
# Toxicity-only emulations: exposure exponential around the published median exposure (months), capped at
# the longest follow-up of the analysis that reported the rate (months).
DISC_ARMS <- list(
  P2_cdk = list(cdk = "palbo", et = "letrozole", eve = NULL, med = 19.8, cap = 36),   # PALOMA-2 (palbociclib label: 19.8 mo)
  P2_ai = list(cdk = NULL, et = "letrozole", eve = NULL, med = 13.8, cap = 36),       # PALOMA-2 placebo arm (13.8 mo)
  PV_let = list(cdk = "palbo", et = "letrozole", eve = NULL, med = 25.2, cap = 60),   # persevERA: median treatment duration
  PV_gir = list(cdk = "palbo", et = "giredestrant", eve = NULL, med = 29.8, cap = 60),#   from the presentation (two reports)
  B4 = list(cdk = NULL, et = "letrozole", eve = 10, med = 75.1 * 7 / MONTH, cap = 42,  # BOLERO-4 (letrozole 75.1 weeks;
            over = NODEX),                                                             #   no mouthwash prophylaxis then)
  EV_ctl = list(cdk = NULL, et = "exemestane", eve = 10, med = 4.6, cap = 24),        # evERA control (4.6 mo)
  EV_exp = list(cdk = NULL, et = "giredestrant", eve = 10, med = 7.1, cap = 24),      # evERA experimental (7.1 mo)
  M3 = list(cdk = "abema", et = "letrozole", eve = NULL, med = 15.1, cap = 27)        # MONARCH 3 (abemaciclib label: 15.1 mo)
)
DISC_TARGETS <- data.frame(
  key = c("P2_cdk_cdk", "P2_cdk_any", "P2_ai_any", "PV_let_et", "PV_gir_et", "B4_any", "EV_ctl_eve", "EV_ctl_et",
          "EV_exp_et", "M3_cdk"),
  arm = c("P2_cdk", "P2_cdk", "P2_ai", "PV_let", "PV_gir", "B4", "EV_ctl", "EV_ctl", "EV_exp", "M3"),
  outcome = c("cdk", "any", "any", "et", "et", "any", "eve", "et", "et", "cdk"),
  target = c(0.074, 0.097, 0.059, 0.055, 0.065, 0.252, 0.118, 0.065, 0.082, 0.196),
  n = c(444, 444, 222, 497, 495, 202, 186, 186, 182, 327),
  label = c("PALOMA-2: palbociclib stopped", "PALOMA-2 palbociclib + letrozole: any study drug stopped",
            "PALOMA-2 letrozole + placebo: any study drug stopped", "persevERA: letrozole stopped (with palbociclib)",
            "persevERA: giredestrant stopped (with palbociclib)", "BOLERO-4 everolimus + letrozole: any drug stopped",
            "evERA control: everolimus stopped", "evERA control: endocrine therapy stopped",
            "evERA experimental: giredestrant stopped", "MONARCH 3: abemaciclib stopped"),
  source = c("NEJMoa1607303 p7", "NEJMoa1607303 p1, p7", "NEJMoa1607303 p1, p7", "persevERA_Turner2026_LBA1006 p1",
             "persevERA_Turner2026_LBA1006 p1", "jamaoncology_royce_2018 p4", "NEJMoa2602457 p10", "NEJMoa2602457 p10",
             "NEJMoa2602457 p10", "goetz2017 p3"),
  stringsAsFactors = FALSE)
# primary target of each hazard level (fitted by bisection for a given decay shape)
DISC_LEVELS <- data.frame(par = c("disc_cdk", "disc_et", "disc_gir", "disc_eve", "disc_et_eve", "disc_abema_mult"),
                          key = c("P2_cdk_cdk", "P2_ai_any", "PV_gir_et", "EV_ctl_eve", "EV_ctl_et", "M3_cdk"),
                          stringsAsFactors = FALSE)

# Trial-like attrition (SONIA: no endocrine-based second line in 19% and 13% of the patients who stopped line 1)
ATTR_TRIAL <- c(base = 0.15, low = 0.10, high = 0.20)
