# Model parameters and dose tables.
#
# Units: time in days, drug exposure as relative concentration (1 = typical patient at reference
# dose at steady state), estradiol (E2) relative to the untreated postmenopausal level.
# Values taken from publications and labels are listed with their sources in verification/SOURCE_VERIFICATION.md.
# The calibration steps write the fitted values to results/params_*.json, and load_params() (R/runner.R)
# applies them on top of default_params().

SIGNALING <- list(
  Kd = 0.10,      # E2 level giving half-maximal activation of wild-type ER
  pm = 0.7,       # share of ER transcriptional output from the mutant allele (dominant)
  fc = 0.90,      # ligand-independent (constitutive) activity of mutant ER
  s_mut = 3.0,    # SERD potency shift against mutant ER (higher = less potent)
  dS = 4.0,       # SERD-induced ER degradation at full occupancy (turnover x (1 + dS))
  aF = 1.0,       # FOXO-driven ESR1 synthesis gain
  KF = 1.0, nF = 2.0,    # FOXO inhibition by AKT
  bS = 0.30,      # S6K (Ser167) contribution to ER transcriptional activity
  cER = 0.50,     # ER-driven IRS1/IGF1R expression
  Kfb = 1.0, hfb = 2.0,  # S6K -> IRS-1 negative feedback
  gP = 1.8,       # PI3K pathway gain in PIK3CA/AKT1/PTEN-altered clones
  theta4E = 0.40, # everolimus effect on 4E-BP1 axis relative to S6K (partial)
  KT = 0.5,       # saturation of cap-dependent translation
  e0 = 0.5,       # E2F-independent cyclin E fraction
  gE_R = 4.0,     # cyclin E / CDK2 gain in CDK4/6i-escape (R) clones
  cE = 0.7,       # ER-independent share of cyclin E/CDK2 activity
  tE = 0.5,       # mTORC1 (translation) dependence of cyclin E/CDK2 activity
  wD = 0.7,       # weight of CDK4/6 vs CDK2 in Rb phosphorylation
  KRb = 0.6, nRb = 3.0,
  eps_g = 0.50    # dependence of cell growth on S6K (ER-independent everolimus effect)
)

TUMOR <- list(
  kmax_med = 0.030,    # median max proliferation rate (1/day) at P = 1
  kmax_sd = 0.40,      # between-patient SD of log(kmax)
  rd = 0.80,           # basal death as fraction of untreated proliferation
  rd_sd = 0.04,
  b0_mean = 0.22,      # ER-independent share of cyclin D1 drive (endocrine independence)
  b0_conc = 8.0,       # beta concentration (higher = less spread)
  pP = 0.38,           # prevalence of truncal PIK3CA/AKT1/PTEN alteration (first line)
  log10_fE0 = -4.0,    # median initial ESR1-mutant subclone fraction
  sd_log10_fE0 = 1.0,
  pE_base = 0.03,      # probability of clinically detectable baseline ESR1 subclone
  log10_fR0 = -5.0,    # median initial CDK4/6i-escape subclone fraction
  sd_log10_fR0 = 1.0,
  mu_E = 1e-9, mu_P = 1e-9, mu_R = 1e-8,   # mutation probability per division
  sigma_R = 2.0e-4,    # adaptive switching to CDK4/6i escape per day at full CDK4/6 inhibition
  kappa = 0.5,         # mTORC1 dependence of adaptive switching (0 = none, 1 = full)
  N0_med = 2e10, N0_sd = 0.8,   # baseline tumor cells (lognormal)
  SLD_ref = 60.0,      # mm, SLD at N0_med
  meas_cv = 0.05,      # SLD measurement error
  lam_NL = 1.5e-4,     # new-lesion hazard (1/day) at nadir burden
  gam_NL = 2.0,
  K_TF = 2e11,         # tumor fraction TF = N/(N + K_TF)
  LOD = 0.0025,        # ctDNA variant allele fraction limit of detection
  h_t = 4.0e-4,        # tumor-related death hazard (1/day) at N = N_ref
  N_ref = 2e11,
  gam_d = 1.5,
  h_post = 0.0,        # extra death hazard (1/day) after endocrine-based therapy has failed
  gomp_a = -10.5, gomp_b = 0.09,   # background mortality: exp(a + b*age) per year
  dc0_med = 0.020, dc0_sd = 0.4,   # post-endocrine therapy (chemo/ADC): initial extra kill rate (1/day)
  tau_c = 150.0,       # decay of post-endocrine effect (days)
  # attrition (no endocrine-based 2nd line after 1st-line progression at >= 6 months): ESME-MBC (Varnier 2025,
  # Suppl. Table A.2), 1st-line PFS > 12 months 75/257 = 29%, 6-12 months 75/131 = 57%
  p_attr = 0.29,       # progression after more than 12 months on line 1
  early_prog_days = 183,  # progression before this on line 1 -> chemotherapy
  p_attr_mid = 0.57,   # progression after 6 to 12 months on line 1
  attr_switch_days = 365,
  # effective detection limit (variant allele fraction) of the PADA-1 ctDNA assay, calibrated with the ESR1 emergence
  # targets (starting value 0.2%); SERENA-6 uses it scaled by 0.26 / 0.8 (R/arms.R)
  LOD_pd = 0.002
)

PHARM <- list(
  Imax_c = 0.95, IC50_c = 0.35, h_c = 1.5,   # CDK4/6 inhibition by relative exposure
  I10_e = 0.70, IC50_e = 0.25,   # everolimus: inhibition at 10 mg exposure, dose-response IC50
  uG_ref = 400.0,      # SERD potency index, giredestrant 30 mg
  uF_ref = 40.0,       # fulvestrant 500 mg
  E2_let = 0.02, E2_exe = 0.03,  # E2 under aromatase inhibitors
  k_palbo = 0.57, k_ribo = 0.52, k_abema = 0.90, k_eve = 0.55, k_gir = 0.40, k_fulv = 0.0173,
  sd_palbo = 0.30, sd_ribo = 0.35, sd_abema = 0.40, sd_eve = 0.40, sd_gir = 0.30, sd_fulv = 0.30,
  # CYP3A4 time-dependent inhibition: kinact_ribo and kinact_palbo are set in default_params() from the
  # midazolam AUC ratios at the reference dose (ribociclib label 12.3: 5.2-fold predicted at 600 mg;
  # palbociclib label 12.3: +61%)
  kdeg_3A4 = 0.5, KI_3A4 = 0.5, kinact_ribo = NA_real_, kinact_palbo = NA_real_,
  fm_eve = 1 - 1 / 15, # ketoconazole raises everolimus AUC 15-fold (everolimus label 12.3), so fm >= 1 - 1/15
  gir_ddi = 1.0        # giredestrant effect on everolimus exposure (unknown, 1 = none)
)

TOX <- list(
  MTT = 5.0, gam_N = 0.17, ANC0_med = 4.0, ANC0_sd = 0.30,          # neutropenia (Friberg)
  slope_palbo = 1.0, slope_ribo = 0.9, slope_abema = 0.35, slope_eve = 0.05,
  aS_eve = 0.010, aS_cdk = 0.0010, tauS = 45.0, p2S = 0.20, p3S = 0.065, nodex_mult = 1.6,  # stomatitis
  aP_eve = 0.0006, aP_cdk = 0.00003,                                  # pneumonitis
  aG_eve = 0.00015,                                                   # hyperglycemia grade >= 3
  aD_abema = 0.006, aD_eve = 0.0004, aD_gir = 0.0001, aD_palbo = 0.0002, aD_ribo = 0.0002,
  tauD_abema = 60.0, q3D_abema = 0.30, q3D_other = 0.20,              # diarrhea grade >= 2
  aL_ribo = 0.0006, aL_abema = 0.0003, aL_palbo = 0.0001, aL_eve = 0.00005, aL_gir = 0.00008,
  aQ_ribo = 0.0012, tauQ = 30.0,                                      # QTc > 480 ms (ribociclib)
  aB_gir = 0.0007, tauB = 60.0, aA_eve = 0.0004,                      # counted only
  disc_cdk = 0.00025, disc_eve = 0.0006, disc_et = 0.0003, disc_gir = 0.0004,
  tox_int = 1.0,                                                      # interaction multiplier (>= 3 drugs)
  # discontinuation for other AEs: the hazard falls from the first day of each line (shape fitted in
  # scripts/02_calibrate_tox.R), everolimus discontinuation is proportional to everolimus exposure, and endocrine
  # therapy has an extra hazard while everolimus is given
  disc_t0 = 0,
  disc_eve_gamma = 1,  # everolimus discontinuation proportional to exposure (assumption; 0 in a sensitivity analysis)
  disc_tau = 120,      # starting values, fitted in 02_calibrate_tox.R
  disc_floor = 0.25,
  disc_et_eve = 3e-4,
  disc_abema_mult = 1  # abemaciclib multiplier of the CDK4/6i hazard (MONARCH 3)
)

STRUCT <- list(s6k_er = TRUE, adaptive_R = TRUE, reversible_R = FALSE)

default_params <- function() {
  p <- list(sig = SIGNALING, tum = TUMOR, ph = PHARM, tox = TOX, st = STRUCT)
  p$ph$kinact_ribo <- kinact_for_aucr(p$ph, 5.2, 1.0)
  p$ph$kinact_palbo <- kinact_for_aucr(p$ph, 1.61, 1.0)
  p
}

# ---------------------------------------------------------------------------- doses
CDK_DOSES <- list(palbo = c(125, 100, 75), ribo = c(600, 400, 200), abema = c(150, 100, 50))
CDK_REF <- c(palbo = 125, ribo = 600, abema = 150)
CDK_ONOFF <- c(palbo = TRUE, ribo = TRUE, abema = FALSE)   # 21 days on / 7 off
CDK_CODE <- c(palbo = 1L, ribo = 2L, abema = 3L)
ET_CODE <- c(letrozole = 1L, exemestane = 2L, fulvestrant = 3L, giredestrant = 4L)
EVE_DOSES <- list("10" = c(10, 5, 2.5), "5" = c(5, 2.5), "2.5" = c(2.5))
