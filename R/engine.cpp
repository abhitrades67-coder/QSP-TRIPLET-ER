// QSP-TRIPLET-ER population engine: signalling grid and daily simulation of virtual patients.
//
// Called from R (R/engine.R) through Rcpp. One daily time step per virtual patient. Patients are
// independent, so they are simulated in parallel (OpenMP). Sums use a fixed pairwise summation order.
//
// Random numbers. Each (patient, day, stream) gets its own uniform from a counter-based
// generator (SplitMix64), so the same patient sees the same random stream under every arm
// (common random numbers) and results do not depend on the number of threads. The engine can instead
// read externally supplied random streams (rng_mode = 1; for testing).
//
// Genotype code g = E + 2P + 4R (E: ESR1 mutant, P: PI3K pathway altered, R: CDK4/6i escape).

// [[Rcpp::plugins(openmp)]]
#include <Rcpp.h>
#include <cmath>
#include <cstdint>
#include <vector>
#include <algorithm>
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace Rcpp;

// ============================================================================ helpers
static inline double num(const List& L, const char* nm) {
  if (!L.containsElementNamed(nm)) Rcpp::stop(std::string("missing parameter: ") + nm);
  return as<double>(L[nm]);
}
static inline double num_or(const List& L, const char* nm, double dflt) {
  if (!L.containsElementNamed(nm)) return dflt;
  SEXP x = L[nm];
  if (Rf_isNull(x)) return dflt;
  double v = as<double>(x);
  return ISNAN(v) ? dflt : v;
}
static inline bool flag_or(const List& L, const char* nm, bool dflt) {
  if (!L.containsElementNamed(nm)) return dflt;
  SEXP x = L[nm];
  if (Rf_isNull(x)) return dflt;
  return as<bool>(x);
}

// Pairwise summation along a contiguous axis: below 8 elements a plain loop, from 8 elements 8 partial sums.
static inline double np_sum4(double a0, double a1, double a2, double a3) {
  double r = 0.;
  r += a0; r += a1; r += a2; r += a3;
  return r;
}
static inline double np_sum8(const double* a) {
  return ((a[0] + a[1]) + (a[2] + a[3])) + ((a[4] + a[5]) + (a[6] + a[7]));
}
static inline double np_sum16(const double* a) {
  double r0 = a[0] + a[8], r1 = a[1] + a[9], r2 = a[2] + a[10], r3 = a[3] + a[11];
  double r4 = a[4] + a[12], r5 = a[5] + a[13], r6 = a[6] + a[14], r7 = a[7] + a[15];
  return ((r0 + r1) + (r2 + r3)) + ((r4 + r5) + (r6 + r7));
}

// counter-based uniform random numbers (SplitMix64 output function)
static inline uint64_t mix64(uint64_t z) {
  z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
  z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
  return z ^ (z >> 31);
}
static inline uint64_t patient_key(uint64_t seed, uint64_t pid) {
  uint64_t h = mix64(seed * 0x9E3779B97F4A7C15ULL + 0x632BE59BD9B4E019ULL);
  return mix64(h + (pid + 1) * 0xD1B54A32D192ED03ULL);
}
static inline double unif_ctr(uint64_t key, uint64_t t, uint64_t k) {
  uint64_t x = mix64(key + (t * 32ULL + k + 1ULL) * 0x9E3779B97F4A7C15ULL);
  return ((double)(x >> 11) + 0.5) * (1.0 / 9007199254740992.0);   // (0, 1)
}

// ============================================================================ signaling
struct Sig {
  double Kd, pm, fc, s_mut, dS, aF, KF, nF, bS, cER, Kfb, hfb, gP, theta4E, KT, e0, gE_R, cE,
         tE, wD, KRb, nRb, eps_g;
};

static Sig make_sig(const List& s, const List& st) {
  Sig q;
  q.Kd = num(s, "Kd"); q.pm = num(s, "pm"); q.fc = num(s, "fc"); q.s_mut = num(s, "s_mut");
  q.dS = num(s, "dS"); q.aF = num(s, "aF"); q.KF = num(s, "KF"); q.nF = num(s, "nF");
  q.bS = flag_or(st, "s6k_er", true) ? num(s, "bS") : 0.0;
  q.cER = num(s, "cER"); q.Kfb = num(s, "Kfb"); q.hfb = num(s, "hfb"); q.gP = num(s, "gP");
  q.theta4E = num(s, "theta4E"); q.KT = num(s, "KT"); q.e0 = num(s, "e0"); q.gE_R = num(s, "gE_R");
  q.cE = num(s, "cE"); q.tE = num(s, "tE"); q.wD = num(s, "wD"); q.KRb = num(s, "KRb");
  q.nRb = num(s, "nRb"); q.eps_g = num(s, "eps_g");
  return q;
}

// Quasi-steady targets of the six signaling states.
// y = (ER, IRS, AKT, CycD, CycE, pRb); tg = targets; aux = (S6K, Tr, A_ER, CDK46, occ)
static inline void sig_targets(const double* y, double E2, double u, double Ie, double Ic, double b0,
                               int isE, int isP, int isR, const Sig& s, double* tg, double* aux) {
  const double ER = y[0], IRS = y[1], AKT = y[2], CycD = y[3], CycE = y[4], pRb = y[5];
  double Kd = s.Kd;
  double fE = E2 / (E2 + Kd);
  double fE_ref = 1.0 / (1.0 + Kd);
  double occW = u / (u + 1.0 + E2 / Kd);
  double occM = u / (u + s.s_mut * (1.0 + E2 / Kd));
  double pmE = isE ? s.pm : 0.0;
  double occ = (1 - pmE) * occW + pmE * occM;
  double L = (1 - pmE) * (1 - occW) * fE + pmE * (1 - occM) * (s.fc + (1 - s.fc) * fE);
  double bS = s.bS;
  double S6K = AKT * (1 - Ie);
  double T4 = AKT * (1 - s.theta4E * Ie);
  double Tr = T4 * (1 + s.KT) / (T4 + s.KT);
  double FOXO = 1.0 / (1.0 + std::pow(AKT / s.KF, s.nF));
  double FOXO0 = 1.0 / (1.0 + std::pow(1.0 / s.KF, s.nF));
  double ER_t = (1 + s.aF * FOXO) / (1 + s.aF * FOXO0) / (1 + s.dS * occ);
  double A_ER = ER * (L / fE_ref) * (1 + bS * S6K) / (1 + bS);
  double IRS_t = (1 + s.cER * A_ER) / (1 + s.cER) * (1 + std::pow(1 / s.Kfb, s.hfb)) /
                 (1 + std::pow(S6K / s.Kfb, s.hfb));
  double gPg = isP ? s.gP : 1.0;
  double AKT_t = gPg * IRS;
  double CycD_t = (b0 + (1 - b0) * A_ER) * Tr;
  double CDK46 = CycD * (1 - Ic);
  double pRb0 = 1.0 / (std::pow(s.KRb, s.nRb) + 1.0);
  double gEg = isR ? s.gE_R : 1.0;
  double CycE_t = gEg * (s.e0 + (1 - s.e0) * pRb / pRb0) * (s.cE + (1 - s.cE) * std::min(A_ER, 1.5)) *
                  (1 - s.tE * (1 - Tr));
  double X = s.wD * CDK46 + (1 - s.wD) * CycE;
  double pRb_t = std::pow(X, s.nRb) / (std::pow(s.KRb, s.nRb) + std::pow(X, s.nRb));
  tg[0] = ER_t; tg[1] = IRS_t; tg[2] = AKT_t; tg[3] = CycD_t; tg[4] = CycE_t; tg[5] = pRb_t;
  aux[0] = S6K; aux[1] = Tr; aux[2] = A_ER; aux[3] = CDK46; aux[4] = occ;
}

// Vectorized damped fixed-point iteration. The stopping rule is global (all points iterate until the
// largest change over all points is below tol).
// [[Rcpp::export]]
List steady_state_cpp(NumericVector E2, NumericVector u, NumericVector Ie, NumericVector Ic,
                      NumericVector b0, IntegerVector g, List sig, List st, int n_iter = 400,
                      double damp = 0.5, double tol = 1e-8, int nthreads = 1) {
  const int n = E2.size();
  if (u.size() != n || Ie.size() != n || Ic.size() != n || b0.size() != n || g.size() != n)
    stop("steady_state_cpp: all inputs must have the same length");
  Sig s = make_sig(sig, st);
  double pRb0 = 1.0 / (std::pow(s.KRb, s.nRb) + 1.0);
  std::vector<double> Y((size_t)n * 6);
  for (int i = 0; i < n; i++) {
    double* y = &Y[(size_t)i * 6];
    y[0] = 1; y[1] = 1; y[2] = 1; y[3] = 1; y[4] = 1; y[5] = pRb0;
  }
  const double *pE2 = E2.begin(), *pu = u.begin(), *pIe = Ie.begin(), *pIc = Ic.begin(), *pb0 = b0.begin();
  const int* pg = g.begin();
  double diff = 0.0;
  int it_done = 0;
  for (int it = 0; it < n_iter; it++) {
    double dmax = 0.0;
#pragma omp parallel for num_threads(nthreads) reduction(max:dmax) schedule(static)
    for (int i = 0; i < n; i++) {
      double* y = &Y[(size_t)i * 6];
      double tg[6], aux[5];
      int gg = pg[i];
      sig_targets(y, pE2[i], pu[i], pIe[i], pIc[i], pb0[i], gg & 1, (gg >> 1) & 1, (gg >> 2) & 1, s, tg, aux);
      double dl = 0.0;
      for (int j = 0; j < 6; j++) {
        double d = std::fabs(tg[j] - y[j]);
        if (d > dl || std::isnan(d)) dl = d;
      }
      if (dl > dmax || std::isnan(dl)) dmax = dl;
      for (int j = 0; j < 6; j++) y[j] = damp * y[j] + (1 - damp) * tg[j];
    }
    diff = dmax;
    it_done = it + 1;
    if (diff < tol) break;
  }
  NumericVector ER(n), IRS(n), AKT(n), CycD(n), CycE(n), pRb(n), S6K(n), A_ER(n), Tr(n), CDK46(n),
      occ(n), P(n);
  for (int i = 0; i < n; i++) {
    double* y = &Y[(size_t)i * 6];
    double tg[6], aux[5];
    int gg = pg[i];
    sig_targets(y, pE2[i], pu[i], pIe[i], pIc[i], pb0[i], gg & 1, (gg >> 1) & 1, (gg >> 2) & 1, s, tg, aux);
    ER[i] = y[0]; IRS[i] = y[1]; AKT[i] = y[2]; CycD[i] = y[3]; CycE[i] = y[4]; pRb[i] = y[5];
    S6K[i] = aux[0]; Tr[i] = aux[1]; A_ER[i] = aux[2]; CDK46[i] = aux[3]; occ[i] = aux[4];
    double pp = y[5] * (1 - s.eps_g * std::max(0.0, 1.0 - aux[0]));
    P[i] = pp < 0.0 ? 0.0 : pp;   // np.clip(., 0, None)
  }
  return List::create(_["ER"] = ER, _["IRS"] = IRS, _["AKT"] = AKT, _["CycD"] = CycD, _["CycE"] = CycE,
                      _["pRb"] = pRb, _["S6K"] = S6K, _["A_ER"] = A_ER, _["Tr"] = Tr, _["CDK46"] = CDK46,
                      _["occ"] = occ, _["P"] = P, _["n_iter"] = it_done, _["max_diff"] = diff);
}

// Right-hand side of the signaling ODE (hours), for the R ODE tier (deSolve).
// [[Rcpp::export]]
NumericVector sig_rhs_cpp(NumericVector y, double E2, double u, double Ie, double Ic, double b0, int g,
                          List sig, List st, NumericVector k) {
  Sig s = make_sig(sig, st);
  double tg[6], aux[5];
  sig_targets(y.begin(), E2, u, Ie, Ic, b0, g & 1, (g >> 1) & 1, (g >> 2) & 1, s, tg, aux);
  NumericVector out(6);
  for (int j = 0; j < 6; j++) out[j] = k[j] * (tg[j] - y[j]);
  return out;
}

// Auxiliary signaling outputs (S6K, Tr, A_ER, CDK46, occ, P) for given states.
// [[Rcpp::export]]
NumericVector sig_aux_cpp(NumericVector y, double E2, double u, double Ie, double Ic, double b0, int g,
                          List sig, List st) {
  Sig s = make_sig(sig, st);
  double tg[6], aux[5];
  sig_targets(y.begin(), E2, u, Ie, Ic, b0, g & 1, (g >> 1) & 1, (g >> 2) & 1, s, tg, aux);
  double pp = y[5] * (1 - s.eps_g * std::max(0.0, 1.0 - aux[0]));
  return NumericVector::create(_["S6K"] = aux[0], _["Tr"] = aux[1], _["A_ER"] = aux[2],
                               _["CDK46"] = aux[3], _["occ"] = aux[4], _["P"] = pp < 0 ? 0 : pp);
}

// ============================================================================ grid
static const double U_GRID[9] = {0, 1, 3, 10, 30, 100, 300, 1000, 3000};
static const double IE_GRID[11] = {0, .1, .2, .3, .4, .5, .6, .7, .8, .9, .98};
static const double IC_GRID[11] = {0, .1, .2, .3, .4, .5, .6, .7, .8, .9, .95};
static const double B0_GRID[10] = {0, .05, .1, .15, .2, .3, .4, .55, .7, .9};
static const int sB = 1, sC = 10, sE = 110, sU = 1210, sG = 10890, sE2 = 87120;

static inline int locate(const double* grid, int len, double x, double& w) {
  int cnt = 0;                                   // np.searchsorted(grid, x, side='right')
  while (cnt < len && grid[cnt] <= x) cnt++;
  int i = cnt - 1;
  if (i < 0) i = 0;
  if (i > len - 2) i = len - 2;
  double ww = (x - grid[i]) / (grid[i + 1] - grid[i]);
  w = ww < 0.0 ? 0.0 : (ww > 1.0 ? 1.0 : ww);
  return i;
}

struct Grid {
  const double* P;
  const double* S6K;
  double lu[9];
  int offs[16];
  int bits[16][4];
  double P0_ref;
};

static inline void grid_lookup(const Grid& G, int e2i, double u, double Ie, double Ic, int ib, double wb,
                               double* P, double* S) {
  double wu, we, wc;
  int iu = locate(G.lu, 9, std::log1p(u), wu);
  int ie = locate(IE_GRID, 11, Ie, we);
  int ic = locate(IC_GRID, 11, Ic, wc);
  int base = e2i * sE2 + iu * sU + ie * sE + ic * sC + ib;
  double W[16];
  for (int k = 0; k < 16; k++) {
    W[k] = (G.bits[k][0] ? wu : 1 - wu) * (G.bits[k][1] ? we : 1 - we) *
           (G.bits[k][2] ? wc : 1 - wc) * (G.bits[k][3] ? wb : 1 - wb);
  }
  double a[16], b[16];
  for (int gg = 0; gg < 8; gg++) {
    int off = base + gg * sG;
    for (int k = 0; k < 16; k++) {
      a[k] = G.P[off + G.offs[k]] * W[k];
      b[k] = G.S6K[off + G.offs[k]] * W[k];
    }
    P[gg] = np_sum16(a);
    S[gg] = np_sum16(b);
  }
}

// [[Rcpp::export]]
List grid_lookup_cpp(NumericVector Pg, NumericVector Sg, IntegerVector e2i, NumericVector u,
                     NumericVector Ie, NumericVector Ic, NumericVector b0) {
  Grid G;
  G.P = Pg.begin(); G.S6K = Sg.begin();
  for (int j = 0; j < 9; j++) G.lu[j] = std::log1p(U_GRID[j]);
  int k = 0;
  for (int du = 0; du < 2; du++) for (int de = 0; de < 2; de++) for (int dc = 0; dc < 2; dc++)
    for (int db = 0; db < 2; db++) {
      G.offs[k] = du * sU + de * sE + dc * sC + db * sB;
      G.bits[k][0] = du; G.bits[k][1] = de; G.bits[k][2] = dc; G.bits[k][3] = db; k++;
    }
  int n = u.size();
  NumericMatrix P(n, 8), S(n, 8);
  for (int i = 0; i < n; i++) {
    double wb;
    int ib = locate(B0_GRID, 10, b0[i], wb);
    double pp[8], ss[8];
    grid_lookup(G, e2i[i], u[i], Ie[i], Ic[i], ib, wb, pp, ss);
    for (int gg = 0; gg < 8; gg++) { P(i, gg) = pp[gg]; S(i, gg) = ss[gg]; }
  }
  return List::create(_["P"] = P, _["S6K"] = S);
}

// ============================================================================ engine
struct Reg {           // one treatment line
  int cdk_code;        // 0 none, 1 palbociclib, 2 ribociclib, 3 abemaciclib
  double cdk_doses[3]; // dose levels (mg), -1 = none
  double cdk_ref;      // reference dose (mg) for relative exposure
  int cdk_onoff;       // 21 days on / 7 off
  double k_cdk, slope_n;
  int expo_idx;        // 0 palbo, 1 ribo, 2 abema
  int et_code;         // 0 none, 1 letrozole, 2 exemestane, 3 fulvestrant, 4 giredestrant
  double eve_doses[3]; // everolimus dose levels (mg), eve_doses[0] <= 0: no everolimus
};

struct Par {
  // tumor
  double kmax_unused, rd_unused, mu_E, mu_P, mu_R, sigma_R, kappa, N0_med, SLD_ref, meas_cv, lam_NL,
      gam_NL, K_TF, LOD, h_t, N_ref, gam_d, h_post, gomp_a, gomp_b, tau_c, p_attr, early_prog_days,
      p_attr_mid, attr_switch_days;
  // pharmacology
  double Imax_c, IC50_c, h_c, I10_e, IC50_e, uG_ref, uF_ref, k_eve, k_gir, k_fulv, kdeg_3A4, KI_3A4,
      kinact_ribo, kinact_palbo, fm_eve, gir_ddi;
  // toxicity
  double MTT, gam_N, slope_eve, aS_eve, aS_cdk, tauS, p2S, p3S, aP_eve, aP_cdk, aG_eve, aD_abema, aD_eve,
      aD_gir, aD_palbo, aD_ribo, tauD_abema, q3D_abema, q3D_other, aL_ribo, aL_abema, aL_palbo, aL_eve,
      aL_gir, aQ_ribo, tauQ, aB_gir, tauB, aA_eve, disc_cdk, disc_eve, disc_et, disc_gir, tox_int,
      dex_mult, dex_sev, disc_tau, disc_floor, disc_t0, disc_eve_gamma, disc_et_eve, disc_abema_mult;
  bool has_disc_tau;
  // structure
  bool adaptive_R, reversible_R;
};

struct ArmCfg {
  int nlines;
  Reg reg[2];
  double p_attr, p_attr_mid, attr_switch_days, early_prog_days, elig_min_days;
  bool d3_switch, d3_monitor, tox_only;
  int ass[2][4];   // tumor assessment schedule per line: first scan, interval, switch day, later interval
  // ctDNA monitoring: first sample and interval (days after the line starts), detection
  // limit, endocrine therapy after a ctDNA-guided switch, delay from detection to the switch (trial:
  // randomization), and an imaging check at detection and on the switch day (no switch, and the patient
  // flagged, if the tumor already meets the progression criteria). Defaults: 56, 56, tum.LOD,
  // giredestrant, no delay, no check.
  int ct_first, ct_int, switch_et, ct_lag;
  double ct_lod;
  bool ct_confirm;
};

struct VP {   // virtual population (pointers into R vectors)
  const double *kmax, *rd, *b0, *fE0, *fR0, *N0, *expo[3], *expo_eve, *expo_gir, *expo_fulv, *ANC0,
      *frail, *neut_sens, *age, *dc0;
  const int* Ptrunc;
};

struct Out {  // per-patient output pointers
  int *pfs1, *pfs2, *edct, *os, *line2_start, *chemo_start, *esr1_detect, *d3_switch, *resp1_day, *d3_excl;
  double *vaf0, *vaf_det;
  int *pfs1_ev, *pfs2_ev, *edct_ev, *os_ev, *resp1, *resp2, *in_line2, *neut_g3, *neut_g4, *stom_any,
      *stom_g2, *stom_g3, *pneu_any, *pneu_g3, *hyperg_g3, *diar_g2, *diar_g3, *alt_g3, *qtc, *brady,
      *anemia_g3, *esr1_at_p1, *esr1_at_edct, *skip2, *prog1, *ae_disc_cdk, *ae_disc_eve, *ae_disc_et;
  double *tox_days, *cum_plan_cdk, *cum_act_cdk, *cum_plan_eve, *cum_act_eve, *days_cdk, *days_eve,
      *days_gir, *vaf_p1, *vaf_edct, *best_sld1, *best_sld2, *clones_p1, *clones_edct, *N_end;
  int n;
  // trajectories
  int m, W;
  double *trN, *treve, *trcdk, *tranc;
  int* trline;
};

// modulo with the sign of the divisor
static inline int floormod(int a, int b) { int r = a % b; return (r != 0 && ((r < 0) != (b < 0))) ? r + b : r; }

struct PatState {
  // treatment state
  int line_start, cdk_code, cdk_lvl, cdk_hold, et_code, gir_hold, eve_lvl, eve_hold, fulv_start, neut_holds,
      stom_count, pneu_g3c, diar_count, line_days, cdk_days_on, eve_days_on, next_assess, pr_streak;
  double cdk_dose, cdk_ref, cdk_planned, expo_cdk, k_cdk, slope_n, eve_dose, eve_planned;
  bool cdk_disc, eve_disc, cdk_onoff, ae_disc_cdk, ae_disc_eve, ae_disc_et, nl_pending, neut_held;
  double cdk_list[3], eve_list[3];
};

static inline void apply_regimen(PatState& S, const Reg& reg, int t, const VP& v, int i) {
  S.line_start = t;
  if (reg.cdk_code == 0) {
    S.cdk_code = 0; S.cdk_dose = 0.0; S.cdk_planned = 0.0;
  } else {
    S.cdk_code = reg.cdk_code;
    int nd = 0;
    for (int j = 0; j < 3; j++) { S.cdk_list[j] = reg.cdk_doses[j]; if (reg.cdk_doses[j] > 0) nd++; }
    S.cdk_lvl = 0;
    S.cdk_dose = reg.cdk_doses[0];
    S.cdk_ref = reg.cdk_ref;
    S.cdk_planned = reg.cdk_doses[0] / reg.cdk_ref;
    S.cdk_onoff = reg.cdk_onoff != 0;
    S.expo_cdk = v.expo[reg.expo_idx][i];
    S.k_cdk = reg.k_cdk;
    S.slope_n = reg.slope_n;
    S.cdk_disc = false; S.cdk_hold = 0; S.neut_holds = 0; S.cdk_days_on = 0;
  }
  S.et_code = reg.et_code;
  S.gir_hold = 0;
  if (reg.et_code == 3) S.fulv_start = t;
  if (!(reg.eve_doses[0] > 0)) {
    S.eve_dose = 0.0; S.eve_planned = 0.0;
  } else {
    for (int j = 0; j < 3; j++) S.eve_list[j] = reg.eve_doses[j];
    S.eve_lvl = 0;
    S.eve_dose = reg.eve_doses[0];
    S.eve_planned = reg.eve_doses[0] / 10.0;
    S.eve_disc = false; S.eve_hold = 0; S.eve_days_on = 0; S.stom_count = 0; S.pneu_g3c = 0;
  }
  S.diar_count = 0;
  S.line_days = 0;
}

// hold drug for hold_days, then resume at the next lower dose; discontinue if none is left
static inline void reduce_cdk(PatState& S, int t, int hold_days) {
  int lvl = S.cdk_lvl + 1;
  double nxt = lvl < 3 ? S.cdk_list[std::min(lvl, 2)] : -1.0;
  bool disc = nxt <= 0;
  if (!disc) { S.cdk_lvl = lvl; S.cdk_dose = nxt; }
  S.cdk_hold = std::max(S.cdk_hold, t + hold_days);
  if (disc) { S.cdk_disc = true; S.ae_disc_cdk = true; }
}
static inline void reduce_eve(PatState& S, int t, int hold_days) {
  int lvl = S.eve_lvl + 1;
  double nxt = lvl < 3 ? S.eve_list[std::min(lvl, 2)] : -1.0;
  bool disc = nxt <= 0;
  if (!disc) { S.eve_lvl = lvl; S.eve_dose = nxt; }
  S.eve_hold = std::max(S.eve_hold, t + hold_days);
  if (disc) { S.eve_disc = true; S.ae_disc_eve = true; }
}

static void sim_patient(int i, const Par& p, const ArmCfg& arm, const VP& v, const Grid& G, int H,
                        const int* fixed_end, int rng_mode, uint64_t seed, int pid, const double* extU,
                        const double* extZ, int n_all, Out& o) {
  const int NSTREAM = 18;
  const bool tox_only = arm.tox_only;
  PatState S;
  // ---- initial state (before the day loop)
  S.line_start = 0; S.cdk_code = 0; S.cdk_lvl = 0; S.cdk_hold = 0; S.et_code = 0; S.gir_hold = 0;
  S.eve_lvl = 0; S.eve_hold = 0; S.fulv_start = 0; S.neut_holds = 0; S.stom_count = 0; S.pneu_g3c = 0;
  S.diar_count = 0; S.line_days = 0; S.cdk_days_on = 0; S.eve_days_on = 0; S.next_assess = 0; S.pr_streak = 0;
  S.cdk_dose = 0; S.cdk_ref = 1.0; S.cdk_planned = 0; S.expo_cdk = 0; S.k_cdk = 0.5; S.slope_n = 0;
  S.eve_dose = 0; S.eve_planned = 0;
  S.cdk_disc = S.eve_disc = S.cdk_onoff = S.ae_disc_cdk = S.ae_disc_eve = S.ae_disc_et = false;
  S.nl_pending = S.neut_held = false;
  for (int j = 0; j < 3; j++) { S.cdk_list[j] = -1; S.eve_list[j] = -1; }

  // records
  int pfs1 = -1, pfs2 = -1, edct = -1, os = -1, line2_start = -1, chemo_start = -1, esr1_detect = -1,
      d3_switch = -1, resp1_day = -1, d3_excl = -1, d3_due = -1;
  double vaf0 = 0.0, vaf_det = 0.0;
  bool pfs1_ev = false, pfs2_ev = false, edct_ev = false, os_ev = false, resp1 = false, resp2 = false,
       in_line2 = false, neut_g3 = false, neut_g4 = false, stom_any = false, stom_g2 = false, stom_g3 = false,
       pneu_any = false, pneu_g3 = false, hyperg_g3 = false, diar_g2 = false, diar_g3 = false, alt_g3 = false,
       qtc = false, brady = false, anemia_g3 = false, esr1_at_p1 = false, esr1_at_edct = false, skip2 = false,
       prog1 = false;
  double tox_days = 0, cum_plan_cdk = 0, cum_act_cdk = 0, cum_plan_eve = 0, cum_act_eve = 0, days_cdk = 0,
         days_eve = 0, days_gir = 0, vaf_p1 = 0, vaf_edct = 0, best_sld1 = 9.0, best_sld2 = 9.0;
  double clones_p1[8] = {0}, clones_edct[8] = {0};

  int line = 0;
  bool alive = true;
  apply_regimen(S, arm.reg[0], 0, v, i);

  double c_cdk = 0, c_eve = 0, c_gir = 0, c_fulv = 0, cyp = 1.0;
  const double ktr = 4.0 / p.MTT;
  double nP = 1, nT1 = 1, nT2 = 1, nT3 = 1, nC = 1;

  // tumor
  double N[8] = {0};
  {
    int gP = v.Ptrunc[i] ? 2 : 0;
    double fE = v.fE0[i], fR = v.fR0[i], N0 = v.N0[i];
    N[gP] = N0 * (1 - fE) * (1 - fR);
    N[gP + 1] = N0 * fE * (1 - fR);
    N[gP + 4] = N0 * (1 - fE) * fR;
    N[gP + 5] = N0 * fE * fR;
  }
  const double third = 1.0 / 3.0;
  const double sld_scale = p.SLD_ref / std::pow(p.N0_med, third);
  double sld_base = sld_scale * std::pow(v.N0[i], third);
  double nadir = sld_base;
  double nadir_cells = v.N0[i];
  S.next_assess = arm.ass[0][0];
  int ib = 0;
  double wb = 0, d0 = 0;
  if (!tox_only) {
    ib = locate(B0_GRID, 10, v.b0[i], wb);
    d0 = v.rd[i] * v.kmax[i] * G.P0_ref;
  }
  const double h_bg = std::exp(p.gomp_a + p.gomp_b * v.age[i]) / 365.0;
  const double Imax_e = std::min(0.98, p.I10_e * (p.IC50_e + 1.0));
  const double frail = v.frail[i];
  const double kmax = v.kmax[i];
  const double mu_b[3] = {p.mu_E, p.mu_P, p.mu_R};
  const double ee = std::exp(-p.k_eve), eg = std::exp(-p.k_gir), ef = std::exp(-p.k_fulv);
  const bool rec_traj = i < o.m;

  double P[8] = {0}, S6K[8] = {0};
  int last_e2i = 0;
  double last_u = 0, last_Ie = 0, last_Ic = 0;
  double vaf = 0, Ntot = 1.0, last_sld = 0;
  const uint64_t key = patient_key(seed, (uint64_t)pid);

  for (int t = 0; t < H; t++) {
    // ---- random numbers for this patient-day
    double U[NSTREAM];
    double Zm;
    if (rng_mode == 1) {
      const double* base = extU + (size_t)t * NSTREAM * n_all;
      for (int k = 0; k < NSTREAM; k++) U[k] = base[(size_t)k * n_all + pid];
      Zm = extZ[(size_t)t * n_all + pid];
    } else {
      for (int k = 0; k < NSTREAM; k++) U[k] = unif_ctr(key, (uint64_t)t, (uint64_t)k);
      double u1 = unif_ctr(key, (uint64_t)t, 20), u2 = unif_ctr(key, (uint64_t)t, 21);
      Zm = std::sqrt(-2.0 * std::log(u1)) * std::cos(6.283185307179586 * u2);
    }
    bool endo = alive && (line <= 1);
    if (tox_only) endo = endo && (t < fixed_end[i]);
    int dline = t - S.line_start;

    // ---------------- dosing
    bool sched_on = S.cdk_onoff ? (floormod(dline, 28) < 21) : true;
    bool cdk_planned_today = endo && (S.cdk_code > 0) && sched_on;
    double cdk_rate = (cdk_planned_today && !S.cdk_disc && (t >= S.cdk_hold)) ? S.cdk_dose / S.cdk_ref : 0.0;
    bool eve_planned_today = endo && (S.eve_planned > 0);
    double eve_rate = (eve_planned_today && !S.eve_disc && (t >= S.eve_hold)) ? S.eve_dose / 10.0 : 0.0;
    bool gir_on = endo && (S.et_code == 4) && (t >= S.gir_hold);
    bool fulv_on = endo && (S.et_code == 3);
    cum_plan_cdk += cdk_planned_today ? S.cdk_planned : 0;
    cum_act_cdk += cdk_rate;
    cum_plan_eve += eve_planned_today ? S.eve_planned : 0;
    cum_act_eve += eve_rate;
    days_cdk += cdk_planned_today;
    days_eve += eve_planned_today;
    days_gir += gir_on;

    // ---------------- PK and DDI
    double kin = (S.cdk_code == 2 ? p.kinact_ribo : 0.0) + (S.cdk_code == 1 ? p.kinact_palbo : 0.0);
    kin = kin * c_cdk / (p.KI_3A4 + c_cdk);
    cyp = (cyp + p.kdeg_3A4) / (1 + p.kdeg_3A4 + kin);
    double aucr = 1.0 / (p.fm_eve * cyp + 1 - p.fm_eve);
    aucr = aucr * (gir_on ? p.gir_ddi : 1.0);
    double ek = std::exp(-S.k_cdk);
    c_cdk = c_cdk * ek + cdk_rate * S.expo_cdk * (1 - ek);
    c_eve = c_eve * ee + eve_rate * v.expo_eve[i] * aucr * (1 - ee);
    c_gir = c_gir * eg + (gir_on ? 1.0 : 0.0) * v.expo_gir[i] * (1 - eg);
    double fr = fulv_on ? ((t - S.fulv_start < 28) ? 2.0 : 1.0) : 0.0;
    c_fulv = c_fulv * ef + fr * v.expo_fulv[i] * (1 - ef);

    double cch = std::pow(c_cdk, p.h_c);
    double Ic = p.Imax_c * cch / (std::pow(p.IC50_c, p.h_c) + cch);
    double Ie = Imax_e * c_eve / (p.IC50_e + c_eve);
    double u = p.uG_ref * c_gir + p.uF_ref * c_fulv;
    int e2i = (endo && S.et_code == 1) ? 2 : ((endo && S.et_code == 2) ? 1 : 0);

    // ---------------- tumor
    bool prog = false;
    if (!tox_only) {
      if (t == 0) {
        grid_lookup(G, e2i, u, Ie, Ic, ib, wb, P, S6K);
        last_e2i = e2i; last_u = u; last_Ie = Ie; last_Ic = Ic;
      } else {
        bool ch = (e2i != last_e2i) || (std::fabs(std::log1p(u) - std::log1p(last_u)) > 0.01) ||
                  (std::fabs(Ie - last_Ie) > 0.004) || (std::fabs(Ic - last_Ic) > 0.004);
        if (ch) {
          grid_lookup(G, e2i, u, Ie, Ic, ib, wb, P, S6K);
          last_e2i = e2i; last_u = u; last_Ie = Ie; last_Ic = Ic;
        }
      }
      double gr[8];
      for (int k = 0; k < 8; k++) gr[k] = kmax * P[k] - d0;
      if (alive && line == 2) {
        double dc = v.dc0[i] * std::exp(-(double)(t - chemo_start) / p.tau_c);
        for (int k = 0; k < 8; k++) gr[k] = gr[k] - dc;
      }
      for (int k = 0; k < 8; k++) N[k] = N[k] * std::exp(gr[k]);
      double divs[8];
      for (int k = 0; k < 8; k++) divs[k] = kmax * P[k];
      for (int bi = 0; bi < 3; bi++) {
        int b = 1 << bi;
        double tr[8];
        for (int s = 0; s < 8; s++) {
          if (s & b) continue;
          double rate = mu_b[bi] * divs[s];
          if (b == 4 && p.adaptive_R)
            rate = rate + p.sigma_R * Ic * (1 - p.kappa + p.kappa * std::min(S6K[s], 1.5));
          tr[s] = std::min(rate, 0.5) * N[s];
        }
        for (int s = 0; s < 8; s++) {
          if (s & b) continue;
          N[s] -= tr[s];
          N[s | b] += tr[s];
        }
      }
      if (p.reversible_R) {
        double back[4];
        for (int r = 0; r < 4; r++) back[r] = 0.01 * (1 - Ic) * N[4 + r];
        for (int r = 0; r < 4; r++) N[4 + r] -= back[r];
        for (int r = 0; r < 4; r++) N[r] += back[r];
      }
      if (!alive) for (int k = 0; k < 8; k++) N[k] = 0.0;
      Ntot = np_sum8(N);

      // new lesions (detected at the next assessment)
      double h_nl = p.lam_NL * std::pow(Ntot / std::max(nadir_cells, 1.0), p.gam_NL);
      S.nl_pending = S.nl_pending || (endo && (U[1] < h_nl));

      // ctDNA monitoring
      double fracE = np_sum4(N[1], N[3], N[5], N[7]) / std::max(Ntot, 1.0);
      double TF = Ntot / (Ntot + p.K_TF);
      vaf = 0.5 * fracE * TF;
      if (t == 0) vaf0 = vaf;
      bool ct_day = endo && (dline > 0) && (dline >= arm.ct_first) && (floormod(dline - arm.ct_first, arm.ct_int) == 0);
      bool newdet = ct_day && (vaf >= arm.ct_lod) && (esr1_detect < 0);
      if (newdet) { esr1_detect = t; vaf_det = vaf; }
      // ctDNA-guided switch (d3_switch) or monitoring only (d3_monitor), line 1 only. The switch (in a trial,
      // the randomization) takes place ct_lag days after detection. With ct_confirm, a tumor that meets the
      // progression criteria (noise-free size, or a pending new lesion) at detection or on the switch day is
      // not switched or randomized (d3_excl); progression before the switch day also prevents it.
      if ((arm.d3_switch || arm.d3_monitor) && line == 0 && d3_switch < 0 && d3_excl < 0) {
        bool due = (newdet && arm.ct_lag == 0) || (d3_due >= 0 && t == d3_due);
        if (newdet && arm.ct_lag > 0) d3_due = t + arm.ct_lag;
        if (arm.ct_confirm && (newdet || due)) {
          double sld_now = sld_scale * std::pow(Ntot, third);
          bool pg = ((sld_now >= 1.2 * nadir) && (sld_now - nadir >= 5.0)) || S.nl_pending;
          if (pg) { d3_excl = t; d3_due = -1; due = false; }
        }
        if (due) {
          if (arm.d3_switch) {
            if (S.et_code != arm.switch_et) {
              d3_switch = t; S.et_code = arm.switch_et; S.gir_hold = 0;
              if (arm.switch_et == 3) S.fulv_start = t;
            }
          } else {
            d3_switch = t;
          }
        }
      }

      // assessments
      bool ass = endo && (t == S.next_assess);
      if (ass) {
        double sld = sld_scale * std::pow(Ntot, third) * (1 + p.meas_cv * Zm);
        double rel = sld / std::max(sld_base, 1e-9);
        if (line == 0) best_sld1 = std::min(best_sld1, rel);
        if (line == 1) best_sld2 = std::min(best_sld2, rel);
        S.pr_streak = rel <= 0.7 ? S.pr_streak + 1 : 0;
        bool conf = S.pr_streak >= 2;
        if (conf && line == 0 && !resp1) resp1_day = t;
        if (conf && line == 0) resp1 = true;
        if (conf && line == 1) resp2 = true;
        prog = ((sld >= 1.2 * nadir) && (sld - nadir >= 5.0)) || S.nl_pending;
        if (!prog) {
          nadir = std::min(nadir, sld);
          nadir_cells = std::min(nadir_cells, Ntot);
        }
        const int* sch = arm.ass[line == 1 ? 1 : 0];
        int nxt = dline < sch[2] ? sch[1] : sch[3];
        S.next_assess = t + nxt;
        last_sld = sld;
      }
    } else {
      Ntot = 1.0;
    }

    // ---------------- death
    bool die = false;
    if (!tox_only) {
      double h_d = h_bg + p.h_t * std::pow(Ntot / p.N_ref, p.gam_d) + p.h_post * (line == 2 ? 1.0 : 0.0);
      die = alive && (U[0] < h_d);
    }

    // ---------------- toxicity
    {
      double E = (S.slope_n * c_cdk + p.slope_eve * c_eve) * v.neut_sens[i];
      E = std::min(E, 0.99);
      for (int sub = 0; sub < 2; sub++) {
        double fb = std::pow(1.0 / std::max(nC, 1e-3), p.gam_N);
        double dP = ktr * nP * (1 - E) * fb - ktr * nP;
        double dT1 = ktr * (nP - nT1);
        double dT2 = ktr * (nT1 - nT2);
        double dT3 = ktr * (nT2 - nT3);
        double dC = ktr * (nT3 - nC);
        nP += 0.5 * dP; nT1 += 0.5 * dT1; nT2 += 0.5 * dT2; nT3 += 0.5 * dT3; nC += 0.5 * dC;
      }
      double anc = v.ANC0[i] * nC;
      bool on_cdk = endo && (S.cdk_code > 0) && !S.cdk_disc;
      bool check = on_cdk && ((floormod(dline, 28) == 0) || (dline == 14) || (dline == 42));
      bool heldchk = on_cdk && S.neut_held && (floormod(t - S.cdk_hold, 7) == 0);
      bool chk = check || heldchk;
      if (chk) {
        bool g3 = anc < 1.0, g4 = anc < 0.5;
        neut_g3 = neut_g3 || g3;
        neut_g4 = neut_g4 || g4;
        bool recover = heldchk && (anc >= 1.0);
        if (recover) { S.neut_held = false; S.cdk_hold = t; }
        bool new3 = g3 && !S.neut_held;
        if (new3) {
          S.neut_holds += 1;
          bool red = g4 || (S.neut_holds >= 3);
          S.neut_held = true;
          S.cdk_hold = t + 400;   // held until recovery check
          if (red) {
            int lvl = S.cdk_lvl + 1;
            double nxt = lvl < 3 ? S.cdk_list[std::min(lvl, 2)] : -1.0;
            bool stay = nxt <= 0;   // at the lowest dose level neutropenia leads to holds only
            if (!stay) { S.cdk_lvl = lvl; S.cdk_dose = nxt; }
            S.neut_holds = 0;
          }
        }
      }

      // event hazards
      double t_int = p.tox_int;
      bool on_eve = endo && (S.eve_planned > 0) && !S.eve_disc;
      bool on_gir = endo && (S.et_code == 4);
      int code = S.cdk_code;
      bool has_cdk = endo && (code > 0) && !S.cdk_disc;
      int ndrug = (int)has_cdk + (int)on_eve + (int)endo;
      double mult = (ndrug >= 3 ? t_int : 1.0) * frail;
      S.eve_days_on += on_eve;
      S.line_days += endo;
      double c_ribo = code == 2 ? c_cdk : 0.0;
      double c_abema = code == 3 ? c_cdk : 0.0;
      double c_palbo = code == 1 ? c_cdk : 0.0;
      double cdk_any = has_cdk ? std::min(c_cdk, 1.5) : 0.0;
      double gir_f = on_gir ? 1.0 : 0.0;

      // stomatitis
      if (endo) {
        double h = (p.aS_eve * c_eve + p.aS_cdk * cdk_any) * std::exp(-S.line_days / p.tauS) * mult * p.dex_mult;
        if (U[2] < h) {
          stom_any = true;
          double p2 = p.p2S * p.dex_sev, p3 = p.p3S * p.dex_sev;
          bool g3 = U[3] < p3;
          bool g2 = !g3 && (U[3] < p3 + p2);
          stom_g2 = stom_g2 || g2 || g3;
          stom_g3 = stom_g3 || g3;
          tox_days += g3 ? 14 : (g2 ? 10 : 0);
          if (g3 && on_eve) reduce_eve(S, t, 14);
          bool e2 = g2 && on_eve;
          bool first = e2 && (S.stom_count == 0);
          if (first) S.eve_hold = std::max(S.eve_hold, t + 7);
          if (e2 && !first) reduce_eve(S, t, 7);
          S.stom_count += ((g2 || g3) && on_eve) ? 1 : 0;
        }
      }
      // pneumonitis
      if (endo) {
        double h = (p.aP_eve * c_eve + p.aP_cdk * cdk_any) * mult;
        if (U[4] < h) {
          pneu_any = true;
          int gr = U[5] < 0.45 ? 1 : (U[5] < 0.80 ? 2 : (U[5] < 0.97 ? 3 : 4));
          pneu_g3 = pneu_g3 || (gr >= 3);
          tox_days += gr >= 3 ? 28 : (gr == 2 ? 21 : 0);
          bool g4 = (gr == 4) && on_eve;
          bool g3rec = (gr == 3) && on_eve && (S.pneu_g3c >= 1);
          if (g4 || g3rec) { S.eve_disc = true; S.ae_disc_eve = true; }
          S.pneu_g3c += (gr == 3) ? 1 : 0;
          if (on_eve && gr == 3 && !g3rec) reduce_eve(S, t, 28);
          if (on_eve && gr == 2) reduce_eve(S, t, 21);
        }
      }
      // hyperglycemia grade >= 3
      if (endo && (U[6] < p.aG_eve * c_eve * frail)) {
        hyperg_g3 = true;
        tox_days += 7;
        if (on_eve) reduce_eve(S, t, 7);
      }
      // diarrhea grade >= 2
      if (endo) {
        double h = (p.aD_abema * c_abema * std::exp(-S.line_days / p.tauD_abema) + p.aD_eve * c_eve +
                    p.aD_gir * c_gir * gir_f + p.aD_palbo * c_palbo + p.aD_ribo * c_ribo) * mult;
        if (U[7] < h) {
          double q3 = code == 3 ? p.q3D_abema : p.q3D_other;
          bool g3 = U[8] < q3;
          diar_g2 = true;
          diar_g3 = diar_g3 || g3;
          tox_days += g3 ? 7 : 5;
          bool ab = (code == 3) && has_cdk;
          bool first = ab && !g3 && (S.diar_count == 0);
          if (first) S.cdk_hold = std::max(S.cdk_hold, t + 5);
          if (ab && !first) reduce_cdk(S, t, 7);
          if (g3 && !ab && on_eve) reduce_eve(S, t, 7);
          S.diar_count += 1;
        }
      }
      // ALT/AST grade >= 3
      if (endo) {
        double h = (p.aL_ribo * c_ribo + p.aL_abema * c_abema + p.aL_palbo * c_palbo + p.aL_eve * c_eve +
                    p.aL_gir * c_gir * gir_f) * frail;
        if (U[9] < h) {
          alt_g3 = true;
          if (has_cdk) reduce_cdk(S, t, 14);
          if (!has_cdk && on_eve) reduce_eve(S, t, 14);
          if (!has_cdk && !on_eve && on_gir) S.gir_hold = t + 14;
        }
      }
      // QTc (ribociclib)
      if (endo && (U[10] < p.aQ_ribo * c_ribo * std::exp(-S.line_days / p.tauQ))) {
        qtc = true;
        if (has_cdk) reduce_cdk(S, t, 7);
      }
      // counted-only events
      if (on_gir && (U[11] < p.aB_gir * c_gir * std::exp(-S.line_days / p.tauB))) brady = true;
      if (endo && (U[12] < p.aA_eve * c_eve * frail)) { anemia_g3 = true; tox_days += 7; }
      // other AE-related discontinuations. Optional: the hazard falls from disc_t0 days
      // after the line starts towards disc_floor of its initial value; everolimus discontinuation scales
      // with everolimus exposure (c_eve ^ disc_eve_gamma); abemaciclib has its own multiplier of the
      // CDK4/6i hazard (disc_abema_mult); an extra endocrine-therapy hazard applies
      // while everolimus is given. Engine defaults if not set: t0 183 d, gamma 0, multiplier 1, extra 0.
      double dmult = mult;
      if (p.has_disc_tau) {
        double fl = p.disc_floor;
        double late = std::max((double)S.line_days - p.disc_t0, 0.0);
        dmult = mult * (fl + (1 - fl) * std::exp(-late / p.disc_tau));
      }
      double eve_f = p.disc_eve_gamma != 0 ? std::pow(std::max(c_eve, 0.0), p.disc_eve_gamma) : 1.0;
      double cdk_dm = (code == 3) ? p.disc_abema_mult : 1.0;
      if (has_cdk && (U[13] < p.disc_cdk * cdk_dm * dmult)) { S.cdk_disc = true; S.ae_disc_cdk = true; }
      if (on_eve && (U[14] < p.disc_eve * eve_f * dmult)) { S.eve_disc = true; S.ae_disc_eve = true; }
      double et_h = (on_gir ? p.disc_gir : p.disc_et) + (on_eve ? p.disc_et_eve : 0.0);
      if (endo && (S.et_code > 0) && (U[15] < et_h * dmult)) {
        S.ae_disc_et = true;
        if (S.et_code == 4 || S.et_code == 3) S.et_code = 1;   // SERD -> letrozole; AI -> other AI
      }

      // ---------------- trajectories
      if (rec_traj && (t % 7 == 0)) {
        int w = t / 7;
        size_t bw = (size_t)w * o.m + i;
        for (int k = 0; k < 8; k++) o.trN[bw * 8 + k] = tox_only ? 0.0 : N[k];
        o.trline[bw] = alive ? line : 9;
        o.treve[bw] = c_eve;
        o.trcdk[bw] = c_cdk;
        o.tranc[bw] = v.ANC0[i] * nC;
      }
    }

    // ---------------- progression / line transitions
    prog = prog && !die;
    if (prog) {
      bool p0 = line == 0, p1 = line == 1;
      bool go2 = false;
      if (p0) {
        pfs1 = t; pfs1_ev = true; prog1 = true;
        vaf_p1 = vaf; esr1_at_p1 = vaf >= arm.ct_lod;
        for (int k = 0; k < 8; k++) clones_p1[k] = N[k] / Ntot;
        if (arm.nlines > 1) {
          bool early = (t - S.line_start) < arm.early_prog_days;
          // attrition (no endocrine-based 2nd line) may depend on time on line 1
          double pa = ((t - S.line_start) < arm.attr_switch_days) ? arm.p_attr_mid : arm.p_attr;
          go2 = !early && (U[16] >= pa);
          if (arm.elig_min_days > 0) go2 = go2 && ((t - S.line_start) >= arm.elig_min_days);
        }
        skip2 = skip2 || !go2;
        if (go2) {
          line2_start = t; in_line2 = true; line = 1;
          apply_regimen(S, arm.reg[1], t, v, i);
          sld_base = last_sld; nadir = last_sld; nadir_cells = Ntot;
          S.next_assess = t + arm.ass[1][0];
          S.pr_streak = 0; S.nl_pending = false;
        }
      }
      if (p1) { pfs2 = t; pfs2_ev = true; }
      bool end = (p0 && !go2) || p1;
      if (end) {
        edct = t; edct_ev = true;
        vaf_edct = vaf; esr1_at_edct = vaf >= arm.ct_lod;
        for (int k = 0; k < 8; k++) clones_edct[k] = N[k] / Ntot;
        line = 2; chemo_start = t;
        S.eve_planned = 0; S.cdk_code = 0; S.et_code = 0;
      }
    }
    if (die) {
      os = t; os_ev = true;
      if (pfs1 < 0) { pfs1 = t; pfs1_ev = true; }
      if (edct < 0) { edct = t; edct_ev = true; }
      if (line == 1 && pfs2 < 0) { pfs2 = t; pfs2_ev = true; }
      alive = false;
      line = 9;
    }
    // nothing changes for a dead patient (or after a toxicity-only run ends), except recorded
    // trajectories
    if (!rec_traj) {
      if (!alive) break;
      if (tox_only && t + 1 >= fixed_end[i]) break;
    }
  }
  // censoring at the horizon
  if (pfs1 < 0) pfs1 = H;
  if (edct < 0) edct = H;
  if (os < 0) os = H;
  if (in_line2 && pfs2 < 0) pfs2 = H;
  if (!alive) for (int k = 0; k < 8; k++) N[k] = 0.0;

  o.pfs1[i] = pfs1; o.pfs2[i] = pfs2; o.edct[i] = edct; o.os[i] = os; o.line2_start[i] = line2_start;
  o.chemo_start[i] = chemo_start; o.esr1_detect[i] = esr1_detect; o.d3_switch[i] = d3_switch;
  o.resp1_day[i] = resp1_day; o.d3_excl[i] = d3_excl; o.vaf0[i] = vaf0; o.vaf_det[i] = vaf_det;
  o.pfs1_ev[i] = pfs1_ev; o.pfs2_ev[i] = pfs2_ev; o.edct_ev[i] = edct_ev; o.os_ev[i] = os_ev;
  o.resp1[i] = resp1; o.resp2[i] = resp2; o.in_line2[i] = in_line2; o.neut_g3[i] = neut_g3;
  o.neut_g4[i] = neut_g4; o.stom_any[i] = stom_any; o.stom_g2[i] = stom_g2; o.stom_g3[i] = stom_g3;
  o.pneu_any[i] = pneu_any; o.pneu_g3[i] = pneu_g3; o.hyperg_g3[i] = hyperg_g3; o.diar_g2[i] = diar_g2;
  o.diar_g3[i] = diar_g3; o.alt_g3[i] = alt_g3; o.qtc[i] = qtc; o.brady[i] = brady;
  o.anemia_g3[i] = anemia_g3; o.esr1_at_p1[i] = esr1_at_p1; o.esr1_at_edct[i] = esr1_at_edct;
  o.skip2[i] = skip2; o.prog1[i] = prog1; o.ae_disc_cdk[i] = S.ae_disc_cdk;
  o.ae_disc_eve[i] = S.ae_disc_eve; o.ae_disc_et[i] = S.ae_disc_et;
  o.tox_days[i] = tox_days; o.cum_plan_cdk[i] = cum_plan_cdk; o.cum_act_cdk[i] = cum_act_cdk;
  o.cum_plan_eve[i] = cum_plan_eve; o.cum_act_eve[i] = cum_act_eve; o.days_cdk[i] = days_cdk;
  o.days_eve[i] = days_eve; o.days_gir[i] = days_gir; o.vaf_p1[i] = vaf_p1; o.vaf_edct[i] = vaf_edct;
  o.best_sld1[i] = best_sld1; o.best_sld2[i] = best_sld2;
  for (int k = 0; k < 8; k++) {   // column-major n x 8 matrices
    o.clones_p1[(size_t)k * o.n + i] = clones_p1[k];
    o.clones_edct[(size_t)k * o.n + i] = clones_edct[k];
    o.N_end[(size_t)k * o.n + i] = N[k];
  }
}

static Reg make_reg(const List& r) {
  Reg g;
  g.cdk_code = as<int>(r["cdk_code"]);
  NumericVector cd = r["cdk_doses"];
  for (int j = 0; j < 3; j++) g.cdk_doses[j] = j < cd.size() ? cd[j] : -1.0;
  g.cdk_ref = as<double>(r["cdk_ref"]);
  g.cdk_onoff = as<int>(r["cdk_onoff"]);
  g.k_cdk = as<double>(r["k_cdk"]);
  g.slope_n = as<double>(r["slope_n"]);
  g.expo_idx = as<int>(r["expo_idx"]);
  g.et_code = as<int>(r["et_code"]);
  NumericVector ed = r["eve_doses"];
  for (int j = 0; j < 3; j++) g.eve_doses[j] = j < ed.size() ? ed[j] : -1.0;
  return g;
}

// [[Rcpp::export]]
List simulate_cpp(List tum, List ph, List tox, List st, List vp, List grid, List arm, int H,
                  int record_traj = 0, int rng_mode = 0, double seed = 1.0,
                  Nullable<IntegerVector> pid_ = R_NilValue,
                  Nullable<NumericVector> extU_ = R_NilValue, Nullable<NumericVector> extZ_ = R_NilValue,
                  int n_ext = 0, int nthreads = 1) {
  // ---- parameters
  Par p;
  p.mu_E = num(tum, "mu_E"); p.mu_P = num(tum, "mu_P"); p.mu_R = num(tum, "mu_R");
  p.sigma_R = num(tum, "sigma_R"); p.kappa = num(tum, "kappa"); p.N0_med = num(tum, "N0_med");
  p.SLD_ref = num(tum, "SLD_ref"); p.meas_cv = num(tum, "meas_cv"); p.lam_NL = num(tum, "lam_NL");
  p.gam_NL = num(tum, "gam_NL"); p.K_TF = num(tum, "K_TF"); p.LOD = num(tum, "LOD"); p.h_t = num(tum, "h_t");
  p.N_ref = num(tum, "N_ref"); p.gam_d = num(tum, "gam_d"); p.h_post = num_or(tum, "h_post", 0.0);
  p.gomp_a = num(tum, "gomp_a"); p.gomp_b = num(tum, "gomp_b"); p.tau_c = num(tum, "tau_c");
  p.p_attr = num(tum, "p_attr"); p.early_prog_days = num(tum, "early_prog_days");
  p.p_attr_mid = num_or(tum, "p_attr_mid", p.p_attr);
  p.attr_switch_days = num_or(tum, "attr_switch_days", 365.0);
  p.Imax_c = num(ph, "Imax_c"); p.IC50_c = num(ph, "IC50_c"); p.h_c = num(ph, "h_c"); p.I10_e = num(ph, "I10_e");
  p.IC50_e = num(ph, "IC50_e"); p.uG_ref = num(ph, "uG_ref"); p.uF_ref = num(ph, "uF_ref");
  p.k_eve = num(ph, "k_eve"); p.k_gir = num(ph, "k_gir"); p.k_fulv = num(ph, "k_fulv");
  p.kdeg_3A4 = num(ph, "kdeg_3A4"); p.KI_3A4 = num(ph, "KI_3A4"); p.kinact_ribo = num(ph, "kinact_ribo");
  p.kinact_palbo = num(ph, "kinact_palbo"); p.fm_eve = num(ph, "fm_eve"); p.gir_ddi = num(ph, "gir_ddi");
  p.MTT = num(tox, "MTT"); p.gam_N = num(tox, "gam_N"); p.slope_eve = num(tox, "slope_eve");
  p.aS_eve = num(tox, "aS_eve"); p.aS_cdk = num(tox, "aS_cdk"); p.tauS = num(tox, "tauS");
  p.p2S = num(tox, "p2S"); p.p3S = num(tox, "p3S"); p.aP_eve = num(tox, "aP_eve"); p.aP_cdk = num(tox, "aP_cdk");
  p.aG_eve = num(tox, "aG_eve"); p.aD_abema = num(tox, "aD_abema"); p.aD_eve = num(tox, "aD_eve");
  p.aD_gir = num(tox, "aD_gir"); p.aD_palbo = num(tox, "aD_palbo"); p.aD_ribo = num(tox, "aD_ribo");
  p.tauD_abema = num(tox, "tauD_abema"); p.q3D_abema = num(tox, "q3D_abema"); p.q3D_other = num(tox, "q3D_other");
  p.aL_ribo = num(tox, "aL_ribo"); p.aL_abema = num(tox, "aL_abema"); p.aL_palbo = num(tox, "aL_palbo");
  p.aL_eve = num(tox, "aL_eve"); p.aL_gir = num(tox, "aL_gir"); p.aQ_ribo = num(tox, "aQ_ribo");
  p.tauQ = num(tox, "tauQ"); p.aB_gir = num(tox, "aB_gir"); p.tauB = num(tox, "tauB"); p.aA_eve = num(tox, "aA_eve");
  p.disc_cdk = num(tox, "disc_cdk"); p.disc_eve = num(tox, "disc_eve"); p.disc_et = num(tox, "disc_et");
  p.disc_gir = num(tox, "disc_gir"); p.tox_int = num(tox, "tox_int");
  p.dex_mult = num_or(tox, "dex_mult", 1.0); p.dex_sev = num_or(tox, "dex_sev", 1.0);
  p.has_disc_tau = tox.containsElementNamed("disc_tau") && !Rf_isNull(tox["disc_tau"]) &&
                   !ISNAN(as<double>(tox["disc_tau"]));
  p.disc_tau = num_or(tox, "disc_tau", 1.0);
  p.disc_floor = num_or(tox, "disc_floor", 0.25);
  p.disc_t0 = num_or(tox, "disc_t0", 183.0);
  p.disc_eve_gamma = num_or(tox, "disc_eve_gamma", 0.0);
  p.disc_et_eve = num_or(tox, "disc_et_eve", 0.0);
  p.disc_abema_mult = num_or(tox, "disc_abema_mult", 1.0);
  p.adaptive_R = flag_or(st, "adaptive_R", true);
  p.reversible_R = flag_or(st, "reversible_R", false);

  // ---- arm
  ArmCfg A;
  List lines = arm["lines"];
  A.nlines = lines.size();
  if (A.nlines < 1 || A.nlines > 2) stop("arm must have 1 or 2 lines");
  for (int j = 0; j < A.nlines; j++) A.reg[j] = make_reg(as<List>(lines[j]));
  A.p_attr = num_or(arm, "p_attr", p.p_attr);
  // an arm that sets p_attr (trial emulations: 0) applies it to every line-1 duration unless it
  // also sets p_attr_mid
  A.p_attr_mid = num_or(arm, "p_attr_mid", arm.containsElementNamed("p_attr") ? A.p_attr : p.p_attr_mid);
  A.attr_switch_days = num_or(arm, "attr_switch_days", p.attr_switch_days);
  A.early_prog_days = num_or(arm, "early_prog_days", p.early_prog_days);
  A.elig_min_days = num_or(arm, "elig_min_days", 0.0);
  A.d3_switch = flag_or(arm, "d3_switch", false);
  A.d3_monitor = flag_or(arm, "d3_monitor", false);
  A.ct_first = (int)num_or(arm, "ct_first", 56.0);
  A.ct_int = (int)num_or(arm, "ct_int", 56.0);
  if (A.ct_first < 1 || A.ct_int < 1) stop("ctDNA schedule: ct_first and ct_int must be >= 1 day");
  A.ct_lod = num_or(arm, "ct_lod", p.LOD);
  A.switch_et = (int)num_or(arm, "switch_et", 4.0);
  if (A.switch_et != 3 && A.switch_et != 4) stop("switch_et must be 3 (fulvestrant) or 4 (giredestrant)");
  A.ct_confirm = flag_or(arm, "ct_confirm", false);
  A.ct_lag = (int)num_or(arm, "ct_lag", 0.0);
  if (A.ct_lag < 0) stop("ct_lag must be >= 0 days");
  // default schedule: first scan at 8 weeks, every 8 weeks for 18 months, then
  // every 12 weeks, restarting at the start of line 2
  for (int j = 0; j < 2; j++) {
    A.ass[j][0] = 56; A.ass[j][1] = 56; A.ass[j][2] = 540; A.ass[j][3] = 84;
    const char* nm = j == 0 ? "ass1" : "ass2";
    if (arm.containsElementNamed(nm) && !Rf_isNull(arm[nm])) {
      IntegerVector a = as<IntegerVector>(arm[nm]);
      if (a.size() != 4) stop("assessment schedule must be c(first, interval, switch_day, later_interval)");
      for (int k = 0; k < 4; k++) A.ass[j][k] = a[k];
    }
  }
  std::vector<int> fixed_end;
  A.tox_only = arm.containsElementNamed("fixed_duration") && !Rf_isNull(arm["fixed_duration"]);

  // ---- population
  VP v;
  NumericVector kmax = vp["kmax"], rd = vp["rd"], b0 = vp["b0"], fE0 = vp["fE0"], fR0 = vp["fR0"],
                N0 = vp["N0"], e_palbo = vp["expo_palbo"], e_ribo = vp["expo_ribo"], e_abema = vp["expo_abema"],
                e_eve = vp["expo_eve"], e_gir = vp["expo_gir"], e_fulv = vp["expo_fulv"], ANC0 = vp["ANC0"],
                frail = vp["tox_frailty"], ns = vp["neut_sens"], age = vp["age"], dc0 = vp["dc0"];
  IntegerVector Ptr = as<IntegerVector>(vp["Ptrunc"]);
  const int n = kmax.size();
  v.kmax = kmax.begin(); v.rd = rd.begin(); v.b0 = b0.begin(); v.fE0 = fE0.begin(); v.fR0 = fR0.begin();
  v.N0 = N0.begin(); v.expo[0] = e_palbo.begin(); v.expo[1] = e_ribo.begin(); v.expo[2] = e_abema.begin();
  v.expo_eve = e_eve.begin(); v.expo_gir = e_gir.begin(); v.expo_fulv = e_fulv.begin(); v.ANC0 = ANC0.begin();
  v.frail = frail.begin(); v.neut_sens = ns.begin(); v.age = age.begin(); v.dc0 = dc0.begin();
  v.Ptrunc = Ptr.begin();
  if (A.tox_only) {
    IntegerVector fd = as<IntegerVector>(arm["fixed_duration"]);
    if (fd.size() != n) stop("fixed_duration must have one value per patient");
    fixed_end.assign(fd.begin(), fd.end());
  } else {
    fixed_end.assign(n, H);
  }
  IntegerVector pid(n);
  if (pid_.isNotNull()) pid = IntegerVector(pid_); else for (int i = 0; i < n; i++) pid[i] = i;

  // ---- grid
  Grid G;
  NumericVector Pg, Sg;
  if (!A.tox_only) {
    Pg = grid["P"]; Sg = grid["S6K"];
    G.P = Pg.begin(); G.S6K = Sg.begin();
    G.P0_ref = as<double>(grid["P0_ref"]);
  } else {
    G.P = nullptr; G.S6K = nullptr; G.P0_ref = 1.0;
  }
  for (int j = 0; j < 9; j++) G.lu[j] = std::log1p(U_GRID[j]);
  {
    int k = 0;
    for (int du = 0; du < 2; du++) for (int de = 0; de < 2; de++) for (int dc = 0; dc < 2; dc++)
      for (int db = 0; db < 2; db++) {
        G.offs[k] = du * sU + de * sE + dc * sC + db * sB;
        G.bits[k][0] = du; G.bits[k][1] = de; G.bits[k][2] = dc; G.bits[k][3] = db; k++;
      }
  }
  // ---- external random numbers (rng_mode = 1)
  const double* extU = nullptr;
  const double* extZ = nullptr;
  NumericVector eU, eZ;
  if (rng_mode == 1) {
    if (extU_.isNull() || extZ_.isNull()) stop("rng_mode = 1 needs extU and extZ");
    eU = NumericVector(extU_); eZ = NumericVector(extZ_);
    extU = eU.begin(); extZ = eZ.begin();
    if ((size_t)eU.size() < (size_t)H * 18 * n_ext) stop("extU too short for the horizon");
  }

  // ---- outputs
  IntegerVector pfs1(n), pfs2(n), edct(n), os(n), line2_start(n), chemo_start(n), esr1_detect(n),
      d3_switch(n), resp1_day(n), d3_excl(n);
  NumericVector vaf0(n), vaf_det(n);
  LogicalVector pfs1_ev(n), pfs2_ev(n), edct_ev(n), os_ev(n), resp1(n), resp2(n), in_line2(n), neut_g3(n),
      neut_g4(n), stom_any(n), stom_g2(n), stom_g3(n), pneu_any(n), pneu_g3(n), hyperg_g3(n), diar_g2(n),
      diar_g3(n), alt_g3(n), qtc(n), brady(n), anemia_g3(n), esr1_at_p1(n), esr1_at_edct(n), skip2(n),
      prog1(n), ae_disc_cdk(n), ae_disc_eve(n), ae_disc_et(n);
  NumericVector tox_days(n), cum_plan_cdk(n), cum_act_cdk(n), cum_plan_eve(n), cum_act_eve(n), days_cdk(n),
      days_eve(n), days_gir(n), vaf_p1(n), vaf_edct(n), best_sld1(n), best_sld2(n);
  NumericMatrix clones_p1(n, 8), clones_edct(n, 8), N_end(n, 8);
  int m = std::min(record_traj, n);
  int W = H / 7 + 1;
  NumericVector trN((size_t)std::max(m, 0) * W * 8), treve((size_t)std::max(m, 0) * W),
      trcdk((size_t)std::max(m, 0) * W), tranc((size_t)std::max(m, 0) * W);
  IntegerVector trline((size_t)std::max(m, 0) * W);
  Out o;
  o.n = n; o.m = m; o.W = W;
  o.pfs1 = pfs1.begin(); o.pfs2 = pfs2.begin(); o.edct = edct.begin(); o.os = os.begin();
  o.line2_start = line2_start.begin(); o.chemo_start = chemo_start.begin(); o.esr1_detect = esr1_detect.begin();
  o.d3_switch = d3_switch.begin(); o.resp1_day = resp1_day.begin(); o.d3_excl = d3_excl.begin();
  o.vaf0 = vaf0.begin(); o.vaf_det = vaf_det.begin();
  o.pfs1_ev = pfs1_ev.begin(); o.pfs2_ev = pfs2_ev.begin(); o.edct_ev = edct_ev.begin(); o.os_ev = os_ev.begin();
  o.resp1 = resp1.begin(); o.resp2 = resp2.begin(); o.in_line2 = in_line2.begin(); o.neut_g3 = neut_g3.begin();
  o.neut_g4 = neut_g4.begin(); o.stom_any = stom_any.begin(); o.stom_g2 = stom_g2.begin();
  o.stom_g3 = stom_g3.begin(); o.pneu_any = pneu_any.begin(); o.pneu_g3 = pneu_g3.begin();
  o.hyperg_g3 = hyperg_g3.begin(); o.diar_g2 = diar_g2.begin(); o.diar_g3 = diar_g3.begin();
  o.alt_g3 = alt_g3.begin(); o.qtc = qtc.begin(); o.brady = brady.begin(); o.anemia_g3 = anemia_g3.begin();
  o.esr1_at_p1 = esr1_at_p1.begin(); o.esr1_at_edct = esr1_at_edct.begin(); o.skip2 = skip2.begin();
  o.prog1 = prog1.begin(); o.ae_disc_cdk = ae_disc_cdk.begin(); o.ae_disc_eve = ae_disc_eve.begin();
  o.ae_disc_et = ae_disc_et.begin();
  o.tox_days = tox_days.begin(); o.cum_plan_cdk = cum_plan_cdk.begin(); o.cum_act_cdk = cum_act_cdk.begin();
  o.cum_plan_eve = cum_plan_eve.begin(); o.cum_act_eve = cum_act_eve.begin(); o.days_cdk = days_cdk.begin();
  o.days_eve = days_eve.begin(); o.days_gir = days_gir.begin(); o.vaf_p1 = vaf_p1.begin();
  o.vaf_edct = vaf_edct.begin(); o.best_sld1 = best_sld1.begin(); o.best_sld2 = best_sld2.begin();
  o.clones_p1 = clones_p1.begin(); o.clones_edct = clones_edct.begin(); o.N_end = N_end.begin();
  o.trN = trN.begin(); o.treve = treve.begin(); o.trcdk = trcdk.begin(); o.tranc = tranc.begin();
  o.trline = trline.begin();

  const int* pfe = fixed_end.data();
  const int* ppid = pid.begin();
  uint64_t useed = (uint64_t)(int64_t)seed;
  int n_all = rng_mode == 1 ? n_ext : n;

#pragma omp parallel for num_threads(nthreads) schedule(dynamic, 8)
  for (int i = 0; i < n; i++) {
    sim_patient(i, p, A, v, G, H, pfe, rng_mode, useed, ppid[i], extU, extZ, n_all, o);
  }

  List out = List::create(
      _["pfs1"] = pfs1, _["pfs1_ev"] = pfs1_ev, _["prog1"] = prog1, _["pfs2"] = pfs2, _["pfs2_ev"] = pfs2_ev,
      _["edct"] = edct, _["edct_ev"] = edct_ev, _["os"] = os, _["os_ev"] = os_ev,
      _["line2_start"] = line2_start, _["chemo_start"] = chemo_start, _["in_line2"] = in_line2,
      _["resp1"] = resp1, _["resp2"] = resp2, _["resp1_day"] = resp1_day,
      _["esr1_detect"] = esr1_detect, _["d3_switch"] = d3_switch, _["d3_excl"] = d3_excl);
  List tox_out = List::create(
      _["neut_g3"] = neut_g3, _["neut_g4"] = neut_g4, _["stom_any"] = stom_any, _["stom_g2"] = stom_g2,
      _["stom_g3"] = stom_g3, _["pneu_any"] = pneu_any, _["pneu_g3"] = pneu_g3, _["hyperg_g3"] = hyperg_g3,
      _["diar_g2"] = diar_g2, _["diar_g3"] = diar_g3, _["alt_g3"] = alt_g3, _["qtc"] = qtc, _["brady"] = brady,
      _["anemia_g3"] = anemia_g3, _["ae_disc_cdk"] = ae_disc_cdk, _["ae_disc_eve"] = ae_disc_eve,
      _["ae_disc_et"] = ae_disc_et, _["tox_days"] = tox_days);
  List more = List::create(
      _["cum_plan_cdk"] = cum_plan_cdk, _["cum_act_cdk"] = cum_act_cdk, _["cum_plan_eve"] = cum_plan_eve,
      _["cum_act_eve"] = cum_act_eve, _["days_cdk"] = days_cdk, _["days_eve"] = days_eve,
      _["days_gir"] = days_gir, _["esr1_at_p1"] = esr1_at_p1, _["esr1_at_edct"] = esr1_at_edct,
      _["vaf_p1"] = vaf_p1, _["vaf_edct"] = vaf_edct, _["vaf0"] = vaf0, _["vaf_det"] = vaf_det, _["best_sld1"] = best_sld1,
      _["best_sld2"] = best_sld2,
      _["skip2"] = skip2, _["clones_p1"] = clones_p1, _["clones_edct"] = clones_edct, _["N_end"] = N_end);
  List traj = R_NilValue;
  if (m > 0) {
    trN.attr("dim") = IntegerVector::create(8, m, W);       // [clone, patient, week]
    trline.attr("dim") = IntegerVector::create(m, W);
    treve.attr("dim") = IntegerVector::create(m, W);
    trcdk.attr("dim") = IntegerVector::create(m, W);
    tranc.attr("dim") = IntegerVector::create(m, W);
    traj = List::create(_["N"] = trN, _["line"] = trline, _["eve"] = treve, _["cdk"] = trcdk, _["anc"] = tranc);
  }
  return List::create(_["out"] = out, _["tox"] = tox_out, _["more"] = more, _["traj"] = traj);
}

// Uniform random numbers of the engine's counter-based generator for given patients, days 0..days-1
// and one stream (for the statistical checks in scripts/01_verify.R).
// [[Rcpp::export]]
NumericMatrix rng_sample_cpp(double seed, IntegerVector pid, int days, int stream) {
  NumericMatrix out(pid.size(), days);
  uint64_t useed = (uint64_t)(int64_t)seed;
  for (int i = 0; i < pid.size(); i++) {
    uint64_t key = patient_key(useed, (uint64_t)pid[i]);
    for (int t = 0; t < days; t++) out(i, t) = unif_ctr(key, (uint64_t)t, (uint64_t)stream);
  }
  return out;
}
