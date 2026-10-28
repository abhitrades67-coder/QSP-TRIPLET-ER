# Trial-emulation metrics computed from simulation outputs.

first_line <- function(R) {
  prog <- R$prog1
  list(mPFS = months(km_median(R$pfs1, R$pfs1_ev)), ORR = mean(R$resp1),
       esr1_prog = if (any(prog)) mean(R$esr1_at_p1[prog]) else NA_real_, t = R$pfs1, e = R$pfs1_ev)
}

# line-2 PFS from the start of line 2, for patients who entered line 2
second_line <- function(R, mask = NULL) {
  m <- if (is.null(mask)) R$in_line2 else (R$in_line2 & mask)
  list(t = (R$pfs2 - R$line2_start)[m], e = R$pfs2_ev[m], m = m)
}

# evERA enrolled ESR1-mutated patients up to about 55.5% of the trial population
evera_weights <- function(esr1, target = 0.555) {
  f <- mean(esr1)
  if (!is.finite(f) || f <= 0 || f >= 1) return(rep(1, length(esr1)))
  ifelse(esr1, target / f, (1 - target) / (1 - f))
}

wmean <- function(x, w = NULL) if (is.null(w)) mean(x) else sum(x * w) / sum(w)

# Emulate evERA: line-2 entrants after >= 6 months of CDK4/6i + AI. ESR1 status = ctDNA ESR1
# detectable at line-1 progression.
evera <- function(Rexp, Rctl) {
  m <- Rexp$in_line2 & Rctl$in_line2   # identical line 1 (paired)
  esr1 <- Rexp$esr1_at_p1[m]
  out <- list(n = sum(m), esr1_frac_natural = mean(esr1))
  ex <- second_line(Rexp, m); ct <- second_line(Rctl, m)
  w <- evera_weights(esr1)
  for (lab in c("esr1", "noesr1", "all")) {
    sel <- switch(lab, esr1 = esr1, noesr1 = !esr1, all = rep(TRUE, length(esr1)))
    wa <- if (lab == "all") w[sel] else NULL
    out[[paste0("mPFS_exp_", lab)]] <- months(km_median(ex$t[sel], ex$e[sel], wa))
    out[[paste0("mPFS_ctl_", lab)]] <- months(km_median(ct$t[sel], ct$e[sel], wa))
    out[[paste0("HR_", lab)]] <- cox_hr(ex$t[sel], ex$e[sel], ct$t[sel], ct$e[sel], wa, wa)
    for (lm in c(6, 12)) {
      out[[sprintf("PFS%d_exp_%s", lm, lab)]] <- 100 * surv_at(ex$t[sel], ex$e[sel], lm * MONTH, wa)
      out[[sprintf("PFS%d_ctl_%s", lm, lab)]] <- 100 * surv_at(ct$t[sel], ct$e[sel], lm * MONTH, wa)
    }
    out[[paste0("ORR_exp_", lab)]] <- wmean(Rexp$resp2[m][sel], wa)
    out[[paste0("ORR_ctl_", lab)]] <- wmean(Rctl$resp2[m][sel], wa)
    ose <- (Rexp$os - Rexp$line2_start)[m][sel]
    osc <- (Rctl$os - Rctl$line2_start)[m][sel]
    out[[paste0("OS18_exp_", lab)]] <- 100 * surv_at(ose, Rexp$os_ev[m][sel], 18 * MONTH, wa)
    out[[paste0("OS18_ctl_", lab)]] <- 100 * surv_at(osc, Rctl$os_ev[m][sel], 18 * MONTH, wa)
  }
  out$curves <- list(te = ex$t, ee = ex$e, tc = ct$t, ec = ct$e, esr1 = esr1, w = w)
  out
}

two_arm_line2 <- function(Ra, Rb) {
  m <- Ra$in_line2 & Rb$in_line2
  a <- second_line(Ra, m); b <- second_line(Rb, m)
  list(n = sum(m), mPFS_a = months(km_median(a$t, a$e)), mPFS_b = months(km_median(b$t, b$e)),
       HR = cox_hr(a$t, a$e, b$t, b$e), esr1 = Ra$esr1_at_p1[m])
}

# acelERA emulation: a share (42%) after CDK4/6i + AI, the rest after AI alone; giredestrant vs
# fulvestrant (physician's-choice proxy)
acelera <- function(res, prior_cdk_mask) {
  tg <- eg <- tc <- ec <- es <- NULL
  for (grp in c("c", "a")) {
    mask <- if (grp == "c") prior_cdk_mask else !prior_cdk_mask
    Rg <- res[[paste0("AC_gir_", grp)]]; Rf <- res[[paste0("AC_ful_", grp)]]
    m <- Rg$in_line2 & Rf$in_line2 & mask
    a <- second_line(Rg, m); b <- second_line(Rf, m)
    tg <- c(tg, a$t); eg <- c(eg, a$e); tc <- c(tc, b$t); ec <- c(ec, b$e); es <- c(es, Rg$esr1_at_p1[m])
  }
  list(mPFS_gir = months(km_median(tg, eg)), mPFS_ctl = months(km_median(tc, ec)),
       HR = cox_hr(tg, eg, tc, ec), HR_esr1 = cox_hr(tg[es], eg[es], tc[es], ec[es]), esr1_frac = mean(es))
}

# PFS from the ctDNA-guided switch (trial emulations: the randomization, ct_lag days after the
# detection sample) to progression, for patients switched (or, in the monitoring arm, eligible for it) after
# >= min_days of first-line therapy. excl_base: leave out patients with a mutation detectable at baseline
# (PADA-1 randomized rising mutations only). Also returns 24-month PFS and OS from randomization.
serena6 <- function(Rsw, Rco, min_days = 183, excl_base = NULL) {
  d <- Rco$d3_switch
  m <- (d >= min_days) & (Rsw$d3_switch == d)
  if (!is.null(excl_base)) m <- m & !(Rco$vaf0 >= excl_base)
  ts <- (Rsw$pfs1 - d)[m]; tc <- (Rco$pfs1 - d)[m]
  es <- Rsw$pfs1_ev[m]; ec <- Rco$pfs1_ev[m]
  ok <- (ts > 0) & (tc > 0)
  os_s <- (Rsw$os - d)[m][ok]; os_c <- (Rco$os - d)[m][ok]
  list(n = sum(ok), mPFS_switch = months(km_median(ts[ok], es[ok])),
       mPFS_cont = months(km_median(tc[ok], ec[ok])), HR = cox_hr(ts[ok], es[ok], tc[ok], ec[ok]),
       PFS24_switch = 100 * surv_at(ts[ok], es[ok], 24 * MONTH), PFS24_cont = 100 * surv_at(tc[ok], ec[ok], 24 * MONTH),
       OS_HR = cox_hr(os_s, Rsw$os_ev[m][ok], os_c, Rco$os_ev[m][ok]),
       t_det = months(stats::median(Rco$esr1_detect[m][ok])), vaf_det = stats::median(Rco$vaf_det[m][ok]),
       curves = list(ts = ts[ok], es = es[ok], tc = tc[ok], ec = ec[ok]))
}

# PADA-1 step 1 (first-line AI + palbociclib with serial ctDNA; calibration targets): baseline ESR1-mutation
# prevalence; among patients with progression (or death) or a rising mutation before monitoring stopped (uniform
# 30 to 52 months after inclusion), the share with a rising mutation; the share of rising mutations with synchronous
# progression (progression within `lag` days after the detection sample, which rules out randomization); the median
# time from inclusion to randomization (detection + lag) of the others; median PFS; median allele fraction at detection.
# A mutation detectable at baseline does not count as rising (PADA-1's 3-fold-increase rule is not emulated).
pada1_step1 <- function(R, lod = attr(R, "ct_lod"), lag = 30L) {
  n <- length(R$pfs1)
  cens <- (30 + 22 * (seq_len(n) - 0.5) / n) * MONTH   # end of monitoring, uniform 30-52 months by patient index
  bpos <- R$vaf0 >= lod
  rise <- !bpos & R$esr1_detect >= 0 & R$esr1_detect <= R$pfs1
  t_rise <- ifelse(rise, R$esr1_detect, Inf)
  t_pd <- ifelse(R$pfs1_ev, R$pfs1, Inf)
  assess <- pmin(t_rise, t_pd) <= cens
  rise_obs <- rise & t_rise <= cens
  sync <- rise_obs & (t_pd <= t_rise + lag)
  rnd <- rise_obs & !sync
  # net cumulative incidence of a rising mutation (1 - Kaplan-Meier, progression without one and the end of
  # monitoring censor), as in PADA-1 (Cabel 2026, Fig 2D)
  t_net <- pmin(t_rise, t_pd, cens); e_net <- rise & t_rise <= pmin(t_pd, cens)
  list(base_prev = mean(bpos), rise_share = sum(rise_obs) / max(sum(assess), 1), rise_all = mean(rise_obs),
       sync_share = sum(sync) / max(sum(rise_obs), 1), t_rand = months(stats::median(t_rise[rnd] + lag)),
       mPFS = months(km_median(R$pfs1, R$pfs1_ev)), vaf_det = stats::median(R$vaf_det[rise_obs]),
       n_rise = sum(rise_obs), n_assess = sum(assess),
       net_inc = sapply(c(12, 24, 36), function(m) 100 * (1 - surv_at(t_net[!bpos], e_net[!bpos], m * MONTH))),
       curves = list(t_net = t_net[!bpos], e_net = e_net[!bpos]))
}
