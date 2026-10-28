# Population simulation engine: R interface to the C++ engine (R/engine.cpp).
#
# Daily time step. For each virtual patient the engine tracks:
#   * treatment line and regimen (doses, holds, discontinuations)
#   * PK (relative concentrations) and CYP3A4 mechanism-based inactivation (everolimus DDI)
#   * neutrophils (Friberg model, 0.5-day sub-steps) and event-based toxicities
#   * 8 tumor clones (ESR1 / PI3K / CDK4/6i-escape) with proliferation from the signaling grid
#   * RECIST-like assessments, new lesions, ctDNA ESR1 variant allele fraction
#   * death (tumor-burden-dependent + background + post-endocrine) and post-endocrine therapy
# Random numbers are a function of (population seed, patient, day, stream), so the same patient
# sees the same random stream under every arm (common random numbers).

regimen <- function(cdk = NULL, et = NULL, eve = NULL, cdk_doses = NULL, cdk_onoff = NULL, label = "") {
  list(cdk = cdk, et = et, eve = eve, cdk_doses = cdk_doses, cdk_onoff = cdk_onoff, label = label)
}

# regimen -> numeric description for the C++ engine
compile_regimen <- function(reg, p) {
  if (is.null(reg$cdk)) {
    cdk_code <- 0L; doses <- c(-1, -1, -1); ref <- 1; onoff <- 0L; k <- 0.5; slope <- 0; ei <- 0L
  } else {
    base <- if (startsWith(reg$cdk, "ribo")) "ribo" else reg$cdk
    d <- if (is.null(reg$cdk_doses)) CDK_DOSES[[base]] else reg$cdk_doses
    doses <- c(d, rep(-1, 3 - length(d)))
    cdk_code <- CDK_CODE[[base]]
    ref <- CDK_REF[[base]]
    onoff <- as.integer(if (is.null(reg$cdk_onoff)) CDK_ONOFF[[base]] else reg$cdk_onoff)
    k <- p$ph[[paste0("k_", base)]]
    slope <- p$tox[[paste0("slope_", base)]]
    ei <- match(base, c("palbo", "ribo", "abema")) - 1L
  }
  et <- if (is.null(reg$et)) 0L else ET_CODE[[reg$et]]
  if (is.null(reg$eve)) {
    ed <- c(-1, -1, -1)
  } else {
    e <- EVE_DOSES[[format(reg$eve)]]
    ed <- c(e, rep(-1, 3 - length(e)))
  }
  list(cdk_code = as.integer(cdk_code), cdk_doses = as.numeric(doses), cdk_ref = as.numeric(ref),
       cdk_onoff = as.integer(onoff), k_cdk = as.numeric(k), slope_n = as.numeric(slope),
       expo_idx = as.integer(ei), et_code = as.integer(et), eve_doses = as.numeric(ed))
}

ENGINE_KEEP_INT <- c("pfs1", "pfs2", "edct", "os", "line2_start", "chemo_start", "esr1_detect", "d3_excl",
                     "d3_switch", "resp1_day")

# Simulate one arm. arm = list(lines = list(regimen, [regimen]), optional: p_attr (attrition after
# > attr_switch_days on line 1), p_attr_mid (attrition between early_prog_days and attr_switch_days),
# early_prog_days, elig_min_days, d3_switch, d3_monitor, fixed_duration (toxicity-only runs), tox_over (named list of
# toxicity-parameter overrides), ass1 / ass2 (tumor assessment schedule per line:
# c(first scan, interval, switch day, later interval), days); ct_first / ct_int (ctDNA sampling, days
# after the line starts; default 56 / 56), ct_lod (detection limit; default tum.LOD), switch_et (endocrine therapy
# after a ctDNA-guided switch: "giredestrant" (default) or "fulvestrant"), ct_lag (days from detection to the switch;
# default 0), ct_confirm (imaging check at detection and on the switch day)).
# rng: NULL (counter-based generator seeded by v$seed_events) or list(U =, Z =, n =, pid =) to
# use externally supplied random streams (for testing).
simulate <- function(p, v, grid, arm, horizon = 2190L, record_traj = 0L, rng = NULL,
                     nthreads = qsp_threads()) {
  tox <- p$tox
  if (!is.null(arm$tox_over)) for (k in names(arm$tox_over)) tox[[k]] <- arm$tox_over[[k]]
  carm <- list(lines = lapply(arm$lines, compile_regimen, p = p))
  for (k in c("p_attr", "p_attr_mid", "attr_switch_days", "early_prog_days", "elig_min_days", "ct_first", "ct_int", "ct_lag",
              "ct_lod"))
    if (!is.null(arm[[k]])) carm[[k]] <- as.numeric(arm[[k]])
  if (!is.null(arm$switch_et)) carm$switch_et <- as.numeric(if (is.character(arm$switch_et)) ET_CODE[[arm$switch_et]] else arm$switch_et)
  # a detection limit taken from a (calibrated) tumor parameter, optionally scaled (the effective limit
  # of the PADA-1 assay, tum.LOD_pd; SERENA-6 uses it scaled by the ratio of the allele fractions at detection)
  if (!is.null(arm$ct_lod_par)) carm$ct_lod <- as.numeric(p$tum[[arm$ct_lod_par]]) * (if (is.null(arm$ct_lod_mult)) 1 else arm$ct_lod_mult)
  for (k in c("d3_switch", "d3_monitor", "ct_confirm")) if (!is.null(arm[[k]])) carm[[k]] <- isTRUE(arm[[k]])
  for (k in c("ass1", "ass2")) if (!is.null(arm[[k]])) carm[[k]] <- as.integer(arm[[k]])
  if (!is.null(arm$fixed_duration)) carm$fixed_duration <- as.integer(arm$fixed_duration)
  vp <- v[VPOP_FIELDS]
  vp$Ptrunc <- as.integer(vp$Ptrunc)
  if (is.null(grid)) grid <- list(P = numeric(0), S6K = numeric(0), P0_ref = 1)
  if (is.null(rng)) {
    res <- simulate_cpp(p$tum, p$ph, tox, p$st, vp, grid, carm, as.integer(horizon),
                        as.integer(record_traj), 0L, as.numeric(v$seed_events), NULL, NULL, NULL, 0L,
                        as.integer(nthreads))
  } else {
    pid <- if (is.null(rng$pid)) seq_len(v$n) - 1L else as.integer(rng$pid)
    res <- simulate_cpp(p$tum, p$ph, tox, p$st, vp, grid, carm, as.integer(horizon),
                        as.integer(record_traj), 1L, 0, pid, rng$U, rng$Z, as.integer(rng$n),
                        as.integer(nthreads))
  }
  R <- c(res$out, res$tox, res$more)
  attr(R, "ct_lod") <- if (!is.null(carm$ct_lod)) carm$ct_lod else p$tum$LOD   # detection limit used
  if (!is.null(res$traj)) R$traj <- res$traj
  R
}

# ---------------------------------------------------------------------------- PK / DDI helpers
# everolimus AUC ratio at steady state under continuous CDK4/6i dosing at relative concentration c
ddi_aucr <- function(ph, agent, c = 1.0) {
  kin <- ph[[paste0("kinact_", agent)]] * c / (ph$KI_3A4 + c)
  E <- ph$kdeg_3A4 / (ph$kdeg_3A4 + kin)
  1 / (ph$fm_eve * E + 1 - ph$fm_eve)
}

# kinact giving a target AUC ratio at relative concentration c (inverse of ddi_aucr)
kinact_for_aucr <- function(ph, aucr, c = 1.0) {
  E <- (1 / aucr - (1 - ph$fm_eve)) / ph$fm_eve
  kin <- ph$kdeg_3A4 * (1 / E - 1)
  kin * (ph$KI_3A4 + c) / c
}

# Daily PK / DDI profile for a typical patient (relative exposure 1), as in the engine: CDK4/6i
# concentration, CYP3A4 activity and everolimus exposure (relative to 10 mg alone).
pk_profile <- function(ph, agent = NULL, cdk_rel = 1, eve_rel = 1, days = 84, onoff = NULL) {
  kin_par <- if (is.null(agent)) 0 else switch(agent, palbo = ph$kinact_palbo, ribo = ph$kinact_ribo, 0)
  k_c <- if (is.null(agent)) 0.5 else ph[[paste0("k_", agent)]]
  if (is.null(onoff)) onoff <- !is.null(agent) && agent %in% c("palbo", "ribo")
  c_c <- 0; c_e <- 0; cyp <- 1
  out <- data.frame(day = 0:(days - 1), cdk = 0, eve = 0, cyp = 0)
  for (d in 0:(days - 1)) {
    on <- if (onoff) (d %% 28) < 21 else TRUE
    rate <- if (!is.null(agent) && on) cdk_rel else 0
    kin <- kin_par * c_c / (ph$KI_3A4 + c_c)
    cyp <- (cyp + ph$kdeg_3A4) / (1 + ph$kdeg_3A4 + kin)
    aucr <- 1 / (ph$fm_eve * cyp + 1 - ph$fm_eve)
    ek <- exp(-k_c); c_c <- c_c * ek + rate * (1 - ek)
    ee <- exp(-ph$k_eve); c_e <- c_e * ee + eve_rel * aucr * (1 - ee)
    out[d + 1, 2:4] <- c(c_c, c_e, cyp)
  }
  out
}
