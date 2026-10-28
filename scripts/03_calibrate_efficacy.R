# Efficacy calibration in four decoupled stages and a joint refit of stages 1 and 3.
#
# stage1: tumor growth, endocrine sensitivity, ESR1 / CDK4/6i-escape clones, CDK4/6i potency
#         targets: PALOMA-2 (both arms, HR, ORR), ESR1 prevalence after CDK4/6i, BOLERO-2 exemestane arm
# stage2: everolimus effect (I10_e, eps_g); targets: BOLERO-2, evERA control arm
# stage3: SERD potency (uG_ref, s_mut, uF_ref); targets: acelERA
# stage4: survival (post-endocrine death hazard h_post and post-endocrine therapy effect dc0_med;
#         h_t fixed); targets: MONALEESA-2 OS (ribociclib arm), evERA control-arm 18-month OS
# stage13: joint refit of stages 1 and 3, also against persevERA (HR, medians, response difference, landmark PFS),
#          PALOMA-3 and PADA-1 step 1 (ESR1 emergence), with the ESR1 spread, baseline-subclone share, PADA-1 assay
#          limit and new-lesion hazard free
# Trial emulations use each trial's published assessment schedule (R/arms.R).
# Usage: Rscript scripts/03_calibrate_efficacy.R stage1|stage2|stage3|stage4 [generations]

if (!exists("ROOT")) ROOT <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), ".."), winslash = "/")
source(file.path(ROOT, "R", "setup.R"))

NCAL <- 800L
SEED <- 21
ARMS <- c(calibration_arms(schedules = TRUE), persevera_arms(), paloma3_arms(), pada1_arms()["PD_cont"])

pars_df <- function(...) {
  x <- rbind(...)
  data.frame(path = x[, 1], lo = as.numeric(x[, 2]), hi = as.numeric(x[, 3]), lg = as.logical(x[, 4]),
             stringsAsFactors = FALSE)
}
STAGES <- list(
  stage1 = list(pars = pars_df(c("tum.kmax_med", 0.003, 0.06, TRUE), c("tum.rd", 0.40, 0.95, FALSE),
                               c("tum.b0_mean", 0.05, 0.6, FALSE), c("tum.b0_conc", 1.0, 25, TRUE),
                               c("tum.log10_fE0", -7.0, -2.0, FALSE), c("tum.log10_fR0", -10.0, -3.0, FALSE),
                               c("tum.sigma_R", 1e-11, 1e-5, TRUE), c("ph.IC50_c", 0.2, 1.5, TRUE),
                               c("tum.kmax_sd", 0.15, 0.9, FALSE),
                               # PADA-1: spread of the initial ESR1-mutant fraction, share of patients
                               # with a baseline ESR1 subclone, effective detection limit of the PADA-1 assay
                               c("tum.sd_log10_fE0", 0.3, 1.5, FALSE), c("tum.pE_base", 0.0, 0.06, FALSE),
                               c("tum.LOD_pd", 5e-4, 1e-2, TRUE),
                               # persevERA landmarks: new-lesion hazard at nadir burden
                               c("tum.lam_NL", 1e-6, 3e-4, TRUE)),
                arms = c("P2_cdk", "P2_ai", "B2_exe"), horizon = 1460L, grid_fixed = TRUE, popsize = 6),
  stage2 = list(pars = pars_df(c("ph.I10_e", 0.3, 0.95, FALSE), c("sig.eps_g", 0.0, 0.95, FALSE)),
                arms = c("B2_eve", "B2_exe", "EV_ctrl"), horizon = 1825L, grid_fixed = FALSE, popsize = 6),
  stage3 = list(pars = pars_df(c("ph.uG_ref", 30, 5000, TRUE), c("sig.s_mut", 1.0, 30.0, TRUE),
                               c("ph.uF_ref", 3, 300, TRUE)),
                arms = c("AC_gir_c", "AC_ful_c", "AC_gir_a", "AC_ful_a"), horizon = 1825L, grid_fixed = FALSE,
                popsize = 4)
)
# stage 3 with wider SERD-potency bounds (the stage-1 optimum of the R run pushed giredestrant potency
# to the lower bound of the original range)
STAGES$stage3w <- STAGES$stage3
STAGES$stage3w$pars <- pars_df(c("ph.uG_ref", 1, 5000, TRUE), c("sig.s_mut", 1.0, 300.0, TRUE), c("ph.uF_ref", 0.5, 300, TRUE))
# joint refit of stages 1 and 3: the tumor-biology parameters of stage 1 and the SERD potencies of
# stage 3 against both target sets, with persevERA (line 1 giredestrant vs letrozole, each with palbociclib)
STAGES$stage13 <- list(pars = rbind(STAGES$stage1$pars, STAGES$stage3w$pars),
                       arms = c(STAGES$stage1$arms, STAGES$stage3$arms, "PV_gir", "PV_let", "P3_pf", "P3_f", "PD_cont"),
                       horizon = 1825L,
                       grid_fixed = FALSE, popsize = 5)
# stage-1 parameters against the stage-1 and acelERA targets with SERD potencies fixed (structural
# alternatives in 07b_structural.R; the signaling grid does not change)
STAGES$stage1j <- list(pars = STAGES$stage1$pars, arms = STAGES$stage13$arms, horizon = 1825L, grid_fixed = TRUE, popsize = 6)

bounds_of <- function(stage) {
  P <- STAGES[[stage]]$pars
  tr <- function(x) ifelse(P$lg, log(pmax(x, 1e-300)), x)
  list(lower = tr(P$lo), upper = tr(P$hi))
}

apply_theta <- function(p, stage, theta) {
  P <- STAGES[[stage]]$pars
  for (j in seq_len(nrow(P))) p <- set_param(p, P$path[j], if (P$lg[j]) exp(theta[j]) else theta[j])
  p
}

zsc <- function(model, target, tol) if (is.null(model) || length(model) == 0 || !is.finite(model)) 10 else (model - target) / tol

prior_cdk_mask <- function(n) { set.seed(99); stats::runif(n) < 0.42 }   # acelERA: 41.9% after CDK4/6i

loss_stage <- function(stage, res, v) {
  T <- CAL_TARGETS
  if (stage == "stage3w") stage <- "stage3"
  if (stage %in% c("stage13", "stage1j")) {
    l1 <- loss_stage("stage1", res, v); l3 <- loss_stage("stage3", res, v); lp <- loss_stage("pv", res, v)
    l3p <- loss_stage("p3", res, v); le <- loss_stage("esr1", res, v)
    return(list(loss = l1$loss + l3$loss + lp$loss + l3p$loss + le$loss,
                info = c(l1$info, l3$info, lp$info, l3p$info, le$info,
                         list(loss_s1 = l1$loss, loss_s3 = l3$loss, loss_pv = lp$loss, loss_p3 = l3p$loss,
                              loss_esr1 = le$loss))))
  }
  if (stage == "esr1") {   # PADA-1 step 1: ESR1-mutation emergence under first-line AI + palbociclib
    s <- pada1_step1(res$PD_cont, lag = ARMS$PD_cont$ct_lag)
    P <- T$PD1
    zz <- c(zsc(s$base_prev, P$base_prev, 0.01), zsc(s$rise_share, P$rise_share, 0.04), zsc(s$sync_share, P$sync_share, 0.05),
            zsc(s$t_rand, P$t_rand, 2.0), zsc(s$mPFS, P$mPFS, 2.0))
    return(list(loss = sum(zz^2), info = list(PD_base = s$base_prev, PD_rise_share = s$rise_share, PD_sync = s$sync_share,
                                              PD_t_rand = s$t_rand, PD_mPFS = s$mPFS, PD_vaf_det = s$vaf_det)))
  }
  if (stage == "p3") {   # PALOMA-3: fulvestrant +/- palbociclib after progression on an AI
    m <- res$P3_pf$in_line2 & res$P3_f$in_line2
    r <- two_arm_line2(res$P3_pf, res$P3_f)
    o_pf <- mean(res$P3_pf$resp2[m]); o_f <- mean(res$P3_f$resp2[m])
    zz <- c(zsc(r$mPFS_a, T$P3_pf_mPFS, 1.0), zsc(r$mPFS_b, T$P3_f_mPFS, 1.0), zsc(r$HR[1], T$P3_HR, 0.05),
            zsc(o_pf, T$P3_ORR[["pf"]], 0.05), zsc(o_f, T$P3_ORR[["f"]], 0.05))
    return(list(loss = sum(zz^2), info = list(P3_pf = r$mPFS_a, P3_f = r$mPFS_b, P3_HR = r$HR[1], P3_ORR_pf = o_pf,
                                              P3_ORR_f = o_f)))
  }
  if (stage == "pv") {   # persevERA: HR, medians, the response-rate difference and landmark PFS
    g <- first_line(res$PV_gir); l <- first_line(res$PV_let)
    hr <- cox_hr(g$t, g$e, l$t, l$e)[1]
    zz <- c(zsc(hr, T$PV_HR, 0.05), zsc(g$mPFS, T$PV_gir_mPFS, 2.0), zsc(l$mPFS, T$PV_let_mPFS, 2.0),
            zsc(g$ORR - l$ORR, T$PV_ORR_diff, 0.05))
    info <- list(PV_HR = hr, PV_gir = g$mPFS, PV_let = l$mPFS, PV_ORR_diff = g$ORR - l$ORR)
    LM <- T$PV_landmarks
    if (!is.null(LM)) for (j in seq_along(LM$m)) {
      sg <- 100 * surv_at(g$t, g$e, LM$m[j] * MONTH); sl <- 100 * surv_at(l$t, l$e, LM$m[j] * MONTH)
      zz <- c(zz, zsc(sg, LM$gir[j], 4), zsc(sl, LM$let[j], 4))
      info[[sprintf("PV%d_gir", LM$m[j])]] <- sg; info[[sprintf("PV%d_let", LM$m[j])]] <- sl
    }
  } else if (stage == "stage1") {
    a <- first_line(res$P2_cdk); b <- first_line(res$P2_ai)
    hr <- cox_hr(a$t, a$e, b$t, b$e)[1]
    s2 <- second_line(res$B2_exe)
    b2 <- months(km_median(s2$t, s2$e))
    # ESR1 targets: PADA-1 (stage "esr1"); the model's ESR1 share at progression is reported only
    zz <- c(zsc(a$mPFS, T$P2_cdk_mPFS, 1.5), zsc(b$mPFS, T$P2_ai_mPFS, 1.5), zsc(hr, T$P2_HR, 0.05),
            zsc(a$ORR, T$P2_ORR_cdk, 0.05), zsc(b$ORR, T$P2_ORR_ai, 0.05), zsc(b2, T$B2_exe_mPFS, 1.0))
    info <- list(P2_cdk = a$mPFS, P2_ai = b$mPFS, HR = hr, ORR_cdk = a$ORR, ORR_ai = b$ORR,
                 ESR1_prog = a$esr1_prog, B2_exe = b2)
  } else if (stage == "stage2") {
    r <- two_arm_line2(res$B2_eve, res$B2_exe)
    s <- second_line(res$EV_ctrl)
    es <- res$EV_ctrl$esr1_at_p1[s$m]
    w <- evera_weights(es)
    ev_all <- months(km_median(s$t, s$e, w))
    ev_esr1 <- months(km_median(s$t[es], s$e[es]))
    zz <- c(zsc(r$mPFS_a, T$B2_eve_mPFS, 1.0), zsc(r$HR[1], T$B2_HR, 0.06), zsc(ev_all, T$EV_ctrl_mPFS_all, 0.7),
            zsc(ev_esr1, T$EV_ctrl_mPFS_esr1, 0.9))
    info <- list(B2_eve = r$mPFS_a, B2_exe = r$mPFS_b, B2_HR = r$HR[1], EV_ctrl_all = ev_all, EV_ctrl_esr1 = ev_esr1)
  } else {
    r <- acelera(res, prior_cdk_mask(v$n))
    zz <- c(zsc(r$mPFS_gir, T$AC_gir_mPFS, 0.7), zsc(r$mPFS_ctl, T$AC_ctl_mPFS, 0.7), zsc(r$HR[1], T$AC_HR, 0.06),
            zsc(r$HR_esr1[1], T$AC_HR_esr1, 0.08))
    info <- list(AC_gir = r$mPFS_gir, AC_ctl = r$mPFS_ctl, AC_HR = r$HR[1], AC_HR_esr1 = r$HR_esr1[1])
  }
  list(loss = sum(zz^2), info = info)
}

evaluate_stage <- function(stage, theta, p0, v0, grid0 = NULL, return_res = FALSE) {
  p <- apply_theta(p0, stage, theta)
  v <- refresh_derived(p, v0)
  cfg <- STAGES[[stage]]
  grid <- if (cfg$grid_fixed) grid0 else make_grid(p)
  res <- lapply(stats::setNames(cfg$arms, cfg$arms), function(a) simulate(p, v, grid, ARMS[[a]], cfg$horizon))
  L <- loss_stage(stage, res, v)
  if (return_res) c(L, list(res = res)) else L$loss
}

run_de <- function(stage, p, maxiter, popsize = STAGES[[stage]]$popsize, tag = "", p_extra = NULL) {
  v <- make_vpop(p, NCAL, SEED)
  grid <- if (STAGES[[stage]]$grid_fixed) make_grid(p) else NULL
  P <- STAGES[[stage]]$pars
  b <- bounds_of(stage)
  x0 <- sapply(seq_len(nrow(P)), function(j) {
    val <- get_param(p, P$path[j])
    x <- if (P$lg[j]) log(val) else val
    min(max(x, b$lower[j] + 1e-9), b$upper[j] - 1e-9)
  })
  if (!is.null(p_extra)) {   # additional starting points (parameter sets)
    x0 <- rbind(x0, do.call(rbind, lapply(p_extra, function(pe) sapply(seq_len(nrow(P)), function(j) {
      val <- get_param(pe, P$path[j]); x <- if (P$lg[j]) log(val) else val
      min(max(x, b$lower[j] + 1e-9), b$upper[j] - 1e-9) }))))
  }
  fn <- function(theta) tryCatch(evaluate_stage(stage, theta, p, v, grid), error = function(e) {
    cat("eval error:", conditionMessage(e), "\n"); 1e6 })
  t0 <- Sys.time()
  out <- de_optimize(fn, b$lower, b$upper, x0 = x0, popsize = popsize, maxiter = maxiter, tol = 1e-3, seed = 3)
  fin <- evaluate_stage(stage, out$par, p, v, grid, return_res = TRUE)
  cat(stage, "loss", fin$loss, "\n"); print(unlist(fin$info))
  pars <- list()
  for (j in seq_len(nrow(P))) {
    sk <- strsplit(P$path[j], ".", fixed = TRUE)[[1]]
    pars[[sk[1]]][[sk[2]]] <- if (P$lg[j]) exp(out$par[j]) else out$par[j]
  }
  save_json(pars, sprintf("params_%s%s.json", stage, tag))
  save_json(list(loss = fin$loss, info = fin$info, x = out$par, nfev = out$nfev, generations = nrow(out$history) - 1,
                 seconds = as.numeric(difftime(Sys.time(), t0, units = "secs"))),
            sprintf("calib_%s%s.json", stage, tag))
  utils::write.csv(out$history, file.path(RES, sprintf("calib_%s%s_history.csv", stage, tag)), row.names = FALSE)
  pars
}

# ---------------------------------------------------------------------------- stage 4
s4_eval <- function(p, v, grid, hp, dc, ht = NA) {
  p$tum$h_post <- hp
  p$tum$dc0_med <- dc
  if (!is.na(ht)) p$tum$h_t <- ht
  v <- refresh_derived(p, v)
  R <- simulate(p, v, grid, ARMS$ML2_os, 2920L)
  mos <- months(km_median(R$os, R$os_ev))
  E <- simulate(p, v, grid, ARMS$EV_ctrl, 2190L)
  m <- E$in_line2
  w <- evera_weights(E$esr1_at_p1[m])
  os18 <- 100 * surv_at((E$os - E$line2_start)[m], E$os_ev[m], 18 * MONTH, w)
  L <- if (is.finite(mos)) ((mos - CAL_TARGETS$ML2_OS_cdk) / 4)^2 + ((os18 - CAL_TARGETS$EV_ctrl_OS18) / 5)^2 else 1e3
  data.frame(h_post = hp, dc0_med = dc, h_t = ht, ML2_OS = mos, EVc_OS18 = os18, loss = L)
}

stage4 <- function(p) {
  v <- make_vpop(p, 1500, SEED)
  grid <- make_grid(p)
  p$tum$h_t <- TUMOR$h_t
  t0 <- Sys.time()
  hps <- c(0.0, 2e-4, 4e-4, 7e-4, 1e-3, 1.5e-3, 2.2e-3)
  dcs <- c(0.03, 0.04, 0.05, 0.065, 0.085, 0.11, 0.15, 0.20)
  coarse <- do.call(rbind, lapply(hps, function(hp) do.call(rbind, lapply(dcs, function(dc) s4_eval(p, v, grid, hp, dc)))))
  ceiling_ <- s4_eval(p, v, grid, 0.0, p$tum$dc0_med, 1e-5)   # tumor-burden hazard alone at its floor
  print(coarse, digits = 4)
  allg <- coarse
  b <- allg[which.min(allg$loss), ]
  cat(sprintf("coarse grid (%.0f s): best h_post %.2e, dc0_med %.3f, loss %.3f\n",
              as.numeric(difftime(Sys.time(), t0, units = "secs")), b$h_post, b$dc0_med, b$loss))
  # refinement around the best point, repeated while the best point lies on the edge of the refinement grid
  for (round in 1:4) {
    hp0 <- max(b$h_post, 1e-4)
    hp_f <- exp(seq(log(hp0 / 1.6), log(hp0 * 1.6), length.out = 5))
    dc_f <- exp(seq(log(b$dc0_med / 1.3), log(b$dc0_med * 1.3), length.out = 5))
    fine <- do.call(rbind, lapply(hp_f, function(hp) do.call(rbind, lapply(dc_f, function(dc) s4_eval(p, v, grid, hp, dc)))))
    allg <- rbind(allg, fine)
    b <- allg[which.min(allg$loss), ]
    edge <- isTRUE(all.equal(b$dc0_med, max(dc_f))) || isTRUE(all.equal(b$dc0_med, min(dc_f))) ||
      isTRUE(all.equal(b$h_post, max(hp_f))) || isTRUE(all.equal(b$h_post, min(hp_f)))
    cat(sprintf("refinement %d: best h_post %.3e, dc0_med %.4f, ML2 OS %.1f, evERA ctl OS18 %.1f, loss %.3f%s\n",
                round, b$h_post, b$dc0_med, b$ML2_OS, b$EVc_OS18, b$loss, if (edge) " (edge: refine again)" else ""))
    if (!edge) break
  }
  utils::write.csv(allg, file.path(RES, "stage4_grid.csv"), row.names = FALSE)
  base <- coarse[coarse$h_post == 0, ]
  base <- base[which.min(base$loss), ]
  save_json(list(tum = list(h_t = p$tum$h_t, h_post = b$h_post, dc0_med = b$dc0_med)), "params_stage4.json")
  save_json(list(ML2_OS = b$ML2_OS, EVc_OS18 = b$EVc_OS18, loss = b$loss, h_t = p$tum$h_t, h_post = b$h_post,
                 dc0_med = b$dc0_med, h_t_only_ceiling_OS = ceiling_$ML2_OS, no_h_post_best = as.list(base),
                 seconds = as.numeric(difftime(Sys.time(), t0, units = "secs"))), "calib_stage4.json")
}

# run as a script (other scripts set CAL_AS_LIBRARY before sourcing this file to reuse its functions)
if (!exists("CAL_AS_LIBRARY")) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) >= 1 && args[1] %in% c("stage1", "stage2", "stage3", "stage3w", "stage13", "stage4")) {
    stage <- args[1]
    p <- load_params(PFILES)
    if (stage == "stage4") {
      stage4(p)
    } else {
      gens <- if (length(args) > 1) as.integer(args[2]) else c(stage1 = 25, stage2 = 15, stage3 = 12, stage3w = 12, stage13 = 25)[[stage]]
      # joint fit: start from the current parameter files (calibration_start/, copied by 03a_start_state.R) and,
      # as further members of the initial population, from the optimum of the SERD-potency scan (giredestrant 11,
      # mutant-ER shift 4, fulvestrant 7) with three settings of the ESR1 parameters.
      p_extra <- NULL
      if (stage == "stage13") {
        b <- function(upd) { q <- p; for (k in names(upd)) q <- set_param(q, k, upd[[k]]); q }
        serd <- list("ph.uG_ref" = 11, "sig.s_mut" = 4, "ph.uF_ref" = 7)
        s_c <- c(serd, list("ph.uG_ref" = 9, "tum.LOD_pd" = 0.002, "tum.pE_base" = 0.01, "tum.log10_fE0" = -3.3))
        p_extra <- list(b(serd),
                        b(c(serd, list("tum.LOD_pd" = 0.001, "tum.pE_base" = 0.01, "tum.sd_log10_fE0" = 0.6))),
                        b(c(serd, list("tum.LOD_pd" = 0.003, "tum.pE_base" = 0.03, "tum.log10_fE0" = -3.5))),
                        b(s_c),
                        # long first-line PFS tail (diagnostics: lower new-lesion hazard and adaptive escape)
                        b(c(s_c, list("tum.lam_NL" = 1e-5, "tum.sigma_R" = 1e-10))),
                        b(c(s_c, list("tum.lam_NL" = 3e-5, "tum.sigma_R" = 1e-10, "tum.b0_conc" = 2.5))))
      }
      pars <- run_de(stage, p, maxiter = gens, p_extra = p_extra)
      # the joint fit replaces the stage-1 and stage-3 parameter files
      if (stage == "stage13") {
        for (f in c("params_stage1", "calib_stage1", "params_stage3", "calib_stage3")) {   # keep the staged fit
          src <- file.path(RES, paste0(f, ".json")); dst <- file.path(RES, paste0(f, "_staged.json"))
          if (file.exists(src) && !file.exists(dst)) file.copy(src, dst)
        }
        P1 <- STAGES$stage1$pars$path; P3 <- STAGES$stage3w$pars$path
        split_pars <- function(paths) { out <- list(); for (pp in paths) { sk <- strsplit(pp, ".", fixed = TRUE)[[1]]
          out[[sk[1]]][[sk[2]]] <- pars[[sk[1]]][[sk[2]]] }; out }
        save_json(split_pars(P1), "params_stage1.json")
        save_json(split_pars(P3), "params_stage3.json")
      }
      if (stage == "stage3w") file.copy(file.path(RES, "params_stage3w.json"), file.path(RES, "params_stage3.json"), overwrite = TRUE)
    }
    progress_note(paste0("03_calibrate_efficacy ", stage), paste("done;", paste(args, collapse = " ")))
  }
}
