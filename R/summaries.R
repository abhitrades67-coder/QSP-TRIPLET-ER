# Per-arm endpoint summaries for the in silico trial. There is no
# censoring before the 6-year horizon, so restricted means are exact sample means.

KEEP <- c("pfs1", "pfs1_ev", "prog1", "edct", "edct_ev", "os", "os_ev", "in_line2", "line2_start",
          "pfs2", "pfs2_ev", "resp1", "resp2", "neut_g3", "neut_g4", "stom_any", "stom_g2", "stom_g3",
          "pneu_any", "pneu_g3", "hyperg_g3", "diar_g2", "diar_g3", "alt_g3", "qtc", "brady",
          "anemia_g3", "ae_disc_cdk", "ae_disc_eve", "ae_disc_et", "tox_days", "cum_plan_cdk",
          "cum_act_cdk", "cum_plan_eve", "cum_act_eve", "esr1_detect", "esr1_at_p1", "esr1_at_edct",
          "clones_p1", "clones_edct", "d3_switch", "best_sld1", "skip2", "chemo_start")
AE_KEYS <- c("neut_g3", "neut_g4", "stom_any", "stom_g2", "stom_g3", "pneu_any", "pneu_g3", "hyperg_g3",
             "diar_g2", "diar_g3", "alt_g3", "qtc", "brady", "anemia_g3", "ae_disc_cdk", "ae_disc_eve",
             "ae_disc_et")

rmean <- function(t, tau_m) mean(pmin(t, tau_m * MONTH)) / MONTH

qtwist <- function(R, u_tox = 0.5, u_rel = 0.5, tau_m = 60) {
  tau <- tau_m * MONTH
  E <- pmin(R$edct, tau); O <- pmin(R$os, tau)
  TOX <- pmin(R$tox_days, E)
  mean(u_tox * TOX + (E - TOX) + u_rel * (O - E)) / MONTH
}

subset_R <- function(R, m) {
  out <- list()
  for (k in intersect(KEEP, names(R))) {
    x <- R[[k]]
    out[[k]] <- if (is.matrix(x)) x[m, , drop = FALSE] else x[m]
  }
  out
}

summarize_arm <- function(R, mask = NULL) {
  r <- if (is.null(mask)) R else subset_R(R, mask)
  s <- list()
  s$EDCT_RMST48 <- rmean(r$edct, 48)
  s$EDCT_RMST36 <- rmean(r$edct, 36)
  s$EDCT_RMST60 <- rmean(r$edct, 60)
  s$EDCT_median <- months(km_median(r$edct, r$edct_ev))
  s$OS_RMST60 <- rmean(r$os, 60)
  s$OS_RMST72 <- rmean(r$os, 72)
  s$OS_median <- months(km_median(r$os, r$os_ev))
  s$OS_36 <- mean(r$os > 36 * MONTH)
  s$OS_60 <- mean(r$os > 60 * MONTH)
  s$PFS1_median <- months(km_median(r$pfs1, r$pfs1_ev))
  s$PFS1_RMST48 <- rmean(r$pfs1, 48)
  s$ORR1 <- mean(r$resp1)
  l2 <- r$in_line2
  s$frac_line2 <- mean(l2)
  s$PFS_line2_median <- if (any(l2)) months(km_median((r$pfs2 - r$line2_start)[l2], r$pfs2_ev[l2])) else NA_real_
  for (k in AE_KEYS) s[[k]] <- mean(r[[k]])
  s$ae_disc_any <- mean(r$ae_disc_cdk | r$ae_disc_eve | r$ae_disc_et)
  s$sympt_g3 <- mean(r$stom_g3 | r$pneu_g3 | r$hyperg_g3 | r$diar_g3 | r$anemia_g3)
  pc <- r$cum_plan_cdk > 0
  pe <- r$cum_plan_eve > 0
  s$RDI_cdk <- if (any(pc)) mean(r$cum_act_cdk[pc] / r$cum_plan_cdk[pc]) else NA_real_
  s$RDI_eve <- if (any(pe)) mean(r$cum_act_eve[pe] / r$cum_plan_eve[pe]) else NA_real_
  s$tox_months <- mean(r$tox_days) / MONTH
  s$QTWiST <- qtwist(r)
  det <- r$esr1_detect
  s$ESR1_detect_24m <- mean(det >= 0 & det <= 24 * MONTH)
  fail <- r$edct_ev & (r$os > r$edct)
  s$ESR1_at_EDCT <- if (any(fail)) mean(r$esr1_at_edct[fail]) else NA_real_
  if (any(fail)) {
    ce <- r$clones_edct[fail, , drop = FALSE]
    s$frac_E_at_EDCT <- mean(rowSums(ce[, c(2, 4, 6, 8), drop = FALSE]))
    s$frac_R_at_EDCT <- mean(rowSums(ce[, c(5, 6, 7, 8), drop = FALSE]))
  } else {
    s$frac_E_at_EDCT <- NA_real_
    s$frac_R_at_EDCT <- NA_real_
  }
  s
}
