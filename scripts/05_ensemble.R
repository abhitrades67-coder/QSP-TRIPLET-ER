# Parameter-uncertainty ensemble by approximate Bayesian computation (ABC rejection with a threshold
# rule) around the calibrated point estimate.
#
# Identifiable (calibrated) parameters are perturbed around the best fit (normal, SD = 5% of the
# calibration range). Parameters the calibration data cannot inform are drawn from prior ranges:
#   kappa      (mTORC1 dependence of adaptive CDK4/6i escape)           U(0, 1)
#   IC50_e     (everolimus dose-response, relative exposure)            logU(0.08, 0.6)
#   gir_ddi    (giredestrant effect on everolimus exposure)             U(1.0, 1.4)
#   p_attr     (no endocrine 2nd line, progression after 12 months)     U(0.20, 0.40)  [ESME-MBC]
#   p_attr_mid (no endocrine 2nd line, progression at 6-12 months)      U(0.45, 0.70)  [ESME-MBC]
#   dc0_med    (post-endocrine therapy kill rate) x                     logU(0.7, 1.4)
# Writes results/ensemble.json and results/abc_candidates.csv
# Usage: Rscript scripts/05_ensemble.R [n_candidates]

if (!exists("ROOT")) ROOT <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), ".."), winslash = "/")
CAL_AS_LIBRARY <- TRUE
source(file.path(ROOT, "scripts", "03_calibrate_efficacy.R"))

args <- commandArgs(trailingOnly = TRUE)
NC <- if (length(args) >= 1) as.integer(args[1]) else 360L
NPAT <- as.integer(Sys.getenv("NPAT", "2000"))
ARMS_C <- c(calibration_arms(schedules = TRUE), persevera_arms(), paloma3_arms(), pada1_arms()["PD_cont"])

# perturbation ranges: those of the joint stage 1 + 3 fit (SERD potencies with the wider stage-3w bounds) and of stage 2
PERTURB <- do.call(rbind, list(STAGES$stage13$pars, STAGES$stage2$pars))
PERTURB <- rbind(PERTURB, data.frame(path = c("tum.h_t", "tum.h_post"), lo = c(1e-6, 1e-5), hi = c(1e-1, 1e-2),
                                     lg = c(TRUE, TRUE)))

sample_candidate <- function(p_best) {
  p <- p_best
  rec <- list()
  for (j in seq_len(nrow(PERTURB))) {
    path <- PERTURB$path[j]; lo <- PERTURB$lo[j]; hi <- PERTURB$hi[j]
    x0 <- get_param(p_best, path)
    if (PERTURB$lg[j] && x0 <= 0) { rec[[path]] <- x0; next }   # e.g. h_post = 0
    if (PERTURB$lg[j]) {
      width <- log(hi) - log(lo)
      x <- exp(min(max(log(x0) + stats::rnorm(1, 0, 0.05 * width), log(lo)), log(hi)))
    } else {
      width <- hi - lo
      x <- min(max(x0 + stats::rnorm(1, 0, 0.05 * width), lo), hi)
    }
    p <- set_param(p, path, x); rec[[path]] <- x
  }
  rec[["tum.kappa"]] <- stats::runif(1, 0, 1)
  rec[["ph.IC50_e"]] <- exp(stats::runif(1, log(0.08), log(0.6)))
  rec[["ph.gir_ddi"]] <- stats::runif(1, 1.0, 1.4)
  rec[["tum.p_attr"]] <- stats::runif(1, ATTR_PRIOR$p_attr[1], ATTR_PRIOR$p_attr[2])
  rec[["tum.p_attr_mid"]] <- stats::runif(1, ATTR_PRIOR$p_attr_mid[1], ATTR_PRIOR$p_attr_mid[2])
  rec[["tum.dc0_med"]] <- get_param(p_best, "tum.dc0_med") * exp(stats::runif(1, log(0.7), log(1.4)))
  for (k in c("tum.kappa", "ph.IC50_e", "ph.gir_ddi", "tum.p_attr", "tum.p_attr_mid", "tum.dc0_med"))
    p <- set_param(p, k, rec[[k]])
  list(p = p, rec = rec)
}

evaluate_candidate <- function(p) {
  v <- make_vpop(p, NPAT, SEED)
  grid <- make_grid(p)
  arms_needed <- sort(unique(c(STAGES$stage13$arms, STAGES$stage2$arms)))
  res <- lapply(stats::setNames(arms_needed, arms_needed), function(a) simulate(p, v, grid, ARMS_C[[a]], 1825L))
  out <- list(); zs <- c()
  for (st in c("stage1", "stage2", "stage3", "pv", "p3", "esr1")) {
    L <- loss_stage(st, res, v)
    out <- c(out, L$info); zs <- c(zs, L$loss)
  }
  R <- simulate(p, v, grid, ARMS_C$ML2_os, 2920L)
  mos <- months(km_median(R$os, R$os_ev))
  out$ML2_OS <- mos
  zos <- if (is.finite(mos)) ((mos - CAL_TARGETS$ML2_OS_cdk) / 4)^2 else 100
  E <- res$EV_ctrl
  m <- E$in_line2
  w <- evera_weights(E$esr1_at_p1[m])
  os18 <- 100 * surv_at((E$os - E$line2_start)[m], E$os_ev[m], 18 * MONTH, w)
  out$EVc_OS18 <- os18
  zos <- zos + ((os18 - CAL_TARGETS$EV_ctrl_OS18) / 5)^2
  out$loss <- sum(zs) + zos
  out$loss_s1 <- zs[1]; out$loss_s2 <- zs[2]; out$loss_s3 <- zs[3]; out$loss_pv <- zs[4]; out$loss_p3 <- zs[5]
  out$loss_esr1 <- zs[6]; out$loss_os <- zos
  out
}

p_best <- load_params(PFILES)
set.seed(2024)
cands <- lapply(seq_len(NC), function(i) sample_candidate(p_best))
t0 <- Sys.time()
rows <- vector("list", NC)
for (i in seq_len(NC)) {
  o <- tryCatch(evaluate_candidate(cands[[i]]$p), error = function(e) list(loss = 1e6, error = conditionMessage(e)))
  rows[[i]] <- c(cands[[i]]$rec, o, list(id = i - 1L))
  if (i %% 20 == 0) cat(sprintf("  %d / %d candidates (%.0f s)\n", i, NC, as.numeric(difftime(Sys.time(), t0, units = "secs"))))
}
df <- do.call(rbind, lapply(rows, function(r) as.data.frame(lapply(r, function(x) if (is.null(x)) NA else x),
                                                             check.names = FALSE, stringsAsFactors = FALSE)))
df <- df[order(df$loss), ]
utils::write.csv(df, file.path(RES, "abc_candidates.csv"), row.names = FALSE)
best_loss <- df$loss[1]
# acceptance: loss within the 95% chi-square bound for the number of targets above the best candidate, at most 60
# members. 36 targets: PALOMA-2 and BOLERO-2 6, acelERA 4, persevERA 10 with the landmarks, PALOMA-3 5, PADA-1 5,
# everolimus stage 4, survival 2
N_TARGETS <- 36L
thr <- best_loss + stats::qchisq(0.95, N_TARGETS)
acc <- utils::head(df[df$loss <= thr, ], 60)
cat(sprintf("evaluated %d candidates in %.0f s; best loss %.2f; accepted %d with loss <= %.2f\n", NC,
            as.numeric(difftime(Sys.time(), t0, units = "secs")), best_loss, nrow(acc), thr))
pcols <- grep(".", names(df), fixed = TRUE, value = TRUE)
members <- lapply(seq_len(nrow(acc)), function(i) {
  d <- list()
  for (k in pcols) {
    sk <- strsplit(k, ".", fixed = TRUE)[[1]]
    d[[sk[1]]][[sk[2]]] <- acc[[k]][i]
  }
  list(id = acc$id[i], loss = acc$loss[i], params = d)
})
save_json(list(threshold = thr, best_loss = best_loss, n_candidates = NC, members = members), "ensemble.json")
progress_note("05_ensemble", sprintf("%d candidates, %d accepted (best loss %.2f)", NC, nrow(acc), best_loss))
