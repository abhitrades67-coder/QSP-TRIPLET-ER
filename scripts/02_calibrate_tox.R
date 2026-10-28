# Stage 0: calibrate toxicity hazards to published incidence (toxicity-only simulations with fixed
# treatment durations drawn around the published median duration of each regimen).
# Targets and durations: R/inputs.R. Discontinuation for other AEs (R/params.R): the hazard falls from the start of
# each line, everolimus discontinuation is proportional to exposure and endocrine therapy has an extra hazard with
# everolimus; it is fitted to short- and long-exposure
# published rates: hazard levels by bisection on their primary targets, decay shape by grid search and
# Nelder-Mead on all ten targets.
# Writes results/params_tox.json, results/tox_calibration.csv, results/disc_calibration.csv and
# results/disc_shape_grid.csv

if (!exists("ROOT")) ROOT <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), ".."), winslash = "/")
source(file.path(ROOT, "R", "setup.R"))
N <- as.integer(Sys.getenv("NTOX", "4000"))

p <- default_params()
v <- make_vpop(p, N, seed = 11)

# durations: exponential around the published median, capped (default 3 years); the same uniforms for
# every regimen (common random numbers)
set.seed(77)
U_DUR <- stats::runif(N)
durations <- function(med, cap = 1095L) pmin(as.integer(-med / log(2) * log(U_DUR)) + 1L, as.integer(cap))

ARMS <- list(
  palbo = list(regimen("palbo", "letrozole"), TOX_DURATIONS[["palbo"]]),
  ribo = list(regimen("ribo", "letrozole"), TOX_DURATIONS[["ribo"]]),
  abema = list(regimen("abema", "letrozole"), TOX_DURATIONS[["abema"]]),
  evectl = list(regimen(NULL, "exemestane", 10), TOX_DURATIONS[["evectl"]]),
  evexp = list(regimen(NULL, "giredestrant", 10), TOX_DURATIONS[["evexp"]])
)

# target key, regimen, outcome, parameter, kind (event toxicities; discontinuation is fitted separately)
MAP <- read.table(header = TRUE, stringsAsFactors = FALSE, text = "
key            arm     outcome      par          kind
neut_g3_palbo  palbo   neut_g3      slope_palbo  haz
neut_g3_ribo   ribo    neut_g3      slope_ribo   haz
neut_g3_abema  abema   neut_g3      slope_abema  haz
stom_any_eve   evectl  stom_any     aS_eve       haz
stom_g3_eve    evectl  stom_g3      p3S          prob
pneu_any_eve   evectl  pneu_any     aP_eve       haz
hyperg_g3_eve  evectl  hyperg_g3    aG_eve       haz
diar_g3_eve    evectl  diar_g3      aD_eve       haz
anemia_g3_eve  evectl  anemia_g3    aA_eve       haz
diar_g2_abema  abema   diar_g2      aD_abema     haz
diar_g3_abema  abema   diar_g3      q3D_abema    prob
alt_g3_ribo    ribo    alt_g3       aL_ribo      haz
qtc_ribo       ribo    qtc          aQ_ribo      haz
brady_gir      evexp   brady        aB_gir       haz
")
DISC_PARS <- c("disc_cdk", "disc_eve", "disc_et", "disc_gir", "disc_et_eve", "disc_abema_mult", "disc_tau", "disc_floor", "disc_t0",
               "disc_eve_gamma")

run_arm <- function(p, k) {
  a <- ARMS[[k]]
  simulate(p, v, NULL, list(lines = list(a[[1]]), fixed_duration = durations(a[[2]])), 1100L)
}

# monotone bisection on log(parameter)
bisect <- function(p, key, arm, outc, par, kind, n_it) {
  tgt <- TOX_TARGETS[[key]]
  x0 <- p$tox[[par]]
  lo <- log(x0) - 3; hi <- log(x0) + 3
  if (kind == "prob") hi <- min(hi, log(0.999))
  for (it in seq_len(n_it)) {
    mid <- 0.5 * (lo + hi)
    p$tox[[par]] <- exp(mid)
    obs <- mean(run_arm(p, arm)[[outc]])
    if (obs < tgt) lo <- mid else hi <- mid
  }
  p$tox[[par]] <- exp(0.5 * (lo + hi))
  p
}

# ---------------------------------------------------------------------------- discontinuation sub-model
run_disc <- function(p, k) {
  a <- DISC_ARMS[[k]]
  cap <- as.integer(round(a$cap * MONTH))
  arm <- list(lines = list(regimen(a$cdk, a$et, a$eve)), fixed_duration = durations(a$med * MONTH, cap))
  if (!is.null(a$over)) arm$tox_over <- a$over
  r <- simulate(p, v, NULL, arm, cap + 5L)
  c(cdk = mean(r$ae_disc_cdk), eve = mean(r$ae_disc_eve), et = mean(r$ae_disc_et),
    any = mean(r$ae_disc_cdk | r$ae_disc_eve | r$ae_disc_et))
}
disc_all <- function(p) {
  arms <- unique(DISC_TARGETS$arm)
  out <- lapply(stats::setNames(arms, arms), function(k) run_disc(p, k))
  vapply(seq_len(nrow(DISC_TARGETS)), function(j) out[[DISC_TARGETS$arm[j]]][[DISC_TARGETS$outcome[j]]], numeric(1))
}
# deviation in standard errors of the published proportion
disc_z <- function(m) (m - DISC_TARGETS$target) / sqrt(DISC_TARGETS$target * (1 - DISC_TARGETS$target) / DISC_TARGETS$n)

# hazard levels for a given decay shape: bisection of each level on its primary target, repeated so that
# the coupled levels (endocrine hazard with and without everolimus) settle
fit_levels <- function(p, rounds = 2, width = 2, n_it = 9) {
  for (r in seq_len(rounds)) for (j in seq_len(nrow(DISC_LEVELS))) {
    par <- DISC_LEVELS$par[j]
    tg <- DISC_TARGETS[DISC_TARGETS$key == DISC_LEVELS$key[j], ]
    x0 <- log(max(p$tox[[par]], 1e-8)); lo <- x0 - width; hi <- x0 + width
    for (it in seq_len(n_it)) {
      mid <- 0.5 * (lo + hi)
      p$tox[[par]] <- exp(mid)
      obs <- run_disc(p, tg$arm)[[tg$outcome]]
      if (obs < tg$target) lo <- mid else hi <- mid
    }
    p$tox[[par]] <- exp(0.5 * (lo + hi))
  }
  p
}
shape_eval <- function(p, tau, fl, width = 2, n_it = 9) {
  p$tox$disc_tau <- tau; p$tox$disc_floor <- fl
  p <- fit_levels(p, rounds = 2, width = width, n_it = n_it)
  m <- disc_all(p)
  list(p = p, loss = sum(disc_z(m)^2), model = m)
}
fit_shape <- function(p, first) {
  grid_rows <- list()
  if (first) {   # grid over the decay shape, plus a constant hazard (floor 1)
    G <- rbind(expand.grid(tau = c(45, 90, 180, 365, 730), fl = c(0.05, 0.15, 0.30, 0.50)), data.frame(tau = 180, fl = 1.0))
    best <- NULL
    for (i in seq_len(nrow(G))) {
      e <- shape_eval(p, G$tau[i], G$fl[i])
      grid_rows[[i]] <- data.frame(tau = G$tau[i], floor = G$fl[i], loss = e$loss, t(stats::setNames(e$model, DISC_TARGETS$key)))
      cat(sprintf("  shape tau %4.0f d, floor %.2f: loss %7.2f\n", G$tau[i], G$fl[i], e$loss))
      # constant hazards refitted to the same targets: structural alternative for 06b
      if (G$fl[i] >= 1) save_json(list(tox = e$p$tox[DISC_PARS], loss = e$loss, model = as.list(stats::setNames(e$model, DISC_TARGETS$key))),
                                  "disc_const.json")
      if (is.null(best) || e$loss < best$loss) best <- e
    }
    utils::write.csv(do.call(rbind, grid_rows), file.path(RES, "disc_shape_grid.csv"), row.names = FALSE)
    p <- best$p
    if (p$tox$disc_floor >= 1) return(best)   # the constant hazard fits best: nothing to refine
  }
  # local refinement (Nelder-Mead on log tau and logit floor), levels re-fitted at every point
  cur <- list(p = p, loss = Inf)
  fn <- function(x) {
    e <- shape_eval(cur$p, exp(x[1]), stats::plogis(x[2]), width = 1, n_it = 7)
    if (e$loss < cur$loss) cur <<- e
    e$loss
  }
  stats::optim(c(log(p$tox$disc_tau), stats::qlogis(p$tox$disc_floor)), fn, method = "Nelder-Mead",
               control = list(maxit = if (first) 24 else 12, reltol = 1e-3))
  cur
}

t0 <- Sys.time()
disc_loss_prev <- Inf
for (pass in 0:4) {
  for (j in seq_len(nrow(MAP))) {
    m <- MAP[j, ]
    p <- bisect(p, m$key, m$arm, m$outcome, m$par, m$kind, if (pass == 0) 12 else 8)
  }
  fs <- fit_shape(p, first = pass == 0)
  p <- fs$p
  res <- lapply(stats::setNames(names(ARMS), names(ARMS)), function(k) run_arm(p, k))
  err <- sapply(seq_len(nrow(MAP)), function(j) abs(mean(res[[MAP$arm[j]]][[MAP$outcome[j]]]) - TOX_TARGETS[[MAP$key[j]]]))
  dz <- disc_z(fs$model)
  cat(sprintf("pass %d: event toxicities max abs error %.4f; discontinuation loss %.2f (max |z| %.2f), tau %.0f d, floor %.3f (%.0f s)\n",
              pass, max(err), fs$loss, max(abs(dz)), p$tox$disc_tau, p$tox$disc_floor,
              as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  if (pass >= 1 && max(err) < 0.004 && abs(disc_loss_prev - fs$loss) < 0.5) break
  disc_loss_prev <- fs$loss
}

rows <- data.frame(target = MAP$key, regimen = MAP$arm, observed_target = unlist(TOX_TARGETS[MAP$key]),
                   model = sapply(seq_len(nrow(MAP)), function(j) mean(res[[MAP$arm[j]]][[MAP$outcome[j]]])),
                   parameter = MAP$par, value = unlist(p$tox[MAP$par]), stringsAsFactors = FALSE)
info <- NULL
for (k in names(ARMS)) for (o in c("neut_g3", "stom_any", "stom_g2", "stom_g3", "pneu_any", "pneu_g3", "diar_g2",
                                   "diar_g3", "alt_g3", "ae_disc_cdk", "ae_disc_eve", "ae_disc_et"))
  info <- rbind(info, data.frame(target = paste0("info:", o), regimen = k, observed_target = NA,
                                 model = mean(res[[k]][[o]]), parameter = "", value = NA))
utils::write.csv(rbind(rows, info), file.path(RES, "tox_calibration.csv"), row.names = FALSE)
print(rows, digits = 4)
dm <- disc_all(p)
dtab <- data.frame(DISC_TARGETS[, c("key", "label", "target", "n", "source")], model = round(dm, 4),
                   z = round(disc_z(dm), 2), stringsAsFactors = FALSE)
utils::write.csv(dtab, file.path(RES, "disc_calibration.csv"), row.names = FALSE)
print(dtab[, c("key", "target", "model", "z")])
print(unlist(p$tox[DISC_PARS]))
save_json(list(tox = c(p$tox[MAP$par], p$tox[DISC_PARS])), "params_tox.json")
progress_note("02_calibrate_tox", sprintf(
  "event toxicities max abs error %.4f; discontinuation sub-model fitted (10 targets, max |z| %.2f; hazard falls with time constant %.0f d to %.0f%% of its start)",
  max(err), max(abs(disc_z(dm))), p$tox$disc_tau, 100 * p$tox$disc_floor))
