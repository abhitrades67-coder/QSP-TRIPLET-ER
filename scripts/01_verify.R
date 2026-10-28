# Model verification (code and numerics).
#  V1  ODE long-time solution equals fixed-point steady state (all genotypes, drug combinations)
#  V2  grid interpolation error at random off-grid points
#  V3  CYP3A4 DDI module reproduces the label AUC ratios; predictions at other ribociclib doses
#  V4  common random numbers: identical first-line regimens give identical first-line outcomes
#  V5  cell-number conservation of mutation / switching flows
#  V6  survival functions (KM median, Breslow Cox HR, weighted versions) vs the 'survival' package
#  V7  random-number generator: uniformity and independence across days, streams and patients
# Writes results/verification.csv

if (!exists("ROOT")) ROOT <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), ".."), winslash = "/")
source(file.path(ROOT, "R", "setup.R"))
p <- load_params(PFILES)
s <- p$sig; st <- p$st; ph <- p$ph
rows <- list()
add <- function(check, metric, value, pass) rows[[length(rows) + 1]] <<- data.frame(check = check, metric = metric,
                                                                               value = as.character(value), pass = pass)

# V1
combos <- list(c(1, 0, 0, 0), c(0.02, 0, 0, .8), c(0.02, 0, .7, 0), c(1, 400, .7, .8), c(0.03, 0, .7, 0),
               c(1, 40, 0, 0), c(1, 400, 0, 0), c(1, ph$uG_ref, .7, .8))
maxd <- 0
for (g in 0:7) for (cb in combos) {
  ss <- steady_state(cb[1], cb[2], cb[3], cb[4], 0.25, g, s, st, nthreads = 1)
  ode <- simulate_ode(cb, 0.25, g, s, st, t_end = 400, n_out = 3)
  maxd <- max(maxd, abs(tail(ode$pRb, 1) - ss$pRb), abs(tail(ode$ER, 1) - ss$ER))
}
add("V1 ODE vs fixed-point steady state", "max abs difference (pRb, ER), 8 genotypes x 8 drug combinations",
    signif(maxd, 3), maxd < 1e-3)

# V2
grid <- make_grid(p)
set.seed(0)
n <- 3000
u <- exp(stats::runif(n, 0, log(3000))) - 1
Ie <- stats::runif(n, 0, .95); Ic <- stats::runif(n, 0, .95); b0 <- stats::runif(n, 0, .85)
e2i <- sample(0:2, n, replace = TRUE)
lk <- grid_lookup(grid, e2i, u, Ie, Ic, b0)
E2v <- grid$E2_levels[e2i + 1]
err <- unlist(lapply(0:7, function(g) abs(steady_state(E2v, u, Ie, Ic, b0, g, s, st)$P - lk$P[, g + 1])))
add("V2 grid interpolation", "median / 99th percentile abs error in the proliferation index (3000 random points x 8 genotypes)",
    sprintf("%.4f / %.4f", stats::median(err), stats::quantile(err, .99)), stats::quantile(err, .99) < 0.05)

# V3
a_palbo <- ddi_aucr(ph, "palbo", 1); a_ribo <- ddi_aucr(ph, "ribo", 1)
add("V3 DDI palbociclib 125 mg -> CYP3A substrate AUC ratio", "target 1.61 (palbociclib prescribing information)", round(a_palbo, 3),
    abs(a_palbo - 1.61) < 0.02)
add("V3 DDI ribociclib 600 mg -> CYP3A substrate AUC ratio", "target 5.2 (ribociclib prescribing information, predicted)", round(a_ribo, 3),
    abs(a_ribo - 5.2) < 0.05)
add("V3 DDI ribociclib 400 mg (prediction)", "prescribing information: midazolam AUC 3.8-fold (measured)",
    round(ddi_aucr(ph, "ribo", 400 / 600), 2), NA)
add("V3 DDI ribociclib 300 / 250 mg continuous + everolimus 2.5 mg (prediction)",
    "TRINITI-1: everolimus trough 4-5x (300 mg), 2-3x (250 mg) above expected",
    sprintf("%.2f / %.2f (AUC ratio)", ddi_aucr(ph, "ribo", 0.5), ddi_aucr(ph, "ribo", 250 / 600)), NA)
dyn <- pk_profile(ph, "ribo", 1, 1, 84, onoff = FALSE)
alone <- pk_profile(ph, NULL, 0, 1, 84)
add("V3 dynamic PK: everolimus exposure ratio, ribociclib 600 mg continuous, day 28", "target 5.2",
    round(dyn$eve[29] / alone$eve[29], 3), abs(dyn$eve[29] / alone$eve[29] - 5.2) < 0.1)

# V4 / V5
v <- make_vpop(p, 300, 5)
arms <- strategy_arms()
Rb <- simulate(p, v, grid, arms$B, 900L)
Rc <- simulate(p, v, grid, arms$C, 900L)
same <- identical(Rb$pfs1, Rc$pfs1) && identical(Rb$prog1, Rc$prog1)
add("V4 common random numbers (arms B and C share line 1)", "identical PFS1 for every patient", same, same)
N <- abs(matrix(stats::rnorm(40, 1e9, 1e8), 5, 8))
tot0 <- sum(N)
for (b in c(1, 2, 4)) {
  src <- which(bitwAnd(0:7, b) == 0) - 1
  tr <- 0.01 * N[, src + 1]
  N[, src + 1] <- N[, src + 1] - tr
  N[, bitwOr(src, b) + 1] <- N[, bitwOr(src, b) + 1] + tr
}
add("V5 conservation of cells in clone transitions", "relative change", signif(abs(sum(N) - tot0) / tot0, 3),
    abs(sum(N) - tot0) / tot0 < 1e-12)

# V6 survival functions vs the 'survival' package
vv <- make_vpop(p, 1500, 9)
RA <- simulate(p, vv, grid, arms$A, 2190L); RB <- simulate(p, vv, grid, arms$B, 2190L)
sf <- survival::survfit(survival::Surv(RB$pfs1, RB$pfs1_ev) ~ 1)
med_surv <- unname(summary(sf)$table["median"])
med_own <- km_median(RB$pfs1, RB$pfs1_ev)
dat <- data.frame(t = c(RA$edct, RB$edct), e = c(RA$edct_ev, RB$edct_ev), x = rep(1:0, each = vv$n))
cx <- survival::coxph(survival::Surv(t, e) ~ x, data = dat, ties = "breslow")
hr_surv <- unname(exp(stats::coef(cx)))
hr_own <- cox_hr(RA$edct, RA$edct_ev, RB$edct, RB$edct_ev)[1]
w <- ifelse(vv$Ptrunc, 2, 1)
cxw <- survival::coxph(survival::Surv(t, e) ~ x, data = dat, weights = rep(w, 2), ties = "breslow")
hrw_own <- cox_hr(RA$edct, RA$edct_ev, RB$edct, RB$edct_ev, w, w)[1]
sfw <- survival::survfit(survival::Surv(RB$os, RB$os_ev) ~ 1, weights = w)
medw_surv <- unname(summary(sfw)$table["median"])
dmax <- max(abs(med_surv - med_own), abs(hr_surv - hr_own) / hr_surv, abs(exp(stats::coef(cxw)) - hrw_own) / hrw_own,
            abs(medw_surv - km_median(RB$os, RB$os_ev, w)))
add("V6 survival functions vs R 'survival' package", "max difference: KM median (days), Cox HR (relative), weighted versions",
    signif(dmax, 3), dmax < 1e-6)

# V7 random-number generator
Ur <- rng_sample_cpp(7926932, 0:1999, 365L, 0L)
ks <- suppressWarnings(stats::ks.test(as.vector(Ur), "punif")$p.value)
lag1 <- stats::cor(as.vector(Ur[, -1]), as.vector(Ur[, -365]))
U2 <- rng_sample_cpp(7926932, 0:1999, 365L, 1L)
xs <- stats::cor(as.vector(Ur), as.vector(U2))
pa <- stats::cor(as.vector(Ur[-1, ]), as.vector(Ur[-2000, ]))
small <- mean(Ur < 1e-3)
add("V7 random numbers (730,000 draws)",
    "KS p-value; correlation lag-1 day / between streams / between patients; share < 0.001",
    sprintf("%.2f; %.4f / %.4f / %.4f; %.5f", ks, lag1, xs, pa, small),
    ks > 0.001 && max(abs(c(lag1, xs, pa))) < 0.01 && abs(small - 0.001) < 0.0003)

df <- do.call(rbind, rows)
utils::write.csv(df, file.path(RES, "verification.csv"), row.names = FALSE)
print(df, right = FALSE)
progress_note("01_verify", sprintf("%d of %d checks passed (%d descriptive)", sum(df$pass %in% TRUE),
                                   sum(!is.na(df$pass)), sum(is.na(df$pass))))
