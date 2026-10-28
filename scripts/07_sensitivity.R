# Global sensitivity (Latin hypercube + partial rank correlation, PRCC) and expected value of partial
# perfect information (EVPPI, regression method) for the A-versus-B decision. The sampled parameters include
# attrition (ranges around the real-world data, R/inputs.R), the exposure dependence (gamma), time constant and
# floor (late hazard) of AE discontinuation and the new-lesion hazard.
# Writes results/gsa_samples.csv, results/gsa_prcc.csv, results/gsa_evppi.csv
# Usage: Rscript scripts/07_sensitivity.R [n_samples]

if (!exists("ROOT")) ROOT <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), ".."), winslash = "/")
source(file.path(ROOT, "R", "setup.R"))
args <- commandArgs(trailingOnly = TRUE)
NS <- if (length(args) >= 1) as.integer(args[1]) else 160L
NPAT <- as.integer(Sys.getenv("NPAT", "1200"))

# path, mode ('abs' uniform, 'logabs' log-uniform, 'mult' log-uniform multiplier, 'add' additive), lo, hi
SPACE <- read.table(header = TRUE, stringsAsFactors = FALSE, text = sprintf("
path                 mode    lo       hi
tum.kappa            abs     0.0      1.0
ph.IC50_e            logabs  0.08     0.6
tox.tox_int          abs     1.0      2.0
tox.disc_eve_gamma   abs     0.0      1.5
tox.disc_tau         mult    0.5      2.0
tox.disc_floor       abs     0.0      0.25
ph.gir_ddi           abs     1.0      1.4
sig.s_mut            mult    0.5      2.0
tum.sigma_R          mult    0.333333 3.0
tum.log10_fE0        add     -0.5     0.5
tum.log10_fR0        add     -0.5     0.5
tum.b0_mean          mult    0.85     1.15
tum.p_attr           abs     %g       %g
tum.p_attr_mid       abs     %g       %g
tum.dc0_med          mult    0.5      2.0
tum.h_t              mult    0.5      2.0
tum.h_post           mult    0.5      2.0
ph.uG_ref            mult    0.333333 3.0
sig.eps_g            mult    0.75     1.25
ph.I10_e             mult    0.85     1.15
tum.early_prog_days  abs     120      270
tum.lam_NL           mult    0.5      2.0
", ATTR_GSA$p_attr[1], ATTR_GSA$p_attr[2], ATTR_GSA$p_attr_mid[1], ATTR_GSA$p_attr_mid[2]))

make_sample <- function(p0, u) {
  p <- p0; rec <- list()
  for (j in seq_len(nrow(SPACE))) {
    path <- SPACE$path[j]; lo <- SPACE$lo[j]; hi <- SPACE$hi[j]; x <- u[j]
    val <- switch(SPACE$mode[j],
                  abs = lo + x * (hi - lo),
                  logabs = exp(log(lo) + x * (log(hi) - log(lo))),
                  mult = get_param(p0, path) * exp(log(lo) + x * (log(hi) - log(lo))),
                  add = get_param(p0, path) + lo + x * (hi - lo))
    if (path %in% c("sig.eps_g", "ph.I10_e")) val <- min(val, 0.95)
    if (path == "tum.early_prog_days") val <- as.integer(val)
    p <- set_param(p, path, val); rec[[path]] <- val
  }
  list(p = p, rec = rec)
}

prcc <- function(X, y) {
  Xr <- apply(X, 2, rank); yr <- rank(y)
  n <- length(yr); k <- ncol(X) - 1
  t(sapply(seq_len(ncol(X)), function(j) {
    A <- cbind(1, Xr[, -j, drop = FALSE])
    rx <- stats::lm.fit(A, Xr[, j])$residuals
    ry <- stats::lm.fit(A, yr)$residuals
    r <- stats::cor(rx, ry)
    tstat <- r * sqrt((n - 2 - k) / max(1 - r * r, 1e-12))
    c(PRCC = r, p_value = 2 * (1 - stats::pt(abs(tstat), n - 2 - k)))
  }))
}

# regression-based EVPPI for a two-option decision with incremental net benefit nb
evppi <- function(x, nb, nbins = 8) {
  o <- order(x)
  bins <- split(o, cut(seq_along(o), nbins, labels = FALSE))
  cond <- numeric(length(nb))
  for (b in bins) cond[b] <- mean(nb[b])
  mean(pmax(cond, 0)) - max(mean(nb), 0)
}

p0 <- load_params(PFILES)
U <- lhs_unit(NS, nrow(SPACE), seed = 7)
arms <- strategy_arms()
t0 <- Sys.time()
out <- vector("list", NS)
for (i in seq_len(NS)) {
  sm <- make_sample(p0, U[i, ])
  p <- sm$p
  v <- make_vpop(p, NPAT, 1001)
  grid <- make_grid(p)
  s <- lapply(c(A = "A", B = "B", C = "C"), function(k) summarize_arm(simulate(p, v, grid, arms[[k]], 2190L)))
  o <- sm$rec
  for (cc in list(c("A", "B"), c("A", "C"), c("B", "C"))) {
    tag <- paste0(cc, collapse = "")
    o[[paste0("dEDCT48_", tag)]] <- s[[cc[1]]]$EDCT_RMST48 - s[[cc[2]]]$EDCT_RMST48
    o[[paste0("dOS60_", tag)]] <- s[[cc[1]]]$OS_RMST60 - s[[cc[2]]]$OS_RMST60
    o[[paste0("dQTW_", tag)]] <- s[[cc[1]]]$QTWiST - s[[cc[2]]]$QTWiST
  }
  o$AEdisc_A <- s$A$ae_disc_any; o$AEdisc_B <- s$B$ae_disc_any; o$id <- i - 1L
  out[[i]] <- as.data.frame(o, check.names = FALSE)
  if (i %% 20 == 0) cat(sprintf("  %d / %d samples (%.0f s)\n", i, NS, as.numeric(difftime(Sys.time(), t0, units = "secs"))))
}
df <- do.call(rbind, out)
utils::write.csv(df, file.path(RES, "gsa_samples.csv"), row.names = FALSE)
names_ <- SPACE$path
X <- as.matrix(df[, names_])
rows <- list()
for (yname in c("dEDCT48_AB", "dOS60_AB", "dQTW_AB", "AEdisc_A", "dEDCT48_AC", "dEDCT48_BC")) {
  r <- prcc(X, df[[yname]])
  rows[[yname]] <- data.frame(output = yname, parameter = names_, PRCC = r[, "PRCC"], p_value = r[, "p_value"])
}
pr <- do.call(rbind, rows)
utils::write.csv(pr, file.path(RES, "gsa_prcc.csv"), row.names = FALSE)
erows <- list()
for (yname in c("dEDCT48_AB", "dQTW_AB", "dOS60_AB")) {
  nb <- df[[yname]]
  erows[[length(erows) + 1]] <- data.frame(output = yname, parameter = "EVPI (all)", value = mean(pmax(nb, 0)) - max(mean(nb), 0))
  for (j in seq_along(names_)) erows[[length(erows) + 1]] <- data.frame(output = yname, parameter = names_[j],
                                                                          value = evppi(X[, j], nb))
}
utils::write.csv(do.call(rbind, erows), file.path(RES, "gsa_evppi.csv"), row.names = FALSE)
x <- pr[pr$output == "dEDCT48_AB", ]
print(x[order(x$PRCC), ], digits = 3)
cat(sprintf("%d samples in %.0f s\n", NS, as.numeric(difftime(Sys.time(), t0, units = "secs"))))
progress_note("07_sensitivity", sprintf("%d Latin hypercube samples", NS))
