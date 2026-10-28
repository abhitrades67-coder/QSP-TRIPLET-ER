# Virtual population (first-line ER+/HER2- advanced breast cancer).
#
# All patient-level randomness is drawn once, as z-scores and uniforms, so that the same patient
# (index) can be simulated under every arm (paired, digital-twin design) and under every
# calibration candidate (common random numbers). refresh_derived() maps these draws to patient
# characteristics for the current population-level parameters.

VPOP_FIELDS <- c("kmax", "rd", "b0", "Ptrunc", "fE0", "fR0", "N0", "expo_palbo", "expo_ribo",
                 "expo_abema", "expo_eve", "expo_gir", "expo_fulv", "ANC0", "tox_frailty",
                 "neut_sens", "age", "dc0")
DRUGS <- c("palbo", "ribo", "abema", "eve", "gir", "fulv")

make_vpop <- function(p, n, seed = 1) {
  old <- if (exists(".Random.seed", envir = globalenv())) get(".Random.seed", envir = globalenv()) else NULL
  on.exit(if (!is.null(old)) assign(".Random.seed", old, envir = globalenv()))
  set.seed(seed, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
  v <- list(n = as.integer(n), seed = seed)
  v$z_fE <- rnorm(n); v$z_fR <- rnorm(n); v$u_base <- runif(n); v$u_base2 <- runif(n)
  v$z_N0 <- rnorm(n)
  for (d in DRUGS) v[[paste0("z_expo_", d)]] <- rnorm(n)
  v$z_ANC0 <- rnorm(n); v$z_frail <- rnorm(n); v$z_neut <- rnorm(n); v$z_age <- rnorm(n)
  v$z_k <- rnorm(n); v$z_rd <- rnorm(n); v$u_b0 <- runif(n); v$u_P <- runif(n); v$z_dc <- rnorm(n)
  v$seed_events <- seed * 7919 + 13
  refresh_derived(p, v)
}

refresh_derived <- function(p, v) {
  t <- p$tum; ph <- p$ph; tx <- p$tox
  fE0 <- 10^(t$log10_fE0 + t$sd_log10_fE0 * v$z_fE)
  baseE <- v$u_base < t$pE_base
  v$fE0 <- pmin(pmax(ifelse(baseE, 0.02 + 0.28 * v$u_base2, fE0), 0), 0.5)
  v$fR0 <- pmin(pmax(10^(t$log10_fR0 + t$sd_log10_fR0 * v$z_fR), 0), 0.2)
  v$kmax <- t$kmax_med * exp(t$kmax_sd * v$z_k)
  v$rd <- pmin(pmax(t$rd + t$rd_sd * v$z_rd, 0.3), 0.98)
  a <- t$b0_mean * t$b0_conc
  b <- (1 - t$b0_mean) * t$b0_conc
  v$b0 <- qbeta(pmin(pmax(v$u_b0, 1e-6), 1 - 1e-6), a, b)
  v$Ptrunc <- v$u_P < t$pP
  v$dc0 <- t$dc0_med * exp(t$dc0_sd * v$z_dc)
  v$N0 <- t$N0_med * exp(t$N0_sd * v$z_N0)
  for (d in DRUGS) {
    sd <- ph[[paste0("sd_", d)]]
    v[[paste0("expo_", d)]] <- exp(-0.5 * sd^2 + sd * v[[paste0("z_expo_", d)]])
  }
  v$ANC0 <- tx$ANC0_med * exp(tx$ANC0_sd * v$z_ANC0)
  v$tox_frailty <- exp(-0.125 + 0.5 * v$z_frail)   # shared susceptibility to adverse events
  v$neut_sens <- exp(-0.08 + 0.4 * v$z_neut)       # myelosuppression sensitivity
  v$age <- pmin(pmax(61 + 11 * v$z_age, 28), 85)
  v
}

build_vpop <- function(p, n, seed = 1) make_vpop(p, n, seed)

# first m patients of a population (for trajectories)
vpop_head <- function(v, m) {
  w <- v
  for (k in names(v)) if (length(v[[k]]) == v$n && k != "n") w[[k]] <- v[[k]][seq_len(m)]
  w$n <- as.integer(m)
  w
}

