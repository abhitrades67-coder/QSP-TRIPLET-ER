# Survival summaries: Kaplan-Meier, medians, landmarks, RMST, Cox hazard ratio
# (checked against the 'survival' package in scripts/01_verify.R).

MONTH <- 30.4375

months <- function(days) days / MONTH

# Kaplan-Meier estimate at the distinct times; optional case weights
km <- function(time, event, weights = NULL) {
  time <- as.numeric(time)
  event <- as.logical(event)
  w <- if (is.null(weights)) rep(1, length(time)) else as.numeric(weights)
  o <- order(time, method = "radix")
  time <- time[o]; event <- event[o]; w <- w[o]
  n <- length(time)
  if (n == 0) return(list(t = numeric(0), s = numeric(0)))
  first <- which(c(TRUE, diff(time) != 0))
  ut <- time[first]
  last <- c(first[-1] - 1L, n)
  cumw <- c(0, cumsum(w))
  at_risk <- cumw[n + 1] - cumw[first]
  ce <- cumsum(w * event)
  d <- ce[last] - c(0, ce[last][-length(last)])
  s <- cumprod(1 - ifelse(at_risk > 0, d / pmax(at_risk, 1e-12), 0))
  list(t = ut, s = s)
}

km_median <- function(time, event, weights = NULL) {
  k <- km(time, event, weights)
  b <- which(k$s <= 0.5)
  if (length(b)) k$t[b[1]] else NA_real_
}

surv_at <- function(time, event, x, weights = NULL) {
  k <- km(time, event, weights)
  i <- findInterval(x, k$t)   # number of distinct times <= x
  if (i < 1) 1.0 else k$s[i]
}

rmst <- function(time, event, tau, weights = NULL) {
  k <- km(time, event, weights)
  tt <- c(0, k$t[k$t <= tau], tau)
  ss <- c(1, k$s)[seq_len(length(tt) - 1)]
  sum(diff(tt) * ss)
}

# HR of group 1 vs group 0: single binary covariate, Breslow ties, Newton-Raphson.
# Returns c(HR, lower 95%, upper 95%).
cox_hr <- function(t1, e1, t0, e0, w1 = NULL, w0 = NULL, iters = 30) {
  t <- c(as.numeric(t1), as.numeric(t0))
  e <- c(as.logical(e1), as.logical(e0))
  x <- c(rep(1, length(t1)), rep(0, length(t0)))
  w <- c(if (is.null(w1)) rep(1, length(t1)) else w1, if (is.null(w0)) rep(1, length(t0)) else w0)
  if (length(t) == 0 || !any(e)) return(c(NA_real_, NA_real_, NA_real_))
  o <- order(-t, method = "radix")   # descending time, stable
  t <- t[o]; e <- e[o]; x <- x[o]; w <- w[o]
  n <- length(t)
  first <- which(c(TRUE, diff(t) != 0))
  last <- c(first[-1] - 1L, n)
  we <- w * e
  cwe <- cumsum(we); cwex <- cumsum(we * x)
  dw <- cwe[last] - c(0, cwe[last][-length(last)])
  dx <- cwex[last] - c(0, cwex[last][-length(last)])
  beta <- 0
  I <- NA_real_
  for (it in seq_len(iters)) {
    r <- w * exp(beta * x)
    S0 <- cumsum(r)[last]
    S1 <- cumsum(r * x)[last]
    m <- S1 / S0
    U <- sum(dx - dw * m)
    I <- sum(dw * (m - m * m))
    if (!is.finite(I) || I <= 0) break
    step <- U / I
    beta <- beta + step
    if (abs(step) < 1e-9) break
  }
  se <- if (is.finite(I) && I > 0) 1 / sqrt(I) else NA_real_
  c(exp(beta), exp(beta - 1.96 * se), exp(beta + 1.96 * se))
}
