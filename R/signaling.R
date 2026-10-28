# ER / PI3K-mTOR / CDK4/6 / Rb signaling module.
#
# Two tiers:
#   * simulate_ode(): dynamic ODE form (hours, deSolve::lsoda), used for time-course figures and
#     for verification;
#   * steady_state(): fixed-point solution of the same equations (C++), used to pre-compute the
#     grid that the population engine interpolates every simulated day (quasi-steady-state
#     reduction: signaling settles in hours, tumor size changes over weeks).
# All quantities are normalized to an untreated wild-type cell (= 1).
# Genotype code g = E + 2*P + 4*R (E: ESR1 mutant, P: PI3K pathway altered, R: CDK4/6i escape).

GENO_NAMES <- c("S", "E", "P", "EP", "R", "ER", "PR", "EPR")
U_GRID <- c(0, 1, 3, 10, 30, 100, 300, 1000, 3000)
IE_GRID <- c(0, .1, .2, .3, .4, .5, .6, .7, .8, .9, .98)
IC_GRID <- c(0, .1, .2, .3, .4, .5, .6, .7, .8, .9, .95)
B0_GRID <- c(0, .05, .1, .15, .2, .3, .4, .55, .7, .9)
K_ODE <- c(ER = 0.12, IRS = 0.25, AKT = 2.0, CycD = 0.35, CycE = 0.25, pRb = 0.6)   # 1/hour

pRb0_of <- function(s) 1 / (s$KRb^s$nRb + 1)

steady_state <- function(E2, u, Ie, Ic, b0, g, s, st, n_iter = 400L, damp = 0.5, tol = 1e-8,
                         nthreads = qsp_threads()) {
  n <- max(length(E2), length(u), length(Ie), length(Ic), length(b0), length(g))
  steady_state_cpp(rep_len(as.numeric(E2), n), rep_len(as.numeric(u), n), rep_len(as.numeric(Ie), n),
                   rep_len(as.numeric(Ic), n), rep_len(as.numeric(b0), n), rep_len(as.integer(g), n),
                   s, st, as.integer(n_iter), damp, tol, as.integer(nthreads))
}

# Grid axes: E2 level (untreated, exemestane, letrozole), genotype, SERD index u, everolimus
# inhibition Ie, CDK4/6 inhibition Ic, ER-independent cyclin D drive b0. Stored in C order
# (last axis fastest).
make_grid <- function(p, nthreads = qsp_threads()) {
  E2_levels <- c(1.0, p$ph$E2_exe, p$ph$E2_let)
  eg <- expand.grid(b = B0_GRID, ic = IC_GRID, ie = IE_GRID, u = U_GRID, g = 0:7, e2 = E2_levels,
                    KEEP.OUT.ATTRS = FALSE)
  ss <- steady_state_cpp(eg$e2, eg$u, eg$ie, eg$ic, eg$b, eg$g, p$sig, p$st, 400L, 0.5, 1e-8,
                         as.integer(nthreads))
  ref <- steady_state_cpp(1.0, 0.0, 0.0, 0.0, 0.2, 0L, p$sig, p$st, 400L, 0.5, 1e-8, 1L)
  list(P = ss$P, S6K = ss$S6K, A_ER = ss$A_ER, E2_levels = E2_levels, P0_ref = ref$P,
       max_diff = ss$max_diff, n_iter = ss$n_iter)
}

grid_lookup <- function(grid, e2i, u, Ie, Ic, b0) {
  grid_lookup_cpp(grid$P, grid$S6K, as.integer(e2i), as.numeric(u), as.numeric(Ie), as.numeric(Ic),
                  as.numeric(b0))
}

# Dynamic ODE tier. inputs: function(t) returning c(E2, u, Ie, Ic) (t in hours) or a constant vector.
simulate_ode <- function(inputs, b0, g, s, st, t_end = 240, y0 = NULL, n_out = 481) {
  if (is.null(y0)) y0 <- c(1, 1, 1, 1, 1, pRb0_of(s))
  fin <- if (is.function(inputs)) inputs else function(t) inputs
  k <- unname(K_ODE)
  rhs <- function(t, y, parms) {
    x <- fin(t)
    list(sig_rhs_cpp(y, x[1], x[2], x[3], x[4], b0, as.integer(g), s, st, k))
  }
  ts <- seq(0, t_end, length.out = n_out)
  sol <- deSolve::lsoda(y0, ts, rhs, parms = NULL, rtol = 1e-8, atol = 1e-10, hmax = 0.5)
  Y <- sol[, 2:7, drop = FALSE]
  colnames(Y) <- names(K_ODE)
  aux <- t(vapply(seq_along(ts), function(i) {
    x <- fin(ts[i])
    sig_aux_cpp(Y[i, ], x[1], x[2], x[3], x[4], b0, as.integer(g), s, st)
  }, numeric(6)))
  c(list(t = ts), as.list(as.data.frame(Y)), list(S6K = aux[, 1], A_ER = aux[, 3], P = aux[, 6]))
}
