# Differential evolution, following scipy.optimize.differential_evolution: strategy 'best1bin',
# mutation dithered in [0.5, 1) once per generation,
# recombination 0.7, 'deferred' updating (the whole generation is evaluated before selection),
# quasi-random initial population (Latin hypercube here; Sobol in scipy) with x0 as the first
# member, stopping when sd(energies) <= atol + tol * |mean(energies)| or after maxiter generations.

de_optimize <- function(fn, lower, upper, x0 = NULL, popsize = 6, maxiter = 20, tol = 1e-3, atol = 0,
                        mutation = c(0.5, 1), recombination = 0.7, seed = 3, trace = TRUE) {
  set.seed(seed)
  d <- length(lower)
  np_ <- max(5L, as.integer(popsize * d))
  # Latin hypercube in the unit cube
  pop <- sapply(seq_len(d), function(j) (sample.int(np_) - stats::runif(np_)) / np_)
  if (!is.matrix(pop)) pop <- matrix(pop, nrow = np_)
  scale_x <- function(u) lower + u * (upper - lower)
  if (!is.null(x0)) {   # one starting point (vector) or several (matrix, one per row)
    x0 <- matrix(x0, ncol = d)
    for (r in seq_len(min(nrow(x0), np_))) pop[r, ] <- pmin(pmax((x0[r, ] - lower) / (upper - lower), 0), 1)
  }
  energies <- vapply(seq_len(np_), function(i) fn(scale_x(pop[i, ])), numeric(1))
  energies[!is.finite(energies)] <- Inf
  nfev <- np_
  best <- which.min(energies)
  hist <- data.frame(gen = 0L, best = energies[best], mean = mean(energies[is.finite(energies)]))
  t0 <- Sys.time()
  if (trace) cat(sprintf("  init: best %.4f (%d evaluations, %.0f s)\n", energies[best], nfev,
                         as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  for (gen in seq_len(maxiter)) {
    Fm <- stats::runif(1, mutation[1], mutation[2])
    trials <- pop
    for (i in seq_len(np_)) {
      r <- sample(setdiff(seq_len(np_), i), 2)
      bprime <- pop[best, ] + Fm * (pop[r[1], ] - pop[r[2], ])
      cross <- stats::runif(d) < recombination
      cross[sample.int(d, 1)] <- TRUE
      tr <- ifelse(cross, bprime, pop[i, ])
      out <- tr < 0 | tr > 1
      tr[out] <- stats::runif(sum(out))
      trials[i, ] <- tr
    }
    te <- vapply(seq_len(np_), function(i) fn(scale_x(trials[i, ])), numeric(1))
    te[!is.finite(te)] <- Inf
    nfev <- nfev + np_
    better <- te <= energies
    pop[better, ] <- trials[better, , drop = FALSE]
    energies[better] <- te[better]
    best <- which.min(energies)
    fin <- energies[is.finite(energies)]
    hist <- rbind(hist, data.frame(gen = gen, best = energies[best], mean = mean(fin)))
    if (trace) cat(sprintf("  gen %d: best %.4f, mean %.3f (%d evaluations, %.0f s)\n", gen, energies[best],
                           mean(fin), nfev, as.numeric(difftime(Sys.time(), t0, units = "secs"))))
    if (all(is.finite(energies)) && stats::sd(energies) <= atol + tol * abs(mean(energies))) break
  }
  list(par = scale_x(pop[best, ]), value = energies[best], nfev = nfev, history = hist,
       population = t(apply(pop, 1, scale_x)), energies = energies)
}

# Latin hypercube sample (n x d) in the unit cube
lhs_unit <- function(n, d, seed = 7) {
  set.seed(seed)
  sapply(seq_len(d), function(j) (sample.int(n) - stats::runif(n)) / n)
}
