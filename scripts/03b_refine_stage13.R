# Local refinement of the joint stage 1 + 3 fit: compass (pattern) search from the differential-evolution result.
# The differential evolution (03_calibrate_efficacy.R stage13, 30 generations) did not improve on its best starting
# point: its population, spread over the full bounds, produced trial points far from the
# optimum. A compass search tries each parameter in turn, one step up and one step down (in the transformed space of
# the fit), keeps any improvement, extends a successful move, and halves the steps when a full sweep improves nothing.
# Same loss, parameters, bounds and calibration population (800 patients, seed 21) as the joint fit.
# Writes results/params_stage13.json, params_stage1.json, params_stage3.json and calib_stage13.json (the
# differential-evolution files are kept as *_de.json) and results/refine_stage13_history.csv.
# Usage: Rscript scripts/03b_refine_stage13.R [max_sweeps]

if (!exists("ROOT")) ROOT <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), ".."), winslash = "/")
CAL_AS_LIBRARY <- TRUE
source(file.path(ROOT, "scripts", "03_calibrate_efficacy.R"))
args <- commandArgs(trailingOnly = TRUE)
MAX_SWEEPS <- if (length(args) >= 1) as.integer(args[1]) else 12L

stage <- "stage13"
P <- STAGES[[stage]]$pars
b <- bounds_of(stage)
p <- load_params(PFILES)
v <- make_vpop(p, NCAL, SEED)
theta <- sapply(seq_len(nrow(P)), function(j) { val <- get_param(p, P$path[j]); if (P$lg[j]) log(val) else val })
clamp <- function(x, j) min(max(x, b$lower[j]), b$upper[j])
nfev <- 0L
fn <- function(th) { nfev <<- nfev + 1L
  tryCatch(evaluate_stage(stage, th, p, v), error = function(e) { cat("eval error:", conditionMessage(e), "\n"); 1e6 }) }
t0 <- Sys.time()
best <- fn(theta)
cat(sprintf("start: loss %.3f\n", best))
hist <- data.frame(sweep = 0L, loss = best, nfev = nfev, step_scale = 1)
step <- 0.10 * (b$upper - b$lower)          # initial step: 10% of each parameter's (transformed) range
min_step <- 0.005 * (b$upper - b$lower)
set.seed(7)
for (sweep in seq_len(MAX_SWEEPS)) {
  improved <- FALSE
  for (j in sample(nrow(P))) {
    for (dir in c(1, -1)) {
      cand <- theta; cand[j] <- clamp(theta[j] + dir * step[j], j)
      if (cand[j] == theta[j]) next
      L <- fn(cand)
      if (L < best - 1e-9) {
        best <- L; theta <- cand; improved <- TRUE
        repeat {   # extend a successful move
          cand <- theta; cand[j] <- clamp(theta[j] + dir * step[j], j)
          if (cand[j] == theta[j]) break
          L <- fn(cand)
          if (L < best - 1e-9) { best <- L; theta <- cand } else break
        }
        cat(sprintf("  sweep %d: %s -> %.4g, loss %.3f (%d evaluations, %.0f s)\n", sweep, P$path[j],
                    if (P$lg[j]) exp(theta[j]) else theta[j], best, nfev, as.numeric(difftime(Sys.time(), t0, units = "secs"))))
        break
      }
    }
  }
  if (!improved) step <- step / 2
  hist <- rbind(hist, data.frame(sweep = sweep, loss = best, nfev = nfev, step_scale = step[1] / (0.10 * (b$upper[1] - b$lower[1]))))
  cat(sprintf("sweep %d done: loss %.3f, %s, %d evaluations\n", sweep, best, if (improved) "improved" else "no improvement (steps halved)", nfev))
  if (all(step < min_step)) break
}
fin <- evaluate_stage(stage, theta, p, v, NULL, return_res = TRUE)
cat(stage, "refined loss", fin$loss, "\n"); print(unlist(fin$info))
# keep the differential-evolution files, then write the refined fit in the same format
for (f in c("params_stage13", "params_stage1", "params_stage3", "calib_stage13")) {
  src <- file.path(RES, paste0(f, ".json")); dst <- file.path(RES, paste0(f, "_de.json"))
  if (file.exists(src) && !file.exists(dst)) file.copy(src, dst)
}
pars <- list()
for (j in seq_len(nrow(P))) { sk <- strsplit(P$path[j], ".", fixed = TRUE)[[1]]; pars[[sk[1]]][[sk[2]]] <- if (P$lg[j]) exp(theta[j]) else theta[j] }
save_json(pars, "params_stage13.json")
split_pars <- function(paths) { out <- list(); for (pp in paths) { sk <- strsplit(pp, ".", fixed = TRUE)[[1]]
  out[[sk[1]]][[sk[2]]] <- pars[[sk[1]]][[sk[2]]] }; out }
save_json(split_pars(STAGES$stage1$pars$path), "params_stage1.json")
save_json(split_pars(STAGES$stage3w$pars$path), "params_stage3.json")
de <- load_json("calib_stage13_de.json")
save_json(list(loss = fin$loss, info = fin$info, x = theta, nfev = de$nfev + nfev, generations = de$generations,
               de_loss = de$loss, refine_nfev = nfev, refine_sweeps = nrow(hist) - 1L,
               seconds = de$seconds + as.numeric(difftime(Sys.time(), t0, units = "secs"))), "calib_stage13.json")
utils::write.csv(hist, file.path(RES, "refine_stage13_history.csv"), row.names = FALSE)
progress_note("03b_refine_stage13", sprintf("compass search: loss %.2f -> %.2f in %d evaluations", de$loss, fin$loss, nfev))
