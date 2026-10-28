# Start state of the efficacy calibration. The joint refit of stages 1 and 3 (03_calibrate_efficacy.R stage13) starts
# from the parameter files in results/: this step copies the starting values of stages 1 to 4 (calibration_start/) into
# results/ (params_tox.json is fitted from the defaults by 02_calibrate_tox.R) and removes results/structural.csv,
# because 07b_structural.R skips structures that are already in that file.
# Usage: Rscript scripts/03a_start_state.R

if (!exists("ROOT")) ROOT <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), ".."), winslash = "/")
RES <- Sys.getenv("QSP_RES", file.path(ROOT, "results"))
START <- file.path(ROOT, "calibration_start")
dir.create(RES, showWarnings = FALSE, recursive = TRUE)
for (f in c("params_stage1.json", "params_stage2.json", "params_stage3.json", "params_stage4.json"))
  if (!file.copy(file.path(START, f), file.path(RES, f), overwrite = TRUE)) stop("cannot copy ", f, " from ", START)
unlink(file.path(RES, "structural.csv"))
cat("start state: starting values copied from", START, "; structural.csv removed\n")
