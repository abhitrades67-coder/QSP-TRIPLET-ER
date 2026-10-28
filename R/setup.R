# Loads the QSP-TRIPLET-ER R model: source this file from any script.
#   source(file.path(<repository root>, "R", "setup.R"))
# Sets ROOT / RES / FIG / TAB, sources the modules and compiles (or loads the cached) C++ engine.

local({
  if (!exists("ROOT", envir = globalenv())) {
    args <- commandArgs(trailingOnly = FALSE)
    f <- sub("^--file=", "", args[grep("^--file=", args)])
    root <- if (length(f)) normalizePath(file.path(dirname(f[1]), ".."), winslash = "/") else getwd()
    if (!file.exists(file.path(root, "R", "engine.cpp"))) root <- normalizePath(getwd(), winslash = "/")
    assign("ROOT", root, envir = globalenv())
  }
})
RES <- Sys.getenv("QSP_RES", file.path(ROOT, "results"))   # override for smoke tests
FIG <- Sys.getenv("QSP_FIG", file.path(ROOT, "figures"))
TAB <- file.path(RES, "tables")
for (d in c(RES, FIG, TAB)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

for (f in c("params.R", "runner.R", "signaling.R", "vpop.R", "engine.R", "analysis.R", "metrics.R",
            "summaries.R", "arms.R", "inputs.R", "de.R")) {
  fn <- file.path(ROOT, "R", f)
  if (file.exists(fn)) source(fn)
}
compile_engine(ROOT)

# parameter files written by the calibration steps (in order)
PFILES <- c("params_tox.json", "params_stage1.json", "params_stage2.json", "params_stage3.json",
            "params_stage4.json")
