# Helpers: paths, parameter I/O, compilation of the C++ engine.

qsp_threads <- function() {
  n <- suppressWarnings(as.integer(Sys.getenv("QSP_THREADS", "")))
  if (is.na(n) || n < 1) n <- max(1L, parallel::detectCores(logical = TRUE))
  n
}

deep_update <- function(p, upd) {
  for (k in names(upd)) {
    if (is.list(upd[[k]]) && !is.null(p[[k]]) && is.list(p[[k]])) {
      p[[k]] <- deep_update(p[[k]], upd[[k]])
    } else {
      p[[k]] <- upd[[k]]
    }
  }
  p
}

save_json <- function(obj, name, dir = RES) {
  jsonlite::write_json(obj, file.path(dir, name), auto_unbox = TRUE, digits = NA, pretty = TRUE,
                       na = "null")
}

load_json <- function(name, dir = RES) {
  jsonlite::fromJSON(file.path(dir, name), simplifyVector = TRUE, simplifyDataFrame = FALSE,
                     simplifyMatrix = FALSE)
}

# default parameters (R/params.R) updated by the parameter files that exist in `dir`, in order
load_params <- function(files = PFILES, dir = RES, base = default_params()) {
  p <- base
  for (nm in files) {
    fn <- file.path(dir, nm)
    if (file.exists(fn)) p <- deep_update(p, load_json(nm, dir))
  }
  p
}

set_param <- function(p, path, value) {
  sk <- strsplit(path, ".", fixed = TRUE)[[1]]
  p[[sk[1]]][[sk[2]]] <- value
  p
}

get_param <- function(p, path) {
  sk <- strsplit(path, ".", fixed = TRUE)[[1]]
  p[[sk[1]]][[sk[2]]]
}

compile_engine <- function(root) {
  cache <- file.path(root, ".cache")
  dir.create(cache, showWarnings = FALSE)
  # no fused multiply-add contraction, for reproducible floating-point results
  Sys.setenv(PKG_CXXFLAGS = "-O2 -ffp-contract=off")
  Rcpp::sourceCpp(file.path(root, "R", "engine.cpp"), cacheDir = cache, rebuild = FALSE,
                  env = globalenv(), showOutput = FALSE)
}

progress_note <- function(step, text) {
  if (nzchar(Sys.getenv("QSP_RES"))) return(invisible())   # smoke tests do not log progress
  fn <- file.path(ROOT, "PROGRESS.md")
  line <- sprintf("- %s  **%s**: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M"), step, text)
  cat(line, file = fn, append = TRUE)
}
