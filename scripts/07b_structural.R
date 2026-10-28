# Structural uncertainty: alternative model structures, each recalibrated on the stage that the
# structural change affects, then arms A, B, C simulated again.
# Structures:
#   base                      : reference model
#   no S6K-ER cross-talk      : bS = 0 (recalibrate stage 2: everolimus effect)
#   no adaptive escape        : CDK4/6i escape only by pre-existing / mutational clones (recalibrate the
#                               stage-1 parameters against the stage-1 and acelERA targets)
#   reversible escape         : escape phenotype reverts when CDK4/6i is off (recalibrated as above)
#   kappa = 0 / kappa = 1     : extremes of the unidentifiable mTORC1 dependence of adaptive escape
#   burden-only death hazard  : survival model without the post-endocrine hazard (h_post = 0,
#                               dc0_med matched to MONALEESA-2 OS alone)
#   constant attrition 10%    : 10% attrition for all line-1 durations
# Writes results/structural.csv
# Usage: Rscript scripts/07b_structural.R [structures ...]   (default: all)

if (!exists("ROOT")) ROOT <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), ".."), winslash = "/")
CAL_AS_LIBRARY <- TRUE
source(file.path(ROOT, "scripts", "03_calibrate_efficacy.R"))
N <- as.integer(Sys.getenv("NSTRUCT", "3000"))
GENS <- as.integer(Sys.getenv("SGENS", "12"))

contrasts_of <- function(p, label, loss = NA) {
  v <- make_vpop(p, N, 1001)
  grid <- make_grid(p)
  arms <- strategy_arms()
  s <- lapply(c(A = "A", B = "B", C = "C"), function(k) summarize_arm(simulate(p, v, grid, arms[[k]], 2190L)))
  do.call(rbind, lapply(list(c("A", "B"), c("A", "C"), c("B", "C")), function(cc)
    data.frame(structure = label, contrast = paste0(cc[1], "-", cc[2]),
               dEDCT48 = s[[cc[1]]]$EDCT_RMST48 - s[[cc[2]]]$EDCT_RMST48,
               dOS60 = s[[cc[1]]]$OS_RMST60 - s[[cc[2]]]$OS_RMST60,
               EDCT48_a = s[[cc[1]]]$EDCT_RMST48, EDCT48_b = s[[cc[2]]]$EDCT_RMST48, calib_loss = loss,
               stringsAsFactors = FALSE)))
}

args <- commandArgs(trailingOnly = TRUE)
which_ <- if (length(args)) args else c("base", "kappa", "nohpost", "attr10", "nos6k", "noadapt", "reversible")
out_fn <- file.path(RES, "structural.csv")
rows <- if (file.exists(out_fn)) utils::read.csv(out_fn, stringsAsFactors = FALSE) else NULL
done <- unique(rows$structure)
add_rows <- function(x) { rows <<- rbind(rows, x); utils::write.csv(rows, out_fn, row.names = FALSE) }
p0 <- load_params(PFILES)

if ("base" %in% which_ && !("base" %in% done)) add_rows(contrasts_of(p0, "base"))
if ("kappa" %in% which_) for (kv in c(0, 1)) {
  lab <- sprintf("kappa = %g", kv)
  if (!(lab %in% done)) { p <- p0; p$tum$kappa <- kv; add_rows(contrasts_of(p, lab)) }
}
if ("nohpost" %in% which_ && !("burden-only death hazard" %in% done)) {
  g <- utils::read.csv(file.path(RES, "stage4_grid.csv"))
  g <- g[g$h_post == 0 & is.na(g$h_t), ]
  g <- g[order(abs(g$ML2_OS - CAL_TARGETS$ML2_OS_cdk)), ]
  p <- p0; p$tum$h_post <- 0; p$tum$dc0_med <- g$dc0_med[1]
  add_rows(contrasts_of(p, "burden-only death hazard"))
}
if ("attr10" %in% which_ && !("constant attrition 10%" %in% done)) {
  p <- p0; p$tum$p_attr <- 0.10; p$tum$p_attr_mid <- 0.10
  add_rows(contrasts_of(p, "constant attrition 10%"))
}
if ("nos6k" %in% which_ && !("no S6K-ER cross-talk" %in% done)) {
  p <- p0; p$sig$bS <- 0; p$st$s6k_er <- FALSE
  pars <- run_de("stage2", p, maxiter = GENS, popsize = 8, tag = "_nos6k")
  p <- deep_update(p, pars)
  add_rows(contrasts_of(p, "no S6K-ER cross-talk", load_json("calib_stage2_nos6k.json")$loss))
}
if ("noadapt" %in% which_ && !("no adaptive escape" %in% done)) {
  p <- p0; p$st$adaptive_R <- FALSE
  pars <- run_de("stage1j", p, maxiter = GENS, popsize = 6, tag = "_noadapt")
  p <- deep_update(p, pars)
  add_rows(contrasts_of(p, "no adaptive escape", load_json("calib_stage1j_noadapt.json")$loss))
}
if ("reversible" %in% which_ && !("reversible escape" %in% done)) {
  p <- p0; p$st$reversible_R <- TRUE
  pars <- run_de("stage1j", p, maxiter = GENS, popsize = 6, tag = "_reversible")
  p <- deep_update(p, pars)
  add_rows(contrasts_of(p, "reversible escape", load_json("calib_stage1j_reversible.json")$loss))
}
print(rows, digits = 3)
progress_note("07b_structural", paste(which_, collapse = ", "))
