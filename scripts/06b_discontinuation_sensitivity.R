# Sensitivity of the outcomes to the AE-discontinuation sub-model.
#
# Base case: the 'other' discontinuation hazards fall from the first day of each line
# (time constant and floor fitted), everolimus discontinuation is proportional to everolimus exposure,
# and endocrine therapy has an extra hazard while everolimus is given; fitted to short- and long-exposure
# published rates (scripts/02_calibrate_tox.R). Alternatives:
#   const   constant hazards refitted to the same targets (results/disc_const.json)
#   gamma0  everolimus discontinuation independent of exposure (levels not refitted)
# Same reference population (n = 6000, seed 1001) as 06_insilico_trial.R ref.
# Writes results/disc_sensitivity.csv

if (!exists("ROOT")) ROOT <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), ".."), winslash = "/")
source(file.path(ROOT, "R", "setup.R"))
p <- load_params(PFILES)
v <- make_vpop(p, as.integer(Sys.getenv("NREF", "6000")), 1001)
grid <- make_grid(p)
s <- strategy_arms("palbo")
s$A_eve5 <- list(lines = list(regimen("palbo", "giredestrant", 5)))

VARIANTS <- list(base = NULL, gamma0 = list(disc_eve_gamma = 0))
if (file.exists(file.path(RES, "disc_const.json"))) {
  dc <- load_json("disc_const.json")$tox
  VARIANTS$const <- c(dc[c("disc_cdk", "disc_eve", "disc_et", "disc_gir", "disc_et_eve", "disc_abema_mult")],
                      list(disc_tau = NA_real_, disc_eve_gamma = dc$disc_eve_gamma))
}
arms <- list()
for (k in c("A", "A_eve5", "B", "C", "D2")) for (vn in names(VARIANTS))
  arms[[paste0(k, "_", vn)]] <- if (is.null(VARIANTS[[vn]])) s[[k]] else c(s[[k]], list(tox_over = VARIANTS[[vn]]))
res <- lapply(arms, function(a) simulate(p, v, grid, a, 2190L))
df <- do.call(rbind, lapply(names(res), function(k) data.frame(arm = k, as.data.frame(summarize_arm(res[[k]])),
                                                               check.names = FALSE)))
utils::write.csv(df, file.path(RES, "disc_sensitivity.csv"), row.names = FALSE)
print(df[, c("arm", "EDCT_RMST48", "OS_RMST60", "ae_disc_cdk", "ae_disc_eve", "ae_disc_et", "ae_disc_any", "RDI_eve",
             "QTWiST")], digits = 3)
progress_note("06b_discontinuation_sensitivity", sprintf("%d variants x %d arms", length(VARIANTS), 5L))
