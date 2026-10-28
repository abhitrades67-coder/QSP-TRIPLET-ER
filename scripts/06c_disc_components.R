# Components of AE-related discontinuation across the parameter ensemble (for the post hoc secondary tolerability
# criterion).
#
# Re-simulates the arms below for every accepted parameter set, with the same virtual population (n = 6000, seed
# 1001), grid and horizon as 06_insilico_trial.R ens, and records per arm and member the share of patients who
# stopped each drug for adverse events over the whole treatment period, the share who stopped the backbone (the
# CDK4/6 inhibitor or the endocrine agent; a stopped endocrine agent is replaced by another one in the model) and the
# share who stopped everolimus only. The any-drug share must equal ensemble_summary.csv exactly (same random
# streams); the script stops if it does not.
# Writes results/ensemble_disc_components.csv
# Usage: Rscript scripts/06c_disc_components.R

if (!exists("ROOT")) ROOT <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), ".."), winslash = "/")
TRIAL_AS_LIBRARY <- TRUE
source(file.path(ROOT, "scripts", "06_insilico_trial.R"))   # setup, all_arms(), NENS, H
DISC_ARMS <- c("A", "A_let", "A_eve5", "A_eve2.5", "A_tox1.5", "B", "B_trial", "C", "D1", "D2", "D3")

disc_components <- function(R) {
  cdk <- R$ae_disc_cdk; eve <- R$ae_disc_eve; et <- R$ae_disc_et
  c(any = mean(cdk | eve | et), cdk = mean(cdk), eve = mean(eve), et = mean(et), backbone = mean(cdk | et),
    eve_only = mean(eve & !cdk & !et))
}

E <- load_json("ensemble.json")$members
p0 <- load_params(PFILES)
arms <- all_arms()
es <- utils::read.csv(file.path(RES, "ensemble_summary.csv"), check.names = FALSE, stringsAsFactors = FALSE)
t0 <- Sys.time()
rows <- list()
for (j in seq_along(E)) {
  p <- deep_update(p0, E[[j]]$params)
  v <- make_vpop(p, NENS, 1001)
  grid <- make_grid(p)
  for (k in DISC_ARMS) {
    dc <- disc_components(simulate(p, v, grid, arms[[k]], H))
    any0 <- es$ae_disc_any[es$arm == k & es$member == j - 1L]
    if (length(any0) != 1 || abs(dc[["any"]] - any0) > 1e-12)
      stop(sprintf("member %d, arm %s: any-drug discontinuation %.6f differs from ensemble_summary.csv (%s)", j - 1L, k,
                   dc[["any"]], paste(any0, collapse = ", ")))
    rows[[length(rows) + 1]] <- data.frame(arm = k, member = j - 1L, t(dc))
  }
  if (j %% 10 == 0) cat(sprintf("  member %d / %d (%.0f s)\n", j, length(E), as.numeric(difftime(Sys.time(), t0, units = "secs"))))
}
utils::write.csv(do.call(rbind, rows), file.path(RES, "ensemble_disc_components.csv"), row.names = FALSE)
cat(sprintf("disc components: %d arms x %d members in %.0f s; any-drug shares identical to ensemble_summary.csv\n",
            length(DISC_ARMS), length(E), as.numeric(difftime(Sys.time(), t0, units = "secs"))))
progress_note("06c_disc_components", sprintf("%d arms x %d ensemble members; any-drug shares identical to ensemble_summary.csv",
                                             length(DISC_ARMS), length(E)))
