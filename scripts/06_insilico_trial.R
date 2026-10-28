# In silico trial.
#
# Part 1 (ref): best-fit parameters, large virtual population, all arms and variants, paired
#               (digital-twin) design. Per-patient outputs -> results/ref_arms.rds
# Part 2 (ens): every accepted ABC parameter set, core arms; summary metrics and paired differences
#               per set -> results/ensemble_summary.csv, ensemble_subgroups.csv, ensemble_curves.rds
# Usage: Rscript scripts/06_insilico_trial.R ref|ens

if (!exists("ROOT")) ROOT <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), ".."), winslash = "/")
source(file.path(ROOT, "R", "setup.R"))
NREF <- as.integer(Sys.getenv("NREF", "6000"))
NENS <- as.integer(Sys.getenv("NENS", "6000"))
H <- 2190L

all_arms <- function() {
  a <- strategy_arms("palbo")
  a[["A_eve5"]] <- list(lines = list(regimen("palbo", "giredestrant", 5)))
  # the triplet with letrozole instead of giredestrant (persevERA found little first-line benefit of giredestrant over
  # letrozole)
  a[["A_let"]] <- list(lines = list(regimen("palbo", "letrozole", 10)))
  a[["A_eve2.5"]] <- list(lines = list(regimen("palbo", "giredestrant", 2.5)))
  a[["A_ribo_eve2.5"]] <- list(lines = list(regimen("ribo", "giredestrant", 2.5)))
  a[["A_ribo_eve10"]] <- list(lines = list(regimen("ribo", "giredestrant", 10)))
  a[["A_ribo300_eve2.5"]] <- list(lines = list(regimen("ribo", "giredestrant", 2.5, c(300, 200))))
  a[["A_abema_eve10"]] <- list(lines = list(regimen("abema", "giredestrant", 10)))
  a[["A_abema_eve5"]] <- list(lines = list(regimen("abema", "giredestrant", 5)))
  for (cdk in c("ribo", "abema")) {
    s <- strategy_arms(cdk)
    a[[paste0("B_", cdk)]] <- s$B
    a[[paste0("C_", cdk)]] <- s$C
  }
  a[["A_tox1.5"]] <- list(lines = list(regimen("palbo", "giredestrant", 10)), tox_over = list(tox_int = 1.5))
  a[["A_tox2.0"]] <- list(lines = list(regimen("palbo", "giredestrant", 10)), tox_over = list(tox_int = 2.0))
  a[["A_eve5_tox2.0"]] <- list(lines = list(regimen("palbo", "giredestrant", 5)), tox_over = list(tox_int = 2.0))
  a[["A_nodex"]] <- list(lines = list(regimen("palbo", "giredestrant", 10)), tox_over = NODEX)
  a[["B_nodex"]] <- c(strategy_arms("palbo")$B, list(tox_over = NODEX))
  a[["A_dexstrict"]] <- list(lines = list(regimen("palbo", "giredestrant", 10)), tox_over = DEX_STRICT)
  a[["B_dexstrict"]] <- c(strategy_arms("palbo")$B, list(tox_over = DEX_STRICT))
  a[["A_eve5_dexstrict"]] <- list(lines = list(regimen("palbo", "giredestrant", 5)), tox_over = DEX_STRICT)
  a[["C_ful"]] <- list(lines = list(regimen("palbo", "letrozole"), regimen(NULL, "fulvestrant", 10)))
  a[["B_ful1L"]] <- list(lines = list(regimen("palbo", "fulvestrant"), regimen(NULL, "giredestrant", 10)))
  a[["C_ful1L"]] <- list(lines = list(regimen("palbo", "fulvestrant"), regimen(NULL, "exemestane", 10)))
  # trial-like attrition: Arm B / C patients who stop line 1 mostly start the protocol line 2
  for (k in names(ATTR_TRIAL)) {
    nm <- if (k == "base") "" else paste0("_", round(100 * ATTR_TRIAL[[k]]))
    a[[paste0("B_trial", nm)]] <- c(strategy_arms("palbo")$B, list(p_attr = ATTR_TRIAL[[k]], p_attr_mid = ATTR_TRIAL[[k]]))
  }
  a[["C_trial"]] <- c(strategy_arms("palbo")$C, list(p_attr = ATTR_TRIAL[["base"]], p_attr_mid = ATTR_TRIAL[["base"]]))
  s2 <- scenario2_arms()
  for (k in names(s2)) a[[paste0("S2_", k)]] <- s2[[k]]
  a
}

ENS_ARMS <- c("A", "B", "C", "D1", "D2", "D3", "A_eve5", "A_eve2.5", "A_tox1.5", "S2_A'(abema)", "S2_B'", "S2_C'",
              "B_trial", "A_let")

features <- function(p, v) {
  t <- p$tum
  TF0 <- v$N0 / (v$N0 + t$K_TF)
  vaf0 <- 0.5 * v$fE0 * TF0
  data.frame(b0 = v$b0, Ptrunc = v$Ptrunc, fE0 = v$fE0, fR0 = v$fR0, kmax = v$kmax, N0 = v$N0,
             expo_eve = v$expo_eve, expo_palbo = v$expo_palbo, frailty = v$tox_frailty, age = v$age, rd = v$rd,
             esr1_base = vaf0 >= t$LOD, TF = TF0)
}

subgroup_masks <- function(f) {
  q <- function(x) stats::quantile(x, c(1 / 3, 2 / 3), names = FALSE)
  b <- q(f$b0); n0 <- q(f$N0); k <- q(f$kmax); fr <- q(f$frailty)
  list(
    "All" = rep(TRUE, nrow(f)),
    "PI3K pathway altered" = f$Ptrunc, "PI3K pathway wild-type" = !f$Ptrunc,
    "Baseline ESR1 detectable" = f$esr1_base, "Baseline ESR1 not detectable" = !f$esr1_base,
    "Endocrine dependence high (b0 low tercile)" = f$b0 <= b[1],
    "Endocrine dependence mid" = f$b0 > b[1] & f$b0 <= b[2],
    "Endocrine dependence low (b0 high tercile)" = f$b0 > b[2],
    "Tumor burden low tercile" = f$N0 <= n0[1], "Tumor burden high tercile" = f$N0 > n0[2],
    "Proliferation low tercile" = f$kmax <= k[1], "Proliferation high tercile" = f$kmax > k[2],
    "Tumor DNA fraction >= 10%" = f$TF >= 0.10, "Tumor DNA fraction < 10%" = f$TF < 0.10,
    "Age >= 65" = f$age >= 65, "Age < 65" = f$age < 65,
    "AE susceptibility high tercile" = f$frailty > fr[2]
  )
}

keep_fields <- function(R) R[intersect(c(KEEP, "days_cdk", "days_eve", "days_gir"), names(R))]

run_ref <- function() {
  p <- load_params(PFILES)
  v <- make_vpop(p, NREF, 1001)
  grid <- make_grid(p)
  arms <- all_arms()
  t0 <- Sys.time()
  res <- lapply(arms, function(a) simulate(p, v, grid, a, H))
  cat(sprintf("reference run: %d arms x %d patients in %.0f s\n", length(arms), NREF,
              as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  f <- features(p, v)
  saveRDS(f, file.path(RES, "ref_features.rds"))
  saveRDS(lapply(res, keep_fields), file.path(RES, "ref_arms.rds"))
  sm <- do.call(rbind, lapply(names(res), function(k) data.frame(arm = k, as.data.frame(summarize_arm(res[[k]])),
                                                                 check.names = FALSE)))
  utils::write.csv(sm, file.path(RES, "ref_summary.csv"), row.names = FALSE)
  print(sm[, c("arm", "EDCT_RMST48", "EDCT_median", "OS_RMST60", "PFS1_median", "ae_disc_any", "QTWiST")], digits = 3)
  # example patients (digital twins) for A, B, C: first 40 patients, same random streams as above
  m <- 40L
  vv <- vpop_head(v, m)
  trajs <- lapply(c(A = "A", B = "B", C = "C"), function(k) {
    R <- simulate(p, vv, grid, arms[[k]], H, record_traj = m)
    list(traj = R$traj, edct = R$edct, os = R$os, pfs1 = R$pfs1)
  })
  saveRDS(trajs, file.path(RES, "ref_traj.rds"))
  progress_note("06_insilico_trial ref", sprintf("%d arms x %d patients", length(arms), NREF))
}

run_ens <- function() {
  E <- load_json("ensemble.json")$members
  p0 <- load_params(PFILES)
  arms <- all_arms()
  grid_m <- (0:72) * MONTH
  t0 <- Sys.time()
  rows <- list(); srows <- list(); curves <- list()
  for (j in seq_along(E)) {
    p <- deep_update(p0, E[[j]]$params)
    v <- make_vpop(p, NENS, 1001)
    grid <- make_grid(p)
    res <- lapply(stats::setNames(ENS_ARMS, ENS_ARMS), function(k) simulate(p, v, grid, arms[[k]], H))
    for (k in ENS_ARMS) rows[[length(rows) + 1]] <- data.frame(arm = k, member = j - 1L,
                                                                as.data.frame(summarize_arm(res[[k]])), check.names = FALSE)
    for (k in c("A", "B", "C", "D1", "D2", "D3")) for (ep in c("edct", "os", "pfs1")) {
      key <- paste0(k, "__", ep)
      curves[[key]] <- rbind(curves[[key]], vapply(grid_m, function(g) mean(res[[k]][[ep]] > g), numeric(1)))
    }
    masks <- subgroup_masks(features(p, v))
    for (sg in names(masks)) {
      mk <- masks[[sg]]
      if (sum(mk) < 30) next
      for (cc in list(c("A", "B"), c("A", "C"), c("B", "C"), c("A", "B_trial"))) {
        a <- res[[cc[1]]]; b <- res[[cc[2]]]
        dE <- mean(pmin(a$edct[mk], 48 * MONTH) - pmin(b$edct[mk], 48 * MONTH)) / MONTH
        dO <- mean(pmin(a$os[mk], 60 * MONTH) - pmin(b$os[mk], 60 * MONTH)) / MONTH
        srows[[length(srows) + 1]] <- data.frame(member = j - 1L, subgroup = sg, contrast = paste0(cc[1], "-", cc[2]),
                                                 dEDCT48 = dE, dOS60 = dO, n = sum(mk))
      }
    }
    if (j %% 10 == 0) cat(sprintf("  member %d / %d (%.0f s)\n", j, length(E), as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  }
  utils::write.csv(do.call(rbind, rows), file.path(RES, "ensemble_summary.csv"), row.names = FALSE)
  utils::write.csv(do.call(rbind, srows), file.path(RES, "ensemble_subgroups.csv"), row.names = FALSE)
  saveRDS(curves, file.path(RES, "ensemble_curves.rds"))
  cat(sprintf("ensemble: %d members in %.0f s\n", length(E), as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  progress_note("06_insilico_trial ens", sprintf("%d ensemble members x %d patients", length(E), NENS))
}

if (!exists("TRIAL_AS_LIBRARY")) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) >= 1 && args[1] == "ref") run_ref()
  if (length(args) >= 1 && args[1] == "ens") run_ens()
}
