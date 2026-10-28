# AE discontinuation in full-model trial emulations (treatment until progression; common follow-up = the
# trial's median follow-up), compared with the published rates. A check of the discontinuation sub-model
# outside its toxicity-only calibration; some rates were fitted (marked), the others were not.
# Usage: Rscript scripts/04d_disc_check.R [params_dir]   (params_dir: the params_*.json to use; default results/)
# Writes results/disc_check.csv and results/tables/disc_check.md

if (!exists("ROOT")) ROOT <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), ".."), winslash = "/")
source(file.path(ROOT, "R", "setup.R"))
args <- commandArgs(trailingOnly = TRUE)
PDIR <- if (length(args) >= 1) args[1] else RES
NP <- as.integer(Sys.getenv("NDISC", "6000"))

p <- load_params(PFILES, dir = PDIR)
v <- make_vpop(p, NP, 71)
grid <- make_grid(p)
pv <- persevera_arms(); so <- sonia_arms()
ev_exp <- trial_arm(list(R_("palbo", "letrozole", NULL), R_(NULL, "giredestrant", 10)))
CHK <- list(
  list("PALOMA-2", "palbociclib + letrozole: any study drug", list(lines = list(R_("palbo", "letrozole")), ass1 = SCHED$q12w), 23, "any", 9.7, "yes"),
  list("PALOMA-2", "palbociclib + letrozole: palbociclib", list(lines = list(R_("palbo", "letrozole")), ass1 = SCHED$q12w), 23, "cdk", 7.4, "yes"),
  list("PALOMA-2", "letrozole + placebo: any study drug", list(lines = list(R_(NULL, "letrozole")), ass1 = SCHED$q12w), 23, "any", 5.9, "yes"),
  list("MONALEESA-2", "ribociclib + letrozole: any study drug (published: both drugs stopped)", list(lines = list(R_("ribo", "letrozole"))), 15.3, "any", 7.5, "no"),
  list("MONARCH 3", "abemaciclib + aromatase inhibitor: abemaciclib", list(lines = list(R_("abema", "letrozole"))), 17.8, "cdk", 19.6, "yes"),
  list("BOLERO-4", "everolimus + letrozole (line 1): any drug", list(lines = list(R_(NULL, "letrozole", 10)), ass1 = SCHED$q8w, tox_over = NODEX), 29.5, "any", 25.2, "yes"),
  list("persevERA", "giredestrant (+ palbociclib)", pv$PV_gir, 52.2, "et", 6.5, "yes"),
  list("persevERA", "letrozole (+ palbociclib)", pv$PV_let, 52.2, "et", 5.5, "yes"),
  list("evERA", "giredestrant + everolimus (line 2, patients who started it): everolimus", ev_exp, NA, "eve", 17.0, "no", "line2"),
  list("SONIA", "CDK4/6i-first group: CDK4/6i (9% stopped it while continuing ET)", so$SO_first, 37.3, "cdk", 9.0, "no")
)
rows <- lapply(CHK, function(k) {
  # line-2 regimens: patients who started line 2, followed until it ends (the rate refers to the regimen)
  l2 <- length(k) >= 8 && identical(k[[8]], "line2")
  R <- simulate(p, v, grid, k[[3]], if (l2) 2190L else as.integer(round(k[[4]] * MONTH)))
  mk <- if (l2) R$in_line2 else rep(TRUE, length(R$os))
  m <- 100 * switch(k[[5]], any = mean((R$ae_disc_cdk | R$ae_disc_eve | R$ae_disc_et)[mk]), cdk = mean(R$ae_disc_cdk[mk]),
                    eve = mean(R$ae_disc_eve[mk]), et = mean(R$ae_disc_et[mk]))
  tol <- max(3, 0.2 * k[[6]])
  fitted <- k[[7]]
  data.frame(trial = k[[1]], item = k[[2]], follow_up_mo = if (l2) "to end of line 2" else k[[4]], published_pct = k[[6]], model_pct = round(m, 1),
             fitted = fitted, criterion = sprintf("within %.0f points", tol), met = if (abs(m - k[[6]]) <= tol) "yes" else "NO",
             stringsAsFactors = FALSE)
})
tab <- do.call(rbind, rows)
utils::write.csv(tab, file.path(RES, "disc_check.csv"), row.names = FALSE)
md <- c("| trial | stopped for adverse events | follow-up (mo) | published % | model % | fitted in toxicity calibration | criterion | met |",
        "|---|---|---|---|---|---|---|---|",
        sprintf("| %s | %s | %s | %s | %s | %s | %s | %s |", tab$trial, tab$item, tab$follow_up_mo, tab$published_pct,
                tab$model_pct, tab$fitted, tab$criterion, tab$met))
writeLines(md, file.path(RES, "tables", "disc_check.md"))
cat(md, sep = "\n")
progress_note("04d_disc_check", sprintf("%d of %d discontinuation checks met (full-model emulations)",
                                                       sum(tab$met == "yes"), nrow(tab)))
