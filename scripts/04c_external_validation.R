# External validation with the calibrated parameters: SONIA (Sonke, Nature 2024) and its ctDNA analysis (Jongbloed,
# Nat Med 2025), CDK4/6i in line 1 vs line 2, and overall survival in persevERA (ASCO 2026, LBA1006), line 1
# giredestrant vs letrozole, each with palbociclib. Both trials were published before the model was built. persevERA PFS
# and discontinuation are calibration targets: they are listed with the criterion "calibrated" and not counted.
# Usage: Rscript scripts/04c_external_validation.R [params_dir]   (params_dir: the params_*.json to use; default results/)
# Writes results/external_validation.csv and results/tables/external_validation.md

if (!exists("ROOT")) ROOT <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), ".."), winslash = "/")
source(file.path(ROOT, "R", "setup.R"))
args <- commandArgs(trailingOnly = TRUE)
PDIR <- if (length(args) >= 1) args[1] else RES
NP <- as.integer(Sys.getenv("NPRO", "6000"))

p <- load_params(PFILES, dir = PDIR)
v <- make_vpop(p, NP, 61)
grid <- make_grid(p)
H <- 2190L
arms <- c(persevera_arms(), sonia_arms())
res <- lapply(arms, function(a) simulate(p, v, grid, a, H))

# administrative censoring as in each trial (uniform follow-up between the shortest and longest)
set.seed(62)
U <- stats::runif(NP)
fu_sonia <- (15 + 46 * U) * MONTH      # SONIA: 46 months of accrual, then 15 months to the data cut-off
fu_pv <- (40 + 24 * U) * MONTH         # persevERA: median follow-up 52.2 months
cens <- function(t, e, fu) list(t = pmin(t, fu), e = as.logical(e) & (t <= fu))
mo <- function(x) x / MONTH
hr_of <- function(a, b, ep, fu, mk = rep(TRUE, NP)) {
  x <- cens(a[[ep]][mk], a[[paste0(ep, "_ev")]][mk], fu[mk]); y <- cens(b[[ep]][mk], b[[paste0(ep, "_ev")]][mk], fu[mk])
  cox_hr(x$t, x$e, y$t, y$e)
}
med_of <- function(a, ep, fu, mk = rep(TRUE, NP)) {
  x <- cens(a[[ep]][mk], a[[paste0(ep, "_ev")]][mk], fu[mk]); mo(km_median(x$t, x$e))
}
lm_of <- function(a, ep, fu, m) { x <- cens(a[[ep]], a[[paste0(ep, "_ev")]], fu); 100 * surv_at(x$t, x$e, m * MONTH) }

rows <- list()
add <- function(trial, item, obs, ci, model, rule, pass) {
  rows[[length(rows) + 1]] <<- data.frame(trial = trial, item = item, published = obs, published_ci = ci,
                                          model = model, criterion = rule, met = pass, stringsAsFactors = FALSE)
}
yn <- function(x) if (is.na(x)) "" else if (x) "yes" else "NO"
f2 <- function(x) sprintf("%.2f", x)
f1 <- function(x) sprintf("%.1f", x)
inci <- function(x, lo, hi) x >= lo && x <= hi
within <- function(x, ref, tol) abs(x / ref - 1) <= tol

# ---- persevERA
g <- res$PV_gir; l <- res$PV_let
h <- hr_of(g, l, "pfs1", fu_pv)
crit <- function(rule) "calibrated"   # persevERA PFS and discontinuation are calibration targets (rule not applied)
met <- function(x) ""
add("persevERA", "PFS HR, giredestrant vs letrozole (+ palbociclib)", "0.89", "0.76-1.05", f2(h[1]),
    crit("within observed 95% CI"), met(yn(inci(h[1], 0.76, 1.05))))
mg <- med_of(g, "pfs1", fu_pv); ml <- med_of(l, "pfs1", fu_pv)
add("persevERA", "median PFS, giredestrant + palbociclib (mo)", "33.1", "30.2-38.3", f1(mg), crit("within +/-20%"), met(yn(within(mg, 33.1, 0.2))))
add("persevERA", "median PFS, letrozole + palbociclib (mo)", "28.2", "25.0-33.1", f1(ml), crit("within +/-20%"), met(yn(within(ml, 28.2, 0.2))))
for (m in c(12, 24, 36)) {
  ob <- c(`12` = "77.8 vs 77.3", `24` = "59.7 vs 57.3", `36` = "45.8 vs 41.9")[[as.character(m)]]
  add("persevERA", sprintf("PFS at %d mo, giredestrant vs letrozole (%%; presentation, two conference reports)", m), ob, "",
      sprintf("%.1f vs %.1f", lm_of(g, "pfs1", fu_pv, m), lm_of(l, "pfs1", fu_pv, m)), crit("descriptive"), "")
}
og <- 100 * mean(g$resp1); ol <- 100 * mean(l$resp1)
add("persevERA", "ORR difference, giredestrant - letrozole (points)", "1.4 (60.2 vs 58.8)", "", sprintf("%.1f (%.1f vs %.1f)", og - ol, og, ol),
    crit("|difference - 1.4| <= 10 points"), met(yn(abs((og - ol) - 1.4) <= 10)))
ho <- hr_of(g, l, "os", fu_pv)
add("persevERA", "OS HR", "1.03", "0.83-1.28", f2(ho[1]), "within observed 95% CI", yn(inci(ho[1], 0.83, 1.28)))
add("persevERA", "endocrine therapy stopped for AEs, giredestrant / letrozole (%)", "6.5 / 5.5", "",
    sprintf("%.1f / %.1f", 100 * mean(g$ae_disc_et), 100 * mean(l$ae_disc_et)), crit("each within 3 points"),
    met(yn(abs(100 * mean(g$ae_disc_et) - 6.5) <= 3 && abs(100 * mean(l$ae_disc_et) - 5.5) <= 3)))

# ---- SONIA
a <- res$SO_first; b <- res$SO_second
h1 <- hr_of(a, b, "pfs1", fu_sonia)
add("SONIA", "PFS1 HR, CDK4/6i first vs second", "0.59", "0.51-0.69", f2(h1[1]), "within observed 95% CI", yn(inci(h1[1], 0.51, 0.69)))
m1a <- med_of(a, "pfs1", fu_sonia); m1b <- med_of(b, "pfs1", fu_sonia)
add("SONIA", "median PFS1, first / second (mo)", "24.7 / 16.1", "", sprintf("%.1f / %.1f", m1a, m1b), "each within +/-20%",
    yn(within(m1a, 24.7, 0.2) && within(m1b, 16.1, 0.2)))
h2 <- hr_of(a, b, "edct", fu_sonia)
add("SONIA", "PFS2 HR (model: end of endocrine-based therapy)", "0.87", "0.74-1.03", f2(h2[1]), "within observed 95% CI", yn(inci(h2[1], 0.74, 1.03)))
m2a <- med_of(a, "edct", fu_sonia); m2b <- med_of(b, "edct", fu_sonia)
add("SONIA", "median PFS2, first / second (mo)", "31.0 / 26.8", "", sprintf("%.1f / %.1f", m2a, m2b), "each within +/-20%",
    yn(within(m2a, 31.0, 0.2) && within(m2b, 26.8, 0.2)))
h3 <- hr_of(a, b, "os", fu_sonia)
add("SONIA", "OS HR", "0.98", "0.80-1.20", f2(h3[1]), "within observed 95% CI", yn(inci(h3[1], 0.80, 1.20)))
dca <- mo(stats::median(a$days_cdk[a$days_cdk > 0])); dcb <- mo(stats::median(b$days_cdk[b$days_cdk > 0]))
add("SONIA", "median CDK4/6i treatment duration, first / second (mo)", "24.6 / 8.1", "", sprintf("%.1f / %.1f", dca, dcb),
    "each within +/-30%", yn(within(dca, 24.6, 0.3) && within(dcb, 8.1, 0.3)))
add("SONIA", "CDK4/6i stopped for AEs, CDK4/6i-first group (%)", "9 (while continuing endocrine therapy)", "",
    f1(100 * mean(a$ae_disc_cdk)), "between 6 and 15%", yn(inci(100 * mean(a$ae_disc_cdk), 6, 15)))
# ctDNA: high aneuploidy score (about >= 10% tumor DNA) in 34%; model tumor fraction at baseline
TF0 <- v$N0 / (v$N0 + p$tum$K_TF)
hi10 <- TF0 >= 0.10
add("SONIA ctDNA", "baseline tumor DNA fraction >= 10% (%)", "34", "", f1(100 * mean(hi10)), "descriptive", "")
for (def in c("TF >= 10%", "top 34% of TF")) {
  hm <- if (def == "TF >= 10%") hi10 else TF0 >= stats::quantile(TF0, 0.66)
  hh <- hr_of(a, b, "edct", fu_sonia, hm); hl <- hr_of(a, b, "edct", fu_sonia, !hm)
  add("SONIA ctDNA", sprintf("PFS2 HR, high tumor DNA (%s)", def), "0.58", "0.38-0.88", f2(hh[1]), "within observed 95% CI", yn(inci(hh[1], 0.38, 0.88)))
  add("SONIA ctDNA", sprintf("PFS2 HR, low tumor DNA (%s)", if (def == "TF >= 10%") "TF < 10%" else "bottom 66% of TF"), "1.36", "0.95-1.96", f2(hl[1]), "within observed 95% CI", yn(inci(hl[1], 0.95, 1.96)))
  add("SONIA ctDNA", sprintf("benefit of CDK4/6i first larger with high tumor DNA (%s)", def), "0.58 vs 1.36 (ratio 0.43)", "",
      sprintf("%s vs %s (ratio %s)", f2(hh[1]), f2(hl[1]), f2(hh[1] / hl[1])), "HR(high) / HR(low) <= 0.80",
      yn(hh[1] / hl[1] <= 0.80))
}

tab <- do.call(rbind, rows)
utils::write.csv(tab, file.path(RES, "external_validation.csv"), row.names = FALSE)
md <- c("| trial | item | published | 95% CI | model | criterion | met |", "|---|---|---|---|---|---|---|",
        sprintf("| %s | %s | %s | %s | %s | %s | %s |", tab$trial, tab$item, tab$published, tab$published_ci, tab$model,
                tab$criterion, tab$met))
dir.create(file.path(RES, "tables"), showWarnings = FALSE)
writeLines(md, file.path(RES, "tables", "external_validation.md"))
cat(md, sep = "\n")
nm <- sum(tab$met == "yes"); nc <- sum(tab$met %in% c("yes", "NO"))
progress_note("04c_external_validation", sprintf("%d of %d external-validation criteria met (persevERA OS, SONIA)", nm, nc))
