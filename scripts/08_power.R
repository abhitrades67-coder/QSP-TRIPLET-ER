# Virtual clinical trials of realistic size drawn from the reference virtual population.
# Design: 1:1 randomization, uniform accrual over 24 months, minimum follow-up
# 36 months (administrative censoring), independent patients per arm. Endpoints: EDCT (log-rank and
# RMST at 36 months), OS (log-rank). 1000 replicates per design. Designs include trial-like attrition in Arm B
# (15%; 10% and 20%; SONIA) and two enriched designs.
# Writes results/power.csv

if (!exists("ROOT")) ROOT <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), ".."), winslash = "/")
source(file.path(ROOT, "R", "setup.R"))
REPS <- as.integer(Sys.getenv("REPS", "1000"))

# two-sample log-rank test, p-value (checked against survival::survdiff below)
logrank_p <- function(ta, ea, tb, eb) {
  t <- c(ta, tb); e <- c(ea, eb); g <- c(rep(1L, length(ta)), rep(0L, length(tb)))
  o <- order(t); t <- t[o]; e <- e[o]; g <- g[o]
  n <- length(t)
  first <- which(c(TRUE, diff(t) != 0)); last <- c(first[-1] - 1L, n)
  n_risk <- n - first + 1
  n1_risk <- rev(cumsum(rev(g)))[first]
  ce <- cumsum(e); ce1 <- cumsum(e * g)
  d <- ce[last] - c(0, ce[last][-length(last)])
  d1 <- ce1[last] - c(0, ce1[last][-length(last)])
  k <- d > 0
  E1 <- sum((n1_risk * d / n_risk)[k])
  V <- sum((n1_risk * (n_risk - n1_risk) * d * (n_risk - d) / (n_risk^2 * pmax(n_risk - 1, 1)))[k])
  stats::pchisq((sum(d1) - E1)^2 / V, 1, lower.tail = FALSE)
}

ref <- readRDS(file.path(RES, "ref_arms.rds"))
arms <- ref[c("A", "B", "C")]
n <- length(arms$A$os)
# check of the log-rank implementation
chk <- local({
  h <- min(300L, n %/% 2)
  ta <- arms$A$edct[1:h]; tb <- arms$B$edct[(h + 1):(2 * h)]
  ea <- arms$A$edct_ev[1:h]; eb <- arms$B$edct_ev[(h + 1):(2 * h)]
  sd_ <- survival::survdiff(survival::Surv(c(ta, tb), c(ea, eb)) ~ rep(1:2, each = h))
  abs(logrank_p(ta, ea, tb, eb) - stats::pchisq(sd_$chisq, 1, lower.tail = FALSE))
})
cat(sprintf("log-rank check vs survival::survdiff: |p difference| = %.2e\n", chk))
save_json(list(logrank_p_difference = chk), "power_check.json")   # reported in Supplementary Table 21
set.seed(12)
SIZES <- c(100L, 150L, 200L, 300L, 400L, 600L, 800L, 1200L, 1600L, 2400L)
f <- readRDS(file.path(RES, "ref_features.rds"))
# designs: contrast label, arm 1, arm 2, patients eligible (indices into the reference population)
DESIGNS <- list(
  list("A vs B", "A", "B", seq_len(n)),
  list("A vs C", "A", "C", seq_len(n)),
  list("B vs C", "B", "C", seq_len(n)),
  # a randomized trial in which Arm B patients who stop line 1 mostly start the protocol line 2
  list("A vs B, trial-like attrition 15%", "A", "B_trial", seq_len(n)),
  list("A vs B, trial-like attrition 10%", "A", "B_trial_10", seq_len(n)),
  list("A vs B, trial-like attrition 20%", "A", "B_trial_20", seq_len(n)),
  # enriched designs (trial-like attrition 15%): high tumor DNA fraction; high proliferation tercile
  list("A vs B, trial-like 15%, tumor DNA >= 10% only", "A", "B_trial", which(f$TF >= 0.10)),
  list("A vs B, trial-like 15%, proliferation high tercile only", "A", "B_trial",
       which(f$kmax > stats::quantile(f$kmax, 2 / 3)))
)
rows <- list()
for (d in DESIGNS) {
  a <- ref[[d[[2]]]]; b <- ref[[d[[3]]]]; pool <- d[[4]]
  for (npa in SIZES) {
    if (2 * npa > length(pool)) next   # not enough virtual patients for this design
    hits <- c(EDCT_logrank = 0, EDCT_RMST36 = 0, OS_logrank = 0)
    for (r in seq_len(REPS)) {
      perm <- pool[sample.int(length(pool))]
      ia <- perm[1:npa]; ib <- perm[(npa + 1):(2 * npa)]
      fu <- (36 + stats::runif(2 * npa, 0, 24)) * MONTH
      fa <- fu[1:npa]; fb <- fu[(npa + 1):(2 * npa)]
      for (ep in c("edct", "os")) {
        ta <- a[[ep]][ia]; tb <- b[[ep]][ib]
        ea <- a[[paste0(ep, "_ev")]][ia] & (ta <= fa)
        eb <- b[[paste0(ep, "_ev")]][ib] & (tb <= fb)
        ta <- pmin(ta, fa); tb <- pmin(tb, fb)
        key <- if (ep == "edct") "EDCT_logrank" else "OS_logrank"
        hits[key] <- hits[key] + (logrank_p(ta, ea, tb, eb) < 0.05)
        if (ep == "edct") {
          xa <- pmin(ta, 36 * MONTH) / MONTH; xb <- pmin(tb, 36 * MONTH) / MONTH
          z <- (mean(xa) - mean(xb)) / sqrt(stats::var(xa) / npa + stats::var(xb) / npa)
          hits["EDCT_RMST36"] <- hits["EDCT_RMST36"] + (abs(z) > stats::qnorm(0.975))
        }
      }
    }
    gain <- mean(pmin(a$edct[pool], 48 * MONTH) - pmin(b$edct[pool], 48 * MONTH)) / MONTH
    for (k in names(hits)) rows[[length(rows) + 1]] <- data.frame(contrast = d[[1]], n_per_arm = npa, endpoint = k,
                                                                  power = hits[[k]] / REPS, gain_EDCT48 = gain,
                                                                  eligible = length(pool))
    cat(d[[1]], npa, round(hits / REPS, 3), "\n")
  }
}
utils::write.csv(do.call(rbind, rows), file.path(RES, "power.csv"), row.names = FALSE)
progress_note("08_power", sprintf("%d designs (real-world and trial-like attrition, enriched)", length(DESIGNS)))
