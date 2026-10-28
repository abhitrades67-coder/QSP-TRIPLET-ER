# Full QSP-TRIPLET-ER pipeline. Approximate wall time on a 6-core / 12-thread PC: about 7 hours.
# Run from a terminal:   Rscript run_all.R            (all steps)
#                        Rscript run_all.R 05 06      (only the steps whose names start with 05 or 06)
# Each step runs in its own R process; output goes to logs/<step>.log and progress to PROGRESS.md.
# Threads: set the environment variable QSP_THREADS (default: all logical cores).
# Calibration order: the joint refit of stages 1 and 3 starts from the values in calibration_start/ (03_start), and
# stages 2 and 4 are fitted once before the local refinement of the joint fit and again after it.
# Outputs: results/ (numbers, tables and Supplementary Data 1), results/article_tables/, figures/ (supplementary
# figures) and figures/article/ (main figures).

ROOT <- normalizePath(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), winslash = "/")
RSCRIPT <- file.path(R.home("bin"), "Rscript")
STEPS <- list(
  c("02_calibrate_tox", "scripts/02_calibrate_tox.R"),           # stage 0: toxicity hazards and AE discontinuation
  c("03_start", "scripts/03a_start_state.R"),                    # start state of the efficacy calibration; clears structural.csv
  c("03_stage13", "scripts/03_calibrate_efficacy.R stage13 30"), # joint refit of stages 1 and 3
  c("03_stage2a", "scripts/03_calibrate_efficacy.R stage2 15"),  # everolimus effect, first pass (used by 03b)
  c("03_stage4a", "scripts/03_calibrate_efficacy.R stage4"),     # survival, first pass
  c("03b_refine", "scripts/03b_refine_stage13.R 12"),            # local (compass) refinement of the joint fit
  c("03_stage2", "scripts/03_calibrate_efficacy.R stage2 15"),   # everolimus effect
  c("03_stage4", "scripts/03_calibrate_efficacy.R stage4"),      # survival (post-endocrine hazard and therapy effect)
  c("01_verify", "scripts/01_verify.R"),                         # code / numerical verification
  c("04b_calibration_check", "scripts/04b_calibration_check.R"), # calibration targets, independent population
  c("04_validate", "scripts/04_validate.R"),                     # held-out validation
  c("04c_external_validation", "scripts/04c_external_validation.R"), # external validation: SONIA, persevERA OS
  c("04d_disc_check", "scripts/04d_disc_check.R"),               # AE discontinuation in full-model trial emulations
  c("05_ensemble", "scripts/05_ensemble.R 360"),                 # parameter ensemble
  c("06_ref", "scripts/06_insilico_trial.R ref"),                # reference in silico trial (n = 6000 per arm)
  c("06_ens", "scripts/06_insilico_trial.R ens"),                # ensemble in silico trial (n = 6000 per arm and member)
  c("06b_disc", "scripts/06b_discontinuation_sensitivity.R"),    # alternative AE-discontinuation models
  c("06c_disc", "scripts/06c_disc_components.R"),                # discontinuation by drug across the ensemble (secondary criterion)
  c("07_sensitivity", "scripts/07_sensitivity.R 160"),           # PRCC + EVPPI
  c("07b_structural", "scripts/07b_structural.R"),               # structural alternatives (recalibrated)
  c("08_power", "scripts/08_power.R"),                           # virtual trial power
  c("09_figures", "scripts/09_figures.R"),                       # figures
  c("10_tables", "scripts/10_tables.R"),                         # tables and headline numbers
  c("11_article_figures", "scripts/11_article_figures.R"),       # main figures of the article (QSP_MS_FIG, QSP_MS_PNG)
  c("12_article_tables", "scripts/12_article_tables.R")          # article tables and Supplementary Data 1 (QSP_MS_TAB, QSP_SUB)
)
sel <- commandArgs(trailingOnly = TRUE)
dir.create(file.path(ROOT, "logs"), showWarnings = FALSE)
for (s in STEPS) {
  if (length(sel) && !any(startsWith(s[1], sel))) next
  parts <- strsplit(s[2], " ")[[1]]
  log <- file.path(ROOT, "logs", paste0(s[1], ".log"))
  cat(sprintf("[%s] %s ... ", format(Sys.time(), "%H:%M:%S"), s[1]))
  t0 <- Sys.time()
  rc <- system2(RSCRIPT, c(shQuote(file.path(ROOT, parts[1])), parts[-1]), stdout = log, stderr = log)
  cat(sprintf("%s (%.1f min)\n", if (rc == 0) "ok" else paste("FAILED, see", log),
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  if (rc != 0) stop("step failed: ", s[1])
}
