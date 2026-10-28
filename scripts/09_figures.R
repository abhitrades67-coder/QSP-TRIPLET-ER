# Figures (PNG; resolution QSP_DPI, default 200 dpi) -> figures/ (or QSP_FIG): model behaviour, calibration, validation,
# outcomes, toxicity, sensitivity and the trial emulations, including the supplementary figures of the article.
# Usage: Rscript scripts/09_figures.R [figure numbers ...]   (default: all)

if (!exists("ROOT")) ROOT <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), ".."), winslash = "/")
source(file.path(ROOT, "R", "setup.R"))
suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })

COL <- c(A = "#c0392b", B = "#2471a3", C = "#7f8c8d", D1 = "#8e44ad", D2 = "#d68910", D3 = "#17a589")
LAB <- c(A = "A: upfront triplet", B = "B: sequential (gir+eve 2L)", C = "C: standard (ET+eve 2L)",
         D1 = "D1: CDK4/6i+gir → gir+eve", D2 = "D2: CDK4/6i+letrozole+eve → gir+eve", D3 = "D3: ctDNA-guided switch")
noai <- function(x) gsub("\\bAI\\b", "aromatase inhibitor", x, perl = TRUE)   # labels spell out the abbreviation
theme_q <- function(base = 9) theme_minimal(base_size = base) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(size = base, face = "bold"),
        legend.key.size = grid::unit(0.8, "lines"), strip.text = element_text(face = "bold"))
save_fig <- function(plot, name, w, h) {
  ggsave(file.path(FIG, name), plot, width = w, height = h, dpi = as.numeric(Sys.getenv("QSP_DPI", "200")), bg = "white")   # QSP_DPI: 300 for the manuscript
  cat("saved", name, "\n")
}
p <- load_params(PFILES)
km_df <- function(t, e, w = NULL, label = "", scale_t = MONTH) {
  k <- km(t, e, w)
  data.frame(time = c(0, k$t / scale_t), surv = 100 * c(1, k$s), label = label)
}
ref_arm <- local({ cache <- NULL; function(k) { if (is.null(cache)) cache <<- readRDS(file.path(RES, "ref_arms.rds")); cache[[k]] } })

# ---------------------------------------------------------------------------- 1
fig01 <- function() {
  boxes <- data.frame(
    x = c(0.1, 3.1, 6.9, 0.1, 4.8), y = c(3.4, 3.4, 3.4, 0.3, 0.3), w = c(2.6, 3.4, 3.0, 4.2, 5.1), h = c(2.4, 2.4, 2.4, 2.6, 2.6),
    fill = c("#fdebd0", "#d6eaf8", "#d5f5e3", "#fadbd8", "#e8daef"),
    title = c("1. PK and DDI", "2. Signaling (per clone)", "3. Tumor clones", "4. Toxicity and dose modification", "5. Outcomes"),
    body = c("palbociclib / ribociclib / abemaciclib\neverolimus, giredestrant\nletrozole, exemestane, fulvestrant\nCYP3A4 mechanism-based inactivation\n(ribociclib, palbociclib -> everolimus)\nestradiol suppression by AI",
             "E2 / SERD -> ER (WT or ESR1-mutant)\nER -> cyclin D1, IRS1; FOXO -> ESR1\nIRS1 -> PI3K/AKT -> mTORC1 -> S6K, 4E-BP1\nS6K -| IRS1 (feedback); S6K -> ER Ser167\ncyclin D1-CDK4/6 + cyclin E-CDK2 -> pRb\noutput: proliferation index P",
             "8 genotypes: ESR1 (E) x PI3K (P)\nx CDK4/6i escape (R)\nmutation + adaptive switching\n(adaptive switch dampened by mTORC1\ninhibition, weight kappa)\nRECIST SLD, new lesions, ctDNA VAF",
             "neutrophils: Friberg model (CDK4/6i)\nstomatitis, pneumonitis, hyperglycemia,\ndiarrhea, ALT/AST, QTc, bradycardia, anemia\nlabel-based holds / reductions / stops\nrealized dose feeds back to PK",
             "PFS1, line-2 PFS, endocrine-based disease\ncontrol time (EDCT), OS\ntumor-burden-dependent death hazard + background\npost-endocrine hazard and therapy (chemo/ADC)\nattrition before 2nd endocrine line (real-world)\nQ-TWiST, ESR1 emergence, clonal makeup"))
  g <- ggplot() + xlim(0, 10) + ylim(0, 6.2) + theme_void() +
    geom_rect(data = boxes, aes(xmin = x, xmax = x + w, ymin = y, ymax = y + h), fill = boxes$fill, colour = "#333333") +
    geom_text(data = boxes, aes(x = x + w / 2, y = y + h - 0.2, label = title), fontface = "bold", size = 3.9, vjust = 1) +
    geom_text(data = boxes, aes(x = x + w / 2, y = y + h - 0.6, label = body), size = 3.2, vjust = 1, lineheight = 1.05) +
    annotate("segment", x = c(2.7, 6.5, 8.4, 4.3), xend = c(3.1, 6.9, 8.4, 4.8), y = c(4.6, 4.6, 3.4, 1.6),
             yend = c(4.6, 4.6, 2.9, 1.6), arrow = grid::arrow(length = grid::unit(0.15, "cm")), linewidth = 0.6) +
    annotate("segment", x = c(1.2, 1.6), xend = c(1.2, 1.6), y = c(3.4, 2.9), yend = c(2.9, 3.4),
             arrow = grid::arrow(length = grid::unit(0.15, "cm")), linewidth = 0.6) +
    annotate("text", x = 5, y = 6.1, label = "QSP-TRIPLET-ER model structure", fontface = "bold", size = 4)
  save_fig(g, "fig01_model_schematic.png", 10, 5.8)
}

# ---------------------------------------------------------------------------- 2
fig02 <- function() {
  s <- p$sig; st <- p$st; ph <- p$ph
  Ie <- min(0.98, ph$I10_e * (ph$IC50_e + 1)) * 1 / (ph$IC50_e + 1)
  Ic <- ph$Imax_c / (1 + ph$IC50_c^ph$h_c)
  conds <- list("letrozole" = c(ph$E2_let, 0, 0, 0), "letrozole + CDK4/6i" = c(ph$E2_let, 0, 0, Ic),
                "exemestane + everolimus" = c(ph$E2_exe, 0, Ie, 0), "giredestrant + everolimus" = c(1, ph$uG_ref, Ie, 0),
                "triplet (CDK4/6i + gir + eve)" = c(1, ph$uG_ref, Ie, Ic))
  cols <- c("#2471a3", "#1abc9c", "#7f8c8d", "#d68910", "#c0392b")
  d <- NULL
  for (g in 0:1) for (nm in names(conds)) {
    r <- simulate_ode(conds[[nm]], 0.25, g, s, st, t_end = 168, n_out = 337)
    for (k in c("A_ER", "AKT", "S6K", "CycD", "P"))
      d <- rbind(d, data.frame(t = r$t, y = r[[k]], var = k, cond = nm, geno = if (g == 0) "wild-type ER" else "ESR1-mutant"))
  }
  d$var <- factor(d$var, levels = c("A_ER", "AKT", "S6K", "CycD", "P"),
                  labels = c("ER activity", "AKT (feedback)", "S6K (mTORC1)", "cyclin D1", "proliferation index"))
  d$cond <- factor(d$cond, levels = names(conds))
  d$geno <- factor(d$geno, levels = c("wild-type ER", "ESR1-mutant"))
  g <- ggplot(d, aes(t, y, colour = cond)) + geom_line(linewidth = 0.7) +
    facet_grid(geno ~ var, scales = "free_y") + scale_colour_manual(values = cols, name = NULL) +
    labs(x = "hours after drug start", y = NULL,
         title = "Signalling responses to each regimen (dynamic model from the untreated steady state; b₀ = 0.25)") +
    theme_q() + theme(legend.position = "bottom")
  save_fig(g, "fig02_signaling_timecourses.png", 13, 5.6)
}

# ---------------------------------------------------------------------------- 3
fig03 <- function() {
  s <- p$sig; st <- p$st; ph <- p$ph
  grid <- make_grid(p)
  thr <- p$tum$rd * grid$P0_ref
  ic <- seq(0, 0.95, length.out = 40); ie <- seq(0, 0.95, length.out = 40)
  gg <- expand.grid(Ic = ic, Ie = ie)
  d <- NULL
  for (row in list(c("letrozole", ph$E2_let, 0), c("giredestrant 30 mg", 1, ph$uG_ref))) for (g in c(0, 1, 4, 5)) {
    r <- steady_state(as.numeric(row[2]), as.numeric(row[3]), gg$Ie, gg$Ic, 0.25, g, s, st)
    d <- rbind(d, data.frame(gg, P = r$P, panel = sprintf("%s clone, %s", GENO_NAMES[g + 1], row[1]),
                             row = row[1], geno = GENO_NAMES[g + 1]))
  }
  d$geno <- factor(d$geno, levels = c("S", "E", "R", "ER"))
  d$row <- factor(d$row, levels = c("letrozole", "giredestrant 30 mg"))
  g <- ggplot(d, aes(Ic, Ie, z = P)) + geom_contour_filled(aes(fill = after_stat(level_mid)), breaks = seq(0, 1.1, 0.05)) +
    geom_contour(breaks = thr, colour = "white", linewidth = 0.8) +
    facet_grid(row ~ geno, labeller = labeller(geno = function(x) paste(x, "clone"))) +
    scale_fill_viridis_c(name = "proliferation\nindex P", limits = c(0, 1.1)) +
    labs(x = "CDK4/6 inhibition", y = "mTORC1 (S6K) inhibition",
         title = "Steady-state proliferation by clone and drug exposure. White line: growth/regression threshold for a median patient (the clone shrinks above and to the right of the line)") +
    theme_q() + coord_fixed()
  save_fig(g, "fig03_proliferation_maps.png", 13, 6.5)
}

# ---------------------------------------------------------------------------- 4
fig04 <- function() {
  ph <- p$ph; tx <- p$tox
  sc <- list(list("everolimus 10 mg alone", NULL, 0, 1, "#7f8c8d"), list("+ palbociclib 125 mg, eve 10 mg", "palbo", 1, 1, "#0072B2"),
             list("+ ribociclib 600 mg, eve 10 mg", "ribo", 1, 1, "#D55E00"), list("+ ribociclib 600 mg, eve 2.5 mg", "ribo", 1, 0.25, "#E69F00"),
             list("+ abemaciclib 150 mg BID, eve 10 mg", "abema", 1, 1, "#009E73"))
  d1 <- NULL; d2 <- NULL
  for (x in sc) {
    pk <- pk_profile(ph, x[[2]], x[[3]], x[[4]], 84)
    d1 <- rbind(d1, data.frame(day = pk$day, y = pk$eve, lab = x[[1]]))
    if (!is.null(x[[2]]) && x[[4]] == 1) d2 <- rbind(d2, data.frame(day = pk$day, y = pk$cyp, lab = strsplit(x[[1]], ",")[[1]][1]))
  }
  cols <- stats::setNames(sapply(sc, `[[`, 5), sapply(sc, `[[`, 1))
  ktr <- 4 / tx$MTT
  d3 <- NULL
  for (x in list(c("palbociclib 125 mg", "palbo", tx$slope_palbo), c("ribociclib 600 mg", "ribo", tx$slope_ribo),
                 c("abemaciclib 150 mg BID", "abema", tx$slope_abema))) {
    cc <- pk_profile(ph, x[2], 1, 0, 84)$cdk
    Pn <- T1 <- T2 <- T3 <- C <- 1; anc <- numeric(84)
    for (dd in 1:84) {
      E <- min(as.numeric(x[3]) * cc[dd], 0.99)
      for (sub in 1:2) {
        fb <- (1 / max(C, 1e-3))^tx$gam_N
        dP <- ktr * Pn * (1 - E) * fb - ktr * Pn
        d1_ <- ktr * (Pn - T1); d2_ <- ktr * (T1 - T2); d3_ <- ktr * (T2 - T3); dC <- ktr * (T3 - C)
        Pn <- Pn + .5 * dP; T1 <- T1 + .5 * d1_; T2 <- T2 + .5 * d2_; T3 <- T3 + .5 * d3_; C <- C + .5 * dC
      }
      anc[dd] <- tx$ANC0_med * C
    }
    d3 <- rbind(d3, data.frame(day = 0:83, y = anc, lab = x[1]))
  }
  g1 <- ggplot(d1, aes(day, y, colour = lab)) + geom_line() + scale_colour_manual(values = cols, name = NULL) +
    labs(title = "Everolimus exposure (relative to 10 mg alone)", x = "day", y = NULL) + theme_q() +
    theme(legend.position = "bottom", legend.direction = "vertical", legend.text = element_text(size = 7))
  g2 <- ggplot(d2, aes(day, y, colour = lab)) + geom_line() +
    scale_colour_manual(values = c("+ palbociclib 125 mg" = "#0072B2", "+ ribociclib 600 mg" = "#D55E00",
                                   "+ abemaciclib 150 mg BID" = "#009E73"), name = NULL) +
    labs(title = "Hepatic/gut CYP3A4 activity", x = "day", y = NULL) + theme_q() +
    theme(legend.position = "bottom", legend.direction = "vertical", legend.text = element_text(size = 7))
  g3 <- ggplot(d3, aes(day, y, colour = lab)) + geom_line() + geom_hline(yintercept = 1, linetype = 2, linewidth = 0.4) +
    annotate("text", x = 4, y = 1.12, label = "grade 3", size = 2.5) +
    scale_colour_manual(values = c("palbociclib 125 mg" = "#0072B2", "ribociclib 600 mg" = "#D55E00",
                                   "abemaciclib 150 mg BID" = "#009E73"), name = NULL) +
    labs(title = "Neutrophils, typical patient (Friberg model)", x = "day", y = "ANC (10^9/L)") + theme_q() +
    theme(legend.position = "bottom", legend.direction = "vertical", legend.text = element_text(size = 7))
  save_fig(g1 + g2 + g3 + plot_annotation(title = "Pharmacokinetics, CYP3A4 drug interaction (ribociclib and palbociclib prescribing information) and myelosuppression"),
           "fig04_pk_ddi_neutrophils.png", 13, 4.8)
}

# ---------------------------------------------------------------------------- 5
fig05 <- function() {
  df <- utils::read.csv(file.path(RES, "calibration_table.csv"), stringsAsFactors = FALSE)
  d <- readRDS(file.path(RES, "calibration_curves.rds"))
  T <- CAL_TARGETS
  mk <- function(curves, pts, ttl, xmax, cols) {
    dd <- do.call(rbind, curves)
    ggplot(dd, aes(time, surv, colour = label)) + geom_step(linewidth = 0.6) +
      geom_point(data = pts, aes(x, y, colour = label), shape = 21, fill = "white", size = 2.6, stroke = 1) +
      scale_colour_manual(values = cols, name = NULL) + coord_cartesian(xlim = c(0, xmax)) +
      labs(title = ttl, x = "months", y = "%") + theme_q() +
      theme(legend.position = "bottom", legend.direction = "vertical")
  }
  g1 <- mk(list(km_df(d$p2c_t, d$p2c_e, NULL, "palbociclib + letrozole"), km_df(d$p2a_t, d$p2a_e, NULL, "letrozole")),
           data.frame(x = c(T$P2_cdk_mPFS, T$P2_ai_mPFS), y = 50, label = c("palbociclib + letrozole", "letrozole")),
           "PALOMA-2: 1L PFS (scans every 12 weeks; extended follow-up)", 60,
           c("palbociclib + letrozole" = "#c0392b", "letrozole" = "#2471a3"))
  g2 <- mk(list(km_df(d$b2e_t, d$b2e_e, NULL, "everolimus + exemestane"), km_df(d$b2x_t, d$b2x_e, NULL, "exemestane")),
           data.frame(x = c(7.8, 3.2), y = 50, label = c("everolimus + exemestane", "exemestane")),
           "BOLERO-2: 2L PFS after an aromatase inhibitor (scans every 6 weeks)", 30,
           c("everolimus + exemestane" = "#d68910", "exemestane" = "#7f8c8d"))
  g3 <- mk(list(km_df(d$evc_t, d$evc_e, d$evc_w, "evERA control 2L PFS (ESR1-enriched)"),
                km_df(d$ml2_t, d$ml2_e, NULL, "MONALEESA-2 OS, ribociclib + letrozole")),
           data.frame(x = c(5.5, 63.9), y = 50, label = c("evERA control 2L PFS (ESR1-enriched)", "MONALEESA-2 OS, ribociclib + letrozole")),
           "evERA control 2L PFS; MONALEESA-2 OS", 96,
           c("evERA control 2L PFS (ESR1-enriched)" = "#7f8c8d", "MONALEESA-2 OS, ribociclib + letrozole" = "#8e44ad"))
  gl <- list(g1, g2, g3)
  if (!is.null(d$pvg_t)) {
    LM <- T$PV_landmarks
    pts <- rbind(data.frame(x = c(T$PV_gir_mPFS, T$PV_let_mPFS), y = 50, label = c("giredestrant + palbociclib", "letrozole + palbociclib")),
                 data.frame(x = rep(LM$m, 2), y = c(LM$gir, LM$let), label = rep(c("giredestrant + palbociclib", "letrozole + palbociclib"), each = 3)))
    gl[[4]] <- mk(list(km_df(d$pvg_t, d$pvg_e, NULL, "giredestrant + palbociclib"), km_df(d$pvl_t, d$pvl_e, NULL, "letrozole + palbociclib")),
                  pts, "persevERA: 1L PFS (medians and 12/24/36-month landmarks)", 60,
                  c("giredestrant + palbociclib" = "#c0392b", "letrozole + palbociclib" = "#2471a3"))
  }
  if (!is.null(d$p3pf_t)) {
    gl[[5]] <- mk(list(km_df(d$p3pf_t, d$p3pf_e, NULL, "palbociclib + fulvestrant"), km_df(d$p3f_t, d$p3f_e, NULL, "fulvestrant")),
                  data.frame(x = c(T$P3_pf_mPFS, T$P3_f_mPFS), y = 50, label = c("palbociclib + fulvestrant", "fulvestrant")),
                  "PALOMA-3: PFS after progression on an aromatase inhibitor", 24,
                  c("palbociclib + fulvestrant" = "#c0392b", "fulvestrant" = "#7f8c8d"))
  }
  if (!is.null(d$pd_net_t)) {
    pub <- data.frame(x = c(12, 24, 36), y = c(12.9, 28.2, 38.6), label = "rising ESR1 mutation (net incidence)")
    k <- km(d$pd_net_t, d$pd_net_e)
    inc <- data.frame(time = c(0, k$t / MONTH), surv = 100 * (1 - c(1, k$s)), label = "rising ESR1 mutation (net incidence)")
    gl[[6]] <- mk(list(inc, km_df(d$pd_pfs_t, d$pd_pfs_e, NULL, "PFS, aromatase inhibitor + palbociclib")),
                  rbind(pub, data.frame(x = T$PD1$mPFS, y = 50, label = "PFS, aromatase inhibitor + palbociclib")),
                  "PADA-1 step 1: ESR1 emergence (points: read from the published curve) and PFS", 48,
                  c("rising ESR1 mutation (net incidence)" = "#17a589", "PFS, aromatase inhibitor + palbociclib" = "#2471a3"))
  }
  df$err <- ifelse(is.na(df$rel_error_pct), df$abs_error_points, df$rel_error_pct)
  df$ok <- ifelse(is.na(df$rel_error_pct), abs(df$abs_error_points) <= 5, abs(df$rel_error_pct) <= 15)
  df$lab <- sprintf("%s: %s (target %g, model %.3g)%s", df$trial, df$target, df$published, df$model,
                    ifelse(is.na(df$rel_error_pct), " [points]", ""))
  df$lab <- gsub("randomiz", "randomis", noai(df$lab), fixed = TRUE)   # British spelling as in the text
  df$lab <- factor(df$lab, levels = rev(df$lab))
  g4 <- ggplot(df, aes(err, lab, fill = ok)) + geom_col() +
    geom_vline(xintercept = 0) + geom_vline(xintercept = c(-15, -5, 5, 15), linetype = 3) +
    scale_fill_manual(values = c(`TRUE` = "#0072B2", `FALSE` = "#E69F00"), guide = "none") +
    labs(x = "error, model vs published: relative (%) for medians and HRs; percentage points for proportions [points]",
         y = NULL) + theme_q(8)
  top <- wrap_plots(gl, nrow = 2)
  save_fig(wrap_elements(full = top) / g4 + plot_layout(heights = c(1.5, 1.6)) +
             plot_annotation(title = "Calibration: simulated curves (lines) vs published values (open circles); all targets on an independent virtual population (n = 3,000)"),
           "fig05_calibration.png", 14, 15)
}

# ---------------------------------------------------------------------------- 6
fig06 <- function() {
  d <- readRDS(file.path(RES, "validation_curves.rds"))
  es <- d$esr1
  obs <- list("ESR1-mutated" = list(lm = rbind(c(6, 66.1, 38.1), c(12, 40.5, 15.2)), med = c(10.0, 5.5)),
              "Overall (ESR1-enriched to 55.5%)" = list(lm = rbind(c(6, 57.4, 39.6), c(12, 34.1, 18.1)), med = c(8.8, 5.5)),
              "No detectable ESR1 mutation" = list(lm = NULL, med = c(5.7, 5.5)))
  gl <- list()
  for (lab in names(obs)) {
    sel <- if (startsWith(lab, "ESR1")) es else if (startsWith(lab, "Overall")) rep(TRUE, length(es)) else !es
    w <- if (startsWith(lab, "Overall")) d$w[sel] else NULL
    dd <- rbind(km_df(d$te[sel], d$ee[sel], w, "model: giredestrant + everolimus"), km_df(d$tc[sel], d$ec[sel], w, "model: ET + everolimus"))
    pts <- data.frame(x = obs[[lab]]$med, y = 50, label = c("model: giredestrant + everolimus", "model: ET + everolimus"), shp = "median")
    if (!is.null(obs[[lab]]$lm)) pts <- rbind(pts, data.frame(x = rep(obs[[lab]]$lm[, 1], 2), y = c(obs[[lab]]$lm[, 2], obs[[lab]]$lm[, 3]),
                                                              label = rep(c("model: giredestrant + everolimus", "model: ET + everolimus"), each = 2), shp = "landmark"))
    gl[[lab]] <- ggplot(dd, aes(time, surv, colour = label)) + geom_step(linewidth = 0.6) +
      geom_point(data = pts, aes(x, y, colour = label, shape = shp), fill = "white", size = 2.4, stroke = 1) +
      scale_shape_manual(values = c(median = 23, landmark = 22), name = "published") +
      scale_colour_manual(values = c("model: giredestrant + everolimus" = "#c0392b", "model: ET + everolimus" = "#7f8c8d"), name = NULL) +
      coord_cartesian(xlim = c(0, 30)) + labs(title = lab, x = "months from start of 2nd line", y = "% progression-free") + theme_q()
  }
  save_fig(wrap_plots(gl, nrow = 1, guides = "collect") & theme(legend.position = "bottom"),
           "fig06_validation_evERA.png", 13, 4.4)
}

# ---------------------------------------------------------------------------- 7
fig07 <- function() {
  df <- utils::read.csv(file.path(RES, "validation.csv"), stringsAsFactors = FALSE)
  df <- df[!is.na(df$passed), ]
  df$o <- suppressWarnings(as.numeric(df$observed)); df$m <- suppressWarnings(as.numeric(df$model))
  df <- df[is.finite(df$o) & is.finite(df$m), ]
  df$rel <- 100 * (df$m - df$o) / df$o
  df$endpoint <- noai(df$endpoint)
  for (pr in list(c("gir\\+eve", "giredestrant + everolimus"), c("ET\\+eve", "ET + everolimus"), c("\\(esr1\\)", "(ESR1-mutated)"),
                  c("\\(all\\)", "(overall)"), c("\\(noesr1\\)", "(no ESR1 mutation)"), c("% exp", "%, experimental arm"), c("% ctl", "%, control arm"),
                  c("^stom_any", "stomatitis, any grade"), c("^stom_g3", "stomatitis, grade 3"), c("^pneu_any", "pneumonitis, any grade"),
                  c("^pneu_g3", "pneumonitis, grade 3 or higher"), c("^ae_disc_eve", "everolimus stopped for AEs"), c("\\b1L\\b", "first line"),
                  c(">=", "≥"), c("<=", "≤")))   # readable endpoint names, as in Supplementary Table 4
    df$endpoint <- gsub(pr[1], pr[2], df$endpoint, perl = TRUE)
  df$lab <- factor(sprintf("%s: %s", df$trial, df$endpoint), levels = rev(unique(sprintf("%s: %s", df$trial, df$endpoint))))
  g <- ggplot(df, aes(rel, lab, fill = passed)) + geom_col() + geom_vline(xintercept = 0) +
    scale_fill_manual(values = c(`TRUE` = "#0072B2", `FALSE` = "#D55E00"), labels = c(`TRUE` = "criterion met", `FALSE` = "not met"), name = NULL) +
    labs(x = "relative difference, model vs published (%)", y = NULL, title = "Validation summary (held-out data only)") +
    theme_q(8) + theme(legend.position = "bottom")
  save_fig(g, "fig07_validation_summary.png", 9, max(4, 0.24 * nrow(df)))
}

# ---------------------------------------------------------------------------- 8
fig08 <- function() {
  curves <- if (file.exists(file.path(RES, "ensemble_curves.rds"))) readRDS(file.path(RES, "ensemble_curves.rds")) else NULL
  gl <- list()
  for (ep in list(c("pfs1", "PFS on line 1"), c("edct", "Endocrine-based disease control time (EDCT)"), c("os", "Overall survival"))) {
    dd <- do.call(rbind, lapply(c("A", "B", "C"), function(k) { R <- ref_arm(k); km_df(R[[ep[1]]], R[[paste0(ep[1], "_ev")]], NULL, LAB[[k]]) }))
    g <- ggplot(dd, aes(time, surv, colour = label)) + geom_step(linewidth = 0.7)
    if (!is.null(curves)) {
      band <- do.call(rbind, lapply(c("A", "B", "C"), function(k) {
        arr <- curves[[paste0(k, "__", ep[1])]]
        # month 72 lies just beyond the 6-year horizon (2190 days), so the band stops at month 71
        data.frame(time = 0:71, lo = 100 * apply(arr[, 1:72], 2, stats::quantile, 0.05),
                   hi = 100 * apply(arr[, 1:72], 2, stats::quantile, 0.95), label = LAB[[k]])
      }))
      g <- g + geom_ribbon(data = band, aes(x = time, ymin = lo, ymax = hi, fill = label), inherit.aes = FALSE, alpha = 0.12)
    }
    gl[[ep[1]]] <- g + scale_colour_manual(values = stats::setNames(COL[c("A", "B", "C")], LAB[c("A", "B", "C")]), name = NULL) +
      scale_fill_manual(values = stats::setNames(COL[c("A", "B", "C")], LAB[c("A", "B", "C")]), guide = "none") +
      scale_x_continuous(breaks = seq(0, 72, 12), limits = c(0, 72)) +
      labs(title = ep[2], x = "months from start of first-line therapy", y = "%") + theme_q()
  }
  save_fig(wrap_plots(gl, nrow = 1, guides = "collect") & theme(legend.position = "bottom"),
           "fig08_scenario1_km.png", 13, 4.6)
}

ens_pivot <- function(metric) {
  es <- utils::read.csv(file.path(RES, "ensemble_summary.csv"), check.names = FALSE, stringsAsFactors = FALSE)
  w <- stats::reshape(es[, c("member", "arm", metric)], idvar = "member", timevar = "arm", direction = "wide")
  names(w) <- sub(paste0("^", metric, "\\."), "", names(w))
  w
}

# ---------------------------------------------------------------------------- 9
fig09 <- function() {
  cons <- list(c("A", "B"), c("A", "B_trial"), c("A", "C"), c("B", "C"), c("D1", "B"), c("D2", "B"), c("D3", "B"))
  gl <- list()
  for (m in list(c("EDCT_RMST48", "Delta RMST of EDCT at 48 months"), c("OS_RMST60", "Delta RMST of OS at 60 months"),
                 c("QTWiST", "Delta Q-TWiST (60 months, u = 0.5/0.5)"))) {
    pv <- ens_pivot(m[1])
    dd <- do.call(rbind, lapply(cons, function(cc) data.frame(contrast = paste0(cc[1], "-", cc[2]), d = pv[[cc[1]]] - pv[[cc[2]]])))
    dd$contrast <- factor(dd$contrast, levels = sapply(cons, paste, collapse = "-"))
    g <- ggplot(dd, aes(contrast, d)) + geom_violin(fill = "#aed6f1", colour = "#2471a3", scale = "width") +
      stat_summary(fun = stats::median, geom = "point", size = 1.6) + geom_hline(yintercept = 0) +
      labs(title = m[2], x = NULL, y = "months") + theme_q()
    if (m[1] == "EDCT_RMST48") g <- g + geom_hline(yintercept = 3, colour = "#c0392b", linetype = 2) +
      annotate("text", x = 1, y = 3.3, label = "decision threshold (3 mo)", colour = "#c0392b", size = 2.5, hjust = 0)
    gl[[m[1]]] <- g
  }
  save_fig(wrap_plots(gl, nrow = 1) + plot_annotation(title = "Paired treatment-strategy differences across the parameter ensemble"),
           "fig09_ensemble_differences.png", 13, 4.3)
}

# ---------------------------------------------------------------------------- 10
fig10 <- function() {
  tr <- readRDS(file.path(RES, "ref_traj.rds"))
  groups <- list("S/P (sensitive)" = c(1, 3), "ESR1-mutant (E, EP)" = c(2, 4), "CDK4/6i escape (R, PR)" = c(5, 7),
                 "ESR1-mutant + escape (ER, EPR)" = c(6, 8))
  gcol <- c("#aed6f1", "#f5b041", "#a569bd", "#c0392b")
  gl <- list()
  for (k in c("A", "B", "C")) {
    N <- tr[[k]]$traj$N     # [clone, patient, week]
    line <- tr[[k]]$traj$line
    tot <- apply(N, c(2, 3), sum)
    ok <- line <= 2 & tot > 0
    W <- dim(N)[3]
    dd <- do.call(rbind, lapply(names(groups), function(gn) {
      x <- apply(N[groups[[gn]], , , drop = FALSE], c(2, 3), sum) / tot
      x[!ok] <- NA
      data.frame(month = (0:(W - 1)) * 7 / MONTH, frac = colMeans(x, na.rm = TRUE), group = gn)
    }))
    dd$frac[!is.finite(dd$frac)] <- 0
    dd$group <- factor(dd$group, levels = rev(names(groups)))
    gl[[k]] <- ggplot(dd, aes(month, frac, fill = group)) + geom_area() + coord_cartesian(xlim = c(0, 60)) +
      scale_fill_manual(values = stats::setNames(rev(gcol), rev(names(groups))), name = NULL) +
      labs(title = paste0(LAB[[k]], ": mean clonal makeup"), x = "months", y = "fraction of tumour cells") + theme_q()
  }
  arms6 <- c("A", "B", "C", "D1", "D2", "D3")
  det <- do.call(rbind, lapply(arms6, function(k) {
    d_ <- ref_arm(k)$esr1_detect
    data.frame(month = 0:60, pct = sapply(0:60, function(m) 100 * mean(d_ >= 0 & d_ <= m * MONTH)), arm = LAB[[k]])
  }))
  g4 <- ggplot(det, aes(month, pct, colour = arm)) + geom_line() + scale_colour_manual(values = stats::setNames(COL, LAB), name = NULL) +
    labs(title = "Cumulative ctDNA ESR1 detection on endocrine therapy", x = "months", y = "% of patients") + theme_q() +
    theme(legend.position = "bottom", legend.text = element_text(size = 6.5)) + guides(colour = guide_legend(ncol = 2))
  cl <- do.call(rbind, lapply(arms6, function(k) {
    R <- ref_arm(k); m <- R$edct_ev & (R$os > R$edct)
    ce <- R$clones_edct[m, , drop = FALSE]
    data.frame(arm = k, type = c("ESR1-mutant cells", "CDK4/6i-escape cells"),
               pct = 100 * c(mean(rowSums(ce[, c(2, 4, 6, 8)])), mean(rowSums(ce[, 5:8]))))
  }))
  g5 <- ggplot(cl, aes(arm, pct, fill = type)) + geom_col(position = "dodge") +
    scale_fill_manual(values = c("ESR1-mutant cells" = "#f5b041", "CDK4/6i-escape cells" = "#a569bd"), name = NULL) +
    labs(title = "Clonal makeup when endocrine-based strategy fails", x = NULL, y = "% of tumour cells") + theme_q() +
    theme(legend.position = "bottom")
  # example patient: one who receives the second endocrine line under B for at least 4 months
  cand <- which(tr$B$edct - tr$B$pfs1 > 120 & tr$A$edct > 0)
  pick <- if (length(cand)) cand[1] else 4
  tw <- do.call(rbind, lapply(c("A", "B", "C"), function(k) {
    N <- tr[[k]]$traj$N; tot <- colSums(N[, pick, ])
    d_ <- data.frame(month = (seq_along(tot) - 1) * 7 / MONTH, sld = 60 * (tot / 2e10)^(1 / 3), arm = LAB[[k]])
    d_[d_$month <= tr[[k]]$edct[pick] / MONTH, ]   # tumour size is tracked until endocrine-based therapy ends
  }))
  ends <- data.frame(arm = LAB[c("A", "B", "C")], e = sapply(c("A", "B", "C"), function(k) tr[[k]]$edct[pick] / MONTH))
  g6 <- ggplot(tw, aes(month, sld, colour = arm)) + geom_line() + geom_vline(data = ends, aes(xintercept = e, colour = arm), linetype = 3) +
    scale_colour_manual(values = stats::setNames(COL[c("A", "B", "C")], LAB[c("A", "B", "C")]), name = NULL) +
    coord_cartesian(xlim = c(0, 72)) + labs(title = "One virtual patient under A, B and C (dotted: end of EDCT)", x = "months",
                                            y = "tumour size (SLD-equivalent, mm)") + theme_q() +
    theme(legend.position = "bottom", legend.direction = "vertical")
  top <- (gl$A | gl$B | gl$C) + plot_layout(guides = "collect")
  save_fig(wrap_elements(full = top) / (g4 | g5 | g6) + plot_annotation(title = "Clonal evolution, ESR1 emergence and an example patient"),
           "fig10_clonal_dynamics.png", 14, 7.5)
}

# ---------------------------------------------------------------------------- 11
fig11 <- function() {
  s <- utils::read.csv(file.path(RES, "ref_summary.csv"), check.names = FALSE, stringsAsFactors = FALSE)
  rownames(s) <- s$arm
  arms <- intersect(c("B", "C", "A", "A_eve5", "A_eve2.5", "A_ribo_eve2.5", "A_ribo300_eve2.5", "A_abema_eve10", "A_abema_eve5",
                      "A_tox1.5", "A_tox2.0", "A_nodex", "A_dexstrict", "A_eve5_dexstrict"), s$arm)
  mets <- c(neut_g3 = "grade ≥3 neutropenia", stom_g3 = "grade 3 stomatitis", pneu_any = "pneumonitis (any)",
            pneu_g3 = "pneumonitis grade ≥3", diar_g3 = "grade ≥3 diarrhoea", hyperg_g3 = "grade ≥3 hyperglycaemia",
            ae_disc_any = "AE-related discontinuation (any drug)")
  dd <- do.call(rbind, lapply(names(mets), function(m) data.frame(arm = arms, met = mets[[m]], pct = 100 * s[arms, m])))
  dd$arm <- factor(dd$arm, levels = arms); dd$met <- factor(dd$met, levels = mets)
  g1 <- ggplot(dd, aes(arm, pct, fill = met)) + geom_col(position = "dodge") +
    scale_fill_manual(values = c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"), name = NULL) +   # colour-blind-safe (Okabe-Ito)
    labs(title = "Adverse events by strategy and triplet variant", x = NULL, y = "% of patients (all endocrine lines)") +
    theme_q() + theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "bottom")
  g2 <- ggplot(data.frame(arm = factor(arms, levels = arms), r = 100 * s[arms, "RDI_eve"]), aes(arm, r)) + geom_col(fill = "#d68910") +
    labs(title = "Relative dose intensity of everolimus (%)", x = NULL, y = NULL) + coord_cartesian(ylim = c(0, 105)) +
    theme_q() + theme(axis.text.x = element_text(angle = 35, hjust = 1))
  save_fig((g1 | g2) + plot_layout(widths = c(3, 1.3)) + plot_annotation(title = "Safety and dose intensity (reference population)"),
           "fig11_toxicity.png", 14, 5.4)
}

# ---------------------------------------------------------------------------- 12
fig12 <- function() {
  A <- ref_arm("A"); B <- ref_arm("B"); tau <- 60 * MONTH
  us <- seq(0, 1, length.out = 21)
  q <- function(R, ut, ur) { E <- pmin(R$edct, tau); O <- pmin(R$os, tau); T <- pmin(R$tox_days, E); mean(ut * T + (E - T) + ur * (O - E)) / MONTH }
  dd <- expand.grid(ut = us, ur = us)
  dd$z <- mapply(function(a, b) q(A, a, b) - q(B, a, b), dd$ut, dd$ur)
  lim <- max(abs(dd$z))
  g <- ggplot(dd, aes(ut, ur, z = z)) + geom_raster(aes(fill = z), interpolate = TRUE) + geom_contour(breaks = 0, colour = "black") +
    scale_fill_distiller(palette = "RdBu", limits = c(-lim, lim), name = "Delta Q-TWiST\nA minus B (mo)") +
    labs(x = "utility of time with toxicity (u_TOX)", y = "utility of time after strategy failure (u_REL)",
         title = "Q-TWiST difference, A minus B, over the utility weights") + theme_q() + coord_fixed()
  save_fig(g, "fig12_qtwist_threshold.png", 6, 4.6)
}

# ---------------------------------------------------------------------------- 13
fig13 <- function() {
  df <- utils::read.csv(file.path(RES, "ensemble_subgroups.csv"), stringsAsFactors = FALSE)
  d <- df[df$contrast %in% c("A-B", "A-B_trial"), ]
  d$key <- paste(d$contrast, d$subgroup, sep = "||")
  ag <- do.call(rbind, lapply(split(d$dEDCT48, d$key), function(x) c(med = stats::median(x), lo = stats::quantile(x, .05, names = FALSE),
                                                                         hi = stats::quantile(x, .95, names = FALSE))))
  order_ <- c("All", "PI3K pathway altered", "PI3K pathway wild-type", "Baseline ESR1 detectable", "Baseline ESR1 not detectable",
              "Endocrine dependence high (b0 low tercile)", "Endocrine dependence mid", "Endocrine dependence low (b0 high tercile)",
              "Tumor burden low tercile", "Tumor burden high tercile", "Tumor DNA fraction < 10%", "Tumor DNA fraction >= 10%",
              "Proliferation low tercile", "Proliferation high tercile",
              "Age < 65", "Age >= 65", "AE susceptibility high tercile")
  order_ <- intersect(order_, unique(d$subgroup))
  kk <- rownames(ag); cc <- sub("||", "\r", kk, fixed = TRUE); cc <- sub("\r.*", "", cc); sg_ <- sub("^.*?[|][|]", "", kk, perl = TRUE)
  dd <- data.frame(sg = factor(sg_, levels = rev(order_)), contrast = ifelse(cc == "A-B", "real-world attrition", "trial-like attrition (15%)"), ag)
  dd <- dd[!is.na(dd$sg), ]
  g <- ggplot(dd, aes(med, sg, colour = contrast)) + geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.3, position = position_dodge(width = 0.6)) +
    geom_point(size = 2, position = position_dodge(width = 0.6)) + geom_vline(xintercept = 0) + geom_vline(xintercept = 3, colour = "#c0392b", linetype = 2) +
    scale_colour_manual(values = c(`real-world attrition` = "#c0392b", `trial-like attrition (15%)` = "#2471a3"), name = "Arm B attrition") +
    labs(x = "Delta RMST of EDCT at 48 months, A minus B (months); median and 90% ensemble interval", y = NULL,
         title = "Subgroup analysis: upfront triplet vs sequential") + theme_q() + theme(legend.position = "bottom")
  save_fig(g, "fig13_subgroups.png", 8, 7.5)
}

# ---------------------------------------------------------------------------- 14
fig14 <- function() {
  gl <- list()
  for (m in list(c("EDCT_RMST48", "EDCT, RMST 48 months"), c("OS_RMST60", "OS, RMST 60 months"))) {
    pv <- ens_pivot(m[1])
    comp <- list("SERD upfront (D1-B)" = pv$D1 - pv$B, "mTOR inhibitor upfront (D2-B)" = pv$D2 - pv$B, "both (A-B)" = pv$A - pv$B,
                 "interaction (A-D1-D2+B)" = pv$A - pv$D1 - pv$D2 + pv$B, "ctDNA-guided switch (D3-B)" = pv$D3 - pv$B)
    dd <- data.frame(lab = factor(names(comp), levels = names(comp)), med = sapply(comp, stats::median),
                     lo = sapply(comp, stats::quantile, 0.05), hi = sapply(comp, stats::quantile, 0.95))
    gl[[m[1]]] <- ggplot(dd, aes(lab, med, fill = lab)) + geom_col() + geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.25) +
      geom_hline(yintercept = 0) + scale_fill_manual(values = c("#8e44ad", "#d68910", "#c0392b", "#95a5a6", "#17a589"), guide = "none") +
      labs(title = m[2], x = NULL, y = "months vs B") + theme_q() + theme(axis.text.x = element_text(angle = 20, hjust = 1))
  }
  save_fig(wrap_plots(gl, nrow = 1) + plot_annotation(title = "Mechanistic decomposition (2 x 2 factorial on upfront SERD and upfront mTOR inhibition; median and 90% ensemble interval)"),
           "fig14_factorial.png", 12, 4.2)
}

# ---------------------------------------------------------------------------- 15
fig15 <- function() {
  s <- utils::read.csv(file.path(RES, "ref_summary.csv"), check.names = FALSE, stringsAsFactors = FALSE)
  s <- s[grepl("^[ABC]", s$arm) & !grepl("^S2", s$arm), ]
  s$col <- substr(s$arm, 1, 1)
  gl <- list()
  for (x in list(c("ae_disc_any", "AE-related discontinuation of any drug (%)"), c("sympt_g3", "selected grade ≥3 AEs (%)"))) {
    gl[[x[1]]] <- ggplot(s, aes(100 * .data[[x[1]]], EDCT_RMST48, colour = col, label = arm)) + geom_point(size = 2) +
      ggrepel::geom_text_repel(size = 2.3, show.legend = FALSE, max.overlaps = Inf, min.segment.length = 0.3, segment.size = 0.2, seed = 1) +
      scale_colour_manual(values = COL[c("A", "B", "C")], guide = "none") +
      labs(x = x[2], y = "RMST of EDCT at 48 months (months)") + theme_q()
  }
  save_fig(wrap_plots(gl, nrow = 1) + plot_annotation(title = "Efficacy–toxicity trade-off across CDK4/6i partner, everolimus dose, toxicity-interaction and prophylaxis variants"),
           "fig15_variants_tradeoff.png", 13, 5)
}

# ---------------------------------------------------------------------------- 16
fig16 <- function() {
  keys <- c("S2_A'(abema)", "S2_A'(palbo cont.)", "S2_B'", "S2_C'")
  S2L <- c("S2_A'(abema)" = "abemaciclib + giredestrant + everolimus", "S2_A'(palbo cont.)" = "palbociclib continued + giredestrant + everolimus",
           "S2_B'" = "giredestrant + everolimus", "S2_C'" = "exemestane + everolimus")
  cols <- c("#c0392b", "#e59866", "#2471a3", "#7f8c8d")
  dd <- do.call(rbind, lapply(keys, function(k) {
    R <- ref_arm(k); m <- R$in_line2
    km_df((R$pfs2 - R$line2_start)[m], R$pfs2_ev[m], NULL, S2L[[k]])
  }))
  g <- ggplot(dd, aes(time, surv, colour = label)) + geom_step(linewidth = 0.7) +
    scale_colour_manual(values = stats::setNames(cols, unname(S2L[keys])), name = NULL) + coord_cartesian(xlim = c(0, 36)) +
    labs(x = "months from start of the second line", y = "% progression-free",
         title = "evERA population: second-line PFS with or without an added CDK4/6 inhibitor") + theme_q() +
    theme(legend.position = c(0.66, 0.78))
  save_fig(g, "fig16_scenario2.png", 7.5, 4.5)
}

# ---------------------------------------------------------------------------- 17
fig17 <- function() {
  PNM <- c(ph.uG_ref = "Giredestrant potency index", sig.eps_g = "Dependence of growth on S6K activity",
           tum.p_attr = "Attrition, progression after more than 12 months", ph.I10_e = "mTORC1 inhibition at everolimus 10 mg",
           tum.log10_fE0 = "Initial ESR1-mutant fraction, median (log10)", tum.p_attr_mid = "Attrition, progression after 6 to 12 months",
           tum.lam_NL = "New-lesion hazard at nadir (per day)", tum.log10_fR0 = "Initial CDK4/6i-escape fraction (log10)",
           tum.early_prog_days = "Early-progression cut-off on line 1 (days)", tum.kappa = "mTORC1 dependence of adaptive escape (κ)",
           tum.h_post = "Death hazard after endocrine therapy (per day)", tum.dc0_med = "Post-endocrine therapy effect (per day)",
           sig.s_mut = "SERD potency shift against mutant ER", tum.sigma_R = "Adaptive escape rate (per day)",
           tum.h_t = "Tumour-burden death hazard coefficient", ph.gir_ddi = "Giredestrant effect on everolimus exposure (fold)",
           tox.disc_eve_gamma = "Exposure dependence of everolimus discontinuation", tum.b0_mean = "Endocrine independence b₀, mean",
           tox.disc_tau = "Time constant of the discontinuation hazard (days)", ph.IC50_e = "Everolimus dose–response IC₅₀",
           tox.disc_floor = "Late discontinuation hazard (share of initial)", tox.tox_int = "Three-drug toxicity interaction factor")
  pnm <- function(x) unname(ifelse(x %in% names(PNM), PNM[x], x))
  pr <- utils::read.csv(file.path(RES, "gsa_prcc.csv"), stringsAsFactors = FALSE)
  ev <- utils::read.csv(file.path(RES, "gsa_evppi.csv"), stringsAsFactors = FALSE)
  gl <- list()
  for (o in list(c("dEDCT48_AB", "PRCC: EDCT gain, A vs B"), c("dOS60_AB", "PRCC: OS gain, A vs B"))) {
    d <- pr[pr$output == o[1], ]; d$parameter <- pnm(d$parameter); d$parameter <- factor(d$parameter, levels = d$parameter[order(d$PRCC)])
    gl[[o[1]]] <- ggplot(d, aes(PRCC, parameter, fill = p_value < 0.05)) + geom_col() + geom_vline(xintercept = 0) +
      scale_fill_manual(values = c(`TRUE` = "#c0392b", `FALSE` = "#bdc3c7"), guide = "none") + coord_cartesian(xlim = c(-1, 1)) +
      labs(title = o[2], y = NULL) + theme_q(8)
  }
  # value of information for the decision whether the gain of A over B reaches 3 months
  gs <- utils::read.csv(file.path(RES, "gsa_samples.csv"), check.names = FALSE)
  pn <- grep(".", names(gs), fixed = TRUE, value = TRUE)
  nb3 <- gs$dEDCT48_AB - 3
  evb <- function(x, nb, nbins = 8) { o <- order(x); bins <- split(o, cut(seq_along(o), nbins, labels = FALSE))
    cond <- numeric(length(nb)); for (b in bins) cond[b] <- mean(nb[b]); mean(pmax(cond, 0)) - max(mean(nb), 0) }
  d <- data.frame(parameter = pn, value = sapply(pn, function(pp) evb(gs[[pp]], nb3)))
  evpi <- mean(pmax(nb3, 0)) - max(mean(nb3), 0)
  d$parameter <- pnm(d$parameter); d$parameter <- factor(d$parameter, levels = d$parameter[order(d$value)])
  g3 <- ggplot(d, aes(value, parameter)) + geom_col(fill = "#2471a3") +
    labs(title = sprintf("EVPPI: does the gain of A over B reach 3 months?\nEVPI (all parameters) = %.2f months", evpi),
         x = "months", y = NULL) + theme_q(8)
  save_fig((gl[[1]] | gl[[2]] | g3) + plot_annotation(title = "Global sensitivity (red = p < 0.05) and value of information"),
           "fig17_sensitivity_voi.png", 18, 5.6)
}

# ---------------------------------------------------------------------------- 18
fig18 <- function() {
  df <- utils::read.csv(file.path(RES, "power.csv"), stringsAsFactors = FALSE)
  gains <- if ("gain_EDCT48" %in% names(df)) tapply(df$gain_EDCT48, df$contrast, function(x) x[1]) else NULL
  if (!is.null(gains)) df$contrast <- sprintf("%s
(EDCT gain %.1f mo)", df$contrast, gains[df$contrast])
  df$contrast <- factor(df$contrast, levels = unique(df$contrast))
  g <- ggplot(df, aes(n_per_arm, 100 * power, colour = endpoint)) + geom_line() + geom_point() + facet_wrap(~contrast, ncol = 4) +
    geom_hline(yintercept = 80, linetype = 3) + scale_colour_manual(values = c(EDCT_logrank = "#c0392b", EDCT_RMST36 = "#d68910", OS_logrank = "#2471a3"), name = NULL) +
    coord_cartesian(ylim = c(0, 100)) + scale_x_log10(breaks = c(100, 200, 400, 800, 1600)) + labs(x = "patients per arm (log scale)", y = "power (%)",
                                             title = "Simulated power of a two-arm randomized trial (1000 virtual trials per point; 24-month accrual, 36-month minimum follow-up, two-sided alpha 0.05)") +
    theme_q() + theme(legend.position = "bottom")
  save_fig(g, "fig18_power.png", 13, 7.2)
}

# ---------------------------------------------------------------------------- 19
fig19 <- function() {
  fn <- file.path(RES, "structural.csv"); if (!file.exists(fn)) return(invisible())
  df <- utils::read.csv(fn, stringsAsFactors = FALSE)
  SL <- c("kappa = 0" = "κ = 0", "kappa = 1" = "κ = 1", "burden-only death hazard" = "tumour-burden-only death hazard",
          "constant attrition 10%" = "constant 10% attrition", "no S6K-ER cross-talk" = "no S6K–ER cross-talk")
  df$structure <- unname(ifelse(df$structure %in% names(SL), SL[df$structure], df$structure))
  df$structure <- factor(df$structure, levels = unique(df$structure))
  g <- ggplot(df, aes(structure, dEDCT48, fill = contrast)) + geom_col(position = "dodge") + geom_hline(yintercept = 0) +
    geom_hline(yintercept = 3, linetype = 2, colour = "#c0392b") +
    scale_fill_manual(values = c("A-B" = "#c0392b", "A-C" = "#e59866", "B-C" = "#2471a3"), name = NULL) +
    labs(x = NULL, y = "Delta RMST EDCT 48 mo (months)",
         title = "Structural uncertainty: strategy contrasts under alternative model structures (S6K–ER cross-talk, adaptive and reversible escape recalibrated)") +
    theme_q() + theme(axis.text.x = element_text(angle = 15, hjust = 1))
  save_fig(g, "fig19_structural.png", 9.5, 4.3)
}

# ---------------------------------------------------------------------------- 20
fig20 <- function() {
  A <- ref_arm("A"); B <- ref_arm("B"); f <- readRDS(file.path(RES, "ref_features.rds"))
  d <- (pmin(A$edct, 48 * MONTH) - pmin(B$edct, 48 * MONTH)) / MONTH
  ben <- 100 * mean(d >= 3); harm <- 100 * mean(d <= -3)
  g1 <- ggplot(data.frame(d = d), aes(d)) + geom_histogram(breaks = seq(-48, 48, 2), fill = "#7f8c8d") +
    geom_vline(xintercept = 3, colour = "#c0392b", linetype = 2) + geom_vline(xintercept = -3, colour = "#2471a3", linetype = 2) +
    labs(title = sprintf("Per-patient Delta EDCT (A minus B), 48-month horizon\nbenefit ≥3 mo: %.0f%%, harm ≤−3 mo: %.0f%%, neutral: %.0f%%",
                         ben, harm, 100 - ben - harm), x = "months", y = "virtual patients") + theme_q()
  dec <- cut(f$b0, stats::quantile(f$b0, seq(0, 1, 0.1)), include.lowest = TRUE, labels = FALSE)
  dd <- do.call(rbind, lapply(c(TRUE, FALSE), function(P) data.frame(dec = 1:10, y = sapply(1:10, function(i) mean(d[f$Ptrunc == P & dec == i])),
                                                                    grp = if (P) "PI3K pathway altered" else "PI3K pathway wild-type")))
  g2 <- ggplot(dd, aes(dec, y, colour = grp)) + geom_line() + geom_point() + geom_hline(yintercept = 0) +
    scale_colour_manual(values = c("PI3K pathway altered" = "#d68910", "PI3K pathway wild-type" = "#2471a3"), name = NULL) +
    scale_x_continuous(breaks = 1:10) + labs(x = "decile of ER-independent cyclin D drive (b0; 10 = least endocrine-dependent)",
                                             y = "mean Delta EDCT (months)", title = "Who benefits from the upfront triplet") +
    theme_q() + theme(legend.position = c(0.25, 0.85))
  save_fig((g1 | g2) + plot_annotation(title = "Per-patient differences in the reference population"), "fig20_individual_effects.png", 12, 4.2)
}

# ---------------------------------------------------------------------------- 21
fig21 <- function() {
  gs <- utils::read.csv(file.path(RES, "gsa_samples.csv"), check.names = FALSE)
  g1 <- ggplot(gs, aes(tum.p_attr, dEDCT48_AB, colour = tum.p_attr_mid)) + geom_point(size = 1.6) +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = "black", linewidth = 0.6) +
    geom_hline(yintercept = 3, colour = "#c0392b", linetype = 2) + scale_colour_viridis_c(name = "attrition,\n6-12 months") +
    geom_vline(xintercept = c(0.10, 0.15, 0.29), linetype = 3) +
    annotate("text", x = c(0.10, 0.15, 0.29), y = max(gs$dEDCT48_AB) * c(1, 0.93, 1), label = c("constant 10%", "SONIA (trial)", "ESME-MBC"), size = 2.6, vjust = 1, hjust = -0.05) +
    labs(x = "attrition after > 12 months on line 1 (no endocrine-based 2nd line)", y = "Delta EDCT RMST48, A minus B (months)",
         title = "Gain of the triplet over sequential therapy vs attrition (global sensitivity samples)") + theme_q()
  st <- utils::read.csv(file.path(RES, "structural.csv"), stringsAsFactors = FALSE)
  st <- st[st$structure %in% c("base", "constant attrition 10%"), ]
  st$structure <- ifelse(st$structure == "base", "real-world attrition (57% / 29%)", "constant 10% attrition")
  g2 <- ggplot(st, aes(contrast, dEDCT48, fill = structure)) + geom_col(position = "dodge") + geom_hline(yintercept = 3, colour = "#c0392b", linetype = 2) +
    scale_fill_manual(values = c("real-world attrition (57% / 29%)" = "#2471a3", "constant 10% attrition" = "#aab7b8"), name = NULL) +
    labs(x = NULL, y = "Delta EDCT RMST48 (months)", title = "Real-world vs constant 10% attrition (n = 3,000)") + theme_q() +
    theme(legend.position = "bottom")
  save_fig((g1 | g2) + plot_layout(widths = c(1.6, 1)), "fig21_attrition.png", 12.5, 4.6)
}

# ---------------------------------------------------------------------------- 22
fig22 <- function() {
  ph <- p$ph
  doses <- seq(100, 600, by = 25)
  d <- data.frame(dose = doses, model = sapply(doses, function(x) ddi_aucr(ph, "ribo", x / 600)))
  obs <- data.frame(dose = c(400, 600, 250, 300), lo = c(3.8, 5.2, 2, 4), hi = c(3.8, 5.2, 3, 5),
                    src = c("ribociclib prescribing information: midazolam AUC (measured)", "ribociclib prescribing information: midazolam AUC (predicted)",
                            "TRINITI-1: everolimus trough vs expected", "TRINITI-1: everolimus trough vs expected"))
  g <- ggplot(d, aes(dose, model)) + geom_line(colour = "#c0392b") +
    geom_pointrange(data = obs, aes(x = dose, y = (lo + hi) / 2, ymin = lo, ymax = hi, shape = src), inherit.aes = FALSE) +
    scale_shape_manual(values = c(16, 1, 17), name = NULL) + coord_cartesian(ylim = c(1, 6)) +
    labs(x = "ribociclib dose (mg/day)", y = "CYP3A substrate exposure ratio",
         title = "Ribociclib–everolimus interaction: model (line) vs prescribing information and TRINITI-1 data",
         caption = "TRINITI-1 reports trough concentrations (continuous ribociclib, everolimus 2.5 mg); trough ratios exceed AUC ratios when clearance is inhibited.") + theme_q() + theme(legend.position = "bottom", legend.direction = "vertical")
  save_fig(g, "fig22_ddi_ribociclib.png", 6.5, 5)
}

# ---------------------------------------------------------------------------- 23: discontinuation
fig23 <- function() {
  b <- utils::read.csv(file.path(RES, "disc_check.csv"), stringsAsFactors = FALSE)
  lab <- paste0(b$trial, ": ", sub(" \\(published: both drugs stopped\\)| \\(9% stopped it while continuing ET\\)", "", b$item))
  lab <- factor(lab, levels = rev(lab))
  d <- rbind(data.frame(item = lab, pct = b$published_pct, what = "published"),
             data.frame(item = lab, pct = b$model_pct, what = "model"))
  g1 <- ggplot(d, aes(pct, item, colour = what, shape = what)) + geom_point(size = 2.4) +
    scale_colour_manual(values = c(published = "black", model = "#c0392b"), name = NULL) +
    scale_shape_manual(values = c(published = 15, model = 16), name = NULL) +
    labs(x = "stopped for adverse events (%)", y = NULL,
         title = "AE-related discontinuation in full-model trial emulations\n(treatment until progression)") +
    theme_q() + theme(legend.position = "bottom")
  # cumulative probability of stopping each drug for adverse events while staying on treatment (no progression),
  # average patient (frailty 1, everolimus exposure 1, no interaction factor)
  tx <- p$tox
  days <- 0:(36 * 30)
  s_ <- if (is.null(tx$disc_tau) || !is.finite(tx$disc_tau)) rep(1, length(days)) else
    tx$disc_floor + (1 - tx$disc_floor) * exp(-pmax(days - (if (is.null(tx$disc_t0)) 183 else tx$disc_t0), 0) / tx$disc_tau)
  cum <- function(h) 100 * (1 - exp(-cumsum(h)))
  dd <- rbind(data.frame(day = days, pct = cum(tx$disc_cdk * s_), drug = "palbociclib"),
              data.frame(day = days, pct = cum(tx$disc_eve * s_), drug = "everolimus 10 mg"),
              data.frame(day = days, pct = cum(tx$disc_et * s_), drug = "letrozole (alone or with a CDK4/6i)"))
  g2 <- ggplot(dd, aes(day / MONTH, pct, colour = drug)) + geom_line() +
    scale_colour_manual(values = c(palbociclib = "#2471a3", `everolimus 10 mg` = "#c0392b", `letrozole (alone or with a CDK4/6i)` = "#7f8c8d"), name = NULL) +
    labs(x = "months on treatment", y = "cumulative % stopped for adverse events",
         title = "Cumulative discontinuation on continued treatment\n(other adverse events; average patient, no progression)") +
    theme_q() + theme(legend.position = "bottom", legend.box = "vertical")
  save_fig((g1 | g2) + plot_layout(widths = c(1.25, 1)), "fig23_discontinuation.png", 14, 5.8)
}

# ---------------------------------------------------------------------------- 25: ctDNA-guided switch trials
fig25 <- function() {
  d <- readRDS(file.path(RES, "validation_curves.rds"))
  if (is.null(d$pd2c) || is.null(d$s6c)) { cat("figure 25: no switch-trial curves\n"); return(invisible()) }
  VT <- VAL_TARGETS
  mk <- function(c_, pts, ttl, labs_, cols, xmax) {
    dd <- rbind(km_df(c_$ts, c_$es, NULL, labs_[1]), km_df(c_$tc, c_$ec, NULL, labs_[2]))
    ggplot(dd, aes(time, surv, colour = label)) + geom_step(linewidth = 0.6) +
      geom_point(data = pts, aes(x, y, colour = label, shape = what), fill = "white", size = 2.6, stroke = 1) +
      scale_shape_manual(values = c(median = 23, landmark = 22), name = "published") +
      scale_colour_manual(values = cols, name = NULL) + coord_cartesian(xlim = c(0, xmax)) +
      labs(title = ttl, x = "months from randomisation", y = "% progression-free") + theme_q() +
      theme(legend.position = "bottom", legend.box = "vertical")
  }
  l1 <- c("switch to fulvestrant + palbociclib", "continue aromatase inhibitor + palbociclib")
  g1 <- mk(d$pd2c, data.frame(x = c(VT$PADA1$switch[1], VT$PADA1$cont[1]), y = 50, label = l1, what = "median"),
           sprintf("PADA-1 step 2 (held out; model n = %d)", d$pd2_n), l1, stats::setNames(c("#17a589", "#7f8c8d"), l1), 30)
  l2 <- c("switch to SERD (giredestrant as proxy) + CDK4/6i", "continue aromatase inhibitor + CDK4/6i")
  S6 <- VT$SERENA6x
  g2 <- mk(d$s6c, rbind(data.frame(x = c(S6$switch, S6$cont), y = 50, label = l2, what = "median"),
                        data.frame(x = 24, y = S6$PFS24, label = l2, what = "landmark")),
           sprintf("SERENA-6 extended analysis (held out; model n = %d)", d$s6_n), l2, stats::setNames(c("#c0392b", "#7f8c8d"), l2), 36)
  save_fig((g1 | g2) + plot_annotation(title = "ctDNA-guided switch on a rising ESR1 mutation: model (lines) vs published (symbols)"),
           "fig25_switch_trials.png", 12.5, 5)
}

# ---------------------------------------------------------------------------- 26: persevERA
fig26 <- function() {
  arms <- persevera_arms()
  v <- make_vpop(p, 4000, 91); g <- make_grid(p)
  dd <- do.call(rbind, lapply(names(arms), function(k) { R <- simulate(p, v, g, arms[[k]], 2190L)
    km_df(R$pfs1, R$pfs1_ev, NULL, if (k == "PV_gir") "giredestrant + palbociclib" else "letrozole + palbociclib") }))
  LM <- CAL_TARGETS$PV_landmarks
  pts <- data.frame(x = c(LM$m, CAL_TARGETS$PV_gir_mPFS, LM$m, CAL_TARGETS$PV_let_mPFS), y = c(LM$gir, 50, LM$let, 50),
                    label = rep(c("giredestrant + palbociclib", "letrozole + palbociclib"), each = 4))
  cols <- c("giredestrant + palbociclib" = "#c0392b", "letrozole + palbociclib" = "#2471a3")
  g1 <- ggplot(dd, aes(time, surv, colour = label)) + geom_step(linewidth = 0.6) +
    geom_point(data = pts, aes(x, y, colour = label), inherit.aes = FALSE, shape = 21, fill = "white", size = 2.6, stroke = 1) +
    scale_colour_manual(values = cols, name = NULL) +
    coord_cartesian(xlim = c(0, 60)) + labs(title = "persevERA: first-line PFS, model (lines) vs published (circles)",
                                            x = "months", y = "% progression-free") + theme_q() + theme(legend.position = "bottom")
  # ER transcriptional activity under each endocrine agent relative to letrozole (typical b0, palbociclib-level CDK4/6 inhibition)
  f <- function(E2, u, gg) steady_state(E2, u, 0, 0.5, p$tum$b0_mean, gg, p$sig, p$st, nthreads = 1)$A_ER
  da <- do.call(rbind, lapply(0:1, function(gg) data.frame(clone = if (gg == 0) "wild-type ER" else "ESR1-mutant",
    agent = c("giredestrant 30 mg", "fulvestrant 500 mg"), ratio = c(f(1, p$ph$uG_ref, gg), f(1, p$ph$uF_ref, gg)) / f(p$ph$E2_let, 0, gg))))
  da$clone <- factor(da$clone, levels = c("wild-type ER", "ESR1-mutant"))
  g2 <- ggplot(da, aes(agent, ratio)) + geom_col(width = 0.6, fill = "#c0392b") +
    geom_hline(yintercept = 1, linetype = 2) + facet_wrap(~clone) +
    labs(title = "ER activity relative to letrozole (1 = same as letrozole)", x = NULL, y = "ER activity, ratio to letrozole") +
    theme_q()
  save_fig((g1 | g2) + plot_layout(widths = c(1.2, 1)), "fig26_persevera.png", 13.5, 5.2)
}

FIGS <- list(fig01, fig02, fig03, fig04, fig05, fig06, fig07, fig08, fig09, fig10, fig11, fig12, fig13, fig14, fig15, fig16,
             fig17, fig18, fig19, fig20, fig21, fig22, fig23, fig25, fig26)
names(FIGS) <- c(1:23, 25, 26)
args <- commandArgs(trailingOnly = TRUE)
which_ <- if (length(args)) args else names(FIGS)
for (i in which_) tryCatch(FIGS[[i]](), error = function(e) cat(sprintf("figure %s failed: %s\n", i, conditionMessage(e))))
