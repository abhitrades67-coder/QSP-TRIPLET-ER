# Main figures of the article, drawn from saved results only.
# Arial throughout, colour-blind-safe palette (Okabe-Ito), lines >= 1 pt, 180 mm wide. Each figure is written as a
# vector PDF, a 600-dpi TIFF (LZW) and a 300-dpi PNG (for embedding in the manuscript file).
# Output: QSP_MS_FIG (default: figures/article); PNG copies for embedding: QSP_MS_PNG (default: figures/article)
# Usage: Rscript scripts/11_article_figures.R [figure numbers ...]   (default: all)

if (!exists("ROOT")) ROOT <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), ".."), winslash = "/")
source(file.path(ROOT, "R", "setup.R"))
suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })
OUT <- Sys.getenv("QSP_MS_FIG", file.path(ROOT, "figures", "article"))
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
PNG <- Sys.getenv("QSP_MS_PNG", file.path(ROOT, "figures", "article")); dir.create(PNG, showWarnings = FALSE, recursive = TRUE)

FONT <- "Arial"
BASE <- 7.5
LW <- 0.5                       # ggplot linewidth 0.5 = about 1.07 pt
ARMCOL <- c(A = "#D55E00", B = "#0072B2", C = "#7F7F7F", D1 = "#CC79A7", D2 = "#E69F00", D3 = "#009E73")
ARMLAB <- c(A = "A, upfront triplet", B = "B, sequential", C = "C, standard", D1 = "D1, SERD first",
            D2 = "D2, everolimus first", D3 = "D3, ctDNA-guided switch")
ATTCOL <- c(`real-world attrition` = "#D55E00", `trial-like attrition (15%)` = "#0072B2")

theme_ms <- function(base = BASE) {
  theme_classic(base_size = base, base_family = FONT) +
    theme(axis.line = element_line(linewidth = LW, colour = "black"), axis.ticks = element_line(linewidth = LW, colour = "black"),
          axis.text = element_text(colour = "black", size = base), axis.title = element_text(size = base),
          legend.text = element_text(size = base), legend.title = element_text(size = base),
          legend.key.size = unit(3.2, "mm"), legend.background = element_blank(),
          plot.title = element_text(size = base, face = "bold", hjust = 0),
          strip.background = element_blank(), strip.text = element_text(size = base, face = "bold"),
          plot.tag = element_text(size = base + 2, face = "bold", family = FONT),
          plot.margin = margin(2, 3, 2, 2, "mm"))
}
tagged <- function(p) p + plot_annotation(tag_levels = "a") & theme(plot.tag = element_text(size = BASE + 2, face = "bold", family = FONT),
                                                                    plot.tag.position = c(0, 1))
GE <- "\u2265"; MINUS <- "\u2212"

save_ms <- function(p, name, w_mm, h_mm) {
  grDevices::cairo_pdf(file.path(OUT, paste0(name, ".pdf")), width = w_mm / 25.4, height = h_mm / 25.4, family = FONT)
  print(p); grDevices::dev.off()
  ragg::agg_tiff(file.path(OUT, paste0(name, ".tif")), width = w_mm, height = h_mm, units = "mm", res = 600, compression = "lzw")
  print(p); grDevices::dev.off()
  ragg::agg_png(file.path(PNG, paste0(name, ".png")), width = w_mm, height = h_mm, units = "mm", res = 300)   # for embedding in the manuscript file
  print(p); grDevices::dev.off()
  cat("saved", name, "\n")
}

q3 <- function(x) { x <- x[is.finite(x)]; c(m = stats::median(x), lo = unname(stats::quantile(x, .05)), hi = unname(stats::quantile(x, .95))) }
es <- utils::read.csv(file.path(RES, "ensemble_summary.csv"), check.names = FALSE, stringsAsFactors = FALSE)
piv <- function(metric) { w <- stats::reshape(es[, c("member", "arm", metric)], idvar = "member", timevar = "arm", direction = "wide")
  names(w) <- sub(paste0("^", metric, "\\."), "", names(w)); w[order(w$member), ] }
E48 <- piv("EDCT_RMST48")

# ---------------------------------------------------------------------------- Figure 1: design and model
fig1 <- function() {
  DRUG <- c(cdk = "#2F4B7C", ai = "#BDBDBD", gir = "#009E73", eve = "#E69F00", exe = "#BDBDBD")
  TXT <- c(cdk = "white", ai = "black", gir = "white", eve = "black", exe = "black")
  rows <- list(
    A = list(l1 = c(cdk = "CDK4/6i", gir = "Giredestrant", eve = "Everolimus"), l2 = NULL),
    B = list(l1 = c(cdk = "CDK4/6i", ai = "Letrozole"), l2 = c(gir = "Giredestrant", eve = "Everolimus")),
    C = list(l1 = c(cdk = "CDK4/6i", ai = "Letrozole"), l2 = c(exe = "Exemestane", eve = "Everolimus")),
    D1 = list(l1 = c(cdk = "CDK4/6i", gir = "Giredestrant"), l2 = c(gir = "Giredestrant", eve = "Everolimus")),
    D2 = list(l1 = c(cdk = "CDK4/6i", ai = "Letrozole", eve = "Everolimus"), l2 = c(gir = "Giredestrant", eve = "Everolimus")),
    D3 = list(l1 = c(cdk = "CDK4/6i", ai = "Letrozole", gir = "Giredestrant"), l2 = c(gir = "Giredestrant", eve = "Everolimus")))
  chips <- list(); y0 <- 6
  for (k in names(rows)) {
    for (ln in c("l1", "l2")) {
      d <- rows[[k]][[ln]]; if (is.null(d)) next
      x_start <- if (ln == "l1") 2.55 else 7.05; wd <- if (ln == "l1") 1.38 else 1.38
      for (i in seq_along(d)) chips[[length(chips) + 1]] <- data.frame(arm = k, xmin = x_start + (i - 1) * (wd + 0.06), xmax = x_start + (i - 1) * (wd + 0.06) + wd,
                                                                       y = y0, drug = names(d)[i], lab = unname(d[i]))
    }
    y0 <- y0 - 1
  }
  ch <- do.call(rbind, chips)
  ch$lab[ch$arm == "D3" & ch$drug == "gir" & ch$xmin < 7] <- "Giredestrant*"
  arms <- data.frame(arm = names(rows), y = 6:1, lab = ARMLAB[names(rows)])
  pa <- ggplot() +
    annotate("rect", xmin = 2.45, xmax = 6.85, ymin = 0.45, ymax = 6.55, fill = "#F3F3F3", colour = NA) +
    annotate("rect", xmin = 6.95, xmax = 9.95, ymin = 0.45, ymax = 6.55, fill = "#F3F3F3", colour = NA) +
    annotate("text", x = 4.65, y = 7.05, label = "First line", fontface = "bold", size = BASE / .pt, family = FONT) +
    annotate("text", x = 8.45, y = 7.05, label = "Second endocrine-based line", fontface = "bold", size = BASE / .pt, family = FONT) +
    geom_rect(data = ch, aes(xmin = xmin, xmax = xmax, ymin = y - 0.32, ymax = y + 0.32, fill = drug), colour = NA) +
    geom_text(data = ch, aes(x = (xmin + xmax) / 2, y = y, label = lab, colour = drug), size = (BASE - 0.5) / .pt, family = FONT) +
    geom_text(data = arms, aes(x = 0.05, y = y, label = lab), hjust = 0, size = BASE / .pt, family = FONT) +
    annotate("text", x = 8.45, y = 6, label = "none (triplet continues)", fontface = "italic", size = (BASE - 0.5) / .pt, family = FONT, colour = "grey30") +
    ggtext::geom_richtext(data = data.frame(x = 0.05, y = 0.05, label = paste0("CDK4/6i: palbociclib 125 mg (base case). \\*D3: letrozole is switched to giredestrant ",
                          "when ctDNA shows an <i>ESR1</i> mutation.<br>After endocrine-based therapy fails, all arms receive the same chemotherapy or ",
                          "antibody\u2013drug conjugate.")), aes(x = x, y = y, label = label), hjust = 0, vjust = 1, size = (BASE - 0.5) / .pt,
                          family = FONT, lineheight = 1.2, fill = NA, label.colour = NA, label.padding = grid::unit(0, "pt")) +
    scale_fill_manual(values = DRUG, guide = "none") + scale_colour_manual(values = TXT, guide = "none") +
    coord_cartesian(xlim = c(0, 10), ylim = c(-0.25, 7.3), clip = "off") + theme_void(base_family = FONT) +
    theme(plot.margin = margin(2, 2, 4, 2, "mm"))
  # b: model modules and workflow
  mods <- data.frame(id = 1:5, x = c(1, 3, 5, 7, 9), y = 2.9,
                     title = c("1 Drug exposure", "2 Signalling", "3 Tumour clones", "4 Toxicity", "5 Outcomes"),
                     body = c("8 drugs; CYP3A4\ninhibition by\nribociclib or\npalbociclib raises\neverolimus exposure",
                              "ER (wild-type or\nESR1-mutant),\nPI3K/AKT/mTOR,\ncyclin D-CDK4/6,\ncyclin E-CDK2",
                              "8 clones: ESR1,\nPI3K and CDK4/6i\nescape; RECIST\nresponse and\nprogression; ctDNA",
                              "neutropenia,\nstomatitis,\npneumonitis and\nothers; label-based\nholds and stops",
                              "PFS, EDCT, OS;\nattrition before\nline 2; chemotherapy\nor ADC after\nendocrine therapy"))
  flow <- data.frame(x = c(0.9, 2.55, 4.2, 5.85, 7.5, 9.15), y = 0.55,
                     lab = c("Inputs\nchecked against\nsources", "Calibration\n14 + 10 + 36\ntargets", "Validation\nheld-out results\nof 10 trials",
                             "Uncertainty\n60 parameter sets,\n7 structures", "In silico trial\n6 strategies \u00d7\n6,000 patients", "Decision\npre-specified\nrule"))
  pb <- ggplot() +
    geom_rect(data = mods, aes(xmin = x - 0.88, xmax = x + 0.88, ymin = y - 1.05, ymax = y + 0.95), fill = "white", colour = "#2F4B7C", linewidth = LW) +
    geom_text(data = mods, aes(x = x, y = y + 0.72, label = title), fontface = "bold", size = BASE / .pt, family = FONT) +
    ggtext::geom_richtext(data = transform(mods, body = gsub("ESR1", "<i>ESR1</i>", gsub("\n", "<br>", body, fixed = TRUE), fixed = TRUE)),
                          aes(x = x, y = y - 0.12, label = body), size = (BASE - 1) / .pt, family = FONT, lineheight = 1.25,
                          fill = NA, label.colour = NA, label.padding = grid::unit(0, "pt")) +
    annotate("segment", x = c(1.9, 3.9, 5.9, 7.9), xend = c(2.1, 4.1, 6.1, 8.1), y = 2.9, yend = 2.9,
             arrow = arrow(length = unit(1.4, "mm"), type = "closed"), linewidth = LW) +
    geom_rect(data = flow, aes(xmin = x - 0.76, xmax = x + 0.76, ymin = y - 0.5, ymax = y + 0.5), fill = "#2F4B7C", colour = NA) +
    geom_text(data = flow, aes(x = x, y = y, label = lab), colour = "white", size = (BASE - 1) / .pt, family = FONT, lineheight = 0.9) +
    annotate("segment", x = flow$x[-6] + 0.77, xend = flow$x[-1] - 0.77, y = 0.55, yend = 0.55,
             arrow = arrow(length = unit(1.2, "mm"), type = "closed"), linewidth = LW) +
    coord_cartesian(xlim = c(0, 10), ylim = c(0, 3.95), clip = "off") + theme_void(base_family = FONT)
  p <- tagged((pa / pb) + plot_layout(heights = c(1.35, 1)))
  save_ms(p, "Figure_1", 180, 150)
}

# ---------------------------------------------------------------------------- Figure 2: calibration and validation
fig2 <- function() {
  cal <- utils::read.csv(file.path(RES, "calibration_table.csv"), stringsAsFactors = FALSE)
  cal$dev <- ifelse(is.na(cal$rel_error_pct), cal$abs_error_points / 5, cal$rel_error_pct / 15)
  lv <- sub("evERA control", "evERA (control arm)", c("PALOMA-2", "BOLERO-2", "MONALEESA-2", "evERA control", "acelERA", "persevERA", "PALOMA-3", "PADA-1"))
  cal$trial <- factor(sub("evERA control", "evERA (control arm)", cal$trial), levels = rev(lv))
  cal$met <- factor(ifelse(abs(cal$dev) <= 1, "within tolerance", "outside tolerance"), levels = c("within tolerance", "outside tolerance"))
  pa <- ggplot(cal, aes(dev, trial)) +
    annotate("rect", xmin = -1, xmax = 1, ymin = -Inf, ymax = Inf, fill = "#E8F1FA") +
    geom_vline(xintercept = 0, linewidth = LW, colour = "grey55") +
    geom_point(aes(shape = met), position = position_jitter(width = 0, height = 0.17, seed = 7), size = 1.5, colour = "#0072B2", stroke = LW) +
    scale_shape_manual(values = c(`within tolerance` = 16, `outside tolerance` = 1), name = NULL) +
    scale_x_continuous(breaks = -2:2) +
    labs(x = "Model minus published\n(in units of the tolerance)", y = NULL) +
    theme_ms() + theme(legend.position = "top", legend.direction = "vertical")
  # held-out hazard ratios
  val <- utils::read.csv(file.path(RES, "validation.csv"), stringsAsFactors = FALSE)
  hr <- val[grepl("HR", val$endpoint) & !is.na(val$obs_CI) & nzchar(val$obs_CI) & grepl("CI", val$criterion), ]
  ci <- do.call(rbind, lapply(strsplit(hr$obs_CI, "-"), as.numeric))
  lab <- c("evERA, <i>ESR1</i>-mutated", "evERA, overall", "evERA, no <i>ESR1</i> mutation", "SERENA-6, PFS", "SERENA-6, OS", "PADA-1, PFS",
           "postMONARCH, PFS", "MAINTAIN, PFS", "MIRACLE, PFS", "PARSIFAL, PFS", "FALCON, PFS")
  key <- c("evERA|HR PFS (esr1)", "evERA|HR PFS (all)", "evERA|HR PFS (noesr1)", "SERENA-6 (SERD proxy)|HR", "SERENA-6 (SERD proxy)|OS HR (immature)",
           "PADA-1 step 2|HR", "postMONARCH|HR abemaciclib + fulvestrant vs fulvestrant", "MAINTAIN|HR ribociclib + fulvestrant vs fulvestrant",
           "MIRACLE|HR everolimus + letrozole vs letrozole (1L)", "PARSIFAL|HR fulvestrant + palbociclib vs letrozole + palbociclib (1L)",
           "FALCON|HR fulvestrant vs AI (1L, endocrine-naive)")
  hk <- paste(hr$trial, hr$endpoint, sep = "|")
  stopifnot(all(key %in% hk), length(key) == nrow(hr))
  i <- match(key, hk)
  dv <- data.frame(lab = factor(lab, levels = rev(lab)), pub = as.numeric(hr$observed[i]), lo = ci[i, 1], hi = ci[i, 2], mod = as.numeric(hr$model[i]))
  forest <- function(d, xlab) ggplot(d, aes(y = lab)) +
    geom_vline(xintercept = 1, linetype = 2, linewidth = LW, colour = "grey55") +
    geom_errorbar(aes(xmin = lo, xmax = hi), width = 0.3, linewidth = LW, orientation = "y") +
    geom_point(aes(x = pub, shape = "published (95% CI)"), size = 1.6) +
    geom_point(aes(x = mod, shape = "model"), size = 2, colour = "#D55E00") +
    scale_shape_manual(values = c(`published (95% CI)` = 15, model = 18), name = NULL) +
    scale_x_log10(breaks = c(0.25, 0.5, 1, 2), labels = function(x) format(x, drop0trailing = TRUE)) + labs(x = xlab, y = NULL) +
    theme_ms() + theme(legend.position = "top", axis.text.y = ggtext::element_markdown())
  pb <- forest(dv, "Hazard ratio (log scale)")
  # external validation (SONIA; persevERA OS)
  pr <- utils::read.csv(file.path(RES, "external_validation.csv"), stringsAsFactors = FALSE, check.names = FALSE)
  ek <- c("^PFS1 HR", "^PFS2 HR \\(model", "^OS HR$", "^PFS2 HR, high tumor DNA \\(TF >= 10%\\)$", "^PFS2 HR, low tumor DNA \\(TF < 10%\\)$")
  el <- c("SONIA, PFS1", "SONIA, PFS2 (primary)", "SONIA, OS", "SONIA, PFS2, high tumour DNA", "SONIA, PFS2, low tumour DNA")
  so <- pr[pr$trial %in% c("SONIA", "SONIA ctDNA"), ]
  rr <- lapply(seq_along(ek), function(j) { x <- so[grepl(ek[j], so$item), ][1, ]; c(x$published, x$published_ci, x$model) })
  pv <- pr[pr$trial == "persevERA" & pr$item == "OS HR", ][1, ]
  rr[[length(rr) + 1]] <- c(pv$published, pv$published_ci, pv$model); el <- c(el, "persevERA, OS")
  m <- do.call(rbind, rr); cie <- do.call(rbind, lapply(strsplit(m[, 2], "-"), as.numeric))
  de <- data.frame(lab = factor(el, levels = rev(el)), pub = as.numeric(m[, 1]), lo = cie[, 1], hi = cie[, 2], mod = as.numeric(m[, 3]))
  pc <- forest(de, "Hazard ratio (log scale)") + theme(legend.position = "none")
  p <- tagged((pa | (pb / pc + plot_layout(heights = c(11, 6)))) + plot_layout(widths = c(1, 1.15)))
  save_ms(p, "Figure_2", 180, 115)
}

# ---------------------------------------------------------------------------- Figure 3: efficacy
fig3 <- function() {
  cv <- readRDS(file.path(RES, "ensemble_curves.rds"))
  months <- 0:60
  dd <- do.call(rbind, lapply(c("A", "B", "C"), function(a) { m <- cv[[paste0(a, "__edct")]]
    m <- m[, seq_along(months), drop = FALSE]   # columns are months 0..72; the last one is the end of the 6-year horizon
    data.frame(arm = a, t = months, med = 100 * apply(m, 2, stats::median), lo = 100 * apply(m, 2, stats::quantile, .05),
               hi = 100 * apply(m, 2, stats::quantile, .95)) }))
  pa <- ggplot(dd, aes(t, med, colour = arm, fill = arm)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.18, colour = NA) + geom_line(linewidth = LW + 0.1) +
    scale_colour_manual(values = ARMCOL, labels = ARMLAB, name = NULL) + scale_fill_manual(values = ARMCOL, labels = ARMLAB, name = NULL) +
    scale_x_continuous(breaks = seq(0, 60, 12)) + coord_cartesian(xlim = c(0, 60), ylim = c(0, 100), expand = FALSE) +
    labs(x = "Months from start of first line", y = "In endocrine-based disease control (%)") +
    theme_ms() + theme(legend.position = "inside", legend.position.inside = c(0.67, 0.86))
  g <- data.frame(member = E48$member, rw = E48$A - E48$B, tr = E48$A - E48$B_trial)
  gl <- rbind(data.frame(att = "real-world attrition", gain = g$rw), data.frame(att = "trial-like attrition (15%)", gain = g$tr))
  gl$att <- factor(gl$att, levels = names(ATTCOL))
  sm <- do.call(rbind, lapply(levels(gl$att), function(a) { x <- q3(gl$gain[gl$att == a]); data.frame(att = a, m = x[["m"]], lo = x[["lo"]], hi = x[["hi"]]) }))
  sm$att <- factor(sm$att, levels = names(ATTCOL))
  pb <- ggplot(gl, aes(att, gain, colour = att)) +
    geom_hline(yintercept = 3, linetype = 2, linewidth = LW, colour = "grey40") + geom_hline(yintercept = 0, linewidth = LW, colour = "grey70") +
    geom_point(position = position_jitter(width = 0.15, height = 0, seed = 3), size = 1, alpha = 0.6) +
    geom_errorbar(data = sm, aes(x = att, ymin = lo, ymax = hi), inherit.aes = FALSE, width = 0.12, linewidth = LW, position = position_nudge(x = 0.32)) +
    geom_point(data = sm, aes(x = att, y = m), inherit.aes = FALSE, size = 2, position = position_nudge(x = 0.32)) +
    scale_colour_manual(values = ATTCOL, guide = "none") + scale_x_discrete(labels = c("real-world\nattrition", "trial-like\nattrition (15%)")) +
    labs(x = NULL, y = "EDCT gain, A minus B (months)") + theme_ms()
  xs <- seq(0, 5, 0.05)
  pc_d <- rbind(data.frame(att = names(ATTCOL)[1], x = xs, p = sapply(xs, function(v) mean(g$rw >= v))),
                data.frame(att = names(ATTCOL)[2], x = xs, p = sapply(xs, function(v) mean(g$tr >= v))))
  pc_d$att <- factor(pc_d$att, levels = names(ATTCOL))
  p3 <- data.frame(att = factor(names(ATTCOL), levels = names(ATTCOL)), x = 3, p = c(mean(g$rw >= 3), mean(g$tr >= 3)))
  pc <- ggplot(pc_d, aes(x, p, colour = att)) +
    annotate("rect", xmin = 3, xmax = Inf, ymin = 0.8, ymax = Inf, fill = "#EDEDED") +
    geom_hline(yintercept = 0.8, linetype = 2, linewidth = LW, colour = "grey40") + geom_vline(xintercept = 3, linetype = 2, linewidth = LW, colour = "grey40") +
    geom_step(linewidth = LW + 0.1) + geom_point(data = p3, size = 2) +
    geom_text(data = p3, aes(label = sprintf("%.2f", p)), hjust = -0.3, vjust = 0.5, size = BASE / .pt, family = FONT, show.legend = FALSE) +
    scale_colour_manual(values = ATTCOL, name = NULL) + scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    labs(x = "Efficacy bar x (months)", y = paste0("P(gain ", GE, " x) across parameter sets")) +
    guides(colour = guide_legend(nrow = 2)) + theme_ms() + theme(legend.position = "top")
  p <- tagged((pa | pb | pc) + plot_layout(widths = c(1.35, 0.85, 1.1)))
  save_ms(p, "Figure_3", 180, 75)
}

# ---------------------------------------------------------------------------- Figure 5: toxicity and the decision
fig4 <- function() {
  DI <- piv("ae_disc_any")
  arms <- c("A", "A_let", "A_eve5", "C", "D1", "D2", "D3")
  alab <- c(A = "A", A_let = "A with letrozole", A_eve5 = "A, everolimus 5 mg", C = "C", D1 = "D1", D2 = "D2", D3 = "D3")
  acol <- c(A = ARMCOL[["A"]], A_let = ARMCOL[["A"]], A_eve5 = ARMCOL[["A"]], C = ARMCOL[["C"]], D1 = ARMCOL[["D1"]], D2 = ARMCOL[["D2"]], D3 = ARMCOL[["D3"]])
  off <- data.frame(arm = arms, dx = c(0.2, -0.25, -0.3, 0, 0, 0.25, 0.35), dy = c(2.6, -2.6, 2.0, -2.8, 2.8, 2.4, -2.6), hj = c(0, 1, 1, 0.5, 0.5, 0, 0))
  td <- do.call(rbind, lapply(arms, function(a) { gx <- q3(E48[[a]] - E48$B); dx <- q3(100 * DI[[a]])
    data.frame(arm = a, g = gx[["m"]], glo = gx[["lo"]], ghi = gx[["hi"]], d = dx[["m"]], dlo = dx[["lo"]], dhi = dx[["hi"]]) }))
  bx <- q3(100 * DI$B)
  pa <- ggplot(td, aes(g, d, colour = arm)) +
    annotate("rect", xmin = 3, xmax = Inf, ymin = -Inf, ymax = 25, fill = "#E8F1FA") +
    annotate("text", x = 5.9, y = 21, label = "target zone", hjust = 1, size = BASE / .pt, family = FONT, colour = "#0072B2") +
    geom_hline(yintercept = 25, linetype = 2, linewidth = LW, colour = "grey40") + geom_vline(xintercept = 3, linetype = 2, linewidth = LW, colour = "grey40") +
    geom_errorbar(aes(ymin = dlo, ymax = dhi), width = 0, linewidth = LW, alpha = 0.55) +
    geom_errorbar(aes(xmin = glo, xmax = ghi), width = 0, linewidth = LW, orientation = "y", alpha = 0.55) +
    geom_point(size = 1.9) +
    annotate("point", x = 0, y = bx[["m"]], shape = 21, size = 2, fill = ARMCOL[["B"]], colour = ARMCOL[["B"]]) +
    annotate("text", x = -0.25, y = bx[["m"]] + 3, label = "B (reference)", hjust = 1, size = (BASE - 0.5) / .pt, family = FONT, colour = ARMCOL[["B"]]) +
    geom_text(data = merge(td, off), aes(x = g + dx, y = d + dy, label = alab[arm], hjust = hj), size = (BASE - 0.5) / .pt, family = FONT, show.legend = FALSE) +
    scale_colour_manual(values = acol, guide = "none") + scale_x_continuous(limits = c(-2.4, 6.2), breaks = seq(-2, 6, 2)) + scale_y_continuous(limits = c(15, 56)) +
    labs(x = "EDCT gain over B (months)", y = "Stopped at least one drug for AEs (%)") + theme_ms()
  # discontinuation by drug
  dc <- utils::read.csv(file.path(RES, "ensemble_disc_components.csv"), stringsAsFactors = FALSE)
  comp <- c(any = "any drug", backbone = "backbone (CDK4/6i\nor endocrine agent)", eve = "everolimus", cdk = "CDK4/6i", et = "endocrine agent")
  ca <- c("A", "B", "C", "D2")
  dd <- do.call(rbind, lapply(ca, function(a) do.call(rbind, lapply(names(comp), function(k) { x <- q3(100 * dc[[k]][dc$arm == a])
    data.frame(arm = a, comp = comp[[k]], m = x[["m"]], lo = x[["lo"]], hi = x[["hi"]]) }))))
  dd$comp <- factor(dd$comp, levels = rev(unname(comp))); dd$arm <- factor(dd$arm, levels = ca)
  pb <- ggplot(dd, aes(m, comp, colour = arm)) +
    geom_vline(xintercept = 25, linetype = 2, linewidth = LW, colour = "grey40") +
    geom_errorbar(aes(xmin = lo, xmax = hi), width = 0, linewidth = LW, orientation = "y", position = position_dodge(width = 0.7, reverse = TRUE)) +
    geom_point(size = 1.5, position = position_dodge(width = 0.7, reverse = TRUE)) +
    scale_colour_manual(values = ARMCOL[ca], labels = c(A = "A", B = "B", C = "C", D2 = "D2"), name = NULL) +
    labs(x = "Stopped for AEs (%)", y = NULL) + theme_ms() + theme(legend.position = "top")
  # decision by efficacy bar and tolerability limit
  dE <- list(`real-world attrition` = E48$A - E48$B, `trial-like attrition (15%)` = E48$A - E48$B_trial)
  O60 <- piv("OS_RMST60"); dO <- list(`real-world attrition` = O60$A - O60$B, `trial-like attrition (15%)` = O60$A - O60$B_trial)
  discA <- stats::median(DI$A)
  bars <- seq(1, 4, 0.5); lims <- seq(25, 55, 5)
  gr <- expand.grid(bar = bars, lim = lims, att = names(dE), stringsAsFactors = FALSE)
  N <- load_json("report_numbers.json")
  gr$dec <- mapply(function(b, l, a) {
    fav <- mean(dE[[a]] >= b) >= 0.80 && mean(dO[[a]] > 0) >= 0.70 && discA <= l / 100 && isTRUE(N$qtwist_frac_positive >= 0.999)
    if (fav) "favoured" else if (mean(dE[[a]] <= 0) >= 0.5 || discA > l / 100) "not favoured" else "uncertain" }, gr$bar, gr$lim, gr$att)
  gr$dec <- factor(gr$dec, levels = c("favoured", "uncertain", "not favoured")); gr$att <- factor(gr$att, levels = names(dE))
  pc <- ggplot(gr, aes(factor(bar), factor(lim), fill = dec)) + geom_tile(colour = "white", linewidth = LW) +
    geom_point(data = data.frame(bar = "3", lim = "25", att = factor(names(dE), levels = names(dE))), aes(x = bar, y = lim), inherit.aes = FALSE,
               shape = 4, size = 2.2, stroke = 0.8) +
    facet_wrap(~att) + scale_fill_manual(values = c(favoured = "#0072B2", uncertain = "#E69F00", `not favoured` = "#D9D9D9"), name = NULL, drop = FALSE) +
    labs(x = "Efficacy bar (months)", y = "Tolerability limit, any drug (%)") + theme_ms() + theme(legend.position = "right", axis.line = element_blank())
  p <- tagged(((pa | pb) + plot_layout(widths = c(1.15, 1))) / pc + plot_layout(heights = c(1.15, 0.85)))
  save_ms(p, "Figure_5", 180, 150)
}

# ---------------------------------------------------------------------------- Figure 4: where the gain comes from
fig5 <- function() {
  g <- function(a, b = "B") E48[[a]] - E48[[b]]
  inter <- (E48$A - E48$B) - (E48$D1 - E48$B) - (E48$D2 - E48$B)
  fd <- list(g("D1"), g("D2"), inter, g("A"), g("A_let"), g("D3"))
  names(fd) <- c(paste0("SERD first (D1 ", MINUS, " B)"), paste0("Everolimus first (D2 ", MINUS, " B)"), "Interaction",
                 paste0("Triplet (A ", MINUS, " B)"), paste0("Triplet with letrozole ", MINUS, " B"), paste0("ctDNA-guided switch (D3 ", MINUS, " B)"))
  fdd <- do.call(rbind, lapply(names(fd), function(k) { x <- q3(fd[[k]]); data.frame(lab = k, m = x[["m"]], lo = x[["lo"]], hi = x[["hi"]]) }))
  fdd$lab <- factor(fdd$lab, levels = rev(names(fd)))
  pa <- ggplot(fdd, aes(m, lab)) + geom_vline(xintercept = 0, linewidth = LW, colour = "grey55") +
    geom_errorbar(aes(xmin = lo, xmax = hi), width = 0.25, linewidth = LW, orientation = "y") + geom_point(size = 1.8, colour = "#D55E00") +
    labs(x = "EDCT difference (months)", y = NULL) + theme_ms()
  # clonal composition at strategy failure
  fe <- piv("frac_E_at_EDCT"); fr <- piv("frac_R_at_EDCT")
  cl <- do.call(rbind, lapply(names(ARMCOL), function(a) rbind(
    { x <- q3(100 * fe[[a]]); data.frame(arm = a, what = "<i>ESR1</i>-mutant cells", m = x[["m"]], lo = x[["lo"]], hi = x[["hi"]]) },
    { x <- q3(100 * fr[[a]]); data.frame(arm = a, what = "CDK4/6i-escape cells", m = x[["m"]], lo = x[["lo"]], hi = x[["hi"]]) })))
  cl$arm <- factor(cl$arm, levels = rev(names(ARMCOL)))
  cl$what <- factor(cl$what, levels = c("CDK4/6i-escape cells", "<i>ESR1</i>-mutant cells"))
  pb <- ggplot(cl, aes(m, arm, colour = arm)) + facet_wrap(~what, scales = "free_x") +
    geom_errorbar(aes(xmin = lo, xmax = hi), width = 0, linewidth = LW, orientation = "y") + geom_point(size = 1.8) +
    scale_colour_manual(values = ARMCOL, guide = "none") + labs(x = "Share of the tumour at strategy failure (%)", y = NULL) +
    theme_ms() + theme(panel.spacing.x = unit(7, "mm"), strip.text = ggtext::element_markdown())
  # attrition
  gs <- utils::read.csv(file.path(RES, "gsa_samples.csv"), check.names = FALSE)
  fit <- stats::lm(dEDCT48_AB ~ tum.p_attr, data = gs)
  pc <- ggplot(gs, aes(100 * tum.p_attr, dEDCT48_AB)) +
    annotate("rect", xmin = 13, xmax = 19, ymin = -Inf, ymax = Inf, fill = "#E8F1FA") +
    annotate("text", x = 16, y = max(gs$dEDCT48_AB) + 0.3, label = "SONIA", size = (BASE - 0.5) / .pt, family = FONT, colour = "#0072B2") +
    geom_vline(xintercept = 29, linetype = 2, linewidth = LW, colour = "#D55E00") +
    annotate("text", x = 30, y = min(gs$dEDCT48_AB), label = "ESME-MBC", hjust = 0, size = (BASE - 0.5) / .pt, family = FONT, colour = "#D55E00") +
    geom_hline(yintercept = 0, linewidth = LW, colour = "grey70") +
    geom_point(size = 0.9, alpha = 0.6, colour = "grey30") + geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = LW + 0.1, colour = "black") +
    labs(x = "Attrition before line 2 (%; progression\nafter more than 12 months on line 1)", y = "EDCT gain, A minus B (months)") + theme_ms()
  p <- tagged(((pc | pa) + plot_layout(widths = c(1, 1.25))) / pb + plot_layout(heights = c(1, 0.85)))   # a attrition, b factorial, c clones
  save_ms(p, "Figure_4", 180, 125)
}

# ---------------------------------------------------------------------------- Figure 6: who benefits and trial design
fig6 <- function() {
  sg <- utils::read.csv(file.path(RES, "ensemble_subgroups.csv"), stringsAsFactors = FALSE)
  keep <- c("All", "Proliferation high tercile", "Proliferation low tercile", "Endocrine dependence high (b0 low tercile)", "Endocrine dependence low (b0 high tercile)",
            "PI3K pathway wild-type", "PI3K pathway altered", "Tumor DNA fraction >= 10%", "Tumor DNA fraction < 10%", "Baseline ESR1 detectable",
            "Age >= 65", "Age < 65", "AE susceptibility high tercile")
  nice <- c("All patients", "Proliferation, fastest third", "Proliferation, slowest third", "Endocrine dependence, highest third",
            "Endocrine dependence, lowest third", "PI3K pathway wild-type", "PI3K pathway altered", paste0("Tumour DNA fraction ", GE, " 10%"),
            "Tumour DNA fraction &lt; 10%", "<i>ESR1</i> mutation detectable at baseline", paste0("Age ", GE, " 65 years"), "Age &lt; 65 years",
            "AE susceptibility, highest third")
  sd <- do.call(rbind, lapply(seq_along(keep), function(i) do.call(rbind, lapply(c("A-B", "A-B_trial"), function(cn) {
    x <- q3(sg$dEDCT48[sg$subgroup == keep[i] & sg$contrast == cn])
    data.frame(sub = nice[i], att = if (cn == "A-B") names(ATTCOL)[1] else names(ATTCOL)[2], m = x[["m"]], lo = x[["lo"]], hi = x[["hi"]]) }))))
  sd$sub <- factor(sd$sub, levels = rev(nice)); sd$att <- factor(sd$att, levels = names(ATTCOL))
  pa <- ggplot(sd, aes(m, sub, colour = att)) + geom_vline(xintercept = 3, linetype = 2, linewidth = LW, colour = "grey40") +
    geom_vline(xintercept = 0, linewidth = LW, colour = "grey70") +
    geom_errorbar(aes(xmin = lo, xmax = hi), width = 0, linewidth = LW, orientation = "y", position = position_dodge(width = 0.6, reverse = TRUE)) +
    geom_point(size = 1.6, position = position_dodge(width = 0.6, reverse = TRUE)) +
    scale_colour_manual(values = ATTCOL, name = NULL) + labs(x = "EDCT gain, A minus B (months)", y = NULL) +
    theme_ms() + theme(legend.position = "top", axis.text.y = ggtext::element_markdown())
  pw <- utils::read.csv(file.path(RES, "power.csv"), stringsAsFactors = FALSE)
  des <- c("A vs B" = "all patients, real-world attrition", "A vs B, trial-like attrition 15%" = "all patients, trial-like attrition",
           "A vs B, trial-like 15%, tumor DNA >= 10% only" = paste0("tumour DNA ", GE, " 10%, trial-like"),
           "A vs B, trial-like 15%, proliferation high tercile only" = "fastest third, trial-like")
  pd <- pw[pw$endpoint == "EDCT_logrank" & pw$contrast %in% names(des), ]
  pd$design <- factor(des[pd$contrast], levels = unname(des))
  pb <- ggplot(pd, aes(n_per_arm, 100 * power, colour = design)) + geom_hline(yintercept = 80, linetype = 2, linewidth = LW, colour = "grey40") +
    geom_line(linewidth = LW + 0.1) + geom_point(size = 1.3) +
    scale_x_log10(breaks = c(100, 200, 400, 800, 1600)) + scale_y_continuous(limits = c(0, 100)) +
    scale_colour_manual(values = c("#D55E00", "#0072B2", "#999999", "#009E73"), name = NULL) +
    labs(x = "Patients per arm (log scale)", y = "Power, EDCT log-rank test (%)") +
    theme_ms() + theme(legend.position = "inside", legend.position.inside = c(0.98, 0.03), legend.justification = c(1, 0))
  p <- tagged((pa | pb) + plot_layout(widths = c(1.3, 1)))
  save_ms(p, "Figure_6", 180, 85)
}

FIGS <- list(fig1, fig2, fig3, fig5, fig4, fig6)   # fig5() draws Figure 4 (mechanism), fig4() draws Figure 5 (toxicity)
args <- commandArgs(trailingOnly = TRUE)
which_ <- if (length(args)) as.integer(args) else seq_along(FIGS)
for (i in which_) tryCatch(FIGS[[i]](), error = function(e) cat(sprintf("figure %d failed: %s\n", i, conditionMessage(e))))
