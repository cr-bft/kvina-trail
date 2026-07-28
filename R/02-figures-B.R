#!/usr/bin/env Rscript
# R/02-figures-B.R -- additional manuscript and supplement figures. Run this in a 
# separate session from 02-figures-A.R. Sources R/01-stats.R
# (which sources R/functions-paper.R and regenerates the processed tables), so
# `Rscript R/paper-02-figures.R` from the repo root reproduces select figures: 
# the weathering violin and toc into figures/manuscript, and supplement
# figures S1-S7 into figures/supporting-information. (Assign the manuscript
# figure number for the violin when renumbering.)

source("R/01-analysis.R")                 # -> R/functions-paper.R + processed tables

library("patchwork")                     # 25 layouts + 33's S8 composite
library("mgcv")                          # 33's S8 seasonal-GAM climatology
library("lubridate")                     # year()/week() in the TOC pH record
# readxl, data.table, ggplot2, ggtext are attached upstream by functions-paper.R.

fig_manuscript <- file.path("figures", "manuscript")
fig_supp       <- file.path("figures", "supporting-information")
dir.create(fig_manuscript, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_supp, recursive = TRUE, showWarnings = FALSE)

# FIGURE 2 — trial-chemistry weathering scatter + Historic-vs-Trial split violin ----
fig_dir <- fig_manuscript

STN_LEV  <- c("NY", "OB", "LA", "KL")
STN_LAB  <- c(NY = "NY (upstream)", OB = "OB (below doser)",
              LA = "LA (Litleåne)", KL = "KL (mouth)")
STN_COLS <- setNames(unname(SITE_COLS[STN_LEV]), STN_LAB[STN_LEV])

# DOC-corrected carbonate alkalinity and 2(Ca+Mg); NA where DOC is absent:
doc_correct <- function(dt) {
  dt[, ph_for := fifelse(is.finite(ph), ph, DOC_CORR_DEFAULT_PH)]
  dt[, carb_alk_doccorr_mmol_L := fifelse(
    is.finite(alk_meas_mmol_L) & is.finite(doc_for_corr_mg_L),
    carb_alk_from_titration_ueq_L(alk_meas_mmol_L * 1000, ph_for,
                                  doc_for_corr_mg_L) / 1000, NA_real_)]
  dt[, camg_mmol_L := 2 * (ca_mg_L / MW_CA + mg_mg_L / MW_MG)]
  dt[]
}

# 2025 trial: 4 stations, measured DOC
tr <- load_trial_stations()
tr[, doc_for_corr_mg_L := doc_mg_L]
tr <- doc_correct(tr)
tr[, stn := factor(STN_LAB[as.character(station)], levels = STN_LAB)]

# season-match window: first grab -> Sep 15 (TRIAL_WINDOW in functions.R):
trial_start <- min(tr$ts, na.rm = TRUE)
trial_end   <- TRIAL_WINDOW$end
tr <- tr[ts <= trial_end]

# Historical: 4 analogues, NO DOC / NO carbonate-alk correction
HIST_STN <- c("81555" = "NY", "58892" = "OB", "58891" = "LA", "58890" = "KL")
HP <- c(pH = "ph", "Syrenøytraliserende kapasitet (ANC)" = "anc_ueq_L",
        "Total alkalitet" = "alk_meas_mmol_L",   # measured Alk = the DOC-correction base
        Kalsium = "ca_mg_L", Magnesium = "mg_mg_L",
        Konduktivitet = "cond_mS_m", Sulfat = "so4_mg_L",
        Nitrat = "no3n_ugN_L",
        "Totalt organisk karbon (TOC)" = "toc_mg_L")
hraw <- read_vannmiljo_long()
hraw <- hraw[as.character(water_location_id) %in% names(HIST_STN) &
               parameter %in% names(HP) & is.finite(value) &
               as.Date(timestamp) < as.Date("2025-01-01") &
               in_trial_date_window(as.Date(timestamp))]
hraw[, `:=`(ts = as.Date(timestamp), col = HP[parameter], stn3 = HIST_STN[as.character(water_location_id)])]
hw <- dcast(hraw, stn3 + ts ~ col, value.var = "value", fun.aggregate = function(z) mean(z, na.rm = TRUE))
for (nm in unname(HP)) if (!nm %in% names(hw)) hw[, (nm) := NA_real_]
hw[, doc_for_corr_mg_L := NA_real_]
hw <- doc_correct(hw)
hw[, station := factor(stn3, levels = STN_LEV)]
hw[, stn := factor(STN_LAB[as.character(stn3)], levels = STN_LAB)]

th <- function(base = 12)
  theme_kvina_md(base_size = base) +
  theme(plot.tag = element_text(face = "bold", size = base + 2),
        legend.text = element_markdown(size = base - 1),
        axis.title = element_markdown(size = base),
        strip.text = element_markdown(face = "bold", size = base - 1))

# (b) historical-vs-trial box-and-whisker, built from tr/hw (all 4 stations) --
tr[, period := "Trial"]
hw[, period := "Historical"]

# Nitrate to mg/L (NO3-N basis) to match the stats tables.
tr[, no3n_mgN_L := no3n_ugN_L / 1000]
hw[, no3n_mgN_L := if ("no3n_ugN_L" %in% names(hw)) no3n_ugN_L / 1000 else NA_real_]

PARAM_LAB <- c(
  ph         = "pH",
  anc_ueq_L  = "ANC (µeq L<sup>-1</sup>)",
  ca_mg_L    = "Calcium (mg L<sup>-1</sup>)",
  mg_mg_L    = "Magnesium (mg L<sup>-1</sup>)",
  cond_mS_m  = "Conductivity (mS m<sup>-1</sup>)",
  so4_mg_L   = "Sulphate (mg L<sup>-1</sup>)",
  no3n_mgN_L = "Nitrate-N (mg L<sup>-1</sup>)")
box_cols <- names(PARAM_LAB)
for (cc in box_cols) {
  if (!cc %in% names(tr)) tr[, (cc) := NA_real_]
  if (!cc %in% names(hw)) hw[, (cc) := NA_real_]
}
box_long <- rbind(tr[, c("station", "period", box_cols), with = FALSE],
                  hw[, c("station", "period", box_cols), with = FALSE])
box_long <- melt(box_long, id.vars = c("station", "period"), measure.vars = box_cols,
                 variable.name = "col", value.name = "value")[is.finite(value)]
box_long[, col := as.character(col)]
# QC: drop one erroneous NY historical ANC grab (2024-07-01, -140 µeq/L); figures only.
box_long <- box_long[!(col == "anc_ueq_L" & value < -50)]

# Fade a station colour for the HISTORICAL boxes by DESATURATING.
fade <- function(col, sat_mult = 0.4, val_mult = 1.08) {
  hsvc <- grDevices::rgb2hsv(grDevices::col2rgb(col))
  grDevices::hsv(hsvc[1], min(1, hsvc[2] * sat_mult), min(1, hsvc[3] * val_mult))
}
site_fill_palette <- function(sites)
  do.call(c, lapply(sites, function(s) setNames(
    c(fade(SITE_COLS[[s]]), unname(SITE_COLS[[s]])), paste(s, c("Historical", "Trial"), sep = "_"))))
interleave_levels <- function(sites, suffixes, sep) as.vector(t(outer(sites, suffixes, paste, sep = sep)))

# (b-violin) Split-violin variant of panel (b): one violin PER STATION, split by
# period -- Historical = left half, Trial = right half.
build_violin_panel <- function(sites, param_keys, ncol = NULL, nrow = NULL,
                               hist_lab = "Hist", trial_lab = "Trial",
                               spatial_bracket = TRUE) {
  plab <- c(Historical = hist_lab, Trial = trial_lab)   # internal -> displayed
  plots <- lapply(seq_along(param_keys), function(i) {
    ck  <- param_keys[i]
    lab <- unname(PARAM_LAB[[ck]])
    dp <- box_long[col == ck & as.character(station) %in% sites]  # ANC -140 QC done at box_long
    dp[, param := factor(lab, levels = lab)]
    avail <- dp[, .(n_hist = sum(period == "Historical"), n_trial = sum(period == "Trial")),
                by = station]
    keep <- sites[sites %in% avail[n_hist > 0 & n_trial > 0, as.character(station)]]
    dp <- dp[as.character(station) %in% keep]
    dp[, station := factor(as.character(station), levels = keep)]
    flevs <- interleave_levels(keep, c("Historical", "Trial"), "_")
    dp[, fill_key    := factor(paste(station, period, sep = "_"), levels = flevs)]
    dp[, period_disp := factor(unname(plab[period]), levels = unname(plab))]
    
    rng  <- diff(range(dp$value, na.rm = TRUE))
    ymax <- max(dp$value, na.rm = TRUE)
    
    # (1) TEMPORAL stars -- Historical vs Trial, per station (grey, over the violin)
    sig <- dp[, {
      tv <- value[period == "Trial"]; hv <- value[period == "Historical"]
      .(stars = star(welch_test(tv, hv)$p_value))   # shared guard: n >= 3 per group
    }, by = station]
    sig_lab <- sig[stars != "" & stars != "ns"]
    sig_lab[, `:=`(xpos = match(as.character(station), keep), ystar = ymax + 0.08 * rng)]
    
    # (2) SPATIAL bracket -- upstream (NY) vs mouth (KL) DURING the trial.
    up <- "NY"; dn <- if ("KL" %in% keep) "KL" else keep[length(keep)]
    brk <- NULL
    if (spatial_bracket && up %in% keep && dn %in% keep && up != dn) {
      uv <- dp[station == up & period == "Trial", value]
      dv <- dp[station == dn & period == "Trial", value]
      sp_star <- star(welch_test(uv, dv)$p_value)   # NA p (guarded) -> "" star
      if (!sp_star %in% c("", "ns")) {
        xu <- match(up, keep); xd <- match(dn, keep)
        yb <- ymax + 0.21 * rng   # above the temporal stars (0.08) so ticks clear them
        brk <- data.frame(x = xu, xend = xd, y = yb, ytick = yb - 0.025 * rng,
                          xmid = (xu + xd) / 2, ystar = yb + 0.013 * rng, lab = sp_star)
      }
    }
    bracket_layers <- if (!is.null(brk)) list(
      geom_segment(data = brk, aes(x = x, xend = xend, y = y, yend = y),
                   inherit.aes = FALSE, linewidth = 0.35, colour = "grey20"),
      geom_segment(data = brk, aes(x = x, xend = x, y = ytick, yend = y),
                   inherit.aes = FALSE, linewidth = 0.35, colour = "grey20"),
      geom_segment(data = brk, aes(x = xend, xend = xend, y = ytick, yend = y),
                   inherit.aes = FALSE, linewidth = 0.35, colour = "grey20"),
      geom_text(data = brk, aes(x = xmid, y = ystar, label = lab), inherit.aes = FALSE,
                colour = "grey10", fontface = "bold", size = 4.5, vjust = 0)) else NULL
    top_expand <- if (!is.null(brk)) 0.27 else 0.17
    
    period_guide <- guide_legend(order = 2, override.aes = list(
      fill = c("grey80", "grey45"), colour = "grey45", linewidth = 0.2))
    
    ggplot(dp, aes(station, value, fill = fill_key, group = as.integer(fill_key))) +
      geom_split_violin(aes(alpha = period_disp), scale = "width", trim = TRUE,
                        colour = "grey45", linewidth = 0.2,
                        med_colour = "grey35", med_linewidth = 0.35) +
      geom_text(data = sig_lab, aes(x = xpos, y = ystar, label = stars), inherit.aes = FALSE,
                colour = "grey50", fontface = "bold", size = 4.2, vjust = 0.3) +
      bracket_layers +
      facet_wrap(~ param) +
      scale_fill_manual(values = site_fill_palette(keep), guide = "none") +
      scale_alpha_manual(values = setNames(c(1, 1), unname(plab)), name = "Period",
                         drop = FALSE, guide = period_guide) +
      scale_y_continuous(expand = expansion(mult = c(0.05, top_expand))) +
      labs(x = NULL, y = NULL, tag = if (i == 1) "b" else NULL) +
      th() +
      theme(axis.text.x = element_text(size = 9, colour = "grey25"),
            axis.ticks.x = element_blank(),
            panel.spacing = grid::unit(4, "mm"))
  })
  patchwork::wrap_plots(plots, nrow = nrow, ncol = ncol)
}

# (a) weathering scatter for a period frame ----------------------------------
weathering_panel <- function(d, xlab = "DOC-corrected trial Alkalinity (meq L<sup>-1</sup>)",
                             wedge = TRUE, colour_title = NULL, point_size = 1.9,
                             centroids = FALSE, y_basis = c("camg", "ca"),
                             point_alpha = NULL) {
  y_basis <- match.arg(y_basis)
  d <- copy(d)
  d[, y_cat := if (y_basis == "ca") 2 * ca_mg_L / MW_CA else camg_mmol_L]
  ylab <- if (y_basis == "ca") "2Ca (meq L<sup>-1</sup>)"
  else "2(Ca+Mg) (meq L<sup>-1</sup>)"
  dd <- d[is.finite(carb_alk_doccorr_mmol_L) & is.finite(y_cat) & !is.na(stn)]
  lim <- max(c(dd$carb_alk_doccorr_mmol_L, dd$y_cat), na.rm = TRUE) * 1.05
  cen <- dd[, .(cx = median(carb_alk_doccorr_mmol_L, na.rm = TRUE),
                cy = median(y_cat, na.rm = TRUE)), by = stn]
  raw_alpha <- if (!is.null(point_alpha)) point_alpha else if (centroids) 0.30 else 0.85
  cen_layers <- if (centroids) list(
    geom_point(data = cen, aes(cx, cy), inherit.aes = FALSE,
               colour = "grey20", size = 5.6),
    geom_point(data = cen, aes(cx, cy, colour = stn), inherit.aes = FALSE,
               size = 4.4, show.legend = FALSE)) else NULL
  ggplot(dd, aes(carb_alk_doccorr_mmol_L, y_cat)) +
    (if (wedge) geom_ribbon(data = data.frame(x = seq(0, lim, length.out = 256)),
                            aes(x, ymin = pmin(x, lim), ymax = pmin(2 * x, lim)),
                            inherit.aes = FALSE, fill = "grey92")) +
    geom_abline(slope = 1, colour = "grey55", linewidth = 0.4) +
    geom_abline(slope = 2, colour = "grey55", linewidth = 0.4, linetype = "22") +
    annotate("text", x = lim * 0.97, y = lim * 0.9,  label = "1:1", colour = "grey45",
             size = 3.4 * (point_size / 1.9)) +
    annotate("text", x = lim * 0.40, y = lim * 0.96, label = "2:1", colour = "grey45",
             size = 3.4 * (point_size / 1.9)) +
    geom_point(aes(colour = stn), size = point_size, alpha = raw_alpha) +
    cen_layers +
    scale_colour_manual(values = STN_COLS, name = colour_title, drop = FALSE) +
    coord_cartesian(xlim = c(0, lim), ylim = c(0, lim)) +
    labs(x = xlab, y = ylab, tag = "a") +
    th() +
    (if (centroids)
      guides(colour = guide_legend(override.aes = list(alpha = 1, size = 3))))
}

# Violin variant of the clean vertical layout
build_vert_violin <- function(d, file, param_keys, sites = c("NY", "KL"),
                              y_basis = "camg", centroids = FALSE, base = 12,
                              xlab = "DOC-corrected trial Alkalinity (meq L<sup>-1</sup>)",
                              spatial_bracket = TRUE) {
  pt_bump <- base / 12
  pa <- weathering_panel(d, wedge = FALSE, colour_title = "Site", point_size = 2.0 * pt_bump,
                         y_basis = y_basis, centroids = centroids, xlab = xlab,
                         point_alpha = 1) +   # match the full-colour violin fills below
    guides(colour = guide_legend(order = 1, override.aes = list(size = 2.55 * pt_bump, alpha = 1)))
  pb <- build_violin_panel(sites = sites, param_keys = param_keys, nrow = 2,
                           hist_lab = "Historic", spatial_bracket = spatial_bracket)
  p <- (pa / pb + plot_layout(heights = c(1, 1.15), guides = "collect")) &
    theme_cdr(base = base) &
    theme(legend.position = "right", legend.justification = "center",
          legend.text = element_markdown(colour = "grey15", size = (base + 1) * 0.8),
          strip.text = element_markdown(face = "plain", colour = "grey10", size = base - 3,
                                        margin = margin(t = 5, b = 5)),
          axis.text.x = element_text(size = base - 3.5, colour = "grey35"),
          plot.tag = element_text(face = "bold", size = base + 6, colour = "grey5"),
          panel.spacing = grid::unit(3, "mm"))
  p[[1]] <- p[[1]] + theme(
    axis.title.x = element_markdown(size = base - 1, colour = "grey10"),
    axis.title.y = element_markdown(size = base - 1, colour = "grey10"))
  ggsave(file.path(fig_dir, file), p, width = 9.4, height = 11.2, dpi = 450, bg = "white")
  cat(sprintf("wrote %s (vertical violin variant; base = %g)\n", file, base))
}

# CDR-reference-style theme
theme_cdr <- function(base = 12) {
  theme_kvina_md(base_size = base) +
    theme(
      panel.grid.major = element_line(colour = "grey92", linewidth = 0.25),
      panel.grid.minor = element_blank(),
      panel.border     = element_rect(colour = "grey20", fill = NA, linewidth = 0.5),
      strip.background = element_rect(fill = "grey88", colour = "grey20", linewidth = 0.5),
      strip.text       = element_markdown(face = "plain", colour = "grey10", size = base,
                                          margin = margin(t = 5, b = 5)),
      axis.title.x = element_markdown(face = "plain", colour = "grey10", size = base + 1),
      axis.title.y = element_markdown(face = "plain", colour = "grey10", size = base + 1),
      axis.text    = element_text(colour = "grey35", size = base - 1),
      plot.tag     = element_text(face = "bold", size = base + 8, colour = "grey5"),
      legend.title = element_markdown(face = "plain", colour = "grey10",
                                      size = base + 1, hjust = 0),
      legend.text  = element_markdown(colour = "grey15", size = base),
      legend.key   = element_blank(),
      legend.title.position = "top",
      legend.key.spacing.y  = grid::unit(1.2, "mm"),
      legend.spacing.y      = grid::unit(7, "mm")
    )
}

# Figure 2: panel (a) 2Ca weathering scatter, panel (b) Historical-vs-trial split
# violin (with the NY->KL spatial bracket).
build_vert_violin(tr, "weathering-violin.png",
                  c("anc_ueq_L", "ca_mg_L", "so4_mg_L", "no3n_mgN_L"),
                  y_basis = "ca", centroids = FALSE, base = 16,
                  xlab = "HCO<sub>3</sub><sup>−</sup> (meq L<sup>-1</sup>)")

cat("done.\n")

# SUPPLEMENT CHEMISTRY TIME SERIES + SALT PULSE (S1-S6, S8) ------------------
Sys.setlocale("LC_CTYPE", "en_US.UTF-8")   # non-UTF-8 locales drop the µ glyph in ragg output

fig_dir <- fig_supp

TRIAL_XLSX <- trial_xlsx_path()
BAND_LAB   <- "Trial Period (Jun 10 - Sep 15)"

# SHARED: read the 2025 trial WaterChem sheet. -------------------------------
read_trial_wc <- function() {
  raw <- as.data.table(read_waterchem_sorted())
  setnames(raw, names(raw), sub("\\n.*$", "", names(raw)))
  raw[, ts := as.Date(`Sample date`)]
  raw[]
}

# SHARED: build one stacked-facet time-series figure and export it. ----------
build_supp_figure <- function(dt, station_cols, x_scale, linewidth,
                              band = FALSE, angle_x = FALSE, height = 11.5, out_name) {
  p <- ggplot(dt, aes(ts, value, colour = station, group = grp))
  if (band) {
    bd <- data.frame(xmin = STUDY0, xmax = STUDY1, ymin = -Inf, ymax = Inf,
                     fill_key = BAND_LAB)
    p <- p + geom_rect(data = bd, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
                                      fill = fill_key),
                       inherit.aes = FALSE, alpha = 0.55)
  }
  p <- p +
    geom_line(linewidth = linewidth, alpha = 0.9) +
    facet_wrap(~ param, ncol = 1, scales = "free_y", strip.position = "top") +
    scale_colour_manual(values = station_cols, name = "Site") +
    x_scale +
    labs(x = NULL, y = NULL) +
    theme_kvina_md(base_size = 12) +
    theme(legend.position = "right", legend.justification = "center",
          legend.key.size = grid::unit(1.0, "lines"),
          panel.spacing = grid::unit(4, "mm"))
  if (band) {
    p <- p +
      scale_fill_manual(values = setNames("grey85", BAND_LAB), name = "Period") +
      guides(colour = guide_legend(override.aes = list(linewidth = 1.3), order = 1),
             fill   = guide_legend(override.aes = list(alpha = 0.55), order = 2))
  } else {
    p <- p + guides(colour = guide_legend(override.aes = list(linewidth = 1.3)))
  }
  if (angle_x) p <- p + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave(file.path(fig_dir, out_name), p, width = 9.4, height = height,
         dpi = 450, bg = "white", device = ragg::agg_png)
  cat(sprintf("wrote figures/%s\n", out_name))
}

# FULL-RECORD family (S1-S3): monitoring record 2011-2025 --------------------
REC_START <- as.Date("2011-01-01")   # monitoring record window (inclusive)
REC_END   <- as.Date("2025-12-31")
STUDY0    <- TRIAL_WINDOW$start       # trial quantification period (shaded band)
STUDY1    <- TRIAL_WINDOW$end          # Jun 10 - Sep 15 2025 (CDR quantification window)
AXIS_END  <- as.Date("2026-01-01")   # a little past the record end
GAP_DAYS  <- 180L                    # break a line across a real data gap
MIN_SERIES_PTS <- 4L                 # drop a station's series in a panel below this

X_SCALE_FULL <- scale_x_date(date_breaks = "1 year", date_labels = "%Y",
                             limits = c(REC_START, AXIS_END))

# Stations
NY_LAB    <- "NY — Nyland (upstream)"
KL_LAB    <- "KL — Kloster (mouth)"
LA_LAB    <- "LA — Litleåne (limed tributary)"
LITUP_LAB <- "Litleåne (upstream of doser)"
STN <- c("81555" = NY_LAB,
         "58890" = KL_LAB,
         "58891" = LA_LAB,
         "45779" = LITUP_LAB)
TR_STN <- c("74911" = NY_LAB, "74914" = KL_LAB, "74913" = LA_LAB)
ST_COLS <- c(setNames(unname(SITE_COLS["NY"]), NY_LAB),
             setNames(unname(SITE_COLS["KL"]), KL_LAB),
             setNames(unname(SITE_COLS["LA"]), LA_LAB),
             setNames("grey60", LITUP_LAB))

# Parameter registry: historical name -> markdown label -> trial token
PARAM_META <- data.table(
  hist  = c("Totalt organisk karbon (TOC)", "Sulfat", "Nitrat", "Totalfosfor",
            "Totalnitrogen",
            "Kalsium", "Magnesium", "Kalium", "Klorid", "Konduktivitet",
            "pH", "Syrenøytraliserende kapasitet (ANC)",
            "Reaktivt aluminium", "Ikke-labilt aluminium", "Labilt aluminium (derived)"),
  label = c("TOC (mg C L<sup>-1</sup>)",
            "Sulphate (mg L<sup>-1</sup>)",
            "Nitrate-N (mg L<sup>-1</sup>)",
            "Total phosphorus (µg P L<sup>-1</sup>)",
            "Total nitrogen (µg N L<sup>-1</sup>)",
            "Calcium (mg L<sup>-1</sup>)",
            "Magnesium (mg L<sup>-1</sup>)",
            "Potassium (mg L<sup>-1</sup>)",
            "Chloride (mg L<sup>-1</sup>)",
            "Conductivity (mS m<sup>-1</sup>)",
            "pH",
            "ANC (µeq L<sup>-1</sup>)",
            "Reactive Al (µg L<sup>-1</sup>)",
            "Non-labile Al (µg L<sup>-1</sup>)",
            "Labile Al (µg L<sup>-1</sup>)"),
  scale = c(1, 1, 0.001, 1, 1,
            1, 1, 1, 1, 1,
            1, 1,
            1, 1, 1),
  trial = c("TOC", "SO4", "NO3-N", "TOTP", NA,
            "Ca-filt", "Mg-filt", "K-filt", "Cl", "Kond.",
            "pH", "ANC", NA, NA, NA))

# Monitoring long record (2011-2025)
h <- read_vannmiljo_long()
h[, ts := as.Date(timestamp)][, timestamp := NULL]
h <- h[as.character(water_location_id) %in% names(STN) & is.finite(value) &
         !is.na(ts) & ts >= REC_START & ts <= REC_END]

# Coalesce "Nitrat + nitritt" as a fallback for "Nitrat" (repo convention).
no3 <- h[parameter %in% c("Nitrat", "Nitrat + nitritt")]
no3[, pref := fifelse(parameter == "Nitrat", 1L, 2L)]
setorder(no3, water_location_id, ts, pref)
no3 <- unique(no3, by = c("water_location_id", "ts"))
no3[, `:=`(parameter = "Nitrat", pref = NULL)]
h <- rbind(h[!parameter %in% c("Nitrat", "Nitrat + nitritt")], no3)

# Labile Al = Reactive Al - Non-labile Al.
al <- dcast(h[parameter %in% c("Reaktivt aluminium", "Ikke-labilt aluminium")],
            water_location_id + ts ~ parameter, value.var = "value",
            fun.aggregate = function(x) median(x, na.rm = TRUE))
al <- al[is.finite(`Reaktivt aluminium`) & is.finite(`Ikke-labilt aluminium`)]
al[, value := `Reaktivt aluminium` - `Ikke-labilt aluminium`]
al_n <- al[, .N, by = water_location_id]
cat(sprintf("derived labile Al (reactive - non-labile): %s (%d negative, kept)\n",
            paste(al_n$water_location_id, al_n$N, sep = "=", collapse = ", "),
            al[value < 0, .N]))
h <- rbind(h[parameter != "Labilt aluminium"],   # drop redundant direct-labile rows
           al[, .(water_location_id, parameter = "Labilt aluminium (derived)", value, ts)])

# 2025 trial grabs (study period): drop coincident registry rows
wc <- read_trial_wc()
wc <- wc[as.character(`Station id`) %in% names(TR_STN)]
grab <- unique(wc[, .(station = TR_STN[as.character(`Station id`)], ts)])
h[, station := STN[as.character(water_location_id)]]
n0 <- nrow(h)
h <- h[!grab, on = c("station", "ts")][, station := NULL]
cat(sprintf("dropped %d registry rows on trial-grab dates\n", n0 - nrow(h)))

# Build one grouped full-record figure
build_record <- function(param_hist, out_name, height, drop_max = NULL) {
  meta <- PARAM_META[match(param_hist, hist)]
  lev  <- unique(meta$label)
  
  hh <- h[parameter %in% param_hist]
  dt <- hh[, .(ts, station = STN[as.character(water_location_id)],
               param  = meta$label[match(parameter, meta$hist)],
               value  = value * meta$scale[match(parameter, meta$hist)])]
  
  if (!is.null(drop_max)) {                     # cull the single Tot-P outlier
    lab <- PARAM_META[hist == drop_max, label]
    idx <- dt[, .I[param == lab]]
    if (length(idx)) dt <- dt[-idx[which.max(dt$value[idx])]]
  }
  
  dt <- dt[, if (.N >= MIN_SERIES_PTS) .SD, by = .(param, station)]
  setorder(dt, param, station, ts)
  dt[, seg := cumsum(c(TRUE, diff(as.integer(ts)) > GAP_DAYS)), by = .(param, station)]
  dt[, grp := paste(station, seg)]
  
  dt[, param   := factor(param,   levels = lev)]
  dt[, station := factor(station, levels = STN)]
  
  cat(sprintf("%s: %d monitoring obs\n", out_name, nrow(dt)))
  
  build_supp_figure(dt, station_cols = ST_COLS, x_scale = X_SCALE_FULL,
                    linewidth = 0.4, band = TRUE, angle_x = TRUE,
                    height = height, out_name = out_name)
}

# TRIAL family (S4-S6): trial-only window 2025-06-10 -> 2025-09-15 -----------
TRIAL_START <- TRIAL_WINDOW$start
TRIAL_END   <- TRIAL_WINDOW$end

X_SCALE_TRIAL <- scale_x_date(date_breaks = "1 month", date_labels = "%b",
                              limits = c(TRIAL_START, TRIAL_END),
                              expand = expansion(mult = 0.02))

station_labels <- c(
  "74911" = "NY — Nyland (upstream)",
  "74912" = "OB — Below doser",
  "74913" = "LA — Litleåne (limed tributary)",
  "74914" = "KL — Kloster (mouth)"
)
station_cols <- c(
  setNames(unname(SITE_COLS["NY"]), station_labels["74911"]),
  setNames(unname(SITE_COLS["OB"]), station_labels["74912"]),
  setNames(unname(SITE_COLS["LA"]), station_labels["74913"]),
  setNames(unname(SITE_COLS["KL"]), station_labels["74914"])
)

param_meta <- data.table(
  pkey = c("TOC", "DOC", "SO4", "NO3N", "TOTP",
           "Ca", "Mg", "K", "Cl", "Conductivity",
           "pH", "ANC", "Alkalinity", "Al"),
  column = c("TOC", "DOC", "SO4", "NO3-N", "TOTP",
             "Ca", "Mg", "K", "Cl", "Kond.",
             "pH", "ANC", "Alk-filt", "Al"),
  label = c("TOC (mg C L<sup>-1</sup>)", "DOC (mg C L<sup>-1</sup>)",
            "Sulphate (mg L<sup>-1</sup>)", "Nitrate-N (mg L<sup>-1</sup>)",
            "Total phosphorus (µg P L<sup>-1</sup>)",
            "Calcium (mg L<sup>-1</sup>)", "Magnesium (mg L<sup>-1</sup>)",
            "Potassium (mg L<sup>-1</sup>)", "Chloride (mg L<sup>-1</sup>)",
            "Conductivity (mS m<sup>-1</sup>)", "pH", "ANC (µeq L<sup>-1</sup>)",
            "Alkalinity (mmol L<sup>-1</sup>)",
            "Total Al (µg L<sup>-1</sup>)"),
  scale = c(1, 1, 1, 0.001, 1,
            1, 1, 1, 1, 1,
            1, 1, 1, 1)
)

trial <- read_trial_wc()
trial <- trial[as.character(`Station id`) %in% names(station_labels) &
                 ts >= TRIAL_START & ts <= TRIAL_END]
trial[, station := station_labels[as.character(`Station id`)]]

for (col in param_meta$column) {
  if (!col %in% names(trial)) stop(sprintf("Missing expected trial column: %s", col))
  set(trial, j = col, value = .num(trial[[col]]))
}

# Build one grouped trial figure
build_trial <- function(keys, out_name) {
  meta <- param_meta[match(keys, pkey)]
  wide <- trial[, c("ts", "station", meta$column), with = FALSE]
  long <- melt(wide, id.vars = c("ts", "station"), measure.vars = meta$column,
               variable.name = "column", value.name = "value")
  long <- long[is.finite(value)]
  long[, value := value * meta$scale[match(as.character(column), meta$column)]]
  long[, param := meta$label[match(column, meta$column)]]
  long[, `:=`(
    station = factor(station, levels = unname(station_labels)),
    param = factor(param, levels = meta$label)
  )]
  long[, grp := station]                            # one line per station
  setorder(long, param, station, ts)
  
  cat(sprintf("%s: %d observations (%s to %s)\n", out_name, nrow(long),
              min(long$ts), max(long$ts)))
  
  build_supp_figure(long, station_cols = station_cols, x_scale = X_SCALE_TRIAL,
                    linewidth = 0.55, band = FALSE, angle_x = FALSE,
                    height = 11.5, out_name = out_name)
}

# Emit the six supplement figures (S1-S6) ------------------------------------
# S1) Full-record acid / aluminium: pH, ANC, Reactive Al, Non-labile Al, Labile Al
build_record(c("pH", "Syrenøytraliserende kapasitet (ANC)",
               "Reaktivt aluminium", "Ikke-labilt aluminium", "Labilt aluminium (derived)"),
             "S1_acid_aluminium_full_record.png", height = 11.5)
# S2) Full-record major ions: Ca, Mg, K, Cl, Conductivity
build_record(c("Kalsium", "Magnesium", "Kalium", "Klorid", "Konduktivitet"),
             "S2_major_ions_full_record.png", height = 11.5)
# S3) Full-record nutrients & organics: TOC, Sulphate, NO3-N, Total P, Total N
build_record(c("Totalt organisk karbon (TOC)", "Sulfat", "Nitrat", "Totalfosfor",
               "Totalnitrogen"),
             "S3_nutrients_organics_full_record.png", height = 11.5,
             drop_max = "Totalfosfor")

# S4) Trial carbonate / acid / aluminium: pH, ANC, Alkalinity, total Al
build_trial(c("pH", "ANC", "Alkalinity", "Al"),
            "S4_carbonate_acid_aluminium_trial.png")
# S5) Trial major ions: total Ca, total Mg, K, Cl, Conductivity
build_trial(c("Ca", "Mg", "K", "Cl", "Conductivity"),
            "S5_major_ions_trial.png")
# S6) Trial nutrients & organics: TOC, DOC, Sulphate, NO3-N, Total P
build_trial(c("TOC", "DOC", "SO4", "NO3N", "TOTP"),
            "S6_nutrients_organics_trial.png")

# S8) Salt pulse: Cl / Na / Mg / Ca marine tracers at NY / LA / KL framed as a ----
# departure from the 2011-2024 seasonal-GAM climatology.
doy_of <- function(dd) as.integer(format(dd, "%j"))

# constants
TRACERS <- c(Klorid = "Chloride (mg L<sup>-1</sup>)",
             Natrium = "Sodium (mg L<sup>-1</sup>)",
             Konduktivitet = "Conductivity (mS m<sup>-1</sup>)",
             Magnesium = "Magnesium (mg L<sup>-1</sup>)",
             Kalsium = "Calcium (mg L<sup>-1</sup>)")
MARINE3    <- TRACERS[c("Klorid", "Natrium", "Magnesium")]              # pure marine tracers
MARINE3_CA <- TRACERS[c("Klorid", "Natrium", "Magnesium", "Kalsium")]  # + calcium (the liming tracer)
TRIAL0 <- doy_of(TRIAL_WINDOW$start); TRIAL1 <- doy_of(TRIAL_WINDOW$end)
TRIAL_LAB  <- "Trial Period (Jun 10 - Sep 15)"
TRIAL_FILL <- "grey75"
QMONTHS    <- c(1, 4, 7, 10)   # quarterly month ticks (Jan, Apr, Jul, Oct)
mon_breaks <- doy_of(as.Date(sprintf("2023-%02d-01", QMONTHS)))
grid    <- data.table(doy = 1:366)
SRC_LAB <- c("Norwegian (monthly)" = "Monthly registry", "Trial (weekly)" = "Trial (weekly)")
SHP     <- c("Norwegian (monthly)" = 16, "Trial (weekly)" = 17)

# 2025 trial grabs (weekly)
wc <- as.data.table(read_waterchem_sorted())
wc <- wc[!is.na(`Project Id`)]; wc[, d := as.Date(`Sample date`)]

# NY / LA / KL 2025 record + 2011-2024 seasonal-GAM climatology
STN3    <- c("81555" = "NY — upstream (above doser)", "45779" = "LA — Litleåne",
             "58890" = "KL — river mouth")
TR_STN3 <- c("74911" = STN3[["81555"]], "74913" = STN3[["45779"]], "74914" = STN3[["58890"]])
SITE3   <- c("81555" = "NY", "45779" = "LA", "58890" = "KL")
stncol3 <- setNames(unname(SITE_COLS[SITE3[names(STN3)]]), unname(STN3))

h3 <- read_vannmiljo_long()
h3[, d := as.Date(timestamp)]
h3 <- h3[parameter %in% names(TRACERS) & is.finite(value) & value > 0 &
           as.character(water_location_id) %in% names(STN3)]
h3[, `:=`(station = STN3[as.character(water_location_id)], tracer = TRACERS[parameter], doy = doy_of(d))]
h3[, era := fifelse(d < as.Date("2025-01-01"), "clim", "y2025")]
cur3 <- h3[era == "y2025", .(station, tracer, doy, d, value, src = "Norwegian (monthly)")]

tr3 <- wc[as.character(`Station id`) %in% names(TR_STN3),
          .(station = TR_STN3[as.character(`Station id`)], d,
            Klorid = .num(`Cl\nmg/L`), Natrium = .num(`Na\nmg/L`),
            Konduktivitet = .num(`Kond.\nmS/m`), Magnesium = .num(`Mg\nmg/L`),
            Kalsium = .num(`Ca\nmg/L`))]
tr3 <- melt(tr3, id.vars = c("station", "d"), variable.name = "pn", value.name = "value")[is.finite(value)]
tr3[, `:=`(tracer = TRACERS[as.character(pn)], doy = doy_of(d), src = "Trial (weekly)")]
tr3 <- tr3[, .(station, tracer, doy, d, value, src)]
cur3 <- rbind(cur3, tr3)[order(station, tracer, d, src != "Trial (weekly)")]
cur3 <- unique(cur3, by = c("station", "tracer", "d"))

clim_pred3 <- rbindlist(lapply(split(h3[era == "clim"], by = c("station", "tracer"), drop = TRUE), function(df) {
  if (nrow(df) < 25) return(NULL)
  m <- mgcv::gam(value ~ s(doy, bs = "cc", k = 8), data = df,
                 knots = list(doy = c(0.5, 366.5)), method = "REML")
  pr <- predict(m, newdata = grid, se.fit = TRUE)
  data.table(station = df$station[1], tracer = df$tracer[1], doy = grid$doy,
             fit = as.numeric(pr$fit), lo = pr$fit - 2 * pr$se.fit, hi = pr$fit + 2 * pr$se.fit)
}))
for (dt in list(clim_pred3, cur3)) {
  dt[, tracer  := factor(tracer,  levels = TRACERS)]
  dt[, station := factor(station, levels = STN3)]
}
cur3[, src := factor(src, levels = c("Norwegian (monthly)", "Trial (weekly)"))]

# labels / anomalies
NY_LAB <- STN3[["81555"]]; LA_LAB <- STN3[["45779"]]; KL_LAB <- STN3[["58890"]]
SHORT3 <- setNames(c("NY (upstream)", "LA (Litleåne)", "KL (mouth)"),
                   c(NY_LAB, LA_LAB, KL_LAB))
LITUP_GREY <- "grey60"
LITUP_LAB  <- "Litleåne (upstream of doser)"
band3 <- data.frame(x1 = TRIAL0, x2 = TRIAL1)

anom3 <- merge(cur3, clim_pred3[, .(station, tracer, doy, fit, lo, hi)],
               by = c("station", "tracer", "doy"), all.x = TRUE)
anom3[, `:=`(anom = value - fit, pct = (value / fit - 1) * 100)]
env3 <- clim_pred3[, .(ylo = -2 * mean((hi - lo) / 4), yhi = 2 * mean((hi - lo) / 4)),
                   by = tracer]

theme_left <- function(md_y = FALSE)
  theme_kvina_md(base_size = 12) +   # supplement scale: axis title 12, axis text 11, legend 11
  theme(axis.title.y = element_markdown(size = 12, face = "bold"),
        plot.tag = element_text(size = 16, face = "bold", colour = "grey5"),
        legend.position = "right", legend.justification = "center",
        legend.text = element_markdown(size = 11),
        axis.text.x = element_text(size = 11), axis.text.y = element_text(size = 11),
        plot.margin = margin(4, 8, 4, 4))

# (a) % above normal (Cl/Na/Mg mean), NY/LA/KL
idx3 <- anom3[tracer %in% as.character(MARINE3) & is.finite(pct),
              .(pct = mean(pct), src = src[1]), by = .(station, d, doy)]
setorder(idx3, station, doy)
p_a <- ggplot(idx3, aes(doy, pct, colour = station)) +
  geom_rect(data = band3, aes(xmin = x1, xmax = x2, ymin = -Inf, ymax = Inf),
            fill = TRIAL_FILL, alpha = 0.5, inherit.aes = FALSE) +
  geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.4) +
  geom_smooth(aes(group = station), method = "loess", span = 0.65, se = FALSE, linewidth = 1) +
  geom_point(aes(shape = src), size = 1.9) +
  scale_colour_manual(values = stncol3, labels = SHORT3, name = "Site") +
  scale_shape_manual(values = SHP, labels = SRC_LAB, name = "Sampling") +
  scale_x_continuous(breaks = mon_breaks, labels = month.abb[QMONTHS],
                     limits = c(1, 366), expand = expansion(mult = c(0.015, 0.015))) +
  guides(colour = guide_legend(order = 1, override.aes = list(shape = 16, linewidth = 1)),
         shape  = guide_legend(order = 2, override.aes = list(colour = "grey25", linetype = 0))) +
  labs(x = NULL, y = "Cl, Na, Mg anomaly (%)", tag = "a") +
  theme_left()

# (b) Magnesium anomaly, NY/LA/KL
aMg   <- anom3[tracer == TRACERS[["Magnesium"]]]
envMg <- env3[tracer == TRACERS[["Magnesium"]]]
p_b <- ggplot() +
  geom_rect(data = band3, aes(xmin = x1, xmax = x2, ymin = -Inf, ymax = Inf),
            fill = TRIAL_FILL, alpha = 0.5, inherit.aes = FALSE) +
  geom_rect(data = envMg, aes(xmin = -Inf, xmax = Inf, ymin = ylo, ymax = yhi),
            fill = "grey55", alpha = 0.22, inherit.aes = FALSE) +
  geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.4) +
  geom_line(data = aMg, aes(doy, anom, colour = station, group = station), linewidth = 0.7) +
  geom_point(data = aMg, aes(doy, anom, colour = station, shape = src), size = 1.7) +
  scale_colour_manual(values = stncol3, labels = SHORT3, name = "Site") +
  scale_shape_manual(values = SHP, labels = SRC_LAB, name = "Sampling") +
  scale_x_continuous(breaks = mon_breaks, labels = month.abb[QMONTHS],
                     limits = c(1, 366), expand = expansion(mult = c(0.015, 0.015))) +
  scale_y_continuous(n.breaks = 5, expand = expansion(mult = c(0.05, 0.07))) +
  guides(colour = guide_legend(order = 1, override.aes = list(shape = 16, linewidth = 1)),
         shape  = guide_legend(order = 2, override.aes = list(colour = "grey25", linetype = 0))) +
  labs(x = NULL, y = "Magnesium anomaly (mg L<sup>-1</sup>)", tag = "b") +
  theme_left(md_y = TRUE)

# (c) the Cl/Na/Mg/Ca record as ONE FACETED plot
tl4 <- as.character(MARINE3_CA)
cd4 <- cur3[tracer %in% tl4];       cd4[, tracer := factor(tracer, levels = tl4)]
cp4 <- clim_pred3[tracer %in% tl4]; cp4[, tracer := factor(tracer, levels = tl4)]
ca_lab  <- TRACERS[["Kalsium"]]
laca_cp <- cp4[station == LA_LAB & tracer == ca_lab]     # LA Ca historical baseline (45779) -> grey
cp4m    <- cp4[!(station == LA_LAB & tracer == ca_lab)]  # all other historical GAMs (colour-mapped)
p_c_facet <- ggplot() +
  geom_rect(data = data.frame(x1 = TRIAL0, x2 = TRIAL1),
            aes(xmin = x1, xmax = x2, ymin = -Inf, ymax = Inf),
            fill = TRIAL_FILL, alpha = 0.55, inherit.aes = FALSE) +
  geom_ribbon(data = cp4m, aes(doy, ymin = lo, ymax = hi, group = station, fill = station), alpha = 0.16) +
  geom_line(data = cp4m, aes(doy, fit, colour = station, group = station),
            linetype = "dashed", linewidth = 0.6) +
  geom_ribbon(data = laca_cp, aes(doy, ymin = lo, ymax = hi),
              fill = LITUP_GREY, alpha = 0.18, inherit.aes = FALSE) +
  geom_line(data = laca_cp, aes(doy, fit), colour = LITUP_GREY,
            linetype = "dashed", linewidth = 0.6, inherit.aes = FALSE) +
  geom_smooth(data = cd4, aes(doy, value, colour = station, group = station),
              method = "loess", span = 0.6, se = FALSE, linewidth = 0.9) +
  facet_wrap(~ tracer, ncol = 1, scales = "free_y", strip.position = "top") +
  scale_colour_manual(values = stncol3, labels = SHORT3, name = "Site") +
  scale_fill_manual(values = stncol3, guide = "none") +
  scale_x_continuous(breaks = mon_breaks, labels = month.abb[QMONTHS],
                     limits = c(1, 366), expand = expansion(mult = c(0.015, 0.015))) +
  scale_y_continuous(n.breaks = 5, expand = expansion(mult = c(0.04, 0.06))) +
  labs(x = NULL, y = NULL) +
  theme_kvina_md(base_size = 12) +
  # supplement scale: strips 12 pt, axis text 11 (consistent with S1-S6, panels a/b):
  theme(strip.text = element_markdown(size = 12, face = "bold", colour = "grey10"),
        panel.spacing = grid::unit(2.5, "mm"),
        axis.text.x = element_text(size = 11), axis.text.y = element_text(size = 11),
        plot.tag = element_text(size = 16, face = "bold", colour = "grey5"),
        plot.margin = margin(4, 6, 4, 4))

# shared legend machinery (deterministic single legend)
no_leg <- theme(legend.position = "none")
extract_legend <- function(p, pos = "bottom") {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  g <- ggplotGrob(p)
  i <- which(g$layout$name == paste0("guide-box-", pos))
  if (!length(i)) i <- which(grepl("guide-box", g$layout$name))[1]
  g$grobs[[i]]
}
make_legend <- function(keys = c("site", "sampling", "gam", "trial"), litup = FALSE) {
  st_lev <- c(NY_LAB, LA_LAB, KL_LAB); st_val <- stncol3; st_lab <- SHORT3
  if (litup) {
    st_lev <- c(st_lev, "LITUP")
    st_val <- c(stncol3, LITUP = LITUP_GREY)
    st_lab <- c(SHORT3,  LITUP = LITUP_LAB)
  }
  ld <- data.frame(x = 1, y = 1, station = factor(st_lev, levels = st_lev))
  d  <- ggplot(ld, aes(x, y)) +
    scale_colour_manual(values = st_val, labels = st_lab, name = "Site")
  gl <- list(colour = if ("site" %in% keys)
    guide_legend(order = 1, override.aes = list(shape = 16, linewidth = 1)) else "none")
  if ("site" %in% keys)
    d <- d + geom_point(aes(colour = station), size = 1.9)
  if ("trial" %in% keys) {
    d <- d + geom_rect(data = data.frame(b = TRIAL_LAB),
                       aes(xmin = 0, xmax = 0, ymin = 0, ymax = 0, fill = b), inherit.aes = FALSE) +
      scale_fill_manual(values = setNames(TRIAL_FILL, TRIAL_LAB),
                        labels = setNames("Trial Period", TRIAL_LAB), name = NULL)
    gl$fill <- guide_legend(order = 4, override.aes = list(alpha = 0.55))
  }
  if ("sampling" %in% keys) {
    ls <- data.frame(x = 1, y = 1,
                     src = factor(c("Norwegian (monthly)", "Trial (weekly)"),
                                  levels = c("Norwegian (monthly)", "Trial (weekly)")))
    d <- d + geom_point(data = ls, aes(x, y, shape = src), size = 1.9, inherit.aes = FALSE) +
      scale_shape_manual(values = SHP, labels = SRC_LAB, name = "Sampling")
    gl$shape <- guide_legend(order = 2, override.aes = list(colour = "grey25", linetype = 0))
  }
  if ("gam" %in% keys) {
    d <- d + geom_line(aes(colour = station, linetype = "Historical GAM", group = station)) +
      scale_linetype_manual(values = c("Historical GAM" = "dashed"), name = NULL)
    gl$linetype <- guide_legend(order = 3, override.aes = list(colour = "grey35"))
  }
  d <- d + do.call(guides, gl) +
    theme_kvina_md(base_size = 12) +
    theme(legend.position = "right", legend.direction = "vertical", legend.box = "vertical",
          legend.title = element_markdown(size = 12, face = "bold"),
          legend.text = element_markdown(size = 11), legend.key.size = grid::unit(1.1, "lines"))
  extract_legend(d, "right")
}
leg <- make_legend(c("site", "sampling", "gam", "trial"), litup = TRUE)

# combined composite (portrait, right-side legend for consistency with S1-S6)
design <- paste(c("AC", "BC"), collapse = "\n")
p_body <- (p_b + labs(tag = "a") + no_leg) +
  free(p_a + labs(tag = "b") + no_leg, type = "space", side = "b") +
  (p_c_facet + labs(tag = "c") + no_leg) +
  plot_layout(design = design, widths = c(1.08, 1))
p_comb <- (p_body | wrap_elements(full = leg)) +
  plot_layout(widths = c(1, 0.34))
ggsave(file.path(fig_dir, "S7_salt_pulse.png"),
       p_comb, width = 9.8, height = 10.6, dpi = 450, bg = "white")
cat(sprintf("wrote S7_salt_pulse.png; anom3 rows: %d\n", nrow(anom3)))

cat("done.\n")

# GRAPHICAL ABSTRACT / TABLE-OF-CONTENTS FIGURE (TOC.png) --------------------
# Pairs the Kvitla downstream-pH record (historical years grey, 2025 trial in
# KL blue) with a compact five-stage CDR waterfall. Independent of the chemistry
# tables above: reads the two doser workbooks for hourly "pH Kvitla" telemetry.
# The italic CDR statement is baked into the figure deliberately.

# ES&T Table-of-Contents graphic: max 3.25 x 1.75 in, >=300 dpi. The physical
# size is capped small (3.25 in), so to match the pixel resolution of the ~4230-px
# supplement figures the dpi is raised: 3.25 x 1.625 in @ 1300 dpi = 4225 x 2112 px.
# Physical size still fits the ES&T cap on both axes and 1300 dpi is far above the
# 300 dpi floor. Type/line sizes were tuned on a 7.5-in reference and are
# multiplied by S so the figure keeps those proportions at TOC size.
PLOT_WIDTH_IN  <- 3.25
PLOT_HEIGHT_IN <- 3.25 / 2
PLOT_DPI       <- 1300
S              <- PLOT_WIDTH_IN / 7.5   # type/line scale from the 7.5-in reference

TOC_DOSER_XLSX <- data_path("Nyland, Kvina, 30.09.2025 00 - 15.05.2025 00 COMPLETE PROJECT DATA.xlsx")
HIST_XLSX  <- data_path("Nyland dosererdata.xlsx")

REF_YEAR_START   <- as.Date("2001-01-01")
REF_YEAR_END     <- as.Date("2002-01-01")
REF_TRIAL_START  <- as.Date("2001-06-10")
REF_TRIAL_END    <- as.Date("2001-09-15")
REF_AXIS_BREAKS  <- as.Date(c("2001-01-01", "2001-07-01", "2002-01-01"))

# Kvitla pH record: weekly traces for 2015-2025 -------------------------------

clean_ph_sensor <- function(x, lo = 4.5, hi = 8.5, stuck_run = 6,
                            dom_frac = 0.4, roll_k = 11, mad_k = 6) {
  x <- as.numeric(x)
  x[x < lo | x > hi] <- NA_real_
  r <- rle(x)
  long <- !is.na(r$values) & r$lengths >= stuck_run
  if (any(long)) {
    tab <- table(r$values[long])
    sentinels <- as.numeric(names(tab)[tab / sum(long) >= dom_frac])
    if (length(sentinels)) x[x %in% sentinels] <- NA_real_
  }
  n <- length(x)
  half <- roll_k %/% 2
  med <- madv <- rep(NA_real_, n)
  for (i in seq_len(n)) {
    if (is.na(x[i])) next
    w <- x[max(1, i - half):min(n, i + half)]
    w <- w[!is.na(w)]
    if (length(w) >= 3) {
      med[i]  <- median(w)
      madv[i] <- mad(w)
    }
  }
  is_spike <- !is.na(med) & !is.na(madv) & madv > 0 & abs(x - med) > mad_k * madv
  x[is_spike] <- NA_real_
  x
}

read_ph_signals <- function(path, upstream_col, kvitla_col) {
  sig <- as.data.table(read_excel(path, sheet = "Signals", skip = 2))
  missing_cols <- setdiff(c("Date", upstream_col, kvitla_col), names(sig))
  if (length(missing_cols))
    stop("Missing required column(s) in ", basename(path), ": ",
         paste(missing_cols, collapse = ", "), call. = FALSE)
  sig[, ts := as.POSIXct(Date, format = "%d.%m.%Y %H", tz = "UTC")]
  if (all(is.na(sig$ts)))
    stop("Could not parse any timestamps in ", basename(path), call. = FALSE)
  setnames(sig, c(upstream_col, kvitla_col), c("upstream", "kvitla"))
  sig[!is.na(ts)][order(ts), .(ts, kvitla)]
}

trial_sig <- read_ph_signals(TOC_DOSER_XLSX, "pH oppstrøm", "pH Kvitla")
hist_sig  <- read_ph_signals(HIST_XLSX,  "pH upstream", "pH Kvitla")[ts < min(trial_sig$ts)]

sig <- rbindlist(list(hist_sig, trial_sig), use.names = TRUE)[order(ts)]
sig[, yr := year(ts)]
sig[, ph_clean := clean_ph_sensor(kvitla), by = yr]
sig[, week := week(ts)]

ph_years <- sig[yr <= 2025, .(pH = mean(ph_clean, na.rm = TRUE)), by = .(yr, week)]
ph_years <- ph_years[is.finite(pH)]
if (!nrow(ph_years) || !2025 %in% ph_years$yr)
  stop("Cleaned pH record is empty or does not contain 2025.", call. = FALSE)
ph_years[, day := REF_YEAR_START + (week - 1) * 7]

KL_BLUE      <- unname(SITE_COLS[["KL"]])
TRIAL_MARKER <- "#8CC8E8"
historical_years <- sort(unique(ph_years[yr < 2025, yr]))
if (!length(historical_years))
  stop("No pre-2025 pH record remains after cleaning.", call. = FALSE)
ph_years[, yr := factor(yr, levels = c(historical_years, 2025))]
YEAR_COLS <- c(
  setNames(gray.colors(length(historical_years), start = 0.78, end = 0.25),
           as.character(historical_years)),
  "2025" = KL_BLUE)

make_ph_panel <- function(trial_onset = REF_TRIAL_START, trial_end = REF_TRIAL_END) {
  p <- ggplot(ph_years, aes(day, pH, group = yr, colour = yr)) +
    annotate("rect", xmin = trial_onset, xmax = trial_end,
             ymin = -Inf, ymax = Inf, fill = TRIAL_MARKER, alpha = 0.16) +
    geom_line(data = ph_years[yr != "2025"], linewidth = 0.35 * S) +
    geom_line(data = ph_years[yr == "2025"], linewidth = 0.9 * S) +
    annotate("text", x = as.Date("2001-02-04"), y = 5.72,
             label = "Historical data", hjust = 0, size = 3.4 * S,
             colour = "grey35", fontface = "italic") +
    annotate("label", x = trial_onset + (trial_end - trial_onset) / 2,
             y = 7.18, label = "Trial Period", size = 3.4 * S,
             colour = "grey10", fontface = "bold",
             fill = scales::alpha("white", 0.58), linewidth = 0,
             label.padding = grid::unit(0.04, "lines")) +
    scale_colour_manual(values = YEAR_COLS, guide = "none") +
    scale_x_date(breaks = REF_AXIS_BREAKS, labels = c("Jan", "Jul", "Jan")) +
    scale_y_continuous(breaks = c(6.0, 6.4, 6.8)) +
    coord_cartesian(xlim = c(REF_YEAR_START, REF_YEAR_END), ylim = c(5.7, 7.2)) +
    labs(x = NULL, y = "Downstream pH") +
    theme_bw(base_size = 9 * S) +
    theme(
      panel.grid = element_line(colour = "grey91", linewidth = 0.25 * S),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(colour = "grey35", fill = NA, linewidth = 0.35 * S),
      axis.text = element_text(size = 8 * S, colour = "grey30"),
      axis.ticks = element_line(colour = "grey35", linewidth = 0.25 * S),
      axis.title.y = element_text(size = 12 * S, colour = "grey15", margin = margin(r = 4 * S)),
      aspect.ratio = 1,
      plot.margin = margin(1, 1, 1, 1))
}

# Five-stage CDR waterfall ----------------------------------------------------

cdr <- list(
  delta_alk = 1020, delta_alk_se = 40,
  feedstock_carbon = 629,
  efficiency_losses = 143 + 100.1,
  lca = 35.9,
  net_cdr = 110, net_cdr_se = 30)

after_feedstock <- cdr$delta_alk - cdr$feedstock_carbon
after_losses    <- cdr$net_cdr + cdr$lca

# Reported rounded values do not close perfectly: geometry uses closed
# accounting levels while labels retain reported values; the rounding residual
# sits in the pooled efficiency-loss step.
fmt_whole <- function(x) format(round(x), trim = TRUE, scientific = FALSE)

wf <- data.frame(
  stage = c("Trial\nΔ Alkalinity", "Feedstock\nCarbon",
            "Efficiency\nlosses", "LCA", "Net CDR"),
  type  = c("subtotal", "deduction", "deduction", "deduction", "final"),
  y0    = c(0, after_feedstock, after_losses, cdr$net_cdr, 0),
  y1    = c(cdr$delta_alk, cdr$delta_alk, after_feedstock, after_losses, cdr$net_cdr),
  label = c(fmt_whole(cdr$delta_alk),
            paste0("−", fmt_whole(cdr$feedstock_carbon)),
            paste0("−", fmt_whole(cdr$efficiency_losses)),
            paste0("−", fmt_whole(cdr$lca)),
            fmt_whole(cdr$net_cdr)),
  se    = c(cdr$delta_alk_se, NA, NA, NA, cdr$net_cdr_se),
  stringsAsFactors = FALSE)
wf$x       <- seq_len(nrow(wf))
wf$ymin    <- pmin(wf$y0, wf$y1)
wf$ymax    <- pmax(wf$y0, wf$y1)
wf$reached <- ifelse(wf$type == "deduction", wf$ymin, wf$y1)
wf$val_y   <- wf$ymax + ifelse(is.na(wf$se), 0, wf$se) + 36
wf$stage_y <- wf$val_y + 96
# the first two bars share a top accounting level; align their two-line headings
wf$stage_y[1:2] <- max(wf$stage_y[1:2])

connectors <- data.frame(x = head(wf$x, -1) + 0.36,
                         xend = tail(wf$x, -1) - 0.36,
                         y = head(wf$reached, -1))
errorbars <- with(wf[!is.na(wf$se), ], data.frame(x = x, y = ymax, se = se))

make_waterfall_panel <- function() {
  wf_plot <- wf
  
  p <- ggplot(wf_plot, aes(x = x)) +
    geom_rect(aes(xmin = x - 0.36, xmax = x + 0.36, ymin = ymin, ymax = ymax,
                  fill = type), colour = "grey30", linewidth = 0.25 * S) +
    geom_segment(data = connectors, inherit.aes = FALSE,
                 aes(x = x, xend = xend, y = y, yend = y),
                 colour = "grey55", linetype = "22", linewidth = 0.3 * S) +
    geom_errorbar(data = errorbars, inherit.aes = FALSE,
                  aes(x = x, ymin = y - se, ymax = y + se),
                  width = 0.11, linewidth = 0.35 * S, colour = "grey10")
  
  p <- p +
    geom_text(data = subset(wf_plot, type != "final"),
              aes(y = val_y, label = label), size = 3.5 * S, vjust = 0) +
    geom_text(data = subset(wf_plot, type != "final"),
              aes(y = stage_y, label = stage), size = 3.0 * S,
              fontface = "italic", colour = "grey30", vjust = 0) +
    geom_text(data = subset(wf_plot, type == "final"),
              aes(y = val_y, label = label), size = 3.5 * S,
              fontface = "bold", colour = "grey10", vjust = 0) +
    geom_text(data = subset(wf_plot, type == "final"),
              aes(y = stage_y, label = stage), size = 3.4 * S,
              fontface = "bold", colour = "grey10", vjust = 0) +
    annotate("text", x = 4.05, y = 975,
             label = "Carbon Dioxide Removal\nattributable to added\nalkalinity",
             size = 3.1 * S, fontface = "italic", colour = "grey25",
             lineheight = 0.95)
  
  p +
    scale_fill_manual(values = c(subtotal = "grey45", deduction = "grey85",
                                 final = KL_BLUE), guide = "none") +
    scale_x_continuous(expand = expansion(mult = 0.035)) +
    scale_y_continuous(breaks = c(0, 500, 1000),
                       expand = expansion(mult = c(0, 0.31))) +
    labs(x = NULL, y = expression("t CO"[2] * "e")) +
    coord_cartesian(clip = "off") +
    theme_bw(base_size = 9 * S) +
    theme(
      panel.grid = element_blank(), panel.border = element_blank(),
      axis.line = element_line(colour = "grey30", linewidth = 0.35 * S),
      axis.ticks.x = element_blank(), axis.text.x = element_blank(),
      axis.title.y = element_text(size = 12 * S, colour = "grey10", margin = margin(r = 3 * S)),
      axis.text.y = element_text(size = 8.5 * S, colour = "grey25"),
      aspect.ratio = 1,
      plot.margin = margin(1, 2, 1, 2))
}

arrow_panel <- ggplot() +
  geom_segment(aes(x = 0.10, xend = 0.90, y = 0.5, yend = 0.5),
               colour = "black", linewidth = 0.65 * S,
               arrow = arrow(angle = 24, length = grid::unit(2.2 * S, "mm"),
                             type = "closed")) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
  theme_void()

toc <- make_ph_panel() + arrow_panel + make_waterfall_panel() +
  plot_layout(widths = c(0.94, 0.14, 1.06))

out <- file.path(fig_manuscript, "toc.png")
save_args <- list(filename = out, plot = toc, width = PLOT_WIDTH_IN,
                  height = PLOT_HEIGHT_IN, dpi = PLOT_DPI, bg = "white")
if (requireNamespace("ragg", quietly = TRUE)) {
  save_args$device <- ragg::agg_png
} else {
  save_args$type <- "cairo"
}
do.call(ggsave, save_args)
message("wrote figures/", basename(out))
