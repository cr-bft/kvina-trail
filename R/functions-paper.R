# R/functions-paper.R -- shared helpers: alkalinity correction, plotting theme,
# split-violin geom, stats helpers, calcite saturation index. Methods are
# described in the manuscript. Sourced with working directory = repo root.

suppressPackageStartupMessages({
  library(readxl)
  library(data.table)
  library(ggplot2)
  library(ggtext)
  library(dplyr)
  library(phreeqc)
})

# data layout ---------------------------------------------------------------
# This repo is self-contained: every input it reads lives under data/. Files
# shared with kvina-trail keep that repo's filenames, so the
# same scripts also run unchanged inside it. Search order:
#   1. KVINA_DATA_DIRS  -- ":"-separated override
#   2. data-raw/        -- when these scripts are run inside the public repo
#   3. this repo's own data/raw/ and data/metadata/

data_dirs <- function() {
  env <- Sys.getenv("KVINA_DATA_DIRS", "")
  c(if (nzchar(env)) strsplit(env, ":", fixed = TRUE)[[1]],
    "data-raw", file.path("data", "raw"), file.path("data", "metadata"))
}

data_path <- function(file) {
  dirs <- data_dirs()
  hit  <- Filter(file.exists, file.path(dirs, file))
  if (!length(hit))
    stop(sprintf("input not found in %s: %s",
                 paste(dirs, collapse = " / "), file), call. = FALSE)
  hit[1]
}

# the 2025 trial workbook ships under two names; the "WaterChem Sorted" sheet is
# identical in both, and no other sheet is read:

TRIAL_XLSX_NAMES <- c("Data report Kvina complete.xlsx", "Kvina_RAE_Trial_Period.xlsx")

trial_xlsx_path <- function() {
  dirs <- data_dirs()
  hit  <- Filter(file.exists, as.vector(t(outer(dirs, TRIAL_XLSX_NAMES, file.path))))
  if (!length(hit))
    stop(sprintf("trial workbook not found; looked for %s under %s",
                 paste(TRIAL_XLSX_NAMES, collapse = " / "),
                 paste(dirs, collapse = " / ")), call. = FALSE)
  hit[1]
}

# the "WaterChem Sorted" sheet as a plain data frame: header row, then data rows
# only (the sheet carries an uncertainty row plus a trailing caption block):

read_waterchem_sorted <- function() {
  raw <- as.data.frame(read_excel(trial_xlsx_path(), sheet = "WaterChem Sorted"),
                       check.names = FALSE)
  # the two workbooks break header lines differently ("pH\n" vs "pH\r\n"); the
  # column references below and in 01 assume "\n":
  names(raw) <- gsub("\r\n", "\n", names(raw), fixed = TRUE)
  raw <- raw[!is.na(raw[["Station code"]]), , drop = FALSE]
  # A column holding a censored "< x" comes back as character, so its numeric
  # cells keep whatever text the workbook stored ("37.0", "1.1000000000000001",
  # "5.7E-3"). Re-render those to canonical form, then apply read.csv()'s type
  # inference, so the processed tables do not depend on workbook serialisation.
  chr <- vapply(raw, is.character, logical(1))
  raw[chr] <- lapply(raw[chr], function(x) {
    n <- suppressWarnings(as.numeric(x))
    x[!is.na(n)] <- as.character(n[!is.na(n)])
    x
  })
  raw[chr] <- type.convert(raw[chr], as.is = TRUE)
  raw
}

# the Vannmiljo monitoring record (5 stations, 2011-2025). The public 2011-2024
# export supplies everything through 2024-12-02; vannmiljo_export_2025.csv adds
# the 2025 rows. Each measurement is one row, so read_vannmiljo_long() is a
# column selection, not a reshape.

# four VK3 (58892) measurements Vannmiljo retracted between the 2024 and 2025
# exports; the public file predates the retraction, so drop them:
VANNMILJO_RETRACTED <- data.table(
  Vannlokalitet_kode = "025-58892",
  Parameter_navn = c("Konduktivitet", "Konduktivitet", "Konduktivitet", "Kalsium"),
  date10 = c("2012-01-03", "2015-02-02", "2015-03-02", "2023-09-04"))

read_vannmiljo <- function() {
  wide <- rbind(
    fread(data_path("Kvina historical_data_2011_2025.csv"), showProgress = FALSE),
    fread(data_path("vannmiljo_export_2025.csv"), showProgress = FALSE),
    fill = TRUE)
  wide[, date10 := substr(Tid_provetak, 1, 10)]
  wide <- wide[!VANNMILJO_RETRACTED, on = c("Vannlokalitet_kode", "Parameter_navn", "date10")]
  wide[, date10 := NULL]
  as.data.frame(wide)
}

read_vannmiljo_long <- function() {
  w <- read_vannmiljo()
  data.table(
    water_location_id = sub("^025-", "", w$Vannlokalitet_kode),
    parameter         = w$Parameter_navn,
    timestamp         = substr(w$Tid_provetak, 1, 10),
    value             = parse_value(w$Verdi))
}

# molar masses (g/mol):

MW_CA <- 40.078
MW_MG <- 24.305

# per-site colour palette (Okabe-Ito):

SITE_COLS <- c(NY = "#E69F00", OB = "#56B4E9", LA = "#009E73", KL = "#0072B2")

# analysis windows (trial and 2011-2024 KL historical baseline; month-day 610-915):

TRIAL_WINDOW <- list(md_min = 610L, md_max = 915L, year_min = 2025L, year_max = 2025L,
                     locality_match = "KL", label = "Trial 2025 (Jun 10-Sep 15)",
                     start = as.Date("2025-06-10"), end = as.Date("2025-09-15"))
HIST_WINDOW  <- list(md_min = 610L, md_max = 915L, year_min = 2011L, year_max = 2024L,
                     locality_match = "Kvina ved Klosterøyna",
                     label = "Site KL Jun 10-Sep 15 (2011-2024)")

# DOC / organic-acid correction to alkalinity (see manuscript):

DOC_CORR_DEFAULT_PH <- 6.3
ALK_ENDPOINT_PH     <- 4.5
CARB_TEMP_C         <- 10

# LA / St. 8 grab on 2025-06-10 is unfiltered in the filtered-alkalinity column:

ALK_UNFILTERED_STATION <- 74913L
ALK_UNFILTERED_DATE    <- as.Date("2025-06-10")

# hruska et al. (2003) triprotic organic-acid model:

HRUSKA_SD  <- 10.2
HRUSKA_PKA <- c(3.04, 4.51, 6.46)

# mean organic-acid charge per molecule at a given pH:

hruska_zbar <- function(pH) {
  H  <- 10^(-pH)
  K1 <- 10^(-HRUSKA_PKA[1]); K2 <- 10^(-HRUSKA_PKA[2]); K3 <- 10^(-HRUSKA_PKA[3])
  (K1 * H^2 + 2 * K1 * K2 * H + 3 * K1 * K2 * K3) /
    (H^3 + K1 * H^2 + K1 * K2 * H + K1 * K2 * K3)
}

# organic charge titrated between the sample pH and the endpoint:

hruska_organic_titrated_ueq_L <- function(doc_mg_L, pH, ph_e = ALK_ENDPOINT_PH)
  doc_mg_L * (HRUSKA_SD / 3) * (hruska_zbar(pH) - hruska_zbar(ph_e))

# water dissociation at CARB_TEMP_C:

CARB_PKW <- 4471 / (CARB_TEMP_C + 273.15) + 0.01706 * (CARB_TEMP_C + 273.15) - 6.0875
oh_ueq_L <- function(pH) 10^(pH - CARB_PKW) * 1e6

# carbonate alkalinity from a raw ISO 9963-1 titration:

carb_alk_from_titration_ueq_L <- function(alk_ueq_L, pH, doc_mg_L,
                                          ph_e = ALK_ENDPOINT_PH) {
  H0  <- 10^(-pH) * 1e6
  OH0 <- oh_ueq_L(pH)
  He  <- 10^(-ph_e) * 1e6
  A_titr <- hruska_organic_titrated_ueq_L(doc_mg_L, pH, ph_e)
  alk_ueq_L - He + H0 - OH0 - A_titr
}

# comma-or-dot decimal safe numeric coercion:

.num <- function(x) suppressWarnings(as.numeric(gsub(",", ".", as.character(x))))

# significance stars from a p-value:

star <- function(p) data.table::fifelse(is.na(p), "", data.table::fifelse(p < 0.001, "***",
                                                                          data.table::fifelse(p < 0.01, "**", data.table::fifelse(p < 0.05, "*", "ns"))))

# TRUE where a date's month-day falls in the trial window:

in_trial_date_window <- function(dates, start = trial_start, end = trial_end) {
  md <- format(as.Date(dates), "%m-%d")
  start_md <- format(start, "%m-%d")
  end_md <- format(end, "%m-%d")
  if (start_md <= end_md) md >= start_md & md <= end_md else md >= start_md | md <= end_md
}

# 2025 trial stations (NY/OB/LA/KL), filtered and total, for flux work:

load_trial_stations <- function(stations = c(NY = 74911L, OB = 74912L,
                                             LA = 74913L, KL = 74914L)) {
  wc <- as.data.table(read_waterchem_sorted())
  wc[, ts := as.Date(`Sample date`)]
  out <- wc[, .(
    station_id = as.integer(`Station id`),
    ts,
    ph        = .num(`pH\n`),
    anc_ueq_L = .num(`ANC\nµEkv/L`),
    alk_meas_mmol_L = .num(`Alk-filt\nmmol/l`),
    ca_mg_L   = .num(`Ca-filt\nmg/l`),
    mg_mg_L   = .num(`Mg-filt\nmg/l`),
    na_mg_L   = .num(`Na-filt\nmg/l`),
    k_mg_L    = .num(`K-filt\nmg/l`),
    ca_tot_mg_L = .num(`Ca\nmg/L`),
    mg_tot_mg_L = .num(`Mg\nmg/L`),
    na_tot_mg_L = .num(`Na\nmg/L`),
    k_tot_mg_L  = .num(`K\nmg/L`),
    toc_mg_L  = .num(`TOC\nmg C/l`),
    doc_mg_L  = .num(`DOC\nmg C/L`),
    dic_mg_L  = .num(`DIC\nmg C/L`),
    so4_mg_L  = .num(`SO4\nmg/L`),
    no3n_ugN_L = .num(`NO3-N\nµg/l`),
    cond_mS_m = .num(`Kond.\nmS/m`)
  )]
  out <- out[station_id %in% stations]
  out[, alk_unfiltered := station_id == ALK_UNFILTERED_STATION &
        ts == ALK_UNFILTERED_DATE]
  out[, station := factor(names(stations)[match(station_id, stations)],
                          levels = names(stations))]
  out[, carb_alk_meq_L := anc_ueq_L / 1000]
  out[, camg_meq_L := 2 * (ca_mg_L / MW_CA + mg_mg_L / MW_MG)]
  out[order(station, ts)]
}


# ggplot theme (theme_bw base + ggtext markdown labels):

theme_kvina_md <- function(base_size = 10) {
  theme_bw(base_size = base_size) +
    theme(
      plot.background  = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA),
      panel.grid.major = element_line(colour = "grey94", linewidth = 0.20),
      panel.grid.minor = element_blank(),
      panel.border     = element_rect(colour = "grey30", fill = NA, linewidth = 0.55),
      panel.spacing    = grid::unit(8, "mm"),
      axis.title.x = element_markdown(face = "bold", colour = "grey10", size = base_size),
      axis.title.y = element_markdown(face = "bold", colour = "grey10", size = base_size),
      axis.text    = element_text(colour = "grey25", size = base_size - 1),
      plot.title   = element_markdown(face = "bold", colour = "grey5", size = base_size + 1),
      plot.margin  = margin(8, 10, 8, 8),
      legend.title = element_markdown(face = "bold", size = base_size - 1),
      legend.text  = element_markdown(size = base_size - 1),
      strip.text       = element_markdown(face = "bold", colour = "grey10",
                                          margin = margin(t = 4, b = 4)),
      strip.background = element_rect(fill = "grey92", colour = "grey30", linewidth = 0.5)
    )
}

# split-violin geom (each group = one half; optional median line per half):

GeomSplitViolin <- ggproto("GeomSplitViolin", GeomViolin,
                           draw_group = function(self, data, ..., med_quantiles = 0.5,
                                                 med_colour = "grey20", med_linewidth = 0.5) {
                             data <- transform(data, xminv = x - violinwidth * (x - xmin),
                                               xmaxv = x + violinwidth * (xmax - x))
                             grp     <- data[1, "group"]
                             centre  <- data[1, "x"]
                             left    <- grp %% 2 == 1
                             poly    <- if (left) transform(data, x = xminv)[order(data$y), ]
                             else      transform(data, x = xmaxv)[order(-data$y), ]
                             n <- nrow(poly)
                             newdata <- rbind(poly[1, ], poly, poly[n, ], poly[1, ])
                             newdata[c(1, n + 2, n + 3), "x"] <- round(centre)
                             poly_grob <- GeomPolygon$draw_panel(newdata, ...)
                             
                             if (length(med_quantiles) > 0 && !isTRUE(all.equal(min(data$y), max(data$y)))) {
                               o    <- order(data$y)
                               yy   <- data$y[o]
                               cdf  <- cumsum(data$density[o]) / sum(data$density[o])
                               qy   <- stats::approx(cdf, yy, xout = med_quantiles, ties = "ordered")$y
                               edge <- if (left) stats::approx(yy, data$xminv[o], xout = qy)$y
                               else      stats::approx(yy, data$xmaxv[o], xout = qy)$y
                               seg  <- data[rep(1, 2 * length(qy)), , drop = FALSE]
                               seg$x         <- as.vector(rbind(rep(centre, length(qy)), edge))
                               seg$y         <- rep(qy, each = 2)
                               seg$group     <- rep(seq_along(qy), each = 2)
                               seg$colour    <- med_colour
                               seg$linewidth <- med_linewidth
                               seg$linetype  <- 1
                               seg$alpha     <- 1
                               seg  <- seg[stats::complete.cases(seg$x, seg$y), , drop = FALSE]
                               if (nrow(seg) >= 2)
                                 return(grid::grobTree(poly_grob, GeomPath$draw_panel(seg, ...)))
                             }
                             poly_grob
                           })
geom_split_violin <- function(mapping = NULL, data = NULL, stat = "ydensity",
                              position = "identity", ..., med_quantiles = 0.5,
                              med_colour = "grey20", med_linewidth = 0.5,
                              trim = TRUE, scale = "width",
                              na.rm = FALSE, show.legend = NA, inherit.aes = TRUE) {
  layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin,
        position = position, show.legend = show.legend, inherit.aes = inherit.aes,
        params = list(med_quantiles = med_quantiles, med_colour = med_colour,
                      med_linewidth = med_linewidth, trim = trim, scale = scale,
                      na.rm = na.rm, ...))
}

# stats helpers ------------------------------------------------------------

read_csv_raw <- function(...) read.csv(file.path(...), check.names = FALSE,
                                       fileEncoding = "UTF-8", stringsAsFactors = FALSE)

# restrict a data frame to the trial window by month-day:

clip_to_trial_window <- function(df, date_col, trial_win) {
  md <- as.integer(format(as.Date(df[[date_col]]), "%m%d"))
  df[md >= trial_win$md_min & md <= trial_win$md_max, , drop = FALSE]
}

# numeric value of a raw cell ("," -> "."; below-detection "< x" -> x):

parse_value <- function(x) {
  x <- trimws(as.character(x))
  x[x == ""] <- NA
  x <- gsub(",", ".", x, fixed = TRUE)
  x <- sub("^<\\s*", "", x)
  suppressWarnings(as.numeric(x))
}

# standardise a NIVA / Vannmiljo long-format export to the tidy schema:

niva_long <- function(df) {
  transmute(df,
            locality   = Vannlokalitetsnavn,
            date       = as.Date(substr(Tid_provetak, 1, 10)),
            year       = as.integer(format(date, "%Y")),
            month      = as.integer(format(date, "%m")),
            param_navn = Parameter_navn,
            fraction   = ifelse(Filtrert_Prove == "Filtrert", "filtered", "unfiltered"),
            operator   = Operator,
            unit       = Enhet,
            value      = parse_value(Verdi)
  )
}

# significance threshold:

alpha <- 0.05

# NA-safe summary statistics:

s_min <- function(x) if (length(x))     min(x)    else NA_real_
s_max <- function(x) if (length(x))     max(x)    else NA_real_
s_med <- function(x) if (length(x))     median(x) else NA_real_
s_mn  <- function(x) if (length(x))     mean(x)   else NA_real_
s_sd  <- function(x) if (length(x) > 1) sd(x)     else NA_real_

# display label composed from the analyte name plus its fraction:

disp_label <- function(label, fraction) {
  mapply(function(l, f) {
    if (is.na(f) || f == "") l
    else if (endsWith(l, ")")) sub("\\)$", paste0(", ", f, ")"), l)
    else sprintf("%s (%s)", l, f)
  }, label, fraction, USE.NAMES = FALSE)
}

# load and validate the analyte registry (data/metadata/analytes.csv):

load_analytes <- function() {
  parm <- read_csv_raw(data_path("analytes.csv"))
  flag_cols <- c("in_ny_ob", "in_hist_vs_trial", "in_routine_2025")
  parm[flag_cols] <- lapply(parm[flag_cols], as.logical)
  for (col in c("trial_column", "historical_names", "fraction", "fraction_match"))
    parm[[col]][is.na(parm[[col]])] <- ""
  parm$conv[is.na(parm$conv)] <- 1
  parm$hist_name_list <- strsplit(parm$historical_names, ";\\s*")
  parm$label <- disp_label(parm$display_label, parm$fraction)
  stopifnot(
    "analytes.csv: duplicate param_id" = !anyDuplicated(parm$param_id),
    "analytes.csv: duplicate ord" = !anyDuplicated(parm$ord),
    "analytes.csv: in_hist_vs_trial rows need a trial_column (except SI_calcite)" =
      all(nzchar(parm$trial_column[parm$in_hist_vs_trial & parm$param_id != "SI_calcite"])),
    "analytes.csv: in_hist_vs_trial / in_routine_2025 rows need historical_names (except SI_calcite)" =
      all(nzchar(parm$historical_names[(parm$in_hist_vs_trial | parm$in_routine_2025) &
                                         parm$param_id != "SI_calcite"])),
    "analytes.csv: CALCITE_PARAM_IDS must all be present" =
      all(CALCITE_PARAM_IDS %in% parm$param_id)
  )
  parm[order(parm$ord), ]
}

# Welch t-test (unpaired, two-sided); "Not tested" if n < min_n or zero variance:

welch_test <- function(x, y, min_n = 3) {
  out <- data.frame(t = NA_real_, df = NA_real_, p_value = NA_real_,
                    Significance = "Not tested", stringsAsFactors = FALSE)
  if (length(x) < min_n || length(y) < min_n) return(out)
  if (sd(x) == 0 || sd(y) == 0)                return(out)
  tt <- tryCatch(t.test(x, y, paired = FALSE, var.equal = FALSE),
                 error = function(e) NULL)
  if (is.null(tt)) return(out)
  out$t <- unname(tt$statistic)
  out$df <- unname(tt$parameter)
  out$p_value <- tt$p.value
  out$Significance <- if (tt$p.value < alpha) "Significant" else "Not significant"
  out
}

# descriptive statistics for one value vector, as one row:

summarise_values <- function(x)
  data.frame(n = length(x), Median = s_med(x), Mean = s_mn(x), SD = s_sd(x),
             Min = s_min(x), Max = s_max(x))

# round numeric columns to 3 significant figures at write time (counts kept):

fmt_stats <- function(tbl) {
  num <- names(tbl)[vapply(tbl, is.numeric, logical(1))]
  num <- num[!grepl("(^|[ _])n$|^n_", num)]
  tbl[num] <- lapply(tbl[num], \(x) signif(x, 3))
  tbl
}

# narrow a selection mask to the wanted fraction; warns on drop / fallback:

match_fraction <- function(sel, fraction, tf, label, wlabel, src, noun) {
  matched <- sel & fraction == tf
  if (any(matched)) {
    dropped <- sel & fraction != tf
    if (any(dropped)) warning(call. = FALSE, sprintf(
      "%s: dropped %d non-%s %s for '%s' (%s); kept %d (other fractions: %s)",
      src, sum(dropped), tf, noun, label, wlabel, sum(matched),
      paste(sort(unique(fraction[dropped])), collapse = "/")))
    return(matched)
  }
  if (any(sel)) warning(call. = FALSE, sprintf(
    "%s: no %s %s for '%s' (%s); using all fractions present (%s)",
    src, tf, noun, label, wlabel, paste(sort(unique(fraction[sel])), collapse = "/")))
  sel
}

# fraction to match a baseline/routine series against:

trial_fraction <- function(p) {
  if (!is.null(p$fraction) && !is.na(p$fraction) && nzchar(p$fraction)) return(p$fraction)
  if (grepl("_filt", p$trial_column)) "filtered" else "unfiltered"
}

# historical baseline values for one parameter, matched to the trial fraction:

hist_values <- function(hist, p, w) {
  md  <- as.integer(format(as.Date(hist$date), "%m%d"))
  sel <- hist$param_navn %in% p$hist_name_list[[1]] &
    hist$year >= w$year_min & hist$year <= w$year_max &
    md        >= w$md_min   & md        <= w$md_max
  if (w$locality_match != "all") sel <- sel & hist$locality == w$locality_match
  
  # fraction_match == "any" pools both fractions (NO3; see manuscript):
  if (!is.null(p$fraction_match) && identical(p$fraction_match, "any"))
    return(hist$value[sel] * p$conv)
  
  tf <- trial_fraction(p)
  sel <- match_fraction(sel, hist$fraction, tf, p$label, w$label,
                        "hist_values", "historical records")
  hist$value[sel] * p$conv
}

# same-sample ion set (by date) for the calcite SI, keyed to the registry:

CALCITE_PARAM_IDS <- c(pH = "pH", Alkalinity_mmolL = "TA", Ca_mgL = "Ca",
                       Mg_mgL = "Mg", Na_mgL = "Na", K_mgL = "K", Cl_mgL = "Cl", SO4_mgL = "SO4",
                       NO3N_ugL = "NO3", TOTP_ugPL = "TP")

hist_wide_ions <- function(hist, parm, w) {
  md <- as.integer(format(as.Date(hist$date), "%m%d"))
  cols <- lapply(names(CALCITE_PARAM_IDS), function(out_col) {
    p <- parm[parm$param_id == CALCITE_PARAM_IDS[[out_col]], ][1, ]
    sel <- hist$param_navn %in% p$hist_name_list[[1]] &
      hist$year >= w$year_min & hist$year <= w$year_max &
      md        >= w$md_min   & md        <= w$md_max &
      hist$locality == w$locality_match
    if (!identical(p$fraction_match, "any")) {
      tf  <- trial_fraction(p)
      sel <- match_fraction(sel, hist$fraction, tf, p$label, w$label,
                            "hist_wide_ions", "historical records")
    }
    agg <- aggregate(value ~ date, data.frame(date = hist$date[sel], value = hist$value[sel]),
                     FUN = median)
    names(agg)[names(agg) == "value"] <- out_col
    agg
  })
  Reduce(function(a, b) merge(a, b, by = "date"), cols)
}

# trial (Site KL) ion set for the calcite SI; alkalinity input is ANC-derived:

trial_wide_ions <- function(trial_win) {
  df <- read_csv_raw("data", "processed", "kvina_grab_samples_2025.csv")
  df <- df[df$station_label == "KL", , drop = FALSE]
  df <- clip_to_trial_window(df, "sample_date", trial_win)
  data.frame(
    pH = df$pH, Alkalinity_mmolL = parse_value(df$ANC_uEkv_L) / 1000,
    Ca_mgL = df$Ca_mgL, Mg_mgL = df$Mg_mgL, Na_mgL = df$Na_mgL, K_mgL = df$K_mgL,
    Cl_mgL = df$Cl_mgL, SO4_mgL = df$SO4_mgL, NO3N_ugL = parse_value(df$NO3N_ugL),
    TOTP_ugPL = df$TOTP_ugPL
  )
}

# calcite saturation index (PHREEQC; one SOLUTION block per sample, Al omitted) --

.calcite_solution_block <- function(i, row) {
  sprintf(paste(
    "SOLUTION %d",
    "    pH    %s",
    "    Alkalinity    %s mmol/l",
    "    K    %s",
    "    Ca    %s",
    "    Mg    %s",
    "    Na    %s",
    "    Cl    %s",
    "    S(6)    %s as SO4",
    "    P    %s ug/l as P",
    "    N    %s ug/l as N",
    "    units    mg/l",
    sep = "\n"
  ), i, row$pH, row$Alkalinity_mmolL, row$K_mgL, row$Ca_mgL, row$Mg_mgL,
  row$Na_mgL, row$Cl_mgL, row$SO4_mgL, row$TOTP_ugPL, row$NO3N_ugL)
}

si_calcite <- function(df) {
  need <- c("pH", "Alkalinity_mmolL", "Ca_mgL", "Mg_mgL", "Na_mgL", "K_mgL",
            "Cl_mgL", "SO4_mgL", "NO3N_ugL", "TOTP_ugPL")
  df <- df[stats::complete.cases(df[, need]), , drop = FALSE]
  if (nrow(df) == 0) return(numeric(0))
  
  blocks <- vapply(seq_len(nrow(df)), \(i) .calcite_solution_block(i, df[i, ]),
                   character(1))
  phr_input <- paste(c(blocks, "SELECTED_OUTPUT\n    saturation_indices    Calcite"),
                     collapse = "\n")
  
  phreeqc::phrLoadDatabaseString(phreeqc::phreeqc.dat)
  phreeqc::phrRunString(phr_input)
  out <- lapply(.Call("getSelOutLst", PACKAGE = "phreeqc"), as.data.frame)$n1
  out$si_Calcite
}
