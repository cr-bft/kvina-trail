#!/usr/bin/env Rscript
# R/01-analysis.R -- statistics pipeline (no figures). Section 01 writes the
# data/processed CSVs that 02-03 read; every output derives from the analyte
# registry (analytes.csv via load_analytes()). Methods: see manuscript.

suppressPackageStartupMessages({
  library("dplyr")
  library("tidyr")
  library("stringr")
  library("tibble")
  library("purrr")
})

source("R/functions-paper.R")

# Section 01: build tidy tables ----------------------------------------------

out_csv <- file.path("data", "processed", "kvina_grab_samples_2025.csv")
stations <- read_csv_raw(data_path("stations.csv"))

clean_names <- c(
  "pH" = "pH", "Kond. mS/m" = "cond_mS_m", "Alk-filt mmol/l" = "alk_filt_mmol_L",
  "TOC mg C/l" = "TOC_mgC_L", "DOC mg C/L" = "DOC_mgC_L", "DIC mg C/L" = "DIC_mgC_L",
  "Farge_F mg Pt/l" = "farge_mgPt_L", "DO mg O/l" = "DO_mgO_L", "ANC µEkv/L" = "ANC_uEkv_L",
  "Al µg/l" = "Al_ugL", "Al-filt µg/l" = "Al_filt_ugL", "Ca mg/L" = "Ca_mgL",
  "Ca-filt mg/l" = "Ca_filt_mgL", "Mg mg/L" = "Mg_mgL", "Mg-filt mg/l" = "Mg_filt_mgL",
  "Na mg/L" = "Na_mgL", "Na-filt mg/l" = "Na_filt_mgL", "K mg/L" = "K_mgL",
  "K-filt mg/l" = "K_filt_mgL", "SO4 mg/L" = "SO4_mgL", "NO3-N µg/l" = "NO3N_ugL",
  "Cl mg/L" = "Cl_mgL", "Fluorid mg/L" = "Fluorid_mgL", "PO4-P µg P/l" = "PO4P_ugPL",
  "TOTP µg P/l" = "TOTP_ugPL", "TSM mg/l" = "TSM_mgL", "Turbiditet FNU" = "turbiditet_FNU",
  "Ag µg/L" = "Ag_ugL", "Ag-filt µg/l" = "Ag_filt_ugL", "As µg/l" = "As_ugL",
  "As-filt µg/l" = "As_filt_ugL", "Cd µg/l" = "Cd_ugL", "Cd-filt µg/l" = "Cd_filt_ugL",
  "Co µg/L" = "Co_ugL", "Co-filt µg/l" = "Co_filt_ugL", "Cr µg/l" = "Cr_ugL",
  "Cr-filt µg/l" = "Cr_filt_ugL", "Cu µg/L" = "Cu_ugL", "Cu-filt µg/l" = "Cu_filt_ugL",
  "Fe µg/l" = "Fe_ugL", "Fe-filt µg/l" = "Fe_filt_ugL", "Mn µg/l" = "Mn_ugL",
  "Mn-filt µg/l" = "Mn_filt_ugL", "Ni µg/L" = "Ni_ugL", "Ni-filt µg/l" = "Ni_filt_ugL",
  "Pb µg/L" = "Pb_ugL", "Pb-filt µg/l" = "Pb_filt_ugL", "U µg/l" = "U_ugL",
  "U-filt µg/l" = "U_filt_ugL", "V µg/L" = "V_ugL", "V-filt µg/l" = "V_filt_ugL",
  "Zn µg/l" = "Zn_ugL"
)
param_cols <- unname(clean_names)

raw <- read_waterchem_sorted()

nm <- str_squish(names(raw))
names(raw) <- ifelse(nm %in% names(clean_names), clean_names[nm], nm)

# warn on any raw station code absent from stations.csv (its rows are dropped):
unmatched <- setdiff(unique(raw[["Station code"]]), stations$station_code)
unmatched <- unmatched[!is.na(unmatched) & nzchar(trimws(unmatched))]
if (length(unmatched)) {
  warning(call. = FALSE, sprintf(
    "01: %d raw 'Station code' value(s) absent from stations.csv; their rows dropped: %s",
    length(unmatched), paste(unmatched, collapse = ", ")
  ))
}

grab <- raw |>
  left_join(stations |> select(station_code, station_label, station_name),
    by = c("Station code" = "station_code")
  ) |>
  filter(!is.na(station_label)) |>
  mutate(sample_date = format(as.Date(`Sample date`), "%Y-%m-%d")) |>
  rename(station_code = `Station code`) |>
  select(station_code, station_label, station_name, sample_date, all_of(param_cols)) |>
  arrange(match(station_label, stations$station_label), sample_date)

if (nrow(grab) == 0L) {
  stop("01: no rows matched a station after the join on 'Station code'. Check that ",
    "the 'WaterChem Sorted' sheet still has a 'Station code' column whose values ",
    "match stations.csv.",
    call. = FALSE
  )
}

dir.create(file.path("data", "processed"), showWarnings = FALSE, recursive = TRUE)
write.csv(grab, out_csv, row.names = FALSE, fileEncoding = "UTF-8")
cat(sprintf("Wrote %d rows x %d cols -> %s\n", nrow(grab), ncol(grab), out_csv))
cat(sprintf(
  "Stations: %s | dates: %d (%s to %s)\n",
  paste(unique(grab$station_label), collapse = ", "),
  length(unique(grab$sample_date)), min(grab$sample_date), max(grab$sample_date)
))

parm <- load_analytes()

hist_names <- parm$hist_name_list |>
  unlist() |>
  str_trim() |>
  unique()
hist_names <- hist_names[nzchar(hist_names)]

historical_long <- read_vannmiljo() |>
  filter(Parameter_navn %in% hist_names) |>
  niva_long() |>
  filter(!is.na(value), !is.na(year))

# reactive aluminium only; no derived "total" Al (see manuscript).

# Site-KL 2025 routine supplement (same-method reactive Al / silica / total N),
# kept separate from the 2011-2024 baseline:
export_kl_2025_long <- read_vannmiljo() |>
  filter(
    Vannlokalitetsnavn == "Kvina ved Klosterøyna",
    substr(Tid_provetak, 1, 4) == "2025"
  ) |>
  niva_long() |>
  filter(!is.na(value), !is.na(date))

# Site-KL grab series for the historical-vs-trial tests:
trial_cols <- parm$trial_column[parm$in_hist_vs_trial & nzchar(parm$trial_column)]
trial_long <- read_csv_raw("data", "processed", "kvina_grab_samples_2025.csv") |>
  filter(station_label == "KL") |>
  mutate(date = as.Date(sample_date)) |>
  select(date, all_of(trial_cols)) |>
  mutate(across(all_of(trial_cols), as.character)) |>
  pivot_longer(
    cols = all_of(trial_cols), names_to = "trial_column",
    values_to = "value_raw"
  ) |>
  mutate(value = parse_value(value_raw)) |>
  filter(!is.na(value)) |>
  select(date, trial_column, value)

write.csv(historical_long, file.path("data", "processed", "historical_long.csv"),
  row.names = FALSE, fileEncoding = "UTF-8"
)
write.csv(trial_long, file.path("data", "processed", "trial_long.csv"),
  row.names = FALSE, fileEncoding = "UTF-8"
)
write.csv(export_kl_2025_long, file.path("data", "processed", "export_kl_2025_long.csv"),
  row.names = FALSE, fileEncoding = "UTF-8"
)

message(sprintf(
  "01: historical_long %d rows (%d source names) | trial_long %d rows (%d dates)",
  nrow(historical_long), length(hist_names),
  nrow(trial_long), length(unique(trial_long$date))
))
message(sprintf(
  "01: export_kl_2025_long %d rows (%d params, %d dates) -> 2025 KL same-method supplement",
  nrow(export_kl_2025_long), length(unique(export_kl_2025_long$param_navn)),
  length(unique(export_kl_2025_long$date))
))

# Section 02: canonical all-site summary -------------------------------------
# One long-form summary per site x analyte x period, from full-precision values.
# Windows come from TRIAL_WINDOW / HIST_WINDOW in functions.R.

parm <- load_analytes()
stn <- read_csv_raw(data_path("stations.csv"))
grab <- read_csv_raw("data", "processed", "kvina_grab_samples_2025.csv")
hist <- read_csv_raw("data", "processed", "historical_long.csv")
export2025 <- read_csv_raw("data", "processed", "export_kl_2025_long.csv")
dir.create(file.path("stats"), showWarnings = FALSE, recursive = TRUE)

trial_win <- TRIAL_WINDOW
w <- HIST_WINDOW
grab <- clip_to_trial_window(grab, "sample_date", trial_win)

md_lab <- function(md) format(as.Date(sprintf("2000%04d", md), "%Y%m%d"), "%b %d")
window_label <- paste(md_lab(trial_win$md_min), "-", md_lab(trial_win$md_max))
hist_period <- sprintf("Historical %d-%d", w$year_min, w$year_max)

SRC_GRAB <- "Kvina Karbon grab campaign"
SRC_NIVA <- "NIVA Vannmiljo routine monitoring"
SRC_DERIVED <- "Derived (PHREEQC)"

# station order upstream -> downstream -> tributary -> outlet, validated against
# the trial stations in stations.csv (rows with a "St." grab code):
STATIONS <- c("NY", "OB", "LA", "KL")
trial_labels <- stn$station_label[grepl("^St\\.", stn$station_code)]
stopifnot(
  "02: STATIONS must match the trial stations in stations.csv" =
    setequal(STATIONS, trial_labels)
)
stopifnot(
  "02: every registry trial_column must exist in the grab table" =
    all(parm$trial_column[nzchar(parm$trial_column)] %in% names(grab))
)

# measured values for one grab column at one station:
station_values <- function(col, conv, st) {
  na.omit(parse_value(grab[[col]][grab$station_label == st])) * conv
}

# 2025 routine values for one parameter (season-clipped, fraction-matched):
export_values <- function(edat, p, w) {
  md <- as.integer(format(as.Date(edat$date), "%m%d"))
  sel <- edat$param_navn %in% p$hist_name_list[[1]] & md >= w$md_min & md <= w$md_max
  if (identical(p$fraction_match, "any")) {
    return(edat$value[sel] * p$conv)
  }
  tf <- trial_fraction(p)
  sel <- match_fraction(
    sel, edat$fraction, tf, p$label, w$label,
    "export_values", "2025 routine records"
  )
  edat$value[sel] * p$conv
}

sum_row <- function(p, period, source, method, site, x) {
  bind_cols(
    tibble(
      ord = p$ord, Period = period, Source = source, Window = window_label,
      Site = site, Parameter = p$display_label, Fraction = p$fraction,
      Method = method, Unit = p$unit
    ),
    summarise_values(x)
  )
}

# trial 2025: every grab series at all four stations:
trial_parm <- parm[nzchar(parm$trial_column), ]
trial_rows <- map_dfr(seq_len(nrow(trial_parm)), function(i) {
  p <- trial_parm[i, ]
  map_dfr(STATIONS, \(st) sum_row(
    p, "Trial 2025", SRC_GRAB, "trial_grab", st,
    station_values(p$trial_column, p$conv, st)
  ))
})

# historical baseline: Site KL only; ICP-MS aluminium excluded (see manuscript):
hist_parm <- parm[nzchar(parm$historical_names) &
  !parm$param_id %in% c("AlTot", "AlDis"), ]
hist_rows <- map_dfr(seq_len(nrow(hist_parm)), function(i) {
  p <- hist_parm[i, ]
  sum_row(p, hist_period, SRC_NIVA, "routine_niva", "KL", hist_values(hist, p, w))
})

# routine 2025: same-method NIVA series at Site KL (5 in-window dates):
routine_parm <- parm[parm$in_routine_2025, ]
routine_rows <- map_dfr(seq_len(nrow(routine_parm)), function(i) {
  p <- routine_parm[i, ]
  sum_row(
    p, "Routine 2025", SRC_NIVA, "routine_niva", "KL",
    export_values(export2025, p, w)
  )
})

# calcite SI (derived), computed once and reused by section 03:
p_si <- parm[parm$param_id == "SI_calcite", ][1, ]
si_hist <- si_calcite(hist_wide_ions(hist, parm, w))
si_trial <- si_calcite(trial_wide_ions(trial_win))
si_rows <- bind_rows(
  sum_row(p_si, hist_period, SRC_DERIVED, "derived_phreeqc", "KL", si_hist),
  sum_row(p_si, "Trial 2025", SRC_DERIVED, "derived_phreeqc", "KL", si_trial)
)

summary_long <- bind_rows(trial_rows, hist_rows, routine_rows, si_rows) |>
  filter(n > 0) |>
  mutate(
    Period = factor(Period, levels = c(hist_period, "Trial 2025", "Routine 2025")),
    Site = factor(Site, levels = STATIONS)
  ) |>
  arrange(ord, Period, Site) |>
  mutate(Period = as.character(Period), Site = as.character(Site))

write.csv(fmt_stats(select(summary_long, -ord)),
  file.path("stats", "summary_stats_all_sites.csv"),
  row.names = FALSE, fileEncoding = "UTF-8"
)

message(sprintf(
  paste0(
    "02: summary_stats_all_sites %d rows ",
    "(%d trial / %d historical / %d routine / %d derived)"
  ),
  nrow(summary_long),
  sum(summary_long$Period == "Trial 2025" & summary_long$Method == "trial_grab"),
  sum(summary_long$Period == hist_period & summary_long$Method != "derived_phreeqc"),
  sum(summary_long$Period == "Routine 2025"),
  sum(summary_long$Method == "derived_phreeqc")
))

# Section 03: Welch t-tests --------------------------------------------------
# NY-vs-OB grabs; Site-KL historical vs 2025 grabs; historical vs 2025 routine
# (kept separate -- 5 routine obs, never merged with the 15-sample grab table).
# All via the shared welch_test(); reuses section 02's objects.

trial <- read_csv_raw("data", "processed", "trial_long.csv")
trial <- clip_to_trial_window(trial, "date", trial_win)

# trial values (Site KL) for one parameter:
trial_values <- function(p) trial$value[trial$trial_column == p$trial_column] * p$conv

# n/Median/Mean/SD block with a group prefix ("NY n", "Historical Median", ...):
stats4 <- function(x, prefix) {
  s <- summarise_values(x)[c("n", "Median", "Mean", "SD")]
  setNames(s, paste(prefix, names(s)))
}

ttest_row <- function(p, x, y, xlab, ylab) {
  bind_cols(
    tibble(Parameter = p$display_label, Fraction = p$fraction, Unit = p$unit),
    stats4(x, xlab), stats4(y, ylab), welch_test(x, y)
  )
}

run_ttests <- function(rows, x_fun, y_fun, xlab, ylab) {
  map_dfr(seq_len(nrow(rows)), function(i) {
    p <- rows[i, ]
    ttest_row(p, x_fun(p), y_fun(p), xlab, ylab)
  })
}

# move "Not tested" rows to the bottom, preserving ord within each group:
sink_not_tested <- function(tbl) {
  tbl[order(tbl$Significance == "Not tested"), , drop = FALSE]
}

ny_ob <- run_ttests(
  parm[parm$in_ny_ob, ],
  \(p) station_values(p$trial_column, p$conv, "NY"),
  \(p) station_values(p$trial_column, p$conv, "OB"),
  "NY", "OB"
)
write.csv(fmt_stats(sink_not_tested(ny_ob)), file.path("stats", "ttest_ny_vs_ob.csv"),
  row.names = FALSE, fileEncoding = "UTF-8"
)

# SI_calcite is derived (values from section 02), appended to the registry-driven rows:
hist_vs_trial <- bind_rows(
  run_ttests(
    parm[parm$in_hist_vs_trial & parm$param_id != "SI_calcite", ],
    \(p) hist_values(hist, p, w), trial_values, "Historical", "Trial"
  ),
  ttest_row(p_si, si_hist, si_trial, "Historical", "Trial")
)
write.csv(fmt_stats(sink_not_tested(hist_vs_trial)),
  file.path("stats", "ttest_kl_historical_vs_trial.csv"),
  row.names = FALSE, fileEncoding = "UTF-8"
)

hist_vs_routine <- run_ttests(
  parm[parm$in_routine_2025, ],
  \(p) hist_values(hist, p, w),
  \(p) export_values(export2025, p, w),
  "Historical", "Routine"
)
write.csv(fmt_stats(sink_not_tested(hist_vs_routine)),
  file.path("stats", "ttest_kl_historical_vs_routine_2025.csv"),
  row.names = FALSE, fileEncoding = "UTF-8"
)

cat("\n====  NY vs OB Welch t-tests (2025 grab campaign)  ====\n")
print(
  as.data.frame(fmt_stats(ny_ob))[, c(
    "Parameter", "Fraction", "Unit", "NY Mean",
    "OB Mean", "p_value", "Significance"
  )],
  row.names = FALSE
)
message(sprintf(
  "03: ttest_ny_vs_ob %d rows | hist_vs_trial %d rows | hist_vs_routine %d rows",
  nrow(ny_ob), nrow(hist_vs_trial), nrow(hist_vs_routine)
))
