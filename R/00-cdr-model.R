# R/00-cdr-model.R -- CDR quantification according to the Isometric 
# methodology: https://registry.isometric.com/protocol/river-alkalinity-enhancement/1.1
# Sources R/functions-cdr-model.R. Sourced with working directory = repo root.

# libraries ---------------------------------------------------------------

library("data.table") 
library("readxl")
library("janitor")
library("ggplot2")
library("ggtext")
library("patchwork")
library("reticulate") # see https://rstudio.github.io/reticulate/ for details on setting up a virtual Python environment
library("furrr")
library("brms") 
library("imputeTS")
library("withr")

source("R/functions-cdr-model.R")

# setup -------------------------------------------------------------------

plan("multisession")

theme_set(theme_bw())

# py_install("PyCO2SYS")
pyco2sys <- import("PyCO2SYS") 

# height-discharge data ---------------------------------------------------

stage_discharge <- fread("data-raw/nyland-flow-curve.csv")

model_discharge <- nls(
  discharge_m3_s ~ exp(x0 + x1 * water_level_cm + x2 * water_level_cm ^ 2 + x3 * water_level_cm ^ 3), 
  data = stage_discharge, 
  start = list(x0 = -3.3, x1 = 8e-2, x2 = -2e-4, x3 = 2e-7)
)

model_discharge_input <- data.table(water_level_cm = seq(0, stage_discharge[, max(water_level_cm)], length.out = 100))
model_discharge_input[, discharge_m3_s := predict(model_discharge, newdata = model_discharge_input)]

# plot:

stage_discharge |> 
  ggplot(aes(water_level_cm, discharge_m3_s)) +
  geom_point() +
  geom_line(data = model_discharge_input)

# parameters --------------------------------------------------------------

# time period:

reporting_period_start <- as.POSIXct("2025-06-10 00:00:00", tz = "Europe/Oslo") # first day of grab sampling
reporting_period_end <- as.POSIXct("2025-09-15 12:00:00", tz = "Europe/Oslo") # dosing returns to baseline

# physical:

mass_CO2 <- 44.009
mass_CaCO3 <- 100.0864
mass_dolomite <- 184.3994
mass_C <- 12.011

# discharge:

withdrawals_m3_s <- 0.3 # magnesium smelter (Eramet) withdraws a constant flow of 0.3 m3/s (as described in the PDD)
discharge_scaling_factor <- 1.744 # from NHC report ("nhc-report/R code and input files/Kvina.R", line 259)

# slurry specs:

slurry_specific_gravity <- 1.83 
slurry_solids_fraction <- 0.715 # see Rock and Mineral Feedstock Characterization Appendix
# calculation of carbon mass fraction, based on XRD data in Table A.27.2 of the Rock and Mineral Feedstock Characterization Appendix:
calcite_wt_pct <- 0.9608
dolomite_wt_pct <- 0.57e-2
carbon_mass_fraction_mg_kg <- calcite_wt_pct * (1e6 * mass_C / mass_CaCO3) + dolomite_wt_pct * (1e6 * 2 * mass_C / mass_dolomite)
# associated variance propagation, assuming independence:
sd_calcite_wt_pct <- 0.33e-2
sd_dolomite_wt_pct <- 0.18e-2
uncertainty_carbon_mass_fraction_mg_kg <- sqrt(sd_calcite_wt_pct ^ 2 * (1e6 * mass_C / mass_CaCO3) ^ 2 + sd_dolomite_wt_pct ^ 2 * (1e6 * 2 * mass_C / mass_dolomite) ^ 2)

# additional uncertainties (expressed as standard deviations or factors):

uncertainty_pH_measured <- 0.3 / 3 # NIVA uncertainty is 0.3 (3 SDs == 0.3) (data-raw/Data report Kvina.xlsx)
uncertainty_TA_fraction <- 0.2 / 3 # NIVA uncertainty is 20% (3 SDs == 20%) (data-raw/Data report Kvina.xlsx)

# renforth-henderson ------------------------------------------------------

# data source: approx. 3 years of data from ICOS OTC SOOP Release from G.O. Sars
# https://www.icos-cp.eu/data-services/about-data-portal

data_north_sea <- list.files("data-raw/ICOS_OTC_SOOP_Release_GO_Sars/", pattern = "^58G2", full.names = TRUE) |> 
  lapply(fread) |> 
  rbindlist(fill = TRUE) |> 
  clean_names()

data_north_sea_subset <- data_north_sea[
  depth_m < 10, 
  .SD, 
  .SDcols = grep("date_time|latitude|longitude|depth_m|p_sal_psu|p_co2_uatm|temp_deg_c", names(data_north_sea), value = TRUE)
]

data_north_sea_flags <- data_north_sea_subset[, .SD, .SDcols = grep("_flag$", names(data_north_sea_subset), value = TRUE)]

qc_pass <- data_north_sea_flags |> 
  apply(2, \(x) x == 2 | is.na(x)) |> 
  rowSums()

data_north_sea_qaqc <- data_north_sea_subset[qc_pass == ncol(data_north_sea_flags)]
data_north_sea_qaqc[, rh_factor := calculate_rh_factor(salinity = p_sal_psu, pco2 = p_co2_uatm, temperature = temp_deg_c)]

rh_factor <- data_north_sea_qaqc[, .(
  rh_factor = mean(rh_factor, na.rm = TRUE), 
  sd_rh_factor = sd(rh_factor, na.rm = TRUE),
  se_rh_factor = sd(rh_factor, na.rm = TRUE),
  n = length(na.omit(rh_factor))
)]

rh_factor[, se_rh_factor := sd_rh_factor / sqrt(n)][]

# carbonrun sensor data ---------------------------------------------------

data_sensor <- fread("data-raw/2025-11-19-carbon-run-sensor-data-temperature.csv")

# delivery data -----------------------------------------------------------

deliveries <- fread("data-raw/2025-11-25-feedstock-deliveries-nyland.csv")
deliveries <- deliveries[month(date_end) < 10]

# doser data --------------------------------------------------------------

data_doser <- setDT(read_excel("data-raw/Nyland, Kvina, 30.09.2025 00 - 15.05.2025 00 COMPLETE PROJECT DATA.xlsx", skip = 2))

data_doser[, Timestamp_Local := as.POSIXct(Date, format = "%d.%m.%Y %H", tz = "Europe/Oslo")]
data_doser <- data_doser[order(Timestamp_Local)]
data_doser[, delta_time_s := difftime(Timestamp_Local, shift(Timestamp_Local), units = "secs")]
data_doser[, delta_time_s := as.numeric(c(delta_time_s[2], delta_time_s[-1]))]

data_doser_historical <- setDT(read_excel("data-raw/Nyland dosererdata.xlsx", skip = 2))

data_doser_historical[, Timestamp_Local := as.POSIXct(Date, format = "%d.%m.%Y %H", tz = "Europe/Oslo")]

model_baseline_dose_input <- data_doser_historical[
  Timestamp_Local < reporting_period_start
][
  , .(
    Timestamp_Local,
    log_dose_l_h = log(`Dosering (liters per h)`), 
    log_discharge_m3_s = log(`Vannføring m3/sec`),
    month = factor(month(Timestamp_Local), labels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))
  )
][log_dose_l_h > 0]

model_baseline_test_data <- model_baseline_dose_input[
  year(Timestamp_Local) == 2025
][
  month(Timestamp_Local) %in% 1:5
]

model_baseline_dose_input <- model_baseline_dose_input[
  (year(Timestamp_Local) == 2024 & month(Timestamp_Local) != 1) | 
    (year(Timestamp_Local) == 2023 & month(Timestamp_Local) == 1)
]

model_baseline_dose <- brm(
  log_dose_l_h ~ log_discharge_m3_s + (1 + log_discharge_m3_s | month),
  data = model_baseline_dose_input,
  family = "student",
  cores = 4,
  control = list(adapt_delta = 0.99, max_treedepth = 14),
  file = "brms/model-baseline-dose",
  prior = c(prior(normal(0, 1), class = "b"), prior(normal(0, 1), class = "Intercept"))
)

# plot baseline dose:

model_baseline_dose_posterior <- posterior_epred(model_baseline_dose)
model_baseline_dose_fitted_values <- model_baseline_dose_posterior |> 
  apply(2, \(x) quantile(x, c(0.025, 0.5, 0.975))) |> 
  apply(2, exp)

model_baseline_dose_input |> 
  ggplot(aes(exp(log_discharge_m3_s), exp(log_dose_l_h))) +
  facet_wrap(vars(month)) +
  geom_point() +
  geom_ribbon(aes(ymin = model_baseline_dose_fitted_values[1, ], ymax = model_baseline_dose_fitted_values[3, ]), fill = "red", alpha = 0.3) +
  geom_line(aes(y = model_baseline_dose_fitted_values[2, ]), col = "red") +
  theme(
    axis.title.x = element_markdown(),
    axis.title.y = element_markdown(),
  ) +
  labs(
    x = "Discharge (m<sup>3</sup> s<sup>-1</sup>)",
    y = "Slurry dose (L h<sup>-1</sup>)",
  )

# ggsave("figures/baseline-dose-model.png", dpi = 600, width = 6.5, height = 4)

y <- exp(model_baseline_dose_input$log_dose_l_h)
yhat <- apply(model_baseline_dose_posterior, 2, \(x) exp(mean(x)))
rsq <- r_sq(y, yhat)

# test predictions (2025):

model_baseline_dose_test_predictions <- posterior_epred(model_baseline_dose, model_baseline_test_data)

model_baseline_test_data |> 
  ggplot(aes(exp(log_dose_l_h))) +
  facet_wrap(vars(month)) +
  geom_point(aes(y = exp(apply(model_baseline_dose_test_predictions, 2, mean)))) +
  scale_x_log10() +
  scale_y_log10() +
  geom_abline() +
  theme(
    axis.title.x = element_markdown(),
    axis.title.y = element_markdown(),
  ) +
  labs(
    x = "Discharge (m<sup>3</sup> s<sup>-1</sup>)",
    y = "Slurry dose (L h<sup>-1</sup>)",
  )

# get feedstock dose uncertainty ------------------------------------------

# combine dosing records:

data_dose_combined <- rbind(
  data_doser[, .(
    Timestamp_Local, 
    dose_ml_m3 = `Kalkdose korr.`,
    discharge_m3_s = Vannføring
  )], 
  data_doser_historical[Timestamp_Local < data_doser[, min(Timestamp_Local)], .(
    Timestamp_Local,
    dose_ml_m3 = `Limedose corrected (ml/m3)`,
    discharge_m3_s = `Vannføring m3/sec`
  )]
)[order(Timestamp_Local)]

data_dose_combined[, delta_time_s := difftime(Timestamp_Local, shift(Timestamp_Local), units = "secs")]
data_dose_combined[, delta_time_s := as.numeric(c(delta_time_s[2], delta_time_s[-1]))]

# merge with deliveries:

data_dose_merged <- data_dose_combined[
  , date_start := as.Date(Timestamp_Local, tz = "Europe/Oslo")
][
  date_start >= deliveries[, min(date_start)]
][
  date_start <= deliveries[, max(date_end)]
] |> 
  merge(deliveries, by = "date_start", all = TRUE)

data_dose_merged <- data_dose_merged[order(Timestamp_Local)]

data_dose_merged[, slurry_delivered_t := nafill(slurry_delivered_t, type = "locf")]
data_dose_merged[, date_end := nafill(date_end, type = "locf")]

# get rrmse:

model_dose_calibration_input <- data_dose_merged[
  , .(
    date_start = min(date_start),
    feedstock_dosed_t = 1e-6 * sum(dose_ml_m3 * discharge_m3_s * delta_time_s) * slurry_solids_fraction * slurry_specific_gravity,
    feedstock_delivered_t = unique(slurry_delivered_t) * slurry_solids_fraction
  ),
  by = "date_end"
]

model_dose_calibration <- lm(feedstock_delivered_t ~ 0 + feedstock_dosed_t, data = model_dose_calibration_input)

# plot calibration:

model_dose_calibration_input |> 
  ggplot(aes(feedstock_dosed_t, feedstock_delivered_t)) +
  geom_abline() +
  geom_point() +
  lims(x = c(50, 650), y = c(50, 650)) +
  labs(
    x = "Feedstock delivered (t)",
    y = "Feedstock dosed (t)",
  )

# ggsave("figures/feedstock-dose-calibration.png", width = 3, height = 2.5)

# uncertainties for dose and discharge:

uncertainty_dosing_fraction <- sqrt(mean(residuals(model_dose_calibration) ^ 2)) / model_dose_calibration_input[, mean(feedstock_delivered_t)]
uncertainty_flow_rate <- sigma(model_discharge) * discharge_scaling_factor

# discharge data ----------------------------------------------------------

data_discharge <- data_doser[
  , .(river_flow_m3_s = mean(Vannføring * discharge_scaling_factor - withdrawals_m3_s, na.rm = TRUE)),
  by = .(date = as.Date(Timestamp_Local, tz = "Europe/Oslo"))
]

# grab sample data --------------------------------------------------------

data_grab_2025 <- setDT(clean_names(read_excel("data-raw/Data report Kvina complete.xlsx"), replace = c("\u00b5" = "u", "pH" = "ph")))

attr(data_grab_2025$sample_date, "tzone") <- "Europe/Oslo"

data_grab_historical <- list.files("data-raw", pattern = "^Vann", full.names = TRUE) |> 
  setNames(list.files("data-raw/", pattern = "^Vann")) |> 
  lapply(fread) |> 
  rbindlist(idcol = "file")

# clean historical data:

data_grab_historical[, site := regmatches(file, m = regexpr("Kloster|Litl", file))]
data_grab_historical[, site := gsub("Litl", "LA", site)]
data_grab_historical[, site := gsub("Kloster", "KL", site)]
data_grab_historical[, param := regmatches(Parameter, m = regexpr("pH|Total alkali", Parameter))]
data_grab_historical[, param := gsub("Total alkali", "alkalinity_mmol_l", param)]
data_grab_historical[, date := as.Date(`Prøvetak. tidspunkt`, format = "%d-%m-%Y", tz = "Europe/Oslo")]

# baseline model input ----------------------------------------------------

model_baseline_input <- data_grab_historical[, .(date, site, param, value = Verdi)] |> 
  dcast(date ~ param + site, value.var = "value")

# model baseline ----------------------------------------------------------

model_baseline_pH <- lm(pH_KL ~ alkalinity_mmol_l_LA + pH_LA, data = model_baseline_input)
model_baseline_TA <- lm(alkalinity_mmol_l_KL ~ alkalinity_mmol_l_LA + pH_LA, data = model_baseline_input)

# error:

predictions_pH_loo_input <- na.omit(model_baseline_input[, .(alkalinity_mmol_l_LA, pH_LA, pH_KL)])
predictions_pH_loo <- loocv(predictions_pH_loo_input, response = "pH_KL")
uncertainty_pH <- get_error(predictions_pH_loo_input, "pH_KL", predictions_pH_loo, use = "complete")

predictions_TA_loo_input <- na.omit(model_baseline_input[, .(pH_LA, alkalinity_mmol_l_LA, alkalinity_mmol_l_KL)])
predictions_TA_loo <- loocv(predictions_TA_loo_input, response = "alkalinity_mmol_l_KL")
uncertainty_TA <- get_error(predictions_TA_loo_input, "alkalinity_mmol_l_KL", predictions_TA_loo)

# plot predictions --------------------------------------------------------

pH_error_labels <- data.table(
  x = -Inf, 
  y = Inf, 
  label = paste(names(uncertainty_pH), "=", lapply(uncertainty_pH, \(x) signif(x, 2)), collapse = "<br>")
)

p1 <- model_baseline_input[, .(alkalinity_mmol_l_LA, pH_LA, pH_KL)] |> 
  na.omit() |> 
  ggplot(aes(pH_KL, predictions_pH_loo)) +
  geom_abline() +
  geom_point() + 
  ggtext::geom_richtext(
    data = pH_error_labels,
    aes(x = x, y = y, label = label),
    hjust = "inward", vjust = "inward",
    label.size = 0, size = 2.5, alpha = 0.7
  ) +
  labs(x = "pH (observed)", y = "pH (predicted)")

TA_error_labels <- data.table(
  x = -Inf, 
  y = Inf, 
  label = paste(names(uncertainty_TA), "=", lapply(uncertainty_TA, \(x) signif(x, 2)), collapse = "<br>")
)

p2 <- model_baseline_input[, .(alkalinity_mmol_l_LA, alkalinity_mmol_l_KL)] |> 
  na.omit() |> 
  ggplot(aes(alkalinity_mmol_l_KL, predictions_TA_loo)) +
  geom_abline() +
  geom_point() +
  ggtext::geom_richtext(
    data = TA_error_labels,
    aes(x = x, y = y, label = label),
    hjust = "inward", vjust = "inward",
    label.size = 0, size = 2.5, alpha = 0.7
  ) +
  labs(x = "Total alkalinity<br>(observed, mmol L<sup>-1</sup>)", y = "Total alkalinity<br>(predicted, mmol L<sup>-1</sup>)")

# residuals plots:

p3_in <- data.table(rowid = seq(52), "Residuals\n(alkalinity model)" = c(residuals(model_baseline_TA), NA_real_), "Residuals\n(pH model)" = residuals(model_baseline_pH)) |> 
  melt(id.vars = "rowid")

setattr(p3_in$variable, "levels", c("Residuals\n(pH model)", "Residuals\n(alkalinity model)"))

p3 <- p3_in |> 
  ggplot(aes(sample = value)) +
  facet_wrap(vars(variable), scales = "free") +
  geom_qq() +
  geom_qq_line() +
  labs(x = "Standard normal quantiles", y = "Sample quantiles")

wrap_plots(p1, p2, p3, design = "AB\nCC") +
  plot_annotation(tag_levels = "a") &
  theme(
    plot.tag = element_text(face = "bold"),
    axis.title.x = element_markdown(),
    axis.title.y = element_markdown()
  )

# ggsave("figures/baseline-model-residuals.png", dpi = 600, width = 3.5, height = 2.5)

# partial residuals plots:

partial_residual <- function(fit, xname) {
  eps  <- resid(fit, "partial")
  data <- fit$model
  data.table(x = data[[xname]], partial_residual = eps[, xname])
}

p1 <- rbindlist(list(
  "Total alkalinity (mmol L<sup>-1</sup>)" = partial_residual(model_baseline_pH, "alkalinity_mmol_l_LA"), 
  "pH" = partial_residual(model_baseline_pH, "pH_LA")
), idcol = "param") |> 
  ggplot(aes(x, partial_residual)) +
  facet_wrap(vars(param), ncol = 2, scales = "free") +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  labs(x = NULL, y = "Partial residual")

p2 <- rbindlist(list(
  "Total alkalinity (mmol L<sup>-1</sup>)" = partial_residual(model_baseline_TA, "alkalinity_mmol_l_LA"), 
  "pH" = partial_residual(model_baseline_TA, "pH_LA")
), idcol = "param") |> 
  ggplot(aes(x, partial_residual)) +
  facet_wrap(vars(param), ncol = 2, scales = "free") +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  labs(x = NULL, y = "Partial residual")

wrap_plots(p1, p2, ncol = 1) +
  plot_annotation(tag_levels = "a") &
  theme(
    plot.tag = element_text(face = "bold"),
    strip.text = element_markdown()
  )

# ggsave("figures/baseline-model-partial-residuals.png", dpi = 600, width = 5.5, height = 4.5)

# clean temperature data --------------------------------------------------

# using temperature at KL:

data_sensor[variable == "Temp_degC"][Site == "KL"][Flag_NOAA == 2] |> 
  ggplot(aes(Timestamp_Local, value)) +
  geom_line()

data_temperature <- data_sensor[variable == "Temp_degC"][Site == "KL"][Flag_NOAA == 2][
  , .(temp_median_deg_c = median(value, na.rm = TRUE)),
  by = .(date = as.Date(Timestamp_Local, tz = "Europe/Oslo"))
]

temp_lab_deg_c <- 23

# subset grab sample data -------------------------------------------------

data_grab_subset_2025 <- data_grab_2025[
  station_code %in% c("St. 4b", "St. 8"),
  .(
    date = as.Date(sample_date, tz = "Europe/Oslo"), 
    station_code = fifelse(station_code == "St. 8", "LA", "KL"),
    pH = as.numeric(ph),
    alkalinity_mmol_l = as.numeric(alk_filt_mmol_l),
    dic_mg_l = dic_mg_c_l
  )
] |> 
  merge(data_temperature[, .(date, temp_deg_c = temp_median_deg_c)], by = "date", all.x = TRUE)

# clean up dic data:

data_grab_subset_2025[, dic_mg_l := gsub("<", "", dic_mg_l)]
data_grab_subset_2025[, dic_mg_l := as.numeric(gsub(",", ".", dic_mg_l))]

# convert units:

data_grab_subset_2025[, alkalinity_umol_kg := 1e3 * alkalinity_mmol_l / rho(temp_lab_deg_c)][]
data_grab_subset_2025[, dic_umol_kg := 1e3 * dic_mg_l / 12.011 / rho(temp_lab_deg_c)][]

# make wide:

data_grab_subset_2025 <- dcast(
  data_grab_subset_2025, 
  date + temp_deg_c ~ station_code, 
  value.var = c("pH", "alkalinity_umol_kg", "alkalinity_mmol_l", "dic_umol_kg")
)

# predict baseline --------------------------------------------------------

# make predictions:

data_grab_subset_2025[, pH_KL_baseline := predict(model_baseline_pH, newdata = data_grab_subset_2025)][]
data_grab_subset_2025[, alkalinity_mmol_l_KL_baseline := predict(model_baseline_TA, newdata = data_grab_subset_2025)][]

# convert units:

data_grab_subset_2025[, alkalinity_umol_kg_KL_baseline := 1e3 * alkalinity_mmol_l_KL_baseline / rho(temp_lab_deg_c)][]

# check that 2025 data are within the range of the training data:

predictions_pH_loo_input_ranges <- predictions_pH_loo_input[, lapply(.SD, range), .SDcols = c("alkalinity_mmol_l_LA", "pH_LA")]
data_grab_subset_2025_ranges <- data_grab_subset_2025[, lapply(.SD, range), .SDcols = c("alkalinity_mmol_l_LA", "pH_LA")]

stopifnot(all((predictions_pH_loo_input_ranges[1, ] <=  data_grab_subset_2025_ranges[1, ])[1, ]))
stopifnot(all((predictions_pH_loo_input_ranges[2, ] >=  data_grab_subset_2025_ranges[2, ])[1, ]))
stopifnot(min(data_grab_subset_2025$alkalinity_mmol_l_KL_baseline) >= min(predictions_TA_loo_input$alkalinity_mmol_l_KL))
stopifnot(max(data_grab_subset_2025$alkalinity_mmol_l_KL_baseline) <= max(predictions_TA_loo_input$alkalinity_mmol_l_KL))

# check that the post-project pH/TA data are within 2 SDs of the baseline model predictions:

stopifnot(data_grab_subset_2025[date > as.Date(reporting_period_end, tz = "Europe/Oslo"), abs(pH_KL - pH_KL_baseline) / uncertainty_pH$RMSE] < 2)
stopifnot(data_grab_subset_2025[date > as.Date(reporting_period_end, tz = "Europe/Oslo"), abs(alkalinity_umol_kg_KL - alkalinity_umol_kg_KL_baseline) / (1e3 * uncertainty_TA$RMSE / rho(temp_lab_deg_c))] < 2)

# calculate DIC -----------------------------------------------------------

# treatment:

pyco2sys_out_treatment <- with(data_grab_subset_2025, pyco2sys$sys(
  par1 = pH_KL, 
  par2 = alkalinity_umol_kg_KL,
  temperature = temp_lab_deg_c,
  temperature_out = temp_deg_c,
  par1_type = 3, # # pH
  par2_type = 1, # alkalinity
  salinity = 0,
  opt_k_carbonic = 8,
  uncertainty_into = c("alkalinity", "dic"),
  uncertainty_from = list(par1 = uncertainty_pH_measured, par2 = uncertainty_TA_fraction *  alkalinity_umol_kg_KL)
))

# baseline:

pyco2sys_out_baseline <- with(data_grab_subset_2025, pyco2sys$sys(
  par1 = pH_KL_baseline, 
  par2 = alkalinity_umol_kg_KL_baseline,
  temperature = temp_lab_deg_c,
  temperature_out = temp_deg_c,
  par1_type = 3, # # pH
  par2_type = 1, # alkalinity
  salinity = 0,
  opt_k_carbonic = 8,
  uncertainty_into = c("alkalinity", "dic"),
  uncertainty_from = list(
    par1 = sqrt(uncertainty_pH_measured ^ 2 + uncertainty_pH$RMSE ^ 2), 
    par2 = sqrt((uncertainty_TA_fraction *  alkalinity_umol_kg_KL) ^ 2 + (1e3 * uncertainty_TA$RMSE / rho(temp_lab_deg_c)) ^ 2)
  )
))

# Site LA:

pyco2sys_out_LA <- with(data_grab_subset_2025, pyco2sys$sys(
  par1 = pH_LA, 
  par2 = alkalinity_umol_kg_LA,
  temperature = temp_lab_deg_c,
  temperature_out = temp_deg_c,
  par1_type = 3, # # pH
  par2_type = 1, # alkalinity
  salinity = 0,
  opt_k_carbonic = 8,
  uncertainty_into = c("alkalinity", "dic"),
  uncertainty_from = list(par1 = uncertainty_pH_measured, par2 = uncertainty_TA_fraction *  alkalinity_umol_kg_LA)
))

# add to dataframe:

data_grab_subset_2025[, dic_umol_kg_KL_predicted := pyco2sys_out_treatment$dic][]
data_grab_subset_2025[, dic_umol_kg_KL_baseline := pyco2sys_out_baseline$dic][]
data_grab_subset_2025[, dic_umol_kg_LA_predicted := pyco2sys_out_LA$dic][]

data_grab_subset_2025[, sd_dic_umol_kg_KL_predicted := pyco2sys_out_treatment$u_dic][]
data_grab_subset_2025[, sd_dic_umol_kg_KL_baseline := pyco2sys_out_baseline$u_dic][]

data_grab_subset_2025[, sd_alkalinity_umol_kg_KL := pyco2sys_out_treatment$u_alkalinity][]
data_grab_subset_2025[, sd_alkalinity_umol_kg_KL_baseline := pyco2sys_out_baseline$u_alkalinity][]

# remove extra columns:

data_co2_removal <- data_grab_subset_2025[
  date <= as.Date(reporting_period_end, tz = "Europe/Oslo"), 
  -grep("mmol_l", names(data_grab_subset_2025), value = TRUE), 
  with = FALSE
]

# do linear interpolation:

data_co2_removal_interpolated <- merge(
  data.table(date = seq(as.Date(reporting_period_start, tz = "Europe/Oslo"), as.Date(reporting_period_end, tz = "Europe/Oslo"), by = "1 day")),
  data_co2_removal, by = "date", all.x = TRUE
)[order(date)]

data_co2_removal_interpolated[, delta_time_s := as.numeric(difftime(date, shift(date), units = "secs"))]
data_co2_removal_interpolated[, delta_time_s := c(delta_time_s[2], delta_time_s[-1])]

interpolate_these <- names(data_co2_removal_interpolated)[sapply(data_co2_removal_interpolated, function(u) sum(is.na(u))) > 0]

data_co2_removal_interpolated[, `:=`((interpolate_these), lapply(.SD, imputeTS::na_interpolation)), .SDcols = interpolate_these]

# aggregate feedstock totals ----------------------------------------------

baseline_dose_predicted_l_h <- posterior_epred(
  model_baseline_dose,
  data.frame(
    log_discharge_m3_s = log(data_doser$Vannføring),
    month = factor(
      month(data_doser$Timestamp_Local), 
      labels = c("May", "Jun", "Jul", "Aug", "Sep")
    )
  )
)

baseline_dose_predicted_l_s <- exp(apply(baseline_dose_predicted_l_h, 2, mean)) / 3600

data_doser[, feedstock_dose_t_baseline := 1e-3 * baseline_dose_predicted_l_s * delta_time_s * slurry_specific_gravity * slurry_solids_fraction]

data_doser[, feedstock_dose_t_total := 1e-6 * `Kalkdose korr.` * 
             `Vannføring` * delta_time_s * slurry_specific_gravity * slurry_solids_fraction]

data_doser[, feedstock_dose_t := feedstock_dose_t_total - feedstock_dose_t_baseline]

data_doser_aggregated <- data_doser[
  Timestamp_Local < reporting_period_end, 
  .(
    feedstock_dose_t_total = sum(feedstock_dose_t_total),
    feedstock_dose_t = sum(feedstock_dose_t), 
    feedstock_dose_t_baseline = sum(feedstock_dose_t_baseline)
  ), 
  by = .(date = as.Date(Timestamp_Local, tz = "Europe/Oslo"))
]

data_co2_removal_interpolated <- merge(data_co2_removal_interpolated, data_doser_aggregated, by = "date")

# add discharge data ------------------------------------------------------

data_co2_removal_interpolated <- merge(data_co2_removal_interpolated, data_discharge, by = "date")

# calculate cdr -----------------------------------------------------------

# inputs:

cdr_input <- copy(data_co2_removal_interpolated)

feedstock_mass_kg <- sum(cdr_input$feedstock_dose_t) * 1e3

water_density_kg_m3 <- cdr_input[, 1e3 * rho(temp_deg_c)]

river_flow_downstream_m3_s <- cdr_input[, river_flow_m3_s]

dic_treatment_mmol_kg <- cdr_input[, 1e-3 * dic_umol_kg_KL_predicted]
# dic_treatment_mmol_kg <- cdr_input[, 1e-3 * dic_umol_kg_KL] # use measured DIC instead (no difference to CDR)
dic_baseline_mmol_kg <- cdr_input[, 1e-3 * dic_umol_kg_KL_baseline]

alk_treatment_mmol_kg <- cdr_input[, 1e-3 * alkalinity_umol_kg_KL]
alk_baseline_mmol_kg <- cdr_input[, 1e-3 * alkalinity_umol_kg_KL_baseline]

time_interval_s <- cdr_input[, delta_time_s]

# calculation:

isometric_co2e_stored <- co2e_stored(feedstock_mass_kg, carbon_mass_fraction_mg_kg, rh_factor$rh_factor, dic_treatment_mmol_kg, dic_baseline_mmol_kg, alk_treatment_mmol_kg, alk_baseline_mmol_kg, water_density_kg_m3, river_flow_downstream_m3_s, time_interval_s)

isometric_co2e_stored

# alternate calculation based on alkalinity should match, as long as ocean retention factors are < 1:

data_co2_removal_interpolated[, .(
  co2e_stored = mass_CO2 * (
    rh_factor$rh_factor * sum((alkalinity_umol_kg_KL - alkalinity_umol_kg_KL_baseline) * 1e3 * rho(temp_deg_c) * river_flow_m3_s * delta_time_s * 1e-12) - 
      sum(feedstock_dose_t * carbon_mass_fraction_mg_kg * 1e-6) / mass_C
  )
)]

# monte carlo -------------------------------------------------------------

with_seed(1245324, {
  cdr_mc <- future_map(seq(1e4), \(x) {
    
    cdr_input <- copy(data_co2_removal_interpolated)
    
    # add feedstock:
    
    baseline_dose_predicted_l_s <- exp(baseline_dose_predicted_l_h[sample(seq(nrow(baseline_dose_predicted_l_h)), 1), ]) / 3600
    
    data_doser_mc <- copy(data_doser)
    
    data_doser_mc[, feedstock_dose_t_baseline := 1e-3 * baseline_dose_predicted_l_s * delta_time_s * slurry_specific_gravity * slurry_solids_fraction]
    
    data_doser_mc[, feedstock_dose_t_total := 1e-6 * `Kalkdose korr.` * `Vannføring` * delta_time_s * slurry_specific_gravity * slurry_solids_fraction]
    
    data_doser_mc[, feedstock_dose_t_total := rnorm(feedstock_dose_t_total, mean = feedstock_dose_t_total, sd = uncertainty_dosing_fraction * feedstock_dose_t_total)]
    
    data_doser_mc[, feedstock_dose_t := feedstock_dose_t_total - feedstock_dose_t_baseline]
    
    data_doser_aggregated_mc <- data_doser_mc[
      Timestamp_Local < reporting_period_end, 
      .(
        feedstock_dose_t_total = sum(feedstock_dose_t_total),
        feedstock_dose_t = sum(feedstock_dose_t), 
        feedstock_dose_t_baseline = sum(feedstock_dose_t_baseline)
      ), 
      by = .(date = as.Date(Timestamp_Local, tz = "Europe/Oslo"))
    ]
    
    cdr_input <- merge(cdr_input[, -c("feedstock_dose_t_total", "feedstock_dose_t", "feedstock_dose_t_baseline")], data_doser_aggregated_mc, by = "date")
    
    # get feedstock totals:
    
    carbon_mass_fraction_mg_kg_rn <- rnorm(1, mean = carbon_mass_fraction_mg_kg, sd = uncertainty_carbon_mass_fraction_mg_kg)
    
    feedstock_mass_kg <- sum(cdr_input$feedstock_dose_t * 1e3)
    
    water_density_kg_m3 <- cdr_input[, 1e3 * rho(temp_deg_c)]
    
    # flow:
    
    data_doser_mc <- copy(data_doser)
    
    data_discharge_mc <- data_doser_mc[
      , river_flow_m3_s := Vannføring * discharge_scaling_factor
    ][
      , river_flow_m3_s := rnorm(river_flow_m3_s, mean = river_flow_m3_s, sd = uncertainty_flow_rate)
    ][
      , .(river_flow_m3_s = mean(river_flow_m3_s - withdrawals_m3_s, na.rm = TRUE)),
      by = .(date = as.Date(Timestamp_Local, tz = "Europe/Oslo"))
    ]
    
    cdr_input <- merge(cdr_input[, !"river_flow_m3_s"], data_discharge_mc, by = "date")
    
    river_flow_downstream_m3_s <- cdr_input[, river_flow_m3_s]
    
    # DIC/TA:
    
    dic_treatment_mmol_kg <- cdr_input[, 1e-3 * dic_umol_kg_KL_predicted]
    dic_treatment_mmol_kg <- rnorm(dic_treatment_mmol_kg, mean = dic_treatment_mmol_kg, sd = 1e-3 * cdr_input$sd_dic_umol_kg_KL_predicted)
    dic_baseline_mmol_kg <- cdr_input[, 1e-3 * dic_umol_kg_KL_baseline]
    dic_baseline_mmol_kg <- rnorm(dic_baseline_mmol_kg, mean = dic_baseline_mmol_kg, sd = 1e-3 * cdr_input$sd_dic_umol_kg_KL_baseline)
    
    alk_treatment_mmol_kg <- cdr_input[, 1e-3 * alkalinity_umol_kg_KL]
    alk_treatment_mmol_kg <- rnorm(alk_treatment_mmol_kg, mean = alk_treatment_mmol_kg, sd = 1e-3 * cdr_input$sd_alkalinity_umol_kg_KL)
    alk_baseline_mmol_kg <- cdr_input[, 1e-3 * alkalinity_umol_kg_KL_baseline]
    alk_baseline_mmol_kg <- rnorm(alk_baseline_mmol_kg, mean = alk_baseline_mmol_kg, sd = 1e-3 * cdr_input$sd_alkalinity_umol_kg_KL_baseline)
    
    time_interval_s <- cdr_input[, delta_time_s]
    
    rh_factor_rn <- rnorm(1, mean = rh_factor$rh_factor, sd = rh_factor$se_rh_factor)
    
    output <- co2e_stored(
      feedstock_mass_kg, 
      carbon_mass_fraction_mg_kg_rn, 
      rh_factor_rn, 
      dic_treatment_mmol_kg, 
      dic_baseline_mmol_kg, 
      alk_treatment_mmol_kg, 
      alk_baseline_mmol_kg, 
      water_density_kg_m3, 
      river_flow_downstream_m3_s, 
      time_interval_s
    )
    
    output$feedstock_mass_baseline_kg <- sum(cdr_input$feedstock_dose_t_baseline * 1e3)
    output$feedstock_mass_total_kg <- sum(cdr_input$feedstock_dose_t_total * 1e3)
    output$rh_factor <- rh_factor_rn
    
    output
    
  }, .progress = TRUE, .options = furrr::furrr_options(seed = TRUE)) |> 
    rbindlist()
})

cdr_mc_sd <- apply(cdr_mc, 2, sd) |> 
  t() |> 
  as.data.frame()

cdr_mc_summary <- cdr_mc[, .(co2e_stored_t = mean(co2e_stored_t), sd_co2e_stored_t = sd(co2e_stored_t), p16 = quantile(co2e_stored_t, 0.16))]

cdr_mc_summary

# isometric inputs --------------------------------------------------------

input_isometric_scalar <- isometric_co2e_stored[, .(
  carbon_mass_fraction_mg_kg,
  feedstock_mass_treatment_kg = 1e3 * data_doser[Timestamp_Local < reporting_period_end & Timestamp_Local >= reporting_period_start, sum(feedstock_dose_t_baseline + feedstock_dose_t)],
  feedstock_mass_kg_counterfactual = 1e3 * data_doser[Timestamp_Local < reporting_period_end & Timestamp_Local >= reporting_period_start, sum(feedstock_dose_t_baseline)],
  ocean_retention_treatment,
  ocean_retention_baseline,
  # final monitoring point is the river mouth, so there are no riverine losses:
  river_retention_treatment = 1,
  river_retention_baseline = 1,
  uncertainty_discount_t = cdr_mc_summary[, sd_co2e_stored_t]
)]

input_isometric_vector <- data_co2_removal_interpolated[, .(
  date,
  water_density_kg_m3 = 1e3 * rho(temp_deg_c),
  dic_treatment_mmol_kg = dic_umol_kg_KL_predicted * 1e-3,
  dic_counterfactual_mmol_kg = dic_umol_kg_KL_baseline * 1e-3,
  river_flow_m3_s = river_flow_m3_s,
  time_interval_s = delta_time_s
)]
