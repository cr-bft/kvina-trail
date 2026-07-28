# R/02-figures-A.R -- key manuscript/supplement figures and analysis. Sources R/00-cdr-model.R
# (which sources R/functions-cdr-model.R)

# code from source --------------------------------------------------------

dir.create("figures/", recursive = TRUE, showWarnings = FALSE)
dir.create("figures/manuscript", recursive = TRUE, showWarnings = FALSE)
dir.create("figures/supporting-information", recursive = TRUE, showWarnings = FALSE)

source("R/00-cdr-model.R")

# libraries ---------------------------------------------------------------

# remotes::install_github("cran/phreeqc")
library("phreeqc")
library("glue")

# setup -------------------------------------------------------------------

# shared site colour palette (Okabe-Ito):

site_colours <- c(NY = "#E69F00", OB = "#56B4E9", LA = "#009E73", KL = "#0072B2")

options(datatable.week = "sequential")

# read data ---------------------------------------------------------------

data_historical <- fread("data-raw/Kvina historical_data_2011_2025.csv")
data_nve <- list.files(path = "data-raw", pattern = "Discharge-day", full.names = TRUE) |>
  setNames(nm = _) |>
  lapply(fread)

data_krageholen_1 <- data_nve$`data-raw/Krågehølen 25.51.0-Discharge-day-v1.csv`
data_krageholen_2 <- data_nve$`data-raw/Krågehølen 25.51.0-Discharge-day-v1-new.csv`
data_krageholen <- rbind(data_krageholen_1, data_krageholen_2[Time > data_krageholen_1[, max(Time)]])
data_roynestad <- data_nve$`data-raw/Røynestad 25.85.0-Discharge-day-v0.csv`

setnames(data_krageholen, "Discharge (m³/s)", "discharge_krageholen_m3_s")
setnames(data_roynestad, "Discharge (m³/s)", "discharge_roynestad_m3_s")

# PHREEQC:

phr_kinetic_model <- readLines("phreeqc/kinetic-model.txt")
phr_calcite_addition <- readLines("phreeqc/calcite-addition.txt")
phr_calcite_saturation_template <- readLines("phreeqc/calcite-saturation.txt")

# parameters --------------------------------------------------------------

dose_cdr_mass_ratio <- 1 / 0.307 # feedstock-limited, does not consider water chemistry
mass_Ca <- 40.0784
mass_CaCO3 <- 100.0864

# plot inputs -------------------------------------------------------------

plot_removal_data_input <- melt(
  data_co2_removal_interpolated,
  id.vars = "date",
  measure.vars = grep("^pH|^dic|^alkalinity|flow", names(data_co2_removal_interpolated), value = TRUE, perl = TRUE)
)[
  , param := regmatches(variable, m = regexpr("[^K^L]+(?=_)", variable, perl = TRUE))
][
  , site := gsub("flow", "KL", regmatches(variable, m = regexpr("[KL][LA]_?[a-z]{0,}|flow", variable)))
][
  , type := gsub("^$|^s$", "Measured", regmatches(variable, m = regexpr("[a-z]{0,}$", variable)))
][
  , site := regmatches(site, m = regexpr("[KL][LA]", site))
][
  param != "dic_umol_kg" | (param == "dic_umol_kg" & type != "Measured")
][
  , type_plot := fcase(
    type == "predicted", "Calculated",
    type == "baseline", "Counterfactual",
    type == "measured", "Measured",
    default = type
  )
][
  # counterfactuals represent the KL alternative and retain its blue colour.
  , colour_site := fifelse(type_plot == "Counterfactual", "KL", site)
][
  # river flow is contextual rather than a site-chemistry series.
  , colour_site := fifelse(param == "river_flow_m3", "flow", colour_site)
][
  , scenario := fifelse(type_plot == "Counterfactual", "Counterfactual", "Trial")
]

plot_input <- data_co2_removal_interpolated[, .(
  river_flow_m3_s,
  pH_KL,
  co2e_stored = mass_CO2 * (
    rh_factor$rh_factor * sum((alkalinity_umol_kg_KL - alkalinity_umol_kg_KL_baseline) * 1e3 * rho(temp_deg_c) * river_flow_m3_s * delta_time_s * 1e-12) -
      sum(feedstock_dose_t * carbon_mass_fraction_mg_kg * 1e-6) / mass_C
  ),
  co2e_stored_upper_bound = feedstock_dose_t / dose_cdr_mass_ratio
),
by = "date"
]

# figure 2 ----------------------------------------------------------------

observed_colour <- "#D55E00" # vermilion (Okabe-Ito)

cdr_point_estimate_t <- plot_input[, sum(co2e_stored)]

p1_by_day <- plot_input |>
  ggplot(aes(date, co2e_stored_upper_bound)) +
  geom_col(aes(date, co2e_stored, fill = "Observed")) +
  geom_line(aes(col = "Potential")) +
  scale_color_manual(values = c(Potential = "grey50")) +
  scale_fill_manual(values = c(Observed = observed_colour)) +
  theme(
    legend.text = element_markdown(),
    strip.text = element_markdown(),
    legend.position = "bottom",
    legend.margin = margin()
  ) +
  labs(x = NULL, y = "t CO<sub>2</sub>e<sub>stored</sub>", fill = NULL, col = NULL)

p2_by_day <- plot_input[
  order(date), .(
    date,
    "Observed" = cumsum(co2e_stored),
    "Potential" = cumsum(co2e_stored_upper_bound)
  )
] |>
  melt(id.vars = "date") |>
  ggplot(aes(date, value, col = variable)) +
  geom_line() +
  scale_color_manual(values = c(Observed = observed_colour, Potential = "grey50")) +
  theme(
    axis.title.y = element_markdown(),
    legend.text = element_markdown(),
    legend.position = "bottom"
  ) +
  labs(x = NULL, y = "t CO<sub>2</sub>e<sub>stored</sub>", col = NULL)

p3 <- ggplot(data.table(x = cdr_mc$co2e_stored_t), aes(x)) +
  geom_histogram(bins = 30, fill = "grey30", colour = NA) +
  geom_vline(xintercept = cdr_point_estimate_t, colour = "white", linewidth = 0.9) +
  theme(axis.title.x = element_markdown()) +
  labs(x = "t CO<sub>2</sub>e<sub>stored</sub>", y = "Count")

p4 <- plot_removal_data_input |>
  ggplot(aes(date, value, col = colour_site, linetype = scenario)) +
  facet_wrap(
    vars(param),
    scales = "free_y", ncol = 1,
    labeller = as_labeller(c(
      alkalinity_umol_kg = "Alkalinity (&mu;mol kg<sup>-1</sup>)",
      dic_umol_kg = "DIC (&mu;mol kg<sup>-1</sup>)",
      pH = "pH",
      river_flow_m3 = "River flow (m<sup>3</sup> s<sup>-1</sup>)"
    ))
  ) +
  geom_line() +
  scale_color_manual(
    values = c(KL = unname(site_colours["KL"]), LA = unname(site_colours["LA"]), flow = "grey50"),
    breaks = c("KL", "LA"),
    labels = c(KL = "KL (river mouth)", LA = "LA (limed tributary)")
  ) +
  scale_linetype_manual(values = c(Trial = "solid", Counterfactual = "dashed")) +
  theme(
    strip.text = element_markdown(),
    legend.position = "right",
    legend.text = element_markdown()
  ) +
  labs(x = NULL, y = NULL, linetype = "Scenario", col = "Site")

wrap_plots(p1_by_day, p2_by_day, p3, p4, design = "AD\nBD\nCD\n") +
  plot_annotation(tag_levels = "a") &
  theme(
    plot.tag = element_text(face = "bold"),
    axis.title.y = element_markdown()
  )

ggsave("figures/manuscript/figure-2.png", width = 6.5, height = 6.5, dpi = 600, bg = "white")

# table 2 -----------------------------------------------------------------

# alkalinity (total):

isometric_co2e_stored[, signif(river_alk_export_treatment_mmol * 1e-9 * mass_CO2, 3)]
cdr_mc[, sd(river_alk_export_treatment_mmol) * 1e-9 * mass_CO2]

# alkalinity (CDR scenario):

isometric_co2e_stored[, signif((river_alk_export_treatment_mmol - river_alk_export_baseline_mmol) * 1e-9 * mass_CO2, 3)]
cdr_mc[, sd(river_alk_export_treatment_mmol - river_alk_export_baseline_mmol) * 1e-9 * mass_CO2]

# alkalinity (BAU liming):

alkalinity_export_bau <- input_isometric_scalar[, 2 * feedstock_mass_kg_counterfactual * carbon_mass_fraction_mg_kg * 1e-9 * mass_CO2 / mass_C]
alkalinity_export_bau
cdr_mc[, 2 * sd(feedstock_mass_baseline_kg * carbon_mass_fraction_mg_kg) * 1e-9 * mass_CO2 / mass_C]

# alkalinity (upstream baseline):

isometric_co2e_stored[, signif(river_alk_export_baseline_mmol * 1e-9 * mass_CO2, 3)] - alkalinity_export_bau
cdr_mc[, mass_CO2 * 1e-9 * sd(river_alk_export_baseline_mmol - 2 * feedstock_mass_baseline_kg * carbon_mass_fraction_mg_kg / mass_C)]

# discharge:

data_doser[
  Timestamp_Local >= reporting_period_start & Timestamp_Local < reporting_period_end,
  .(
    range_river_flow_m3_s = range(Vannføring * discharge_scaling_factor - withdrawals_m3_s, na.rm = TRUE),
    mean_river_flow_m3_s = mean(Vannføring * discharge_scaling_factor - withdrawals_m3_s, na.rm = TRUE)
  )
]

# feedstock (total):

input_isometric_scalar[, feedstock_mass_treatment_kg * 1e-3]
cdr_mc[, sd(feedstock_mass_total_kg) * 1e-3] # SD

# feedstock (CDR scenario):

input_isometric_scalar[, (feedstock_mass_treatment_kg - feedstock_mass_kg_counterfactual) * 1e-3]
cdr_mc[, sd(feedstock_dose_kg) * 1e-3] # SD

# feedstock (BAU liming):

input_isometric_scalar[, (feedstock_mass_kg_counterfactual) * 1e-3]
cdr_mc[, sd(feedstock_mass_baseline_kg) * 1e-3] # SD

# C_feedstock (total):

input_isometric_scalar[, feedstock_mass_treatment_kg * carbon_mass_fraction_mg_kg * 1e-9 * mass_CO2 / mass_C]
cdr_mc[, sd(feedstock_mass_total_kg * carbon_mass_fraction_mg_kg * 1e-9 * mass_CO2 / mass_C)] # SD

# C_feedstock (CDR scenario):

input_isometric_scalar[, (feedstock_mass_treatment_kg - feedstock_mass_kg_counterfactual) * carbon_mass_fraction_mg_kg * 1e-9 * mass_CO2 / mass_C]
cdr_mc[, sd(feedstock_dose_kg * carbon_mass_fraction_mg_kg * 1e-9 * mass_CO2 / mass_C)] # SD

# C_feedstock (BAU liming):

input_isometric_scalar[, feedstock_mass_kg_counterfactual * carbon_mass_fraction_mg_kg * 1e-9 * mass_CO2 / mass_C]
cdr_mc[, sd(feedstock_mass_baseline_kg * carbon_mass_fraction_mg_kg * 1e-9 * mass_CO2 / mass_C)] # SD

# CDR (total):

input_isometric_scalar[
  ,
  1e-9 * mass_CO2 *
    (ocean_retention_treatment * input_isometric_vector[, sum(dic_treatment_mmol_kg * water_density_kg_m3 * river_flow_m3_s * time_interval_s)] -
       ocean_retention_baseline * input_isometric_vector[, sum(dic_counterfactual_mmol_kg * water_density_kg_m3 * river_flow_m3_s * time_interval_s)] -
       (feedstock_mass_treatment_kg - feedstock_mass_kg_counterfactual) * carbon_mass_fraction_mg_kg / mass_C) +
    1e-3 * feedstock_mass_kg_counterfactual / dose_cdr_mass_ratio
]

cdr_mc[, sd(
  1e-9 * mass_CO2 * (ocean_retention_treatment * river_dic_export_treatment_mmol - ocean_retention_baseline * river_dic_export_baseline_mmol - feedstock_dose_kg * carbon_mass_fraction_mg_kg / mass_C) +
    1e-3 * feedstock_mass_baseline_kg / dose_cdr_mass_ratio
)]

# CDR (project):

input_isometric_scalar[
  ,
  1e-9 * mass_CO2 *
    (ocean_retention_treatment * input_isometric_vector[, sum(dic_treatment_mmol_kg * water_density_kg_m3 * river_flow_m3_s * time_interval_s)] -
       ocean_retention_baseline * input_isometric_vector[, sum(dic_counterfactual_mmol_kg * water_density_kg_m3 * river_flow_m3_s * time_interval_s)] -
       (feedstock_mass_treatment_kg - feedstock_mass_kg_counterfactual) * carbon_mass_fraction_mg_kg / mass_C)
]

cdr_mc[, sd(ocean_retention_treatment * river_dic_export_treatment_mmol - ocean_retention_baseline * river_dic_export_baseline_mmol - feedstock_dose_kg * carbon_mass_fraction_mg_kg / mass_C) * 1e-9 * mass_CO2]

# CDR (BAU):

input_isometric_scalar[, 1e-3 * feedstock_mass_kg_counterfactual / dose_cdr_mass_ratio]
cdr_mc[, 1e-3 * sd(feedstock_mass_baseline_kg) / dose_cdr_mass_ratio]

# gross efficiency ratio (total):

input_isometric_scalar[
  ,
  (1e-9 * mass_CO2 *
     (ocean_retention_treatment * input_isometric_vector[, sum(dic_treatment_mmol_kg * water_density_kg_m3 * river_flow_m3_s * time_interval_s)] -
        ocean_retention_baseline * input_isometric_vector[, sum(dic_counterfactual_mmol_kg * water_density_kg_m3 * river_flow_m3_s * time_interval_s)] -
        (feedstock_mass_treatment_kg - feedstock_mass_kg_counterfactual) * carbon_mass_fraction_mg_kg / mass_C) +
     1e-3 * feedstock_mass_kg_counterfactual / dose_cdr_mass_ratio) /
    (1e-3 * feedstock_mass_treatment_kg)
]

cdr_mc[, sd(
  (1e-9 * mass_CO2 * (ocean_retention_treatment * river_dic_export_treatment_mmol - ocean_retention_baseline * river_dic_export_baseline_mmol - feedstock_dose_kg * carbon_mass_fraction_mg_kg / mass_C) +
     1e-3 * feedstock_mass_baseline_kg / dose_cdr_mass_ratio) /
    (1e-3 * feedstock_mass_total_kg)
)]

# gross efficiency ratio (project):

input_isometric_scalar[
  ,
  (1e-9 * mass_CO2 * (ocean_retention_treatment * input_isometric_vector[, sum(dic_treatment_mmol_kg * water_density_kg_m3 * river_flow_m3_s * time_interval_s)] -
                        ocean_retention_baseline * input_isometric_vector[, sum(dic_counterfactual_mmol_kg * water_density_kg_m3 * river_flow_m3_s * time_interval_s)] -
                        (feedstock_mass_treatment_kg - feedstock_mass_kg_counterfactual) * carbon_mass_fraction_mg_kg / mass_C)) /
    (1e-3 * (feedstock_mass_treatment_kg - feedstock_mass_kg_counterfactual))
]

cdr_mc[, sd(((ocean_retention_treatment * river_dic_export_treatment_mmol - ocean_retention_baseline * river_dic_export_baseline_mmol - feedstock_dose_kg * carbon_mass_fraction_mg_kg / mass_C) * 1e-9 * mass_CO2) / (1e-3 * feedstock_dose_kg))]

# gross efficiency ratio (BAU): zero uncertainty because ratio is fixed

input_isometric_scalar[, (1e-3 * feedstock_mass_kg_counterfactual / dose_cdr_mass_ratio) / (1e-3 * feedstock_mass_kg_counterfactual)]
cdr_mc[, sd((1e-3 * feedstock_mass_baseline_kg / dose_cdr_mass_ratio) / (1e-3 * feedstock_mass_baseline_kg))]

# net efficiency ratio (total):

input_isometric_scalar[
  ,
  (1e-9 * mass_CO2 *
     (ocean_retention_treatment * input_isometric_vector[, sum(dic_treatment_mmol_kg * water_density_kg_m3 * river_flow_m3_s * time_interval_s)] -
        ocean_retention_baseline * input_isometric_vector[, sum(dic_counterfactual_mmol_kg * water_density_kg_m3 * river_flow_m3_s * time_interval_s)] -
        (feedstock_mass_treatment_kg - feedstock_mass_kg_counterfactual) * carbon_mass_fraction_mg_kg / mass_C) +
     1e-3 * feedstock_mass_kg_counterfactual / dose_cdr_mass_ratio - 42.4 - 118.3) /
    (1e-3 * feedstock_mass_treatment_kg)
]

cdr_mc[, sd(
  (1e-9 * mass_CO2 * (ocean_retention_treatment * river_dic_export_treatment_mmol - ocean_retention_baseline * river_dic_export_baseline_mmol - feedstock_dose_kg * carbon_mass_fraction_mg_kg / mass_C) +
     1e-3 * feedstock_mass_baseline_kg / dose_cdr_mass_ratio - 42.4 - 118.3) /
    (1e-3 * feedstock_mass_total_kg)
)]

# gross efficiency ratio (project):

input_isometric_scalar[
  ,
  (1e-9 * mass_CO2 * (ocean_retention_treatment * input_isometric_vector[, sum(dic_treatment_mmol_kg * water_density_kg_m3 * river_flow_m3_s * time_interval_s)] -
                        ocean_retention_baseline * input_isometric_vector[, sum(dic_counterfactual_mmol_kg * water_density_kg_m3 * river_flow_m3_s * time_interval_s)] -
                        (feedstock_mass_treatment_kg - feedstock_mass_kg_counterfactual) * carbon_mass_fraction_mg_kg / mass_C) - 35.9 - 100.1) /
    (1e-3 * (feedstock_mass_treatment_kg - feedstock_mass_kg_counterfactual))
]

cdr_mc[, sd(((ocean_retention_treatment * river_dic_export_treatment_mmol - ocean_retention_baseline * river_dic_export_baseline_mmol - feedstock_dose_kg * carbon_mass_fraction_mg_kg / mass_C) * 1e-9 * mass_CO2 - 35.9 - 100.1) / (1e-3 * feedstock_dose_kg))]

# gross efficiency ratio (BAU): zero uncertainty because ratio is fixed

input_isometric_scalar[, (1e-3 * feedstock_mass_kg_counterfactual / dose_cdr_mass_ratio - 6.52 - 18.2) / (1e-3 * feedstock_mass_kg_counterfactual)]
cdr_mc[, sd((1e-3 * feedstock_mass_baseline_kg / dose_cdr_mass_ratio - 6.52 - 18.2) / (1e-3 * feedstock_mass_baseline_kg))]

# sense check on the BAU alkalinity export: -------------------------------

# first, alkalinity (upstream baseline):

alkalinity_upstream <- data_grab_2025[
  station_code == "St. 1",
  .(
    date = as.Date(sample_date, tz = "Europe/Oslo"),
    alkalinity_mmol_l = as.numeric(alk_filt_mmol_l)
  )
]

alkalinity_upstream_interpolated <- merge(
  data.table(date = seq(as.Date(reporting_period_start, tz = "Europe/Oslo"), as.Date(reporting_period_end, tz = "Europe/Oslo"), by = "1 day")),
  alkalinity_upstream,
  by = "date", all.x = TRUE
)[order(date)]

alkalinity_upstream_interpolated[, delta_time_s := as.numeric(difftime(date, shift(date), units = "secs"))]
alkalinity_upstream_interpolated[, delta_time_s := c(delta_time_s[2], delta_time_s[-1])]

alkalinity_upstream_interpolated[, alkalinity_mmol_l := imputeTS::na_interpolation(alkalinity_mmol_l)]

alkalinity_export_upstream <- merge(alkalinity_upstream_interpolated, data_discharge, by = "date")[, sum(alkalinity_mmol_l * river_flow_m3_s * delta_time_s) * mass_CO2 * 1e-6]

isometric_co2e_stored[, signif(river_alk_export_baseline_mmol * 1e-9 * mass_CO2, 3)] - alkalinity_export_upstream

# figure 4 ----------------------------------------------------------------

make_color_ramp <- colorRampPalette(c("grey80", "grey10"))

p_historical_input <- data_dose_combined[
  order(Timestamp_Local),
  .(
    Timestamp_Local,
    value = cumsum(nafill(1e-6 * dose_ml_m3 * discharge_m3_s * delta_time_s, fill = 0)),
    param = "Cumulative slurry dose (m<sup>3</sup>)"
  ),
  by = .(year = year(Timestamp_Local))
]

p_historical_input <- rbind(
  data_doser[, .(Timestamp_Local, `pH Kvitla`)],
  data_doser_historical[Timestamp_Local < data_doser[, min(Timestamp_Local)], .(Timestamp_Local, `pH Kvitla`)],
  fill = TRUE
)[
  `pH Kvitla` > 4,
  .(Timestamp_Local = mean(Timestamp_Local), value = mean(`pH Kvitla`, na.rm = TRUE), param = "pH at Kvitla"),
  by = .(year = as.factor(year(Timestamp_Local)), week = week(Timestamp_Local))
] |>
  rbind(p_historical_input, fill = TRUE)

p_historical_input |>
  ggplot(aes(yday(Timestamp_Local), value, col = year, linewidth = year)) +
  facet_wrap(vars(param), ncol = 1, scales = "free_y") +
  scale_color_manual(values = c(make_color_ramp(10), "red")) +
  scale_linewidth_manual(values = c(rep(0.3, 10), 0.6)) +
  geom_vline(xintercept = yday(as.Date("2025-05-15")), linetype = 3) +
  geom_label(
    data = data.table(
      param = "Cumulative slurry dose (m<sup>3</sup>)",
      x = yday(as.Date("2025-05-15")) + 28,
      y = 1700,
      label = "CarbonRun RAE\ntrial begins (2025)"
    ),
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    hjust = "right", linewidth = 0,
    size = 2, alpha = 1, label.padding = unit(0.5, "lines")
  ) +
  theme(strip.text = element_markdown()) +
  geom_line(data = \(x) x[!year %in% c(2015, 2022) | param == "pH at Kvitla"]) +
  guides(linewidth = "none") +
  scale_x_date(date_labels = "%b") +
  labs(x = NULL, y = NULL, col = NULL)

ggsave("figures/manuscript/figure-4.png", dpi = 600, width = 3.33, height = 4)

# figure 5 ----------------------------------------------------------------

melt(data_grab_2025[station_code %in% c("St. 1", "St. 2")], measure.vars = c("turbiditet_fnu", "al_filt_ug_l", "fe_filt_ug_l"), variable.name = "param")[
  , .(date = as.Date(sample_date), param, station_code, value)
] |>
  rbind(data_co2_removal_interpolated[, .(date, param = "flow", value = river_flow_m3_s, station_code = "flow")], fill = TRUE) |>
  ggplot(aes(date, value, col = station_code)) +
  facet_wrap(
    vars(param),
    scales = "free_y", ncol = 1,
    labeller = as_labeller(c(
      "flow" = "Discharge (m<sup>3</sup> s<sup>-1</sup>)",
      "turbiditet_fnu" = "Turbidity (FNU)",
      "al_filt_ug_l" = "[Al]<sub>dissolved</sub> (&mu;g L<sup>-1</sup>)",
      "fe_filt_ug_l" = "[Fe]<sub>dissolved</sub> (&mu;g L<sup>-1</sup>)"
    ))
  ) +
  geom_line(data = \(x) x[date <= data_co2_removal_interpolated[, max(date)]]) +
  scale_color_manual(
    breaks = c("St. 1", "St. 2"),
    labels = c("Site NY", "Site OB"),
    values = c("St. 1" = unname(site_colours["NY"]), "St. 2" = unname(site_colours["OB"]), "flow" = "grey50")
  ) +
  theme(strip.text = element_markdown()) +
  labs(x = NULL, y = NULL, col = NULL)

ggsave("figures/manuscript/figure-5.png", width = 3.33, height = 4.5, dpi = 600)

# supplementary figure ----------------------------------------------------

n_replicates <- 1000

specific_surface_area_cm2_mol <- runif(n_replicates, 1e6, 5e6)
saturation_inhibition <- runif(n_replicates, 0.5, 1)

model_ensemble <- future_map2(specific_surface_area_cm2_mol, saturation_inhibition, \(specific_surface_area, exponent) {
  phr_kinetic_model_modified <- sub("parms    4e\\+06    1", paste0("parms    ", specific_surface_area, "    ", exponent), phr_kinetic_model)
  phreeqc::phrLoadDatabaseString(phreeqc::phreeqc.dat)
  phreeqc::phrRunString(phr_kinetic_model_modified)
  lapply(.Call("getSelOutLst", PACKAGE = "phreeqc"), as.data.table)$n1
}, .progress = TRUE) |>
  rbindlist(idcol = "rep")

p_ensemble <- model_ensemble[
  , .(pH = mean(pH), pH_lower = quantile(pH, 0.025), pH_upper = quantile(pH, 0.975)),
  by = .(time_h = time / 3600)
] |>
  ggplot(aes(time_h, pH)) +
  geom_ribbon(aes(ymin = pH_lower, ymax = pH_upper), alpha = 0.4) +
  geom_line(data = model_ensemble[rep %in% 1:100], aes(time / 3600, pH, group = rep), linewidth = 0.1) +
  geom_line(col = "white") +
  labs(x = "Time (h)") +
  guides(col = "none") +
  theme(plot.tag = element_text(face = "bold")) +
  lims(x = c(0, NA), y = c(5.4, NA))

sensitivity_input <- expand.grid(
  specific_surface_area_cm2_mol = seq(1e6, 5e6, length.out = 10),
  saturation_inhibition = seq(0.5, 1, length.out = 10)
)

model_sensitivity <- future_map2(sensitivity_input$specific_surface_area_cm2_mol, sensitivity_input$saturation_inhibition, \(specific_surface_area, exponent) {
  phr_kinetic_model_modified <- sub("parms    4e\\+06    1", paste0("parms    ", specific_surface_area, "    ", exponent), phr_kinetic_model)
  phreeqc::phrLoadDatabaseString(phreeqc::phreeqc.dat)
  phreeqc::phrRunString(phr_kinetic_model_modified)
  lapply(.Call("getSelOutLst", PACKAGE = "phreeqc"), as.data.table)$n1
}, .progress = TRUE) |>
  rbindlist(idcol = "rep")

initial_dose_mol_kgw <- as.numeric(sub("-m0\\s+([^\\s][0-9.eE+-]+).*", "\\1", phr_kinetic_model[grepl("-m0", phr_kinetic_model)]))
p_sensitivity <- model_sensitivity[, .(p_Calcite = 1e2 * (1 - k_Calcite[which.max(step)] / initial_dose_mol_kgw)), by = "rep"] |>
  ggplot(aes(sensitivity_input$specific_surface_area_cm2_mol * 1e-4, sensitivity_input$saturation_inhibition)) +
  geom_raster(aes(fill = p_Calcite)) +
  scale_fill_viridis_c() +
  theme(
    axis.title.x = element_markdown(),
    axis.title.y = element_markdown(),
    legend.position = "bottom",
    plot.tag = element_text(face = "bold"),
    margins = margin(r = 10)
  ) +
  labs(x = "Specific surface area (m<sup>2</sup> mol<sup>-1</sup>)", y = "Surface area<br>exponent &gamma;", fill = "Calcite\n(% dissolved)")

wrap_plots(p_ensemble, p_sensitivity, ncol = 1) +
  plot_annotation(tag_level = "a")

ggsave("figures/supporting-information/kinetic-dissolution.png", width = 3, height = 5, dpi = 600)

# supplementary figure ----------------------------------------------------

p_chem <- data_historical[Vannlokalitetsnavn %in% c("Litlåne ved Åmot", "Kvina ved Klosterøyna") & !is.na(Verdi)][
  , `:=`(
    date = as.Date(Tid_provetak, tz = "UTC"),
    site = fifelse(Vannlokalitetsnavn == "Kvina ved Klosterøyna", "KL", "LA"),
    param = fcase(
      grepl("alkalitet", Parameter_navn), "alkalinity_mmol_l",
      grepl("alum", Parameter_navn), "Al",
      default = Parameter_navn
    )
  )
][param != "Kalium"] |>
  dcast(date + param ~ site, value.var = "Verdi", fun.aggregate = mean) |>
  na.omit() |>
  rbind(dcast(data_grab_historical, date + param ~ site, value.var = "Verdi")) |> # with(unique(param))
  ggplot(aes(LA, KL)) +
  facet_wrap(vars(param), scales = "free", labeller = as_labeller(c(
    "alkalinity_mmol_l" = "Total alkalinity<br>(mmol L<sup>-1</sup>)",
    "pH" = "pH",
    "Al" = "Al (&mu;g L<sup>-1</sup>)",
    "Kalsium" = "Ca (mg L<sup>-1</sup>)",
    "Konduktivitet" = "Conductivity<br>(mS m<sup>-1</sup>)",
    "Totalt organisk karbon (TOC)" = "Total organic<br>carbon (mg L<sup>-1</sup>)"
  )), ncol = 2) +
  theme(strip.text = element_markdown()) +
  geom_abline() +
  geom_point(shape = 16, alpha = 0.6) +
  labs(x = "LA (limed control)", y = "KL (river mouth)")

p_flow <- merge(data_roynestad, data_krageholen, by = "Time", all = TRUE) |>
  ggplot(aes(discharge_krageholen_m3_s, discharge_roynestad_m3_s)) +
  theme(
    axis.title.x = element_markdown(),
    axis.title.y = element_markdown()
  ) +
  scale_x_log10() +
  scale_y_log10() +
  geom_abline() +
  geom_point(shape = 16, alpha = 0.6) +
  labs(
    x = "Discharge at Kr&aring;geh&oslash;len (m<sup>3</sup> s<sup>-1</sup>)",
    y = "Discharge at<br>R&oslash;ynestad (m<sup>3</sup> s<sup>-1</sup>)"
  )

wrap_plots(p_chem, p_flow, ncol = 1, heights = c(2.5, 1)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold"))

ggsave("figures/supporting-information/baseline-data.png", width = 3.33, height = 6, dpi = 600)

# ocean retention factor --------------------------------------------------
#
# Seasonal climatology of the Renforth-Henderson ocean retention factor,
# salinity, pCO2, and temperature, built from the QC-filtered ICOS North Sea
# underway data already loaded and retention-factor-scored in R/01-run-model.R.

rh_plot_input <- melt(data_north_sea_qaqc, measure.vars = c("rh_factor", "p_sal_psu", "p_co2_uatm", "temp_deg_c"))[
  , .(
    mean = mean(value, na.rm = TRUE),
    sd = sd(value, na.rm = TRUE)
  ),
  by = .(variable, doy = yday(date_time))
]

merge(rh_plot_input, expand.grid(variable = c("rh_factor", "p_sal_psu", "p_co2_uatm", "temp_deg_c"), doy = 1:366), all.y = TRUE)[
  , date_plot := as.Date(paste("2024", doy), format = "%Y %j")
][] |>
  ggplot(aes(date_plot, mean)) +
  facet_wrap(
    vars(variable),
    nrow = 1, scales = "free_y",
    labeller = as_labeller(c(
      rh_factor = "Ocean retention factor",
      p_sal_psu = "Salinity (PSU)",
      p_co2_uatm = "pCO<sub>2</sub> (&mu;atm)",
      temp_deg_c = "Temp. (&deg;C)"
    ))
  ) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0, col = "grey") +
  geom_line() +
  theme(strip.text = element_markdown()) +
  scale_x_date(date_labels = "%b") +
  labs(x = NULL, y = NULL)

ggsave("figures/supporting-information/ocean-retention-data.png", width = 7, height = 1.8, dpi = 600)

# extras ------------------------------------------------------------------

# get cumulative observed vs potential CDR:

plot_input[, .(date, ratio = cumsum(co2e_stored) / cumsum(co2e_stored_upper_bound), river_flow_m3_s)]

# extrapolate CDR to the entire year using discharge weighting

t_co2e_m3 <- plot_input[, .(t_co2e_m3 = sum(co2e_stored) / sum(river_flow_m3_s * 3600 * 24))]

stopifnot("data_dose_combined is not in chronological order" = data_dose_combined[, all(diff(Timestamp_Local) > 0)])

data_dose_combined[
  , .(
    m3_yr = sum((discharge_m3_s * discharge_scaling_factor - withdrawals_m3_s) * delta_time_s),
    days = uniqueN(yday(Timestamp_Local)),
    date_min = min(as.Date(Timestamp_Local, tz = "Europe/Oslo")),
    date_max = max(as.Date(Timestamp_Local, tz = "Europe/Oslo"))
  ),
  by = .(year = year(Timestamp_Local))
][
  days > 364
][
  , m3_yr_fraction := m3_yr / plot_input[, sum(river_flow_m3_s * 3600 * 24)]
][
  , t_co2e_yr := m3_yr * t_co2e_m3[, t_co2e_m3]
][
  , t_co2e_yr_NET := m3_yr * plot_input[, (76.69 + input_isometric_scalar$uncertainty_discount_t) / sum(river_flow_m3_s * 3600 * 24)]
][
  , .(t_co2e_yr = mean(t_co2e_yr), t_co2e_yr_NET = mean(t_co2e_yr_NET))
]

# calculate saturation indices (historical data):

calcite_saturation_input <- data_historical[Vannlokalitetsnavn == "Kvina ved Klosterøyna" & !is.na(Verdi)][
  , param_name := paste(Parameter_navn, Enhet)
] |>
  dcast(Tid_provetak ~ param_name, value.var = "Verdi", fun.aggregate = median)

calcite_saturation_input <- na.omit(calcite_saturation_input[, c(
  "Tid_provetak", "Kalium mg/l", "Kalsium mg/l", "Klorid mg/l", "Magnesium mg/l", "Natrium mg/l", "Reaktivt aluminium µg/l Al",
  "Sulfat mg/l", "Total alkalitet mmol/l", "Totalfosfor µg/l P", "Totalnitrogen µg/l N", "pH <ubenevnt>"
)])

phr_calcite_saturation <- lapply(seq(nrow(calcite_saturation_input)), \(i) {
  with(calcite_saturation_input[i], glue_data(
    list(
      number = i,
      pH = `pH <ubenevnt>`,
      Alkalinity = paste(`Total alkalitet mmol/l`, "mmol/l"),
      # major cations:
      K = `Kalium mg/l`,
      Ca = `Kalsium mg/l`,
      Mg = `Magnesium mg/l`,
      Na = `Natrium mg/l`,
      # major anions:
      Cl = `Klorid mg/l`,
      S = paste(`Sulfat mg/l`, "as SO4"),
      P = paste(`Totalfosfor µg/l P`, "ug/l as P"),
      N = paste(`Totalnitrogen µg/l N`, "ug/l as N"),
      # metals:
      Al = paste(`Reaktivt aluminium µg/l Al`, "ug/l")
    ),
    paste(phr_calcite_saturation_template, collapse = "\n")
  ))
}) |>
  c("SELECTED_OUTPUT\n    saturation_indices    Calcite") |>
  paste(collapse = "\n")

phreeqc::phrLoadDatabaseString(phreeqc::phreeqc.dat)
phreeqc::phrRunString(phr_calcite_saturation)
calcite_saturation_output <- lapply(.Call("getSelOutLst", PACKAGE = "phreeqc"), as.data.table)$n1

calcite_saturation_output[
  !is.na(si_Calcite),
  unlist(lapply(
    .SD,
    \(x) {
      list(min = min(x), max = max(x), median = median(x), mean = mean(x), sd = sd(x))
    }
  )),
  .SDcols = "si_Calcite"
] |>
  lapply(\(x) signif(x, 3))

# calculate saturation indices (during trial):

calcite_saturation_input_trial <- data_grab_2025[station_code == "St. 4b" & as.Date(sample_date) %between% c("2025-06-10", "2025-09-15")]

phr_calcite_saturation_trial <- lapply(seq(nrow(calcite_saturation_input_trial)), \(i) {
  with(calcite_saturation_input_trial[i], glue_data(
    list(
      number = i,
      pH = ph,
      Alkalinity = paste(1e-3 * anc_u_ekv_l, "mmol/l"),
      # major cations:
      K = k_filt_mg_l,
      Ca = ca_filt_mg_l,
      Mg = mg_filt_mg_l,
      Na = na_filt_mg_l,
      # major anions:
      Cl = cl_mg_l,
      S = paste(so4_mg_l, "as SO4"),
      P = paste(totp_ug_p_l, "ug/l"),
      N = paste(no3_n_ug_l, "ug/l as N"),
      # metals:
      Al = paste(al_filt_ug_l, "ug/l")
    ),
    paste(phr_calcite_saturation_template, collapse = "\n")
  ))
}) |>
  c("SELECTED_OUTPUT\n    saturation_indices    Calcite") |>
  paste(collapse = "\n")

phreeqc::phrLoadDatabaseString(phreeqc::phreeqc.dat)
phreeqc::phrRunString(phr_calcite_saturation_trial)
calcite_saturation_output_trial <- lapply(.Call("getSelOutLst", PACKAGE = "phreeqc"), as.data.table)$n1

calcite_saturation_output_trial[
  !is.na(si_Calcite),
  unlist(lapply(
    .SD,
    \(x) {
      list(min = min(x), max = max(x), median = median(x), mean = mean(x), sd = sd(x))
    }
  )),
  .SDcols = "si_Calcite"
] |>
  lapply(\(x) signif(x, 3))

# sense check on counterfactual vs additional dose (simplistic carbonate-system-only case):

phreeqc::phrLoadDatabaseString(phreeqc::phreeqc.dat)
phreeqc::phrRunString(phr_calcite_addition)
counterfactual <- lapply(.Call("getSelOutLst", PACKAGE = "phreeqc"), as.data.table)$n1

phreeqc::phrRunString(sub("Fix_pH    -6", "Fix_pH    -7", phr_calcite_addition))
observed <- lapply(.Call("getSelOutLst", PACKAGE = "phreeqc"), as.data.table)$n1

(-observed$d_Fix_pH[2] - -counterfactual$d_Fix_pH[2]) / -counterfactual$d_Fix_pH[2]

# Ca mass balance during the trial:

calcium_concentrations <- data_grab_2025[
  !is.na(station_name) & !is.na(sample_date)
][
  station_code %in% c("St. 1", "St. 4b", "St. 8")
][
  , `:=`(
    date = as.Date(sample_date, tz = "Europe/Oslo"),
    station_code = fcase(station_code == "St. 1", "NY", station_code == "St. 4b", "KL", station_code == "St. 8", "LA")
  )
] |>
  dcast(date ~ station_code, value.var = "ca_mg_l")

calcium_exports <- merge(data_co2_removal_interpolated, calcium_concentrations, by = "date", all.x = TRUE)[
  , `:=`(
    calcium_mg_l_KL = na_interpolation(KL),
    calcium_mg_l_NY = na_interpolation(NY),
    calcium_mg_l_LA = na_interpolation(LA)
  )
][
  , `:=`(
    calcium_t_KL = 1e-6 * calcium_mg_l_KL * (river_flow_m3_s + withdrawals_m3_s) * delta_time_s,
    calcium_t_NY = 1e-6 * calcium_mg_l_NY * (river_flow_m3_s + withdrawals_m3_s) * delta_time_s / discharge_scaling_factor
  )
] |>
  merge(data_krageholen[, date := as.Date(Time)][date %between% c("2025-06-10", "2025-09-15")], by = "date", all.x = TRUE)

calcium_exports[, calcium_t_LA := 1e-6 * calcium_mg_l_LA * discharge_krageholen_m3_s * delta_time_s][
  , lapply(.SD, sum),
  .SDcols = patterns("^calcium_t_")
][
  # convert to CaCO3:
  , (calcium_t_KL - calcium_t_NY - calcium_t_LA) * mass_CaCO3 / mass_Ca / calcite_wt_pct
]

data_doser_aggregated[date >= "2025-06-10", lapply(.SD, sum), .SDcols = is.numeric] # compare to dosing records
