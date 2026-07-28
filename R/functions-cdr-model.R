
# R/functions-cdr-model.R -- Functions to support CDR quantification according to the Isometric 
# methodology: https://registry.isometric.com/protocol/river-alkalinity-enhancement/1.1
# Sourced with working directory = repo root.

# functions ---------------------------------------------------------------

# calculate error metrics:

get_error <- function(data, response, predictions, ...) {
  error <- data[[response]] - predictions
  list(
    RMSE = sqrt(mean(error ^ 2, na.rm = TRUE)),
    Bias = mean(error, na.rm = TRUE),
    `R<sup>2</sup>` = r_sq(data[[response]], predictions)
  )
}

# r-squared:

r_sq <- function (observations, predictions, mu = NULL) 
{
  if (!identical(length(observations), length(predictions))) 
    stop("length mismatch between 'observations' and 'predictions'")
  if (identical(length(observations), 1L)) 
    stop("r-squared is not defined for N=1")
  if (length(mu) > 1) 
    stop("mu must be scalar")
  if (!(is.numeric(observations) & is.numeric(predictions))) 
    stop("'observations' and 'predictions' must be numeric")
  data <- na.omit(data.frame(observations, predictions))
  if (identical(nrow(data), 0L)) 
    stop("no complete cases")
  if (is.null(mu)) 
    mu <- mean(data$observations)
  1 - sum((data$observations - data$predictions)^2)/sum((data$observations - mu)^2)
}

# calculate renforth-henderson:

calculate_rh_factor <- function (salinity, pco2, temperature) {
  (salinity * 10^-3.009 + 10^-1.519) * log(pco2) - 
    salinity * 10^-2.1 - 
    temperature * pco2 * (salinity * 10^-7.501 - 10^-5.598) - 
    temperature * 10^-2.337 + 
    10^-0.102
}

# cross-validation:

loocv <- function (data, response) {
  n <- nrow(data)
  these_rows <- seq(n)
  predictions <- vector(mode = "numeric", length = n)
  this_formula <- stats::as.formula(paste0(response, " ~ ."))
  for (row in these_rows) {
    data_train <- data[-row, ]
    data_test <- data[row, ]
    data_train <- mean_imputation(data_train)
    data_test <- mean_imputation(data_train, data_test)
    model <- stats::lm(this_formula, data = data_train)
    predictions[row] <- predict(model, data_test)
  }
  predictions
}

mean_imputation <- function (data_train, data_test = NULL) {
  column_means <- colMeans(data_train, na.rm = TRUE)
  column_means_names <- setNames(names(column_means), names(column_means))
  data_out <- if (is.null(data_test)) 
    data_train
  else data_test
  data.table::as.data.table(lapply(column_means_names, function(x) {
    data_out[[x]] <- data.table::fifelse(is.na(data_out[[x]]), column_means[[x]], data_out[[x]])
    data_out[[x]]
  }))
}

# calculate water density (Formula source: doi:10.6028/jres.097.013):

rho <- function (temperature) {
  a <- 999.83952
  b <- 16.945176
  c <- -0.00798704
  d <- -4.6170461e-05
  e <- 1.0556302e-07
  f <- -2.8054253e-10
  g <- 0.01689785
  0.001 * (a + b * temperature + c * temperature^2 + d * temperature^3 + 
             e * temperature^4 + f * temperature^5)/(1 + g * temperature)
}

# calculate co2e stored:

co2e_stored <- function(
    feedstock_mass_kg, 
    carbon_mass_fraction_mg_kg, 
    ocean_retention, 
    dic_treatment_mmol_kg, 
    dic_baseline_mmol_kg,
    alk_treatment_mmol_kg,
    alk_baseline_mmol_kg,
    water_density_kg_m3,
    river_flow_downstream_m3_s,
    time_interval_s,
    C = mass_C,
    CO2 = mass_CO2) {
  
  water_mass_kg <- water_density_kg_m3 * river_flow_downstream_m3_s * time_interval_s
  
  river_dic_export_treatment_mmol <- sum(dic_treatment_mmol_kg * water_mass_kg)
  river_dic_export_baseline_mmol <- sum(dic_baseline_mmol_kg * water_mass_kg)
  
  river_alk_export_treatment_mmol <- sum(alk_treatment_mmol_kg * water_mass_kg)
  river_alk_export_baseline_mmol <- sum(alk_baseline_mmol_kg * water_mass_kg)
  
  ocean_retention <- data.table(
    treatment = pmin(ocean_retention * river_alk_export_treatment_mmol / river_dic_export_treatment_mmol, 1),
    baseline = ocean_retention * river_alk_export_baseline_mmol / river_dic_export_baseline_mmol
  )
  
  co2e_feedstock_t <- (CO2 / C) * feedstock_mass_kg * carbon_mass_fraction_mg_kg * 1e-9
  river_dic_export_treatment_t <- CO2 * ocean_retention$treatment * river_dic_export_treatment_mmol * 1e-9
  river_dic_export_baseline_t <- CO2 * ocean_retention$baseline * river_dic_export_baseline_mmol * 1e-9
  
  data.table(
    co2e_stored_t = river_dic_export_treatment_t - river_dic_export_baseline_t - co2e_feedstock_t,
    feedstock_dose_kg = feedstock_mass_kg,
    ocean_retention_treatment = ocean_retention$treatment,
    ocean_retention_baseline = ocean_retention$baseline,
    carbon_mass_fraction_mg_kg = carbon_mass_fraction_mg_kg,
    river_dic_export_treatment_mmol,
    river_dic_export_baseline_mmol,
    river_alk_export_treatment_mmol,
    river_alk_export_baseline_mmol,
    water_mass_kg = sum(water_mass_kg)
  )
}