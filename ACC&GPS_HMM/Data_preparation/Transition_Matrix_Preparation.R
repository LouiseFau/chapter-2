#' -----------------------------------------------------------------------------
#' Title: Aerial-origin transition datasets for two GPS and two ACC datasets ----
#' Author: Louise Faure
#' Date: 03.08.26, cleaned on 25.11.26
#'
#' **Info:** follows 'ACC_data_preparation' and 'GPS_data_preparation'.
#' **Output:** one standardized dataset of aerial-origin transitions per sampling
#' resolution, ready for the GLMMs fitted in Covariates_HFI_Selection.R.
#'
#' **Main steps:**
#' (1) compute age since emigration (ACC only);
#' (2) prepare behavioural states and elapsed state durations;
#' (3) build consecutive transitions at the expected time step;
#' (4) retain transitions originating from the aerial state;
#' (5) remove individuals with fewer than 30 aerial transitions, then
#'     incomplete observations;
#' (6) standardize covariates within each dataset;
#' (7) export.
#' -----------------------------------------------------------------------------

# Libraries ----
library(dplyr)
library(tidyr)
library(tibble)
library(sf)

# Data obtained from ACC and GPS preparation scripts ----
input_datasets <- list(
  GPS_20 = readRDS("/Users/louisefaure/Desktop/dossier sans titre/donnees filtree/GE_gps_20_covariates_hfi.rds"),
  GPS_60 = readRDS("/Users/louisefaure/Desktop/dossier sans titre/donnees filtree/GE_gps_60_covariates_hfi.rds"),
  ACC_20 = readRDS("/Users/louisefaure/Desktop/dossier sans titre/donnees filtree/GE_acc_20_covariates_hfi.rds"),
  ACC_60 = readRDS("/Users/louisefaure/Desktop/dossier sans titre/donnees filtree/GE_acc_60_covariates_hfi.rds"))

emig_dates_raw <- readRDS("/Users/louisefaure/Library/CloudStorage/OneDrive-Personnel/THESE/CHAPITRE 2/git/chapter-2/DONNEES AIGLES/emigration dates/emigration_dates_20250417.rds")

# Parameters ----
state_levels <- list(
  GPS = c("aerial","terrestrial"),
  ACC = c("aerial","resting","feeding"))

dataset_parameters <- tibble::tribble(
  ~dataset, ~data_type, ~expected_dt_min, ~dt_tolerance_min, ~output_file,
  "GPS_20", "GPS", 20,  5, "gps_20.rds",
  "GPS_60", "GPS", 60, 10, "gps_60.rds",
  "ACC_20", "ACC", 20,  5, "acc_20.rds",
  "ACC_60", "ACC", 60, 10, "acc_60.rds")

dataset_names     <- dataset_parameters$dataset
acc_dataset_names <- dataset_parameters$dataset[dataset_parameters$data_type == "ACC"]

minimum_aerial_transitions <- 30L
minimum_age_days <- 0
maximum_age_days <- 105

output_directory <- "/Users/louisefaure/Desktop/dossier sans titre/donnees filtree"
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)

required_numeric_variables <- c(
  "age_since_emig_days","age_since_emig_weeks","aerial_duration_min",
  "cos_diel","sin_time","ruggedness_100m","slope_100m",
  "distance_to_ridgeline_100m","elevation_100m","prop_forest_5cells",
  "prop_low_vegetation_5cells","prop_rocky_terrain_5cells",
  "prop_other_5cells","settlement_density","population_density",
  "settlement_density_raw")

environmental_hfi_variables <- setdiff(required_numeric_variables,
  c("age_since_emig_days","age_since_emig_weeks","aerial_duration_min",
    "cos_diel","sin_time"))


# Functions ----
# Emigration dates function ----
#' Clean the emigration date dataset by removing space and keeping one info per individuals
emig_dates <- emig_dates_raw %>%
  dplyr::filter(did_disperse, !is.na(dispersal_date)) %>%
  dplyr::mutate(
    individual.local.identifier = trimws(as.character(individual.local.identifier)),
    dispersal_date = if (inherits(dispersal_date, "POSIXct")) dispersal_date
    else as.POSIXct(as.character(dispersal_date), tz = "UTC")) %>%
  dplyr::select(individual.local.identifier, dispersal_date) %>%
  dplyr::distinct(individual.local.identifier, .keep_all = TRUE)

# Age since emigration ----
#' add age since emigration to each individuals
calculate_acc_age <- function(data) {
  data %>%
    dplyr::select(-dplyr::any_of(c("dispersal_date","age_since_emig_days",
                                   "age_since_emig_weeks"))) %>%
    dplyr::mutate(
      individual.local.identifier = trimws(as.character(individual.local.identifier)),
      timestamp = as.POSIXct(timestamp, tz = "UTC")) %>%
    dplyr::left_join(emig_dates, by = "individual.local.identifier") %>%
    dplyr::mutate(
      age_since_emig_days  = as.numeric(difftime(timestamp, dispersal_date, units = "days")),
      age_since_emig_weeks = age_since_emig_days / 7) %>%
    dplyr::filter(
      age_since_emig_days >= minimum_age_days,
      age_since_emig_days <= maximum_age_days)}

# Behavioural states and elapsed duration within the current state function ----
#' Create a binary state category for ACC and three state category for GPS
#' Calculate the time spent in each state
prepare_states <- function(data, data_type) {
  
  behavior_values <- if (data_type == "GPS") {
    data$behavior_binary %>% as.character() %>% trimws() %>% tolower() %>%
      dplyr::recode("aerian" = "aerial", "flight" = "aerial",
                    "flying" = "aerial", "ground" = "terrestrial")} else {
    data$behavior_reclassified %>% as.character() %>% trimws() %>% tolower() %>%
      dplyr::recode("aerian" = "aerial", "flight" = "aerial",
                    "flying" = "aerial", "rest" = "resting",
                    "feed" = "feeding", "foraging" = "feeding")}
  
  data %>%
    dplyr::mutate(
      timestamp = as.POSIXct(timestamp, tz = "UTC"),
      individual.local.identifier = trimws(as.character(individual.local.identifier)),
      burst_id = as.character(burst_id),
      behavior_state = factor(behavior_values, levels = state_levels[[data_type]])) %>%
    dplyr::arrange(individual.local.identifier, burst_id, timestamp) %>%
    dplyr::group_by(individual.local.identifier, burst_id) %>%
    dplyr::mutate(
      state_bout_n = cumsum(
        dplyr::row_number() == 1L |
          dplyr::coalesce(behavior_state != dplyr::lag(behavior_state), TRUE))) %>%
    dplyr::group_by(individual.local.identifier, burst_id, state_bout_n) %>%
    dplyr::mutate(
      state_duration_min = as.numeric(
        difftime(timestamp, dplyr::first(timestamp), units = "mins"))) %>%
    dplyr::ungroup()}

# Consecutive transitions at the expected time step ----
#' To compute transition, now each line of the dataset correspond to a transition (and not a location)
#' Keep transition that are not too gapy 
construct_transitions <- function(data, expected_dt_min, dt_tolerance_min) {
  data %>%
    dplyr::arrange(individual.local.identifier, burst_id, timestamp) %>%
    dplyr::group_by(individual.local.identifier, burst_id) %>%
    dplyr::mutate(
      behavior_to = dplyr::lead(behavior_state),
      dt_min = as.numeric(difftime(dplyr::lead(timestamp), timestamp, units = "mins"))) %>%
    dplyr::ungroup() %>%
    dplyr::filter(
      !is.na(behavior_state),
      !is.na(behavior_to),
      is.finite(dt_min),
      dt_min > 0,
      abs(dt_min - expected_dt_min) <= dt_tolerance_min)}

# Standardization function ----
standardize_dataset <- function(data, dataset_name) {
  variables_to_standardize <- c("age_since_emig_weeks","aerial_duration_min",
                                environmental_hfi_variables)
  parameters <- tibble::tibble(
    dataset        = dataset_name,
    variable       = variables_to_standardize,
    transformation = "center and standardize",
    center = vapply(variables_to_standardize,
                    function(variable) mean(data[[variable]]), numeric(1)),
    scale  = vapply(variables_to_standardize,
                    function(variable) stats::sd(data[[variable]]), numeric(1)))
  
  if (any(!is.finite(parameters$scale) | parameters$scale <= 0))
    stop("At least one variable has a non-finite or zero SD in ", dataset_name, ".")
  for (i in seq_len(nrow(parameters))) {
    variable <- parameters$variable[[i]]
    data[[paste0(variable, "_z")]] <-
      (data[[variable]] - parameters$center[[i]]) / parameters$scale[[i]]}
  
  cos_center <- mean(data$cos_diel)
  sin_center <- mean(data$sin_time)
  
  data <- data %>%
    dplyr::mutate(
      age_z       = age_since_emig_weeks_z,
      duration_z  = aerial_duration_min_z,
      age_z2      = age_z^2,
      duration_z2 = duration_z^2,
      cos_diel_c  = cos_diel - cos_center,
      sin_diel_c  = sin_time - sin_center)
  
  diel_parameters <- tibble::tibble(
    dataset        = dataset_name,
    variable       = c("cos_diel","sin_time"),
    transformation = "center only",
    center         = c(cos_center, sin_center),
    scale          = 1)
  list(data = data, parameters = dplyr::bind_rows(parameters, diel_parameters))}


#------------------------------------------------------------------------------- STEP 1: loop through the datasets ----
# Information to keep in mind ----
aerial_transitions         <- stats::setNames(vector("list", length(dataset_names)), dataset_names)
standardization_parameters <- stats::setNames(vector("list", length(dataset_names)), dataset_names)
removed_individuals        <- stats::setNames(vector("list", length(dataset_names)), dataset_names)

# 1.0 Loop through the dataset
for (dataset_name in dataset_names) {
  
  parameters_i <- dataset_parameters[dataset_parameters$dataset == dataset_name, ]
  data_type_i  <- parameters_i$data_type
  levels_i     <- state_levels[[data_type_i]]
  
  message("--- ", dataset_name, " ---")
  
  # 1.1 Age since emigration (ACC only) ----
  data_i <- input_datasets[[dataset_name]]
  if (data_type_i == "ACC") data_i <- calculate_acc_age(data_i)
  
  # 1.2 Behavioural states and elapsed state duration ----
  data_i <- prepare_states(data_i, data_type_i)
  
  # 1.3 Transitions at the expected time step ----
  data_i <- construct_transitions(data_i,
                                  parameters_i$expected_dt_min,
                                  parameters_i$dt_tolerance_min)
  
  # 1.4 Transitions originating from the aerial state ----
  aerial_i <- data_i %>%
    dplyr::filter(behavior_state == "aerial") %>%
    dplyr::mutate(
      transition_destination = stats::relevel(
        factor(behavior_to, levels = levels_i), ref = "aerial"),
      remain_aerial       = as.integer(transition_destination == "aerial"),
      aerial_duration_min = state_duration_min)
  
  # 1.5 Aerial transitions per individual, then support filter ----
  support_i <- aerial_i
  if (inherits(support_i, "sf")) support_i <- sf::st_drop_geometry(support_i)
  
  support_i <- support_i %>%
    dplyr::group_by(individual.local.identifier) %>%
    dplyr::summarise(n_aerial_transitions = dplyr::n(), .groups = "drop")
  
  removed_individuals[[dataset_name]] <- support_i %>%
    dplyr::filter(n_aerial_transitions < minimum_aerial_transitions) %>%
    dplyr::mutate(dataset = dataset_name, .before = 1)
  
  retained_i <- support_i %>%
    dplyr::filter(n_aerial_transitions >= minimum_aerial_transitions) %>%
    dplyr::pull(individual.local.identifier)
  
  # 1.6 Complete observations on the modelling variables ----
  missing_variables <- setdiff(required_numeric_variables, names(aerial_i))
  if (length(missing_variables) > 0L)
    stop("Missing numeric variables in ", dataset_name, ": ",
         paste(missing_variables, collapse = ", "))
  
  aerial_i <- aerial_i %>%
    dplyr::filter(individual.local.identifier %in% retained_i) %>%
    dplyr::mutate(individual_id = factor(individual.local.identifier)) %>%
    tidyr::drop_na(individual_id, burst_id, transition_destination, remain_aerial,
                   dplyr::all_of(required_numeric_variables)) %>%
    dplyr::filter(dplyr::if_all(dplyr::all_of(required_numeric_variables), is.finite)) %>%
    dplyr::mutate(individual_id = droplevels(individual_id))
  
  # 1.7 Within-dataset standardization ----
  standardized_i <- standardize_dataset(aerial_i, dataset_name)
  
  aerial_transitions[[dataset_name]]         <- standardized_i$data
  standardization_parameters[[dataset_name]] <- standardized_i$parameters
  
  # 1.8 Export ----
  saveRDS(standardized_i$data,
          file = file.path(output_directory, parameters_i$output_file))
  
  message("    ", nrow(standardized_i$data), " aerial transitions | ",
          dplyr::n_distinct(standardized_i$data$individual_id), " individuals | ",
          nrow(removed_individuals[[dataset_name]]), " individuals removed")
}

standardization_parameters <- dplyr::bind_rows(standardization_parameters)


#------------------------------------------------------------------------------- STEP 2: summaries ----
# 2.1 Individuals removed for insufficient aerial transitions ----
removed_individuals_summary <- dplyr::bind_rows(removed_individuals) %>%
  dplyr::rename(individual = individual.local.identifier) %>%
  dplyr::arrange(dataset, n_aerial_transitions)

# 2.2 Aerial-to-terrestrial transitions, all four datasets ----
#' For ACC, "terrestrial" is any non-aerial destination, i.e. resting + feeding.
aerial_to_terrestrial_summary <- dplyr::bind_rows(
  lapply(dataset_names, function(dataset_name) {
    destination <- aerial_transitions[[dataset_name]]$transition_destination
    tibble::tibble(
      dataset                 = dataset_name,
      n_aerial_origin         = length(destination),
      n_aerial_to_aerial      = sum(destination == "aerial"),
      n_aerial_to_terrestrial = sum(destination != "aerial"))}))

# 2.3 Aerial-to-feeding and aerial-to-resting transitions, ACC datasets ----
aerial_to_acc_states_summary <- dplyr::bind_rows(
  lapply(acc_dataset_names, function(dataset_name) {
    destination <- aerial_transitions[[dataset_name]]$transition_destination
    tibble::tibble(
      dataset             = dataset_name,
      n_aerial_to_resting = sum(destination == "resting"),
      n_aerial_to_feeding = sum(destination == "feeding"))}))


cat("\n== Individuals removed (< ", minimum_aerial_transitions, " aerial transitions) ==\n", sep = "")
if (nrow(removed_individuals_summary) == 0L) {cat("None.\n")} else {print(as.data.frame(removed_individuals_summary), row.names = FALSE)}

cat("\n== Aerial-to-terrestrial transitions (ACC: resting + feeding) ==\n")
print(as.data.frame(aerial_to_terrestrial_summary), row.names = FALSE)
#dataset n_aerial_origin n_aerial_to_aerial n_aerial_to_terrestrial
#GPS_20           22456              11621                   10835
#GPS_60           10872               4244                    6628
#ACC_20           19966               9703                   10263
#ACC_60           10042               3843                    6199

cat("\n== Aerial-to-resting and aerial-to-feeding transitions (ACC) ==\n")
print(as.data.frame(aerial_to_acc_states_summary), row.names = FALSE)
# dataset n_aerial_to_resting n_aerial_to_feeding
# ACC_20                9499                 764
# ACC_60                5650                 549