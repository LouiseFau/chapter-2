#' -----------------------------------------------------------------------------
#' Title: Transition matrices preparation for two GPS and one ACC dataset ----
#' Author: Louise Faure
#' Date: 03.08.26
#'
#' Datasets:
#' (i) GPS thinned at 60 minutes;
#' (ii) GPS thinned at 20 minutes;
#' (iii) ACC aggregated at 60 minutes.
#'
#' Main steps:
#' (1) prepare behavioural states and elapsed state durations;
#' (2) construct transition matrices separately for the three datasets;
#' (3) retain transitions originating from the aerial state;
#' (4) remove individuals with fewer than 30 aerial transitions;
#' (5) standardize covariates separately for each dataset;
#' (6) assign equal total weight to each individual;
#' (7) summarize and export the three datasets.
#' -----------------------------------------------------------------------------

# Libraries ----
library(dplyr)
library(tidyr)
library(tibble)
library(sf)

# Input data ----
input_datasets <- list(
  GPS_20 = readRDS("/Users/louisefaure/Desktop/dossier sans titre/donnees filtree/GE_gps_20_covariates_hfi.rds"),
  GPS_60 = readRDS("/Users/louisefaure/Desktop/dossier sans titre/donnees filtree/GE_gps_60_min_covariates_hfi(2).rds"),
  ACC_20 = readRDS("/Users/louisefaure/Desktop/dossier sans titre/donnees filtree/GE_acc_20_covariates_hfi.rds"),
  ACC_60 = readRDS("/Users/louisefaure/Desktop/dossier sans titre/donnees filtree/GE_acc_60_min_covariates_hfi.rds"))

emig_dates_raw <- readRDS( "/Users/louisefaure/Library/CloudStorage/OneDrive-Personnel/THESE/CHAPITRE 2/git/chapter-2/DONNEES AIGLES/emigration dates/emigration_dates_20250417.rds")

# Parameters ----
state_levels <- list(GPS = c("aerial","terrestrial"), ACC = c("aerial","resting","feeding"))

dataset_parameters <- tibble::tribble(
  ~dataset, ~data_type, ~expected_dt_min, ~dt_tolerance_min, ~output_file,
  "GPS_20", "GPS", 20,  5, "gps_20_weighted.rds",
  "GPS_60", "GPS", 60, 10, "gps_60_weighted.rds",
  "ACC_20", "ACC", 20, 20, "acc_20_weighted.rds",
  "ACC_60", "ACC", 60, 60, "acc_60_weighted.rds")

minimum_aerial_transitions <- 30L
minimum_age_days <- 0
maximum_age_days <- 105

output_directory <- "/Users/louisefaure/Desktop/dossier sans titre/donnees filtree"
dir.create(output_directory,recursive = TRUE,showWarnings = FALSE)

required_numeric_variables <- c(
  "age_since_emig_days","age_since_emig_weeks","aerial_duration_min",
  "cos_diel","sin_time","ruggedness_100m","slope_100m",
  "distance_to_ridgeline_100m","elevation_100m","prop_forest_5cells",
  "prop_low_vegetation_5cells","prop_rocky_terrain_5cells",
  "prop_other_5cells","settlement_density","population_density")

environmental_hfi_variables <- setdiff(
  required_numeric_variables,
  c("age_since_emig_days","age_since_emig_weeks","aerial_duration_min","cos_diel","sin_time"))

gps_dataset_names <- dataset_parameters$dataset[dataset_parameters$data_type == "GPS"]
acc_dataset_names <- dataset_parameters$dataset[dataset_parameters$data_type == "ACC"]








#------------------------------------------------------------------------------- STEP 1: prepare datasets ----
# 1.1 Prepare emigration dates ----
emig_dates <- emig_dates_raw %>%
  dplyr::transmute(
    individual.local.identifier = trimws(as.character(individual.local.identifier)),
    dispersal_date = as.POSIXct(as.character(dispersal_date),tz = "UTC")) %>%
  dplyr::filter(!is.na(individual.local.identifier),!is.na(dispersal_date)) %>%
  dplyr::distinct(individual.local.identifier,.keep_all = TRUE)

# 1.2 Calculate age for ACC data ----
calculate_acc_age <- function(data) {
  data %>%
    dplyr::select(-dplyr::any_of(c(
      "dispersal_date","age_since_emig_days","age_since_emig_weeks"))) %>%
    dplyr::mutate(
      individual.local.identifier = trimws(as.character(individual.local.identifier)),
      timestamp = as.POSIXct(timestamp,tz = "UTC")) %>%
    dplyr::left_join(emig_dates,by = "individual.local.identifier") %>%
    dplyr::mutate(
      age_since_emig_days = as.numeric(difftime(timestamp,dispersal_date,units = "days")),
      age_since_emig_weeks = age_since_emig_days / 7) %>%
    dplyr::filter(
      age_since_emig_days >= minimum_age_days,
      age_since_emig_days <= maximum_age_days)}

datasets_with_age <- input_datasets
datasets_with_age[acc_dataset_names] <- lapply(input_datasets[acc_dataset_names],calculate_acc_age)

# 1.3 Prepare GPS behavioural states ----
prepare_locations <- function(data,data_type) {
  behavior_values <- if(data_type == "GPS") {
    data$behavior_binary %>%
      as.character() %>%
      trimws() %>%
      tolower() %>%
      dplyr::recode(
        "aerian" = "aerial",
        "flight" = "aerial",
        "flying" = "aerial",
        "ground" = "terrestrial"
      )
  } else {
    data$behavior_reclassified %>%
      as.character() %>%
      trimws() %>%
      tolower() %>%
      dplyr::recode(
        "aerian" = "aerial",
        "flight" = "aerial",
        "flying" = "aerial",
        "rest" = "resting",
        "feed" = "feeding",
        "foraging" = "feeding"
      )
  }
  
  data %>%
    dplyr::mutate(
      timestamp = as.POSIXct(timestamp,tz = "UTC"),
      individual.local.identifier = trimws(as.character(individual.local.identifier)),
      burst_id = as.character(burst_id),
      behavior_state = factor(behavior_values,levels = state_levels[[data_type]])
    ) %>%
    dplyr::arrange(individual.local.identifier,burst_id,timestamp) %>%
    dplyr::group_by(individual.local.identifier,burst_id) %>%
    dplyr::mutate(
      state_bout_n = cumsum(
        dplyr::row_number() == 1L |
          dplyr::coalesce(behavior_state != dplyr::lag(behavior_state),TRUE)
      )
    ) %>%
    dplyr::group_by(individual.local.identifier,burst_id,state_bout_n) %>%
    dplyr::mutate(
      state_duration_min = as.numeric(
        difftime(timestamp,dplyr::first(timestamp),units = "mins")
      )
    ) %>%
    dplyr::ungroup()
}

locations <- vector("list",length(input_datasets))
names(locations) <- names(input_datasets)

locations[gps_dataset_names] <- lapply(
  gps_dataset_names,
  function(dataset_name) prepare_locations(datasets_with_age[[dataset_name]],"GPS"))

# 1.4 Prepare ACC behavioural states ----
locations[acc_dataset_names] <- lapply(acc_dataset_names,  function(dataset_name) prepare_locations(datasets_with_age[[dataset_name]],"ACC"))

# 1.5 Standardization function ----
standardize_aerial_dataset <- function(data,dataset_name) {
  variables_to_standardize <- c("age_since_emig_weeks","aerial_duration_min",environmental_hfi_variables)
  
  parameters <- tibble::tibble(
    dataset = dataset_name,
    variable = variables_to_standardize,
    transformation = "center and standardize",
    center = vapply(
      variables_to_standardize,
      function(variable) mean(data[[variable]]),
      numeric(1)),
    scale = vapply(
      variables_to_standardize,
      function(variable) stats::sd(data[[variable]]),
      numeric(1)))
  
  if(any(!is.finite(parameters$scale) | parameters$scale <= 0)) {
    stop("At least one variable has a non-finite or zero SD in ",dataset_name,".")}
  
  for(i in seq_len(nrow(parameters))) {
    variable <- parameters$variable[[i]]
    data[[paste0(variable,"_z")]] <-
      (data[[variable]] - parameters$center[[i]]) / parameters$scale[[i]]
  }
  
  cos_center <- mean(data$cos_diel)
  sin_center <- mean(data$sin_time)
  
  data <- data %>%
    dplyr::mutate(
      age_z = age_since_emig_weeks_z,
      duration_z = aerial_duration_min_z,
      age_z2 = age_z^2,
      duration_z2 = duration_z^2,
      cos_diel_c = cos_diel - cos_center,
      sin_diel_c = sin_time - sin_center)
  
  diel_parameters <- tibble::tibble(
    dataset = dataset_name,
    variable = c("cos_diel","sin_time"),
    transformation = "center only",
    center = c(cos_center,sin_center),
    scale = 1)
  
  list(data = data,parameters = dplyr::bind_rows(parameters,diel_parameters))}




#------------------------------------------------------------------------------- STEP 2: construct transition matrices ----
# 2.1 Construct transitions ----
construct_transitions <- function(data,expected_dt_min,dt_tolerance_min) {
  data %>%
    dplyr::arrange(individual.local.identifier,burst_id,timestamp) %>%
    dplyr::group_by(individual.local.identifier,burst_id) %>%
    dplyr::mutate(
      behavior_from = behavior_state,
      behavior_to = dplyr::lead(behavior_state),
      timestamp_next = dplyr::lead(timestamp),
      dt_min = as.numeric(difftime(timestamp_next,timestamp,units = "mins"))
    ) %>%
    dplyr::ungroup() %>%
    dplyr::filter(
      !is.na(behavior_from),
      !is.na(behavior_to),
      is.finite(dt_min),
      dt_min > 0,
      abs(dt_min - expected_dt_min) <= dt_tolerance_min
    )
}

transitions <- setNames(
  lapply(dataset_parameters$dataset,function(dataset_name) {
    parameters_i <- dataset_parameters %>%
      dplyr::filter(dataset == dataset_name)
    
    construct_transitions(
      data = locations[[dataset_name]],
      expected_dt_min = parameters_i$expected_dt_min,
      dt_tolerance_min = parameters_i$dt_tolerance_min
    )
  }),
  dataset_parameters$dataset)

# 2.2 Calculate empirical transition matrices ----
calculate_transition_matrices <- function(data,dataset_name) {
  count_matrix <- xtabs(~ behavior_from + behavior_to,data = data)
  probability_matrix <- prop.table(count_matrix,margin = 1)
  
  count_table <- as.data.frame(count_matrix,responseName = "n_transitions")
  probability_table <- as.data.frame(
    probability_matrix,
    responseName = "transition_probability")
  
  summary <- dplyr::left_join(
    count_table,
    probability_table,
    by = c("behavior_from","behavior_to")
  ) %>%
    dplyr::mutate(dataset = dataset_name,.before = 1)
  
  list(
    count_matrix = count_matrix,
    probability_matrix = probability_matrix,
    summary = summary
  )
}

transition_matrices <- Map(calculate_transition_matrices,transitions,names(transitions))
transition_count_matrices <- lapply(transition_matrices,function(result) result$count_matrix)
transition_probability_matrices <- lapply(transition_matrices,function(result) result$probability_matrix)

transition_matrix_summary <- dplyr::bind_rows(
  lapply(transition_matrices,function(result) result$summary)) %>%
  dplyr::arrange(dataset,behavior_from,behavior_to)

# 2.3 Retain aerial-origin transitions ----
retain_aerial_transitions <- function(data,state_levels_i) {
  data %>%
    dplyr::filter(behavior_from == "aerial") %>%
    dplyr::mutate(
      transition_destination = stats::relevel(
        factor(behavior_to,levels = state_levels_i),
        ref = "aerial"
      ),
      remain_aerial = as.integer(transition_destination == "aerial"),
      aerial_duration_min = state_duration_min
    )}

aerial_transitions_raw <- setNames(
  lapply(dataset_parameters$dataset,function(dataset_name) {
    data_type_i <- dataset_parameters$data_type[
      dataset_parameters$dataset == dataset_name]
    
    retain_aerial_transitions(
      transitions[[dataset_name]],
      state_levels[[data_type_i]])}),
  dataset_parameters$dataset)

# 2.4 Count aerial transitions per individual ----
summarise_individual_transitions <- function(data,dataset_name) {
  data %>%
    dplyr::group_by(individual.local.identifier) %>%
    dplyr::summarise(
      dataset = dataset_name,
      n_aerial_transitions = dplyr::n(),
      n_bursts = dplyr::n_distinct(burst_id),
      .groups = "drop"
    )}

transitions_by_individual_list <- Map(summarise_individual_transitions,
  aerial_transitions_raw,names(aerial_transitions_raw))

transitions_by_individual <- dplyr::bind_rows(
  transitions_by_individual_list) %>%
  dplyr::arrange(dataset,n_aerial_transitions)

# 2.5 Remove low-support individuals and incomplete observations ----
prepare_final_aerial_dataset <- function(data,transition_summary,state_levels_i,dataset_name) {
  missing_variables <- setdiff(required_numeric_variables,names(data))
  
  if(length(missing_variables) > 0L) {
    stop(
      "Missing numeric variables in ",dataset_name,": ",
      paste(missing_variables,collapse = ", ")
    )
  }
  
  retained_individuals <- transition_summary %>%
    dplyr::filter(n_aerial_transitions >= minimum_aerial_transitions) %>%
    dplyr::select(individual.local.identifier)
  
  data_after_support <- data %>%
    dplyr::semi_join(retained_individuals,by = "individual.local.identifier")
  
  final_data <- data_after_support %>%
    dplyr::mutate(
      individual_id = factor(individual.local.identifier),
      transition_destination = stats::relevel(
        factor(transition_destination,levels = state_levels_i),
        ref = "aerial"
      )
    ) %>%
    tidyr::drop_na(
      individual_id,
      burst_id,
      transition_destination,
      remain_aerial,
      dplyr::all_of(required_numeric_variables)
    ) %>%
    dplyr::filter(
      dplyr::if_all(dplyr::all_of(required_numeric_variables),is.finite)
    ) %>%
    dplyr::mutate(
      individual_id = droplevels(individual_id),
      transition_destination = droplevels(transition_destination)
    )
  
  list(
    data = final_data,
    n_removed_low_support = nrow(data) - nrow(data_after_support),
    n_removed_incomplete = nrow(data_after_support) - nrow(final_data)
  )
}

prepared_datasets <- setNames(
  lapply(dataset_parameters$dataset,function(dataset_name) {
    data_type_i <- dataset_parameters$data_type[
      dataset_parameters$dataset == dataset_name
    ]
    
    prepare_final_aerial_dataset(
      data = aerial_transitions_raw[[dataset_name]],
      transition_summary = transitions_by_individual_list[[dataset_name]],
      state_levels_i = state_levels[[data_type_i]],
      dataset_name = dataset_name
    )
  }),
  dataset_parameters$dataset
)

# 2.6 Standardize each dataset separately ----
standardized_datasets <- setNames(
  lapply(dataset_parameters$dataset,function(dataset_name) {
    standardize_aerial_dataset(
      prepared_datasets[[dataset_name]]$data,
      dataset_name
    )
  }),dataset_parameters$dataset)

aerial_transitions <- lapply(
  standardized_datasets,
  function(result) result$data
)

standardization_parameters <- dplyr::bind_rows(
  lapply(standardized_datasets,function(result) result$parameters)
)

# 2.7 Give each individual the same total weight ----
add_equal_individual_weights <- function(data) {
  data %>%
    dplyr::select(-dplyr::any_of(c(
      "n_transitions_individual","individual_weight_raw","individual_weight"
    ))) %>%
    dplyr::group_by(individual_id) %>%
    dplyr::mutate(
      n_transitions_individual = dplyr::n(),
      individual_weight_raw = 1 / n_transitions_individual
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      individual_weight = individual_weight_raw / mean(individual_weight_raw)
    )
}

weighted_datasets <- lapply(aerial_transitions,add_equal_individual_weights)

#-------------------------------------------------------------------------------
# STEP 3: summaries and export ----

# 3.1 Summarise aerial-origin destinations ----
summarise_aerial_destinations <- function(data,dataset_name,state_levels_i) {
  data %>%
    dplyr::count(transition_destination,name = "n_transitions") %>%
    tidyr::complete(
      transition_destination = factor(state_levels_i,levels = state_levels_i),
      fill = list(n_transitions = 0L)
    ) %>%
    dplyr::mutate(
      dataset = dataset_name,
      transition_destination = as.character(transition_destination),
      transition_probability = n_transitions / sum(n_transitions)
    ) %>%
    dplyr::select(
      dataset,
      transition_destination,
      n_transitions,
      transition_probability
    )
}

aerial_transition_summary <- dplyr::bind_rows(
  lapply(dataset_parameters$dataset,function(dataset_name) {
    data_type_i <- dataset_parameters$data_type[
      dataset_parameters$dataset == dataset_name
    ]
    
    summarise_aerial_destinations(
      weighted_datasets[[dataset_name]],
      dataset_name,
      state_levels[[data_type_i]]
    )
  })
)

# 3.2 Summarise final datasets ----
summarise_final_dataset <- function(data,dataset_name,prepared_data) {
  data %>%
    dplyr::summarise(
      dataset = dataset_name,
      n_observations = dplyr::n(),
      n_individuals = dplyr::n_distinct(individual_id),
      n_bursts = dplyr::n_distinct(
        interaction(individual_id,burst_id,drop = TRUE)
      ),
      n_aerial_to_aerial = sum(transition_destination == "aerial"),
      n_aerial_to_terrestrial = sum(transition_destination == "terrestrial"),
      n_aerial_to_resting = sum(transition_destination == "resting"),
      n_aerial_to_feeding = sum(transition_destination == "feeding"),
      proportion_remain_aerial = mean(remain_aerial),
      n_removed_low_support = prepared_data$n_removed_low_support,
      n_removed_incomplete = prepared_data$n_removed_incomplete
    )
}

final_dataset_summary <- dplyr::bind_rows(
  lapply(dataset_parameters$dataset,function(dataset_name) {
    summarise_final_dataset(
      weighted_datasets[[dataset_name]],
      dataset_name,
      prepared_datasets[[dataset_name]]
    )
  })
)

# summaries
base::print(as.data.frame(aerial_transition_summary),row.names = FALSE)
base::print(as.data.frame(final_dataset_summary),row.names = FALSE)
base::print(as.data.frame(transitions_by_individual),row.names = FALSE)

# 3.4 Export the four weighted datasets ----
saveRDS(weighted_datasets$GPS_60,file = file.path(output_directory,"gps_60_weighted.rds"))
saveRDS(weighted_datasets$GPS_20,file = file.path(output_directory,"gps_20_weighted.rds"))
saveRDS(weighted_datasets$ACC_60,file = file.path(output_directory,"acc_60_weighted.rds"))
saveRDS(weighted_datasets$ACC_20,file = file.path(output_directory,"acc_20_weighted.rds"))


#-------------------------------------------------------------------------------
# STEP 4: test model and residual spatio-temporal autocorrelation ----
library(lme4)
library(DHARMa)
library(spdep)

# 4.1 Prepare model data ----
weight_variable_20 <- intersect(c("weight","individual_weight"),names(gps_20_weighted))[1]
if(is.na(weight_variable_20)) stop("No weight or individual_weight column found.")

gps_20_weighted <- gps_20_weighted %>%
  dplyr::mutate(
    cos_time_c = dplyr::coalesce(
      if("cos_time_c" %in% names(.)) cos_time_c else NA_real_,
      cos_diel_c
    ),
    sin_time_c = dplyr::coalesce(
      if("sin_time_c" %in% names(.)) sin_time_c else NA_real_,
      sin_diel_c
    ),
    model_weight = .data[[weight_variable_20]],
    timestamp = as.POSIXct(timestamp,tz = "UTC")
  )

required_model_variables_20 <- c(
  "remain_aerial","cos_time_c","sin_time_c","duration_z",
  "settlement_density_z","elevation_100m_z","ruggedness_100m_z",
  "individual_id","burst_id","timestamp","lon","lat","model_weight"
)

missing_model_variables_20 <- setdiff(required_model_variables_20,names(gps_20_weighted))
if(length(missing_model_variables_20) > 0L) {
  stop("Missing variables: ",paste(missing_model_variables_20,collapse = ", "))
}

model_data_20 <- gps_20_weighted %>%
  dplyr::select(dplyr::all_of(required_model_variables_20)) %>%
  tidyr::drop_na() %>%
  dplyr::filter(
    model_weight > 0,
    dplyr::if_all(
      c(
        remain_aerial,cos_time_c,sin_time_c,duration_z,
        settlement_density_z,elevation_100m_z,ruggedness_100m_z,
       lon,lat,model_weight
      ),
      is.finite
    )
  ) %>%
  dplyr::arrange(individual_id,burst_id,timestamp) %>%
  dplyr::mutate(individual_id = droplevels(factor(individual_id)))

# 4.2 Fit binomial GLMM ----
model_formula_20 <- remain_aerial ~ cos_time_c + sin_time_c + duration_z +
  settlement_density_z + elevation_100m_z + ruggedness_100m_z +
  (1 | individual_id)

model_20 <- suppressWarnings(
  lme4::glmer(
    formula = model_formula_20,
    data = model_data_20,
    weights = model_weight,
    family = stats::binomial(link = "logit"),
    control = lme4::glmerControl(
      optimizer = "bobyqa",
      optCtrl = list(maxfun = 2e5)
    )
  )
)

print(summary(model_20))

# 4.3 Generate DHARMa pseudo-residuals ----
# 4.3 Generate conditional DHARMa pseudo-residuals ----
set.seed(123)

fitted_probability_20 <- stats::predict(
  model_20,
  type = "response",
  re.form = NULL
)

simulated_response_20 <- replicate(
  500,
  stats::rbinom(
    n = nrow(model_data_20),
    size = 1,
    prob = fitted_probability_20
  )
)

pseudo_residuals_20 <- DHARMa::createDHARMa(
  simulatedResponse = simulated_response_20,
  observedResponse = model_data_20$remain_aerial,
  fittedPredictedResponse = fitted_probability_20,
  integerResponse = TRUE
)

plot(pseudo_residuals_20)

residual_data_20 <- model_data_20 %>%
  dplyr::mutate(
    pseudo_residual = pseudo_residuals_20$scaledResiduals
  )

residual_points_20 <- residual_data_20 %>%
  sf::st_as_sf(
    coords = c("lon","lat"),
    crs = 4326,
    remove = FALSE
  ) %>%
  sf::st_transform(3035)

coordinates_3035_20 <- sf::st_coordinates(residual_points_20)

residual_data_20 <- residual_data_20 %>%
  dplyr::mutate(
    x_3035 = coordinates_3035_20[,1],
    y_3035 = coordinates_3035_20[,2]
  )

# 4.4 Temporal autocorrelation within bursts ----
temporal_autocorrelation_by_burst_20 <- residual_data_20 %>%
  dplyr::group_by(individual_id,burst_id) %>%
  dplyr::filter(dplyr::n() >= 5L) %>%
  dplyr::arrange(timestamp,.by_group = TRUE) %>%
  dplyr::group_modify(~{
    residual_vector <- .x$pseudo_residual
    
    correlation_test <- tryCatch(
      stats::cor.test(
        residual_vector[-length(residual_vector)],
        residual_vector[-1],
        method = "pearson"
      ),
      error = function(e) NULL
    )
    
    if(is.null(correlation_test)) {
      return(tibble::tibble(
        n = length(residual_vector),
        lag1_correlation = NA_real_,
        p_value = NA_real_
      ))
    }
    
    tibble::tibble(
      n = length(residual_vector),
      lag1_correlation = unname(correlation_test$estimate),
      p_value = correlation_test$p.value
    )
  }) %>%
  dplyr::ungroup() %>%
  dplyr::filter(is.finite(lag1_correlation),is.finite(p_value)) %>%
  dplyr::mutate(p_adjusted = stats::p.adjust(p_value,method = "BH"))

temporal_summary_20 <- temporal_autocorrelation_by_burst_20 %>%
  dplyr::summarise(
    statistic = stats::weighted.mean(
      lag1_correlation,
      w = pmax(n - 1,1)
    ),
    fisher_statistic = -2 * sum(
      log(pmax(p_value,.Machine$double.xmin))
    ),
    p_value = stats::pchisq(
      fisher_statistic,
      df = 2 * dplyr::n(),
      lower.tail = FALSE
    ),
    n_series = dplyr::n(),
    n_significant = sum(p_adjusted < 0.05)
  )

# 4.5 Spatial autocorrelation of pseudo-residuals ----
# 4.5 Spatial autocorrelation of pseudo-residuals ----
spatial_grid_m_20 <- 1000

spatial_residuals_20 <- residual_data_20 %>%
  dplyr::mutate(
    grid_x = floor(x_3035 / spatial_grid_m_20),
    grid_y = floor(y_3035 / spatial_grid_m_20)
  ) %>%
  dplyr::group_by(grid_x,grid_y) %>%
  dplyr::summarise(
    x_3035 = mean(x_3035),
    y_3035 = mean(y_3035),
    pseudo_residual = mean(pseudo_residual),
    .groups = "drop"
  )

spatial_k_20 <- min(8L,nrow(spatial_residuals_20) - 1L)

spatial_neighbours_20 <- spdep::knn2nb(
  spdep::knearneigh(
    cbind(
      spatial_residuals_20$x_3035,
      spatial_residuals_20$y_3035
    ),
    k = spatial_k_20
  )
)

spatial_weights_20 <- spdep::nb2listw(
  spatial_neighbours_20,
  style = "W",
  zero.policy = TRUE
)

spatial_test_20 <- spdep::moran.test(
  spatial_residuals_20$pseudo_residual,
  spatial_weights_20,
  alternative = "two.sided",
  zero.policy = TRUE
)

print(spatial_test_20)
# 4.6 Compile results ----
autocorrelation_summary_20 <- dplyr::bind_rows(
  temporal_summary_20 %>%
    dplyr::transmute(
      test = "Temporal lag-1 within bursts",
      statistic,
      p_value,
      n_units = n_series
    ),
  tibble::tibble(
    test = paste0("Spatial Moran's I, ",spatial_grid_m_20," m cells"),
    statistic = unname(spatial_test_20$estimate[1]),
    p_value = spatial_test_20$p.value,
    n_units = nrow(spatial_residuals_20)
  )
)

print(autocorrelation_summary_20)
print(temporal_autocorrelation_by_burst_20,n = Inf)
