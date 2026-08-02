#' -----------------------------------------------------------------------------
# Title: Matrices of transition preparation for ACC and GPS behaviorally classified datasets ----
#' Authors : Louise Faure
#' Date : 30.07.26
#' 
#' Info : this script follow the Extract_covariates.R script where covariates are
#' extracted below each location and within two buffers. This script can be applied
#' to both GPS and ACC data. 
#' 
# Main steps:
#' (1) prepare the dataset:
#'     (i) calculate elapsed duration in the current behavioural state for both
#'     acc and gps data;
#'     (ii) for acc data, calculate age since emigration
#'     (iii) centre and standardize age and elapsed aerial duration
#'     (iv) standardize environmental covariates and HFI metrics separately for 
#'     GPS and ACC data;
#'     
#' (2) construct two transition matrix transition matrix and retain 
#'     aerial-origin transitions;
#'     (i) construct transition matrix
#'     (ii) give a weight to each individuals
#'     (iii) remove individuals with less than 30 transition from aerial
#'     (iv) summarize the dataset information 
#' 
#' (3) summary statistics and export
#'    (i) print the number of transition from aerial to feeding, resting and 
#'    terrestrial as well as the probabilities 
#'    (ii) export the acc and gps datasets
#'------------------------------------------------------------------------------


# library
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)


# acc and gps data, emigration dates ----
GE_60_gps_covariates_hfi <- readRDS("/Users/louisefaure/Desktop/dossier sans titre/donnees filtree/GE_gps_60_min_covariates_hfi(2).rds")
GE_60_acc_covariates_hfi <- readRDS("/Users/louisefaure/Desktop/dossier sans titre/donnees filtree/GE_acc_60_min_covariates_hfi(2).rds")
emig_dates_60_raw <- readRDS("/Users/louisefaure/Library/CloudStorage/OneDrive-Personnel/THESE/CHAPITRE 2/git/chapter-2/DONNEES AIGLES/emigration dates/emigration_dates_20250417.rds")

# Parameters ----
state_levels_gps_60 <- c("aerial","terrestrial")
state_levels_acc_60 <- c("aerial","resting","feeding")

required_numeric_variables_60 <- c(
  "age_since_emig_days","age_since_emig_weeks","aerial_duration_min","cos_diel","sin_time",
  "ruggedness_100m","slope_100m","distance_to_ridgeline_100m","elevation_100m",
  "prop_forest_5cells","prop_low_vegetation_5cells","prop_rocky_terrain_5cells",
  "prop_other_5cells", "settlement_density", "population_density"
)

# Additional parameters ----
expected_dt_min_60 <- 60
dt_tolerance_min_60 <- 60
minimum_aerial_transitions_60 <- 30L
minimum_age_days_60 <- 0
maximum_age_days_60 <- 105

output_directory_60 <- "/Users/louisefaure/Desktop/dossier sans titre/donnees filtree"
dir.create(output_directory_60,recursive = TRUE,showWarnings = FALSE)

environmental_hfi_variables_60 <- setdiff(required_numeric_variables_60,c("age_since_emig_days","age_since_emig_weeks","aerial_duration_min","cos_diel","sin_time"))

#------------------------------------------------------------------------------ STEP 1: prepare GPS and ACC datasets ----
#' **Steps:**
#' (i) calculate elapsed duration in the current behavioural state;
#' (ii) calculate age since emigration for ACC data;
#' (iii) prepare age, duration and diel variables;
#' (iv) define dataset-specific standardization of environmental covariates
#'      and HFI metrics.

# 1.1 Prepare emigration dates and calculate ACC age ----
emig_dates_60 <- emig_dates_60_raw %>%
  dplyr::transmute(
    individual.local.identifier = trimws(as.character(individual.local.identifier)),
    dispersal_date = as.POSIXct(as.character(dispersal_date),tz = "UTC")
  ) %>%
  dplyr::filter(!is.na(individual.local.identifier),!is.na(dispersal_date)) %>%
  dplyr::distinct()

acc_data_with_age_60 <- GE_60_acc_covariates_hfi %>%
  dplyr::select(
    -dplyr::any_of(c(
      "dispersal_date","age_since_emig_days","age_since_emig_weeks"
    ))
  ) %>%
  dplyr::mutate(
    individual.local.identifier = trimws(as.character(individual.local.identifier)),
    timestamp = as.POSIXct(timestamp,tz = "UTC")
  ) %>%
  dplyr::left_join(emig_dates_60,by = "individual.local.identifier") %>%
  dplyr::mutate(
    age_since_emig_days = as.numeric(
      difftime(timestamp,dispersal_date,units = "days")
    ),
    age_since_emig_weeks = age_since_emig_days / 7
  ) %>%
  dplyr::filter(
    age_since_emig_days >= minimum_age_days_60,
    age_since_emig_days <= maximum_age_days_60
  )

# 1.2 Prepare GPS behavioural states and state durations ----
locations_gps_60 <- GE_60_gps_covariates_hfi %>%
  dplyr::mutate(
    timestamp = as.POSIXct(timestamp,tz = "UTC"),
    individual.local.identifier = trimws(as.character(individual.local.identifier)),
    behavior_state = behavior_binary %>%
      as.character() %>%
      trimws() %>%
      tolower() %>%
      dplyr::recode(
        "aerian" = "aerial",
        "flight" = "aerial",
        "flying" = "aerial",
        "ground" = "terrestrial"
      ) %>%
      factor(levels = state_levels_gps_60)
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

# 1.3 Prepare ACC behavioural states and state durations ----
locations_acc_60 <- acc_data_with_age_60 %>%
  dplyr::mutate(
    behavior_state = behavior_reclassified %>%
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
      ) %>%
      factor(levels = state_levels_acc_60)
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

# 1.4 Define dataset-specific standardization ----
standardize_aerial_dataset_60 <- function(data,dataset_name) {
  variables_to_standardize <- c(
    "age_since_emig_weeks",
    "aerial_duration_min",
    environmental_hfi_variables_60
  )
  
  parameters <- tibble::tibble(
    dataset = dataset_name,
    variable = variables_to_standardize,
    transformation = "center and standardize",
    center = vapply(
      variables_to_standardize,
      function(variable) mean(data[[variable]]),
      numeric(1)
    ),
    scale = vapply(
      variables_to_standardize,
      function(variable) stats::sd(data[[variable]]),
      numeric(1)
    )
  )
  
  if(any(!is.finite(parameters$scale) | parameters$scale <= 0)) {
    stop("At least one variable has a non-finite or zero SD in ",dataset_name,".")
  }
  
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
      sin_diel_c = sin_time - sin_center
    )
  
  diel_parameters <- tibble::tibble(
    dataset = dataset_name,
    variable = c("cos_diel","sin_time"),
    transformation = "center only",
    center = c(cos_center,sin_center),
    scale = 1
  )
  
  list(
    data = data,
    parameters = dplyr::bind_rows(parameters,diel_parameters)
  )
}

#------------------------------------------------------------------------------ STEP 2: construct transitions and retain aerial origins ----
#' **Steps:**
#' (i) construct GPS and ACC transition matrices;
#' (ii) retain aerial-origin transitions;
#' (iii) remove individuals with fewer than 30 aerial-origin transitions;
#' (iv) remove incomplete observations;
#' (v) standardize variables and give equal total weight to each individual.

# 2.1 Construct transitions between consecutive locations ----
construct_transitions_60 <- function(data) {
  data %>%
    dplyr::group_by(individual.local.identifier,burst_id) %>%
    dplyr::mutate(
      behavior_from = behavior_state,
      behavior_to = dplyr::lead(behavior_state),
      timestamp_next = dplyr::lead(timestamp),
      dt_min = as.numeric(
        difftime(timestamp_next,timestamp,units = "mins")
      )
    ) %>%
    dplyr::ungroup() %>%
    dplyr::filter(
      !is.na(behavior_from),
      !is.na(behavior_to),
      is.finite(dt_min),
      dt_min > 0,
      abs(dt_min - expected_dt_min_60) <= dt_tolerance_min_60
    )
}

transitions_gps_60 <- construct_transitions_60(locations_gps_60)
transitions_acc_60 <- construct_transitions_60(locations_acc_60)

# 2.2 Calculate empirical transition matrices ----
transition_count_matrix_gps_60 <- xtabs(
  ~ behavior_from + behavior_to,
  data = transitions_gps_60)

transition_count_matrix_acc_60 <- xtabs(
  ~ behavior_from + behavior_to,
  data = transitions_acc_60)

transition_probability_matrix_gps_60 <- prop.table(
  transition_count_matrix_gps_60,
  margin = 1)

transition_probability_matrix_acc_60 <- prop.table(
  transition_count_matrix_acc_60,
  margin = 1)

compile_transition_matrix_60 <- function(count_matrix,probability_matrix,dataset_name) {
  count_table <- as.data.frame(
    count_matrix,
    responseName = "n_transitions")
  
  probability_table <- as.data.frame(
    probability_matrix,
    responseName = "transition_probability")
  
  dplyr::left_join(
    count_table,
    probability_table,
    by = c("behavior_from","behavior_to")
  ) %>%
    dplyr::mutate(dataset = dataset_name,.before = 1)
}

transition_matrix_summary_60 <- dplyr::bind_rows(
  compile_transition_matrix_60(
    transition_count_matrix_gps_60,
    transition_probability_matrix_gps_60,
    "GPS"
  ),
  compile_transition_matrix_60(
    transition_count_matrix_acc_60,
    transition_probability_matrix_acc_60,
    "ACC"
  )
) %>%
  dplyr::arrange(dataset,behavior_from,behavior_to)

# 2.3 Retain transitions originating from the aerial state ----
aerial_transitions_raw_gps_60 <- transitions_gps_60 %>%
  dplyr::filter(behavior_from == "aerial") %>%
  dplyr::mutate(
    transition_destination = stats::relevel(
      factor(behavior_to,levels = state_levels_gps_60),
      ref = "aerial"
    ),
    remain_aerial = as.integer(transition_destination == "aerial"),
    aerial_duration_min = state_duration_min
  )

aerial_transitions_raw_acc_60 <- transitions_acc_60 %>%
  dplyr::filter(behavior_from == "aerial") %>%
  dplyr::mutate(
    transition_destination = stats::relevel(
      factor(behavior_to,levels = state_levels_acc_60),
      ref = "aerial"
    ),
    remain_aerial = as.integer(transition_destination == "aerial"),
    aerial_duration_min = state_duration_min)

# 2.4 Count aerial-origin transitions per individual ----
summarise_individual_transitions_60 <- function(data,dataset_name) {
  data %>%
    dplyr::group_by(individual.local.identifier) %>%
    dplyr::summarise(
      dataset = dataset_name,
      n_aerial_transitions = dplyr::n(),
      n_bursts = dplyr::n_distinct(burst_id),
      .groups = "drop"
    )}

transitions_by_individual_gps_60 <- summarise_individual_transitions_60(aerial_transitions_raw_gps_60,"GPS")
transitions_by_individual_acc_60 <- summarise_individual_transitions_60(aerial_transitions_raw_acc_60,"ACC")

transitions_by_individual_60 <- dplyr::bind_rows(
  transitions_by_individual_gps_60,
  transitions_by_individual_acc_60
) %>%
  dplyr::arrange(dataset,n_aerial_transitions)

# 2.5 Remove low-support individuals and incomplete observations ----
prepare_final_aerial_dataset_60 <- function(data,transition_summary,state_levels) {
  retained_individuals <- transition_summary %>%
    dplyr::filter(n_aerial_transitions >= minimum_aerial_transitions_60) %>%
    dplyr::select(individual.local.identifier)
  
  data_after_support <- data %>%
    dplyr::semi_join(
      retained_individuals,
      by = "individual.local.identifier"
    )
  
  final_data <- data_after_support %>%
    dplyr::mutate(
      individual_id = factor(individual.local.identifier),
      transition_destination = stats::relevel(
        factor(transition_destination,levels = state_levels),
        ref = "aerial"
      )
    ) %>%
    tidyr::drop_na(
      individual_id,
      burst_id,
      transition_destination,
      remain_aerial,
      dplyr::all_of(required_numeric_variables_60)
    ) %>%
    dplyr::filter(
      dplyr::if_all(
        dplyr::all_of(required_numeric_variables_60),
        is.finite
      )
    ) %>%
    dplyr::mutate(
      individual_id = droplevels(individual_id),
      transition_destination = droplevels(transition_destination)
    )
  
  list(
    data = final_data,
    n_removed_low_support = nrow(data) - nrow(data_after_support),
    n_removed_incomplete = nrow(data_after_support) - nrow(final_data)
  )}

prepared_gps_60 <- prepare_final_aerial_dataset_60(
  data = aerial_transitions_raw_gps_60,
  transition_summary = transitions_by_individual_gps_60,
  state_levels = state_levels_gps_60)

prepared_acc_60 <- prepare_final_aerial_dataset_60(
  data = aerial_transitions_raw_acc_60,
  transition_summary = transitions_by_individual_acc_60,
  state_levels = state_levels_acc_60)

# 2.6 Standardize age, duration, environmental covariates and HFI ----
standardized_gps_60 <- standardize_aerial_dataset_60(prepared_gps_60$data,"GPS")
standardized_acc_60 <- standardize_aerial_dataset_60(prepared_acc_60$data,"ACC")

aerial_transitions_gps_60 <- standardized_gps_60$data
aerial_transitions_acc_60 <- standardized_acc_60$data

standardization_parameters_60 <- dplyr::bind_rows(
  standardized_gps_60$parameters,
  standardized_acc_60$parameters)

# 2.7 Give each retained individual the same total weight ----
add_equal_individual_weights_60 <- function(data) {
  data %>%
    dplyr::select(
      -dplyr::any_of(
        c(
          "n_transitions_individual",
          "individual_weight_raw",
          "individual_weight"
        )
      )
    ) %>%
    dplyr::group_by(individual_id) %>%
    dplyr::mutate(
      n_transitions_individual = dplyr::n(),
      individual_weight_raw = 1 / n_transitions_individual
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      individual_weight =
        individual_weight_raw / mean(individual_weight_raw))}

gps_weighted <- add_equal_individual_weights_60(aerial_transitions_gps_60)

acc_weighted <- add_equal_individual_weights_60(aerial_transitions_acc_60)


#------------------------------------------------------------------------------ STEP 3: summary statistics and export ----
#' **Steps:**
#' (i) summarise transition counts and probabilities;
#' (ii) summarise retained observations, individuals and bursts;
#' (iii) export weighted ACC and GPS datasets.

# 3.1 Summarise aerial-origin destinations ----
summarise_aerial_destinations_60 <- function(data,dataset_name,state_levels) {
  data %>%
    dplyr::count(transition_destination,name = "n_transitions") %>%
    tidyr::complete(
      transition_destination = factor(
        state_levels,
        levels = state_levels
      ),
      fill = list(n_transitions = 0L)
    ) %>%
    dplyr::mutate(
      dataset = dataset_name,
      transition_destination = as.character(transition_destination),
      transition_probability =
        n_transitions / sum(n_transitions)
    ) %>%
    dplyr::select(
      dataset,
      transition_destination,
      n_transitions,
      transition_probability
    )}

aerial_transition_summary_60 <- dplyr::bind_rows(
  summarise_aerial_destinations_60(
    gps_weighted,
    "GPS",
    state_levels_gps_60
  ),
  summarise_aerial_destinations_60(
    acc_weighted,
    "ACC",
    state_levels_acc_60))

# 3.2 Summarise final weighted datasets ----
final_dataset_summary_60 <- dplyr::bind_rows(
  gps_weighted %>%
    dplyr::summarise(
      dataset = "GPS",
      n_observations = dplyr::n(),
      n_individuals = dplyr::n_distinct(individual_id),
      n_bursts = dplyr::n_distinct(
        interaction(individual_id,burst_id,drop = TRUE)
      ),
      n_aerial_to_aerial = sum(transition_destination == "aerial"),
      n_aerial_to_terrestrial = sum(transition_destination == "terrestrial"),
      n_aerial_to_resting = 0L,
      n_aerial_to_feeding = 0L,
      proportion_remain_aerial = mean(remain_aerial),
      n_removed_low_support = prepared_gps_60$n_removed_low_support,
      n_removed_incomplete = prepared_gps_60$n_removed_incomplete
    ),
  acc_weighted %>%
    dplyr::summarise(
      dataset = "ACC",
      n_observations = dplyr::n(),
      n_individuals = dplyr::n_distinct(individual_id),
      n_bursts = dplyr::n_distinct(
        interaction(individual_id,burst_id,drop = TRUE)
      ),
      n_aerial_to_aerial = sum(transition_destination == "aerial"),
      n_aerial_to_terrestrial = 0L,
      n_aerial_to_resting = sum(transition_destination == "resting"),
      n_aerial_to_feeding = sum(transition_destination == "feeding"),
      proportion_remain_aerial = mean(remain_aerial),
      n_removed_low_support = prepared_acc_60$n_removed_low_support,
      n_removed_incomplete = prepared_acc_60$n_removed_incomplete))

# 3.3 Display summaries ----
base::print(as.data.frame(transition_matrix_summary_60),row.names = FALSE)
base::print( as.data.frame(aerial_transition_summary_60),row.names = FALSE)
base::print(as.data.frame(final_dataset_summary_60),row.names = FALSE)
# dataset transition_destination n_transitions transition_probability
# GPS                 aerial          4244             0.39036056
# GPS            terrestrial          6628             0.60963944
# ACC                 aerial          3822             0.38193265
# ACC                resting          5637             0.56330569
# ACC                feeding           548             0.05476167

# 3.4 Export weighted GPS and ACC datasets ----
saveRDS(gps_weighted,file = file.path(output_directory_60,"gps_weighted(2).rds"))
saveRDS(acc_weighted,file = file.path(output_directory_60,"acc_weighted(2).rds"))
