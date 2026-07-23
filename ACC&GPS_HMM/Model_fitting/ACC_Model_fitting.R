#' ----------------------------------------------------------------------------- 
# Title: Covariate and HFI variable selection for ACC-classified behaviours ----
#' Author: Louise Faure
#' Date: 22.07.26
#'
#' Info:
#' This script follows "prepare gps-acc data.R", in which environmental
#' covariates and HFI variables were extracted at each location and within
#' spatial buffers.
#'
#' Main steps:
#' (1) prepare the three-state transition dataset;
#' (2) fit a population-level three-state transition model;
#' (3) fit a model with individual variation in HFI responses.
#' -----------------------------------------------------------------------------


# Libraries ----
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(mgcv)
library(corrplot)

# Data ----
GE_60_min_covariates_hfi <- readRDS("/Users/louisefaure/Desktop/dossier sans titre/donnees filtree/GE_acc_class_60_min_covariates_hfi.rds")
emig_date <- readRDS("/Users/louisefaure/Library/CloudStorage/OneDrive-Personnel/THESE/CHAPITRE 2/git/chapter-2/DONNEES AIGLES/emigration dates/emigration_dates_20250417.rds")

# Parameters ----
state_levels_60 <- c("aerial","resting","feeding")
transition_levels_60 <- tidyr::expand_grid(
  behavior_from =
    state_levels_60,
  
  behavior_to =
    state_levels_60
) %>%
  dplyr::transmute(
    transition_type =
      paste(
        behavior_from,
        behavior_to,
        sep = "_to_"
      )
  ) %>%
  dplyr::pull(
    transition_type)




#------------------------------------------------------------------------------ STEP 1: prepare the three-state transition dataset ----
#' **Steps:**
#' (i) clean the ACC behavioural classification;
#' (ii) calculate elapsed time in the current behavioural state;
#' (iii) construct transitions between consecutive 60-minute locations;
#' (iv) calculate the empirical three-state transition matrix;
#' (v) retain transitions originating from the aerial state;
#' (vi) document and remove individuals with insufficient information;
#' (vii) remove incomplete observations.


# 1.1 Prepare behavioural states and identify behavioural bouts ----
locations_60 <- GE_60_min_covariates_hfi %>%
  dplyr::mutate(
    timestamp =
      as.POSIXct(
        timestamp,
        tz = "UTC"
      ),
    
    behavior_state =
      behavior_reclassified %>%
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
      factor(
        levels =
          state_levels_60
      ),
    
    individual_id =
      factor(
        individual.local.identifier
      )
  ) %>%
  dplyr::arrange(
    individual.local.identifier,
    burst_id,
    timestamp
  ) %>%
  dplyr::group_by(
    individual.local.identifier,
    burst_id
  ) %>%
  dplyr::mutate(
    # Start a new behavioural bout at the beginning of a burst
    # or whenever the behavioural state changes.
    state_bout_n =
      cumsum(
        dplyr::row_number() == 1L |
          dplyr::coalesce(
            behavior_state !=
              dplyr::lag(
                behavior_state
              ),
            TRUE
          )
      )
  ) %>%
  dplyr::group_by(
    individual.local.identifier,
    burst_id,
    state_bout_n
  ) %>%
  dplyr::mutate(
    # Elapsed time already spent in the current behavioural state.
    state_duration_min =
      as.numeric(
        difftime(
          timestamp,
          dplyr::first(
            timestamp
          ),
          units = "mins"
        )
      )
  ) %>%
  dplyr::ungroup()

# 1.2 Construct transitions between consecutive 60-minute locations ----
transitions_60 <- locations_60 %>%
  dplyr::group_by(
    individual.local.identifier,
    burst_id
  ) %>%
  dplyr::mutate(
    behavior_from =
      behavior_state,
    
    behavior_to =
      dplyr::lead(
        behavior_state
      ),
    
    timestamp_next =
      dplyr::lead(
        timestamp
      ),
    
    dt_min =
      as.numeric(
        difftime(
          timestamp_next,
          timestamp,
          units = "mins"
        )
      )
  ) %>%
  dplyr::ungroup() %>%
  dplyr::filter(
    !is.na(
      behavior_from
    ),
    
    !is.na(
      behavior_to
    ),
    
    is.finite(
      dt_min
    ),
    
    abs(
      dt_min -
        expected_dt_min_60
    ) <=
      dt_tolerance_min_60
  ) %>%
  dplyr::mutate(
    transition_type =
      factor(
        paste(
          behavior_from,
          behavior_to,
          sep = "_to_"
        ),
        levels =
          transition_levels_60
      )
  )


# 1.3 Calculate the empirical three-state transition matrix ----
transition_count_matrix_60 <- xtabs(
  ~ behavior_from + behavior_to,
  data =
    transitions_60)

transition_probability_matrix_60 <- prop.table(
  transition_count_matrix_60,
  margin = 1)

print(transition_count_matrix_60)
print(round(transition_probability_matrix_60,digits = 4))
print(rowSums(transition_probability_matrix_60))

# 1.4 Retain transitions originating from the aerial state ----
aerial_transitions_raw_60 <- transitions_60 %>%
  dplyr::filter(
    behavior_from == "aerial"
  ) %>%
  dplyr::mutate(
    # Main three-state response.
    # The first factor level, aerial, is the reference category.
    transition_destination =
      stats::relevel(
        factor(
          behavior_to,
          levels =
            state_levels_60
        ),
        ref = "aerial"
      ),
    
    # Numeric coding required by mgcv::multinom():
    # aerial = 0, resting = 1, feeding = 2.
    transition_destination_code =
      as.integer(
        transition_destination
      ) -
      1L,
    
    # Auxiliary binary response retained only for comparison with the
    # previous GPS aerial-versus-terrestrial model.
    remain_aerial =
      as.integer(
        transition_destination == "aerial"
      ),
    
    aerial_duration_min =
      state_duration_min)


# 1.5 Inspect aerial-origin transitions by individual ----
transitions_by_individual_60 <- aerial_transitions_raw_60 %>%
  dplyr::group_by(
    individual.local.identifier
  ) %>%
  dplyr::summarise(
    n_aerial_transitions =
      dplyr::n(),
    
    n_aerial_to_aerial =
      sum(
        transition_destination == "aerial"
      ),
    
    n_aerial_to_resting =
      sum(
        transition_destination == "resting"
      ),
    
    n_aerial_to_feeding =
      sum(
        transition_destination == "feeding"
      ),
    
    proportion_aerial =
      mean(
        transition_destination == "aerial"
      ),
    
    proportion_resting =
      mean(
        transition_destination == "resting"
      ),
    
    proportion_feeding =
      mean(
        transition_destination == "feeding"
      ),
    
    n_bursts =
      dplyr::n_distinct(
        burst_id
      ),
    
    .groups = "drop"
  ) %>%
  dplyr::arrange(
    n_aerial_transitions
  )

print(transitions_by_individual_60,n = Inf)

# 1.6 Individuals to remove ----
# We exclude the individuals with less that 30 transitions, at the end of this stage, over the 66 birds, 2 having no acc data for the period under study
# 6 individual having less than 30 transitions, we therefore have 58 individuals left. 
individuals_to_remove_60 <- c("Mals2_20 (eobs 7579)","Langgries21 (eobs 7586)","Almen18 (eobs 5861)","Untersberg21 (eobs 7501)","Schreital22 (eobs 10537)","Schlanders18 (eobs 6226)")
aerial_transitions_60 <- aerial_transitions_raw_60 %>%
  dplyr::filter(
    !individual.local.identifier %in%
      individuals_to_remove_60
  ) %>%
  dplyr::mutate(
    individual_id =
      droplevels(
        factor(
          individual.local.identifier
        )
      ),
    
    transition_destination =
      stats::relevel(
        droplevels(
          factor(
            transition_destination,
            levels =
              state_levels_60
          )
        ),
        ref = "aerial"
      ),
    
    transition_destination_code =
      as.integer(
        transition_destination
      ) -
      1L
  )

# 1.8 Summarise the final three-state transition dataset ----
final_dataset_summary_60 <- aerial_transitions_60 %>%
  dplyr::summarise(
    n_observations =
      dplyr::n(),
    
    n_individuals =
      dplyr::n_distinct(
        individual_id
      ),
    
    n_bursts =
      dplyr::n_distinct(
        burst_id
      ),
    
    n_aerial_to_aerial =
      sum(
        transition_destination == "aerial"
      ),
    
    n_aerial_to_resting =
      sum(
        transition_destination == "resting"
      ),
    
    n_aerial_to_feeding =
      sum(
        transition_destination == "feeding"
      ),
    
    proportion_aerial =
      mean(
        transition_destination == "aerial"
      ),
    
    proportion_resting =
      mean(
        transition_destination == "resting"
      ),
    
    proportion_feeding =
      mean(
        transition_destination == "feeding"
      )
  )

print(final_dataset_summary_60)




# ------------------------------------------------------------------------------ STEP 2 : model fitting ----
#' **Philisophy:** the reference model is : remain_aerial ~ cos_diel + sin_time + age_z + duration_z + duration_z2 + 
#' s(individual_id, bs = "re") + hfi_mean_1000m_z + dem_elevation_z + prop_low_vegetation_5cells_z
#' We fitted this model with the GPS-based classification of behaviors for two behavioral state.
#'
#' **Steps:**
#' (i) prepare diel and standardized covariates;
#' (ii) fit the three-state population-level model;
#' (iii) extract the HFI effects for resting and feeding relative to aerial;
#' (iv) calculate Q05-Q95 HFI log-odds contrasts.

# 2.0. Calculate age since emigration ----
# 2.0.1 Correct the individual Stürfis20 (eobs 7049) name ----
emig_date <- emig_date %>%
  dplyr::mutate(
    individual.local.identifier =
      trimws(
        as.character(
          individual.local.identifier)))

# 2.0.2 Prepare dispersal dates ----
emig_date_60 <- emig_date %>%
  dplyr::transmute(
    individual.local.identifier =
      as.character(
        individual.local.identifier
      ),
    
    dispersal_date =
      as.Date(
        dispersal_date
      )
  ) %>%
  dplyr::filter(
    !is.na(
      individual.local.identifier),
    !is.na(
      dispersal_date
    )
  ) %>%
  dplyr::distinct()

# 2.0.2 Join dispersal dates and calculate age (day since emigration) ----
aerial_transitions_60 <- aerial_transitions_60 %>%
  dplyr::select(
    -dplyr::any_of(
      c(
        "dispersal_date",
        "dispersal_date.x",
        "dispersal_date.y",
        "age_days",
        "age_weeks"
      )
    )
  ) %>%
  dplyr::mutate(
    individual.local.identifier =
      as.character(
        individual.local.identifier
      ),
    
    timestamp =
      as.POSIXct(
        timestamp,
        tz = "UTC"
      )
  ) %>%
  dplyr::left_join(
    emig_date_60,
    by =
      "individual.local.identifier"
  ) %>%
  dplyr::mutate(
    age_days =
      as.numeric(
        as.Date(
          timestamp,
          tz = "UTC"
        ) -
          dispersal_date))


# 2.1 Raw models parameters ----
# 2.1.0 Parameters and variable names ----
age_variable_60 <- "age_days"
hfi_variable_60 <- "hfi_mean_1000m"
elevation_variable_60 <- "elevation_100m"
low_vegetation_variable_60 <- "prop_low_vegetation_5cells"
confidence_level_60 <- 0.95
confidence_multiplier_60 <- stats::qnorm(1 - (1 - confidence_level_60) / 2)

# 2.1.1 Prepare the raw model dataset ----
three_state_model_data_raw_60 <- aerial_transitions_60 %>%
  dplyr::mutate(
    aerial_duration_min =
      state_duration_min,
    
    individual_id =
      droplevels(
        factor(
          individual.local.identifier
        )
      ),
    
    transition_destination =
      factor(
        transition_destination,
        levels = c(
          "aerial",
          "resting",
          "feeding"
        )
      ),
    
    transition_destination_code =
      as.integer(
        transition_destination
      ) -
      1L)

# 2.2 Calculate standardization parameters ----
# 2.2.1 Mean and ecart type for each variables ----
standardization_parameters_60 <- tibble::tibble(
  variable = c(
    age_variable_60,
    "aerial_duration_min",
    hfi_variable_60,
    elevation_variable_60,
    low_vegetation_variable_60
  ),
  
  mean = c(
    mean(
      three_state_model_data_raw_60[[
        age_variable_60
      ]],
      na.rm = TRUE
    ),
    
    mean(
      three_state_model_data_raw_60$
        aerial_duration_min,
      na.rm = TRUE
    ),
    
    mean(
      three_state_model_data_raw_60[[
        hfi_variable_60
      ]],
      na.rm = TRUE
    ),
    
    mean(
      three_state_model_data_raw_60[[
        elevation_variable_60
      ]],
      na.rm = TRUE
    ),
    
    mean(
      three_state_model_data_raw_60[[
        low_vegetation_variable_60
      ]],
      na.rm = TRUE
    )
  ),
  
  standard_deviation = c(
    stats::sd(
      three_state_model_data_raw_60[[
        age_variable_60
      ]],
      na.rm = TRUE
    ),
    
    stats::sd(
      three_state_model_data_raw_60$
        aerial_duration_min,
      na.rm = TRUE
    ),
    
    stats::sd(
      three_state_model_data_raw_60[[
        hfi_variable_60
      ]],
      na.rm = TRUE
    ),
    
    stats::sd(
      three_state_model_data_raw_60[[
        elevation_variable_60
      ]],
      na.rm = TRUE
    ),
    
    stats::sd(
      three_state_model_data_raw_60[[
        low_vegetation_variable_60
      ]],
      na.rm = TRUE)))

print(standardization_parameters_60,n = Inf)

# 2.2.2 Extract the standardization parameters for each variables ----
age_mean_60 <-
  standardization_parameters_60$
  mean[
    standardization_parameters_60$
      variable ==
      age_variable_60
  ]

age_sd_60 <-
  standardization_parameters_60$
  standard_deviation[
    standardization_parameters_60$
      variable ==
      age_variable_60
  ]

duration_mean_60 <-
  standardization_parameters_60$
  mean[
    standardization_parameters_60$
      variable ==
      "aerial_duration_min"
  ]

duration_sd_60 <-
  standardization_parameters_60$
  standard_deviation[
    standardization_parameters_60$
      variable ==
      "aerial_duration_min"
  ]

hfi_mean_60 <-
  standardization_parameters_60$
  mean[
    standardization_parameters_60$
      variable ==
      hfi_variable_60
  ]

hfi_sd_60 <-
  standardization_parameters_60$
  standard_deviation[
    standardization_parameters_60$
      variable ==
      hfi_variable_60
  ]

elevation_mean_60 <-
  standardization_parameters_60$
  mean[
    standardization_parameters_60$
      variable ==
      elevation_variable_60
  ]

elevation_sd_60 <-
  standardization_parameters_60$
  standard_deviation[
    standardization_parameters_60$
      variable ==
      elevation_variable_60
  ]

low_vegetation_mean_60 <-
  standardization_parameters_60$
  mean[
    standardization_parameters_60$
      variable ==
      low_vegetation_variable_60
  ]

low_vegetation_sd_60 <-
  standardization_parameters_60$
  standard_deviation[
    standardization_parameters_60$
      variable ==
      low_vegetation_variable_60
  ]


# 2.2.3 Applied standaridization for each variables ----
# duration_z2 is the square of the standardized duration, consistent with
# the previous GPS-based model.
three_state_model_data_60 <- three_state_model_data_raw_60 %>%
  dplyr::mutate(
    age_z =
      (
        .data[[
          age_variable_60
        ]] -
          age_mean_60
      ) /
      age_sd_60,
    
    duration_z =
      (
        aerial_duration_min -
          duration_mean_60
      ) /
      duration_sd_60,
    
    duration_z2 =
      duration_z^2,
    
    hfi_mean_1000m_z =
      (
        .data[[
          hfi_variable_60
        ]] -
          hfi_mean_60
      ) /
      hfi_sd_60,
    
    dem_elevation_z =
      (
        .data[[
          elevation_variable_60
        ]] -
          elevation_mean_60
      ) /
      elevation_sd_60,
    
    prop_low_vegetation_5cells_z =
      (
        .data[[
          low_vegetation_variable_60
        ]] -
          low_vegetation_mean_60
      ) /
      low_vegetation_sd_60,
    
    individual_id =
      droplevels(
        factor(
          individual_id
        )
      ),
    
    transition_destination =
      factor(
        transition_destination,
        levels = c(
          "aerial",
          "resting",
          "feeding"
        )
      ),
    
    transition_destination_code =
      as.integer(
        transition_destination
      ) -
      1L
  ) %>%
  tidyr::drop_na(
    transition_destination_code,
    individual_id,
    cos_diel,
    sin_time,
    age_z,
    duration_z,
    duration_z2,
    hfi_mean_1000m_z,
    dem_elevation_z,
    prop_low_vegetation_5cells_z
  ) %>%
  dplyr::filter(
    dplyr::if_all(
      c(
        cos_diel,
        sin_time,
        age_z,
        duration_z,
        duration_z2,
        hfi_mean_1000m_z,
        dem_elevation_z,
        prop_low_vegetation_5cells_z
      ),
      is.finite))

# 2.3 Define the population-level three-state model ----
# 2.3.0 Give an equal weight to each individuals ----
#' Without this weight, eagle with 300 transition contributes more than eagles
#' with only 30 transitions, therefore the system is dominated by some individuals
#' only.
three_state_model_data_60 <- three_state_model_data_60 %>%
  dplyr::group_by(
    individual_id
  ) %>%
  dplyr::mutate(
    n_transitions_individual =
      dplyr::n(),
    
    individual_equal_weight =
      1 /
      n_transitions_individual
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    individual_equal_weight =
      individual_equal_weight /
      mean(
        individual_equal_weight
      )
  )


# 2.3.1 Model formula ----
# Formula 1: resting relative to aerial.
# Formula 2: feeding relative to aerial.
# Each destination has its own coefficients and its own individual random
# intercept, but HFI slopes are initially common among individuals.
three_state_formula_60 <- list(
  transition_destination_code ~
    cos_diel +
    sin_time +
    age_z +
    duration_z +
    duration_z2 +
    hfi_mean_1000m_z +
    dem_elevation_z +
    prop_low_vegetation_5cells_z +
    s(
      individual_id,
      bs = "re"
    ),
  
  ~
    cos_diel +
    sin_time +
    age_z +
    duration_z +
    duration_z2 +
    hfi_mean_1000m_z +
    dem_elevation_z +
    prop_low_vegetation_5cells_z +
    s(
      individual_id,
      bs = "re"
    )
)


# 2.3.2 Fit the population-level three-state model ----
three_state_model_60 <- mgcv::gam(
  formula =
    three_state_formula_60,
  
  family =
    mgcv::multinom(
      K = 2
    ),
  
  data =
    three_state_model_data_60,
  
  method =
    "REML")

three_state_model_summary_60 <- summary(three_state_model_60,re.test = TRUE)
print(three_state_model_summary_60)


# 2.4 Calculate empirical HFI Q05 and Q95 thresholds ----
# For direct comparison with the GPS analysis, replace them later with the
# fixed GPS Q05 and Q95 thresholds transformed using the same HFI scale.
hfi_thresholds_60 <- three_state_model_data_60 %>%
  dplyr::summarise(
    hfi_q05_raw =
      as.numeric(
        stats::quantile(
          .data[[
            hfi_variable_60
          ]],
          probs = 0.05,
          names = FALSE
        )
      ),
    
    hfi_q95_raw =
      as.numeric(
        stats::quantile(
          .data[[
            hfi_variable_60
          ]],
          probs = 0.95,
          names = FALSE
        )
      )
  ) %>%
  dplyr::mutate(
    hfi_q05_z =
      (
        hfi_q05_raw -
          hfi_mean_60
      ) /
      hfi_sd_60,
    
    hfi_q95_z =
      (
        hfi_q95_raw -
          hfi_mean_60
      ) /
      hfi_sd_60,
    
    hfi_q95_q05_difference_z =
      hfi_q95_z -
      hfi_q05_z
  )

print(hfi_thresholds_60)


# 2.5 Extract destination-specific HFI coefficients ----
# The multinomial model contains two parametric HFI coefficients:
# one for resting relative to aerial and one for feeding relative to aerial.
parametric_coefficient_table_60 <- three_state_model_summary_60$
  p.table %>%
  as.data.frame() %>%
  tibble::rownames_to_column(
    var = "term"
  ) %>%
  tibble::as_tibble()

hfi_coefficient_rows_60 <- grep(
  pattern =
    "^hfi_mean_1000m_z",
  
  x =
    parametric_coefficient_table_60$
    term
)

stopifnot(
  length(
    hfi_coefficient_rows_60
  ) == 2L
)

hfi_destination_effects_60 <- parametric_coefficient_table_60 %>%
  dplyr::slice(
    hfi_coefficient_rows_60
  ) %>%
  dplyr::mutate(
    destination_contrast = c(
      "resting versus aerial",
      "feeding versus aerial"
    )
  ) %>%
  dplyr::transmute(
    destination_contrast,
    
    coefficient =
      Estimate,
    
    coefficient_standard_error =
      `Std. Error`,
    
    coefficient_confidence_low =
      coefficient -
      confidence_multiplier_60 *
      coefficient_standard_error,
    
    coefficient_confidence_high =
      coefficient +
      confidence_multiplier_60 *
      coefficient_standard_error,
    
    coefficient_p_value =
      `Pr(>|z|)`
  )

print(hfi_destination_effects_60,n = Inf)


# 2.13 Calculate destination-specific Q05-Q95 log-odds contrasts ----
# A positive value means that increasing HFI from Q05 to Q95 increases
# the log-odds of the destination relative to remaining aerial.
#
# A negative value means that increasing HFI decreases the log-odds of
# the destination relative to remaining aerial.

hfi_q05_q95_log_odds_60 <- hfi_destination_effects_60 %>%
  dplyr::mutate(
    hfi_q05_raw =
      hfi_thresholds_60$
      hfi_q05_raw,
    
    hfi_q95_raw =
      hfi_thresholds_60$
      hfi_q95_raw,
    
    hfi_q05_z =
      hfi_thresholds_60$
      hfi_q05_z,
    
    hfi_q95_z =
      hfi_thresholds_60$
      hfi_q95_z,
    
    hfi_q95_q05_difference_z =
      hfi_thresholds_60$
      hfi_q95_q05_difference_z,
    
    Q95_Q05_log_odds_contrast =
      coefficient *
      hfi_q95_q05_difference_z,
    
    Q95_Q05_log_odds_standard_error =
      coefficient_standard_error *
      abs(
        hfi_q95_q05_difference_z
      ),
    
    Q95_Q05_log_odds_confidence_low =
      Q95_Q05_log_odds_contrast -
      confidence_multiplier_60 *
      Q95_Q05_log_odds_standard_error,
    
    Q95_Q05_log_odds_confidence_high =
      Q95_Q05_log_odds_contrast +
      confidence_multiplier_60 *
      Q95_Q05_log_odds_standard_error,
    
    Q95_Q05_odds_ratio =
      exp(
        Q95_Q05_log_odds_contrast
      ),
    
    Q95_Q05_odds_ratio_confidence_low =
      exp(
        Q95_Q05_log_odds_confidence_low
      ),
    
    Q95_Q05_odds_ratio_confidence_high =
      exp(
        Q95_Q05_log_odds_confidence_high
      )
  ) %>%
  dplyr::select(
    destination_contrast,
    hfi_q05_raw,
    hfi_q95_raw,
    coefficient,
    coefficient_confidence_low,
    coefficient_confidence_high,
    coefficient_p_value,
    Q95_Q05_log_odds_contrast,
    Q95_Q05_log_odds_confidence_low,
    Q95_Q05_log_odds_confidence_high,
    Q95_Q05_odds_ratio,
    Q95_Q05_odds_ratio_confidence_low,
    Q95_Q05_odds_ratio_confidence_high
  )

print(hfi_q05_q95_log_odds_60,n = Inf)

# control hfi distribution per transition types 
hfi_support_by_destination_60 <- three_state_model_data_60 %>%
  dplyr::group_by(
    transition_destination
  ) %>%
  dplyr::summarise(
    n_transitions =
      dplyr::n(),
    
    hfi_min =
      min(
        hfi_mean_1000m,
        na.rm = TRUE
      ),
    
    hfi_q05 =
      stats::quantile(
        hfi_mean_1000m,
        0.05,
        na.rm = TRUE,
        names = FALSE
      ),
    
    hfi_q25 =
      stats::quantile(
        hfi_mean_1000m,
        0.25,
        na.rm = TRUE,
        names = FALSE
      ),
    
    hfi_median =
      stats::median(
        hfi_mean_1000m,
        na.rm = TRUE
      ),
    
    hfi_q75 =
      stats::quantile(
        hfi_mean_1000m,
        0.75,
        na.rm = TRUE,
        names = FALSE
      ),
    
    hfi_q95 =
      stats::quantile(
        hfi_mean_1000m,
        0.95,
        na.rm = TRUE,
        names = FALSE
      ),
    
    hfi_max =
      max(
        hfi_mean_1000m,
        na.rm = TRUE
      ),
    
    .groups = "drop"
  )

print(
  hfi_support_by_destination_60,
  n = Inf
)


#------------------------------------------------------------------------------
# 2.14 Calculate aerial-persistence log-odds from the three-state model ----

# Create prediction datasets at low and high HFI.
prediction_data_q05_60 <- three_state_model_data_60 %>%
  dplyr::mutate(
    hfi_mean_1000m_z =
      hfi_thresholds_60$hfi_q05_z
  )

prediction_data_q95_60 <- three_state_model_data_60 %>%
  dplyr::mutate(
    hfi_mean_1000m_z =
      hfi_thresholds_60$hfi_q95_z
  )


# Predict the three destination probabilities.
probabilities_q05_60 <- stats::predict(
  three_state_model_60,
  newdata =
    prediction_data_q05_60,
  type =
    "response"
)

probabilities_q95_60 <- stats::predict(
  three_state_model_60,
  newdata =
    prediction_data_q95_60,
  type =
    "response"
)

colnames(
  probabilities_q05_60
) <- c(
  "aerial",
  "resting",
  "feeding"
)

colnames(
  probabilities_q95_60
) <- c(
  "aerial",
  "resting",
  "feeding"
)


# Calculate mean predicted probabilities for each individual.
individual_aerial_persistence_60 <- tibble::tibble(
  individual_id =
    three_state_model_data_60$
    individual_id,
  
  aerial_q05 =
    probabilities_q05_60[
      ,
      "aerial"
    ],
  
  resting_q05 =
    probabilities_q05_60[
      ,
      "resting"
    ],
  
  feeding_q05 =
    probabilities_q05_60[
      ,
      "feeding"
    ],
  
  aerial_q95 =
    probabilities_q95_60[
      ,
      "aerial"
    ],
  
  resting_q95 =
    probabilities_q95_60[
      ,
      "resting"
    ],
  
  feeding_q95 =
    probabilities_q95_60[
      ,
      "feeding"
    ]
) %>%
  dplyr::group_by(
    individual_id
  ) %>%
  dplyr::summarise(
    dplyr::across(
      c(
        aerial_q05,
        resting_q05,
        feeding_q05,
        aerial_q95,
        resting_q95,
        feeding_q95
      ),
      mean
    ),
    
    .groups =
      "drop"
  )


# Give every individual the same contribution.
aerial_persistence_log_odds_60 <-
  individual_aerial_persistence_60 %>%
  dplyr::summarise(
    hfi_q05 =
      hfi_thresholds_60$
      hfi_q05_raw,
    
    hfi_q95 =
      hfi_thresholds_60$
      hfi_q95_raw,
    
    probability_aerial_q05 =
      mean(
        aerial_q05
      ),
    
    probability_terrestrial_q05 =
      mean(
        resting_q05 +
          feeding_q05
      ),
    
    probability_aerial_q95 =
      mean(
        aerial_q95
      ),
    
    probability_terrestrial_q95 =
      mean(
        resting_q95 +
          feeding_q95
      ),
    
    aerial_log_odds_q05 =
      log(
        probability_aerial_q05 /
          probability_terrestrial_q05
      ),
    
    aerial_log_odds_q95 =
      log(
        probability_aerial_q95 /
          probability_terrestrial_q95
      ),
    
    Q95_Q05_aerial_log_odds_contrast =
      aerial_log_odds_q95 -
      aerial_log_odds_q05,
    
    Q95_Q05_aerial_odds_ratio =
      exp(
        Q95_Q05_aerial_log_odds_contrast
      ),
    
    absolute_aerial_probability_difference =
      probability_aerial_q95 -
      probability_aerial_q05,
    
    absolute_aerial_difference_percentage_points =
      100 *
      absolute_aerial_probability_difference
  )

print(
  aerial_persistence_log_odds_60,
  n = Inf
)

#------------------------------------------------------------------------------
# 2.15 Inspect individual representation along the high-HFI gradient ----
#' **Steps:**
#' (i) calculate common global HFI thresholds;
#' (ii) count transitions and individuals above each threshold by destination;
#' (iii) inspect the contribution of each individual to high-HFI transitions.


# 2.15.1 Calculate common HFI thresholds ----

high_hfi_thresholds_60 <- stats::quantile(
  three_state_model_data_60$hfi_mean_1000m,
  probs = c(
    0.80,
    0.90,
    0.95
  ),
  na.rm = TRUE,
  names = FALSE
)

names(
  high_hfi_thresholds_60
) <- c(
  "Q80",
  "Q90",
  "Q95"
)

print(
  high_hfi_thresholds_60
)


# 2.15.2 Create one row per transition and threshold ----

high_hfi_transition_data_60 <- three_state_model_data_60 %>%
  dplyr::select(
    individual_id,
    transition_destination,
    hfi_mean_1000m
  ) %>%
  tidyr::crossing(
    hfi_threshold = names(
      high_hfi_thresholds_60
    )
  ) %>%
  dplyr::mutate(
    threshold_value =
      high_hfi_thresholds_60[
        hfi_threshold
      ],
    
    above_threshold =
      hfi_mean_1000m >=
      threshold_value
  )


# 2.15.3 Summarise high-HFI representation by destination ----

high_hfi_support_by_destination_60 <-
  high_hfi_transition_data_60 %>%
  dplyr::group_by(
    hfi_threshold,
    threshold_value,
    transition_destination
  ) %>%
  dplyr::summarise(
    n_transitions_total =
      dplyr::n(),
    
    n_transitions_high_hfi =
      sum(
        above_threshold
      ),
    
    proportion_transitions_high_hfi =
      mean(
        above_threshold
      ),
    
    n_individuals_total =
      dplyr::n_distinct(
        individual_id
      ),
    
    n_individuals_high_hfi =
      dplyr::n_distinct(
        individual_id[
          above_threshold
        ]
      ),
    
    proportion_individuals_high_hfi =
      n_individuals_high_hfi /
      n_individuals_total,
    
    .groups =
      "drop"
  )

print(
  high_hfi_support_by_destination_60,
  n = Inf
)


# 2.15.4 Count high-HFI transitions per individual and destination ----

high_hfi_support_by_individual_60 <-
  high_hfi_transition_data_60 %>%
  dplyr::filter(
    above_threshold
  ) %>%
  dplyr::count(
    hfi_threshold,
    threshold_value,
    transition_destination,
    individual_id,
    name =
      "n_high_hfi_transitions"
  ) %>%
  dplyr::arrange(
    hfi_threshold,
    transition_destination,
    dplyr::desc(
      n_high_hfi_transitions
    )
  )

print(
  high_hfi_support_by_individual_60,
  n = Inf
)


# 2.15.5 Determine whether high-HFI support is concentrated in few individuals ----

high_hfi_individual_concentration_60 <-
  high_hfi_support_by_individual_60 %>%
  dplyr::group_by(
    hfi_threshold,
    threshold_value,
    transition_destination
  ) %>%
  dplyr::summarise(
    n_individuals_with_at_least_1 =
      dplyr::n(),
    
    n_individuals_with_at_least_3 =
      sum(
        n_high_hfi_transitions >= 3L
      ),
    
    n_individuals_with_at_least_5 =
      sum(
        n_high_hfi_transitions >= 5L
      ),
    
    maximum_individual_contribution =
      max(
        n_high_hfi_transitions
      ),
    
    proportion_from_top_individual =
      max(
        n_high_hfi_transitions
      ) /
      sum(
        n_high_hfi_transitions
      ),
    
    proportion_from_top_5_individuals =
      sum(
        sort(
          n_high_hfi_transitions,
          decreasing = TRUE
        )[
          seq_len(
            min(
              5L,
              dplyr::n()
            )
          )
        ]
      ) /
      sum(
        n_high_hfi_transitions
      ),
    
    .groups =
      "drop"
  )

print(
  high_hfi_individual_concentration_60,
  n = Inf
)

#------------------------------------------------------------------------------- STEP 3 : model fitting with a slope response to each individuals

#------------------------------------------------------------------------------
# STEP 3: test individual heterogeneity in HFI responses ----
#' **Philosophy:**
#' Test whether the within-individual HFI response varies among individuals
#' separately for:
#'   (i) resting versus aerial;
#'   (ii) feeding versus aerial.
#'
#' **Steps:**
#' (i) separate within- and between-individual HFI variation;
#' (ii) compare common- and individual-slope models;
#' (iii) test random-slope variance for resting and feeding;
#' (iv) extract individual HFI slopes and uncertainty;
#' (v) calculate individual Q05-Q95 probability changes and uncertainty;
#' (vi) compile ordered individual response tables.


# 3.1 Parameters ----

n_parameter_simulations_60 <- 500L
parameter_simulation_seed_60 <- 20260723L


# 3.2 Separate within- and between-individual HFI variation ----
# hfi_between_z:
#   mean HFI exposure of each individual.
#
# hfi_within_z:
#   deviation of each observation from that individual's mean HFI.
#
# The individual random slopes are estimated for hfi_within_z.

three_state_random_slope_data_60 <- three_state_model_data_60 %>%
  dplyr::group_by(
    individual_id
  ) %>%
  dplyr::mutate(
    hfi_between_z =
      mean(
        hfi_mean_1000m_z
      ),
    
    hfi_within_z =
      hfi_mean_1000m_z -
      hfi_between_z
  ) %>%
  dplyr::ungroup()


# 3.3 Summarise empirical support by individual ----

individual_model_support_60 <- three_state_random_slope_data_60 %>%
  dplyr::group_by(
    individual_id
  ) %>%
  dplyr::summarise(
    n_transitions =
      dplyr::n(),
    
    n_aerial =
      sum(
        transition_destination == "aerial"
      ),
    
    n_resting =
      sum(
        transition_destination == "resting"
      ),
    
    n_feeding =
      sum(
        transition_destination == "feeding"
      ),
    
    n_bursts =
      dplyr::n_distinct(
        burst_id
      ),
    
    hfi_min =
      min(
        hfi_mean_1000m
      ),
    
    hfi_max =
      max(
        hfi_mean_1000m
      ),
    
    hfi_range =
      hfi_max -
      hfi_min,
    
    hfi_within_sd_z =
      stats::sd(
        hfi_within_z
      ),
    
    .groups =
      "drop"
  )

print(
  individual_model_support_60,
  n = Inf
)


# 3.4 Define common- and individual-slope models ----
# Formula 1 always describes resting versus aerial.
# Formula 2 always describes feeding versus aerial.

three_state_slope_formulas_60 <- list(
  
  # Common within-individual HFI slopes
  common_HFI_slopes =
    list(
      transition_destination_code ~
        cos_diel +
        sin_time +
        age_z +
        duration_z +
        duration_z2 +
        hfi_within_z +
        hfi_between_z +
        dem_elevation_z +
        prop_low_vegetation_5cells_z +
        s(
          individual_id,
          bs = "re"
        ),
      
      ~
        cos_diel +
        sin_time +
        age_z +
        duration_z +
        duration_z2 +
        hfi_within_z +
        hfi_between_z +
        dem_elevation_z +
        prop_low_vegetation_5cells_z +
        s(
          individual_id,
          bs = "re"
        )
    ),
  
  # Individual HFI slopes for resting only
  individual_resting_slope =
    list(
      transition_destination_code ~
        cos_diel +
        sin_time +
        age_z +
        duration_z +
        duration_z2 +
        hfi_within_z +
        hfi_between_z +
        dem_elevation_z +
        prop_low_vegetation_5cells_z +
        s(
          individual_id,
          bs = "re"
        ) +
        s(
          hfi_within_z,
          individual_id,
          bs = "re"
        ),
      
      ~
        cos_diel +
        sin_time +
        age_z +
        duration_z +
        duration_z2 +
        hfi_within_z +
        hfi_between_z +
        dem_elevation_z +
        prop_low_vegetation_5cells_z +
        s(
          individual_id,
          bs = "re"
        )
    ),
  
  # Individual HFI slopes for feeding only
  individual_feeding_slope =
    list(
      transition_destination_code ~
        cos_diel +
        sin_time +
        age_z +
        duration_z +
        duration_z2 +
        hfi_within_z +
        hfi_between_z +
        dem_elevation_z +
        prop_low_vegetation_5cells_z +
        s(
          individual_id,
          bs = "re"
        ),
      
      ~
        cos_diel +
        sin_time +
        age_z +
        duration_z +
        duration_z2 +
        hfi_within_z +
        hfi_between_z +
        dem_elevation_z +
        prop_low_vegetation_5cells_z +
        s(
          individual_id,
          bs = "re"
        ) +
        s(
          hfi_within_z,
          individual_id,
          bs = "re"
        )
    ),
  
  # Individual HFI slopes for both destinations
  individual_resting_and_feeding_slopes =
    list(
      transition_destination_code ~
        cos_diel +
        sin_time +
        age_z +
        duration_z +
        duration_z2 +
        hfi_within_z +
        hfi_between_z +
        dem_elevation_z +
        prop_low_vegetation_5cells_z +
        s(
          individual_id,
          bs = "re"
        ) +
        s(
          hfi_within_z,
          individual_id,
          bs = "re"
        ),
      
      ~
        cos_diel +
        sin_time +
        age_z +
        duration_z +
        duration_z2 +
        hfi_within_z +
        hfi_between_z +
        dem_elevation_z +
        prop_low_vegetation_5cells_z +
        s(
          individual_id,
          bs = "re"
        ) +
        s(
          hfi_within_z,
          individual_id,
          bs = "re"
        )
    )
)


# 3.5 Fit the four candidate models ----

three_state_slope_models_60 <- lapply(
  three_state_slope_formulas_60,
  function(model_formula) {
    
    mgcv::gam(
      formula =
        model_formula,
      
      family =
        mgcv::multinom(
          K = 2
        ),
      
      data =
        three_state_random_slope_data_60,
      
      method =
        "REML"
    )
  }
)


three_state_common_slope_model_60 <-
  three_state_slope_models_60$
  common_HFI_slopes

three_state_resting_slope_model_60 <-
  three_state_slope_models_60$
  individual_resting_slope

three_state_feeding_slope_model_60 <-
  three_state_slope_models_60$
  individual_feeding_slope

three_state_random_slope_model_60 <-
  three_state_slope_models_60$
  individual_resting_and_feeding_slopes


# 3.6 Compare common- and individual-slope models ----
# 3.6 Compare common- and individual-slope models ----

three_state_slope_model_comparison_60 <- tibble::tibble(
  model =
    names(
      three_state_slope_models_60
    ),
  
  AIC =
    vapply(
      three_state_slope_models_60,
      stats::AIC,
      numeric(1)
    ),
  
  convergence_status =
    vapply(
      three_state_slope_models_60,
      function(model) {
        model$outer.info$conv
      },
      character(1)
    ),
  
  converged =
    vapply(
      three_state_slope_models_60,
      function(model) {
        identical(
          model$outer.info$conv,
          "full convergence"
        )
      },
      logical(1)
    ),
  
  maximum_absolute_gradient =
    vapply(
      three_state_slope_models_60,
      function(model) {
        max(
          abs(
            model$outer.info$grad
          )
        )
      },
      numeric(1)
    )
) %>%
  dplyr::mutate(
    delta_AIC =
      AIC -
      min(
        AIC
      ),
    
    AIC_change_from_common =
      AIC -
      AIC[
        model ==
          "common_HFI_slopes"
      ]
  ) %>%
  dplyr::arrange(
    AIC
  )

print(
  three_state_slope_model_comparison_60,
  n = Inf
)
# model                                   AIC convergence_status converged maximum_absolute_gradient delta_AIC AIC_change_from_common
# <chr>                                 <dbl> <chr>              <lgl>                         <dbl>     <dbl>                  <dbl>
#   1 individual_resting_slope              8945. full convergence   TRUE                        0.00307     0                     -1.96 
# 2 individual_resting_and_feeding_slopes 8945. full convergence   TRUE                        0.00159     0.137                 -1.82 
# 3 common_HFI_slopes                     8947. full convergence   TRUE                        0.00442     1.96                   0    
# 4 individual_feeding_slope              8947. full convergence   TRUE                        0.00481     2.29                   0.324
#' This shows that the response to feeding is common accross individuals because integrating individual variation in the model does not
#' imporove the model performance. Their is slight improvement for the model when we adjust the slope of resting to individual but it is 
#' not highly significant (0.324)

# 3.7 Test resting and feeding random-slope variances ----

three_state_random_slope_summary_60 <- summary(
  three_state_random_slope_model_60,
  re.test = TRUE
)

random_slope_smooth_table_60 <-
  three_state_random_slope_summary_60$
  s.table %>%
  as.data.frame() %>%
  tibble::rownames_to_column(
    var = "term"
  ) %>%
  tibble::as_tibble() %>%
  dplyr::filter(
    grepl(
      "hfi_within_z.*individual_id",
      term
    )
  )

stopifnot(
  nrow(
    random_slope_smooth_table_60
  ) == 2L
)

random_slope_variance_tests_60 <-
  random_slope_smooth_table_60 %>%
  dplyr::mutate(
    destination_contrast = c(
      "resting versus aerial",
      "feeding versus aerial"
    )
  ) %>%
  dplyr::transmute(
    destination_contrast,
    
    term,
    
    effective_degrees_of_freedom =
      .data[[
        names(
          random_slope_smooth_table_60
        )[2]
      ]],
    
    reference_degrees_of_freedom =
      .data[[
        names(
          random_slope_smooth_table_60
        )[3]
      ]],
    
    test_statistic =
      .data[[
        names(
          random_slope_smooth_table_60
        )[4]
      ]],
    
    p_value =
      .data[[
        names(
          random_slope_smooth_table_60
        )[
          ncol(
            random_slope_smooth_table_60
          )
        ]
      ]]
  )

print(
  random_slope_variance_tests_60,
  n = Inf
)

# 
# destination_contrast  term                            effective_degrees_of_freedom reference_degrees_of_freedom test_statistic p_value
# <chr>                 <chr>                                                  <dbl>                        <dbl>          <dbl>   <dbl>
#   1 resting versus aerial s(hfi_within_z,individual_id)                          12.4                            57           16.9  0.0521
# 2 feeding versus aerial s.1(hfi_within_z,individual_id)                         9.62                           57           12.6  0.0884
#' il y a un retrecissement partiel de la variation individuelle autour de feeding. 

# 3.8 Extract random-slope standard deviations ----
# 3.8 Extract random-slope standard deviations ----

random_slope_variance_components_raw_60 <-
  mgcv::gam.vcomp(
    three_state_random_slope_model_60
  )


# Inspect the object returned by gam.vcomp().
str(
  random_slope_variance_components_raw_60
)


# Convert one gam.vcomp() component into a consistent table.
extract_variance_component_60 <- function(
    variance_component,
    destination_contrast
) {
  
  if (
    is.matrix(
      variance_component
    ) ||
    is.data.frame(
      variance_component
    )
  ) {
    
    variance_component <-
      as.matrix(
        variance_component
      )
    
    tibble::tibble(
      destination_contrast =
        destination_contrast,
      
      term =
        rownames(
          variance_component
        ),
      
      standard_deviation =
        as.numeric(
          variance_component[
            ,
            1
          ]
        ),
      
      confidence_low =
        if (
          ncol(
            variance_component
          ) >= 2L
        ) {
          as.numeric(
            variance_component[
              ,
              2
            ]
          )
        } else {
          rep(
            NA_real_,
            nrow(
              variance_component
            )
          )
        },
      
      confidence_high =
        if (
          ncol(
            variance_component
          ) >= 3L
        ) {
          as.numeric(
            variance_component[
              ,
              3
            ]
          )
        } else {
          rep(
            NA_real_,
            nrow(
              variance_component
            )
          )
        }
    )
    
  } else {
    
    # Case where one component is a named numeric vector.
    variance_component <-
      unlist(
        variance_component,
        recursive = TRUE,
        use.names = TRUE
      )
    
    tibble::tibble(
      destination_contrast =
        destination_contrast,
      
      term =
        names(
          variance_component
        ),
      
      standard_deviation =
        as.numeric(
          variance_component
        ),
      
      confidence_low =
        NA_real_,
      
      confidence_high =
        NA_real_
    )
  }
}


# gam.vcomp() returns one component for each multinomial linear predictor:
# first = resting versus aerial;
# second = feeding versus aerial.

random_slope_variance_components_60 <-
  dplyr::bind_rows(
    extract_variance_component_60(
      random_slope_variance_components_raw_60[[1]],
      "resting versus aerial"
    ),
    
    extract_variance_component_60(
      random_slope_variance_components_raw_60[[2]],
      "feeding versus aerial"
    )
  ) %>%
  dplyr::filter(
    grepl(
      "hfi_within_z.*individual_id",
      term
    )
  ) %>%
  dplyr::select(
    destination_contrast,
    term,
    standard_deviation,
    confidence_low,
    confidence_high
  )


print(
  random_slope_variance_components_60,
  n = Inf
)




# 3.9 Draw parameter vectors for uncertainty estimation ----
# These simulations use the coefficient covariance matrix corrected for
# smoothing-parameter uncertainty when this correction is available.

set.seed(
  parameter_simulation_seed_60
)

model_coefficients_60 <- stats::coef(
  three_state_random_slope_model_60
)

model_covariance_60 <- stats::vcov(
  three_state_random_slope_model_60,
  unconditional = TRUE
)

coefficient_draws_60 <- mgcv::rmvn(
  n =
    n_parameter_simulations_60,
  
  mu =
    model_coefficients_60,
  
  V =
    model_covariance_60
)


# 3.10 Extract individual HFI slopes ----
# The difference between hfi_within_z = 1 and hfi_within_z = 0 isolates:
# population within-HFI slope + individual random-slope deviation.

individual_slope_template_60 <-
  three_state_random_slope_data_60 %>%
  dplyr::arrange(
    individual_id
  ) %>%
  dplyr::group_by(
    individual_id
  ) %>%
  dplyr::slice(
    1
  ) %>%
  dplyr::ungroup()

individual_slope_data_0_60 <-
  individual_slope_template_60 %>%
  dplyr::mutate(
    hfi_within_z = 0
  )

individual_slope_data_1_60 <-
  individual_slope_template_60 %>%
  dplyr::mutate(
    hfi_within_z = 1
  )

slope_matrix_0_60 <- stats::predict(
  three_state_random_slope_model_60,
  newdata =
    individual_slope_data_0_60,
  type =
    "lpmatrix"
)

slope_matrix_1_60 <- stats::predict(
  three_state_random_slope_model_60,
  newdata =
    individual_slope_data_1_60,
  type =
    "lpmatrix"
)

linear_predictor_indices_60 <- attr(
  slope_matrix_0_60,
  "lpi"
)

stopifnot(
  length(
    linear_predictor_indices_60
  ) == 2L
)

slope_difference_matrix_60 <-
  slope_matrix_1_60 -
  slope_matrix_0_60

resting_slope_matrix_60 <-
  slope_difference_matrix_60 *
  0

resting_slope_matrix_60[
  ,
  linear_predictor_indices_60[[1]]
] <-
  slope_difference_matrix_60[
    ,
    linear_predictor_indices_60[[1]]
  ]

feeding_slope_matrix_60 <-
  slope_difference_matrix_60 *
  0

feeding_slope_matrix_60[
  ,
  linear_predictor_indices_60[[2]]
] <-
  slope_difference_matrix_60[
    ,
    linear_predictor_indices_60[[2]]
  ]

resting_slope_draws_60 <-
  resting_slope_matrix_60 %*%
  t(
    coefficient_draws_60
  )

feeding_slope_draws_60 <-
  feeding_slope_matrix_60 %*%
  t(
    coefficient_draws_60
  )

resting_slope_confidence_60 <-
  t(
    apply(
      resting_slope_draws_60,
      1,
      stats::quantile,
      probs = c(
        0.025,
        0.975
      ),
      names = FALSE
    )
  )

feeding_slope_confidence_60 <-
  t(
    apply(
      feeding_slope_draws_60,
      1,
      stats::quantile,
      probs = c(
        0.025,
        0.975
      ),
      names = FALSE
    )
  )


individual_hfi_slopes_60 <- tibble::tibble(
  individual_id =
    as.character(
      individual_slope_template_60$
        individual_id
    ),
  
  resting_HFI_slope =
    as.numeric(
      resting_slope_matrix_60 %*%
        model_coefficients_60
    ),
  
  resting_HFI_slope_confidence_low =
    resting_slope_confidence_60[
      ,
      1
    ],
  
  resting_HFI_slope_confidence_high =
    resting_slope_confidence_60[
      ,
      2
    ],
  
  feeding_HFI_slope =
    as.numeric(
      feeding_slope_matrix_60 %*%
        model_coefficients_60
    ),
  
  feeding_HFI_slope_confidence_low =
    feeding_slope_confidence_60[
      ,
      1
    ],
  
  feeding_HFI_slope_confidence_high =
    feeding_slope_confidence_60[
      ,
      2
    ]
) %>%
  dplyr::mutate(
    resting_slope_sign =
      dplyr::case_when(
        resting_HFI_slope < 0 ~ "negative",
        resting_HFI_slope > 0 ~ "positive",
        TRUE ~ "zero"
      ),
    
    feeding_slope_sign =
      dplyr::case_when(
        feeding_HFI_slope < 0 ~ "negative",
        feeding_HFI_slope > 0 ~ "positive",
        TRUE ~ "zero"
      ),
    
    resting_supported_direction =
      dplyr::case_when(
        resting_HFI_slope_confidence_low > 0 ~
          "positive",
        
        resting_HFI_slope_confidence_high < 0 ~
          "negative",
        
        TRUE ~
          "uncertain"
      ),
    
    feeding_supported_direction =
      dplyr::case_when(
        feeding_HFI_slope_confidence_low > 0 ~
          "positive",
        
        feeding_HFI_slope_confidence_high < 0 ~
          "negative",
        
        TRUE ~
          "uncertain"
      )
  )


# 3.11 Calculate individual Q05-Q95 probability changes ----
# The same global Q05 and Q95 HFI values are used for all individuals.
# hfi_between_z remains fixed at each individual's mean exposure.

prediction_data_q05_60 <-
  three_state_random_slope_data_60 %>%
  dplyr::mutate(
    hfi_within_z =
      hfi_thresholds_60$
      hfi_q05_z -
      hfi_between_z
  )

prediction_data_q95_60 <-
  three_state_random_slope_data_60 %>%
  dplyr::mutate(
    hfi_within_z =
      hfi_thresholds_60$
      hfi_q95_z -
      hfi_between_z
  )

probability_q05_60 <- stats::predict(
  three_state_random_slope_model_60,
  newdata =
    prediction_data_q05_60,
  type =
    "response"
)

probability_q95_60 <- stats::predict(
  three_state_random_slope_model_60,
  newdata =
    prediction_data_q95_60,
  type =
    "response"
)

colnames(
  probability_q05_60
) <- c(
  "aerial",
  "resting",
  "feeding"
)

colnames(
  probability_q95_60
) <- c(
  "aerial",
  "resting",
  "feeding"
)


individual_probability_changes_point_60 <- tibble::tibble(
  individual_id =
    as.character(
      three_state_random_slope_data_60$
        individual_id
    ),
  
  aerial_q05 =
    probability_q05_60[
      ,
      "aerial"
    ],
  
  resting_q05 =
    probability_q05_60[
      ,
      "resting"
    ],
  
  feeding_q05 =
    probability_q05_60[
      ,
      "feeding"
    ],
  
  aerial_q95 =
    probability_q95_60[
      ,
      "aerial"
    ],
  
  resting_q95 =
    probability_q95_60[
      ,
      "resting"
    ],
  
  feeding_q95 =
    probability_q95_60[
      ,
      "feeding"
    ]
) %>%
  dplyr::group_by(
    individual_id
  ) %>%
  dplyr::summarise(
    dplyr::across(
      c(
        aerial_q05,
        resting_q05,
        feeding_q05,
        aerial_q95,
        resting_q95,
        feeding_q95
      ),
      mean
    ),
    
    .groups =
      "drop"
  ) %>%
  dplyr::mutate(
    aerial_probability_change =
      aerial_q95 -
      aerial_q05,
    
    resting_probability_change =
      resting_q95 -
      resting_q05,
    
    feeding_probability_change =
      feeding_q95 -
      feeding_q05
  )


# 3.12 Simulate uncertainty for individual probability changes ----

prediction_matrix_q05_60 <- stats::predict(
  three_state_random_slope_model_60,
  newdata =
    prediction_data_q05_60,
  type =
    "lpmatrix"
)

prediction_matrix_q95_60 <- stats::predict(
  three_state_random_slope_model_60,
  newdata =
    prediction_data_q95_60,
  type =
    "lpmatrix"
)

linear_predictor_indices_60 <- attr(
  prediction_matrix_q05_60,
  "lpi"
)

resting_eta_q05_60 <-
  prediction_matrix_q05_60[
    ,
    linear_predictor_indices_60[[1]],
    drop = FALSE
  ] %*%
  t(
    coefficient_draws_60[
      ,
      linear_predictor_indices_60[[1]],
      drop = FALSE
    ]
  )

feeding_eta_q05_60 <-
  prediction_matrix_q05_60[
    ,
    linear_predictor_indices_60[[2]],
    drop = FALSE
  ] %*%
  t(
    coefficient_draws_60[
      ,
      linear_predictor_indices_60[[2]],
      drop = FALSE
    ]
  )

resting_eta_q95_60 <-
  prediction_matrix_q95_60[
    ,
    linear_predictor_indices_60[[1]],
    drop = FALSE
  ] %*%
  t(
    coefficient_draws_60[
      ,
      linear_predictor_indices_60[[1]],
      drop = FALSE
    ]
  )

feeding_eta_q95_60 <-
  prediction_matrix_q95_60[
    ,
    linear_predictor_indices_60[[2]],
    drop = FALSE
  ] %*%
  t(
    coefficient_draws_60[
      ,
      linear_predictor_indices_60[[2]],
      drop = FALSE
    ]
  )


multinomial_probabilities_60 <- function(
    resting_eta,
    feeding_eta
) {
  
  maximum_eta <-
    pmax(
      0,
      resting_eta,
      feeding_eta
    )
  
  denominator <-
    exp(
      -maximum_eta
    ) +
    exp(
      resting_eta -
        maximum_eta
    ) +
    exp(
      feeding_eta -
        maximum_eta
    )
  
  list(
    aerial =
      exp(
        -maximum_eta
      ) /
      denominator,
    
    resting =
      exp(
        resting_eta -
          maximum_eta
      ) /
      denominator,
    
    feeding =
      exp(
        feeding_eta -
          maximum_eta
      ) /
      denominator
  )
}


simulated_probabilities_q05_60 <-
  multinomial_probabilities_60(
    resting_eta =
      resting_eta_q05_60,
    
    feeding_eta =
      feeding_eta_q05_60
  )

simulated_probabilities_q95_60 <-
  multinomial_probabilities_60(
    resting_eta =
      resting_eta_q95_60,
    
    feeding_eta =
      feeding_eta_q95_60
  )


individual_index_60 <- as.integer(
  three_state_random_slope_data_60$
    individual_id
)

individual_sizes_60 <- tabulate(
  individual_index_60
)

mean_matrix_by_individual_60 <- function(
    values
) {
  
  sweep(
    rowsum(
      values,
      individual_index_60,
      reorder = FALSE
    ),
    MARGIN = 1,
    STATS =
      individual_sizes_60,
    FUN = "/"
  )
}


aerial_probability_change_draws_60 <-
  mean_matrix_by_individual_60(
    simulated_probabilities_q95_60$
      aerial -
      simulated_probabilities_q05_60$
      aerial
  )

resting_probability_change_draws_60 <-
  mean_matrix_by_individual_60(
    simulated_probabilities_q95_60$
      resting -
      simulated_probabilities_q05_60$
      resting
  )

feeding_probability_change_draws_60 <-
  mean_matrix_by_individual_60(
    simulated_probabilities_q95_60$
      feeding -
      simulated_probabilities_q05_60$
      feeding
  )


aerial_probability_change_confidence_60 <-
  t(
    apply(
      aerial_probability_change_draws_60,
      1,
      stats::quantile,
      probs = c(
        0.025,
        0.975
      ),
      names = FALSE
    )
  )

resting_probability_change_confidence_60 <-
  t(
    apply(
      resting_probability_change_draws_60,
      1,
      stats::quantile,
      probs = c(
        0.025,
        0.975
      ),
      names = FALSE
    )
  )

feeding_probability_change_confidence_60 <-
  t(
    apply(
      feeding_probability_change_draws_60,
      1,
      stats::quantile,
      probs = c(
        0.025,
        0.975
      ),
      names = FALSE
    )
  )


individual_probability_change_uncertainty_60 <- tibble::tibble(
  individual_id =
    levels(
      three_state_random_slope_data_60$
        individual_id
    ),
  
  aerial_probability_change_confidence_low =
    aerial_probability_change_confidence_60[
      ,
      1
    ],
  
  aerial_probability_change_confidence_high =
    aerial_probability_change_confidence_60[
      ,
      2
    ],
  
  resting_probability_change_confidence_low =
    resting_probability_change_confidence_60[
      ,
      1
    ],
  
  resting_probability_change_confidence_high =
    resting_probability_change_confidence_60[
      ,
      2
    ],
  
  feeding_probability_change_confidence_low =
    feeding_probability_change_confidence_60[
      ,
      1
    ],
  
  feeding_probability_change_confidence_high =
    feeding_probability_change_confidence_60[
      ,
      2
    ]
)


# 3.13 Compile the complete individual response table ----

individual_hfi_response_table_60 <-
  individual_hfi_slopes_60 %>%
  dplyr::left_join(
    individual_probability_changes_point_60,
    by =
      "individual_id"
  ) %>%
  dplyr::left_join(
    individual_probability_change_uncertainty_60,
    by =
      "individual_id"
  ) %>%
  dplyr::left_join(
    individual_model_support_60 %>%
      dplyr::mutate(
        individual_id =
          as.character(
            individual_id
          )
      ),
    by =
      "individual_id"
  ) %>%
  dplyr::arrange(
    resting_HFI_slope,
    feeding_HFI_slope
  )

print(
  individual_hfi_response_table_60,
  n = Inf
)


# 3.14 Compile the requested ordered sign table ----
# Individuals are ordered from the most negative to the most positive
# resting slope. The feeding slope and sign are shown in parallel.

individual_hfi_slope_sign_table_60 <-
  individual_hfi_response_table_60 %>%
  dplyr::select(
    individual_id,
    
    resting_HFI_slope,
    resting_slope_sign,
    resting_supported_direction,
    
    feeding_HFI_slope,
    feeding_slope_sign,
    feeding_supported_direction
  ) %>%
  dplyr::arrange(
    resting_HFI_slope,
    feeding_HFI_slope
  )

print(
  individual_hfi_slope_sign_table_60,
  n = Inf
)


# 3.15 Summarise combinations of estimated slope signs ----

individual_slope_sign_distribution_60 <-
  individual_hfi_slope_sign_table_60 %>%
  dplyr::count(
    resting_slope_sign,
    feeding_slope_sign,
    name =
      "n_individuals"
  ) %>%
  dplyr::mutate(
    proportion_individuals =
      n_individuals /
      sum(
        n_individuals
      )
  ) %>%
  dplyr::arrange(
    resting_slope_sign,
    feeding_slope_sign
  )

print(
  individual_slope_sign_distribution_60,
  n = Inf
)


# 3.16 Calculate equal-weight population summaries ----
# Every individual contributes exactly one value to these summaries.

equal_weight_individual_response_summary_60 <-
  individual_hfi_response_table_60 %>%
  dplyr::summarise(
    n_individuals =
      dplyr::n(),
    
    mean_aerial_probability_change =
      mean(
        aerial_probability_change
      ),
    
    mean_resting_probability_change =
      mean(
        resting_probability_change
      ),
    
    mean_feeding_probability_change =
      mean(
        feeding_probability_change
      ),
    
    median_resting_HFI_slope =
      stats::median(
        resting_HFI_slope
      ),
    
    median_feeding_HFI_slope =
      stats::median(
        feeding_HFI_slope
      ),
    
    proportion_negative_resting_slopes =
      mean(
        resting_HFI_slope < 0
      ),
    
    proportion_negative_feeding_slopes =
      mean(
        feeding_HFI_slope < 0
      )
  )

print(
  equal_weight_individual_response_summary_60,
  n = Inf
)

# pente variable pour resting seulement 
selected_three_state_slope_model_60 <-
  three_state_resting_slope_model_60

selected_three_state_slope_summary_60 <-
  summary(
    selected_three_state_slope_model_60,
    re.test = TRUE
  )

print(
  selected_three_state_slope_summary_60
)
#------------------------------------------------------------------------------
# Quantify individual intensity and between-individual heterogeneity ----

selected_resting_slope_model_60 <-
  three_state_resting_slope_model_60


# Population mean HFI slope for resting versus aerial.
population_resting_HFI_slope_60 <- unname(
  stats::coef(
    selected_resting_slope_model_60
  )[["hfi_within_z"]]
)


# Extract variance components from the first multinomial predictor:
# resting versus aerial.
resting_variance_component_raw_60 <-
  mgcv::gam.vcomp(
    selected_resting_slope_model_60
  )[[1]]


resting_variance_component_table_60 <-
  resting_variance_component_raw_60 %>%
  as.data.frame() %>%
  tibble::rownames_to_column(
    var = "term"
  )


resting_HFI_slope_SD_60 <-
  resting_variance_component_table_60 %>%
  dplyr::filter(
    grepl(
      "hfi_within_z.*individual_id",
      term
    )
  ) %>%
  dplyr::pull(
    2
  )


resting_HFI_slope_heterogeneity_60 <- tibble::tibble(
  population_mean_slope =
    population_resting_HFI_slope_60,
  
  between_individual_standard_deviation =
    resting_HFI_slope_SD_60,
  
  between_individual_variance =
    resting_HFI_slope_SD_60^2,
  
  approximate_population_range_68_low =
    population_resting_HFI_slope_60 -
    resting_HFI_slope_SD_60,
  
  approximate_population_range_68_high =
    population_resting_HFI_slope_60 +
    resting_HFI_slope_SD_60,
  
  approximate_population_range_95_low =
    population_resting_HFI_slope_60 -
    1.96 *
    resting_HFI_slope_SD_60,
  
  approximate_population_range_95_high =
    population_resting_HFI_slope_60 +
    1.96 *
    resting_HFI_slope_SD_60
)

print(
  resting_HFI_slope_heterogeneity_60,
  n = Inf
)

#------------------------------------------------------------------------------
# Extract resting HFI slope for every individual ----

individual_resting_slope_template_60 <-
  three_state_random_slope_data_60 %>%
  dplyr::arrange(
    individual_id
  ) %>%
  dplyr::group_by(
    individual_id
  ) %>%
  dplyr::slice(
    1
  ) %>%
  dplyr::ungroup()


individual_resting_slope_data_0_60 <-
  individual_resting_slope_template_60 %>%
  dplyr::mutate(
    hfi_within_z = 0
  )


individual_resting_slope_data_1_60 <-
  individual_resting_slope_template_60 %>%
  dplyr::mutate(
    hfi_within_z = 1
  )


resting_slope_matrix_0_60 <- stats::predict(
  selected_resting_slope_model_60,
  newdata =
    individual_resting_slope_data_0_60,
  type =
    "lpmatrix"
)


resting_slope_matrix_1_60 <- stats::predict(
  selected_resting_slope_model_60,
  newdata =
    individual_resting_slope_data_1_60,
  type =
    "lpmatrix"
)


resting_linear_predictor_columns_60 <- attr(
  resting_slope_matrix_0_60,
  "lpi"
)[[1]]


resting_slope_contrast_matrix_60 <-
  resting_slope_matrix_1_60 -
  resting_slope_matrix_0_60


resting_slope_contrast_matrix_60[
  ,
  -resting_linear_predictor_columns_60
] <- 0


resting_model_coefficients_60 <- stats::coef(
  selected_resting_slope_model_60
)


resting_model_covariance_60 <- stats::vcov(
  selected_resting_slope_model_60,
  unconditional = TRUE
)


individual_resting_HFI_slope_60 <- as.numeric(
  resting_slope_contrast_matrix_60 %*%
    resting_model_coefficients_60
)


individual_resting_HFI_slope_variance_60 <- rowSums(
  (
    resting_slope_contrast_matrix_60 %*%
      resting_model_covariance_60
  ) *
    resting_slope_contrast_matrix_60
)


individual_resting_HFI_slope_SE_60 <- sqrt(
  pmax(
    individual_resting_HFI_slope_variance_60,
    0
  )
)


individual_resting_HFI_intensity_60 <- tibble::tibble(
  individual_id =
    as.character(
      individual_resting_slope_template_60$
        individual_id
    ),
  
  resting_HFI_slope =
    individual_resting_HFI_slope_60,
  
  resting_HFI_slope_magnitude =
    abs(
      resting_HFI_slope
    ),
  
  resting_HFI_slope_SE =
    individual_resting_HFI_slope_SE_60,
  
  confidence_low =
    resting_HFI_slope -
    1.96 *
    resting_HFI_slope_SE,
  
  confidence_high =
    resting_HFI_slope +
    1.96 *
    resting_HFI_slope_SE,
  
  estimated_direction =
    dplyr::case_when(
      resting_HFI_slope < 0 ~
        "negative",
      
      resting_HFI_slope > 0 ~
        "positive",
      
      TRUE ~
        "zero"
    ),
  
  supported_direction =
    dplyr::case_when(
      confidence_high < 0 ~
        "negative",
      
      confidence_low > 0 ~
        "positive",
      
      TRUE ~
        "uncertain"
    )
) %>%
  dplyr::arrange(
    resting_HFI_slope
  )


print(
  individual_resting_HFI_intensity_60,
  n = Inf
)


sd(
  individual_resting_HFI_intensity_60$
    resting_HFI_slope
)
