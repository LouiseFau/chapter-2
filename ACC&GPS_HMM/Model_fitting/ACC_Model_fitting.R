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
#' (2) fit a model that assume individuals have the same response to HFI 
#' (3) fit several models that assume varying slope per individuals to HFI 
#'     (either feeding and resting varies, or just feeding, or just resting) and
#'     select the one with the lowest AIC. 
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
expected_dt_min_60 <- 60
dt_tolerance_min_60 <- 60
transition_levels_60 <- tidyr::expand_grid(behavior_from =state_levels_60,behavior_to =state_levels_60) %>%
  dplyr::transmute(transition_type =paste(behavior_from,behavior_to,sep = "_to_")) %>%
  dplyr::pull(transition_type)




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
  dplyr::mutate(timestamp = as.POSIXct( timestamp, tz = "UTC"),
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
        "foraging" = "feeding") %>%
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
    individual.local.identifier) %>%
  dplyr::summarise(
    n_aerial_transitions =
      dplyr::n(),
    n_aerial_to_aerial =
      sum(
        transition_destination == "aerial"),
    n_aerial_to_resting =
      sum(
        transition_destination == "resting"),
    n_aerial_to_feeding =
      sum(
        transition_destination == "feeding"),
    n_bursts =
      dplyr::n_distinct(
        burst_id
      ),
    .groups = "drop" ) %>% dplyr::arrange( n_aerial_transitions)

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
age_mean_60 <-standardization_parameters_60$mean[standardization_parameters_60$variable == age_variable_60]
age_sd_60 <-standardization_parameters_60$standard_deviation[standardization_parameters_60$variable ==age_variable_60]
duration_mean_60 <- standardization_parameters_60$mean[standardization_parameters_60$variable == "aerial_duration_min"]
duration_sd_60 <- standardization_parameters_60$standard_deviation[standardization_parameters_60$variable == "aerial_duration_min"]
hfi_mean_60 <- standardization_parameters_60$mean[standardization_parameters_60$variable == hfi_variable_60]
hfi_sd_60 <- standardization_parameters_60$standard_deviation[standardization_parameters_60$variable == hfi_variable_60]
elevation_mean_60 <- standardization_parameters_60$mean[standardization_parameters_60$variable == elevation_variable_60]
elevation_sd_60 <- standardization_parameters_60$standard_deviation[standardization_parameters_60$variable == elevation_variable_60]
low_vegetation_mean_60 <- standardization_parameters_60$mean[standardization_parameters_60$variable == low_vegetation_variable_60]
low_vegetation_sd_60 <- standardization_parameters_60$standard_deviation[standardization_parameters_60$variable == low_vegetation_variable_60]

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
hfi_q95_q05_difference_z_60 <- three_state_model_data_60 %>%
  dplyr::summarise(
    hfi_q95_q05_difference_z =
      (stats::quantile(
          .data[[hfi_variable_60]],
          probs = 0.95,
          names = FALSE) -
          stats::quantile(
            .data[[hfi_variable_60]],
            probs = 0.05,
            names = FALSE)
      ) /
      hfi_sd_60
  ) %>%
  dplyr::pull(
    hfi_q95_q05_difference_z)


# Extract HFI coefficients and Q05-Q95 log-odds contrasts
hfi_destination_effects_60 <- three_state_model_summary_60$p.table %>%
  as.data.frame() %>%
  tibble::rownames_to_column("term") %>%
  dplyr::filter(
    grepl("^hfi_mean_1000m_z", term)
  ) %>%
  dplyr::mutate(
    destination_contrast = c(
      "resting versus aerial",
      "feeding versus aerial"),
    
    coefficient = Estimate,
    
    coefficient_CI_low =
      coefficient -
      confidence_multiplier_60 * `Std. Error`,
    
    coefficient_CI_high =
      coefficient +
      confidence_multiplier_60 * `Std. Error`,
    
    Q95_Q05_log_odds =
      coefficient *
      hfi_thresholds_60$hfi_q95_q05_difference_z,
    
    Q95_Q05_log_odds_CI_low =
      Q95_Q05_log_odds -
      confidence_multiplier_60 *
      (`Std. Error` *
         abs(hfi_thresholds_60$hfi_q95_q05_difference_z)),
    
    Q95_Q05_log_odds_CI_high =
      Q95_Q05_log_odds +
      confidence_multiplier_60 *
      (`Std. Error` *
         abs(hfi_thresholds_60$hfi_q95_q05_difference_z))
  ) %>%
  dplyr::select(
    destination_contrast,
    coefficient,
    coefficient_CI_low,
    coefficient_CI_high,
    Q95_Q05_log_odds,
    Q95_Q05_log_odds_CI_low,
    Q95_Q05_log_odds_CI_high
  )

hfi_destination_effects_60 <- three_state_model_summary_60$p.table %>%
    as.data.frame() %>%
    tibble::rownames_to_column("term") %>%
    dplyr::filter(
      grepl("^hfi_mean_1000m_z", term)
    ) %>%
    dplyr::mutate(
      destination_contrast = c(
        "resting versus aerial",
        "feeding versus aerial"
      ),
      
      coefficient = Estimate,
      
      coefficient_CI_low =
        coefficient -
        confidence_multiplier_60 * `Std. Error`,
      
      coefficient_CI_high =
        coefficient +
        confidence_multiplier_60 * `Std. Error`,
      
      Q95_Q05_log_odds =
        coefficient *
        hfi_q95_q05_difference_z_60,
      
      Q95_Q05_log_odds_CI_low =
        Q95_Q05_log_odds -
        confidence_multiplier_60 *
        (`Std. Error` *
           abs(hfi_q95_q05_difference_z_60)),
      
      Q95_Q05_log_odds_CI_high =
        Q95_Q05_log_odds +
        confidence_multiplier_60 *
        (`Std. Error` *
           abs(hfi_q95_q05_difference_z_60))
    ) %>%
    dplyr::select(
      destination_contrast,
      coefficient,
      coefficient_CI_low,
      coefficient_CI_high,
      Q95_Q05_log_odds,
      Q95_Q05_log_odds_CI_low,
      Q95_Q05_log_odds_CI_high)
  
print(hfi_destination_effects_60)
#                       coefficient            CI coefficent            Q95_Q05_log_odds                CI log odds
# resting versus aerial -0.02134393        [-0.08095036;0.03826251]       -0.06439953              [-0.2442458;0.1154468] 
# feeding versus aerial -0.07422699        [-0.20805357;0.05959960]       -0.22395987              [-0.6277454;0.1798257] 

# 2.5 Control HFI support by transition destination ----
hfi_q10_absolute_60 <- 0.10
hfi_q80_absolute_60 <- 0.80

hfi_support_by_destination_60 <- three_state_model_data_60 %>%
  dplyr::filter(
    !is.na(transition_destination),
    !is.na(hfi_mean_1000m)
  ) %>%
  dplyr::mutate(
    # Ensures that bursts are unique between individuals
    burst_uid = interaction(
      individual_id,
      burst_id,
      drop = TRUE,
      lex.order = TRUE
    )
  ) %>%
  dplyr::group_by(
    transition_destination
  ) %>%
  dplyr::summarise(
    # Empirical HFI range for this transition
    hfi_min =
      min(
        hfi_mean_1000m,
        na.rm = TRUE
      ),
    
    hfi_max =
      max(
        hfi_mean_1000m,
        na.rm = TRUE
      ),
    
    # Number of individuals represented in each extreme class
    n_individuals_HFI_q10 =
      dplyr::n_distinct(
        individual_id[
          hfi_mean_1000m <= hfi_q10_absolute_60
        ],
        na.rm = TRUE
      ),
    
    n_individuals_HFI_q80 =
      dplyr::n_distinct(
        individual_id[
          hfi_mean_1000m > hfi_q80_absolute_60
        ],
        na.rm = TRUE
      ),
    
    # Number of bursts represented in each extreme class
    n_bursts_HFI_q10 =
      dplyr::n_distinct(
        burst_uid[
          hfi_mean_1000m <= hfi_q10_absolute_60
        ],
        na.rm = TRUE
      ),
    
    n_bursts_HFI_q80 =
      dplyr::n_distinct(
        burst_uid[
          hfi_mean_1000m > hfi_q80_absolute_60
        ],
        na.rm = TRUE
      ),
    
    .groups = "drop"
  ) %>%
  dplyr::arrange(
    transition_destination
  )

print(hfi_support_by_destination_60,n = Inf)
# transition_destination hfi_min hfi_max n_individuals_HFI_q10 n_individuals_HFI_q80 n_bursts_HFI_q10 n_bursts_HFI_q80
# aerial                       0   0.997                    57                     5             1464                5
# resting                      0   0.844                    58                     3             2575                3
# feeding                      0   0.881                    53                     1              350                1



# ------------------------------------------------------------------------------- STEP 3 : model fitting with common and individual HFI slopes ----
#' **Philosophy:**
#' Test whether within-individual HFI responses vary among individuals
#' for resting versus aerial and feeding versus aerial.
#'
#' **Steps:**
#' (i) fit two models:
#'     a) common HFI slope model;
#'     b) individual HFI slope model for resting and feeding.
#' (ii) compare AIC and random-slope variance.
#' (iii) extract individual HFI coefficients and rank individuals.


# 3.1 Separate within- and between-individual HFI variation ----
three_state_random_slope_data_60 <-
  three_state_model_data_60 %>%
  dplyr::mutate(
    individual_id =
      droplevels(
        factor(individual_id)
      )
  ) %>%
  dplyr::group_by(
    individual_id
  ) %>%
  dplyr::mutate(
    hfi_between_z =
      mean(hfi_mean_1000m_z),
    
    hfi_within_z =
      hfi_mean_1000m_z -
      hfi_between_z
  ) %>%
  dplyr::ungroup()


# 3.2 Define common and individual slope models ----

common_HFI_formula_60 <- list(
  
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
)


individual_HFI_formula_60 <- list(
  
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


three_state_slope_formulas_60 <- list(
  common_HFI_slopes =
    common_HFI_formula_60,
  
  individual_HFI_slopes =
    individual_HFI_formula_60
)


# 3.3 Fit both models using ML ----
three_state_slope_models_ml_60 <-
  lapply(
    three_state_slope_formulas_60,
    function(formula_set) {
      
      mgcv::gam(
        formula =
          formula_set,
        
        family =
          mgcv::multinom(
            K = 2
          ),
        
        data =
          three_state_random_slope_data_60,
        
        method =
          "ML"
      )
    }
  )


# 3.4 Compare models using AIC ----
three_state_slope_model_comparison_60 <-
  tibble::tibble(
    model = names(three_state_slope_models_ml_60),
    AIC = vapply(
        three_state_slope_models_ml_60,
        stats::AIC,
        numeric(1))) %>%
  dplyr::mutate(
    delta_AIC =
      AIC -
      min(AIC))


print(as.data.frame(three_state_slope_model_comparison_60),row.names = FALSE)
# model                   AIC      delta_AIC
# common_HFI_slopes     8958.297  2.129646
# individual_HFI_slopes 8956.167  0.000000


# 3.5 Select best model and refit using REML ----
selected_model_name_60 <-
  three_state_slope_model_comparison_60 %>%
  dplyr::slice_min(
    AIC,
    n = 1
  ) %>%
  dplyr::pull(
    model
  )


selected_three_state_slope_model_60 <-
  mgcv::gam(formula =
      three_state_slope_formulas_60[[selected_model_name_60]],
    family =
      mgcv::multinom(
        K = 2
      ),
    data =
      three_state_random_slope_data_60,
    method =
      "REML")


# 3.6 Extract individual HFI coefficients ----
individual_slope_template_60 <-
  three_state_random_slope_data_60 %>%
  dplyr::arrange(
    individual_id
  ) %>%
  dplyr::distinct(
    individual_id,
    .keep_all = TRUE
  )


slope_matrix_0_60 <-
  stats::predict(
    selected_three_state_slope_model_60,
    newdata =
      individual_slope_template_60 %>%
      dplyr::mutate(
        hfi_within_z = 0
      ),
    type = "lpmatrix"
  )


slope_matrix_1_60 <-
  stats::predict(
    selected_three_state_slope_model_60,
    newdata =
      individual_slope_template_60 %>%
      dplyr::mutate(
        hfi_within_z = 1
      ),
    type = "lpmatrix"
  )


slope_difference_matrix_60 <-
  slope_matrix_1_60 -
  slope_matrix_0_60


linear_predictor_indices_60 <-
  attr(
    slope_matrix_0_60,
    "lpi"
  )


model_coefficients_60 <-
  stats::coef(
    selected_three_state_slope_model_60
  )


individual_hfi_slopes_60 <-
  tibble::tibble(
    
    individual_id =
      as.character(
        individual_slope_template_60$individual_id
      ),
    
    resting_HFI_coefficient =
      as.numeric(
        slope_difference_matrix_60[
          ,
          linear_predictor_indices_60[[1]],
          drop = FALSE
        ] %*%
          model_coefficients_60[
            linear_predictor_indices_60[[1]]
          ]
      ),
    
    feeding_HFI_coefficient =
      as.numeric(
        slope_difference_matrix_60[
          ,
          linear_predictor_indices_60[[2]],
          drop = FALSE
        ] %*%
          model_coefficients_60[
            linear_predictor_indices_60[[2]]
          ]
      )
  )


# 3.8 Rank individuals by HFI response combinations ----

individual_hfi_grouped_table_60 <-
  individual_hfi_slopes_60 %>%
  
  dplyr::mutate(
    
    response_group =
      dplyr::case_when(
        
        feeding_HFI_coefficient > 0 &
          resting_HFI_coefficient > 0 ~
          
          "Positive feeding and positive resting",
        
        feeding_HFI_coefficient > 0 &
          resting_HFI_coefficient <= 0 ~
          
          "Positive feeding only",
        
        feeding_HFI_coefficient <= 0 &
          resting_HFI_coefficient > 0 ~
          
          "Positive resting only",
        
        TRUE ~
          
          "Negative feeding and negative resting"
      ),
    
    response_group_order =
      dplyr::case_when(
        
        response_group ==
          "Positive feeding and positive resting" ~ 1L,
        
        response_group ==
          "Positive feeding only" ~ 2L,
        
        response_group ==
          "Positive resting only" ~ 3L,
        
        TRUE ~ 4L
      ),
    
    mean_HFI_coefficient =
      (
        feeding_HFI_coefficient +
          resting_HFI_coefficient
      ) / 2
  ) %>%
  
  dplyr::arrange(
    
    response_group_order,
    
    dplyr::desc(
      dplyr::if_else(
        response_group_order == 4L,
        mean_HFI_coefficient,
        NA_real_
      )
    ),
    
    individual_id
  ) %>%
  
  dplyr::select(
    individual_id,
    feeding_HFI_coefficient,
    resting_HFI_coefficient
  )


base::print(
  as.data.frame(
    individual_hfi_grouped_table_60
  ),
  row.names = FALSE
)


# individual_id                 feeding_HFI_coefficient resting_HFI_coefficient
# Ettenberg22 (eobs 10539)            0.0269657106            0.0314341768
# Fahrntal19 (eobs 7014)            0.0111615815            0.0035551866
# Mals2_23 (eobs 11920)            0.0001213819           -0.0016769867
# Punteglias20 (eobs 6483)            0.0541185759           -0.0262912029
# Reschen21 (eobs 7503)            0.0875493847           -0.0819583185
# Almen19 (eobs 7001)           -0.0104048982            0.0560266339
# Dischma1 19 (eobs 7006)           -0.0408479071            0.0134970845
# Flüela1 21 (eobs 6995)           -0.0850378469            0.0830143770
# Grabernock21 (eobs 7506)           -0.1945659582            0.0005980758
# Grabernock23 (eobs 11916)           -0.0841109723            0.0205776725
# Mals1_23 (eobs 11919)           -0.1941829706            0.0731668082
# Siat20 (eobs 7037)           -0.1039544268            0.0432002325
# Trimmis20 (eobs 7041)           -0.0296064147            0.0237016391
# Vernuga22 (eobs 10538)           -0.1091579224            0.0050119932
# Flüela19 (eobs 7007)           -0.0059689917           -0.0072103624
# Dischma2 19 (eobs 7009)           -0.0188923587           -0.0015575114
# Sampuoir2 19 (eobs 6462)           -0.0011186934           -0.0615493464
# Sampuoir1 19 (eobs 5943)           -0.0091042380           -0.0536618812
# Matsch19 (eobs 7035)           -0.0385038080           -0.0268853753
# Tabland22 (eobs 10534)           -0.0327932871           -0.0348893356
# Tuors1 19 (eobs 7010)           -0.0016529109           -0.0753120796
# Umbrail18 (eobs 5859)           -0.0672733452           -0.0149207593
# Schlappin22 (eobs 5944)           -0.0242459937           -0.0741773876
# Laas2_23 (eobs 11918)           -0.0624631059           -0.0391458250
# ValGrande19 (eobs 7033)           -0.0562137559           -0.0496303991
# Trenzeira19 (eobs 5858)           -0.0910444291           -0.0211689345
# Vrata20 (eobs 7551)           -0.0757188286           -0.0377498556
# Krn20 (eobs 7549)           -0.0910499956           -0.0239480797
# ValSozzine21 (eobs 7500)           -0.1078958794           -0.0107895760
# Johnsbach2_23 (eobs 11913)           -0.0944378893           -0.0284618993
# Burgum21 (eobs 7502)           -0.0550975654           -0.0695247881
# Lischana22 (eobs 5941)           -0.0913107806           -0.0361038128
# Flüela20 (eobs 7040)           -0.1291184950           -0.0036212429
# Tuors2 19 (eobs 7011)           -0.1034797746           -0.0293772581
# Schlappin1 18 (eobs 5858)           -0.1289949922           -0.0072022382
# Lassingbach24 (eobs 11915)           -0.1096252758           -0.0275563762
# Sinestra1 19 (eobs 7003)           -0.1373278721           -0.0061369272
# Kastelbell19 (eobs 7034)           -0.1166658414           -0.0271954364
# Grosio 19 (eobs 7000)           -0.0949653784           -0.0531855882
# Laas1_23 (eobs 11917)           -0.1519108717           -0.0080830406
# Flüela2 21 (eobs 7043)           -0.1200103543           -0.0421684529
# Gaming24 (eobs 11914)           -0.1588517304           -0.0050114683
# Seta19 (eobs 5796)           -0.0552775763           -0.1093898530
# Güstizia18 (eobs 5942)           -0.1071469504           -0.0658690709
# Sinestra2 19 (eobs 7005)           -0.1708937152           -0.0021438582
# Reschen20 (eobs 7556)           -0.1287326616           -0.0471091996
# Nalps18 (eobs 5860)           -0.1658944641           -0.0123143796
# Cornasc20 (eobs 7039)           -0.1278909630           -0.0557342075
# Johnsbach1_23 (eobs 11645)           -0.1429239973           -0.0529544185
# Adamello20 (eobs 7548)           -0.1852416267           -0.0126967549
# Art San Romerio18 (eobs 5941)           -0.1776761005           -0.0377304681
# Tasna18 (eobs 5940)           -0.1720759557           -0.0457719126
# Torta19 (eobs 7002)           -0.1781254002           -0.0554974110
# Sils20 (eobs 7038)           -0.1347077444           -0.1037554982
# Avers20 (eobs 7101)           -0.0858493654           -0.1532864892
# Nalps19 (eobs 5861)           -0.1733887418           -0.1012669749
# Stürfis20 (eobs 7049)           -0.2773412001           -0.0035957838
# Valdidentro_Braulio20 (eobs 7581)           -0.2615049800           -0.0268970116