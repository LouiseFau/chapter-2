#' -----------------------------------------------------------------------------
# Title: Model Validation on GPS-bassed classification of behaviors ----
#' Authors : Louise Faure
#' Date : 24.07.26
#' 
#' Info : this script follow the Covariates&HFI_Selection.R script where covariates
#' and one model representation is selected.  
#' 
# Main steps:
#' (1) prepare the dataset:
#'     (i) calculate elapsed duration in the current behavioural state;
#'     (ii) construct the transition matrix and retain aerial-origin transitions;
#'     (iii) remove incomplete observations and document individual exclusions
#'     (iv) give a weight to individuals 
#'
#' (2) fit two open habitat models Open Habitat models 
#'     (i) standardize covariates
#'     (ii) the first model assume individuals have the same response to HFI 
#'     (iii) the second model assume varying slope per individuals to HFI
#'     (iv) resume results : (a) a table for model comparison (stability of HFI, 
#'     AIC and EDF) and (b) retain the names of the individuals with a negative 
#'     response to HFI.
#'
#' (3) Model controls on autocorrelation within pseudo residuals
#'
#' (4) Confirmation of the results: bootstraping on individuals to confirm the 
#' confidence interval for HFI



# library
library(dplyr)
library(tidyr)
library(tibble)
library(corrplot)
library(ggplot2)
library(mgcv)


# data
GE_60_min_covariates_hfi <- readRDS("/Users/louisefaure/Desktop/dossier sans titre/donnees filtree/GE_60_min_covariates_hfi.rds")

# parameters
state_levels <- c("aerial", "terrestrial")
transition_levels <- c("aerial_to_aerial", "aerial_to_terrestrial", "terrestrial_to_aerial", "terrestrial_to_terrestrial")
gam_selection_method <- "ML"
gam_final_method <- "REML"


#------------------------------------------------------------------------------ STEP 1: prepare dataset ----
#' **Steps:**
#' (i) calculate elapsed time in current behavioural state;
#' (ii) build transition dataset and transition matrix;
#' (iii) retain transitions originating from aerial state;
#' (iv) inspect and remove individuals with insufficient aerial transitions;
#' (v) give a weight to each individuals.


# 1.1 Prepare behavioural states and calculate time spent in current state ----
locations_60 <- GE_60_min_covariates_hfi %>%
  dplyr::mutate(
    timestamp = as.POSIXct(
      timestamp,
      tz = "UTC"
    ),
    
    behavior_state = behavior_binary %>%
      as.character() %>%
      trimws() %>%
      tolower() %>%
      dplyr::recode(
        "aerian" = "aerial",
        "flight" = "aerial"
      ) %>%
      factor(
        levels = state_levels
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
    
    # Identify the beginning of each behavioural bout.
    state_bout_n = cumsum(
      dplyr::row_number() == 1L |
        dplyr::coalesce(
          behavior_state != dplyr::lag(behavior_state),
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
    
    # Time already spent in current behavioural state.
    # This is elapsed duration at the departure point, not total future bout length.
    state_duration_min = as.numeric(
      difftime(
        timestamp,
        dplyr::first(timestamp),
        units = "mins"
      )
    )
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    individual_id = factor(
      individual.local.identifier
    )
  )


# 1.2 Build transitions between consecutive locations ----
transitions_60 <- locations_60 %>%
  dplyr::group_by(
    individual.local.identifier,
    burst_id
  ) %>%
  dplyr::mutate(
    behavior_from = behavior_state,
    
    behavior_to = dplyr::lead(
      behavior_state
    ),
    
    timestamp_next = dplyr::lead(
      timestamp
    ),
    
    dt_min = as.numeric(
      difftime(
        timestamp_next,
        timestamp,
        units = "mins"
      )
    )
  ) %>%
  dplyr::ungroup() %>%
  dplyr::filter(
    !is.na(behavior_from),
    !is.na(behavior_to),
    is.finite(dt_min),
    dt_min > 0
  ) %>%
  dplyr::mutate(
    transition_type = factor(
      paste(
        behavior_from,
        behavior_to,
        sep = "_to_"
      ),
      levels = transition_levels))


# 1.3 Empirical transition matrix (all transitions) ----
transition_count_matrix_60 <- xtabs(
  ~ behavior_from + behavior_to,
  data = transitions_60
)

transition_probability_matrix_60 <- prop.table(
  transition_count_matrix_60,
  margin = 1
)

print(transition_count_matrix_60)
print(round(transition_probability_matrix_60,digits = 4))
#                              Behavior_to
# behavior_from aerial terrestrial       aerial terrestrial
# aerial        4271        6658         0.3908      0.6092
# terrestrial   6520       31274         0.1725      0.8275

# 1.4 Retain only transitions originating from aerial state ----
aerial_transitions_raw_60 <- transitions_60 %>%
  dplyr::filter(
    behavior_from == "aerial") %>%
  dplyr::mutate(
    # 1 = remain aerial
    # 0 = transition to terrestrial
    remain_aerial =
      as.integer(
        behavior_to == "aerial"),
    aerial_duration_min =state_duration_min)

# 1.5 Inspect number of aerial transitions per individual ----
transitions_by_individual_60 <- aerial_transitions_raw_60 %>%
  dplyr::group_by(individual.local.identifier) %>%
  dplyr::summarise(
    n_aerial_transitions = dplyr::n(),
    n_aerial_to_aerial = sum(remain_aerial == 1L),
    n_aerial_to_terrestrial = sum(remain_aerial == 0L),
    n_bursts = dplyr::n_distinct(burst_id),
    .groups = "drop") %>%
  dplyr::arrange(n_aerial_transitions)

print(transitions_by_individual_60,n = Inf)

# 1.6 Remove selected individuals ----
# We remove individuals with less than 30 transitions, we only have 62 individuals.
individuals_to_remove <- c("Langgries21 (eobs 7586)","Almen18 (eobs 5861)","Untersberg21 (eobs 7501)","Schreital22 (eobs 10537)") # for 60 min dataset

# Give an equal weight to each individuals
aerial_transitions_60 <- aerial_transitions_raw_60 %>% dplyr::filter(!individual.local.identifier %in%individuals_to_remove
) %>%
  dplyr::mutate(
    individual_id =
      factor(
        individual.local.identifier
      )) %>%
  tidyr::drop_na(
    individual_id,
    burst_id,
    remain_aerial,
    dplyr::all_of(
      required_numeric_variables_60)) %>%
  dplyr::filter(
    dplyr::if_all(
      dplyr::all_of(
        required_numeric_variables_60
      ),
      is.finite)
  ) %>%
  dplyr::mutate(
    individual_id =
      droplevels(
        individual_id)) %>%
  
  # Number of retained transitions for each individual
  dplyr::group_by(
    individual_id
  ) %>%
  dplyr::mutate(
    n_transitions_individual =
      dplyr::n(),
    
    individual_weight_raw =
      1 /
      n_transitions_individual
  ) %>%
  dplyr::ungroup() %>%
  
  # Normalize weights so their overall mean equals 1
  dplyr::mutate(
    individual_weight =
      individual_weight_raw /
      mean(
        individual_weight_raw))

#------------------------------------------------------------------------------- STEP 2 : model fitting ----
#' **Philosophy**: the reference model is : remain_aerial ~ cos_diel + sin_time 
#' + age_z + duration_z + duration_z2 + s(individual_id, bs = "re") + 
#' hfi_mean_1000m_z + dem_elevation_z + prop_low_vegetation_5cells_z
#' 
#' **Steps**:
#' (i) standardize covariates 
#' (ii) fit two models
#' 