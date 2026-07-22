#' -----------------------------------------------------------------------------
# Title: Covariate and HFI variables selection ----
#' Authors : Louise Faure
#' Date : 20.07.26
#' 
#' Info : this script follow the Extract_covariates.R script where covariates are
#' extracted below each location and within two buffers. 
#' 
# Main steps: ----
#' (1) prepare the dataset:
#'     (i) calculate elapsed duration in the current behavioural state;
#'     (ii) construct the transition matrix and retain aerial-origin transitions;
#'     (iii) remove incomplete observations and document individual exclusions.
#'
#' (2) define the behavioural backbone model:
#'     (i) retain the diel cosinor and individual random intercept;
#'     (ii) compare biologically plausible linear, quadratic or smooth forms for
#'          age and elapsed aerial duration;
#'     (iii) retain the most parsimonious supported structure.
#'
#' (3) identify potential confounders:
#'     (i) standardize covariates
#'     (ii) inspect correlations using Pearson correlation coefficient;
#'     (ii) reduce strongly correlated topographic or habitat variables using
#'           biological justification rather than automatic selection.
#'
#' (4) identify several HFI metrics:
#'     (i) compile pearson correlation coefficient btw hfi metrics and retained covariates
#'     (ii) compile vif amongst group of covariates for the models defined in step 5
#'     
#' (5) fit candidates HFI models
#' (6) diagnose pseudo residuals
#' (7) fit model with variable HFI slope per individuals   


# library ----
library(dplyr)
library(tidyr)
library(tibble)
library(corrplot)
library(ggplot2)
library(mgcv)


# data ----
GE_60_min_covariates_hfi <- readRDS("/Users/louisefaure/Desktop/dossier sans titre/donnees filtree/GE_60_min_covariates_hfi.rds")


# parameters ----
state_levels <- c("aerial", "terrestrial")
transition_levels <- c("aerial_to_aerial", "aerial_to_terrestrial", "terrestrial_to_aerial", "terrestrial_to_terrestrial")
hfi_candidates_60 <- c("hfi_point", "hfi_mean_500m", "hfi_q75_500m", "hfi_mean_1000m", "hfi_q75_1000m", "hfi_q90_500m","hfi_q90_1000m")
topographic_covariates_60 <- c("dem_elevation","ruggedness_100m","slope_100m","distance_to_ridgeline_100m")
habitat_covariates_60 <- c("prop_forest_5cells","prop_low_vegetation_5cells","prop_rocky_terrain_5cells")
environmental_covariates_60 <- c(topographic_covariates_60,habitat_covariates_60)
backbone_variables_60 <- c("age_since_emig_weeks","aerial_duration_min","cos_diel","sin_time")
required_numeric_variables_60 <- c(backbone_variables_60,environmental_covariates_60,hfi_candidates_60)
delta_aic_threshold <- 4

gam_selection_method <- "ML"
gam_final_method <- "REML"
age_k <- 5
duration_k <- 6


#------------------------------------------------------------------------------
# STEP 1: prepare dataset ----
#' **Steps:**
#' (i) calculate elapsed time in current behavioural state;
#' (ii) build transition dataset and transition matrix;
#' (iii) retain transitions originating from aerial state;
#' (iv) inspect and remove individuals with insufficient aerial transitions;
#' (v) remove incomplete observations.


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
      levels = transition_levels
    )
  )


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
print(
  round(
    transition_probability_matrix_60,
    digits = 4
  )
)

# Rows represent departure states and columns represent arrival states.
# Each row of transition_probability_matrix_60 should sum to one.


# 1.4 Retain only transitions originating from aerial state ----
aerial_transitions_raw_60 <- transitions_60 %>%
  dplyr::filter(
    behavior_from == "aerial"
  ) %>%
  dplyr::mutate(
    
    # Response variable:
    # 1 = remain aerial
    # 0 = transition to terrestrial
    remain_aerial = as.integer(
      behavior_to == "aerial"
    ),
    
    aerial_duration_min = state_duration_min
  )


# 1.5 Inspect number of aerial transitions per individual ----
transitions_by_individual_60 <- aerial_transitions_raw_60 %>%
  dplyr::group_by(
    individual.local.identifier
  ) %>%
  dplyr::summarise(
    n_aerial_transitions = dplyr::n(),
    
    n_aerial_to_aerial = sum(
      remain_aerial == 1L
    ),
    
    n_aerial_to_terrestrial = sum(
      remain_aerial == 0L
    ),
    
    n_bursts = dplyr::n_distinct(
      burst_id
    ),
    
    .groups = "drop"
  ) %>%
  dplyr::arrange(
    n_aerial_transitions
  )

print(transitions_by_individual_60,n = Inf)


# 1.6 Remove selected individuals ----
# We remove individuals with less than 30 transitions, we only have 62 individuals.
individuals_to_remove <- c("Langgries21 (eobs 7586)","Almen18 (eobs 5861)","Untersberg21 (eobs 7501)","Schreital22 (eobs 10537)") # for 60 min dataset

aerial_transitions_60 <- aerial_transitions_raw_60 %>%
  dplyr::filter(
    !individual.local.identifier %in%
      individuals_to_remove
  ) %>%
  dplyr::mutate(
    individual_id = factor(
      individual.local.identifier
    )
  ) %>%
  tidyr::drop_na(
    individual_id,
    burst_id,
    remain_aerial,
    dplyr::all_of(
      required_numeric_variables_60
    )
  ) %>%
  dplyr::filter(
    dplyr::if_all(
      dplyr::all_of(
        required_numeric_variables_60
      ),
      is.finite
    )
  ) %>%
  dplyr::mutate(
    individual_id = droplevels(
      individual_id
    )
  )


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
    
    n_remain_aerial =
      sum(
        remain_aerial == 1L
      ),
    
    n_transition_terrestrial =
      sum(
        remain_aerial == 0L
      ),
    
    proportion_remain_aerial =
      mean(
        remain_aerial
      )
  )

print(final_dataset_summary_60)






#-------------------------------------------------------------------------------STEP 2: define the behavioural backbone model ----
#' **Steps:**
#' (i) retain the diel cosinor and individual random intercept;
#' (ii) compare biologically plausible linear, quadratic or smooth forms for
#'      age and elapsed aerial duration;
#' (iii) retain the most parsimonious supported structure, preferring a simpler
#'       form when its delta AIC is below 4.


# 2.1 Prepare the common backbone dataset ----
# All candidate models are fitted to exactly the same observations.
backbone_data_60 <- aerial_transitions_60 %>%
  dplyr::filter(
    !is.na(individual_id),
    dplyr::if_all(
      dplyr::all_of(backbone_variables_60),
      is.finite
    )
  ) %>%
  dplyr::mutate(
    individual_id = droplevels(
      factor(individual_id)
    ),
    
    aerial_duration_h =
      aerial_duration_min / 60
  )


# 2.2 Centre and standardize age and elapsed aerial duration ----
age_center_60 <- mean(backbone_data_60$age_since_emig_weeks)
age_scale_60 <- stats::sd(backbone_data_60$age_since_emig_weeks)
duration_center_60 <- mean(backbone_data_60$aerial_duration_h)
duration_scale_60 <- stats::sd(backbone_data_60$aerial_duration_h)

backbone_data_60 <- backbone_data_60 %>%
  dplyr::mutate(
    age_z = (
      age_since_emig_weeks -
        age_center_60
    ) / age_scale_60,
    
    duration_z = (
      aerial_duration_h -
        duration_center_60
    ) / duration_scale_60,
    
    age_z2 = age_z^2,
    duration_z2 = duration_z^2
  )

backbone_scaling_60 <- tibble::tibble(
  variable = c(
    "age_since_emig_weeks",
    "aerial_duration_h"
  ),
  
  center = c(
    age_center_60,
    duration_center_60
  ),
  
  scale = c(
    age_scale_60,
    duration_scale_60
  )
)

print(backbone_scaling_60)


# 2.3 Define candidate functional forms ----
# Quadratic forms always contain both the linear and squared terms.
# Shrinkage cubic splines ("cs") can be penalized towards a simple or
# effectively absent relationship.
age_terms_60 <- list(
  linear =
    "age_z",
  
  quadratic = c(
    "age_z",
    "age_z2"
  ),
  
  smooth = paste0(
    "s(age_z, bs = 'cs', k = ",
    age_k,
    ")"
  )
)

duration_terms_60 <- list(
  linear =
    "duration_z",
  
  quadratic = c(
    "duration_z",
    "duration_z2"
  ),
  
  smooth = paste0(
    "s(duration_z, bs = 'cs', k = ",
    duration_k,
    ")"
  )
)


# 2.4 Define the candidate backbone model set ----
backbone_model_metadata_60 <- tidyr::expand_grid(
  age_form = c(
    "linear",
    "quadratic",
    "smooth"
  ),
  
  duration_form = c(
    "linear",
    "quadratic",
    "smooth"
  )
) %>%
  dplyr::mutate(
    model = paste(
      "age",
      age_form,
      "duration",
      duration_form,
      sep = "_"
    ),
    
    # This ranking implements the a priori preference for simpler structures.
    age_complexity = dplyr::recode(
      age_form,
      linear = 1L,
      quadratic = 2L,
      smooth = 3L
    ),
    
    duration_complexity = dplyr::recode(
      duration_form,
      linear = 1L,
      quadratic = 2L,
      smooth = 3L
    ),
    
    complexity_score =
      age_complexity +
      duration_complexity
  )


# 2.5 Build candidate model formulas ----
# cos_diel and sin_time are retained together in every model.
# s(individual_id, bs = "re") represents the individual random intercept.
backbone_formulas_60 <- stats::setNames(
  lapply(
    seq_len(
      nrow(backbone_model_metadata_60)
    ),
    function(i) {
      
      model_information <-
        backbone_model_metadata_60[i, ]
      
      model_terms <- c(
        "cos_diel",
        "sin_time",
        
        age_terms_60[[
          model_information$age_form
        ]],
        
        duration_terms_60[[
          model_information$duration_form
        ]],
        
        "s(individual_id, bs = 're')"
      )
      
      stats::as.formula(
        paste(
          "remain_aerial ~",
          paste(
            model_terms,
            collapse = " + "
          )
        )
      )
    }
  ),
  backbone_model_metadata_60$model
)


# 2.6 Fit candidate backbone models ----
backbone_models_60 <- stats::setNames(
  lapply(
    names(backbone_formulas_60),
    function(model_name) {
      
      mgcv::gam(
        formula =
          backbone_formulas_60[[
            model_name
          ]],
        
        data =
          backbone_data_60,
        
        family =
          stats::binomial(
            link = "logit"
          ),
        
        method =
          gam_selection_method
      )
    }
  ),
  names(backbone_formulas_60)
)


# 2.7 Compare candidate backbone models ----
backbone_model_comparison_60 <- dplyr::bind_rows(
  lapply(
    names(backbone_models_60),
    function(model_name) {
      
      fitted_model <-
        backbone_models_60[[
          model_name
        ]]
      
      model_information <-
        backbone_model_metadata_60 %>%
        dplyr::filter(
          model ==
            model_name
        )
      
      model_summary <-
        summary(fitted_model)
      
      loglik_model <-
        stats::logLik(
          fitted_model
        )
      
      tibble::tibble(
        model =
          model_name,
        
        age_form =
          model_information$age_form,
        
        duration_form =
          model_information$duration_form,
        
        complexity_score =
          model_information$complexity_score,
        
        n_observations =
          nrow(
            stats::model.frame(
              fitted_model
            )
          ),
        
        df =
          attr(
            loglik_model,
            "df"
          ),
        
        log_likelihood =
          as.numeric(
            loglik_model
          ),
        
        AIC =
          stats::AIC(
            fitted_model
          ),
        
        deviance_explained =
          model_summary$dev.expl,
        
        converged =
          fitted_model$converged,
        
        convergence_message =
          paste(
            fitted_model$outer.info$conv,
            collapse = "; "
          )
      )
    }
  )
) %>%
  dplyr::arrange(AIC) %>%
  dplyr::mutate(
    delta_AIC =
      AIC - min(AIC),
    
    akaike_weight =
      exp(-0.5 * delta_AIC) /
      sum(
        exp(-0.5 * delta_AIC)
      ),
    
    similar_support =
      delta_AIC <=
      delta_aic_threshold
  )

print(backbone_model_comparison_60,n = Inf)

# Inspect:
# - n_observations: must be identical among all models;
# - AIC: lower values indicate stronger relative support;
# - delta_AIC <= 4: models are retained as reasonably supported;
# - df: effective model complexity, including smooth complexity;
# - deviance_explained: descriptive fit of the complete model;
# - converged: expected TRUE;
# - convergence_message: expected to indicate full convergence.


# 2.8 Retain the most parsimonious supported structure ----
# Models with delta_AIC <= 4 are considered supported.
# Within this set, the simplest predefined structure is selected first.
# AIC is used to order models with the same complexity score.
supported_backbone_models_60 <-
  backbone_model_comparison_60 %>%
  dplyr::filter(
    similar_support
  ) %>%
  dplyr::arrange(
    complexity_score,
    df,
    AIC
  )

print(supported_backbone_models_60,n = Inf)
# we take the first model even though the AIC is higher

provisional_backbone_name_60 <-
  supported_backbone_models_60 %>%
  dplyr::slice(1) %>%
  dplyr::pull(model)

provisional_backbone_information_60 <-
  backbone_model_metadata_60 %>%
  dplyr::filter(
    model ==
      provisional_backbone_name_60)

selected_age_form_60 <-
  provisional_backbone_information_60$
  age_form

selected_duration_form_60 <-
  provisional_backbone_information_60$
  duration_form

provisional_backbone_model_ml_60 <-
  backbone_models_60[[
    provisional_backbone_name_60
  ]]

print(provisional_backbone_name_60)
print(summary(provisional_backbone_model_ml_60))


# ------------------------------------------------------------------------------STEP 3 : identify potential cofounders ----
#'**Steps**:
#' (i) define topographic and habitat covariates from biological hypotheses;
#' (ii) inspect pairwise Pearson correlations;
#' (iii) reduce strongly correlated variables using biological justification
#'       rather than automatic selection.


# Calculate Pearson correlations among environmental covariates
correlation_data_60 <- aerial_transitions_60 %>%
  dplyr::select(
    dplyr::all_of(
      environmental_covariates_60
    )
  )

correlation_matrix_60 <- stats::cor(correlation_data_60,method = "pearson",use = "complete.obs")

#------------------------------------------------------------------------------- Visualisation n°1 : Pearson correlation matrix
corrplot::corrplot(
  correlation_matrix_60,
  
  method = "color",
  type = "lower",
  order = "original",
  
  col = corrplot::COL2(
    "RdBu",
    200
  ),
  
  col.lim = c(
    -1,
    1
  ),
  
  diag = FALSE,
  
  addCoef.col = "black",
  number.digits = 2,
  number.cex = 0.6,
  
  tl.col = "black",
  tl.cex = 0.6,
  tl.srt = 45,
  
  cl.cex = 0.8,
  addgrid.col = "white",
  
  title =
    "Pearson correlations among candidate environmental confounders",
  
  mar = c(
    0,
    0,
    3,
    0
  )
)
#------------------------------------------------------------------------------- (end) Visualisation n°1 : Pearson correlation matrix

#' we remove slope due to its high correlation with ruggedness and lower 
#' explanation of the topographic complexity 
#' we remove as well proportion of rocky terrain for its association with 
#' elevation and limited biological interpret ability











# ------------------------------------------------------------------------------ Step 4 : identify several HFI metrics ----
#' **Steps:**
#' (i) compile support per hfi class;
#' (ii) compile VIF for retained HFI metrics and their covariate groups.
#' 
#' Pearson correlations between HFI metrics and retained environmental
#' covariates were all below |r| = 0.70. HFI metric discrimination therefore
#' relies on biological interpretation, empirical support and VIF.
#'                dem_elevation ruggedness_100m distance_to_ridgeline_100m prop_forest_5cells prop_low_vegetation_5cells
#'                hfi_point              -0.58           -0.23                       0.29               0.30                      -0.05
#'                hfi_mean_500m          -0.63           -0.22                       0.28               0.34                      -0.06
#'                hfi_q75_500m           -0.64           -0.20                       0.26               0.37                      -0.07
#'                hfi_mean_1000m         -0.68           -0.18                       0.24               0.41                      -0.10
#'                hfi_q75_1000m          -0.68           -0.16                       0.23               0.43                      -0.10
#'                hfi_q90_500m           -0.65           -0.17                       0.25               0.40                      -0.08
#'                hfi_q90_1000m          -0.68           -0.13                       0.20               0.46                      -0.11



# 4.1 Describe the distribution and empirical support of HFI metrics ----
# 4.1.1 Define absolute HFI classes and reshape candidate metrics ----
hfi_classes_60 <- c(0,0.10,0.20,0.60,0.80,Inf)
hfi_labels_60 <- c("0-0.10","0.10-0.20","0.20-0.60","0.60-0.80",">0.80")

hfi_support_long_60 <- aerial_transitions_60 %>%
  dplyr::select(
    individual_id,
    burst_id,
    remain_aerial,
    dplyr::all_of(
      hfi_candidates_60
    )
  ) %>%
  tidyr::pivot_longer(
    cols =
      dplyr::all_of(
        hfi_candidates_60
      ),
    
    names_to =
      "hfi_variable",
    
    values_to =
      "hfi_value"
  ) %>%
  dplyr::mutate(
    hfi_variable = factor(
      hfi_variable,
      levels =
        hfi_candidates_60
    )
  )


# 4.1.2 Compile the empirical range of each HFI metric ----
hfi_thresholds_60 <- hfi_support_long_60 %>%
  dplyr::group_by(
    hfi_variable
  ) %>%
  dplyr::summarise(
    hfi_min =
      min(hfi_value),
    
    hfi_q05 =
      as.numeric(
        stats::quantile(
          hfi_value,
          probs = 0.05
        )
      ),
    
    hfi_median =
      stats::median(
        hfi_value
      ),
    
    hfi_q95 =
      as.numeric(
        stats::quantile(
          hfi_value,
          probs = 0.95
        )
      ),
    
    hfi_max =
      max(hfi_value),
    
    .groups =
      "drop"
  )

print(hfi_thresholds_60,n = Inf)


# 4.1.3 Count observations, individuals and bursts across HFI classes ----
hfi_class_support_60 <- hfi_support_long_60 %>%
  dplyr::mutate(
    hfi_class = cut(
      hfi_value,
      breaks =
        hfi_classes_60,
      labels =
        hfi_labels_60,
      include.lowest =
        TRUE,
      right =
        TRUE
    )
  ) %>%
  dplyr::group_by(
    hfi_variable,
    hfi_class
  ) %>%
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
    
    n_remain_aerial =
      sum(
        remain_aerial == 1L
      ),
    
    n_transition_terrestrial =
      sum(
        remain_aerial == 0L
      ),
    
    proportion_remain_aerial =
      mean(
        remain_aerial
      ),
    
    .groups =
      "drop"
  ) %>%
  dplyr::arrange(
    hfi_variable,
    hfi_class
  )

print(hfi_class_support_60,n = Inf)


# 4.1.4 Prepare individual and transition counts for one plot ----
hfi_class_plot_data_60 <- hfi_class_support_60 %>%
  dplyr::select(
    hfi_variable,
    hfi_class,
    n_individuals,
    n_observations
  ) %>%
  tidyr::pivot_longer(
    cols = c(
      n_individuals,
      n_observations
    ),
    
    names_to =
      "support_measure",
    
    values_to =
      "count"
  ) %>%
  dplyr::mutate(
    support_measure = dplyr::recode(
      support_measure,
      n_individuals =
        "Number of individuals",
      n_observations =
        "Number of aerial-origin transitions"
    )
  )

#------------------------------------------------------------------------------- Visualisation n°2 : empirical support
hfi_support_comparison_plot_60 <- hfi_class_plot_data_60 %>%
  ggplot2::ggplot(
    ggplot2::aes(
      x =
        hfi_class,
      y =
        count,
      group =
        hfi_variable,
      color =
        hfi_variable
    )
  ) +
  ggplot2::geom_line(
    linewidth = 0.9
  ) +
  ggplot2::geom_point(
    size = 2.5
  ) +
  ggplot2::facet_wrap(
    ~ support_measure,
    ncol = 1,
    scales = "free_y"
  ) +
  ggplot2::labs(
    x =
      "Absolute HFI class",
    y =
      "Count",
    color =
      "HFI metric",
    title =
      "Empirical support across HFI classes",
    subtitle =
      "Comparison of individual and transition support among candidate HFI metrics"
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    axis.text.x =
      ggplot2::element_text(
        angle = 45,
        hjust = 1
      ),
    
    legend.position =
      "right"
  )

print(
  hfi_support_comparison_plot_60
)
#------------------------------------------------------------------------------- (end) Visualisation n°2 : empirical support
selected_hfi_candidates_60 <- c("hfi_mean_1000m","hfi_q75_1000m","hfi_q75_500m")


# 4.2 Compile VIF for HFI candidate models ----
#' **Steps:**
#' (i) calculate VIF for each retained HFI metric within each environmental model;
#' (ii) compile a detailed table and a summary table;
#' (iii) plot both tables as images.

# 4.2.1 Define the candidate environmental model structures ----
vif_model_structures_60 <- list(
  
  altitudinal_gradient = c(
    "dem_elevation"
  ),
  
  topographic_complexity = c(
    "ruggedness_100m",
    "distance_to_ridgeline_100m", 
    "dem_elevation"
  ),
  
  habitat_refuge = c(
    "prop_forest_5cells",
    "ruggedness_100m",
    "dem_elevation"
  ),
  
  open_habitat = c(
    "prop_low_vegetation_5cells",
    "dem_elevation"
  )
)


# 4.2.2 Define the VIF calculation function ----
# For each predictor:
# VIF = 1 / (1 - R²),
# where R² is obtained by regressing the predictor on the other predictors
# included in the same candidate model.
calculate_vif_60 <- function(data, predictors) {
  
  dplyr::bind_rows(
    lapply(
      predictors,
      function(variable_name) {
        
        other_predictors <- setdiff(
          predictors,
          variable_name
        )
        
        auxiliary_model <- stats::lm(
          formula = stats::reformulate(
            termlabels = other_predictors,
            response = variable_name
          ),
          data = data
        )
        
        tibble::tibble(
          variable = variable_name,
          vif = 1 / (
            1 -
              summary(auxiliary_model)$r.squared
          )
        )
      }
    )
  )
}


# 4.2.3 Calculate detailed VIF values for each HFI metric x model structure ----
vif_detailed_results_60 <- dplyr::bind_rows(
  lapply(
    selected_hfi_candidates_60,
    function(hfi_name) {
      
      dplyr::bind_rows(
        lapply(
          names(vif_model_structures_60),
          function(model_name) {
            
            model_predictors <- c(
              hfi_name,
              vif_model_structures_60[[model_name]]
            )
            
            calculate_vif_60(
              data = aerial_transitions_60,
              predictors = model_predictors
            ) %>%
              dplyr::mutate(
                hfi_variable = hfi_name,
                model_structure = model_name,
                predictor_role = dplyr::if_else(
                  variable == hfi_name,
                  "HFI",
                  "control"
                )
              )
          }
        )
      )
    }
  )
) %>%
  dplyr::select(
    hfi_variable,
    model_structure,
    variable,
    predictor_role,
    vif
  ) %>%
  dplyr::arrange(
    hfi_variable,
    model_structure,
    dplyr::desc(vif)
  )

# 4.2.4 Summarise VIF at the model level ----
vif_summary_results_60 <- vif_detailed_results_60 %>%
  dplyr::group_by(
    hfi_variable,
    model_structure
  ) %>%
  dplyr::summarise(
    hfi_vif = vif[
      predictor_role == "HFI"
    ][1],
    
    maximum_vif = max(vif),
    
    variable_with_maximum_vif = variable[
      which.max(vif)
    ],
    
    .groups = "drop"
  ) %>%
  dplyr::arrange(
    hfi_variable,
    model_structure
  )

print(vif_summary_results_60,n = Inf)


# 4.2.5 Define readable labels for image plots ----
vif_variable_labels_60 <- c(
  hfi_mean_1000m = "HFI mean\n1000 m",
  hfi_q75_1000m = "HFI Q75\n1000 m",
  hfi_q75_500m = "HFI Q75\n500 m",
  dem_elevation = "Elevation",
  ruggedness_100m = "Ruggedness",
  distance_to_ridgeline_100m = "Distance to\nridgeline",
  prop_forest_5cells = "Forest\nproportion",
  prop_low_vegetation_5cells = "Low vegetation\nproportion"
)

vif_structure_labels_60 <- c(
  altitudinal_gradient = "Altitudinal\ngradient",
  topographic_complexity = "Topographic\ncomplexity",
  habitat_refuge = "Habitat\nrefuge",
  open_habitat = "Open\nhabitat"
)

vif_hfi_labels_60 <- c(
  hfi_mean_1000m = "HFI mean 1000 m",
  hfi_q75_1000m = "HFI Q75 1000 m",
  hfi_q75_500m = "HFI Q75 500 m"
)


# 4.2.6 Prepare detailed VIF table for plotting ----
vif_detailed_plot_data_60 <- vif_detailed_results_60 %>%
  dplyr::mutate(
    variable_label = unname(
      vif_variable_labels_60[
        variable
      ]
    ),
    
    model_label = paste(
      unname(
        vif_hfi_labels_60[
          hfi_variable
        ]
      ),
      unname(
        vif_structure_labels_60[
          model_structure
        ]
      ),
      sep = "\n"
    ),
    
    variable_label = factor(
      variable_label,
      levels = rev(
        unique(
          unname(
            vif_variable_labels_60
          )
        )
      )
    ),
    
    model_label = factor(
      model_label,
      levels = unique(model_label)
    )
  )


#------------------------------------------------------------------------------- Visualisation n°3: the detailed VIF table
vif_detailed_table_plot_60 <- ggplot2::ggplot(
  vif_detailed_plot_data_60,
  ggplot2::aes(
    x = model_label,
    y = variable_label,
    fill = vif
  )
) +
  ggplot2::geom_tile(
    color = "white"
  ) +
  ggplot2::geom_text(
    ggplot2::aes(
      label = sprintf(
        "%.2f",
        vif
      )
    ),
    size = 3.2
  ) +
  ggplot2::scale_fill_gradient(
    low = "white",
    high = "darkseagreen3"
  ) +
  ggplot2::labs(
    x = NULL,
    y = NULL,
    fill = "VIF",
    title = "Variance inflation factors across candidate HFI models",
    subtitle = "Detailed VIF values for each HFI metric and environmental model structure"
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(
      angle = 45,
      hjust = 1,
      vjust = 1
    ),
    panel.grid = ggplot2::element_blank()
  )

print(vif_detailed_table_plot_60)
# White cells indicate variables absent from a given model structure.
# Darker cells indicate larger VIF values.
# all HFI metric have an IVF value below 2, and only three models have a maximum 
# VIF value below 3. Their is therefore no major concern regarding collinearity.
#------------------------------------------------------------------------------- (end) Visualisation n°3: the detailed VIF table




#------------------------------------------------------------------------------ STEP 5: fit candidate HFI models ----
#' **Steps:**
#' (i) centre and standardize HFI and environmental covariates;
#' (ii) fit four environmental models for each retained HFI metric;
#' (iii) estimate Q05-Q95 contrasts using randomly sampled observed covariate
#'       combinations rather than fixing covariates at their means;
#' (iv) compile and plot the model results in one table.



# 5.1 Centre and standardize model covariates ----
# 5.1.1 Define covariates requiring standardization ----
# age_z and duration_z were already standardized in STEP 2.
# cos_diel and sin_time are retained on their original [-1, 1] scale.

environmental_model_covariates_60 <- c(
  "dem_elevation",
  "ruggedness_100m",
  "distance_to_ridgeline_100m",
  "prop_forest_5cells",
  "prop_low_vegetation_5cells"
)

variables_to_standardize_60 <- c(
  selected_hfi_candidates_60,
  environmental_model_covariates_60
)


# 5.1.2 Centre and standardize HFI and environmental covariates ----
# scale() subtracts the empirical mean and divides by the empirical SD.

hfi_model_data_60 <- backbone_data_60 %>%
  dplyr::mutate(
    dplyr::across(
      dplyr::all_of(
        variables_to_standardize_60
      ),
      ~ as.numeric(
        scale(.x)
      ),
      .names = "{.col}_z"
    )
  )

standardized_model_covariates_60 <- paste0(
  variables_to_standardize_60,
  "_z"
)


# 5.1.3 Verify centering and standardization ----
standardization_check_60 <- dplyr::bind_rows(
  lapply(
    standardized_model_covariates_60,
    function(variable_name) {
      
      tibble::tibble(
        variable = variable_name,
        
        mean = mean(
          hfi_model_data_60[[
            variable_name
          ]]
        ),
        
        standard_deviation = stats::sd(
          hfi_model_data_60[[
            variable_name
          ]]
        )
      )
    }
  )
)

print(
  standardization_check_60,
  n = Inf
)

# Inspect:
# - mean should be approximately 0;
# - standard_deviation should be approximately 1.


# 5.1.4 Compile raw and standardized HFI thresholds ----
hfi_model_thresholds_60 <- dplyr::bind_rows(
  lapply(
    selected_hfi_candidates_60,
    function(hfi_name) {
      
      hfi_z_name <- paste0(
        hfi_name,
        "_z"
      )
      
      tibble::tibble(
        hfi_variable = hfi_name,
        
        hfi_q05 = as.numeric(
          stats::quantile(
            hfi_model_data_60[[
              hfi_name
            ]],
            probs = 0.05
          )
        ),
        
        hfi_q95 = as.numeric(
          stats::quantile(
            hfi_model_data_60[[
              hfi_name
            ]],
            probs = 0.95
          )
        ),
        
        hfi_q05_z = as.numeric(
          stats::quantile(
            hfi_model_data_60[[
              hfi_z_name
            ]],
            probs = 0.05
          )
        ),
        
        hfi_q95_z = as.numeric(
          stats::quantile(
            hfi_model_data_60[[
              hfi_z_name
            ]],
            probs = 0.95
          )
        )
      )
    }
  )
)

print(
  hfi_model_thresholds_60,
  n = Inf
)


#------------------------------------------------------------------------------
# 5.2 Define and fit the 12 candidate models ----


# 5.2.1 Define environmental model structures ----
# All terms refer to centred and standardized covariates.

hfi_model_structures_60 <- list(
  
  elevation = c(
    "dem_elevation_z"
  ),
  
  topographic_complexity = c(
    "dem_elevation_z",
    "ruggedness_100m_z",
    "distance_to_ridgeline_100m_z"
  ),
  
  habitat_refuge = c(
    "dem_elevation_z",
    "prop_forest_5cells_z",
    "ruggedness_100m_z"
  ),
  
  open_habitat = c(
    "dem_elevation_z",
    "prop_low_vegetation_5cells_z"
  )
)


# 5.2.2 Define the 12 model specifications ----
hfi_model_metadata_60 <- tidyr::expand_grid(
  hfi_variable =
    selected_hfi_candidates_60,
  
  model_type =
    names(
      hfi_model_structures_60
    )
) %>%
  dplyr::mutate(
    model = paste(
      hfi_variable,
      model_type,
      sep = "__"
    )
  )


# 5.2.3 Add each HFI and environmental structure to the backbone ----
selected_backbone_formula_60 <- stats::formula(
  provisional_backbone_model_ml_60
)

hfi_candidate_formulas_60 <- stats::setNames(
  lapply(
    seq_len(
      nrow(
        hfi_model_metadata_60
      )
    ),
    function(i) {
      
      model_information <-
        hfi_model_metadata_60[i, ]
      
      hfi_term <- paste0(
        model_information$hfi_variable,
        "_z"
      )
      
      additional_terms <- c(
        hfi_term,
        hfi_model_structures_60[[
          model_information$model_type
        ]]
      )
      
      stats::update.formula(
        selected_backbone_formula_60,
        
        stats::as.formula(
          paste(
            ". ~ . +",
            paste(
              additional_terms,
              collapse = " + "
            )
          )
        )
      )
    }
  ),
  hfi_model_metadata_60$model
)


# 5.2.4 Fit all candidate models using ML ----
# The same observations and estimation method are used for every model.

hfi_candidate_models_ml_60 <- stats::setNames(
  lapply(
    names(
      hfi_candidate_formulas_60
    ),
    function(model_name) {
      
      mgcv::gam(
        formula =
          hfi_candidate_formulas_60[[
            model_name
          ]],
        
        data =
          hfi_model_data_60,
        
        family =
          stats::binomial(
            link = "logit"
          ),
        
        method =
          gam_selection_method
      )
    }
  ),
  names(
    hfi_candidate_formulas_60
  )
)


#------------------------------------------------------------------------------
# 5.3 Draw observed covariate combinations for Q05-Q95 predictions ----


# 5.3.1 Draw the same number of observed covariate combinations per individual ----
# Sampling complete rows preserves the observed associations among age, diel
# time, aerial duration, topography and habitat.
#
# Each individual contributes the same number of sampled rows and therefore
# receives the same weight in population-level predictions.

prediction_seed_60 <- 20260720
n_prediction_draws_per_individual_60 <- 200L

set.seed(
  prediction_seed_60
)

prediction_reference_data_60 <- hfi_model_data_60 %>%
  dplyr::group_by(
    individual_id
  ) %>%
  dplyr::slice_sample(
    n =
      n_prediction_draws_per_individual_60,
    
    replace =
      TRUE
  ) %>%
  dplyr::ungroup()

prediction_sample_summary_60 <- prediction_reference_data_60 %>%
  dplyr::count(
    individual_id,
    name =
      "n_prediction_rows"
  )

print(
  prediction_sample_summary_60,
  n = Inf
)

# Every individual should have exactly n_prediction_draws_per_individual_60 rows.

#------------------------------------------------------------------------------
# 5.4 Extract fit statistics and Q05-Q95 contrasts ----


# 5.4.1 Calculate model results ----
hfi_model_results_60 <- dplyr::bind_rows(
  lapply(
    names(
      hfi_candidate_models_ml_60
    ),
    function(model_name) {
      
      fitted_model <-
        hfi_candidate_models_ml_60[[
          model_name
        ]]
      
      model_information <-
        hfi_model_metadata_60 %>%
        dplyr::filter(
          model ==
            model_name
        )
      
      hfi_name <-
        model_information$hfi_variable
      
      hfi_term <- paste0(
        hfi_name,
        "_z"
      )
      
      hfi_thresholds <-
        hfi_model_thresholds_60 %>%
        dplyr::filter(
          hfi_variable ==
            hfi_name
        )
      
      prediction_data_q05 <-
        prediction_reference_data_60
      
      prediction_data_q95 <-
        prediction_reference_data_60
      
      prediction_data_q05[[
        hfi_term
      ]] <- hfi_thresholds$hfi_q05_z
      
      prediction_data_q95[[
        hfi_term
      ]] <- hfi_thresholds$hfi_q95_z
      
      # Population-level linear predictors:
      # the individual random intercept is fixed at zero.
      linear_predictor_q05 <- stats::predict(
        fitted_model,
        newdata =
          prediction_data_q05,
        type =
          "link",
        exclude =
          "s(individual_id)"
      )
      
      linear_predictor_q95 <- stats::predict(
        fitted_model,
        newdata =
          prediction_data_q95,
        type =
          "link",
        exclude =
          "s(individual_id)"
      )
      
      probability_q05 <- stats::predict(
        fitted_model,
        newdata =
          prediction_data_q05,
        type =
          "response",
        exclude =
          "s(individual_id)"
      )
      
      probability_q95 <- stats::predict(
        fitted_model,
        newdata =
          prediction_data_q95,
        type =
          "response",
        exclude =
          "s(individual_id)"
      )
      
      coefficient_table <-
        summary(
          fitted_model
        )$p.table
      
      hfi_estimate <-
        coefficient_table[
          hfi_term,
          1
        ]
      
      hfi_standard_error <-
        coefficient_table[
          hfi_term,
          2
        ]
      
      hfi_standardized_range <-
        hfi_thresholds$hfi_q95_z -
        hfi_thresholds$hfi_q05_z
      
      log_odds_difference <-
        hfi_estimate *
        hfi_standardized_range
      
      log_odds_difference_se <-
        abs(
          hfi_standardized_range
        ) *
        hfi_standard_error
      
      # Calculate predictions separately for each individual.
      individual_prediction_summary_60 <- tibble::tibble(
        individual_id =
          prediction_reference_data_60$
          individual_id,
        
        linear_predictor_q05 =
          as.numeric(
            linear_predictor_q05
          ),
        
        linear_predictor_q95 =
          as.numeric(
            linear_predictor_q95
          ),
        
        probability_q05 =
          as.numeric(
            probability_q05
          ),
        
        probability_q95 =
          as.numeric(
            probability_q95
          )
      ) %>%
        dplyr::group_by(
          individual_id
        ) %>%
        dplyr::summarise(
          mean_linear_predictor_q05 =
            mean(
              linear_predictor_q05
            ),
          
          mean_linear_predictor_q95 =
            mean(
              linear_predictor_q95
            ),
          
          mean_probability_q05 =
            mean(
              probability_q05
            ),
          
          mean_probability_q95 =
            mean(
              probability_q95
            ),
          
          .groups =
            "drop"
        ) %>%
        dplyr::mutate(
          log_odds_difference_q95_q05 =
            mean_linear_predictor_q95 -
            mean_linear_predictor_q05,
          
          absolute_probability_difference =
            mean_probability_q95 -
            mean_probability_q05,
          
          relative_probability_difference_percent =
            100 *
            absolute_probability_difference /
            mean_probability_q05
        )
      
      # Equal-weight population summaries:
      # each individual contributes one value to each mean.
      mean_probability_q05 <-
        mean(
          individual_prediction_summary_60$
            mean_probability_q05
        )
      
      mean_probability_q95 <-
        mean(
          individual_prediction_summary_60$
            mean_probability_q95
        )
      
      log_odds_difference <-
        mean(
          individual_prediction_summary_60$
            log_odds_difference_q95_q05
        )
      
      absolute_probability_difference <-
        mean(
          individual_prediction_summary_60$
            absolute_probability_difference
        )
      
      relative_probability_difference_percent <-
        mean(
          individual_prediction_summary_60$
            relative_probability_difference_percent
        )
      
      tibble::tibble(
        model =
          model_name,
        
        hfi_variable =
          hfi_name,
        
        model_type =
          model_information$model_type,
        
        n_observations =
          nrow(
            stats::model.frame(
              fitted_model
            )
          ),
        
        n_prediction_individuals =
          nrow(
            individual_prediction_summary_60
          ),
        
        hfi_q05 =
          hfi_thresholds$hfi_q05,
        
        hfi_q95 =
          hfi_thresholds$hfi_q95,
        
        mean_probability_q05 =
          mean_probability_q05,
        
        mean_probability_q95 =
          mean_probability_q95,
        
        log_odds_difference_q95_q05 =
          log_odds_difference,
        
        log_odds_difference_se =
          log_odds_difference_se,
        
        absolute_probability_difference =
          absolute_probability_difference,
        
        relative_probability_difference_percent =
          relative_probability_difference_percent,
        
        model_edf =
          sum(
            fitted_model$edf
          ),
        
        AIC =
          stats::AIC(
            fitted_model
          ),
        
        converged =
          fitted_model$converged
      )
    }
  )
) %>%
  dplyr::mutate(
    delta_AIC =
      AIC -
      min(AIC)
  ) %>%
  dplyr::arrange(
    factor(
      hfi_variable,
      levels =
        selected_hfi_candidates_60
    ),
    factor(
      model_type,
      levels =
        names(
          hfi_model_structures_60
        )
    )
  )

print(
  hfi_model_results_60,
  n = Inf
)


# Inspect:
# - converged: expected TRUE;
# - n_observations: must be identical for all 12 models;
# - log_odds_difference_q95_q05: signed Q95-Q05 HFI contrast;
# - relative_probability_difference_percent: relative change in predicted
#   probability, averaged across sampled observed covariate combinations;
# - delta_AIC: relative model support among the 12 candidate models.


#------------------------------------------------------------------------------
# 5.5 Compile the requested summary table ----


# 5.5.1 Define readable model labels ----
hfi_result_labels_60 <- c(
  hfi_mean_1000m =
    "HFI mean 1000 m",
  
  hfi_q75_1000m =
    "HFI Q75 1000 m",
  
  hfi_q75_500m =
    "HFI Q75 500 m"
)

model_type_labels_60 <- c(
  elevation =
    "Elevation",
  
  topographic_complexity =
    "Topographic complexity",
  
  habitat_refuge =
    "Refuge",
  
  open_habitat =
    "Open habitat"
)


# 5.5.2 Create the formatted results table ----
hfi_model_summary_table_60 <- hfi_model_results_60 %>%
  dplyr::transmute(
    `HFI metric` =
      unname(
        hfi_result_labels_60[
          hfi_variable
        ]
      ),
    
    `Model type` =
      unname(
        model_type_labels_60[
          model_type
        ]
      ),
    
    `Q95 - Q05\nlog-odds` =
      sprintf(
        "%.3f",
        log_odds_difference_q95_q05
      ),
    
    `Relative Q95 vs Q05\nprobability difference (%)` =
      sprintf(
        "%.2f",
        relative_probability_difference_percent
      ),
    
    `Model EDF` =
      sprintf(
        "%.2f",
        model_edf
      ),
    
    `Contrast SE` =
      sprintf(
        "%.3f",
        log_odds_difference_se
      ),
    
    `AIC` =
      sprintf(
        "%.1f",
        AIC
      )
  ) %>%
  dplyr::group_by(
    `HFI metric`
  ) %>%
  dplyr::mutate(
    `HFI metric` =
      dplyr::if_else(
        dplyr::row_number() == 1L,
        `HFI metric`,
        ""
      )
  ) %>%
  dplyr::ungroup()

print(
  hfi_model_summary_table_60,
  n = Inf
)


#------------------------------------------------------------------------------
# 5.6 Plot the results table as an image ----


# 5.6.1 Reshape the formatted table for plotting ----
hfi_model_table_plot_data_60 <-
  hfi_model_summary_table_60 %>%
  dplyr::mutate(
    row_id =
      dplyr::row_number(),
    
    row_fill =
      dplyr::if_else(
        row_id %% 2L == 0L,
        "grey96",
        "white"
      )
  ) %>%
  tidyr::pivot_longer(
    cols =
      -c(
        row_id,
        row_fill
      ),
    
    names_to =
      "column",
    
    values_to =
      "value"
  )

table_column_order_60 <- names(
  hfi_model_summary_table_60
)

hfi_model_table_plot_data_60 <-
  hfi_model_table_plot_data_60 %>%
  dplyr::mutate(
    column = factor(
      column,
      levels =
        table_column_order_60
    ),
    
    row_id = factor(
      row_id,
      levels =
        rev(
          seq_len(
            nrow(
              hfi_model_summary_table_60
            )
          )
        )
    )
  )


# 5.6.2 Plot table ----
hfi_model_summary_table_plot_60 <- ggplot2::ggplot(
  hfi_model_table_plot_data_60,
  ggplot2::aes(
    x =
      column,
    y =
      row_id
  )
) +
  ggplot2::geom_tile(
    ggplot2::aes(
      fill =
        row_fill
    ),
    color =
      "grey80"
  ) +
  ggplot2::scale_fill_identity() +
  ggplot2::geom_text(
    ggplot2::aes(
      label =
        value
    ),
    size =
      3.2
  ) +
  ggplot2::scale_x_discrete(
    position =
      "top"
  ) +
  ggplot2::labs(
    x =
      NULL,
    y =
      NULL,
    title =
      "Comparison of candidate HFI transition models",
    subtitle =
      "Q05-Q95 contrasts use randomly sampled observed covariate combinations"
  ) +
  ggplot2::theme_void() +
  ggplot2::theme(
    axis.text.x =
      ggplot2::element_text(
        size = 9,
        face = "bold",
        margin = ggplot2::margin(
          b = 8
        )
      ),
    
    plot.title =
      ggplot2::element_text(
        face = "bold",
        size = 14
      ),
    
    plot.subtitle =
      ggplot2::element_text(
        size = 10
      ),
    
    plot.margin =
      ggplot2::margin(
        15,
        15,
        15,
        15
      )
  )

print(hfi_model_summary_table_plot_60)







#------------------------------------------------------------------------------- STEP 6: control spatio temporal autocorrtion
#' **Steps:**
#' (i) extract HFI estimates and Q05-Q95 contrasts;
#' (ii) compare model support using AIC;
#' (iii) print one final formatted comparison table.


# 6.1 Define the four candidate models ----
# for the 60 minutes dataset
selected_models_60 <- c("hfi_mean_1000m__elevation","hfi_mean_1000m__open_habitat","hfi_q75_1000m__open_habitat","hfi_q75_500m__open_habitat")

# for the 20 minutes dataset
selected_models_60 <- c("hfi_mean_1000m__elevation","hfi_mean_1000m__open_habitat","hfi_mean_1000m__habitat_refuge","hfi_mean_1000m__topographic_complexity")

# 6.2 Extract information required for the final table ----
model_comparison_60 <- dplyr::bind_rows(
  lapply(
    selected_models_60,
    function(model_name) {
      
      fitted_model <-
        hfi_candidate_models_ml_60[[
          model_name
        ]]
      
      model_info <-
        hfi_model_results_60 %>%
        dplyr::filter(
          model == model_name
        ) %>%
        dplyr::slice(1)
      
      if (
        nrow(
          model_info
        ) != 1L
      ) {
        stop(
          paste0(
            "Exactly one row was expected in ",
            "hfi_model_results_60 for model: ",
            model_name
          )
        )
      }
      
      model_summary <-
        summary(
          fitted_model
        )
      
      hfi_variable <-
        model_info$hfi_variable[[1]]
      
      hfi_term <-
        paste0(
          hfi_variable,
          "_z"
        )
      
      if (
        !hfi_term %in%
        rownames(
          model_summary$p.table
        )
      ) {
        stop(
          paste0(
            "Coefficient ",
            hfi_term,
            " was not found in model ",
            model_name
          )
        )
      }
      
      hfi_estimate <-
        model_summary$p.table[
          hfi_term,
          "Estimate"
        ]
      
      hfi_se <-
        model_summary$p.table[
          hfi_term,
          "Std. Error"
        ]
      
      contrast_estimate <-
        model_info$
        log_odds_difference_q95_q05[[1]]
      
      contrast_se <-
        model_info$
        log_odds_difference_se[[1]]
      
      tibble::tibble(
        
        # This column was missing in your previous version.
        # It is required for all subsequent joins.
        model =
          model_name,
        
        hfi_variable =
          hfi_variable,
        
        model_type =
          model_info$model_type[[1]],
        
        AIC =
          stats::AIC(
            fitted_model
          ),
        
        hfi_estimate =
          hfi_estimate,
        
        hfi_standard_error =
          hfi_se,
        
        hfi_confidence_low =
          hfi_estimate -
          1.96 *
          hfi_se,
        
        hfi_confidence_high =
          hfi_estimate +
          1.96 *
          hfi_se,
        
        contrast_estimate =
          contrast_estimate,
        
        contrast_standard_error =
          contrast_se,
        
        contrast_confidence_low =
          contrast_estimate -
          1.96 *
          contrast_se,
        
        contrast_confidence_high =
          contrast_estimate +
          1.96 *
          contrast_se,
        
        absolute_difference_pp =
          100 *
          model_info$
          absolute_probability_difference[[1]],
        
        relative_difference_percent =
          model_info$
          relative_probability_difference_percent[[1]]
      )
    }
  )
) %>%
  dplyr::arrange(
    AIC
  ) %>%
  dplyr::mutate(
    delta_AIC =
      AIC -
      min(
        AIC
      )
  )


print(
  model_comparison_60,
  n = Inf
)

first_model_name_60 <- model_comparison_60$model[1]

first_model_60 <- hfi_candidate_models_ml_60[[first_model_name_60]]

summary(
  first_model_60
)

# 6.3 Format and print the final comparison table ----
model_comparison_formatted_60 <- model_comparison_60 %>%
  dplyr::mutate(
    hfi_metric_label =
      dplyr::recode(
        hfi_variable,
        hfi_mean_1000m =
          "HFI mean 1000 m",
        hfi_q75_1000m =
          "HFI Q75 1000 m",
        hfi_q75_500m =
          "HFI Q75 500 m"
      ),
    
    model_type_label =
      dplyr::recode(
        model_type,
        elevation =
          "Elevation",
        open_habitat =
          "Open habitat"
      ),
    
    hfi_effect_label =
      sprintf(
        "%.3f [%.3f, %.3f]",
        hfi_estimate,
        hfi_confidence_low,
        hfi_confidence_high
      ),
    
    contrast_label =
      sprintf(
        "%.3f [%.3f, %.3f]",
        contrast_estimate,
        contrast_confidence_low,
        contrast_confidence_high
      )
  ) %>%
  dplyr::transmute(
    `HFI metric` =
      hfi_metric_label,
    
    `Model type` =
      model_type_label,
    
    `AIC` =
      round(
        AIC,
        1
      ),
    
    `Delta AIC` =
      round(
        delta_AIC,
        2
      ),
    
    `HFI coefficient [95% CI]` =
      hfi_effect_label,
    
    `Q95-Q05 log-odds [95% CI]` =
      contrast_label,
    
    `Absolute difference (percentage points)` =
      round(
        absolute_difference_pp,
        2
      ),
    
    `Relative difference (%)` =
      round(
        relative_difference_percent,
        2
      )
  )

print(model_comparison_formatted_60,n = Inf)
# HFI metric`    `Model type`    AIC `Delta AIC` `HFI coefficient [95% CI]` `Q95-Q05 log-odds [95% CI]` `Absolute difference (percentage points)` `Relative difference (%)`
# HFI mean 1000 m Open habitat 14115.        0    0.035 [-0.021, 0.091]      0.104 [-0.063, 0.272]                                            2.41                      6.33
# HFI Q75 1000 m  Open habitat 14116.        1.9  0.009 [-0.047, 0.065]      0.027 [-0.145, 0.200]                                            0.63                      1.63
# HFI Q75 500 m   Open habitat 14117.        1.96 0.008 [-0.045, 0.061]      0.023 [-0.131, 0.178]                                            0.54                      1.39
# HFI mean 1000 m Elevation    14120.        5.8  0.037 [-0.019, 0.093]      0.110 [-0.058, 0.277]                                            2.53                      6.67


#' For 20 minutes dataset 
#' FI metric`    `Model type`    AIC `Delta AIC` HFI coefficient [95%…¹ Q95-Q05 log-odds [95…² Absolute difference …³ Relative difference …⁴
#' HFI mean 1000 m topographic… 29987.        0    0.083 [0.042, 0.124]   0.246 [0.125, 0.367]                     5.91                   11.9
#' HFI mean 1000 m Open habitat 29989.        2.34 0.085 [0.044, 0.125]   0.251 [0.131, 0.370]                     6.02                   12.2
#' HFI mean 1000 m habitat_ref… 29993.        6.23 0.082 [0.041, 0.123]   0.241 [0.120, 0.362]                     5.8                    11.7
#' HFI mean 1000 m Elevation    29995.        7.78 0.089 [0.048, 0.129]   0.262 [0.143, 0.382]                     6.3                    12.8
#' HFI mean 1000m and q75 1000m have the lowest AIC, HFI coefficient are also stable (0.037 and 0.035) across elevation and open habitat models. Given that low vegetation 
#' increase the AIC by 5.8 units while not modifying the HFI coefficient substantially compare to the model with elevation only, we should prefer open habitat hfi 1000m mean. 


#------------------------------------------------------------------------------ STEP 7: diagnose residual spatio-temporal autocorrelation ----
#' **Philosophy**: we investigate whether the residual of the models have an 
#' underlying spatio temporal structure. 
#' **Steps:**
#' (i) simulate DHARMa residuals for the four candidate models. DHARMa compare real
#' observation to simulated observation generated by the four models. 
#' (ii) temporal series are rebuild and calculate the correlation at the different
#' lag (60, 120 and 180 minutes). Bootstrap across individuals (on tire les aigles
#' avec remise et on recalcule les correlation)
#' (iii) test spatial autocorrelation after aggregation in 1-km cells. We do the Moran's I 
#' test which test whether spatially closed residuals are more simalar than expected 
#' under a randomely generated distribution.
#' (iv) add the diagnostics to the candidate-model comparison table.
#' **Results**: their is no spatio temporal correlation detectable after modelisation.

#------------------------------------------------------------------------------
# 7.1 Parameters
#------------------------------------------------------------------------------

n_dharma_simulations_60 <- 1000L

max_temporal_lag_60 <- 5L

nominal_interval_min_60 <- 60

temporal_tolerance_min_60 <- 15

n_temporal_bootstrap_60 <- 999L

spatial_grid_size_m_60 <- 1000


#------------------------------------------------------------------------------
# 7.2 Simulate DHARMa residuals
#------------------------------------------------------------------------------

dharma_residuals_60 <- stats::setNames(
  
  lapply(
    selected_models_60,
    
    function(model_name){
      
      DHARMa::simulateResiduals(
        
        fittedModel =
          hfi_candidate_models_ml_60[[model_name]],
        
        n =
          n_dharma_simulations_60,
        
        refit =
          FALSE
      )
    }
  ),
  
  selected_models_60
)



#------------------------------------------------------------------------------
# 7.3 Temporal autocorrelation
#------------------------------------------------------------------------------


build_temporal_pairs_60 <- function(
    dharma_output,
    model_data,
    max_lag = 5L,
    nominal_interval_min = 60,
    tolerance_min = 15
){
  
  stopifnot(
    length(dharma_output$scaledResiduals) ==
      nrow(model_data)
  )
  
  
  residual_data <- model_data %>%
    
    dplyr::mutate(
      
      residual_normal =
        stats::qnorm(
          pmin(
            pmax(
              dharma_output$scaledResiduals,
              1e-6
            ),
            1 - 1e-6
          )
        )
    ) %>%
    
    dplyr::arrange(
      individual_id,
      timestamp
    ) %>%
    
    dplyr::group_by(
      individual_id
    )
  
  
  dplyr::bind_rows(
    
    lapply(
      seq_len(max_lag),
      
      function(lag_value){
        
        
        residual_data %>%
          
          dplyr::mutate(
            
            residual_previous =
              dplyr::lag(
                residual_normal,
                lag_value
              ),
            
            timestamp_previous =
              dplyr::lag(
                timestamp,
                lag_value
              ),
            
            time_lag_min =
              as.numeric(
                difftime(
                  timestamp,
                  timestamp_previous,
                  units="mins"
                )
              )
          ) %>%
          
          dplyr::ungroup() %>%
          
          dplyr::filter(
            !is.na(residual_previous),
            
            abs(
              time_lag_min -
                lag_value * nominal_interval_min
            ) <=
              tolerance_min *
              sqrt(lag_value)
          ) %>%
          
          dplyr::transmute(
            
            individual_id,
            
            lag =
              lag_value,
            
            time_lag_min,
            
            residual_previous,
            
            residual_current =
              residual_normal
          )
      }
    )
  )
}



summarise_temporal_lags_60 <- function(
    temporal_pairs,
    n_bootstrap = 999L
){
  
  dplyr::bind_rows(
    
    lapply(
      sort(unique(temporal_pairs$lag)),
      
      function(lag_value){
        
        lag_data <-
          temporal_pairs %>%
          dplyr::filter(
            lag == lag_value
          )
        
        
        observed_correlation <-
          stats::cor(
            lag_data$residual_previous,
            lag_data$residual_current
          )
        
        
        individual_clusters <-
          split(
            lag_data,
            lag_data$individual_id
          )
        
        
        set.seed(
          20260721 + lag_value
        )
        
        
        bootstrap_correlations <-
          replicate(
            
            n_bootstrap,
            
            {
              sampled_clusters <-
                sample(
                  seq_along(individual_clusters),
                  length(individual_clusters),
                  replace=TRUE
                )
              
              
              bootstrap_data <-
                dplyr::bind_rows(
                  individual_clusters[
                    sampled_clusters
                  ]
                )
              
              
              stats::cor(
                bootstrap_data$residual_previous,
                bootstrap_data$residual_current
              )
            }
          )
        
        
        tibble::tibble(
          
          lag =
            lag_value,
          
          mean_time_lag_min =
            mean(
              lag_data$time_lag_min
            ),
          
          residual_correlation =
            observed_correlation,
          
          confidence_low =
            quantile(
              bootstrap_correlations,
              0.025
            ),
          
          confidence_high =
            quantile(
              bootstrap_correlations,
              0.975
            ),
          
          n_pairs =
            nrow(lag_data),
          
          n_individuals =
            dplyr::n_distinct(
              lag_data$individual_id
            )
        )
      }
    )
  )
}



#------------------------------------------------------------------------------
# 7.4 Run temporal diagnostics for each model
#------------------------------------------------------------------------------


temporal_lag_diagnostics_60 <-
  
  dplyr::bind_rows(
    
    lapply(
      selected_models_60,
      
      function(model_name){
        
        
        build_temporal_pairs_60(
          
          dharma_output =
            dharma_residuals_60[[model_name]],
          
          model_data =
            hfi_model_data_60,
          
          max_lag =
            max_temporal_lag_60,
          
          nominal_interval_min =
            nominal_interval_min_60,
          
          tolerance_min =
            temporal_tolerance_min_60
        ) %>%
          
          summarise_temporal_lags_60(
            n_bootstrap =
              n_temporal_bootstrap_60
          ) %>%
          
          dplyr::mutate(
            model =
              model_name,
            .before=1
          )
      }
    )
  )


print(
  temporal_lag_diagnostics_60,
  n=Inf
)



#------------------------------------------------------------------------------
# 7.5 Spatial autocorrelation
#------------------------------------------------------------------------------


test_spatial_autocorrelation_60 <- function(
    dharma_output,
    model_data,
    grid_size_m = 1000
){
  
  
  spatial_data <-
    model_data %>%
    
    dplyr::mutate(
      
      residual =
        dharma_output$scaledResiduals,
      
      cell_x =
        floor(
          x_3035 / grid_size_m
        ),
      
      cell_y =
        floor(
          y_3035 / grid_size_m
        ),
      
      spatial_cell =
        paste(
          cell_x,
          cell_y,
          sep="_"
        )
    )
  
  
  spatial_cells <-
    spatial_data %>%
    
    dplyr::group_by(
      spatial_cell
    ) %>%
    
    dplyr::summarise(
      
      x =
        mean(x_3035),
      
      y =
        mean(y_3035),
      
      .groups="drop"
    )
  
  
  grouped_residuals <-
    DHARMa::recalculateResiduals(
      dharma_output,
      group =
        factor(
          spatial_data$spatial_cell
        )
    )
  
  
  test <-
    DHARMa::testSpatialAutocorrelation(
      grouped_residuals,
      x =
        spatial_cells$x,
      y =
        spatial_cells$y,
      plot =
        FALSE
    )
  
  
  tibble::tibble(
    
    spatial_p_value =
      test$p.value,
    
    spatial_moran_I =
      test$observed,
    
    n_spatial_cells =
      nrow(spatial_cells)
  )
}



#------------------------------------------------------------------------------
# 7.6 Run spatial diagnostics
#------------------------------------------------------------------------------


spatial_diagnostics_60 <-
  
  dplyr::bind_rows(
    
    lapply(
      selected_models_60,
      
      function(model_name){
        
        
        test_spatial_autocorrelation_60(
          
          dharma_output =
            dharma_residuals_60[[model_name]],
          
          model_data =
            hfi_model_data_60,
          
          grid_size_m =
            spatial_grid_size_m_60
        ) %>%
          
          dplyr::mutate(
            model =
              model_name,
            .before=1
          )
      }
    )
  )


print(
  spatial_diagnostics_60,
  n=Inf
)



#------------------------------------------------------------------------------
# 7.7 Final diagnostic table
#------------------------------------------------------------------------------


model_comparison_diagnostics_60 <-
  
  model_comparison_60 %>%
  
  dplyr::left_join(
    
    temporal_lag_diagnostics_60 %>%
      dplyr::filter(
        lag==1
      ) %>%
      dplyr::select(
        model,
        temporal_lag1 =
          residual_correlation,
        temporal_lag1_low =
          confidence_low,
        temporal_lag1_high =
          confidence_high
      ),
    
    by="model"
  ) %>%
  
  dplyr::left_join(
    
    spatial_diagnostics_60,
    
    by="model"
  )



print(
  model_comparison_diagnostics_60,
  n=Inf
)


#------------------------------------------------------------------------------
# 7.7 Summary table of temporal and spatial autocorrelation diagnostics
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
# Summary table of temporal and spatial autocorrelation diagnostics
#------------------------------------------------------------------------------

autocorrelation_summary_60 <-
  
  temporal_lag_diagnostics_60 %>%
  
  filter(
    lag == 1L
  ) %>%
  
  select(
    
    model,
    
    temporal_autocorrelation =
      residual_correlation,
    
    temporal_confidence_low =
      confidence_low,
    
    temporal_confidence_high =
      confidence_high
  ) %>%
  
  left_join(
    
    spatial_diagnostics_60,
    
    by = "model"
  ) %>%
  
  mutate(
    
    temporal_autocorrelation_detected =
      if_else(
        temporal_confidence_low > 0 |
          temporal_confidence_high < 0,
        "yes",
        "no"
      ),
    
    spatial_autocorrelation_detected =
      if_else(
        spatial_p_value < 0.05,
        "yes",
        "no"
      )
  )


print(
  autocorrelation_summary_60,
  n = Inf,
  width = Inf
)



#------------------------------------------------------------------------------ STEP 8: estimate equal-weight HFI contrasts with individual bootstrap ----
#' **Steps:**
#' (i) fit the selected HFI mean 1000-m open-habitat model using fast REML;
#' (ii) calculate Q05-Q95 predictions and contrasts for every individual;
#' (iii) calculate equal-weight population summaries;
#' (iv) estimate uncertainty using a faster individual-level cluster bootstrap;
#' (v) report the individual median and proportion of positive contrasts.


# 8.1 Parameters ----
selected_final_hfi_60 <- "hfi_mean_1000m"
selected_final_hfi_z_60 <- "hfi_mean_1000m_z"

# 300 replicates provide a relatively fast exploratory bootstrap.
# Increase to 500-1000 only for the final publication analysis.
n_individual_bootstrap_60 <- 300L
bootstrap_seed_60 <- 20260721L

# Each bootstrap model uses one thread because bootstrap replicates are
# themselves run in parallel.
n_bootstrap_cores_60 <- max(
  1L,
  min(
    4L,
    parallel::detectCores() - 1L
  )
)


# 8.2 Extract fixed Q05 and Q95 HFI thresholds ----
selected_hfi_thresholds_60 <- hfi_model_thresholds_60 %>%
  dplyr::filter(
    hfi_variable == selected_final_hfi_60
  ) %>%
  dplyr::slice(1)

hfi_q05_raw_60 <-
  selected_hfi_thresholds_60$hfi_q05[[1]]

hfi_q95_raw_60 <-
  selected_hfi_thresholds_60$hfi_q95[[1]]

hfi_q05_z_60 <-
  selected_hfi_thresholds_60$hfi_q05_z[[1]]

hfi_q95_z_60 <-
  selected_hfi_thresholds_60$hfi_q95_z[[1]]


# 8.3 Define and fit the selected final model once ----
final_hfi_formula_60 <- remain_aerial ~
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

final_hfi_model_60 <- mgcv::bam(
  formula =
    final_hfi_formula_60,
  
  data =
    hfi_model_data_60,
  
  family =
    stats::binomial(
      link = "logit"
    ),
  
  method =
    "fREML",
  
  discrete =
    TRUE,
  
  nthreads =
    n_bootstrap_cores_60
)

print(
  summary(
    final_hfi_model_60
  )
)


# 8.4 Pre-compute objects reused throughout the bootstrap ----
# The HFI column is fixed at zero in this design matrix.
# Q05 and Q95 are then introduced analytically through the HFI coefficient.
# This avoids creating two prediction datasets and calling predict() twice
# during every bootstrap iteration.

prediction_data_60 <- hfi_model_data_60
prediction_data_60$hfi_mean_1000m_z <- 0

fixed_prediction_matrix_60 <- stats::model.matrix(
  ~ cos_diel +
    sin_time +
    age_z +
    duration_z +
    duration_z2 +
    hfi_mean_1000m_z +
    dem_elevation_z +
    prop_low_vegetation_5cells_z,
  
  data =
    prediction_data_60
)

individual_id_prediction_60 <- droplevels(
  factor(
    hfi_model_data_60$individual_id
  )
)

individual_levels_60 <- levels(
  individual_id_prediction_60
)

individual_index_60 <- as.integer(
  individual_id_prediction_60
)

n_individuals_60 <- length(
  individual_levels_60
)

individual_prediction_sizes_60 <- tabulate(
  individual_index_60,
  nbins =
    n_individuals_60
)

# Split the fitting data only once, rather than rebuilding these subsets
# during every bootstrap iteration.
individual_data_split_60 <- split(
  hfi_model_data_60,
  individual_id_prediction_60,
  drop = TRUE
)

individual_data_sizes_60 <- vapply(
  individual_data_split_60,
  nrow,
  integer(1)
)


# 8.5 Define fast individual Q05-Q95 contrast calculation ----
# Only fixed effects are used for population-level predictions.
# The individual random intercept is therefore not included in these equations.

calculate_fast_individual_contrasts_60 <- function(
    fitted_model
) {
  
  fixed_coefficients <- stats::coef(
    fitted_model
  )[
    colnames(
      fixed_prediction_matrix_60
    )
  ]
  
  hfi_coefficient <-
    fixed_coefficients[[
      selected_final_hfi_z_60
    ]]
  
  # HFI is zero in fixed_prediction_matrix_60.
  baseline_linear_predictor <- as.numeric(
    fixed_prediction_matrix_60 %*%
      fixed_coefficients
  )
  
  probability_q05 <- stats::plogis(
    baseline_linear_predictor +
      hfi_coefficient *
      hfi_q05_z_60
  )
  
  probability_q95 <- stats::plogis(
    baseline_linear_predictor +
      hfi_coefficient *
      hfi_q95_z_60
  )
  
  mean_by_individual <- function(values) {
    
    as.numeric(
      rowsum(
        values,
        individual_index_60,
        reorder = FALSE
      )
    ) /
      individual_prediction_sizes_60
  }
  
  individual_probability_q05 <-
    mean_by_individual(
      probability_q05
    )
  
  individual_probability_q95 <-
    mean_by_individual(
      probability_q95
    )
  
  individual_absolute_difference <-
    individual_probability_q95 -
    individual_probability_q05
  
  list(
    hfi_coefficient =
      unname(
        hfi_coefficient
      ),
    
    probability_q05 =
      individual_probability_q05,
    
    probability_q95 =
      individual_probability_q95,
    
    absolute_difference =
      individual_absolute_difference,
    
    absolute_difference_percentage_points =
      100 *
      individual_absolute_difference,
    
    relative_difference_percent =
      100 *
      individual_absolute_difference /
      individual_probability_q05
  )
}


# 8.6 Calculate observed equal-weight individual contrasts ----
observed_individual_effects_60 <-
  calculate_fast_individual_contrasts_60(
    final_hfi_model_60
  )

individual_hfi_contrasts_60 <- tibble::tibble(
  individual_id =
    individual_levels_60,
  
  mean_probability_q05 =
    observed_individual_effects_60$
    probability_q05,
  
  mean_probability_q95 =
    observed_individual_effects_60$
    probability_q95,
  
  absolute_probability_difference =
    observed_individual_effects_60$
    absolute_difference,
  
  absolute_difference_percentage_points =
    observed_individual_effects_60$
    absolute_difference_percentage_points,
  
  relative_probability_difference_percent =
    observed_individual_effects_60$
    relative_difference_percent,
  
  positive_effect =
    absolute_probability_difference > 0
)

print(
  individual_hfi_contrasts_60,
  n = Inf
)


# 8.7 Calculate equal-weight population estimates ----
equal_weight_hfi_summary_60 <- individual_hfi_contrasts_60 %>%
  dplyr::summarise(
    hfi_q05 =
      hfi_q05_raw_60,
    
    hfi_q95 =
      hfi_q95_raw_60,
    
    n_individuals =
      dplyr::n(),
    
    mean_probability_q05 =
      mean(
        mean_probability_q05
      ),
    
    mean_probability_q95 =
      mean(
        mean_probability_q95
      ),
    
    mean_absolute_difference =
      mean(
        absolute_probability_difference
      ),
    
    mean_absolute_difference_percentage_points =
      mean(
        absolute_difference_percentage_points
      ),
    
    mean_relative_difference_percent =
      mean(
        relative_probability_difference_percent
      ),
    
    median_individual_difference =
      stats::median(
        absolute_probability_difference
      ),
    
    median_individual_difference_percentage_points =
      stats::median(
        absolute_difference_percentage_points
      ),
    
    proportion_individuals_positive =
      mean(
        positive_effect
      ),
    
    hfi_coefficient =
      observed_individual_effects_60$
      hfi_coefficient
  )

print(
  equal_weight_hfi_summary_60
)


# 8.8 Prepare individual cluster-bootstrap draws once ----
# Every replicate still contains n_individuals_60 sampled clusters.
# Some individuals can be selected more than once and others not at all.

set.seed(
  bootstrap_seed_60
)

bootstrap_draws_60 <- replicate(
  n =
    n_individual_bootstrap_60,
  
  expr =
    sample.int(
      n =
        n_individuals_60,
      
      size =
        n_individuals_60,
      
      replace =
        TRUE
    ),
  
  simplify =
    FALSE
)


# 8.9 Define the bootstrap model ----
bootstrap_hfi_formula_60 <- remain_aerial ~
  cos_diel +
  sin_time +
  age_z +
  duration_z +
  duration_z2 +
  hfi_mean_1000m_z +
  dem_elevation_z +
  prop_low_vegetation_5cells_z +
  s(
    bootstrap_individual_id,
    bs = "re"
  )


# 8.10 Define one optimized bootstrap iteration ----
# 8.10 Define one optimized bootstrap iteration ----
run_fast_individual_bootstrap_60 <- function(
    bootstrap_iteration
) {
  
  sampled_indices <-
    bootstrap_draws_60[[
      bootstrap_iteration
    ]]
  
  
  # These weights give one contribution to each sampled eagle copy.
  sampled_individual_weights <- tabulate(
    sampled_indices,
    nbins =
      n_individuals_60
  )
  
  
  # Reconstruct the fitting dataset once per replicate from pre-split clusters.
  bootstrap_data <- dplyr::bind_rows(
    individual_data_split_60[
      sampled_indices
    ]
  )
  
  
  # Repeated selections of the same eagle receive separate random-intercept IDs.
  bootstrap_data$bootstrap_individual_id <- factor(
    rep(
      seq_len(
        n_individuals_60
      ),
      times =
        individual_data_sizes_60[
          sampled_indices
        ]
    )
  )
  
  
  bootstrap_result <- tryCatch(
    {
      
      bootstrap_model <- mgcv::bam(
        formula =
          bootstrap_hfi_formula_60,
        
        data =
          bootstrap_data,
        
        family =
          stats::binomial(
            link = "logit"
          ),
        
        method =
          "fREML",
        
        discrete =
          TRUE,
        
        nthreads =
          1L
      )
      
      
      individual_effects <-
        calculate_fast_individual_contrasts_60(
          bootstrap_model
        )
      
      
      # ----------------------------------------------------------------------
      # Log-odds Q95-Q05 contrast
      #
      # This represents the model-scale difference between the two HFI levels:
      #
      # logit(P(Q95)) - logit(P(Q05))
      #
      # It is calculated for every bootstrap model.
      # ----------------------------------------------------------------------
      
      log_odds_contrast <-
        individual_effects$hfi_coefficient *
        (
          hfi_q95_z_60 -
            hfi_q05_z_60
        )
      
      
      sampled_absolute_differences <- rep(
        individual_effects$
          absolute_difference,
        
        times =
          sampled_individual_weights
      )
      
      
      sampled_absolute_differences_pp <- rep(
        individual_effects$
          absolute_difference_percentage_points,
        
        times =
          sampled_individual_weights
      )
      
      
      tibble::tibble(
        
        bootstrap_iteration =
          bootstrap_iteration,
        
        
        converged =
          bootstrap_model$converged,
        
        
        # HFI coefficient per bootstrap model
        hfi_coefficient =
          individual_effects$
          hfi_coefficient,
        
        
        # NEW: Q95-Q05 contrast on log-odds scale
        log_odds_contrast =
          log_odds_contrast,
        
        
        # Q95-Q05 contrast on probability scale
        mean_absolute_difference =
          stats::weighted.mean(
            individual_effects$
              absolute_difference,
            
            sampled_individual_weights
          ),
        
        
        mean_absolute_difference_percentage_points =
          stats::weighted.mean(
            individual_effects$
              absolute_difference_percentage_points,
            
            sampled_individual_weights
          ),
        
        
        mean_relative_difference_percent =
          stats::weighted.mean(
            individual_effects$
              relative_difference_percent,
            
            sampled_individual_weights
          ),
        
        
        median_individual_difference =
          stats::median(
            sampled_absolute_differences
          ),
        
        
        median_individual_difference_percentage_points =
          stats::median(
            sampled_absolute_differences_pp
          ),
        
        
        proportion_individuals_positive =
          stats::weighted.mean(
            individual_effects$
              absolute_difference > 0,
            
            sampled_individual_weights
          )
      )
    },
    
    
    error = function(error_message) {
      
      tibble::tibble(
        
        bootstrap_iteration =
          bootstrap_iteration,
        
        
        converged =
          FALSE,
        
        
        hfi_coefficient =
          NA_real_,
        
        
        # NEW: keep same structure if a bootstrap fails
        log_odds_contrast =
          NA_real_,
        
        
        mean_absolute_difference =
          NA_real_,
        
        
        mean_absolute_difference_percentage_points =
          NA_real_,
        
        
        mean_relative_difference_percent =
          NA_real_,
        
        
        median_individual_difference =
          NA_real_,
        
        
        median_individual_difference_percentage_points =
          NA_real_,
        
        
        proportion_individuals_positive =
          NA_real_
      )
    }
  )
  
  
  bootstrap_result
}


# 8.11 Run bootstrap replicates in parallel ----
individual_bootstrap_results_60 <- dplyr::bind_rows(
  parallel::mclapply(
    X =
      seq_len(
        n_individual_bootstrap_60
      ),
    
    FUN =
      run_fast_individual_bootstrap_60,
    
    mc.cores =
      n_bootstrap_cores_60,
    
    mc.preschedule =
      FALSE
  )
)

print(
  individual_bootstrap_results_60,
  n = 20
)


# 8.12 Retain successful bootstrap iterations ----
successful_individual_bootstrap_60 <-
  individual_bootstrap_results_60 %>%
  dplyr::filter(
    converged,
    dplyr::if_all(
      c(
        hfi_coefficient,
        mean_absolute_difference,
        log_odds_contrast,
        mean_absolute_difference_percentage_points,
        mean_relative_difference_percent,
        median_individual_difference,
        median_individual_difference_percentage_points,
        proportion_individuals_positive
      ),
      is.finite
    )
  )


# 8.13 Compile percentile bootstrap confidence intervals ----
bootstrap_hfi_summary_60 <-
  successful_individual_bootstrap_60 %>%
  dplyr::summarise(
    n_bootstrap_requested =
      n_individual_bootstrap_60,
    
    n_bootstrap_successful =
      dplyr::n(),
    
    hfi_coefficient_low =
      stats::quantile(
        hfi_coefficient,
        0.025,
        names = FALSE
      ),
    
    hfi_coefficient_high =
      stats::quantile(
        hfi_coefficient,
        0.975,
        names = FALSE
      ),
    
    mean_difference_low =
      stats::quantile(
        mean_absolute_difference,
        0.025,
        names = FALSE
      ),
    
    mean_difference_high =
      stats::quantile(
        mean_absolute_difference,
        0.975,
        names = FALSE
      ),
    
    mean_difference_percentage_points_low =
      stats::quantile(
        mean_absolute_difference_percentage_points,
        0.025,
        names = FALSE
      ),
    
    mean_difference_percentage_points_high =
      stats::quantile(
        mean_absolute_difference_percentage_points,
        0.975,
        names = FALSE
      ),
    
    mean_relative_difference_percent_low =
      stats::quantile(
        mean_relative_difference_percent,
        0.025,
        names = FALSE
      ),
    
    mean_relative_difference_percent_high =
      stats::quantile(
        mean_relative_difference_percent,
        0.975,
        names = FALSE
      ),
    
    median_individual_difference_low =
      stats::quantile(
        median_individual_difference,
        0.025,
        names = FALSE
      ),
    
    median_individual_difference_high =
      stats::quantile(
        median_individual_difference,
        0.975,
        names = FALSE
      ),
    
    median_individual_difference_percentage_points_low =
      stats::quantile(
        median_individual_difference_percentage_points,
        0.025,
        names = FALSE
      ),
    
    median_individual_difference_percentage_points_high =
      stats::quantile(
        median_individual_difference_percentage_points,
        0.975,
        names = FALSE
      ),
    
    proportion_individuals_positive_low =
      stats::quantile(
        proportion_individuals_positive,
        0.025,
        names = FALSE
      ),
    
    proportion_individuals_positive_high =
      stats::quantile(
        proportion_individuals_positive,
        0.975,
        names = FALSE
      ),
    
    log_odds_contrast_low =
      stats::quantile(
        log_odds_contrast,
        0.025,
        names = FALSE
      ),
    
    log_odds_contrast_high =
      stats::quantile(
        log_odds_contrast,
        0.975,
        names = FALSE
      ), 
    
    log_odds_contrast_mean =
      mean(
        log_odds_contrast
      ),
    
    proportion_bootstrap_mean_effect_positive =
      mean(
        mean_absolute_difference > 0
      )
  )

print(
  bootstrap_hfi_summary_60,
  n = Inf
)


# 8.14 Compile final report ----
# 8.14 Compile final report ----

final_bootstrap_report_60 <- tibble::tibble(
  
  measure = c(
    
    "HFI coefficient",
    
    "Q95-Q05 log-odds contrast",
    
    "Equal-weight mean Q95-Q05 difference",
    
    "Equal-weight mean Q95-Q05 difference (percentage points)",
    
    "Equal-weight relative Q95-Q05 difference (%)",
    
    "Median individual Q95-Q05 difference",
    
    "Median individual Q95-Q05 difference (percentage points)",
    
    "Proportion of individuals with positive contrast"
  ),
  
  
  estimate = c(
    
    equal_weight_hfi_summary_60$
      hfi_coefficient,
    
    bootstrap_hfi_summary_60$
      log_odds_contrast_mean,
    
    equal_weight_hfi_summary_60$
      mean_absolute_difference,
    
    equal_weight_hfi_summary_60$
      mean_absolute_difference_percentage_points,
    
    equal_weight_hfi_summary_60$
      mean_relative_difference_percent,
    
    equal_weight_hfi_summary_60$
      median_individual_difference,
    
    equal_weight_hfi_summary_60$
      median_individual_difference_percentage_points,
    
    equal_weight_hfi_summary_60$
      proportion_individuals_positive
  ),
  
  
  confidence_low = c(
    
    bootstrap_hfi_summary_60$
      hfi_coefficient_low,
    
    bootstrap_hfi_summary_60$
      log_odds_contrast_low,
    
    bootstrap_hfi_summary_60$
      mean_difference_low,
    
    bootstrap_hfi_summary_60$
      mean_difference_percentage_points_low,
    
    bootstrap_hfi_summary_60$
      mean_relative_difference_percent_low,
    
    bootstrap_hfi_summary_60$
      median_individual_difference_low,
    
    bootstrap_hfi_summary_60$
      median_individual_difference_percentage_points_low,
    
    bootstrap_hfi_summary_60$
      proportion_individuals_positive_low
  ),
  
  
  confidence_high = c(
    
    bootstrap_hfi_summary_60$
      hfi_coefficient_high,
    
    bootstrap_hfi_summary_60$
      log_odds_contrast_high,
    
    bootstrap_hfi_summary_60$
      mean_difference_high,
    
    bootstrap_hfi_summary_60$
      mean_difference_percentage_points_high,
    
    bootstrap_hfi_summary_60$
      mean_relative_difference_percent_high,
    
    bootstrap_hfi_summary_60$
      median_individual_difference_high,
    
    bootstrap_hfi_summary_60$
      median_individual_difference_percentage_points_high,
    
    bootstrap_hfi_summary_60$
      proportion_individuals_positive_high
  )
)


print(
  final_bootstrap_report_60,
  n = Inf
)

# measure (60 minutes dataset)                            estimate confidence_low confidence_high
# HFI coefficient                                            0.0351        -0.0149          0.0830
# Q95-Q05 log-odds contrast                                  0.115         -0.0445          0.248 
# Equal-weight mean Q95-Q05 difference                       0.0245        -0.0101          0.0589
# Equal-weight mean Q95-Q05 difference (percentage points)   2.45          -1.01            5.89  
# Equal-weight relative Q95-Q05 difference (%)               6.38          -2.65           15.7   
# Median individual Q95-Q05 difference                       0.0243        -0.0102          0.0581
# Median individual Q95-Q05 difference (percentage points)   2.43          -1.02            5.81  
# Proportion of individuals with positive contrast           1              0               1 

# measure (20 minutes dataset)                              estimate confidence_low confidence_high
# HFI coefficient                                            0.0856         0.0337           0.144
# Q95-Q05 log-odds contrast                                  0.264          0.0996           0.425
# Equal-weight mean Q95-Q05 difference                       0.0609         0.0241           0.101
# Equal-weight mean Q95-Q05 difference (percentage points)   6.09           2.41            10.1  
# Equal-weight relative Q95-Q05 difference (%)              12.3            4.80            20.8  
# Median individual Q95-Q05 difference                       0.0607         0.0240           0.101
# Median individual Q95-Q05 difference (percentage points)   6.07           2.40            10.1  
# Proportion of individuals with positive contrast           1              1                1 


#------------------------------------------------------------------------------
# STEP 9: test individual heterogeneity in HFI responses ----
#' **Steps:**
#' (i) separate within-individual HFI variation from between-individual exposure;
#' (ii) inspect whether each individual has sufficient HFI variation;
#' (iii) compare a common within-individual HFI slope with individual random slopes;
#' (iv) extract individual HFI slopes and Q05-Q95 log-odds contrasts;
#' (v) plot the distribution of individual responses.


# 9.1 Decompose HFI into within- and between-individual components ----
# hfi_between_z represents each individual's mean exposure.
# hfi_within_z represents deviations from that individual's mean exposure.
#
# hfi_within_z is not divided by the individual SD, so one unit remains equal
# to one population SD of HFI and slopes remain comparable among individuals.

hfi_random_slope_data_60 <- hfi_model_data_60 %>%
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
  dplyr::ungroup() %>%
  dplyr::mutate(
    individual_id =
      droplevels(
        factor(
          individual_id
        )
      )
  )


# 9.2 Inspect empirical HFI support within each individual ----
# This is an identifiability check, not an analysis of geographic context.
# An individual slope cannot be estimated reliably when HFI varies little
# within that individual.

global_q05_q95_range_z_60 <-
  hfi_q95_z_60 -
  hfi_q05_z_60

individual_hfi_slope_support_60 <- hfi_random_slope_data_60 %>%
  dplyr::group_by(
    individual_id
  ) %>%
  dplyr::summarise(
    n_transitions =
      dplyr::n(),
    
    n_bursts =
      dplyr::n_distinct(
        burst_id
      ),
    
    n_remain_aerial =
      sum(
        remain_aerial == 1L
      ),
    
    n_transition_terrestrial =
      sum(
        remain_aerial == 0L
      ),
    
    hfi_mean =
      mean(
        hfi_mean_1000m
      ),
    
    hfi_sd =
      stats::sd(
        hfi_mean_1000m
      ),
    
    hfi_min =
      min(
        hfi_mean_1000m
      ),
    
    hfi_q05 =
      as.numeric(
        stats::quantile(
          hfi_mean_1000m,
          0.05
        )
      ),
    
    hfi_q95 =
      as.numeric(
        stats::quantile(
          hfi_mean_1000m,
          0.95
        )
      ),
    
    hfi_max =
      max(
        hfi_mean_1000m
      ),
    
    hfi_within_sd_z =
      stats::sd(
        hfi_within_z
      ),
    
    hfi_within_range_z =
      max(
        hfi_within_z
      ) -
      min(
        hfi_within_z
      ),
    
    observed_range_relative_to_global_q05_q95 =
      hfi_within_range_z /
      global_q05_q95_range_z_60,
    
    .groups =
      "drop"
  ) %>%
  dplyr::arrange(
    hfi_within_range_z
  )

print(
  individual_hfi_slope_support_60,
  n = Inf
)

# Inspect:
# - hfi_within_sd_z and hfi_within_range_z: within-individual HFI information;
# - n_remain_aerial and n_transition_terrestrial: both outcomes are needed;
# - observed_range_relative_to_global_q05_q95:
#   values below 1 indicate that the individual did not observe the complete
#   population Q05-Q95 gradient.


# 9.3 Define models with a common slope and individual random slopes ----
# Both models contain identical fixed effects.
# The second model adds individual deviations from the population HFI slope.

common_hfi_slope_formula_60 <- remain_aerial ~
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

individual_hfi_slope_formula_60 <- remain_aerial ~
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


# 9.4 Fit both models using fast REML ----
common_hfi_slope_model_60 <- mgcv::bam(
  formula =
    common_hfi_slope_formula_60,
  
  data =
    hfi_random_slope_data_60,
  
  family =
    stats::binomial(
      link = "logit"
    ),
  
  method =
    "fREML",
  
  discrete =
    TRUE,
  
  nthreads =
    n_bootstrap_cores_60
)

individual_hfi_slope_model_60 <- mgcv::bam(
  formula =
    individual_hfi_slope_formula_60,
  
  data =
    hfi_random_slope_data_60,
  
  family =
    stats::binomial(
      link = "logit"
    ),
  
  method =
    "fREML",
  
  discrete =
    TRUE,
  
  nthreads =
    n_bootstrap_cores_60
)

common_hfi_slope_summary_60 <- summary(
  common_hfi_slope_model_60,
  re.test = TRUE
)

individual_hfi_slope_summary_60 <- summary(
  individual_hfi_slope_model_60,
  re.test = TRUE
)

print(
  individual_hfi_slope_summary_60
)


# 9.5 Compare common- and random-slope models ----
extract_hfi_slope_model_results_60 <- function(
    fitted_model,
    model_structure
) {
  
  model_summary <- summary(
    fitted_model,
    re.test = TRUE
  )
  
  coefficient_table <-
    model_summary$p.table
  
  tibble::tibble(
    model_structure =
      model_structure,
    
    n_observations =
      model_summary$n,
    
    AIC =
      stats::AIC(
        fitted_model
      ),
    
    total_EDF =
      sum(
        fitted_model$edf
      ),
    
    deviance_explained =
      model_summary$dev.expl,
    
    within_HFI_estimate =
      coefficient_table[
        "hfi_within_z",
        "Estimate"
      ],
    
    within_HFI_SE =
      coefficient_table[
        "hfi_within_z",
        "Std. Error"
      ],
    
    within_HFI_confidence_low =
      within_HFI_estimate -
      1.96 *
      within_HFI_SE,
    
    within_HFI_confidence_high =
      within_HFI_estimate +
      1.96 *
      within_HFI_SE,
    
    between_HFI_estimate =
      coefficient_table[
        "hfi_between_z",
        "Estimate"
      ],
    
    between_HFI_SE =
      coefficient_table[
        "hfi_between_z",
        "Std. Error"
      ],
    
    converged =
      fitted_model$converged
  )
}

individual_hfi_slope_model_comparison_60 <- dplyr::bind_rows(
  extract_hfi_slope_model_results_60(
    common_hfi_slope_model_60,
    "Common HFI slope"
  ),
  
  extract_hfi_slope_model_results_60(
    individual_hfi_slope_model_60,
    "Individual random HFI slopes"
  )
) %>%
  dplyr::arrange(
    AIC
  ) %>%
  dplyr::mutate(
    delta_AIC =
      AIC -
      min(AIC)
  )

print(
  individual_hfi_slope_model_comparison_60,
  n = Inf
)

# Compare:
# - AIC and delta_AIC;
# - within_HFI_estimate and its confidence interval;
# - whether within_HFI_SE increases after allowing individual slopes;
# - deviance_explained and convergence.


# 9.6 Extract the random-slope test ----
random_slope_term_60 <- grep(
  pattern =
    "hfi_within_z.*individual_id|individual_id.*hfi_within_z",
  
  x =
    rownames(
      individual_hfi_slope_summary_60$
        s.table
    ),
  
  value =
    TRUE
)

stopifnot(
  length(
    random_slope_term_60
  ) == 1L
)

random_slope_test_60 <- tibble::tibble(
  term =
    random_slope_term_60,
  
  edf =
    individual_hfi_slope_summary_60$
    s.table[
      random_slope_term_60,
      1
    ],
  
  reference_df =
    individual_hfi_slope_summary_60$
    s.table[
      random_slope_term_60,
      2
    ],
  
  statistic =
    individual_hfi_slope_summary_60$
    s.table[
      random_slope_term_60,
      3
    ],
  
  p_value =
    individual_hfi_slope_summary_60$
    s.table[
      random_slope_term_60,
      ncol(
        individual_hfi_slope_summary_60$
          s.table
      )
    ]
)

print(
  random_slope_test_60
)

# The random-effect p-value tests whether the random-slope variance is
# compatible with zero. It is an approximate boundary likelihood-ratio test.


# 9.7 Extract random-effect variance components ----
# 9.7 Extract random-effect variance components ----

random_slope_variance_components_raw_60 <-
  mgcv::gam.vcomp(
    individual_hfi_slope_model_60
  )


# gam.vcomp() can return either:
# - a matrix directly;
# - or a list containing the matrix in $vc.
random_slope_variance_components_values_60 <-
  if (
    is.list(
      random_slope_variance_components_raw_60
    ) &&
    !is.null(
      random_slope_variance_components_raw_60$vc
    )
  ) {
    random_slope_variance_components_raw_60$vc
  } else {
    random_slope_variance_components_raw_60
  }


# Convert the result into a consistent table.
random_slope_variance_components_60 <-
  if (
    is.null(
      dim(
        random_slope_variance_components_values_60
      )
    )
  ) {
    
    # Case where gam.vcomp() returns a named vector.
    tibble::tibble(
      term =
        names(
          random_slope_variance_components_values_60
        ),
      
      standard_deviation =
        as.numeric(
          random_slope_variance_components_values_60
        ),
      
      confidence_low =
        NA_real_,
      
      confidence_high =
        NA_real_
    )
    
  } else {
    
    # Case where gam.vcomp() returns a matrix.
    n_variance_components_60 <-
      nrow(
        random_slope_variance_components_values_60
      )
    
    variance_component_confidence_low_60 <-
      if (
        ncol(
          random_slope_variance_components_values_60
        ) >= 2L
      ) {
        as.numeric(
          random_slope_variance_components_values_60[
            ,
            2
          ]
        )
      } else {
        rep(
          NA_real_,
          n_variance_components_60
        )
      }
    
    variance_component_confidence_high_60 <-
      if (
        ncol(
          random_slope_variance_components_values_60
        ) >= 3L
      ) {
        as.numeric(
          random_slope_variance_components_values_60[
            ,
            3
          ]
        )
      } else {
        rep(
          NA_real_,
          n_variance_components_60
        )
      }
    
    tibble::tibble(
      term =
        rownames(
          random_slope_variance_components_values_60
        ),
      
      standard_deviation =
        as.numeric(
          random_slope_variance_components_values_60[
            ,
            1
          ]
        ),
      
      confidence_low =
        variance_component_confidence_low_60,
      
      confidence_high =
        variance_component_confidence_high_60
    )
  }


print(
  random_slope_variance_components_60,
  n = Inf
)

# 9.8 Extract individual HFI slopes using the model matrix ----
# For each individual, the linear predictor is evaluated once at
# hfi_within_z = 0 and once at hfi_within_z = 1.
#
# The difference isolates:
# population within-HFI slope + individual random-slope deviation.
#
# All environmental covariates and the random intercept cancel from this
# difference, so the resulting variation reflects the fitted HFI slopes.

individual_slope_template_60 <- hfi_random_slope_data_60 %>%
  dplyr::arrange(
    individual_id
  ) %>%
  dplyr::group_by(
    individual_id
  ) %>%
  dplyr::slice(1) %>%
  dplyr::ungroup()

individual_slope_template_60$individual_id <- factor(
  individual_slope_template_60$individual_id,
  levels =
    levels(
      hfi_random_slope_data_60$
        individual_id
    )
)

individual_slope_data_0_60 <-
  individual_slope_template_60

individual_slope_data_1_60 <-
  individual_slope_template_60

individual_slope_data_0_60$hfi_within_z <- 0
individual_slope_data_1_60$hfi_within_z <- 1

slope_matrix_0_60 <- stats::predict(
  individual_hfi_slope_model_60,
  newdata =
    individual_slope_data_0_60,
  type =
    "lpmatrix"
)

slope_matrix_1_60 <- stats::predict(
  individual_hfi_slope_model_60,
  newdata =
    individual_slope_data_1_60,
  type =
    "lpmatrix"
)

individual_slope_difference_matrix_60 <-
  slope_matrix_1_60 -
  slope_matrix_0_60

individual_slope_coefficients_60 <-
  stats::coef(
    individual_hfi_slope_model_60
  )

individual_slope_covariance_60 <- stats::vcov(
  individual_hfi_slope_model_60,
  unconditional =
    TRUE
)

individual_slope_estimate_60 <- as.numeric(
  individual_slope_difference_matrix_60 %*%
    individual_slope_coefficients_60
)

individual_slope_variance_60 <- rowSums(
  (
    individual_slope_difference_matrix_60 %*%
      individual_slope_covariance_60
  ) *
    individual_slope_difference_matrix_60
)

individual_slope_SE_60 <- sqrt(
  pmax(
    individual_slope_variance_60,
    0
  )
)


# 9.9 Calculate individual Q05-Q95 log-odds contrasts ----
individual_hfi_slopes_60 <- tibble::tibble(
  individual_id =
    individual_slope_template_60$
    individual_id,
  
  HFI_slope =
    individual_slope_estimate_60,
  
  HFI_slope_SE =
    individual_slope_SE_60,
  
  HFI_slope_confidence_low =
    HFI_slope -
    1.96 *
    HFI_slope_SE,
  
  HFI_slope_confidence_high =
    HFI_slope +
    1.96 *
    HFI_slope_SE,
  
  Q95_Q05_log_odds_contrast =
    HFI_slope *
    global_q05_q95_range_z_60,
  
  Q95_Q05_contrast_SE =
    HFI_slope_SE *
    abs(
      global_q05_q95_range_z_60
    ),
  
  Q95_Q05_confidence_low =
    Q95_Q05_log_odds_contrast -
    1.96 *
    Q95_Q05_contrast_SE,
  
  Q95_Q05_confidence_high =
    Q95_Q05_log_odds_contrast +
    1.96 *
    Q95_Q05_contrast_SE,
  
  Q95_Q05_odds_ratio =
    exp(
      Q95_Q05_log_odds_contrast
    ),
  
  positive_estimate =
    Q95_Q05_log_odds_contrast > 0,
  
  confidence_interval_entirely_positive =
    Q95_Q05_confidence_low > 0,
  
  confidence_interval_entirely_negative =
    Q95_Q05_confidence_high < 0
) %>%
  dplyr::left_join(
    individual_hfi_slope_support_60,
    by =
      "individual_id"
  ) %>%
  dplyr::arrange(
    Q95_Q05_log_odds_contrast
  )

print(
  individual_hfi_slopes_60,
  n = Inf
)


# 9.10 Summarise individual HFI-response heterogeneity ----
population_within_HFI_slope_60 <-
  stats::coef(
    individual_hfi_slope_model_60
  )[
    "hfi_within_z"
  ]

population_between_HFI_slope_60 <-
  stats::coef(
    individual_hfi_slope_model_60
  )[
    "hfi_between_z"
  ]

individual_hfi_slope_distribution_60 <- individual_hfi_slopes_60 %>%
  dplyr::summarise(
    n_individuals =
      dplyr::n(),
    
    population_within_HFI_slope =
      population_within_HFI_slope_60,
    
    population_between_HFI_slope =
      population_between_HFI_slope_60,
    
    population_Q95_Q05_log_odds_contrast =
      population_within_HFI_slope_60 *
      global_q05_q95_range_z_60,
    
    median_individual_HFI_slope =
      stats::median(
        HFI_slope
      ),
    
    individual_HFI_slope_Q25 =
      stats::quantile(
        HFI_slope,
        0.25
      ),
    
    individual_HFI_slope_Q75 =
      stats::quantile(
        HFI_slope,
        0.75
      ),
    
    minimum_individual_HFI_slope =
      min(
        HFI_slope
      ),
    
    maximum_individual_HFI_slope =
      max(
        HFI_slope
      ),
    
    median_Q95_Q05_log_odds_contrast =
      stats::median(
        Q95_Q05_log_odds_contrast
      ),
    
    proportion_positive_estimates =
      mean(
        positive_estimate
      ),
    
    proportion_confidently_positive =
      mean(
        confidence_interval_entirely_positive
      ),
    
    proportion_confidently_negative =
      mean(
        confidence_interval_entirely_negative
      )
  )

print(
  individual_hfi_slope_distribution_60,
  n = Inf
)


#------------------------------------------------------------------------------
# Visualisation n°1: individual Q05-Q95 HFI log-odds contrasts ----

individual_hfi_slope_plot_60 <- individual_hfi_slopes_60 %>%
  ggplot2::ggplot(
    ggplot2::aes(
      x =
        stats::reorder(
          individual_id,
          Q95_Q05_log_odds_contrast
        ),
      
      y =
        Q95_Q05_log_odds_contrast,
      
      ymin =
        Q95_Q05_confidence_low,
      
      ymax =
        Q95_Q05_confidence_high
    )
  ) +
  ggplot2::geom_hline(
    yintercept =
      0,
    linetype =
      "dashed"
  ) +
  ggplot2::geom_errorbar(
    width =
      0
  ) +
  ggplot2::geom_point(
    size =
      2
  ) +
  ggplot2::coord_flip() +
  ggplot2::labs(
    x =
      "Individual",
    
    y =
      "Individual Q95-Q05 HFI contrast (log-odds)",
    
    title =
      "Individual heterogeneity in the response to HFI",
    
    subtitle =
      "Points are population slope plus individual random-slope deviations"
  ) +
  ggplot2::theme_bw()

print(
  individual_hfi_slope_plot_60
)


#------------------------------------------------------------------------------
# Visualisation n°2: distribution of individual HFI slopes ----

individual_hfi_slope_distribution_plot_60 <- individual_hfi_slopes_60 %>%
  ggplot2::ggplot(
    ggplot2::aes(
      x =
        Q95_Q05_log_odds_contrast
    )
  ) +
  ggplot2::geom_histogram(
    bins =
      15
  ) +
  ggplot2::geom_vline(
    xintercept =
      0,
    linetype =
      "dashed"
  ) +
  ggplot2::geom_vline(
    xintercept =
      population_within_HFI_slope_60 *
      global_q05_q95_range_z_60,
    linewidth =
      0.8
  ) +
  ggplot2::labs(
    x =
      "Individual Q95-Q05 HFI contrast (log-odds)",
    
    y =
      "Number of individuals",
    
    title =
      "Distribution of individual HFI responses",
    
    subtitle =
      "Solid line: population within-individual HFI effect"
  ) +
  ggplot2::theme_bw()

print(
  individual_hfi_slope_distribution_plot_60
)


# evaluate model 
summary(
  individual_hfi_slope_model_60
)$p.table


#------------------------------------------------------------------------------
# 9.11 Identify individuals with negative HFI responses and empirical support ----
#' Identify whether negative individual slopes are supported by sufficient data.


negative_hfi_individuals_60 <- individual_hfi_slopes_60 %>%
  dplyr::filter(
    
    # Estimated negative response
    Q95_Q05_log_odds_contrast < 0
  ) %>%
  dplyr::select(
    
    individual_id,
    
    HFI_slope,
    
    Q95_Q05_log_odds_contrast,
    
    Q95_Q05_confidence_low,
    
    Q95_Q05_confidence_high,
    
    n_transitions,
    
    n_bursts,
    
    hfi_sd,
    
    hfi_within_range_z,
    
    observed_range_relative_to_global_q05_q95
  ) %>%
  dplyr::arrange(
    Q95_Q05_log_odds_contrast
  )


print(
  negative_hfi_individuals_60,
  n = Inf
)

individual_hfi_slopes_60 %>%
  dplyr::mutate(
    response_group =
      dplyr::if_else(
        Q95_Q05_log_odds_contrast < 0,
        "negative",
        "positive"
      )
  ) %>%
  dplyr::group_by(
    response_group
  ) %>%
  dplyr::summarise(
    n_individuals =
      dplyr::n(),
    
    median_transitions =
      median(
        n_transitions
      ),
    
    median_hfi_range =
      median(
        hfi_within_range_z
      ),
    
    median_hfi_sd =
      median(
        hfi_sd
      )
  )
