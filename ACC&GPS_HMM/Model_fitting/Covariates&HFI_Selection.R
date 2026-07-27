#' -----------------------------------------------------------------------------
# Title: Covariate and HFI variables selection ----
#' Authors : Louise Faure
#' Date : 20.07.26
#' 
#' Info : this script follow the Extract_covariates.R script where covariates are
#' extracted below each location and within two buffers. 
#' 
# Main steps:
#' (1) prepare the dataset:
#'     (i) calculate elapsed duration in the current behavioural state;
#'     (ii) construct the transition matrix and retain aerial-origin transitions;
#'     (iii) remove incomplete observations and document individual exclusions
#'     (iv) attribute a weight to individuals based on their nbr of transitions.
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


#------------------------------------------------------------------------------ STEP 1: prepare dataset ----
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

# Summary of the final dataset
final_dataset_summary_60 <- aerial_transitions_60 %>%
  dplyr::summarise(n_observations = dplyr::n(),
    n_individuals =dplyr::n_distinct(individual_id),
    n_bursts =dplyr::n_distinct(burst_id),
    n_remain_aerial = sum(remain_aerial == 1L),
    n_transition_terrestrial = sum( remain_aerial == 0L),
    proportion_remain_aerial = mean(remain_aerial))

print(final_dataset_summary_60)
# n_observations n_individuals n_bursts n_remain_aerial n_transition_terrestrial proportion_remain_aerial
# 10872            62            4875            4244                     6628                  0.390





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
  dplyr::filter( !is.na(individual_id),
    dplyr::if_all( dplyr::all_of(backbone_variables_60), is.finite)) %>%
  dplyr::mutate(individual_id = droplevels(factor(individual_id)),
    aerial_duration_h = aerial_duration_min / 60)


# 2.2 Centre and standardize age and elapsed aerial duration ----
age_center_60 <- mean(backbone_data_60$age_since_emig_weeks)
age_scale_60 <- stats::sd(backbone_data_60$age_since_emig_weeks)
duration_center_60 <- mean(backbone_data_60$aerial_duration_h)
duration_scale_60 <- stats::sd(backbone_data_60$aerial_duration_h)

backbone_data_60 <- backbone_data_60 %>%
  dplyr::mutate(
    age_z = (age_since_emig_weeks - age_center_60) / age_scale_60,
    duration_z = (aerial_duration_h - duration_center_60) / duration_scale_60,
    age_z2 = age_z^2,
    duration_z2 = duration_z^2)

backbone_scaling_60 <- tibble::tibble(
  variable = c("age_since_emig_weeks","aerial_duration_h"),
  center = c(age_center_60,duration_center_60),
  scale = c(age_scale_60,duration_scale_60))

print(backbone_scaling_60)


# 2.3 Define candidate functional forms ----
# Quadratic forms always contain both the linear and squared terms.
# Shrinkage cubic splines ("cs") can be penalized towards a simple or
# effectively absent relationship.
age_terms_60 <- list(linear ="age_z",
  quadratic = c("age_z","age_z2"),
  smooth = paste0("s(age_z, bs = 'cs', k = ",age_k,")"))
duration_terms_60 <- list(
  linear ="duration_z",
  quadratic = c("duration_z","duration_z2"),
  smooth = paste0("s(duration_z, bs = 'cs', k = ",duration_k,")"))


# 2.4 Define the candidate backbone model set ----
backbone_model_metadata_60 <- tidyr::expand_grid(
  age_form = c("linear","quadratic","smooth"),
  duration_form = c("linear","quadratic","smooth")) %>%
  dplyr::mutate( model = paste("age",age_form,"duration",duration_form,sep = "_"),
    
    # This ranking implements the a priori preference for simpler structures.
    age_complexity = dplyr::recode(
      age_form,
      linear = 1L,
      quadratic = 2L,
      smooth = 3L),
    
    duration_complexity = dplyr::recode(
      duration_form,
      linear = 1L,
      quadratic = 2L,
      smooth = 3L),
    
    complexity_score =
      age_complexity +
      duration_complexity)


# 2.5 Build candidate model formulas ----
# cos_diel and sin_time are retained together in every model.
# s(individual_id, bs = "re") represents the individual random intercept.
backbone_formulas_60 <- stats::setNames(
  lapply(
    seq_len(
      nrow(backbone_model_metadata_60)),
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
  backbone_model_metadata_60$model)


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
  names(backbone_formulas_60))


# 2.7 Compare candidate backbone models ----
backbone_model_comparison_60 <- dplyr::bind_rows(
  lapply(names(backbone_models_60),
    function(model_name) {
      
      fitted_model <-
        backbone_models_60[[model_name]]
      
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
        exp(-0.5 * delta_AIC)),
    
    similar_support =
      delta_AIC <=
      delta_aic_threshold)

print(backbone_model_comparison_60,n = Inf)


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

provisional_backbone_name_60 <- supported_backbone_models_60 %>%
  dplyr::slice(1) %>%
  dplyr::pull(model)

provisional_backbone_information_60 <-
  backbone_model_metadata_60 %>%
  dplyr::filter(
    model == provisional_backbone_name_60)

selected_age_form_60 <- provisional_backbone_information_60$age_form

selected_duration_form_60 <- provisional_backbone_information_60$duration_form

provisional_backbone_model_ml_60 <-backbone_models_60[[provisional_backbone_name_60]]

print(provisional_backbone_name_60) # "age_linear_duration_quadratic"
print(summary(provisional_backbone_model_ml_60))
#' Formula: remain_aerial ~ cos_diel + sin_time + age_z + duration_z + duration_z2 + 
#' s(individual_id, bs = "re")




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
      environmental_covariates_60))

correlation_matrix_60 <- stats::cor(correlation_data_60,method = "pearson",use = "complete.obs")

#------------------------------------------------------------------------------- VISUALISATION n°1 : Pearson correlation matrix ----
corrplot::corrplot(correlation_matrix_60,
  method = "color",
  type = "lower",
  order = "original",
  col = corrplot::COL2("RdBu",200),
  col.lim = c(-1,1),
  diag = FALSE,
  addCoef.col = "black",
  number.digits = 2,
  number.cex = 0.6,
  tl.col = "black",
  tl.cex = 0.6,
  tl.srt = 45,
  cl.cex = 0.8,
  addgrid.col = "white",
  title ="Pearson correlations among candidate environmental confounders",
  mar = c(0,0,3,0))
#------------------------------------------------------------------------------- (end) Visualisation n°1 : Pearson correlation matrix

#' we remove slope due to its high correlation with ruggedness and lower 
#' explanation of the topographic complexity 
#' we remove as well proportion of rocky terrain for its association with 
#' elevation and limited biological interpretability











# ------------------------------------------------------------------------------ STEP 4 : identify several HFI metrics ----
#' **Steps:**
#' (i) evaluate the support per hfi class (range of HFI used by individuals per
#' selected metric of HFI and number of individuals in the extreme range of HFI)
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
      hfi_candidates_60)) %>%
  tidyr::pivot_longer(
    cols =dplyr::all_of(hfi_candidates_60),
    names_to ="hfi_variable",
    values_to ="hfi_value") %>%
  dplyr::mutate(hfi_variable = factor(hfi_variable,
      levels =hfi_candidates_60))

# 4.1.2 Compile the empirical range of each HFI metric ----
hfi_thresholds_60 <- hfi_support_long_60 %>%
  dplyr::group_by(hfi_variable) %>%
  dplyr::summarise(hfi_min =min(hfi_value),
    hfi_q05 =as.numeric(stats::quantile(hfi_value,probs = 0.05)),
    hfi_median =stats::median(hfi_value),
    hfi_q95 =as.numeric(stats::quantile(hfi_value,probs = 0.95)),
    hfi_max =max(hfi_value),.groups ="drop")

print(hfi_thresholds_60,n = Inf)

# 4.1.3 Count observations, individuals and bursts across HFI classes ----
hfi_class_support_60 <- hfi_support_long_60 %>%
  dplyr::mutate(
    hfi_class = cut(
      hfi_value,
      breaks =hfi_classes_60,
      labels =hfi_labels_60,
      include.lowest =TRUE,
      right =TRUE)) %>%
  dplyr::group_by(hfi_variable,hfi_class) %>%
  dplyr::summarise(n_observations =dplyr::n(),
    n_individuals =dplyr::n_distinct(individual_id),
    n_bursts =dplyr::n_distinct(burst_id),
    n_remain_aerial =sum(remain_aerial == 1L),
    n_transition_terrestrial =sum(remain_aerial == 0L),
    proportion_remain_aerial = mean(remain_aerial),
    .groups ="drop") %>%
  dplyr::arrange(hfi_variable,hfi_class)

# 4.1.4 Prepare individual and transition counts for one plot ----
hfi_class_plot_data_60 <- hfi_class_support_60 %>%
  dplyr::select(hfi_variable,hfi_class,n_individuals,n_observations) %>%
  tidyr::pivot_longer(cols = c(n_individuals,n_observations),
    names_to ="support_measure",
    values_to ="count") %>%
  dplyr::mutate(
    support_measure = dplyr::recode(
      support_measure,
      n_individuals =
        "Number of individuals",
      n_observations =
        "Number of aerial-origin transitions"))

#------------------------------------------------------------------------------- VISUALISATION n°2 : empirical support ----
hfi_support_comparison_plot_60 <- hfi_class_plot_data_60 %>%
  ggplot2::ggplot(ggplot2::aes(x =hfi_class,y =count,group =hfi_variable,color =hfi_variable)) +
  ggplot2::geom_line(linewidth = 0.9) +
  ggplot2::geom_point(size = 2.5) +
  ggplot2::facet_wrap(~ support_measure,ncol = 1,scales = "free_y") +
  ggplot2::labs(x ="Absolute HFI class",y ="Count",color ="HFI metric", title = "Empirical support across HFI classes",
    subtitle ="Comparison of individual and transition support among candidate HFI metrics") +
  ggplot2::theme_bw() +
  ggplot2::theme(axis.text.x =ggplot2::element_text(angle = 45,hjust = 1),
    legend.position ="right")

print(hfi_support_comparison_plot_60)

# Summarise support in the lowest and highest HFI classes ----
hfi_extreme_support_60 <- hfi_class_support_60 %>%
  dplyr::filter(
    as.character(hfi_class) %in% c("0-0.10", ">0.80")
  ) %>%
  dplyr::mutate(
    hfi_group = dplyr::recode(
      as.character(hfi_class),
      "0-0.10" = "0_10",
      ">0.80"  = "above_0_80"
    )
  ) %>%
  dplyr::select(
    hfi_variable,
    hfi_group,
    n_individuals,
    n_bursts
  ) %>%
  tidyr::pivot_wider(
    names_from = hfi_group,
    values_from = c(n_individuals, n_bursts),
    names_glue = "{.value}_HFI_{hfi_group}"
  ) %>%
  dplyr::left_join(
    hfi_thresholds_60 %>%
      dplyr::select(
        hfi_variable,
        hfi_min,
        hfi_max
      ),
    by = "hfi_variable"
  ) %>%
  dplyr::select(
    hfi_variable,
    hfi_min,
    hfi_max,
    n_individuals_HFI_0_10,
    n_individuals_HFI_above_0_80,
    n_bursts_HFI_0_10,
    n_bursts_HFI_above_0_80
  )

print(hfi_extreme_support_60,n = Inf)
#------------------------------------------------------------------------------- (end) Visualisation n°2 : empirical support
selected_hfi_candidates_60 <- c("hfi_mean_1000m","hfi_q90_1000m","hfi_q75_500m")


# 4.2 Compile VIF for HFI candidate models ----
#' **Steps:**
#' (i) calculate VIF for each retained HFI metric within each environmental model;
#' (ii) compile a detailed table and a summary table;
#' (iii) plot both tables as images.

# 4.2.1 Define the candidate environmental model structures ----
vif_model_structures_60 <- list(
  altitudinal_gradient = c("dem_elevation"),
  topographic_complexity = c("ruggedness_100m","distance_to_ridgeline_100m", "dem_elevation"),
  habitat_refuge = c("prop_forest_5cells","ruggedness_100m","dem_elevation"),
  open_habitat = c("prop_low_vegetation_5cells","dem_elevation"))

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
        
        other_predictors <- setdiff(predictors,variable_name)
        
        auxiliary_model <- stats::lm(
          formula = stats::reformulate(termlabels = other_predictors,response = variable_name),
          data = data)
        
        tibble::tibble(
          variable = variable_name,
          vif = 1 / (1 -summary(auxiliary_model)$r.squared))}))}

# 4.2.3 Calculate detailed VIF values for each HFI metric x model structure ----
vif_detailed_results_60 <- dplyr::bind_rows(
  lapply(
    selected_hfi_candidates_60,
    function(hfi_name) {
      
      dplyr::bind_rows(
        lapply(names(vif_model_structures_60), function(model_name) {
            
            model_predictors <- c(hfi_name,vif_model_structures_60[[model_name]])
            
            calculate_vif_60(
              data = aerial_transitions_60,
              predictors = model_predictors) %>%
              dplyr::mutate(
                hfi_variable = hfi_name,
                model_structure = model_name,
                predictor_role = dplyr::if_else(
                  variable == hfi_name,
                  "HFI",
                  "control"))}))})) %>%
  dplyr::select(hfi_variable,model_structure,variable,predictor_role,vif) %>%
  dplyr::arrange(hfi_variable,model_structure,dplyr::desc(vif))

# 4.2.4 Summarise VIF at the model level ----
vif_summary_results_60 <- vif_detailed_results_60 %>%
  dplyr::group_by(
    hfi_variable,
    model_structure) %>%
  dplyr::summarise(
    hfi_vif = vif[predictor_role == "HFI"][1],
    maximum_vif = max(vif),
    variable_with_maximum_vif = variable[which.max(vif)],.groups = "drop") %>%
  dplyr::arrange(
    hfi_variable,
    model_structure)


# 4.2.5 Define readable labels for image plots ----
vif_variable_labels_60 <- c(hfi_mean_1000m = "HFI mean\n1000 m",hfi_q75_1000m = "HFI Q75\n1000 m",
  hfi_q75_500m = "HFI Q75\n500 m",dem_elevation = "Elevation",ruggedness_100m = "Ruggedness",
  distance_to_ridgeline_100m = "Distance to\nridgeline",prop_forest_5cells = "Forest\nproportion",
  prop_low_vegetation_5cells = "Low vegetation\nproportion")

vif_structure_labels_60 <- c(altitudinal_gradient = "Altitudinal\ngradient",
  topographic_complexity = "Topographic\ncomplexity",habitat_refuge = "Habitat\nrefuge",
  open_habitat = "Open\nhabitat")

vif_hfi_labels_60 <- c(hfi_mean_1000m = "HFI mean 1000 m",hfi_q75_1000m = "HFI Q75 1000 m",hfi_q75_500m = "HFI Q75 500 m")

# 4.2.6 Prepare detailed VIF table for plotting ----
vif_detailed_plot_data_60 <- vif_detailed_results_60 %>%
  dplyr::mutate(
    variable_label = unname(
      vif_variable_labels_60[
        variable]),
    model_label = paste(
      unname(
        vif_hfi_labels_60[
          hfi_variable]),
      unname(
        vif_structure_labels_60[
          model_structure]),
      sep = "\n" ),
    
    variable_label = factor(
      variable_label,
      levels = rev(unique(unname(vif_variable_labels_60)))),
    
    model_label = factor(model_label,levels = unique(model_label)))


#------------------------------------------------------------------------------- VISUALISATION n°3: the detailed VIF table ----
vif_detailed_table_plot_60 <- ggplot2::ggplot(
  vif_detailed_plot_data_60,
  ggplot2::aes(x = model_label,y = variable_label,fill = vif)) +
  ggplot2::geom_tile(color = "white") +
  ggplot2::geom_text(
    ggplot2::aes(label = sprintf("%.2f",vif)),size = 3.2) +
  ggplot2::scale_fill_gradient(low = "white",high = "darkseagreen3") +
  ggplot2::labs(x = NULL,y = NULL,
    fill = "VIF",
    title = "Variance inflation factors across candidate HFI models",
    subtitle = "Detailed VIF values for each HFI metric and environmental model structure") +
  ggplot2::theme_bw() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45,hjust = 1,vjust = 1),panel.grid = ggplot2::element_blank())

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
environmental_model_covariates_60 <- c("dem_elevation","ruggedness_100m","distance_to_ridgeline_100m","prop_forest_5cells","prop_low_vegetation_5cells")
variables_to_standardize_60 <- c(selected_hfi_candidates_60,environmental_model_covariates_60)

# 5.1.2 Centre and standardize HFI and environmental covariates ----
# scale() subtracts the empirical mean and divides by the empirical SD.
hfi_model_data_60 <- backbone_data_60 %>%
  dplyr::mutate(
    dplyr::across(
      dplyr::all_of(
        variables_to_standardize_60),
      ~ as.numeric(
        scale(.x)),
      .names = "{.col}_z"))

standardized_model_covariates_60 <- paste0(variables_to_standardize_60,"_z")

# 5.1.3 Compile raw and standardized HFI thresholds ----
hfi_model_thresholds_60 <- dplyr::bind_rows(
  lapply(
    selected_hfi_candidates_60,
    function(hfi_name) {
      
      hfi_z_name <- paste0(
        hfi_name,
        "_z")
      
      tibble::tibble(
        hfi_variable = hfi_name,
        
        hfi_q05 = as.numeric(
          stats::quantile(
            hfi_model_data_60[[
              hfi_name
            ]],
            probs = 0.05)),
        
        hfi_q95 = as.numeric(
          stats::quantile(
            hfi_model_data_60[[
              hfi_name
            ]],
            probs = 0.95)),
        
        hfi_q05_z = as.numeric(
          stats::quantile(
            hfi_model_data_60[[
              hfi_z_name
            ]],
            probs = 0.05)),
        
        hfi_q95_z = as.numeric(
          stats::quantile(
            hfi_model_data_60[[
              hfi_z_name
            ]],
            probs = 0.95)))}))

print(hfi_model_thresholds_60,n = Inf)

# 5.2 Define and fit the 12 candidate models ----
# 5.2.1 Define environmental model structures ----
# All terms refer to centred and standardized covariates.
hfi_model_structures_60 <- list(
  elevation = c("dem_elevation_z"),
  topographic_complexity = c("dem_elevation_z","ruggedness_100m_z","distance_to_ridgeline_100m_z"),
  habitat_refuge = c("dem_elevation_z","prop_forest_5cells_z","ruggedness_100m_z"),
  open_habitat = c("dem_elevation_z","prop_low_vegetation_5cells_z"))

# 5.2.2 Define the 12 model specifications ----
hfi_model_metadata_60 <- tidyr::expand_grid(hfi_variable =selected_hfi_candidates_60,
  model_type = names(hfi_model_structures_60)) %>%
  dplyr::mutate(model = paste(hfi_variable,model_type,sep = "__"))

# 5.2.3 Add each HFI and environmental structure to the backbone ----
selected_backbone_formula_60 <- stats::formula(provisional_backbone_model_ml_60)

hfi_candidate_formulas_60 <- stats::setNames(
  lapply(
    seq_len(
      nrow(
        hfi_model_metadata_60)),
    function(i) {
      
      model_information <-
        hfi_model_metadata_60[i, ]
      
      hfi_term <- paste0(
        model_information$hfi_variable,
        "_z")
      
      additional_terms <- c(
        hfi_term,
        hfi_model_structures_60[[model_information$model_type]])
      
      stats::update.formula(
        selected_backbone_formula_60,
        
        stats::as.formula(
          paste(
            ". ~ . +",
            paste(
              additional_terms,
              collapse = " + ")))
      )}),
  hfi_model_metadata_60$model)

# 5.2.4 Fit all candidate models using ML ----
# The same observations and estimation method are used for every model.
hfi_candidate_models_ml_60 <- stats::setNames(
  lapply(
    names(hfi_candidate_formulas_60),
    function(model_name) {
      
      mgcv::gam(
        formula =hfi_candidate_formulas_60[[model_name]],
        data =hfi_model_data_60,
        weights =individual_weight,
        family =stats::binomial(link = "logit"),
        method =gam_selection_method)}),
  names(hfi_candidate_formulas_60))

# 5.3 Prepare observed covariate combinations for Q05-Q95 predictions ----
# Predictions are averaged first within each individual and then across
# individuals. Consequently, each individual contributes one value to the
# population-level summaries, independently of its number of transitions.
prediction_reference_data_60 <- hfi_model_data_60

# 5.4 Extract fit statistics and Q05-Q95 contrasts ----
# 5.4.1 Calculate model results ----
hfi_model_results_60 <- dplyr::bind_rows(
  lapply( names(hfi_candidate_models_ml_60),
    function(model_name) {
      
      fitted_model <- hfi_candidate_models_ml_60[[model_name]]
      
      model_information <-hfi_model_metadata_60 %>%dplyr::filter(
          model ==model_name)
      
      hfi_name <-model_information$hfi_variable
      
      hfi_term <- paste0(hfi_name,"_z")
      
      hfi_thresholds <-hfi_model_thresholds_60 %>%
        dplyr::filter(hfi_variable ==hfi_name)
      
      prediction_data_q05 <-prediction_reference_data_60
      
      prediction_data_q95 <-prediction_reference_data_60
      
      prediction_data_q05[[hfi_term]] <- hfi_thresholds$hfi_q05_z
      
      prediction_data_q95[[hfi_term]] <- hfi_thresholds$hfi_q95_z
      
      # Population-level linear predictors:
      # the individual random intercept is fixed at zero.
      linear_predictor_q05 <- stats::predict(
        fitted_model,
        newdata =prediction_data_q05,
        type ="link",
        exclude ="s(individual_id)")
      
      linear_predictor_q95 <- stats::predict(
        fitted_model,
        newdata =prediction_data_q95,
        type ="link",
        exclude ="s(individual_id)")
      
      probability_q05 <- stats::predict(
        fitted_model,
        newdata =prediction_data_q05,
        type ="response",
        exclude ="s(individual_id)")
      
      probability_q95 <- stats::predict(
        fitted_model,
        newdata =prediction_data_q95,
        type ="response",
        exclude ="s(individual_id)")
      
      coefficient_table <-
        summary(fitted_model)$p.table
      
      hfi_estimate <-coefficient_table[hfi_term,1]
      
      hfi_standard_error <-coefficient_table[hfi_term,2]
      
      hfi_standardized_range <-hfi_thresholds$hfi_q95_z - hfi_thresholds$hfi_q05_z
      
      log_odds_difference <-hfi_estimate *hfi_standardized_range
      
      log_odds_difference_se <-abs( hfi_standardized_range) *hfi_standard_error
      
      # Calculate predictions separately for each individual.
      individual_prediction_summary_60 <- tibble::tibble(
        individual_id =prediction_reference_data_60$individual_id,
        linear_predictor_q05 = as.numeric(linear_predictor_q05),
        linear_predictor_q95 =as.numeric(linear_predictor_q95),
        probability_q05 =as.numeric(probability_q05),
        probability_q95 =as.numeric(probability_q95)) %>%
        dplyr::group_by(individual_id) %>%
        dplyr::summarise(
          mean_linear_predictor_q05 =mean(linear_predictor_q05),
          mean_linear_predictor_q95 =mean(linear_predictor_q95),
          mean_probability_q05 =mean(probability_q05),
          mean_probability_q95 = mean(probability_q95),
          .groups = "drop") %>%
        dplyr::mutate(
          log_odds_difference_q95_q05 = mean_linear_predictor_q95 -mean_linear_predictor_q05,
          absolute_probability_difference = mean_probability_q95 -mean_probability_q05,
          relative_probability_difference_percent = 100 *absolute_probability_difference /mean_probability_q05)
      
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
        
        hfi_coefficient =
          hfi_estimate,
        hfi_coefficient_se =
          hfi_standard_error,
        
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
        )))

#------------------------------------------------------------------------------ VISUALISATION n°4: summary table of models ----
hfi_result_labels_60 <- c(hfi_mean_1000m ="HFI mean 1000 m",
  hfi_q90_1000m ="HFI Q90 1000 m",
  hfi_q75_500m ="HFI Q75 500 m")

model_type_labels_60 <- c(elevation ="Elevation",
  topographic_complexity ="Topographic complexity",
  habitat_refuge ="Refuge",
  open_habitat ="Open habitat")

# Create the formatted results table
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
    
    `HFI coefficient\n(log-odds per SD)` =
      sprintf(
        "%.3f",
        hfi_coefficient
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

# Reshape the formatted table for plotting
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

table_column_order_60 <- names(hfi_model_summary_table_60)

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


# plot table
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
        15))

print(hfi_model_summary_table_plot_60)

# 4.1.4 Summarise support at the 10th and 80th HFI quantiles ----

hfi_q10_q80_support_60 <- hfi_support_long_60 %>%
  dplyr::filter(
    !is.na(hfi_value),
    !is.na(individual_id),
    !is.na(burst_id)
  ) %>%
  dplyr::group_by(hfi_variable) %>%
  dplyr::mutate(
    HFI_q10 = as.numeric(
      stats::quantile(hfi_value, probs = 0.10, na.rm = TRUE)
    ),
    HFI_q80 = as.numeric(
      stats::quantile(hfi_value, probs = 0.80, na.rm = TRUE)
    )
  ) %>%
  dplyr::summarise(
    HFI_q10 = dplyr::first(HFI_q10),
    HFI_q80 = dplyr::first(HFI_q80),
    
    n_individuals_HFI_q10 =
      dplyr::n_distinct(individual_id[hfi_value <= HFI_q10]),
    
    n_individuals_HFI_q80 =
      dplyr::n_distinct(individual_id[hfi_value >= HFI_q80]),
    
    n_bursts_HFI_q10 =
      dplyr::n_distinct(burst_id[hfi_value <= HFI_q10]),
    
    n_bursts_HFI_q80 =
      dplyr::n_distinct(burst_id[hfi_value >= HFI_q80]),
    
    .groups = "drop"
  )

print(hfi_q10_q80_support_60,n = Inf)
#------------------------------------------------------------------------------- (end visualisation n°4)
#' we selected the Open Habitat HFI mean 1000m because the HFI coefficient is stable in comparison to 
#' the model that only include elevation. The AIC is also the smallest compared to the other models. 

# Print the summary of the selected model:
# HFI mean within 1000 m + open-habitat covariates
selected_model_name_60 <-
  "hfi_mean_1000m__open_habitat"

selected_open_habitat_model_60 <-hfi_candidate_models_ml_60[[selected_model_name_60]]

print(summary(selected_open_habitat_model_60))

