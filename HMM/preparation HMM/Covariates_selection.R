#' -----------------------------------------------------------------------------
# Title: Selection of model variables ----
#' Authors : Louise Faure
#' Date : 16.07.26
#' 
#' Info : this script follow the Extract_covariates.R script where covariates are
#' extracted below each location and within two buffers. 
#' 
#' Purpose : 
#' (1) adjust the basic model and inspect individual variation
#' (2) select variable forms using AIC;
#' (3) for highly correlated variables (|r| > 0.5), retain the variable
#'     associated with the lowest AIC;
#' (4) control the biological interpretation of the variables;
#' (5) conduct forward and backward model selection;
#' (6) diagnose the selected model.
#' 
#' Ref. the statistical selection of the variables and model is based on Togunov
#' et al., Mov Ecol, 2022. 
#' -----------------------------------------------------------------------------


# library 
library(dplyr)
library(tidyr)
library(tibble)
library(lme4)

# data 
GE_60_min_covariates_hfi <- readRDS("/Users/louisefaure/Desktop/dossier sans titre/donnees filtree/GE_60_min_covariates_hfi.rds")

# parameters
state_levels <- c("aerial", "terrestrial")
transition_levels <- c("aerial_to_aerial", "aerial_to_terrestrial", "terrestrial_to_aerial", "terrestrial_to_terrestrial")
glmer_control <- lme4::glmerControl(
  optimizer = "bobyqa",
  optCtrl = list(
    maxfun = 200000))



#'------------------------------------------------------------------------------
# STEP 1: adjust basic model and inspect individual variation ----
#'
#' (1) prepare the transition dataset and transition matrix;
#' (2) fit the null transition model;
#' (3) fit a model with an individual random intercept;
#' (4) inspect AIC between two models.


# 1.1 Prepare transitions between consecutive locations ----
# A transition is created only between consecutive locations belonging
# to the same individual and the same burst_n.
transitions_60 <- GE_60_min_covariates_hfi %>%
  arrange(
    individual.local.identifier,
    burst_n,
    timestamp
  ) %>%
  mutate(
    behavior_binary_clean = tolower(
      trimws(
        as.character(behavior_binary)
      )
    ),
    
    behavior_binary_clean = dplyr::recode(
      behavior_binary_clean,
      "aerian" = "aerial",
      "flight" = "aerial"
    ),
    
    behavior_from = factor(
      behavior_binary_clean,
      levels = state_levels
    ),
    
    landcover_from = case_when(
      landcover_100m %in% c(2, 3, 4) ~ "forest",
      landcover_100m %in% c(5, 6, 7) ~ "low_vegetation",
      landcover_100m %in% c(8, 9, 11) ~ "rocky_terrain",
      TRUE ~ NA_character_
    )
  ) %>%
  group_by(
    individual.local.identifier,
    burst_n
  ) %>%
  mutate(
    behavior_to = dplyr::lead(
      behavior_from
    ),
    
    landcover_to = dplyr::lead(
      landcover_from
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
  ungroup() %>%
  filter(
    !is.na(behavior_from),
    !is.na(behavior_to)
  ) %>%
  mutate(
    landcover_from = factor(
      landcover_from,
      levels = c(
        "rocky_terrain",
        "low_vegetation",
        "forest"
      )
    ),
    
    landcover_to = factor(
      landcover_to,
      levels = c(
        "rocky_terrain",
        "low_vegetation",
        "forest"
      )
    ),
    
    individual_id = factor(
      individual.local.identifier
    ),
    
    next_terrestrial = as.integer(
      behavior_to == "terrestrial"
    ),
    
    transition_type = paste(
      behavior_from,
      behavior_to,
      sep = "_to_"
    )
  )


# 1.2 Nbr of transitions & bursts per individual and thinning ----
transitions_by_individual <- transitions_60 %>%
  group_by(
    individual.local.identifier
  ) %>%
  summarise(
    n_transitions = n(),
    n_bursts = n_distinct(burst_n),
    .groups = "drop"
  ) %>%
  arrange(
    n_transitions
  )

print(transitions_by_individual, n = Inf)
# Langgries21 (eobs 7586) and Almen18 (eobs 5861) have respectively
# 3 and 34 transitions and are therefore removed.

individuals_to_remove <- c("Langgries21 (eobs 7586)","Almen18 (eobs 5861)")
transitions_60 <- transitions_60 %>% filter(!individual.local.identifier %in% individuals_to_remove)

# 1.3 Empirical transition-count matrix ----
transition_count_matrix_60 <- xtabs(
  ~ behavior_from + behavior_to,
  data = transitions_60)

print(transition_count_matrix_60)
#                  behavior_to
# behavior_from aerial terrestrial
# aerial        4259        6652
# terrestrial   6518       31257

# 1.4 Empirical transition-probability matrix ----
transition_probability_matrix_60 <- prop.table(
  transition_count_matrix_60,
  margin = 1)

print(round(transition_probability_matrix_60,digits = 4))
#                 behavior_to
# behavior_from aerial terrestrial
# aerial      0.3903      0.6097
# terrestrial 0.1725      0.8275

# 1.5 Fit the null Markov transition model ----
# With aerial as the reference:
model_null_60 <- stats::glm(
  next_terrestrial ~ behavior_from,
  data = transitions_60,
  family = stats::binomial(
    link = "logit"))

print(summary(model_null_60))

# 1.6 Transition probabilities estimated by the null model ----
null_prediction_data <- data.frame(
  behavior_from = factor(
    state_levels,
    levels = state_levels))

p_next_terrestrial_null <- stats::predict(
  model_null_60,
  newdata = null_prediction_data,
  type = "response")

null_transition_matrix_60 <- cbind(
  aerial = 1 - p_next_terrestrial_null,
  terrestrial = p_next_terrestrial_null)

rownames(null_transition_matrix_60) <- state_levels

print(round(null_transition_matrix_60,digits = 4))
#             aerial terrestrial
# aerial      0.3903      0.6097
# terrestrial 0.1725      0.8275

# 1.7 Fit a model with an individual random intercept ----
# The random intercept represents an individual's general tendency to
# be terrestrial at the following location, after accounting for the
# current behavioural state.
model_individual_60 <- lme4::glmer(
  next_terrestrial ~
    behavior_from +
    (1 | individual_id),
  
  data = transitions_60,
  
  family = stats::binomial(
    link = "logit"
  ),
  
  control = glmer_control,
  
  nAGQ = 1)

print(summary(model_individual_60))


# 1.8 Inspect convergence ----
convergence_messages_60 <-
  model_individual_60@optinfo$conv$lme4$messages

if (is.null(convergence_messages_60)) {cat(
    "\nNo lme4 convergence warning was returned.\n")
  } else {
  print(convergence_messages_60)}

cat("\nSingular model:",lme4::isSingular(model_individual_60,tol = 1e-4),"\n")


# 1.9 Descriptive model comparison ----
# The AIC is reported as supporting information. The random-effect variance
# and singularity result must also be examined because zero variance lies on
# the boundary of the parameter space.
basic_model_comparison_60 <- stats::AIC(
  model_null_60,
  model_individual_60
) %>%
  as.data.frame() %>%
  tibble::rownames_to_column(
    "model"
  ) %>%
  arrange(
    AIC
  ) %>%
  mutate(
    delta_AIC = AIC - min(AIC)
  )

print(basic_model_comparison_60)
# model df              AIC        delta_AIC
# model_individual_60  3 48896.57    0.0000
#       model_null_60  2 49346.44  449.8712
# on retient donc le modèle avec effet individuel





#'------------------------------------------------------------------------------
# STEP 2: select covariate format ----
#'
#' (1) standardise continuous covariates;
#' (2) fit each covariate alone in linear and quadratic form;
#' (3) select the form with the lowest AIC;
#' (4) fit diel time as a cosinor model;
#' (5) prepare land cover as a categorical covariate.

# 2.1 Prepare covariates
# 2.1.1 Define continuous candidate covariates ----
control_covariates <- c("dem_elevation", "days_since_emig", "distance_to_nest_km", "ruggedness_100m", "slope_100m", "distance_to_ridgeline_100m")

hfi_covariates <- c("hfi_point", "hfi_mean_500m", "hfi_max_500m", "hfi_q75_500m", "hfi_mean_1000m", "hfi_max_1000m", "hfi_q75_1000m")

continuous_covariates <- c(control_covariates, hfi_covariates)

# 2.1.2 Standardise continuous covariates ----
# Standardisation improves numerical stability, particularly for quadratic
# terms. Original variables are retained.
standardise_variable <- function(x) {
  
  x <- as.numeric(x)
  
  valid_values <- is.finite(x)
  
  x_standardised <- rep(
    NA_real_,
    length(x)
  )
  
  variable_mean <- mean(
    x[valid_values]
  )
  
  variable_sd <- stats::sd(
    x[valid_values]
  )
  
  x_standardised[valid_values] <-
    (
      x[valid_values] -
        variable_mean
    ) /
    variable_sd
  
  x_standardised
}


transitions_60 <- transitions_60 %>%
  mutate(
    across(
      all_of(continuous_covariates),
      standardise_variable,
      .names = "{.col}_z"))

# 2.2 Functions used to inspect model fitting ----
model_converged <- function(model) {
  
  no_lme4_message <- is.null(
    model@optinfo$conv$lme4$messages
  )
  
  optimiser_converged <- isTRUE(
    model@optinfo$conv$opt == 0
  )
  
  no_lme4_message &&
    optimiser_converged}


# 2.4 Function to compare linear and quadratic forms ----
# Linear model: behavior_from * covariate
# Quadratic model: behavior_from * (covariate + covariate²)
#
# The interaction allows covariate effects to differ between transitions
# originating from aerial and terrestrial states.
fit_covariate_forms <- function(covariate, data) {
  
  standardised_covariate <- paste0(
    covariate,
    "_z")
  
  cat(
    "Fitting",
    covariate,
    "\n")
  
  flush.console()
  
  # 2.4.1 Use exactly the same observations for the two forms
  model_data <- data %>%
    filter(
      is.finite(
        .data[[standardised_covariate]]
      )
    ) %>%
    mutate(
      individual_id = droplevels(
        factor(individual_id)
      ),
      
      behavior_from = factor(
        behavior_from,
        levels = state_levels)
    )
  
  # 2.4.2 Linear form
  linear_formula <- stats::as.formula(
    paste0(
      "next_terrestrial ~ ",
      "behavior_from * ",
      standardised_covariate,
      " + (1 | individual_id)"
    )
  )
  
  model_linear <- lme4::glmer(
    formula = linear_formula,
    
    data = model_data,
    
    family = stats::binomial(
      link = "logit"),
    
    control = glmer_control,
    
    nAGQ = 1)
  
  # 2.4.3 Quadratic form
  quadratic_formula <- stats::as.formula(
    paste0(
      "next_terrestrial ~ ",
      "behavior_from * (",
      standardised_covariate,
      " + I(",
      standardised_covariate,
      "^2))",
      " + (1 | individual_id)"
    )
  )
  
  
  model_quadratic <- lme4::glmer(
    formula = quadratic_formula,
    
    data = model_data,
    
    family = stats::binomial(
      link = "logit"
    ),
    
    control = glmer_control,
    
    nAGQ = 1)
  
  # 2.4.4 Model comparison
  comparison <- tibble::tibble(
    covariate = covariate,
    
    form = c(
      "linear",
      "quadratic"
    ),
    
    n_observations = nrow(
      model_data
    ),
    
    n_individuals = nlevels(
      model_data$individual_id
    ),
    
    n_parameters = c(
      attr(
        stats::logLik(model_linear),
        "df"
      ),
      
      attr(
        stats::logLik(model_quadratic),
        "df"
      )
    ),
    
    AIC = c(
      stats::AIC(model_linear),
      stats::AIC(model_quadratic)
    ),
    
    converged = c(
      model_converged(model_linear),
      model_converged(model_quadratic)
    ),
    
    singular = c(
      lme4::isSingular(
        model_linear,
        tol = 1e-4
      ),
      
      lme4::isSingular(
        model_quadratic,
        tol = 1e-4
      )
    )
  )
  
  
  list(
    linear = model_linear,
    quadratic = model_quadratic,
    comparison = comparison
  )
}

# 2.5 Continuous covariates selection ----
# 2.5.1 Fit linear and quadratic forms for all continuous covariates ----
covariate_form_models_60 <- lapply(
  continuous_covariates,
  fit_covariate_forms,
  data = transitions_60)

names(covariate_form_models_60) <- continuous_covariates

# 2.5.2 Compile the AIC results ----
covariate_form_comparison_60 <- dplyr::bind_rows(
  lapply(
    covariate_form_models_60,
    function(model_results) {
      model_results$comparison
    }
  )
) %>%
  group_by(
    covariate
  ) %>%
  arrange(
    AIC,
    .by_group = TRUE
  ) %>%
  mutate(
    AIC_min = min(AIC),
    delta_AIC = AIC - AIC_min
  ) %>%
  ungroup()

print(covariate_form_comparison_60,n = Inf)


# 2.5.3 Retain the most parsimonious form ----
# The linear form is retained when its AIC is within 2 units
# of the minimum AIC. Otherwise, the quadratic form is retained.
covariate_form_comparison_60 <-
  covariate_form_comparison_60 %>%
  group_by(
    covariate
  ) %>%
  mutate(
    selected_form_parsimonious = case_when(
      
      form == "linear" &
        delta_AIC <= 2 ~ TRUE,
      
      form == "quadratic" &
        !any(
          form == "linear" &
            delta_AIC <= 2
        ) ~ TRUE,
      
      TRUE ~ FALSE
    )
  ) %>%
  ungroup()


selected_covariate_forms_60 <-
  covariate_form_comparison_60 %>%
  filter(
    selected_form_parsimonious
  ) %>%
  arrange(
    covariate
  )

print(selected_covariate_forms_60, n = Inf)


# 2.6 Fit diel time as a cosinor model ----
# cos_diel and sin_time are fitted together and must not be separated.
diel_data_60 <- transitions_60 %>%
  filter(
    is.finite(cos_diel),
    is.finite(sin_time)
  ) %>%
  mutate(
    individual_id = droplevels(
      factor(individual_id) ))


# Null model using exactly the same observations
model_diel_null_60 <- lme4::glmer(
  next_terrestrial ~
    behavior_from +
    (1 | individual_id),
  
  data = diel_data_60,
  
  family = stats::binomial(
    link = "logit"
  ),
  control = glmer_control,
  nAGQ = 1)


# Cosinor model
model_diel_cosinor_60 <- lme4::glmer(
  next_terrestrial ~
    behavior_from *
    (
      cos_diel +
        sin_time
    ) +
    (1 | individual_id),
  
  data = diel_data_60,
  
  family = stats::binomial(
    link = "logit"
  ),
  
  control = glmer_control,
  
  nAGQ = 1
)

diel_model_comparison_60 <- stats::AIC(
  model_diel_null_60,
  model_diel_cosinor_60
) %>%
  as.data.frame() %>%
  tibble::rownames_to_column(
    "model"
  ) %>%
  arrange(
    AIC
  ) %>%
  mutate(
    delta_AIC = AIC - min(AIC)
  )

print(diel_model_comparison_60)


# 2.7 Land-cover categorical model ----
if ("landcover_to" %in% names(transitions_60)) {
  
  landcover_data_60 <- transitions_60 %>%
    filter(
      !is.na(landcover_to)
    ) %>%
    mutate(
      landcover_to = factor(
        landcover_to,
        levels = c(
          "rocky_terrain",
          "low_vegetation",
          "forest"
        )
      ),
      
      individual_id = droplevels(
        factor(individual_id)
      )
    )
  
  
  model_landcover_null_60 <- lme4::glmer(
    next_terrestrial ~
      behavior_from +
      (1 | individual_id),
    
    data = landcover_data_60,
    
    family = stats::binomial(
      link = "logit"
    ),
    
    control = glmer_control,
    
    nAGQ = 1
  )
  
  
  model_landcover_factor_60 <- lme4::glmer(
    next_terrestrial ~
      behavior_from *
      landcover_to +
      (1 | individual_id),
    
    data = landcover_data_60,
    
    family = stats::binomial(
      link = "logit"
    ),
    
    control = glmer_control,
    
    nAGQ = 1
  )
  
  
  landcover_model_comparison_60 <- stats::AIC(
    model_landcover_null_60,
    model_landcover_factor_60
  ) %>%
    as.data.frame() %>%
    tibble::rownames_to_column(
      "model"
    ) %>%
    arrange(
      AIC
    ) %>%
    mutate(
      delta_AIC = AIC - min(AIC)
    )
  
  
  print(
    landcover_model_comparison_60
  )
  
} else {
  
  message(
    paste0(
      "landcover_to has not yet been created. ",
      "Reclassify landcover_100m before fitting the categorical model."
    )
  )
}

# synthèse
selected_variable_forms_60 <- selected_covariate_forms_60 %>%
  select(
    covariate,
    form,
    AIC,
    delta_AIC,
    converged,
    singular
  ) %>%
  bind_rows(
    tibble::tibble(
      covariate = "diel_cycle",
      form = "cosinor",
      AIC = stats::AIC(model_diel_cosinor_60),
      delta_AIC = 0,
      converged = model_converged(
        model_diel_cosinor_60
      ),
      singular = lme4::isSingular(
        model_diel_cosinor_60,
        tol = 1e-4
      )
    ),
    
    tibble::tibble(
      covariate = "landcover_to",
      form = "categorical",
      AIC = stats::AIC(
        model_landcover_factor_60
      ),
      delta_AIC = 0,
      converged = model_converged(
        model_landcover_factor_60
      ),
      singular = lme4::isSingular(
        model_landcover_factor_60,
        tol = 1e-4
      )
    )
  )

print(
  selected_variable_forms_60,
  n = Inf
)





#'------------------------------------------------------------------------------
# Step 3 : correlation amongst covariates
#' 
#' **Steps:** 
#' (1) assess correlation amongst covariates 
#' (2) for highly correlated variables (|r| > 0.5), retain the variable
#'     associated with the lowest AIC. 
#' (3) remove covariate for which the biological interpretation is less clear

