#' -----------------------------------------------------------------------------
# Title: Selection of variables ----
#' Authors : Louise Faure
#' Date : 16.07.26
#' 
#' Info : this script follow the Extract_covariates.R script where covariates are
#' extracted below each location and within two buffers. 
#' 
# Main steps: ----
#' (1) adjust the basic model and inspect individual variation
#' (2) prepare covariates under linear, cosinor and quadriatique forms
#' (3) select which form of HFI has the best AIC
#' (4) filter the control covariates by first identifying amongst correlated 
#' variables (|r| > 0.7) (Dormann et al., 2013) the one with highest AIC and 
#' better biological explanation
#' (5) conduct forward and backward model selection on the selected variables
#' (6) export dataset with retained variables and form
#' 
#' Ref. the statistical selection of the variables and model is based on Togunov
#' et al., Mov Ecol, 2022. 
#' -----------------------------------------------------------------------------


# library 
library(dplyr)
library(tidyr)
library(tibble)
library(lme4)
library(corrplot)

# data 
GE_60_min_covariates_hfi <- readRDS("/Users/louisefaure/Desktop/dossier sans titre/donnees filtree/GE_60_min_covariates_hfi.rds")

# parameters
state_levels <- c("aerial", "terrestrial")
transition_levels <- c("aerial_to_aerial", "aerial_to_terrestrial", "terrestrial_to_aerial", "terrestrial_to_terrestrial")
correlation_threshold <- 0.5
glmer_control <- lme4::glmerControl(optimizer = "bobyqa",optCtrl = list(maxfun = 200000))



#'------------------------------------------------------------------------------
# STEP 1: adjust basic model with HFI ----
#'
#' (1) prepare the transition dataset and transition matrix;
#' (2) fit the null transition model that include HFI;
#' (3) fit a model with HFI + random intercept
#' (4) inspect AIC between two models, and chose the one with highest AIC


# 1.1 Prepare transitions between consecutive locations ----
# A transition is created only between consecutive locations belonging
# to the same individual and the same burst_n.
transitions_60 <- GE_60_min_covariates_hfi %>%
  dplyr::arrange(
    individual.local.identifier,
    burst_n,
    timestamp
  ) %>%
  dplyr::mutate(
    timestamp = as.POSIXct(
      timestamp,
      tz = "UTC"
    ),
    
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
    )
  ) %>%
  dplyr::group_by(
    individual.local.identifier,
    burst_n
  ) %>%
  dplyr::mutate(
    # Behaviour at the arrival point
    behavior_to = dplyr::lead(
      behavior_from
    ),
    
    # Land-cover proportions at the arrival point
    prop_forest_5cells_to = dplyr::lead(
      prop_forest_5cells
    ),
    
    prop_low_vegetation_5cells_to = dplyr::lead(
      prop_low_vegetation_5cells
    ),
    
    prop_rocky_terrain_5cells_to = dplyr::lead(
      prop_rocky_terrain_5cells
    ),
    
    prop_other_5cells_to = dplyr::lead(
      prop_other_5cells
    ),
    
    # Time at the arrival point
    timestamp_next = dplyr::lead(
      timestamp
    ),
    
    # Time elapsed between the two locations
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
    !is.na(behavior_to)
  ) %>%
  dplyr::mutate(
    individual_id = factor(
      individual.local.identifier
    ),
    
    # Binary response:
    # 1 = arrival state is terrestrial
    # 0 = arrival state is aerial
    next_terrestrial = as.integer(
      behavior_to == "terrestrial"
    ),
    
    transition_type = paste(
      behavior_from,
      behavior_to,
      sep = "_to_"
    )
  )

# 1.1.2 Nbr of transitions & bursts per individual and thinning ----
transitions_by_individual <- transitions_60 %>%
  group_by(
    individual.local.identifier
  ) %>%
  summarise(
    n_transitions = n(),
    n_bursts = n_distinct(burst_n),
    .groups = "drop") %>%
  arrange(
    n_transitions)

print(transitions_by_individual, n = Inf)
# Langgries21 (eobs 7586) and Almen18 (eobs 5861) have respectively
# 3 and 34 transitions and are therefore removed.

individuals_to_remove <- c("Langgries21 (eobs 7586)","Almen18 (eobs 5861)")
transitions_60 <- transitions_60 %>% filter(!individual.local.identifier %in% individuals_to_remove)

# 1.1.3 Empirical transition-count matrix ----
transition_count_matrix_60 <- xtabs(
  ~ behavior_from + behavior_to,
  data = transitions_60)

print(transition_count_matrix_60)
#                  behavior_to
# behavior_from aerial terrestrial
# aerial        4259        6652
# terrestrial   6518       31257

# 1.1.4 Empirical transition-probability matrix ----
transition_probability_matrix_60 <- prop.table(
  transition_count_matrix_60,
  margin = 1)

print(round(transition_probability_matrix_60,digits = 4))
#                 behavior_to
# behavior_from aerial terrestrial
# aerial      0.3903      0.6097
# terrestrial 0.1725      0.8275

# 1.2 Fit the null Markov chains ----
# With aerial as the reference:
model_null_60 <- stats::glm(
  next_terrestrial ~ behavior_from,
  data = transitions_60,
  family = stats::binomial(
    link = "logit"))

print(summary(model_null_60))

# 1.2.1 Transition probabilities estimated by the null model ----
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

# 1.3 Fit a model with an individual random intercept ----
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


# 1.3.1 Inspect convergence ----
convergence_messages_60 <-
  model_individual_60@optinfo$conv$lme4$messages

if (is.null(convergence_messages_60)) {cat(
    "\nNo lme4 convergence warning was returned.\n")
  } else {
  print(convergence_messages_60)}

cat("\nSingular model:",lme4::isSingular(model_individual_60,tol = 1e-4),"\n")


# 1.4 Descriptive model comparison ----
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
# we retain the model with 





#'------------------------------------------------------------------------------
# STEP 2: select covariate format ----
#'
#' (1) standardise continuous covariates;
#' (2) fit each covariate alone in linear and quadratic form;
#' (3) select the form with the lowest AIC;
#' (4) fit diel time as a cosinor model;


# 2.1 Prepare covariates
# 2.1.0 Add a column which for each location calculate the time spend in the 
# behavioral state (terrestrial or aerian, restart time count when their is a burst)

# 2.1.1 Define continuous candidate covariates ----
control_covariates <- c("dem_elevation", "days_since_emig", "distance_to_nest_km", "ruggedness_100m", "slope_100m", "distance_to_ridgeline_100m", "prop_forest_5cells_to","prop_low_vegetation_5cells_to","prop_rocky_terrain_5cells_to")
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
            delta_AIC <= 4
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
    delta_AIC = AIC - min(AIC))

print(diel_model_comparison_60)

# 2.7 Synthesis of retained variable forms ----
selected_variable_forms_60 <- selected_covariate_forms_60 %>%
  dplyr::select(
    covariate,
    form,
    AIC,
    delta_AIC,
    converged,
    singular
  ) %>%
  dplyr::bind_rows(
    tibble::tibble(
      covariate = "diel_cycle",
      form = "cosinor",
      AIC = stats::AIC(
        model_diel_cosinor_60
      ),
      delta_AIC = 0,
      converged = model_converged(
        model_diel_cosinor_60
      ),
      singular = lme4::isSingular(
        model_diel_cosinor_60,
        tol = 1e-4
      )
    )
  )

print(selected_variable_forms_60,n = Inf)


#'------------------------------------------------------------------------------
# STEP 3: identify the best HFI formulation using AIC ----
#'
#' (1) retrieve the selected linear or quadratic form of each HFI variable;
#' (2) refit all HFI models on exactly the same observations;
#' (3) compare point, 500-m and 1000-m HFI formulations using AIC;
#' (4) retain the HFI formulation with the lowest AIC.


# 3.1 Define HFI candidate formulations ----

hfi_candidate_information_60 <- tibble::tibble(
  covariate = c(
    "hfi_point",
    "hfi_mean_500m",
    "hfi_max_500m",
    "hfi_q75_500m",
    "hfi_mean_1000m",
    "hfi_max_1000m",
    "hfi_q75_1000m"
  ),
  
  spatial_scale = c(
    "point",
    "500m",
    "500m",
    "500m",
    "1000m",
    "1000m",
    "1000m"
  ),
  
  buffer_statistic = c(
    "point_value",
    "mean",
    "maximum",
    "q75",
    "mean",
    "maximum",
    "q75"
  )
)


# 3.2 Retrieve the selected form of each HFI covariate ----

selected_hfi_forms_60 <- selected_covariate_forms_60 %>%
  dplyr::filter(
    covariate %in% hfi_covariates
  ) %>%
  dplyr::select(
    covariate,
    selected_form = form
  ) %>%
  dplyr::left_join(
    hfi_candidate_information_60,
    by = "covariate"
  )


# Check that one form was selected for every HFI covariate

missing_hfi_forms_60 <- setdiff(
  hfi_covariates,
  selected_hfi_forms_60$covariate
)

if (length(missing_hfi_forms_60) > 0) {
  stop(
    paste(
      "No selected form was found for:",
      paste(
        missing_hfi_forms_60,
        collapse = ", "
      )
    )
  )
}


duplicated_hfi_forms_60 <- selected_hfi_forms_60 %>%
  dplyr::count(
    covariate
  ) %>%
  dplyr::filter(
    n != 1
  )

if (nrow(duplicated_hfi_forms_60) > 0) {
  stop(
    "Each HFI covariate must have exactly one selected form."
  )
}

print(
  selected_hfi_forms_60,
  n = Inf
)


# 3.3 Create a common complete-case dataset ----
#
# Every HFI model must use exactly the same transitions for the AIC comparison.

required_hfi_columns_60 <- paste0(
  hfi_covariates,
  "_z"
)

hfi_common_data_60 <- transitions_60 %>%
  dplyr::filter(
    dplyr::if_all(
      dplyr::all_of(required_hfi_columns_60),
      ~ is.finite(as.numeric(.x))
    )
  ) %>%
  dplyr::mutate(
    individual_id = droplevels(
      factor(individual_id)
    ),
    
    behavior_from = factor(
      behavior_from,
      levels = state_levels
    )
  )


cat(
  "\nObservations used for all HFI models:",
  nrow(hfi_common_data_60),
  "\n"
)

cat(
  "Individuals used for all HFI models:",
  nlevels(hfi_common_data_60$individual_id),
  "\n"
)


# 3.4 Fit the null mixed model on the common HFI dataset ----

model_hfi_null_60 <- lme4::glmer(
  next_terrestrial ~
    behavior_from +
    (1 | individual_id),
  
  data = hfi_common_data_60,
  
  family = stats::binomial(
    link = "logit"
  ),
  
  control = glmer_control,
  
  nAGQ = 1
)


# 3.5 Function to fit the selected form of one HFI covariate ----

fit_selected_hfi_model <- function(
    covariate,
    data,
    selected_forms
) {
  
  selected_form <- selected_forms %>%
    dplyr::filter(
      .data$covariate == .env$covariate
    ) %>%
    dplyr::pull(
      selected_form
    )
  
  
  if (length(selected_form) != 1) {
    stop(
      paste(
        "No unique selected form was found for",
        covariate
      )
    )
  }
  
  
  hfi_z <- paste0(
    covariate,
    "_z"
  )
  
  
  # Linear HFI model
  
  if (selected_form == "linear") {
    
    model_formula <- stats::as.formula(
      paste0(
        "next_terrestrial ~ ",
        "behavior_from * ",
        hfi_z,
        " + (1 | individual_id)"
      )
    )
  }
  
  
  # Quadratic HFI model
  
  if (selected_form == "quadratic") {
    
    model_formula <- stats::as.formula(
      paste0(
        "next_terrestrial ~ ",
        "behavior_from * (",
        hfi_z,
        " + I(",
        hfi_z,
        "^2))",
        " + (1 | individual_id)"
      )
    )
  }
  
  
  if (!selected_form %in% c("linear", "quadratic")) {
    stop(
      paste(
        "Unsupported form for",
        covariate,
        ":",
        selected_form
      )
    )
  }
  
  
  fitted_model <- lme4::glmer(
    formula = model_formula,
    
    data = data,
    
    family = stats::binomial(
      link = "logit"
    ),
    
    control = glmer_control,
    
    nAGQ = 1
  )
  
  
  list(
    covariate = covariate,
    selected_form = selected_form,
    formula = model_formula,
    model = fitted_model,
    AIC = stats::AIC(fitted_model),
    logLik = as.numeric(
      stats::logLik(fitted_model)
    ),
    n_parameters = attr(
      stats::logLik(fitted_model),
      "df"
    ),
    converged = model_converged(
      fitted_model
    ),
    singular = lme4::isSingular(
      fitted_model,
      tol = 1e-4
    )
  )
}


# 3.6 Fit all seven HFI candidate models ----

hfi_candidate_models_60 <- lapply(
  hfi_covariates,
  fit_selected_hfi_model,
  data = hfi_common_data_60,
  selected_forms = selected_hfi_forms_60
)

names(hfi_candidate_models_60) <-
  hfi_covariates


# 3.7 Compile HFI model-comparison results ----

hfi_model_comparison_60 <- dplyr::bind_rows(
  lapply(
    hfi_candidate_models_60,
    function(x) {
      
      tibble::tibble(
        covariate = x$covariate,
        selected_form = x$selected_form,
        n_observations = nrow(
          hfi_common_data_60
        ),
        n_individuals = nlevels(
          hfi_common_data_60$individual_id
        ),
        n_parameters = x$n_parameters,
        logLik = x$logLik,
        AIC = x$AIC,
        converged = x$converged,
        singular = x$singular
      )
    }
  )
) %>%
  dplyr::left_join(
    hfi_candidate_information_60,
    by = "covariate"
  ) %>%
  dplyr::arrange(
    AIC
  ) %>%
  dplyr::mutate(
    delta_AIC = AIC - min(AIC),
    
    Akaike_weight = exp(
      -0.5 * delta_AIC
    ) /
      sum(
        exp(
          -0.5 * delta_AIC
        )
      ),
    
    selected_HFI = dplyr::row_number() == 1
  ) %>%
  dplyr::select(
    covariate,
    spatial_scale,
    buffer_statistic,
    selected_form,
    n_observations,
    n_individuals,
    n_parameters,
    logLik,
    AIC,
    delta_AIC,
    Akaike_weight,
    converged,
    singular,
    selected_HFI
  )

print(
  hfi_model_comparison_60,
  n = Inf
)


# 3.8 Compare each HFI model with the common null model ----

null_hfi_AIC_60 <- stats::AIC(
  model_hfi_null_60
)

hfi_model_comparison_60 <- hfi_model_comparison_60 %>%
  dplyr::mutate(
    AIC_null = null_hfi_AIC_60,
    
    AIC_improvement_over_null =
      AIC_null - AIC
  )

print(
  hfi_model_comparison_60,
  n = Inf
)


# 3.9 Retain the best HFI formulation ----

selected_hfi_covariate_60 <-
  hfi_model_comparison_60$covariate[1]

selected_hfi_form_60 <-
  hfi_model_comparison_60$selected_form[1]

selected_hfi_scale_60 <-
  hfi_model_comparison_60$spatial_scale[1]

selected_hfi_statistic_60 <-
  hfi_model_comparison_60$buffer_statistic[1]


model_hfi_selected_60 <-
  hfi_candidate_models_60[[selected_hfi_covariate_60]]$model


cat(
  "\nSelected HFI covariate:",
  selected_hfi_covariate_60,
  "\n"
)

cat(
  "Selected form:",
  selected_hfi_form_60,
  "\n"
)

cat(
  "Selected spatial scale:",
  selected_hfi_scale_60,
  "\n"
)

cat(
  "Selected buffer statistic:",
  selected_hfi_statistic_60,
  "\n"
)

cat(
  "Selected-model AIC:",
  round(
    stats::AIC(model_hfi_selected_60),
    2
  ),
  "\n"
)


# 3.10 Inspect the selected HFI model ----

print(
  summary(
    model_hfi_selected_60
  )
)













#'------------------------------------------------------------------------------
# Step 3 : correlation amongst covariates ----
#' 
#' **Steps:** 
#' (1) calculate correlations among continuous covariates;
#' (2) identify highly correlated pairs (|r| > 0.5);
#' (3) select one competing land-cover proportion;
#' (4) retain non-redundant control covariates using AIC;
#' (5) inspect alternative HFI formulations separately;
#' (6) assess associations with the diel cycle;
#' (7) apply biological screening.


# 3.1 Parameters and covariate groups ----
# Topographic, temporal and spatial control variables
topographic_control_covariates <- control_covariates

# Land-cover proportions are competing descriptions of the same composition.
# They should not all be included together in the final model.
landcover_control_covariates <- landcover_covariates


# 3.2 Correlation matrix among all continuous candidate covariates ----
# Correlations are calculated on the original variables.
correlation_data_60 <- transitions_60 %>%
  dplyr::select(
    dplyr::all_of(continuous_covariates)
  ) %>%
  dplyr::mutate(
    dplyr::across(
      dplyr::everything(),
      as.numeric
    )
  )

correlation_matrix_60 <- stats::cor(
  correlation_data_60,
  use = "pairwise.complete.obs",
  method = "pearson"
)

#' ----------------------------------------------------------------------------- Visualisation n°1: correlation matrix (beginning)
corrplot::corrplot(
  correlation_matrix_60,
  
  # Display correlations as coloured cells
  method = "color",
  
  # Display only the lower triangle
  type = "lower",
  
  # Preserve the order defined in continuous_covariates
  order = "original",
  
  # Diverging palette for correlations from -1 to 1
  col = corrplot::COL2(
    "RdBu",
    200
  ),
  
  col.lim = c(
    -1,
    1
  ),
  
  # Do not display the diagonal
  diag = FALSE,
  
  # Add correlation values inside cells
  addCoef.col = "black",
  number.digits = 2,
  number.cex = 0.55,
  
  # Variable labels
  tl.col = "black",
  tl.cex = 0.65,
  tl.srt = 45,
  
  # Legend
  cl.cex = 0.8,
  
  # Cell borders
  addgrid.col = "white",
  
  # Representation of unavailable correlations
  na.label = "×",
  na.label.col = "grey50",
  
  title = "Pearson correlations among continuous covariates",
  
  mar = c(
    0,
    0,
    2,
    0
  )
)
#'------------------------------------------------------------------------------ Visualisation n°1 : correlation matrix (end)


# 3.3 Convert the correlation matrix to a pairwise table ----
correlation_pairs_60 <- as.data.frame(
  as.table(correlation_matrix_60),
  stringsAsFactors = FALSE
) %>%
  tibble::as_tibble() %>%
  dplyr::rename(
    covariate_1 = Var1,
    covariate_2 = Var2,
    correlation = Freq
  ) %>%
  dplyr::mutate(
    covariate_1 = as.character(covariate_1),
    covariate_2 = as.character(covariate_2)
  ) %>%
  dplyr::filter(
    covariate_1 < covariate_2
  ) %>%
  dplyr::mutate(
    absolute_correlation = abs(correlation),
    
    group_1 = dplyr::case_when(
      covariate_1 %in% control_covariates ~ "control",
      covariate_1 %in% landcover_covariates ~ "landcover",
      covariate_1 %in% hfi_covariates ~ "hfi",
      TRUE ~ "other"
    ),
    
    group_2 = dplyr::case_when(
      covariate_2 %in% control_covariates ~ "control",
      covariate_2 %in% landcover_covariates ~ "landcover",
      covariate_2 %in% hfi_covariates ~ "hfi",
      TRUE ~ "other"
    ),
    
    pair_type = paste(
      pmin(group_1, group_2),
      pmax(group_1, group_2),
      sep = "_"
    )
  ) %>%
  dplyr::arrange(
    dplyr::desc(absolute_correlation)
  )


# 3.4 Identify highly correlated pairs ----
# 3.4 Identify highly correlated pairs among control covariates ----
# 3.4 Identify highly correlated pairs ----
#
# HFI-HFI pairs are excluded because HFI formulations are competing
# alternatives and will never be included together in the same model.

highly_correlated_pairs_60 <- correlation_pairs_60 %>%
  dplyr::filter(
    absolute_correlation > correlation_threshold,
    pair_type != "hfi_hfi"
  )

print(
  highly_correlated_pairs_60,
  n = Inf
)

# 3.5 Retrieve the selected form and AIC of each covariate ----
selected_aic_60 <- selected_covariate_forms_60 %>%
  dplyr::select(
    covariate,
    form,
    n_observations,
    AIC
  ) %>%
  dplyr::rename(
    selected_form = form
  )

# 3.6 Retain the lowest-AIC covariate within correlated control pairs ----
#
# HFI-control pairs are ignored at this stage.
# Only correlations among control and land-cover covariates are considered.

correlated_control_pairs_60 <- highly_correlated_pairs_60 %>%
  dplyr::filter(
    !covariate_1 %in% hfi_covariates,
    !covariate_2 %in% hfi_covariates
  )

print(
  correlated_control_pairs_60,
  n = Inf
)


# Add the selected form and AIC of both covariates
correlated_control_pairs_AIC_60 <-
  correlated_control_pairs_60 %>%
  dplyr::left_join(
    selected_aic_60 %>%
      dplyr::rename(
        covariate_1 = covariate,
        selected_form_1 = selected_form,
        n_observations_1 = n_observations,
        AIC_1 = AIC
      ),
    by = "covariate_1"
  ) %>%
  dplyr::left_join(
    selected_aic_60 %>%
      dplyr::rename(
        covariate_2 = covariate,
        selected_form_2 = selected_form,
        n_observations_2 = n_observations,
        AIC_2 = AIC
      ),
    by = "covariate_2"
  ) %>%
  dplyr::mutate(
    same_dataset =
      n_observations_1 == n_observations_2,
    
    retained_covariate = dplyr::case_when(
      !same_dataset ~ NA_character_,
      AIC_1 < AIC_2 ~ covariate_1,
      AIC_2 < AIC_1 ~ covariate_2,
      TRUE ~ "equal_AIC"
    ),
    
    removed_covariate = dplyr::case_when(
      !same_dataset ~ NA_character_,
      AIC_1 < AIC_2 ~ covariate_2,
      AIC_2 < AIC_1 ~ covariate_1,
      TRUE ~ NA_character_
    ),
    
    delta_AIC_between_covariates =
      abs(AIC_1 - AIC_2)
  ) %>%
  dplyr::arrange(
    dplyr::desc(absolute_correlation)
  )

print(
  correlated_control_pairs_AIC_60,
  n = Inf
)
#' covariate_1           covariate_2      correlation
#' ruggedness_100m       slope_100m             0.899                
#' prop_forest_5cells_to prop_low_vegeta…      -0.655             
#' dem_elevation         prop_forest_5ce…      -0.549
#' we decide to keep the dem as it was also correlated with HFI and remove forest, 
#' we choose the best form between slope and ruggedness using AIC. 
#' 


# 3.X Apply correlation and biological screening decisions ----
#
# Candidate covariates represented distinct but partly overlapping ecological
# mechanisms. Highly correlated covariates were not retained simultaneously
# when this would impair model interpretation.
#
# Decisions:
# (1) ruggedness_100m is retained over slope_100m because the two variables
#     are strongly correlated and ruggedness had the lower selected-model AIC;
# (2) dem_elevation is retained as an essential control for the altitudinal
#     gradient in HFI;
# (3) prop_forest_5cells_to is removed because it is correlated with elevation;
# (4) the remaining land-cover proportions are removed from the main control
#     model because they represent competing components of local land-cover
#     composition;
# (5) distance_to_nest_km is removed because its ecological interpretation
#     overlaps with progression through dispersal, already represented by
#     days_since_emig.


# Controls explicitly retained for forward selection

retained_control_covariates_60 <- c(
  "dem_elevation",
  "distance_to_ridgeline_100m",
  "days_since_emig",
  "ruggedness_100m"
)


# Controls explicitly removed

removed_control_covariates_60 <- c(
  "slope_100m",
  "prop_forest_5cells_to",
  "prop_low_vegetation_5cells_to",
  "prop_rocky_terrain_5cells_to",
  "distance_to_nest_km"
)


# Ensure that all named variables belong to the candidate set

all_control_candidates_60 <- c(
  control_covariates,
  landcover_covariates
)

unknown_retained_covariates_60 <- setdiff(
  retained_control_covariates_60,
  all_control_candidates_60
)

unknown_removed_covariates_60 <- setdiff(
  removed_control_covariates_60,
  all_control_candidates_60
)

if (length(unknown_retained_covariates_60) > 0) {
  stop(
    paste(
      "Unknown retained control covariates:",
      paste(
        unknown_retained_covariates_60,
        collapse = ", "
      )
    )
  )
}

if (length(unknown_removed_covariates_60) > 0) {
  stop(
    paste(
      "Unknown removed control covariates:",
      paste(
        unknown_removed_covariates_60,
        collapse = ", "
      )
    )
  )
}


cat(
  "\nControl covariates removed after correlation and biological screening:\n"
)

print(
  removed_control_covariates_60
)


cat(
  "\nControl covariates retained for forward selection:\n"
)

print(
  retained_control_covariates_60
)


# Summary table of biological decisions
control_screening_summary_60 <- selected_aic_60 %>%
  dplyr::filter(
    covariate %in% all_control_candidates_60
  ) %>%
  dplyr::mutate(
    decision = dplyr::case_when(
      covariate %in% retained_control_covariates_60 ~
        "retained_for_forward_selection",
      
      covariate %in% removed_control_covariates_60 ~
        "removed_after_screening",
      
      TRUE ~
        "not_classified"
    ),
    
    reason = dplyr::case_when(
      covariate == "dem_elevation" ~
        "essential control for the altitudinal HFI gradient",
      
      covariate == "distance_to_ridgeline_100m" ~
        "proxy for potential access to orographic uplift",
      
      covariate == "days_since_emig" ~
        "proxy for behavioural development during dispersal",
      
      covariate == "ruggedness_100m" ~
        "retained over slope based on correlation and AIC",
      
      covariate == "slope_100m" ~
        "strongly correlated with ruggedness",
      
      covariate == "prop_forest_5cells_to" ~
        "correlated with elevation",
      
      covariate %in% c(
        "prop_low_vegetation_5cells_to",
        "prop_rocky_terrain_5cells_to"
      ) ~
        "competing land-cover composition variable",
      
      covariate == "distance_to_nest_km" ~
        "removed in favour of days since emigration",
      
      TRUE ~
        NA_character_
    )
  ) %>%
  dplyr::arrange(
    decision,
    AIC
  )

print(
  control_screening_summary_60,
  n = Inf
)





#'------------------------------------------------------------------------------
# STEP 4: forward selection of control covariates ----
#'
#' (1) start from the null model with an individual random intercept;
#' (2) add each retained control covariate separately;
#' (3) retain the covariate producing the lowest AIC;
#' (4) repeat until no candidate decreases AIC by at least 2 units.



# 4.1 Parameters ----

aic_improvement_threshold <- 4

selected_hfi_covariate_60 <- "hfi_mean_1000m"
selected_hfi_form_60 <- "quadratic"

# Only control covariates enter the forward-selection loop.
forward_candidate_covariates_60 <- c(
  retained_control_covariates_60,
  "diel_cycle"
)


# 4.2 Create a common complete-case dataset ----
#
# HFI is forced into every model.
# All candidate models use exactly the same observations.

continuous_forward_covariates_60 <- setdiff(
  forward_candidate_covariates_60,
  "diel_cycle"
)

required_forward_columns_60 <- c(
  paste0(
    continuous_forward_covariates_60,
    "_z"
  ),
  paste0(
    selected_hfi_covariate_60,
    "_z"
  ),
  "cos_diel",
  "sin_time"
)

forward_data_60 <- transitions_60 %>%
  dplyr::filter(
    dplyr::if_all(
      dplyr::all_of(required_forward_columns_60),
      ~ is.finite(as.numeric(.x))
    )
  ) %>%
  dplyr::mutate(
    individual_id = droplevels(
      factor(individual_id)
    ),
    
    behavior_from = factor(
      behavior_from,
      levels = state_levels
    )
  )

cat(
  "\nObservations used for HFI-conditioned forward selection:",
  nrow(forward_data_60),
  "\n"
)

cat(
  "Individuals used for HFI-conditioned forward selection:",
  nlevels(forward_data_60$individual_id),
  "\n"
)

# 4.3 Construct the model term for a candidate covariate ----

# 4.3 Construct the model term for a candidate covariate ----

create_candidate_term <- function(
    covariate,
    selected_forms
) {
  
  # Diel cycle: cos_diel and sin_time must remain together
  if (covariate == "diel_cycle") {
    
    return(
      paste(
        "cos_diel",
        "sin_time",
        "behavior_from:cos_diel",
        "behavior_from:sin_time",
        sep = " + "
      )
    )
  }
  
  
  # Retrieve the selected linear or quadratic form
  selected_form <- selected_forms %>%
    dplyr::filter(
      .data$covariate == .env$covariate
    ) %>%
    dplyr::pull(
      selected_form
    )
  
  
  if (length(selected_form) == 0) {
    
    stop(
      paste(
        "No selected form was found for",
        covariate
      )
    )
  }
  
  
  if (length(selected_form) > 1) {
    
    stop(
      paste(
        "More than one selected form was found for",
        covariate,
        ":",
        paste(selected_form, collapse = ", ")
      )
    )
  }
  
  
  covariate_z <- paste0(
    covariate,
    "_z"
  )
  
  
  # Linear form
  if (selected_form == "linear") {
    
    return(
      paste(
        covariate_z,
        paste0(
          "behavior_from:",
          covariate_z
        ),
        sep = " + "
      )
    )
  }
  
  
  # Quadratic form
  if (selected_form == "quadratic") {
    
    quadratic_term <- paste0(
      "I(",
      covariate_z,
      "^2)"
    )
    
    return(
      paste(
        covariate_z,
        quadratic_term,
        paste0(
          "behavior_from:",
          covariate_z
        ),
        paste0(
          "behavior_from:",
          quadratic_term
        ),
        sep = " + "
      )
    )
  }
  
  
  stop(
    paste(
      "Unsupported form for",
      covariate,
      ":",
      selected_form
    )
  )
}

# 4.3.1 Construct the mandatory HFI term ----

forced_hfi_term_60 <- create_candidate_term(
  covariate = selected_hfi_covariate_60,
  selected_forms = selected_aic_60
)

cat(
  "\nMandatory HFI term:\n",
  forced_hfi_term_60,
  "\n"
)
# 4.4 Forward AIC selection ----

selected_forward_covariates_60 <- character(0)

remaining_forward_covariates_60 <-
  forward_candidate_covariates_60

current_fixed_terms_60 <- paste(
  "behavior_from",
  forced_hfi_term_60,
  sep = " + "
)

current_fixed_terms_60 <- paste(
  "behavior_from",
  forced_hfi_term_60,
  sep = " + "
)

current_AIC_60 <- stats::AIC(
  current_model_60
)


forward_selection_log_60 <- tibble::tibble(
  step = 0L,
  added_covariate = "HFI_base_model",
  AIC = current_AIC_60,
  AIC_improvement = NA_real_,
  selected_covariates = selected_hfi_covariate_60
)


step_number <- 1L

# 4.4 Initialise HFI-conditioned forward selection ----

selected_forward_covariates_60 <- character(0)

remaining_forward_covariates_60 <-
  forward_candidate_covariates_60


# Mandatory HFI term retained in every model

forced_hfi_term_60 <- create_candidate_term(
  covariate = selected_hfi_covariate_60,
  selected_forms = selected_aic_60
)


# Initial fixed effects:
# current behavioural state + selected HFI formulation

current_fixed_terms_60 <- paste(
  "behavior_from",
  forced_hfi_term_60,
  sep = " + "
)


# Build the initial HFI model formula

current_formula_60 <- stats::as.formula(
  paste0(
    "next_terrestrial ~ ",
    current_fixed_terms_60,
    " + (1 | individual_id)"
  )
)


# Fit the initial HFI model

current_model_60 <- lme4::glmer(
  formula = current_formula_60,
  data = forward_data_60,
  family = stats::binomial(
    link = "logit"
  ),
  control = glmer_control,
  nAGQ = 1
)


# AIC of the HFI base model

current_AIC_60 <- stats::AIC(
  current_model_60
)


cat(
  "\nInitial HFI model AIC:",
  round(current_AIC_60, 2),
  "\n"
)

print(
  current_formula_60
)


# Initialise selection log

forward_selection_log_60 <- tibble::tibble(
  step = 0L,
  added_covariate = "HFI_base_model",
  AIC = current_AIC_60,
  AIC_improvement = NA_real_,
  selected_covariates = selected_hfi_covariate_60
)


step_number <- 1L

while (
  length(remaining_forward_covariates_60) > 0
) {
  
  candidate_models <- vector(
    mode = "list",
    length = length(
      remaining_forward_covariates_60
    )
  )
  
  names(candidate_models) <-
    remaining_forward_covariates_60
  
  
  candidate_results <- lapply(
    remaining_forward_covariates_60,
    function(candidate_covariate) {
      
      candidate_term <- create_candidate_term(
        covariate = candidate_covariate,
        selected_forms = selected_aic_60
      )
      
      candidate_fixed_terms <- paste(
        current_fixed_terms_60,
        candidate_term,
        sep = " + "
      )
      
      candidate_formula <- stats::as.formula(
        paste0(
          "next_terrestrial ~ ",
          candidate_fixed_terms,
          " + (1 | individual_id)"
        )
      )
      
      candidate_model <- lme4::glmer(
        formula = candidate_formula,
        data = forward_data_60,
        family = stats::binomial(
          link = "logit"
        ),
        control = glmer_control,
        nAGQ = 1
      )
      
      list(
        covariate = candidate_covariate,
        term = candidate_term,
        model = candidate_model,
        AIC = stats::AIC(candidate_model),
        converged = model_converged(candidate_model),
        singular = lme4::isSingular(
          candidate_model,
          tol = 1e-4
        )
      )
    }
  )
  
  names(candidate_results) <- remaining_forward_covariates_60
  
  candidate_comparison <- dplyr::bind_rows(
    lapply(
      candidate_results,
      function(x) {
        tibble::tibble(
          covariate = x$covariate,
          AIC = x$AIC,
          converged = x$converged,
          singular = x$singular
        )
      }
    )
  ) %>%
    dplyr::arrange(
      AIC
    )
  
  
  print(
    candidate_comparison,
    n = Inf
  )
  
  
  best_candidate <- candidate_comparison$covariate[1]
  
  best_candidate_AIC <- candidate_comparison$AIC[1]
  
  AIC_improvement <-
    current_AIC_60 -
    best_candidate_AIC
  
  
  if (
    AIC_improvement <
    aic_improvement_threshold
  ) {
    
    cat(
      "\nForward selection stopped: no remaining variable improved AIC by at least",
      aic_improvement_threshold,
      "units.\n"
    )
    
    break
  }
  
  
  best_result <- candidate_results[[best_candidate]]
  
  
  selected_forward_covariates_60 <- c(
    selected_forward_covariates_60,
    best_candidate
  )
  
  remaining_forward_covariates_60 <- setdiff(
    remaining_forward_covariates_60,
    best_candidate
  )
  
  current_fixed_terms_60 <- paste(
    current_fixed_terms_60,
    best_result$term,
    sep = " + "
  )
  
  current_model_60 <- best_result$model
  
  current_AIC_60 <- best_candidate_AIC
  
  
  forward_selection_log_60 <- dplyr::bind_rows(
    forward_selection_log_60,
    
    tibble::tibble(
      step = step_number,
      added_covariate = best_candidate,
      AIC = current_AIC_60,
      AIC_improvement = AIC_improvement,
      selected_covariates = paste(
        selected_forward_covariates_60,
        collapse = " + "
      )
    )
  )
  
  
  cat(
    "\nStep",
    step_number,
    "- retained:",
    best_candidate,
    "- AIC:",
    round(
      current_AIC_60,
      2
    ),
    "- improvement:",
    round(
      AIC_improvement,
      2
    ),
    "\n"
  )
  
  
  step_number <- step_number + 1L
}

# 4.5 Forward-selection results ----

print(
  forward_selection_log_60,
  n = Inf
)


# step added_covariate                 AIC AIC_improvement selected_covariates                                                                                                             
# step added_covariate               AIC AIC_improvement selected_covariates                                                      
# <int> <chr>                       <dbl>           <dbl> <chr>                                                                    
#   1     0 HFI_base_model             48859.            NA   hfi_mean_1000m                                                           
# 2     1 diel_cycle                 47422.          1438.  diel_cycle                                                               
# 3     2 distance_to_ridgeline_100m 47271.           150.  diel_cycle + distance_to_ridgeline_100m                                  
# 4     3 dem_elevation              47195.            76.1 diel_cycle + distance_to_ridgeline_100m + dem_elevation                  
# 5     4 days_since_emig            47176.            19.1 diel_cycle + distance_to_ridgeline_100m + dem_elevation + days_since_emig
cat(
  "\nControl covariates retained by forward selection:\n"
)

print(
  selected_forward_covariates_60
)

model_controls_forward_60 <- current_model_60

print(
  summary(
    model_controls_forward_60
  )
)

#'------------------------------------------------------------------------------
# STEP 5: assess the selected HFI transition model ----
#'
#' (1) calculate likelihood-based deviance reduction;
#' (2) report model and residual degrees of freedom;
#' (3) calculate marginal and conditional R2.


# 5.1 Final model retained by forward selection ----

model_final_60 <- model_controls_forward_60


# 5.2 Null model fitted on exactly the same dataset ----

model_null_final_data_60 <- lme4::glmer(
  next_terrestrial ~
    behavior_from +
    (1 | individual_id),
  data = forward_data_60,
  family = stats::binomial(link = "logit"),
  control = glmer_control,
  nAGQ = 1
)


# 5.3 Deviance reduction and degrees of freedom ----

deviance_null_60 <- -2 * as.numeric(
  stats::logLik(model_null_final_data_60)
)

deviance_final_60 <- -2 * as.numeric(
  stats::logLik(model_final_60)
)

model_fit_summary_60 <- tibble::tibble(
  deviance_null = deviance_null_60,
  deviance_final = deviance_final_60,
  
  deviance_reduction_percent =
    100 * (
      1 - deviance_final_60 / deviance_null_60
    ),
  
  n_parameters = attr(
    stats::logLik(model_final_60),
    "df"
  ),
  
  residual_df = stats::df.residual(
    model_final_60
  ),
  
  AIC = stats::AIC(
    model_final_60
  )
)

print(model_fit_summary_60)


# 5.4 Marginal and conditional R2 ----
library(performance)
model_R2_60 <- performance::r2_nakagawa(
  model_final_60
)

print(model_R2_60)



library(DHARMa)
#'------------------------------------------------------------------------------
# STEP 6: inspect residual temporal autocorrelation ----

simulation_residuals_60 <- DHARMa::simulateResiduals(
  fittedModel = model_final_60,
  n = 1000
)

plot(simulation_residuals_60)

DHARMa::testUniformity(simulation_residuals_60)
DHARMa::testDispersion(simulation_residuals_60)
DHARMa::testOutliers(simulation_residuals_60)




residual_data_60 <- forward_data_60 %>%
  dplyr::mutate(
    dharma_residual =
      simulation_residuals_60$scaledResiduals
  ) %>%
  dplyr::arrange(
    individual.local.identifier,
    burst_n,
    timestamp
  )


acf_by_burst_60 <- residual_data_60 %>%
  dplyr::group_by(
    individual.local.identifier,
    burst_n
  ) %>%
  dplyr::filter(
    dplyr::n() >= 5
  ) %>%
  dplyr::summarise(
    n_transitions = dplyr::n(),
    
    lag1_autocorrelation = stats::acf(
      dharma_residual,
      lag.max = 1,
      plot = FALSE,
      na.action = na.pass
    )$acf[2],
    
    .groups = "drop"
  )

summary(
  acf_by_burst_60$lag1_autocorrelation
)

longest_bursts_60 <- residual_data_60 %>%
  dplyr::count(
    individual.local.identifier,
    burst_n,
    name = "n"
  ) %>%
  dplyr::arrange(
    dplyr::desc(n)
  ) %>%
  dplyr::slice_head(
    n = 10
  )


# longueur burst
# Calculate lag-1 residual autocorrelation within sufficiently long bursts ----

