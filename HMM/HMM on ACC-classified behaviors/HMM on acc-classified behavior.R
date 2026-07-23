#' ----------------------------------------------------------------------------- 
#' Title: Discrete Markov Model with multiple behavioral states and covariates 
#' Authors : Louise Faure
#' Date : 07.07.26
#' Purpose : 
#' (1) Identify transition probabilities between behavioral states 
#' (2) Model fitting
#' (3) Plot the results 
#' Data : this code is run for ACC classified dataset. 
#' Hypothesis : we assume that when eagles are in highly human dominated 
#' environment, the probability of keep flying is higher, and the probability of
#' transitioning to a resting, feeding, walking and peering behavior is reduced. 
#' -----------------------------------------------------------------------------



library(dplyr)
library(tidyr)
library(lubridate)
library(stringr)
library(tibble)
library(mgcv)
library(ggplot2)


ge_20min <- readRDS("/Users/louisefaure/Library/CloudStorage/OneDrive-Personnel/THESE/CHAPITRE 2/git/chapter-2/HMM/HMM on ACC-classified behaviors/donnees intermediaire (2)/GE_20_min_thinned_behavior_assigned_hfi.rds")
emig_dates <- readRDS("/Users/louisefaure/Library/CloudStorage/OneDrive-Personnel/THESE/CHAPITRE 2/git/chapter-2/DONNEES AIGLES/emigration dates/emigration_dates_20250417.rds")
nest_site <- readRDS("/Users/louisefaure/Library/CloudStorage/OneDrive-Personnel/THESE/CHAPITRE 2/git/chapter-2/DONNEES AIGLES/nest site location/nest_site_location/nest_site_location.rds")



#'------------------------------------------------------------------------------
#' ### Step 1 : Identification of the probabilities of transition
#' **Philosophy**: we focus on the transition from flight to terrestrial and 
#' from flight to flight. 
#' 
#' **Covariates**: we consider as covariates the human footprint index, the age 
#' since emigration, the distance to nest and the time of the day. 
#' 
#' **Steps**:
#' (1) compute temporal covariate and distance to nest
#' (2) standardize continuous variables
#' (3) inspect the number of behavior observed per behavioral categories at the
#' population level (and eventually group certain behavioral categories)

ge_20min <- ge_20min |>
  dplyr::mutate(
    individual.local.identifier = as.character(individual.local.identifier),
    time_UTC = as.POSIXct(ge_20min$gps_timestamp, tz = "UTC"),
    time_local = lubridate::with_tz(time_UTC, tzone = "Europe/Zurich"),
    behavior_state = as.character(behavior_reclassified)
  ) |>
  dplyr::filter(
    !is.na(individual.local.identifier),
    !is.na(time_UTC),
    !is.na(behavior_state))

# 1.1 Temporal variation (VonBank et al., 2023) ----
# Calcule continuous covariates representing time of day : the variable cos(Diel) 
# represented diurnal (negative values) and nocturnal (positive values) periods,
# and sin(Time) represented midnight until 11:59 am (positive values) and noon 
# until the following 11:59 pm (negative values).

ge_20min <- ge_20min |>
  dplyr::mutate(
    decimal_hour = lubridate::hour(time_local) +
      lubridate::minute(time_local) / 60 +
      lubridate::second(time_local) / 3600,
    
    cos_Diel = cos(2 * pi * decimal_hour / 24),
    sin_Time = sin(2 * pi * decimal_hour / 24),
    date_local = as.Date(time_local))

#' 1.2 Compute distance to nest ----
# Remove old nest-distance columns if this block has already been run
ge_20min <- ge_20min |>
  dplyr::select(
    -dplyr::any_of(
      c(
        "x_3035",
        "y_3035",
        "nest_x_3035",
        "nest_y_3035",
        "nest_x_3035.x",
        "nest_y_3035.x",
        "nest_x_3035.y",
        "nest_y_3035.y",
        "distance_to_nest_km",
        "distance_to_nest_km.x",
        "distance_to_nest_km.y")))

# GPS points in EPSG:3035
pts_3035 <- ge_20min |>
  sf::st_as_sf(
    coords = c("location.long", "location.lat"),
    crs = 4326,
    remove = FALSE
  ) |>
  sf::st_transform(3035)

gps_xy_3035 <- sf::st_coordinates(pts_3035)

ge_20min <- ge_20min |>
  dplyr::mutate(
    x_3035 = gps_xy_3035[, 1],
    y_3035 = gps_xy_3035[, 2]
  )

# Nest sites in EPSG:3035
nest_site_3035 <- nest_site |>
  dplyr::mutate(
    individual.local.identifier = as.character(individual.local.identifier)
  ) |>
  sf::st_transform(3035)

nest_xy_3035 <- sf::st_coordinates(nest_site_3035)

nest_tbl <- nest_site_3035 |>
  sf::st_drop_geometry() |>
  dplyr::mutate(
    nest_x_3035 = nest_xy_3035[, 1],
    nest_y_3035 = nest_xy_3035[, 2]
  ) |>
  dplyr::select(
    individual.local.identifier,
    nest_x_3035,
    nest_y_3035
  ) |>
  dplyr::distinct(
    individual.local.identifier,
    .keep_all = TRUE
  )

# Join nest coordinates and compute distance
ge_20min <- ge_20min |>
  dplyr::left_join(
    nest_tbl,
    by = "individual.local.identifier"
  ) |>
  dplyr::mutate(
    distance_to_nest_km = sqrt(
      (x_3035 - nest_x_3035)^2 +
        (y_3035 - nest_y_3035)^2
    ) / 1000
  )

# 1.3 Standardize continuous covariates ----
# We standardize continuous covariates to make model coefficients comparable.
# Temporal covariates are not standardized because they are already cyclic
# variables bounded between -1 and 1.
ge_20min_stand <- ge_20min |>
  dplyr::mutate(
    age_days_z = as.numeric(scale(age_days)),
    distance_to_nest_km_z = as.numeric(scale(distance_to_nest_km)),
    
    hfi_point_z = as.numeric(scale(hfi_point)),
    hfi_mean_500m_z = as.numeric(scale(hfi_mean_500m)),
    hfi_mean_1000m_z = as.numeric(scale(hfi_mean_1000m))
  )

#' 1.4 Population-level behavior counts ----
behavior_counts_population <- ge_20min_stand |>
  count(behavior_reclassified, name = "n") |>
  mutate(
    prop = n / sum(n),
    prop_percent = 100 * prop
  ) |>
  arrange(desc(n))

print(behavior_counts_population)
#' behavior    n        prop prop_percent
#' resting 74868 0.72122999    72.122999
#' flight 21886 0.21083560    21.083560
#' feeding  7052 0.06793442     6.793442


# 1.5 Prepare transition table ----
# The transition model requires pairs of consecutive locations within each burst:
# behavior at time t     = behavior_grouped
# behavior at time t + 1 = behavior_grouped_next
state_levels_3 <- c("flight","resting","feeding")

transitions_20 <- ge_20min_stand |>
  dplyr::mutate(
    behavior_grouped = factor(
      behavior_reclassified,
      levels = state_levels_3)
  ) |>
  dplyr::arrange(
    individual.local.identifier,
    burst_id,
    time_UTC
  ) |>
  dplyr::group_by(
    individual.local.identifier,
    burst_id
  ) |>
  dplyr::mutate(
    behavior_grouped_next = dplyr::lead(behavior_grouped),
    time_UTC_next = dplyr::lead(time_UTC),
    gps_row_id_next = dplyr::lead(gps_row_id),
    
    transition_dt_min = as.numeric(
      difftime(time_UTC_next, time_UTC, units = "mins")
    ),
    
    transition_type = paste(
      behavior_grouped,
      behavior_grouped_next,
      sep = "_to_"
    )
  ) |>
  dplyr::ungroup() |>
  dplyr::filter(
    has_next_regular == TRUE,
    !is.na(behavior_grouped),
    !is.na(behavior_grouped_next))


# 1.6 Empirical transition counts ----
transition_count_matrix_20 <- table(
  from = transitions_20$behavior_grouped,
  to = transitions_20$behavior_grouped_next
)

print(transition_count_matrix_20)
#'         to
#'from      flight resting feeding
#'flight    9728    9531     767
#'resting   9020   55925    2601
#'feeding    628    2619    3095

# 1.7 Empirical transition probability matrix ----
row_totals_20 <- rowSums(transition_count_matrix_20)

transition_matrix_20 <- sweep(
  transition_count_matrix_20,
  MARGIN = 1,
  STATS = row_totals_20,
  FUN = "/"
)

transition_matrix_20[row_totals_20 == 0, ] <- NA

print(round(transition_matrix_20, 3))
#'            to
#'  from      flight resting feeding
#'  flight   0.486   0.476   0.038
#'  resting  0.134   0.828   0.039
#'  feeding  0.099   0.413   0.488


# Long-format transition table to inspect rare transitions
transition_counts_20 <- as.data.frame(
  transition_count_matrix_20,
  stringsAsFactors = FALSE
) |>
  rename(
    behavior_grouped = from,
    behavior_grouped_next = to,
    n = Freq
  ) |>
  group_by(behavior_grouped) |>
  mutate(
    n_from_state = sum(n),
    transition_probability = n / n_from_state,
    transition_type = paste(
      behavior_grouped,
      behavior_grouped_next,
      sep = "_to_"
    )
  ) |>
  ungroup() |>
  arrange(behavior_grouped, behavior_grouped_next)

print(transition_counts_20)

rare_transitions_20 <- transition_counts_20 |>
  filter(n > 0, n < 30) |>
  arrange(n)

print(rare_transitions_20)


#'##############################################################################
#' ### Step 3: model fitting
#'
#' We fit one multinomial GAM per HFI spatial scale.
#'
#' Response variable:
#'   behavior_grouped_next
#'
#' Departure state:
#'   behavior_grouped
#'
#' The 3 states are:
#'   1. flight
#'   2. resting
#'   3. feeding
#'
#' For mgcv::multinom(K = 2), the response must be coded as:
#'   0 = flight
#'   1 = resting
#'   2 = feeding
#'
#' The first state, flight, is the reference arrival state.
#'##############################################################################

state_levels_3 <- c(
  "flight",
  "resting",
  "feeding"
)


# 3.1 Prepare one common dataset for all candidate models ----

transition_model_all_20 <- transitions_20 |>
  dplyr::filter(
    !is.na(behavior_grouped),
    !is.na(behavior_grouped_next),
    !is.na(hfi_point_z),
    !is.na(hfi_mean_500m_z),
    !is.na(hfi_mean_1000m_z),
    !is.na(age_days_z),
    !is.na(distance_to_nest_km_z),
    !is.na(cos_Diel),
    !is.na(sin_Time),
    !is.na(individual.local.identifier)
  ) |>
  dplyr::mutate(
    behavior_grouped = factor(
      behavior_grouped,
      levels = state_levels_3
    ),
    behavior_grouped_next = factor(
      behavior_grouped_next,
      levels = state_levels_3
    ),
    individual.local.identifier = factor(individual.local.identifier),
    
    arrival_state_id = as.integer(behavior_grouped_next) - 1
  )


# 3.2 Check model dataset ----

model_data_summary_all_20 <- transition_model_all_20 |>
  dplyr::count(
    behavior_grouped,
    behavior_grouped_next,
    transition_type,
    name = "n"
  ) |>
  dplyr::group_by(behavior_grouped) |>
  dplyr::mutate(
    n_from_state = sum(n),
    transition_probability_raw = n / n_from_state
  ) |>
  dplyr::ungroup() |>
  dplyr::arrange(behavior_grouped, behavior_grouped_next)

print(model_data_summary_all_20)

arrival_state_check <- transition_model_all_20 |>
  dplyr::count(behavior_grouped_next, arrival_state_id, name = "n") |>
  dplyr::arrange(arrival_state_id)

print(arrival_state_check)


# 3.3 Helper function to build multinomial GAM formulas ----
# mgcv::multinom(K = 2) estimates 2 logits relative to the reference category.
# Because we have 3 states, we need a list of 2 formulas.

make_multinom_formula <- function(hfi_var = NULL) {
  
  base_covariates <- c(
    "age_days_z",
    "distance_to_nest_km_z",
    "cos_Diel",
    "sin_Time"
  )
  
  if (is.null(hfi_var)) {
    covariates <- base_covariates
  } else {
    covariates <- c(hfi_var, base_covariates)
  }
  
  rhs <- paste0(
    "behavior_grouped * (",
    paste(covariates, collapse = " + "),
    ") + s(individual.local.identifier, bs = 're')"
  )
  
  list(
    as.formula(paste("arrival_state_id ~", rhs)),
    as.formula(paste("~", rhs))
  )
}


# 3.4 Null model: no HFI ----
m_null_20 <- mgcv::gam(
  formula = make_multinom_formula(hfi_var = NULL),
  family = mgcv::multinom(K = 2),
  data = transition_model_all_20,
  method = "ML"
)

summary(m_null_20)


# 3.5 Model 1: HFI under GPS point ----
m_hfi_point_all_20 <- mgcv::gam(
  formula = make_multinom_formula(hfi_var = "hfi_point_z"),
  family = mgcv::multinom(K = 2),
  data = transition_model_all_20,
  method = "ML"
)

summary(m_hfi_point_all_20)

# 3.6 Model 2: mean HFI within 500 m radius ----
m_hfi_500m_all_20 <- mgcv::gam(
  formula = make_multinom_formula(hfi_var = "hfi_mean_500m_z"),
  family = mgcv::multinom(K = 2),
  data = transition_model_all_20,
  method = "ML"
)

summary(m_hfi_500m_all_20)

# 3.7 Model 3: mean HFI within 1000 m radius ----
m_hfi_1000m_all_20 <- mgcv::gam(
  formula = make_multinom_formula(hfi_var = "hfi_mean_1000m_z"),
  family = mgcv::multinom(K = 2),
  data = transition_model_all_20,
  method = "ML"
)

summary(m_hfi_1000m_all_20)


# 3.8 Compare candidate models with AIC ----
model_comparison_all_20 <- AIC(
  m_null_20,
  m_hfi_point_all_20,
  m_hfi_500m_all_20,
  m_hfi_1000m_all_20
) |>
  as.data.frame()

model_comparison_all_20$model <- rownames(model_comparison_all_20)
rownames(model_comparison_all_20) <- NULL

model_comparison_all_20 <- model_comparison_all_20 |>
  dplyr::mutate(
    delta_AIC = AIC - min(AIC),
    AIC_weight = exp(-0.5 * delta_AIC) / sum(exp(-0.5 * delta_AIC))
  ) |>
  dplyr::arrange(AIC) |>
  dplyr::select(model, df, AIC, delta_AIC, AIC_weight)

print(model_comparison_all_20)
#'              model       df      AIC delta_AIC   AIC_weight
#'m_hfi_1000m_all_20 129.7369 77493.23   0.00000 1.000000e+00
#' m_hfi_500m_all_20 129.8399 77527.37  34.14160 3.856962e-08
#' m_hfi_point_all_20 129.8588 77556.72  63.49247 1.632246e-14
#' m_null_20 123.8969 77725.61 232.37977 3.462379e-51



#'##############################################################################
#' ### Step 4: plot the results
#'
#' For the 3-state multinomial transition model, coefficient-level plots are hard
#' to interpret because mgcv::multinom() estimates log-odds relative to a reference
#' arrival state.
#'
#' We therefore plot the predicted effect of HFI as:
#'
#'   delta probability = P(to_state | high HFI, from_state)
#'                     - P(to_state | low HFI, from_state)
#'
#' For each HFI scale, low and high HFI correspond to the 5% and 95% quantiles of
#' the standardized HFI covariate used in that model.
#'##############################################################################


# 4.1 Define model list and HFI variables ----

state_levels_3 <- c(
  "flight",
  "resting",
  "feeding"
)

model_info_20 <- tibble::tibble(
  model_name = c(
    "m_null_20",
    "m_hfi_point_all_20",
    "m_hfi_500m_all_20",
    "m_hfi_1000m_all_20"
  ),
  panel = c(
    "Null model",
    "HFI point",
    "HFI 500m",
    "HFI 1000m"
  ),
  hfi_var = c(
    NA_character_,
    "hfi_point_z",
    "hfi_mean_500m_z",
    "hfi_mean_1000m_z"
  )
)

model_objects_20 <- list(
  m_null_20 = m_null_20,
  m_hfi_point_all_20 = m_hfi_point_all_20,
  m_hfi_500m_all_20 = m_hfi_500m_all_20,
  m_hfi_1000m_all_20 = m_hfi_1000m_all_20
)


# 4.2 Helper function: predict transition probabilities ----
predict_transition_probabilities <- function(model,
                                             hfi_var,
                                             hfi_values,
                                             panel_name,
                                             model_data = transition_model_all_20) {
  
  newdata <- expand.grid(
    behavior_grouped = state_levels_3,
    hfi_z = hfi_values,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  
  newdata <- newdata |>
    dplyr::mutate(
      behavior_grouped = factor(
        behavior_grouped,
        levels = state_levels_3
      ),
      age_days_z = 0,
      distance_to_nest_km_z = 0,
      cos_Diel = 0,
      sin_Time = 0,
      individual.local.identifier = factor(
        levels(model_data$individual.local.identifier)[1],
        levels = levels(model_data$individual.local.identifier)
      ),
      
      hfi_point_z = 0,
      hfi_mean_500m_z = 0,
      hfi_mean_1000m_z = 0
    )
  
  if (!is.na(hfi_var)) {
    newdata[[hfi_var]] <- newdata$hfi_z
  }
  
  pred_prob <- tryCatch(
    {
      predict(
        model,
        newdata = newdata,
        type = "response",
        exclude = "s(individual.local.identifier)"
      )
    },
    error = function(e) {
      predict(
        model,
        newdata = newdata,
        type = "response"
      )
    }
  )
  
  pred_prob <- as.data.frame(pred_prob)
  
  if (ncol(pred_prob) != length(state_levels_3)) {
    stop(
      "The prediction output has ",
      ncol(pred_prob),
      " columns, but ",
      length(state_levels_3),
      " behavioral states were expected."
    )
  }
  
  names(pred_prob) <- state_levels_3
  
  pred_long <- dplyr::bind_cols(
    newdata |>
      dplyr::select(behavior_grouped, hfi_z),
    pred_prob
  ) |>
    tidyr::pivot_longer(
      cols = dplyr::all_of(state_levels_3),
      names_to = "to_state",
      values_to = "transition_probability"
    ) |>
    dplyr::mutate(
      panel = panel_name,
      from_state = behavior_grouped,
      to_state = factor(
        to_state,
        levels = state_levels_3
      ),
      from_state = factor(
        from_state,
        levels = state_levels_3
      )
    )
  
  pred_long
}


# 4.3 Helper function: extract high-minus-low HFI effect ----

extract_hfi_delta_probability <- function(model,
                                          hfi_var,
                                          panel_name,
                                          model_data = transition_model_all_20) {
  
  if (is.na(hfi_var)) {
    
    delta_null <- expand.grid(
      from_state = state_levels_3,
      to_state = state_levels_3,
      KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE
    ) |>
      dplyr::mutate(
        panel = panel_name,
        low_HFI = NA_real_,
        high_HFI = NA_real_,
        low_probability = NA_real_,
        high_probability = NA_real_,
        delta_probability = 0
      )
    
    return(delta_null)
  }
  
  hfi_low_high <- stats::quantile(
    model_data[[hfi_var]],
    probs = c(0.05, 0.95),
    na.rm = TRUE
  )
  
  pred <- predict_transition_probabilities(
    model = model,
    hfi_var = hfi_var,
    hfi_values = as.numeric(hfi_low_high),
    panel_name = panel_name,
    model_data = model_data
  ) |>
    dplyr::mutate(
      hfi_level = dplyr::case_when(
        hfi_z == as.numeric(hfi_low_high[1]) ~ "low_probability",
        hfi_z == as.numeric(hfi_low_high[2]) ~ "high_probability",
        TRUE ~ NA_character_
      )
    ) |>
    dplyr::select(
      panel,
      from_state,
      to_state,
      hfi_level,
      transition_probability
    ) |>
    tidyr::pivot_wider(
      names_from = hfi_level,
      values_from = transition_probability
    ) |>
    dplyr::mutate(
      low_HFI = as.numeric(hfi_low_high[1]),
      high_HFI = as.numeric(hfi_low_high[2]),
      delta_probability = high_probability - low_probability
    )
  
  pred
}


# 4.4 Extract HFI effects for all candidate models ----

transition_delta_list_20 <- lapply(
  seq_len(nrow(model_info_20)),
  function(i) {
    
    model_name_i <- model_info_20$model_name[i]
    
    extract_hfi_delta_probability(
      model = model_objects_20[[model_name_i]],
      hfi_var = model_info_20$hfi_var[i],
      panel_name = model_info_20$panel[i],
      model_data = transition_model_all_20
    )
  }
)

transition_delta_all_20 <- dplyr::bind_rows(transition_delta_list_20) |>
  dplyr::mutate(
    panel = factor(
      panel,
      levels = model_info_20$panel
    ),
    from_state = factor(
      from_state,
      levels = state_levels_3
    ),
    to_state = factor(
      to_state,
      levels = state_levels_3
    ),
    label = sprintf("%.2f", delta_probability)
  )

print(transition_delta_all_20)


# 4.5 Plot HFI effect on transition probabilities ----

max_abs_delta_20 <- max(
  abs(transition_delta_all_20$delta_probability),
  na.rm = TRUE
)

p_hfi_transition_delta_20 <- ggplot2::ggplot(
  transition_delta_all_20,
  ggplot2::aes(
    x = from_state,
    y = to_state,
    fill = delta_probability
  )
) +
  ggplot2::geom_tile(color = "grey70", linewidth = 0.5) +
  ggplot2::geom_text(ggplot2::aes(label = label), size = 3) +
  ggplot2::facet_wrap(~ panel, nrow = 1) +
  ggplot2::scale_fill_gradient2(
    low = "#3B5BDB",
    mid = "grey95",
    high = "#FF6B4A",
    midpoint = 0,
    limits = c(-max_abs_delta_20, max_abs_delta_20),
    name = expression(Delta*" transition probability")
  ) +
  ggplot2::labs(
    x = "Departure state",
    y = "Arrival state",
    title = "Change in predicted transition probabilities along the HFI gradient",
    subtitle = "P(transition | high HFI) - P(transition | low HFI)"
  ) +
  ggplot2::coord_equal() +
  ggplot2::theme_bw(base_size = 13) +
  ggplot2::theme(
    strip.background = ggplot2::element_rect(fill = "white", color = "white"),
    strip.text = ggplot2::element_text(face = "bold", size = 12),
    panel.grid = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
    legend.title = ggplot2::element_text(size = 11),
    legend.text = ggplot2::element_text(size = 10)
  )

print(p_hfi_transition_delta_20)


# 4.6 Save figure ----

ggplot2::ggsave(
  filename = file.path(
    "/Users/louisefaure/Library/CloudStorage/OneDrive-Personnel/THESE/CHAPITRE 2/git/chapter-2/HMM/HMM on ACC-classified behaviors/donnees intermediaire (2)/transition_hfi_delta_probability_60min_3state.png"
  ),
  plot = p_hfi_transition_delta_20,
  width = 11,
  height = 4.5,
  dpi = 300
)
