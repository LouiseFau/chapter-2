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
  mutate(
    individual.local.identifier = as.character(individual.local.identifier),
    time_UTC = as.POSIXct(gps_timestamp, tz = "UTC"),
    time_local = lubridate::with_tz(time_UTC, tzone = "Europe/Zurich"),
    rf8fitted = as.character(rf8fitted)
  ) |>
  filter(
    !is.na(individual.local.identifier),
    !is.na(time_UTC),
    !is.na(rf8fitted)
  )

# 1.1 Temporal variation (VonBank et al., 2023) ----
# Calcule continuous covariates representing time of day : the variable cos(Diel) 
# represented diurnal (negative values) and nocturnal (positive values) periods,
# and sin(Time) represented midnight until 11:59 am (positive values) and noon 
# until the following 11:59 pm (negative values).

ge_20min <- ge_20min |>
  mutate(
    decimal_hour = lubridate::hour(time_local) +
      lubridate::minute(time_local) / 60 +
      lubridate::second(time_local) / 3600,
    
    cos_Diel = cos(2 * pi * decimal_hour / 24),
    sin_Time = sin(2 * pi * decimal_hour / 24)
  )

#' 1.2 Identify daylight versus night locations ----
tz_loc <- "Europe/Zurich"

ge_20min <- ge_20min |>
  mutate(
    date_local = as.Date(time_local)
  )

sun_times <- suncalc::getSunlightTimes(
  data = data.frame(
    date = ge_20min$date_local,
    lat = ge_20min$location.lat,
    lon = ge_20min$location.long
  ),
  keep = c("sunrise", "sunset"),
  tz = tz_loc
)

ge_20min <- ge_20min |>
  mutate(
    sunrise = sun_times$sunrise,
    sunset = sun_times$sunset,
    is_daylight = time_local >= sunrise & time_local <= sunset,
    is_night = !is_daylight
  )

#' 1.3 Compute distance to nest ----
# Prepare GPS points in EPSG:3035
pts_3035 <- ge_20min |>
  sf::st_as_sf(
    coords = c("location.long", "location.lat"),
    crs = 4326,
    remove = FALSE
  ) |>
  sf::st_transform(3035)

gps_xy_3035 <- sf::st_coordinates(pts_3035)

ge_20min <- ge_20min |>
  mutate(
    x_3035 = gps_xy_3035[, 1],
    y_3035 = gps_xy_3035[, 2]
  )

# Prepare nest sites in EPSG:3035
nest_site_3035 <- nest_site |>
  mutate(
    individual.local.identifier = as.character(individual.local.identifier)
  ) |>
  sf::st_transform(3035)

nest_xy_3035 <- sf::st_coordinates(nest_site_3035)

nest_tbl <- nest_site_3035 |>
  sf::st_drop_geometry() |>
  mutate(
    nest_x_3035 = nest_xy_3035[, 1],
    nest_y_3035 = nest_xy_3035[, 2]
  ) |>
  select(
    individual.local.identifier,
    nest_x_3035,
    nest_y_3035
  ) |>
  distinct(individual.local.identifier, .keep_all = TRUE)

# Join nest coordinates to GPS data and compute distance
ge_20min <- ge_20min |>
  left_join(
    nest_tbl,
    by = "individual.local.identifier"
  ) |>
  mutate(
    distance_to_nest_km = sqrt(
      (x_3035 - nest_x_3035)^2 +
        (y_3035 - nest_y_3035)^2
    ) / 1000
  )


# 1.5 Standardize continuous covariates ----
# We standardize continuous covariates to make model coefficients comparable.
# Temporal covariates are not standardized because they are already cyclic
# variables bounded between -1 and 1.

ge_20min_stand <- ge_20min |>
  mutate(
    age_days_z = as.numeric(scale(age_days)),
    distance_to_nest_km_z = as.numeric(scale(distance_to_nest_km)),
    
    hfi_point_z = as.numeric(scale(hfi_point)),
    hfi_mean_500m_z = as.numeric(scale(hfi_mean_500m)),
    hfi_mean_1000m_z = as.numeric(scale(hfi_mean_1000m))
  )


#' 1.5.1 Population-level behavior counts ----
behavior_counts_population <- ge_20min_stand |>
  count(rf8fitted, name = "n") |>
  mutate(
    prop = n / sum(n),
    prop_percent = 100 * prop
  ) |>
  arrange(desc(n))

print(behavior_counts_population)
#' rf8fitted     n        prop prop_percent
#' Standing 71207 0.651583503   65.1583503
#' Passive  18732 0.171408179   17.1408179
#' Feeding  7081 0.064795073    6.4795073
#' Bodycare 6604 0.060430259    6.0430259
#' Active   3892 0.035613956    3.5613956
#' Walking  1349 0.012344097    1.2344097
#' Undulating 418 0.003824932    0.3824932
#' 
#' We decided to group 'undulating' and 'active' in the category flight
#' We group standing, passive, bodycare under the 'resting' category and then we split 
#' resting into two categories : 'overnight resting' vs 'diurnal resting'


#' 1.5.2 Group behavioral categories ----
ge_20min_stand <- ge_20min_stand |>
  mutate(
    behavior_grouped = case_when(
      rf8fitted %in% c("Active", "Undulating") ~ "flight",
      
      rf8fitted %in% c("Standing", "Passive", "Bodycare") & is_daylight == TRUE ~ "diurnal_resting",
      rf8fitted %in% c("Standing", "Passive", "Bodycare") & is_night == TRUE ~ "overnight_resting",
      
      rf8fitted == "Feeding" ~ "feeding",
      rf8fitted == "Walking" ~ "walking",
      
      rf8fitted %in% c("Standing", "Passive", "Bodycare") & is.na(is_daylight) ~ "resting_unknown_light",
      
      TRUE ~ NA_character_
    ),
    
    behavior_grouped = factor(
      behavior_grouped,
      levels = c(
        "flight",
        "overnight_resting",
        "diurnal_resting",
        "feeding",
        "walking",
        "resting_unknown_light"
      )
    )
  )


# 1.6 Prepare transition table ----
# The transition model requires pairs of consecutive locations within each burst:
# behavior at time t     = behavior_grouped
# behavior at time t + 1 = behavior_grouped_next

state_levels_5 <- c(
  "flight",
  "overnight_resting",
  "diurnal_resting",
  "feeding",
  "walking"
)

ge_20min_stand <- ge_20min_stand |>
  mutate(
    behavior_grouped = factor(
      behavior_grouped,
      levels = state_levels_5
    )
  )

burst_col <- if ("burst_id" %in% names(ge_20min_stand)) "burst_id" else "burst_n"

transitions_20 <- ge_20min_stand |>
  filter(!is.na(behavior_grouped)) |>
  arrange(
    `individual.local.identifier`,
    .data[[burst_col]],
    time_local
  ) |>
  group_by(
    `individual.local.identifier`,
    .data[[burst_col]]
  ) |>
  mutate(
    behavior_grouped_next = dplyr::lead(behavior_grouped),
    time_local_next = dplyr::lead(time_local),
    
    dt_min = as.numeric(
      difftime(time_local_next, time_local, units = "mins")
    ),
    
    transition_type = paste(
      behavior_grouped,
      behavior_grouped_next,
      sep = "_to_"
    )
  ) |>
  ungroup() |>
  filter(!is.na(behavior_grouped_next)) |>
  mutate(
    behavior_grouped = factor(
      behavior_grouped,
      levels = state_levels_5
    ),
    behavior_grouped_next = factor(
      behavior_grouped_next,
      levels = state_levels_5
    ),
    burst_uid = interaction(
      `individual.local.identifier`,
      .data[[burst_col]],
      drop = TRUE
    )
  )

# 1.7 Empirical transition counts ----
transition_count_matrix_20 <- table(
  from = transitions_20$behavior_grouped,
  to = transitions_20$behavior_grouped_next
)

print(transition_count_matrix_20)


# 1.8 Empirical transition probability matrix ----
transition_matrix_20 <- prop.table(
  transition_count_matrix_20,
  margin = 1
)

print(round(transition_matrix_20, 3))
#' we found :
#' from                flight overnight_resting diurnal_resting feeding walking
#' flight             0.079             0.004           0.856   0.042   0.019
#' overnight_resting  0.005             0.879           0.105   0.009   0.003
#' diurnal_resting    0.038             0.003           0.911   0.037   0.011
#' feeding            0.023             0.005           0.483   0.467   0.022
#' walking            0.049             0.006           0.807   0.111   0.026


# Long-format transition table
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

# Rare transitions
rare_transitions_20 <- transition_counts_20 |>
  filter(n > 0, n < 30) |>
  arrange(n)

print(rare_transitions_20)


#'##############################################################################
#' ### Step 2 : model fitting
#'
#' **Philosophy**:
#' We fit one multinomial GAM per HFI spatial scale.
#'
#' The response variable is the arrival behavioral state:
#'   behavior_grouped_next
#'
#' The departure state is:
#'   behavior_grouped
#'
#' The 5 states are:
#'   1. flight
#'   2. overnight_resting
#'   3. diurnal_resting
#'   4. feeding
#'   5. walking
#'
#' We compare:
#'   - a null model without HFI
#'   - a model with HFI under the GPS point
#'   - a model with mean HFI within 500 m
#'   - a model with mean HFI within 1000 m
#'
#' For mgcv::multinom(K = 4), the response must be coded as:
#'   0 = flight
#'   1 = overnight_resting
#'   2 = diurnal_resting
#'   3 = feeding
#'   4 = walking
#'
#' The first state, flight, is the reference arrival state.
#'##############################################################################


# 2.1 Prepare one common dataset for all candidate models ----
# We use the same rows for the null model and the three HFI-scale models.
# This is important because AIC comparisons are only meaningful when models are
# fitted to the same response data.

transition_model_all_20 <- transitions_20 |>
  filter(
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
  mutate(
    behavior_grouped = factor(
      behavior_grouped,
      levels = state_levels_5
    ),
    behavior_grouped_next = factor(
      behavior_grouped_next,
      levels = state_levels_5
    ),
    individual.local.identifier = factor(individual.local.identifier),
    
    arrival_state_id = as.integer(behavior_grouped_next) - 1
  )

# 2.2 Check model dataset ----

model_data_summary_all_20 <- transition_model_all_20 |>
  count(
    behavior_grouped,
    behavior_grouped_next,
    transition_type,
    name = "n"
  ) |>
  group_by(behavior_grouped) |>
  mutate(
    n_from_state = sum(n),
    transition_probability_raw = n / n_from_state
  ) |>
  ungroup() |>
  arrange(behavior_grouped, behavior_grouped_next)

print(model_data_summary_all_20)

arrival_state_check <- transition_model_all_20 |>
  count(behavior_grouped_next, arrival_state_id, name = "n") |>
  arrange(arrival_state_id)

print(arrival_state_check)


# 2.3 Helper function to build multinomial GAM formulas ----
# mgcv::multinom(K = 4) estimates 4 logits relative to the reference category.
# Because we have 5 states, we need a list of 4 formulas.

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
    as.formula(paste("~", rhs)),
    as.formula(paste("~", rhs)),
    as.formula(paste("~", rhs))
  )
}


# 2.4 Null model: no HFI ----
# This model includes departure state, age, distance to nest, time of day,
# and individual identity, but no human footprint index.
m_null_20 <- mgcv::gam(
  formula = make_multinom_formula(hfi_var = NULL),
  family = mgcv::multinom(K = 4),
  data = transition_model_all_20,
  method = "ML"
)

summary(m_null_20)


# 2.5 Model 1: HFI under GPS point ----
m_hfi_point_all_20 <- mgcv::gam(
  formula = make_multinom_formula(hfi_var = "hfi_point_z"),
  family = mgcv::multinom(K = 4),
  data = transition_model_all_20,
  method = "ML"
)

summary(m_hfi_point_all_20)


# 2.6 Model 2: mean HFI within 500 m radius ----
m_hfi_500m_all_20 <- mgcv::gam(
  formula = make_multinom_formula(hfi_var = "hfi_mean_500m_z"),
  family = mgcv::multinom(K = 4),
  data = transition_model_all_20,
  method = "ML"
)

summary(m_hfi_500m_all_20)


# 2.7 Model 3: mean HFI within 1000 m radius ----
m_hfi_1000m_all_20 <- mgcv::gam(
  formula = make_multinom_formula(hfi_var = "hfi_mean_1000m_z"),
  family = mgcv::multinom(K = 4),
  data = transition_model_all_20,
  method = "ML"
)

summary(m_hfi_1000m_all_20)


# 2.8 Compare candidate models with AIC ----
# The null model tells us whether HFI improves the transition model.
# The three HFI models tell us which spatial scale is best supported.

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
  mutate(
    delta_AIC = AIC - min(AIC),
    AIC_weight = exp(-0.5 * delta_AIC) / sum(exp(-0.5 * delta_AIC))
  ) |>
  arrange(AIC) |>
  select(model, df, AIC, delta_AIC, AIC_weight)

print(model_comparison_all_20)
#' we find that the best model is at 1000m resolution
#' model                 df     AIC       delta_AIC   AIC_weight
#' m_hfi_1000m_all_20  279.7587 75737.38  0.00000  9.961537e-01
#' m_hfi_500m_all_20   279.5177 75749.79  12.40755 2.014005e-03
#' m_hfi_point_all_20  279.3300 75749.98  12.59664 1.832318e-03
#' m_null_20           260.5271 75772.64  35.26032 2.196064e-08

#'##############################################################################
#' ### Step 3 : plot the results
#'
#' For the 5-state multinomial transition model, coefficient-level plots are hard
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


# 3.1 Define model list and HFI variables ----

state_levels_5 <- c(
  "flight",
  "overnight_resting",
  "diurnal_resting",
  "feeding",
  "walking"
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


# 3.2 Helper function: predict transition probabilities ----

predict_transition_probabilities <- function(model,
                                             hfi_var,
                                             hfi_values,
                                             panel_name,
                                             model_data = transition_model_all_20) {
  
  newdata <- expand.grid(
    behavior_grouped = state_levels_5,
    hfi_z = hfi_values,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  
  newdata <- newdata |>
    mutate(
      behavior_grouped = factor(
        behavior_grouped,
        levels = state_levels_5
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
  
  if (ncol(pred_prob) != length(state_levels_5)) {
    stop(
      "The prediction output has ",
      ncol(pred_prob),
      " columns, but ",
      length(state_levels_5),
      " behavioral states were expected."
    )
  }
  
  names(pred_prob) <- state_levels_5
  
  pred_long <- bind_cols(
    newdata |>
      select(behavior_grouped, hfi_z),
    pred_prob
  ) |>
    tidyr::pivot_longer(
      cols = all_of(state_levels_5),
      names_to = "to_state",
      values_to = "transition_probability"
    ) |>
    mutate(
      panel = panel_name,
      from_state = behavior_grouped,
      to_state = factor(
        to_state,
        levels = state_levels_5
      ),
      from_state = factor(
        from_state,
        levels = state_levels_5
      )
    )
  
  pred_long
}


# 3.3 Helper function: extract high-minus-low HFI effect ----

extract_hfi_delta_probability <- function(model,
                                          hfi_var,
                                          panel_name,
                                          model_data = transition_model_all_20) {
  
  if (is.na(hfi_var)) {
    
    delta_null <- expand.grid(
      from_state = state_levels_5,
      to_state = state_levels_5,
      KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE
    ) |>
      mutate(
        panel = panel_name,
        low_HFI = NA_real_,
        high_HFI = NA_real_,
        low_probability = NA_real_,
        high_probability = NA_real_,
        delta_probability = 0
      )
    
    return(delta_null)
  }
  
  hfi_low_high <- quantile(
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
    mutate(
      hfi_level = case_when(
        hfi_z == as.numeric(hfi_low_high[1]) ~ "low_probability",
        hfi_z == as.numeric(hfi_low_high[2]) ~ "high_probability",
        TRUE ~ NA_character_
      )
    ) |>
    select(
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
    mutate(
      low_HFI = as.numeric(hfi_low_high[1]),
      high_HFI = as.numeric(hfi_low_high[2]),
      delta_probability = high_probability - low_probability
    )
  
  pred
}


# 3.4 Extract HFI effects for all candidate models ----

transition_delta_all_20 <- purrr::map2_dfr(
  model_info_20$model_name,
  seq_len(nrow(model_info_20)),
  function(model_name, i) {
    
    extract_hfi_delta_probability(
      model = model_objects_20[[model_name]],
      hfi_var = model_info_20$hfi_var[i],
      panel_name = model_info_20$panel[i],
      model_data = transition_model_all_20
    )
  }
) |>
  mutate(
    panel = factor(
      panel,
      levels = model_info_20$panel
    ),
    from_state = factor(
      from_state,
      levels = state_levels_5
    ),
    to_state = factor(
      to_state,
      levels = state_levels_5
    ),
    label = sprintf("%.2f", delta_probability)
  )

print(transition_delta_all_20)


# 3.5 Plot HFI effect on transition probabilities ----

max_abs_delta_20 <- max(
  abs(transition_delta_all_20$delta_probability),
  na.rm = TRUE
)

p_hfi_transition_delta_20 <- ggplot(
  transition_delta_all_20,
  aes(
    x = from_state,
    y = to_state,
    fill = delta_probability
  )
) +
  geom_tile(color = "grey70", linewidth = 0.5) +
  geom_text(aes(label = label), size = 3) +
  facet_wrap(~ panel, nrow = 1) +
  scale_fill_gradient2(
    low = "#3B5BDB",
    mid = "grey95",
    high = "#FF6B4A",
    midpoint = 0,
    limits = c(-max_abs_delta_20, max_abs_delta_20),
    name = expression(Delta*" transition probability")
  ) +
  labs(
    x = "Departure state",
    y = "Arrival state",
    title = "Change in predicted transition probabilities along the HFI gradient",
    subtitle = "P(transition | HFI élevé) - P(transition | HFI faible)"
  ) +
  coord_equal() +
  theme_bw(base_size = 13) +
  theme(
    strip.background = element_rect(fill = "white", color = "white"),
    strip.text = element_text(face = "bold", size = 12),
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10)
  )

print(p_hfi_transition_delta_20)





# 3.7 Save figures ----
output_dir <- "C:/Users/lfaure7/OneDrive/THESE/CHAPITRE 2/git/chapter-2/HMM/HMM on ACC-classified behaviors/donnes filtree intermediaire"

ggsave(
  filename = file.path(
    output_dir,
    "transition_hfi_delta_probability_20min_5state.png"
  ),
  plot = p_hfi_transition_delta_20,
  width = 15,
  height = 5,
  dpi = 300
)
