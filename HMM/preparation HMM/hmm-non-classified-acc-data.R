#' ############################################################################# 
#' Title: Discrete Markov Model with two behavioral states and covariates 
#' Authors : Louise Faure
#' Date : 25.06.26
#' Purpose : 
#' (1) Identify transition probabilities 
#' (2) Fit model focusing on the gps point of departure only
#' (3) Plot the results 
#' Data : this code is run for the non-ACC classified dataset. 
#' Hypothesis : we assume that when eagles are in highly human dominated 
#' environment, the probability of keep flying is higher, and the probability of
#' transitioning to a terrestrial behavior is reduced. 
#' #############################################################################



# 0. Setup ----
# 0.1 library 
library(dplyr)
library(tidyr)
library(lubridate)
library(stringr)
library(tibble)
library(mgcv)
library(ggplot2)

# 0.2 dataset
ge_20min <- readRDS(
  "C:/Users/lfaure7/Documents/git/chapter-2/HMM/preparation HMM/donnees intermediaire/GE_20_min_thinned_hfi.rds"
)


#'##############################################################################
#' ### Step 1 : Identification of the probabilities of transition
#' **Philosophy**: we focus on the transition from flight to terrestrial and 
#' from flight to flight. 
#' 
#' **Covariates**: we consider as covariates the human footprint index, the age 
#' since emigration, the distance to nest and the time of the day. 
#'##############################################################################

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


# 1.2 Standardize continuous covariates ----
# We standardize continuous covariates to make model coefficients comparable.
# Temporal covariates are not standardized because they are already cyclic
# variables bounded between -1 and 1.

ge_20min_stand <- ge_20min |>
  mutate(
    days_since_emig_z = as.numeric(scale(days_since_emig)),
    distance_to_nest_km_z = as.numeric(scale(distance_to_nest_km)),
    
    hfi_point_z = as.numeric(scale(hfi_point)),
    hfi_mean_500m_z = as.numeric(scale(hfi_mean_500m)),
    hfi_mean_1000m_z = as.numeric(scale(hfi_mean_1000m))
  )

# 1.3 Prepare transition table ----
# The transition model requires pairs of consecutive locations within each burst:
# behavior at time t     = behavior_broad
# behavior at time t + 1 = behavior_broad_next

burst_col <- if ("burst_id" %in% names(ge_20min_stand)) "burst_id" else "burst_n"

transitions_20 <- ge_20min_stand |>
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
    behavior_broad_next = dplyr::lead(behavior_broad),
    time_local_next = dplyr::lead(time_local),
    
    dt_min = as.numeric(
      difftime(time_local_next, time_local, units = "mins")
    ),
    
    transition_type = paste(
      behavior_broad,
      behavior_broad_next,
      sep = "_to_"
    ),
    
    next_terrestrial = case_when(
      behavior_broad_next == "terrestrial" ~ 1L,
      behavior_broad_next == "flight" ~ 0L,
      TRUE ~ NA_integer_
    )
  ) |>
  ungroup() |>
  filter(!is.na(behavior_broad_next)) |>
  filter(!is.na(next_terrestrial)) |>
  mutate(
    burst_uid = interaction(
      `individual.local.identifier`,
      .data[[burst_col]],
      drop = TRUE
    )
  )

# 1.4 Empirical transition probabilities ----

transition_counts_20 <- transitions_20 |>
  count(
    behavior_broad,
    behavior_broad_next,
    transition_type,
    name = "n"
  ) |>
  group_by(behavior_broad) |>
  mutate(
    n_from_state = sum(n),
    transition_probability = n / n_from_state
  ) |>
  ungroup() |>
  arrange(behavior_broad, behavior_broad_next)

print(transition_counts_20)

# 1.5 Transition probability matrix ----

transition_matrix_20 <- transition_counts_20 |>
  dplyr::select(
    behavior_broad,
    behavior_broad_next,
    transition_probability
  ) |>
  tidyr::pivot_wider(
    names_from = behavior_broad_next,
    values_from = transition_probability,
    values_fill = 0
  ) |>
  tibble::column_to_rownames("behavior_broad") |>
  as.matrix()

print(transition_matrix_20)

#' we found :
#'              flight    terrestrial
#' flight      0.5064611   0.4935389
#' terrestrial 0.1620199   0.8379801

#'##############################################################################
#' ### Step 2 : model fitting
#' **Philosophy**: we fit one combined binomial GAMM per HFI spatial scale.
#' 
#' The response variable is `next_terrestrial`:
#'   1 = arrival state is terrestrial
#'   0 = arrival state is flight
#'
#' The departure state is given by `behavior_broad`.
#' Therefore:
#'   if behavior_broad == "flight":
#'      next_terrestrial = 1 means flight -> terrestrial
#'      next_terrestrial = 0 means flight -> flight
#'
#'   if behavior_broad == "terrestrial":
#'      next_terrestrial = 1 means terrestrial -> terrestrial
#'      next_terrestrial = 0 means terrestrial -> flight
#'
#' We include interactions between departure state and covariates so that
#' covariate effects can differ depending on whether the animal starts from
#' flight or terrestrial behavior.
#'
#' Individual identity is included as a random intercept using:
#'   s(individual.local.identifier, bs = "re")
#'##############################################################################


# 2.1 Prepare one common dataset for all candidate models ----
# We use the same rows for the null model and the three HFI-scale models.
# This is important because AIC comparisons are only meaningful when models are
# fitted to the same response data.

transition_model_all_20 <- transitions_20 |>
  drop_na(
    next_terrestrial,
    behavior_broad,
    hfi_point_z,
    hfi_mean_500m_z,
    hfi_mean_1000m_z,
    days_since_emig_z,
    distance_to_nest_km_z,
    cos_Diel,
    sin_Time,
    `individual.local.identifier`
  ) |>
  mutate(
    behavior_broad = factor(
      behavior_broad,
      levels = c("flight", "terrestrial")
    ),
    `individual.local.identifier` = factor(`individual.local.identifier`)
  )


# 2.2 Check model dataset ----

model_data_summary_all_20 <- transition_model_all_20 |>
  count(
    behavior_broad,
    behavior_broad_next,
    transition_type,
    name = "n"
  ) |>
  group_by(behavior_broad) |>
  mutate(
    n_from_state = sum(n),
    transition_probability_raw = n / n_from_state
  ) |>
  ungroup()

print(model_data_summary_all_20)


# 2.3 Null model: no HFI ----
# This model tests the transition process with age, distance to nest, time of day,
# and individual identity, but without human footprint.

m_null_20 <- mgcv::bam(
  next_terrestrial ~
    behavior_broad * (
      days_since_emig_z +
        distance_to_nest_km_z +
        cos_Diel +
        sin_Time
    ) +
    s(`individual.local.identifier`, bs = "re"),
  
  family = binomial(link = "logit"),
  data = transition_model_all_20,
  method = "fREML",
  discrete = TRUE
)

summary(m_null_20)


# 2.4 Model 1: HFI under GPS point ----

m_hfi_point_all_20 <- mgcv::bam(
  next_terrestrial ~
    behavior_broad * (
      hfi_point_z +
        days_since_emig_z +
        distance_to_nest_km_z +
        cos_Diel +
        sin_Time
    ) +
    s(`individual.local.identifier`, bs = "re"),
  
  family = binomial(link = "logit"),
  data = transition_model_all_20,
  method = "fREML",
  discrete = TRUE
)

summary(m_hfi_point_all_20)


# 2.5 Model 2: mean HFI within 500 m radius ----

m_hfi_500m_all_20 <- mgcv::bam(
  next_terrestrial ~
    behavior_broad * (
      hfi_mean_500m_z +
        days_since_emig_z +
        distance_to_nest_km_z +
        cos_Diel +
        sin_Time
    ) +
    s(`individual.local.identifier`, bs = "re"),
  
  family = binomial(link = "logit"),
  data = transition_model_all_20,
  method = "fREML",
  discrete = TRUE
)

summary(m_hfi_500m_all_20)


# 2.6 Model 3: mean HFI within 1000 m radius ----

m_hfi_1000m_all_20 <- mgcv::bam(
  next_terrestrial ~
    behavior_broad * (
      hfi_mean_1000m_z +
        days_since_emig_z +
        distance_to_nest_km_z +
        cos_Diel +
        sin_Time
    ) +
    s(`individual.local.identifier`, bs = "re"),
  
  family = binomial(link = "logit"),
  data = transition_model_all_20,
  method = "fREML",
  discrete = TRUE
)

summary(m_hfi_1000m_all_20)


# 2.7 Compare candidate models with AIC ----
# The null model tells us whether HFI improves the transition model.
# The three HFI models tell us which spatial scale is best supported.

model_comparison_all_20 <- AIC(
  m_null_20,
  m_hfi_point_all_20,
  m_hfi_500m_all_20,
  m_hfi_1000m_all_20
) |>
  as.data.frame() |>
  tibble::rownames_to_column("model") |>
  mutate(
    delta_AIC = AIC - min(AIC),
    AIC_weight = exp(-0.5 * delta_AIC) / sum(exp(-0.5 * delta_AIC))
  ) |>
  arrange(AIC)

print(model_comparison_all_20)

#'##############################################################################
#' ### Step 3 : plot the results
#' This plot is inspired by Figure 4 in VonBank et al. (2023).
#'
#' Panels:
#'   - Null model
#'   - HFI point
#'   - HFI 500m
#'   - HFI 1000m
#'
#' Rows:
#'   - To Aerial
#'   - To Terrestrial
#'
#' Columns:
#'   - From Aerial
#'   - From Terrestrial
#'
#' Values:
#'   HFI effect on transition log-odds.
#'
#' In the null model, all HFI effects are zero by definition.
#'##############################################################################


# 3.1 Helper function to identify interaction terms ----

find_interaction_term <- function(model, hfi_term) {
  
  coef_names <- names(coef(model))
  
  possible_terms <- c(
    paste0("behavior_broadterrestrial:", hfi_term),
    paste0(hfi_term, ":behavior_broadterrestrial")
  )
  
  interaction_term <- possible_terms[possible_terms %in% coef_names]
  
  if (length(interaction_term) != 1) {
    stop(
      "Could not uniquely identify the interaction term for ",
      hfi_term,
      ". Check coefficient names with names(coef(model))."
    )
  }
  
  interaction_term
}


# 3.2 Helper function to extract HFI transition effects ----

extract_transition_effects_combined <- function(model, hfi_term, panel_name) {
  
  b <- coef(model)
  
  interaction_term <- find_interaction_term(
    model = model,
    hfi_term = hfi_term
  )
  
  # Effect of HFI when departure state is flight:
  # effect on flight -> terrestrial
  beta_flight_to_terrestrial <- b[[hfi_term]]
  
  # Effect of HFI when departure state is terrestrial:
  # effect on terrestrial -> terrestrial
  beta_terrestrial_to_terrestrial <- b[[hfi_term]] + b[[interaction_term]]
  
  tibble(
    panel = panel_name,
    from_state = c("Aerial", "Aerial", "Terrestrial", "Terrestrial"),
    to_state = c("Aerial", "Terrestrial", "Aerial", "Terrestrial"),
    
    log_odds_estimate = c(
      -beta_flight_to_terrestrial,
      beta_flight_to_terrestrial,
      -beta_terrestrial_to_terrestrial,
      beta_terrestrial_to_terrestrial
    )
  )
}


# 3.3 Null model transition effects ----
# The null model contains no HFI term.
# Therefore, the HFI effect is zero for all transitions.

transition_plot_null_20 <- tibble(
  panel = "Null model",
  from_state = c("Aerial", "Aerial", "Terrestrial", "Terrestrial"),
  to_state = c("Aerial", "Terrestrial", "Aerial", "Terrestrial"),
  log_odds_estimate = c(0, 0, 0, 0)
)


# 3.4 Extract transition effects from HFI models ----

transition_plot_hfi_point_20 <- extract_transition_effects_combined(
  model = m_hfi_point_all_20,
  hfi_term = "hfi_point_z",
  panel_name = "HFI point"
)

transition_plot_hfi_500m_20 <- extract_transition_effects_combined(
  model = m_hfi_500m_all_20,
  hfi_term = "hfi_mean_500m_z",
  panel_name = "HFI 500m"
)

transition_plot_hfi_1000m_20 <- extract_transition_effects_combined(
  model = m_hfi_1000m_all_20,
  hfi_term = "hfi_mean_1000m_z",
  panel_name = "HFI 1000m"
)


# 3.5 Combine plot data ----

transition_plot_all_20 <- bind_rows(
  transition_plot_null_20,
  transition_plot_hfi_point_20,
  transition_plot_hfi_500m_20,
  transition_plot_hfi_1000m_20
) |>
  mutate(
    label = sprintf("%.2f", log_odds_estimate),
    
    from_state = factor(
      from_state,
      levels = c("Aerial", "Terrestrial")
    ),
    
    to_state = factor(
      to_state,
      levels = c("Aerial", "Terrestrial")
    ),
    
    panel = factor(
      panel,
      levels = c("Null model", "HFI point", "HFI 500m", "HFI 1000m")
    )
  )

print(transition_plot_all_20)


# 3.6 Plot complete transition-effect matrix ----

max_abs_hfi_effect_all_20 <- max(
  abs(transition_plot_all_20$log_odds_estimate),
  na.rm = TRUE
)

p_hfi_transition_all_20 <- ggplot(
  transition_plot_all_20,
  aes(x = from_state, y = to_state, fill = log_odds_estimate)
) +
  geom_tile(color = "grey70", linewidth = 0.6) +
  geom_text(aes(label = label), size = 4) +
  facet_wrap(~ panel, nrow = 1) +
  scale_fill_gradient2(
    low = "#3B5BDB",
    mid = "grey92",
    high = "#FF6B4A",
    midpoint = 0,
    limits = c(-max_abs_hfi_effect_all_20, max_abs_hfi_effect_all_20),
    name = "HFI effect\n(log-odds)"
  ) +
  labs(
    x = "From",
    y = "To"
  ) +
  coord_equal() +
  theme_bw(base_size = 14) +
  theme(
    strip.background = element_rect(fill = "white", color = "white"),
    strip.text = element_text(face = "bold", size = 14),
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 11)
  )

print(p_hfi_transition_all_20)


# 3.7 Save figure ----

ggsave(
  filename = "transition_hfi_effects_combined_20min.png",
  plot = p_hfi_transition_all_20,
  width = 12,
  height = 4.5,
  dpi = 300
)

