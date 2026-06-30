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

# 1.6 Focus on transitions departing from flight ----
# This is the dataset needed for the departure-point model.
#
# Response:
# 0 = flight -> flight
# 1 = flight -> terrestrial

transition_flight_departure_20 <- transitions_20 |>
  filter(behavior_broad == "flight")

# 1.7 Prepare the three departure-point model datasets ----
# Each object corresponds to one HFI spatial scale.
# The rows are the same type of transition:
# flight -> flight
# flight -> terrestrial
#
# They differ only in the HFI covariate that will be used in Step 2.

transition_flight_hfi_point_20 <- transition_flight_departure_20 |>
  drop_na(
    next_terrestrial,
    hfi_point_z,
    days_since_emig_z,
    distance_to_nest_km_z,
    cos_Diel,
    sin_Time,
    `individual.local.identifier`
  )

transition_flight_hfi_500m_20 <- transition_flight_departure_20 |>
  drop_na(
    next_terrestrial,
    hfi_mean_500m_z,
    days_since_emig_z,
    distance_to_nest_km_z,
    cos_Diel,
    sin_Time,
    `individual.local.identifier`
  )

transition_flight_hfi_1000m_20 <- transition_flight_departure_20 |>
  drop_na(
    next_terrestrial,
    hfi_mean_1000m_z,
    days_since_emig_z,
    distance_to_nest_km_z,
    cos_Diel,
    sin_Time,
    `individual.local.identifier`
  )


#'##############################################################################
#' ### Step 2 : model fitting
#' **Philosophy**: we fit one binomial GAMM per HFI spatial scale.
#' The response variable is `next_terrestrial`:
#'   0 = flight -> flight
#'   1 = flight -> terrestrial
#'
#' We include individual identity as a random intercept using:
#'   s(individual.local.identifier, bs = "re")
#'
#' This accounts for repeated observations from the same individual and allows
#' each individual to deviate from the population-level intercept.
#'##############################################################################


# 2.1 Ensure individual ID is treated as a factor ----
# mgcv random effects require a factor variable.

transition_flight_hfi_point_20 <- transition_flight_hfi_point_20 |>
  mutate(
    `individual.local.identifier` = factor(`individual.local.identifier`)
  )

transition_flight_hfi_500m_20 <- transition_flight_hfi_500m_20 |>
  mutate(
    `individual.local.identifier` = factor(`individual.local.identifier`)
  )

transition_flight_hfi_1000m_20 <- transition_flight_hfi_1000m_20 |>
  mutate(
    `individual.local.identifier` = factor(`individual.local.identifier`)
  )

# 2.2 Model 1: HFI under GPS point ----

m_hfi_point_20 <- mgcv::bam(
  next_terrestrial ~
    hfi_point_z +
    days_since_emig_z +
    distance_to_nest_km_z +
    cos_Diel +
    sin_Time +
    s(`individual.local.identifier`, bs = "re"),
  
  family = binomial(link = "logit"),
  data = transition_flight_hfi_point_20,
  method = "fREML",
  discrete = TRUE
)

summary(m_hfi_point_20)

#' Family: binomial 
#' Link function: logit 
#' Formula:
#'  next_terrestrial ~ hfi_point_z + days_since_emig_z + distance_to_nest_km_z + 
#'  cos_Diel + sin_Time + s(individual.local.identifier, bs = "re")
#'Parametric coefficients:
#'  Estimate Std. Error z value Pr(>|z|)    
#'  (Intercept)            0.81515    0.05462  14.925  < 2e-16 ***
#'  hfi_point_z           -0.08486    0.01112  -7.632 2.32e-14 ***
#'  days_since_emig_z      0.10466    0.01384   7.563 3.93e-14 ***
#'  distance_to_nest_km_z -0.08374    0.01680  -4.986 6.17e-07 ***
#'  cos_Diel               0.98806    0.05212  18.957  < 2e-16 ***
#'  sin_Time               0.37318    0.02727  13.685  < 2e-16 ***
#'  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#'  
#'  Approximate significance of smooth terms:
#'   edf Ref.df Chi.sq p-value    
#'   s(individual.local.identifier) 46.52     61  266.5  <2e-16 ***
#'   
#'   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#'   
#'   R-sq.(adj) =  0.0364   Deviance explained = 2.83%
#'   fREML =  40883  Scale est. = 1         n = 25615


# 2.3 Model 2: mean HFI within 500 m radius ----

m_hfi_500m_20 <- mgcv::bam(
  next_terrestrial ~
    hfi_mean_500m_z +
    days_since_emig_z +
    distance_to_nest_km_z +
    cos_Diel +
    sin_Time +
    s(`individual.local.identifier`, bs = "re"),
  
  family = binomial(link = "logit"),
  data = transition_flight_hfi_500m_20,
  method = "fREML",
  discrete = TRUE
)

summary(m_hfi_500m_20)

# Family: binomial 
# Link function: logit 
# 
# Formula:
#   next_terrestrial ~ hfi_mean_500m_z + days_since_emig_z + distance_to_nest_km_z + 
#   cos_Diel + sin_Time + s(individual.local.identifier, bs = "re")
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)            0.81599    0.05462  14.939  < 2e-16 ***
#   hfi_mean_500m_z       -0.09064    0.01115  -8.126 4.45e-16 ***
#   days_since_emig_z      0.10128    0.01389   7.293 3.04e-13 ***
#   distance_to_nest_km_z -0.08270    0.01680  -4.922 8.55e-07 ***
#   cos_Diel               0.98702    0.05213  18.935  < 2e-16 ***
#   sin_Time               0.37426    0.02727  13.722  < 2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value    
# s(individual.local.identifier) 46.51     61  267.1  <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.0367   Deviance explained = 2.85%
# fREML =  40879  Scale est. = 1         n = 25615

# 2.4 Model 3: mean HFI within 1000 m radius ----

m_hfi_1000m_20 <- mgcv::bam(
  next_terrestrial ~
    hfi_mean_1000m_z +
    days_since_emig_z +
    distance_to_nest_km_z +
    cos_Diel +
    sin_Time +
    s(`individual.local.identifier`, bs = "re"),
  
  family = binomial(link = "logit"),
  data = transition_flight_hfi_1000m_20,
  method = "fREML",
  discrete = TRUE
)

summary(m_hfi_1000m_20)

# Family: binomial 
# Link function: logit 
# 
# Formula:
#   next_terrestrial ~ hfi_mean_1000m_z + days_since_emig_z + distance_to_nest_km_z + 
#   cos_Diel + sin_Time + s(individual.local.identifier, bs = "re")
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)            0.81564    0.05458  14.943  < 2e-16 ***
#   hfi_mean_1000m_z      -0.10205    0.01166  -8.751  < 2e-16 ***
#   days_since_emig_z      0.09624    0.01397   6.891 5.53e-12 ***
#   distance_to_nest_km_z -0.08092    0.01681  -4.815 1.47e-06 ***
#   cos_Diel               0.98473    0.05214  18.888  < 2e-16 ***
#   sin_Time               0.37564    0.02728  13.769  < 2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value    
# s(individual.local.identifier) 46.47     61  266.9  <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.0371   Deviance explained = 2.88%
# fREML =  40874  Scale est. = 1         n = 25615

# 2.5 Compare models ----

model_comparison_20 <- AIC(
  m_hfi_point_20,
  m_hfi_500m_20,
  m_hfi_1000m_20
)

print(model_comparison_20)

# df      AIC
# m_hfi_point_20 53.00632 34608.48
# m_hfi_500m_20  52.99884 34600.41
# m_hfi_1000m_20 52.95309 34589.84

# 2.6 Extract HFI effects ----

extract_hfi_effect <- function(model, hfi_term, model_name) {
  
  param_table <- summary(model)$p.table
  
  out <- as.data.frame(param_table) |>
    tibble::rownames_to_column("term") |>
    filter(term == hfi_term) |>
    mutate(
      model = model_name,
      odds_ratio = exp(Estimate)
    ) |>
    dplyr::select(
      model,
      term,
      Estimate,
      `Std. Error`,
      `z value`,
      `Pr(>|z|)`,
      odds_ratio
    )
  
  return(out)
}

hfi_effects_20 <- bind_rows(
  extract_hfi_effect(
    model = m_hfi_point_20,
    hfi_term = "hfi_point_z",
    model_name = "hfi_point"
  ),
  extract_hfi_effect(
    model = m_hfi_500m_20,
    hfi_term = "hfi_mean_500m_z",
    model_name = "hfi_500m"
  ),
  extract_hfi_effect(
    model = m_hfi_1000m_20,
    hfi_term = "hfi_mean_1000m_z",
    model_name = "hfi_1000m"
  )
)

print(hfi_effects_20)
# model             term    Estimate Std. Error   z value     Pr(>|z|) odds_ratio
# 1 hfi_point      hfi_point_z -0.08486482 0.01111994 -7.631773 2.315476e-14  0.9186365
# 2  hfi_500m  hfi_mean_500m_z -0.09063786 0.01115443 -8.125726 4.446937e-16  0.9133484
# 3 hfi_1000m hfi_mean_1000m_z -0.10205227 0.01166152 -8.751197 2.111017e-18  0.9029823


#'##############################################################################
#' ### Step 2bis : model fitting for transitions departing from terrestrial
#' **Philosophy**: we now fit the same type of binomial GAMM, but for transitions
#' departing from terrestrial behavior.
#'
#' Response variable:
#'   next_terrestrial = 1 means terrestrial -> terrestrial
#'   next_terrestrial = 0 means terrestrial -> flight
#'
#' Therefore, the HFI coefficient estimates the effect of HFI on the probability
#' of remaining terrestrial, conditional on being terrestrial at the departure point.
#'##############################################################################


# 2bis.1 Focus on transitions departing from terrestrial ----

transition_terrestrial_departure_20 <- transitions_20 |>
  filter(behavior_broad == "terrestrial")


# 2bis.2 Prepare the three terrestrial-departure model datasets ----

transition_terrestrial_hfi_point_20 <- transition_terrestrial_departure_20 |>
  drop_na(
    next_terrestrial,
    hfi_point_z,
    days_since_emig_z,
    distance_to_nest_km_z,
    cos_Diel,
    sin_Time,
    `individual.local.identifier`
  ) |>
  mutate(
    `individual.local.identifier` = factor(`individual.local.identifier`)
  )

transition_terrestrial_hfi_500m_20 <- transition_terrestrial_departure_20 |>
  drop_na(
    next_terrestrial,
    hfi_mean_500m_z,
    days_since_emig_z,
    distance_to_nest_km_z,
    cos_Diel,
    sin_Time,
    `individual.local.identifier`
  ) |>
  mutate(
    `individual.local.identifier` = factor(`individual.local.identifier`)
  )

transition_terrestrial_hfi_1000m_20 <- transition_terrestrial_departure_20 |>
  drop_na(
    next_terrestrial,
    hfi_mean_1000m_z,
    days_since_emig_z,
    distance_to_nest_km_z,
    cos_Diel,
    sin_Time,
    `individual.local.identifier`
  ) |>
  mutate(
    `individual.local.identifier` = factor(`individual.local.identifier`)
  )


# 2bis.3 Check terrestrial-departure model datasets ----

model_data_summary_terrestrial_20 <- bind_rows(
  hfi_point = transition_terrestrial_hfi_point_20,
  hfi_500m = transition_terrestrial_hfi_500m_20,
  hfi_1000m = transition_terrestrial_hfi_1000m_20,
  .id = "hfi_model"
) |>
  group_by(hfi_model) |>
  summarise(
    n_transitions = n(),
    n_individuals = n_distinct(`individual.local.identifier`),
    n_bursts = n_distinct(burst_uid),
    n_terrestrial_to_flight = sum(next_terrestrial == 0),
    n_terrestrial_to_terrestrial = sum(next_terrestrial == 1),
    p_terrestrial_to_terrestrial_raw = mean(next_terrestrial),
    .groups = "drop"
  )

print(model_data_summary_terrestrial_20)


# 2bis.4 Model 1: terrestrial departure, HFI under GPS point ----

m_hfi_point_20_terrestrial <- mgcv::bam(
  next_terrestrial ~
    hfi_point_z +
    days_since_emig_z +
    distance_to_nest_km_z +
    cos_Diel +
    sin_Time +
    s(`individual.local.identifier`, bs = "re"),
  
  family = binomial(link = "logit"),
  data = transition_terrestrial_hfi_point_20,
  method = "fREML",
  discrete = TRUE
)

summary(m_hfi_point_20_terrestrial)


# 2bis.5 Model 2: terrestrial departure, mean HFI within 500 m radius ----

m_hfi_500m_20_terrestrial <- mgcv::bam(
  next_terrestrial ~
    hfi_mean_500m_z +
    days_since_emig_z +
    distance_to_nest_km_z +
    cos_Diel +
    sin_Time +
    s(`individual.local.identifier`, bs = "re"),
  
  family = binomial(link = "logit"),
  data = transition_terrestrial_hfi_500m_20,
  method = "fREML",
  discrete = TRUE
)

summary(m_hfi_500m_20_terrestrial)


# 2bis.6 Model 3: terrestrial departure, mean HFI within 1000 m radius ----

m_hfi_1000m_20_terrestrial <- mgcv::bam(
  next_terrestrial ~
    hfi_mean_1000m_z +
    days_since_emig_z +
    distance_to_nest_km_z +
    cos_Diel +
    sin_Time +
    s(`individual.local.identifier`, bs = "re"),
  
  family = binomial(link = "logit"),
  data = transition_terrestrial_hfi_1000m_20,
  method = "fREML",
  discrete = TRUE
)

summary(m_hfi_1000m_20_terrestrial)


# 2bis.7 Compare terrestrial-departure models ----

model_comparison_terrestrial_20 <- AIC(
  m_hfi_point_20_terrestrial,
  m_hfi_500m_20_terrestrial,
  m_hfi_1000m_20_terrestrial
)

print(model_comparison_terrestrial_20)


# 2bis.8 Extract HFI effects for terrestrial-departure models ----

hfi_effects_terrestrial_20 <- bind_rows(
  extract_hfi_effect(
    model = m_hfi_point_20_terrestrial,
    hfi_term = "hfi_point_z",
    model_name = "hfi_point_terrestrial_departure"
  ),
  extract_hfi_effect(
    model = m_hfi_500m_20_terrestrial,
    hfi_term = "hfi_mean_500m_z",
    model_name = "hfi_500m_terrestrial_departure"
  ),
  extract_hfi_effect(
    model = m_hfi_1000m_20_terrestrial,
    hfi_term = "hfi_mean_1000m_z",
    model_name = "hfi_1000m_terrestrial_departure"
  )
)

print(hfi_effects_terrestrial_20)

#'##############################################################################
#' ### Step 3 : plot the results
#' This plot is inspired by Figure 4 in VonBank et al. (2023).
#'
#' Panels:
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
#'##############################################################################


# 3.1 Extract HFI coefficients from flight-departure models ----
# For flight-departure models:
#   estimate = effect of HFI on flight -> terrestrial
#   -estimate = effect of HFI on flight -> flight

extract_hfi_coef <- function(model, term, panel_name, from_state_name) {
  
  coef_table <- summary(model)$p.table
  
  beta <- coef_table[term, "Estimate"]
  se   <- coef_table[term, "Std. Error"]
  pval <- coef_table[term, "Pr(>|z|)"]
  
  tibble(
    panel = panel_name,
    from_state = from_state_name,
    term = term,
    estimate_to_terrestrial = beta,
    std_error = se,
    p_value = pval
  )
}


hfi_coef_from_aerial_20 <- bind_rows(
  extract_hfi_coef(
    model = m_hfi_point_20,
    term = "hfi_point_z",
    panel_name = "HFI point",
    from_state_name = "Aerial"
  ),
  extract_hfi_coef(
    model = m_hfi_500m_20,
    term = "hfi_mean_500m_z",
    panel_name = "HFI 500m",
    from_state_name = "Aerial"
  ),
  extract_hfi_coef(
    model = m_hfi_1000m_20,
    term = "hfi_mean_1000m_z",
    panel_name = "HFI 1000m",
    from_state_name = "Aerial"
  )
)


# 3.2 Extract HFI coefficients from terrestrial-departure models ----
# For terrestrial-departure models:
#   estimate = effect of HFI on terrestrial -> terrestrial
#   -estimate = effect of HFI on terrestrial -> aerial

hfi_coef_from_terrestrial_20 <- bind_rows(
  extract_hfi_coef(
    model = m_hfi_point_20_terrestrial,
    term = "hfi_point_z",
    panel_name = "HFI point",
    from_state_name = "Terrestrial"
  ),
  extract_hfi_coef(
    model = m_hfi_500m_20_terrestrial,
    term = "hfi_mean_500m_z",
    panel_name = "HFI 500m",
    from_state_name = "Terrestrial"
  ),
  extract_hfi_coef(
    model = m_hfi_1000m_20_terrestrial,
    term = "hfi_mean_1000m_z",
    panel_name = "HFI 1000m",
    from_state_name = "Terrestrial"
  )
)


# 3.3 Build the complete 2 x 2 transition grid ----

hfi_coef_all_20 <- bind_rows(
  hfi_coef_from_aerial_20,
  hfi_coef_from_terrestrial_20
)

transition_plot_20 <- hfi_coef_all_20 |>
  rowwise() |>
  do({
    beta_val <- .$estimate_to_terrestrial
    
    tibble(
      panel = .$panel,
      from_state = .$from_state,
      to_state = c("Aerial", "Terrestrial"),
      log_odds_estimate = c(-beta_val, beta_val),
      std_error = .$std_error,
      p_value = .$p_value
    )
  }) |>
  ungroup() |>
  mutate(
    label = sprintf("%.2f", log_odds_estimate),
    
    from_state = factor(
      from_state,
      levels = c("Aerial", "Terrestrial")
    ),
    
    # This order puts Terrestrial above Aerial in the plot.
    to_state = factor(
      to_state,
      levels = c("Aerial", "Terrestrial")
    ),
    
    panel = factor(
      panel,
      levels = c("HFI point", "HFI 500m", "HFI 1000m")
    )
  )


print(transition_plot_20)


# 3.4 Plot complete transition-effect matrix ----

max_abs_hfi_effect_20 <- max(
  abs(transition_plot_20$log_odds_estimate),
  na.rm = TRUE
)

p_hfi_transition_20 <- ggplot(
  transition_plot_20,
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
    limits = c(-max_abs_hfi_effect_20, max_abs_hfi_effect_20),
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

print(p_hfi_transition_20)
