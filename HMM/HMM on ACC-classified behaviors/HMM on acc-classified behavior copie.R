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




# libraries
library(dplyr)
library(tidyr)
library(lubridate)
library(stringr)
library(tibble)
library(mgcv)
library(ggplot2)
library(terra)


# golden eagle data
ge_20min <- readRDS("/Users/louisefaure/Library/CloudStorage/OneDrive-Personnel/THESE/CHAPITRE 2/git/chapter-2/HMM/HMM on ACC-classified behaviors/donnees intermediaire (2)/GE_20_min_thinned_behavior_assigned_hfi.rds")
emig_dates <- readRDS("/Users/louisefaure/Library/CloudStorage/OneDrive-Personnel/THESE/CHAPITRE 2/git/chapter-2/DONNEES AIGLES/emigration dates/emigration_dates_20250417.rds")
nest_site <- readRDS("/Users/louisefaure/Library/CloudStorage/OneDrive-Personnel/THESE/CHAPITRE 2/git/chapter-2/DONNEES AIGLES/nest site location/nest_site_location/nest_site_location.rds")


# covariate data
covar_1_distance_ridgeline <- terra::rast("/Users/louisefaure/Library/CloudStorage/OneDrive-Personnel/MEMOIRE M2/donnees/raster/topography/distance_to_ridge_line_complete_version.tif")
covar_2_TRI <- terra::rast("/Users/louisefaure/Library/CloudStorage/OneDrive-Personnel/MEMOIRE M2/eagle_projet/data/pretraitements/DEM_25m/TRI/TRI.tif")


# parameter 
percentile <- 0.5




#'------------------------------------------------------------------------------
##### Step 1 : Identification of the probabilities of transition ----
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

# 1.2 Compute distance to nest ----
# Remove old nest-distance columns if this block has already been run
ge_20min <- ge_20min |>
  dplyr::select(
    -dplyr::any_of(
      c("x_3035", "y_3035", "nest_x_3035", "nest_y_3035", "nest_x_3035.x", "nest_y_3035.x",
        "nest_x_3035.y", "nest_y_3035.y", "distance_to_nest_km", "distance_to_nest_km.x",
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
    .keep_all = TRUE)

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
    ) / 1000)

# 1.3 Extract distance to ridgeline and ruggedness for each GPS point ----
# Convert GPS points to a terra vector
pts_terra <- terra::vect(pts_3035)

# Project GPS points to the CRS of each raster
pts_ridgeline <- terra::project(
  pts_terra,
  terra::crs(covar_1_distance_ridgeline)
)

pts_TRI <- terra::project(
  pts_terra,
  terra::crs(covar_2_TRI)
)

# Extract raster values at GPS locations
ge_20min <- ge_20min |>
  dplyr::mutate(
    covar_1_distance_ridgeline = terra::extract(
      covar_1_distance_ridgeline,
      pts_ridgeline,
      ID = FALSE
    )[[1]],
    
    covar_2_TRI = terra::extract(
      covar_2_TRI,
      pts_TRI,
      ID = FALSE
    )[[1]]
  )

# 1.4 Standardize continuous covariates ----
# We standardize continuous covariates to make model coefficients comparable.
# Temporal covariates are not standardized because they are already cyclic
# variables bounded between -1 and 1.
ge_20min_stand <- ge_20min |>
  dplyr::mutate(
    age_days_z = as.numeric(
      scale(age_days)
    ),
    
    distance_to_nest_km_z = as.numeric(
      scale(distance_to_nest_km)
    ),
    
    distance_ridgeline_z = as.numeric(
      scale(covar_1_distance_ridgeline)
    ),
    
    ruggedness_z = as.numeric(
      scale(covar_2_TRI)
    ),
    
    hfi_point_z = as.numeric(
      scale(hfi_point)
    ),
    
    hfi_mean_500m_z = as.numeric(
      scale(hfi_mean_500m)
    ),
    
    hfi_mean_1000m_z = as.numeric(
      scale(hfi_mean_1000m)
    ),
    
    hfi_q75_500m_z = as.numeric(
      scale(hfi_q75_500m)
    ),
    
    hfi_max_500m_z = as.numeric(
      scale(hfi_max_500m)
    ),
    
    hfi_q75_1000m_z = as.numeric(
      scale(hfi_q75_1000m)
    ),
    
    hfi_max_1000m_z = as.numeric(
      scale(hfi_max_1000m)
    ))

# 1.5 Population-level behavior counts ----
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


# 1.6 Prepare transition table ----
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


# 1.7 Empirical transition counts ----
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

# 1.8 Empirical transition probability matrix ----
row_totals_20 <- rowSums(transition_count_matrix_20)

transition_matrix_20 <- sweep(
  transition_count_matrix_20,
  MARGIN = 1,
  STATS = row_totals_20,
  FUN = "/")

transition_matrix_20[row_totals_20 == 0, ] <- NA

print(round(transition_matrix_20, 3))
#'            to
#'  from      flight resting feeding
#'  flight   0.486   0.476   0.038
#'  resting  0.134   0.828   0.039
#'  feeding  0.099   0.413   0.488


# 1.9 Inspect individuals ----
# Long-format transition table to inspect rare transitions and remove individuals 
# with little transitions
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

# Nombre total de transitions observées par individu
individual_transition_counts_20 <- transitions_20 |>
  dplyr::count(
    individual.local.identifier,
    name = "n_transitions"
  ) |>
  dplyr::arrange(
    n_transitions,
    individual.local.identifier
  )

# Idividuals with the minimul of transitions
individual_transition_counts_20 |>
  dplyr::slice_head(n = 10) |>
  print(n = Inf)
# Untersberg21 (eobs 7501), Krn20 (eobs 7549) and Lassingbach24 (eobs 11915) have
# less than 80 transition for the first fifteen weeks of the dispersal period. 
# We remove then from the study. 
transitions_20 <- transitions_20 |>
  dplyr::filter(
    !individual.local.identifier %in% c(
      "Untersberg21 (eobs 7501)",
      "Krn20 (eobs 7549)",
      "Lassingbach24 (eobs 11915)"
    )
  ) |>
  droplevels()


# 1.10 Population distribution of HFI values ----
# Convert all HFI variables to long format
hfi_violin_data_20 <- transitions_20 |>
  dplyr::select(
    hfi_point,
    hfi_mean_500m,
    hfi_q75_500m,
    hfi_max_500m,
    hfi_mean_1000m,
    hfi_q75_1000m,
    hfi_max_1000m
  ) |>
  tidyr::pivot_longer(
    cols = dplyr::everything(),
    names_to = "hfi_measure",
    values_to = "hfi"
  ) |>
  dplyr::filter(
    is.finite(hfi)
  ) |>
  dplyr::mutate(
    hfi_measure = factor(
      hfi_measure,
      levels = c("hfi_point","hfi_mean_500m","hfi_q75_500m","hfi_max_500m",
                 "hfi_mean_1000m","hfi_q75_1000m","hfi_max_1000m"),
      labels = c("GPS point","Mean\n500 m","Q75\n500 m","Maximum\n500 m", 
                 "Mean\n1000 m","Q75\n1000 m","Maximum\n1000 m")))


# Summary statistics for each HFI measure
hfi_summary_20 <- hfi_violin_data_20 |>
  dplyr::group_by(hfi_measure) |>
  dplyr::summarise(
    n = dplyr::n(),
    median_hfi = stats::median(
      hfi,
      na.rm = TRUE),
    q25 = stats::quantile(
      hfi,
      0.25,
      na.rm = TRUE),
    q75 = stats::quantile(
      hfi,
      0.75,
      na.rm = TRUE),
    .groups = "drop")

print(hfi_summary_20)


# Plot one violin per HFI measure
p_hfi_violin_20 <- ggplot2::ggplot(
  hfi_violin_data_20,
  ggplot2::aes(
    x = hfi_measure,
    y = hfi
  )
) +
  ggplot2::geom_violin(
    scale = "width",
    trim = TRUE,
    linewidth = 0.4
  ) +
  ggplot2::geom_point(
    data = hfi_summary_20,
    ggplot2::aes(
      x = hfi_measure,
      y = median_hfi
    ),
    inherit.aes = FALSE,
    size = 2.5
  ) +
  ggplot2::labs(
    x = "HFI measure and spatial scale",
    y = "Human Footprint Index",
    title = "Distribution of HFI values by measure and spatial scale",
    subtitle = paste(
      "Distributions combine observations from all included individuals;",
      "points indicate medians"
    )
  ) +
  ggplot2::theme_bw(base_size = 12) +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    axis.title = ggplot2::element_text(
      face = "bold"
    ),
    axis.text.x = ggplot2::element_text(
      angle = 45,
      hjust = 1
    )
  )

print(p_hfi_violin_20)




#'------------------------------------------------------------------------------
##### Step 2: model fitting ----
#'
#' We fit:
#' - one null model without HFI;
#' - linear and nonlinear models for point HFI;
#' - linear and nonlinear models for mean, Q75 and maximum HFI at 500 m;
#' - linear and nonlinear models for mean, Q75 and maximum HFI at 1000 m.
#'
#' All models include age, distance to nest, distance to ridgeline,
#' terrain ruggedness, time of day, departure state interactions,
#' and an individual random intercept.
#'
#' The observed covariate values in each row are used during fitting.


state_levels_3 <- c("flight","resting","feeding")
spline_k <- 6L


# 2.1 Prepare one common dataset for all candidate models ----
transition_model_all_20 <- transitions_20 |>
  tidyr::drop_na(
    behavior_grouped,
    behavior_grouped_next,
    individual.local.identifier,
    
    hfi_point_z,
    
    hfi_mean_500m_z,
    hfi_q75_500m_z,
    hfi_max_500m_z,
    
    hfi_mean_1000m_z,
    hfi_q75_1000m_z,
    hfi_max_1000m_z,
    
    age_days_z,
    distance_to_nest_km_z,
    distance_ridgeline_z,
    ruggedness_z,
    cos_Diel,
    sin_Time
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
    
    individual.local.identifier = factor(
      individual.local.identifier
    ),
    
    arrival_state_id =
      as.integer(behavior_grouped_next) - 1L
  ) |>
  droplevels()

# Keep only columns required for model fitting to reduce memory use
model_columns_20 <- c(
  "arrival_state_id",
  "behavior_grouped",
  "individual.local.identifier",
  
  "hfi_point_z",
  
  "hfi_mean_500m_z",
  "hfi_q75_500m_z",
  "hfi_max_500m_z",
  
  "hfi_mean_1000m_z",
  "hfi_q75_1000m_z",
  "hfi_max_1000m_z",
  
  "age_days_z",
  "distance_to_nest_km_z",
  "distance_ridgeline_z",
  "ruggedness_z",
  "cos_Diel",
  "sin_Time"
)

transition_model_fit_20 <- transition_model_all_20 |>
  dplyr::select(
    dplyr::all_of(model_columns_20)
  )


# 2.2 Build linear and nonlinear multinomial formulas ----
base_covariates_20 <- c("age_days_z", "distance_to_nest_km_z", "distance_ridgeline_z",
                        "ruggedness_z", "cos_Diel", "sin_Time")

make_multinom_formula <- function(
    hfi_var = NULL,
    nonlinear = FALSE,
    k = 6L
) {
  base_terms <- paste(
    base_covariates_20,
    collapse = " + "
  )
  
  if (is.null(hfi_var)) {
    
    rhs <- paste0(
      "behavior_grouped * (",
      base_terms,
      ") + ",
      "s(individual.local.identifier, bs = 're')"
    )
    
  } else if (!nonlinear) {
    
    rhs <- paste0(
      "behavior_grouped * (",
      paste(
        c(
          hfi_var,
          base_covariates_20
        ),
        collapse = " + "
      ),
      ") + ",
      "s(individual.local.identifier, bs = 're')"
    )
    
  } else {
    
    rhs <- paste0(
      "behavior_grouped * (",
      base_terms,
      ") + ",
      
      "s(",
      hfi_var,
      ", by = behavior_grouped, ",
      "bs = 'cr', k = ",
      k,
      ") + ",
      
      "s(individual.local.identifier, bs = 're')"
    )
  }
  
  list(
    stats::as.formula(
      paste(
        "arrival_state_id ~",
        rhs
      )
    ),
    
    stats::as.formula(
      paste(
        "~",
        rhs
      )
    )
  )
}


# 2.3 Define all candidate models ----
model_specs_20 <- tibble::tribble(
  ~model,                              ~hfi_variable,         ~hfi_statistic, ~hfi_scale,  ~functional_form,
  
  "m_null_20",                         NA_character_,          "none",         "none",      "none",
  
  "m_hfi_point_linear_20",             "hfi_point_z",          "point",        "point",     "linear",
  "m_hfi_point_nonlinear_20",          "hfi_point_z",          "point",        "point",     "nonlinear",
  
  "m_hfi_mean_500m_linear_20",         "hfi_mean_500m_z",      "mean",         "500 m",     "linear",
  "m_hfi_mean_500m_nonlinear_20",      "hfi_mean_500m_z",      "mean",         "500 m",     "nonlinear",
  
  "m_hfi_q75_500m_linear_20",          "hfi_q75_500m_z",       "Q75",          "500 m",     "linear",
  "m_hfi_q75_500m_nonlinear_20",       "hfi_q75_500m_z",       "Q75",          "500 m",     "nonlinear",
  
  "m_hfi_max_500m_linear_20",          "hfi_max_500m_z",       "maximum",      "500 m",     "linear",
  "m_hfi_max_500m_nonlinear_20",       "hfi_max_500m_z",       "maximum",      "500 m",     "nonlinear",
  
  "m_hfi_mean_1000m_linear_20",        "hfi_mean_1000m_z",     "mean",         "1000 m",    "linear",
  "m_hfi_mean_1000m_nonlinear_20",     "hfi_mean_1000m_z",     "mean",         "1000 m",    "nonlinear",
  
  "m_hfi_q75_1000m_linear_20",         "hfi_q75_1000m_z",      "Q75",          "1000 m",    "linear",
  "m_hfi_q75_1000m_nonlinear_20",      "hfi_q75_1000m_z",      "Q75",          "1000 m",    "nonlinear",
  
  "m_hfi_max_1000m_linear_20",         "hfi_max_1000m_z",      "maximum",      "1000 m",    "linear",
  "m_hfi_max_1000m_nonlinear_20",      "hfi_max_1000m_z",      "maximum",      "1000 m",    "nonlinear"
)


# 2.4 Fit all candidate models in parallel ----
physical_cores <- parallel::detectCores(
  logical = FALSE
)

if (is.na(physical_cores)) {
  physical_cores <- parallel::detectCores(
    logical = TRUE
  )
}

n_model_workers <- min(
  4L,
  max(
    1L,
    physical_cores - 1L
  )
)


fit_candidate_model_20 <- function(i) {
  
  specification <- model_specs_20[i, ]
  
  hfi_variable <- if (
    is.na(specification$hfi_variable)
  ) {
    NULL
  } else {
    specification$hfi_variable
  }
  
  nonlinear <- identical(
    specification$functional_form,
    "nonlinear"
  )
  
  model_formula <- make_multinom_formula(
    hfi_var = hfi_variable,
    nonlinear = nonlinear,
    k = spline_k
  )
  
  mgcv::gam(
    formula = model_formula,
    family = mgcv::multinom(K = 2),
    data = transition_model_fit_20,
    method = "ML",
    control = mgcv::gam.control(
      nthreads = 1L,
      trace = FALSE))}


model_objects_20 <- parallel::mclapply(
  X = seq_len(
    nrow(model_specs_20)),
  FUN = fit_candidate_model_20,
  mc.cores = n_model_workers,
  mc.preschedule = FALSE)

names(model_objects_20) <- model_specs_20$model


# 2.5 Compare all candidate models using AIC ----
model_comparison_all_20 <- model_specs_20 |>
  dplyr::mutate(
    df = vapply(
      model_objects_20,
      function(model) {
        attr(
          stats::logLik(model),
          "df"
        )
      },
      numeric(1)
    ),
    
    AIC = vapply(
      model_objects_20,
      stats::AIC,
      numeric(1)
    ),
    
    delta_AIC =
      AIC - min(AIC),
    
    AIC_weight =
      exp(-0.5 * delta_AIC) /
      sum(
        exp(-0.5 * delta_AIC)
      )
  ) |>
  dplyr::arrange(AIC) |>
  dplyr::select(
    model,
    hfi_statistic,
    hfi_scale,
    functional_form,
    df,
    AIC,
    delta_AIC,
    AIC_weight
  )

print(
  model_comparison_all_20,
  n = Inf
)

#' model                         hfi_statistic hfi_scale functional_form    df    AIC delta_AIC AIC_weight
#' m_hfi_mean_1000m_nonlinear_20 mean          1000 m    nonlinear        152. 77067.       0    1.000e+ 0
#' m_hfi_q75_1000m_nonlinear_20  Q75           1000 m    nonlinear        157. 77085.      17.6  1.50 e- 4
#' m_hfi_mean_1000m_linear_20    mean          1000 m    linear           140. 77098.      30.6  2.25 e- 7
#' m_hfi_max_500m_nonlinear_20   maximum       500 m     nonlinear        161. 77102.      34.5  3.17 e- 8
#' m_hfi_q75_1000m_linear_20     Q75           1000 m    linear           140. 77103.      35.7  1.78 e- 8
#' m_hfi_max_1000m_linear_20     maximum       1000 m    linear           140. 77103.      35.8  1.68 e- 8
#' m_hfi_max_1000m_nonlinear_20  maximum       1000 m    nonlinear        151. 77108.      40.8  1.39 e- 9
#' m_hfi_max_500m_linear_20      maximum       500 m     linear           140. 77113.      45.8  1.13 e-10
#' m_hfi_mean_500m_nonlinear_20  mean          500 m     nonlinear        151. 77115.      48.4  3.14 e-11
#' m_hfi_q75_500m_nonlinear_20   Q75           500 m     nonlinear        154. 77118.      50.8  9.38 e-12
#' m_hfi_q75_500m_linear_20      Q75           500 m     linear           140. 77129.      61.9  3.63 e-14
#' m_hfi_mean_500m_linear_20     mean          500 m     linear           140. 77141.      73.4  1.15 e-16
#' m_hfi_point_nonlinear_20      point         point     nonlinear        151. 77162.      94.5  2.95 e-21
#' m_hfi_point_linear_20         point         point     linear           140. 77171.     104.   2.75 e-23
#' m_null_20                     none          none      none             134. 77305.     238.   1.82 e-52

# 2.6 Save models 
# model_bundle_all_20 <- list(
#   models = model_objects_20,
#   model_specifications = model_specs_20,
#   model_comparison = model_comparison_all_20,
#   model_data = transition_model_all_20,
#   state_levels = state_levels_3,
#   spline_k = spline_k,
#   saved_at = Sys.time(),
#   session_info = utils::sessionInfo()
# )
# 
# step3_output_directory <- paste0("/Users/louisefaure/Library/CloudStorage/OneDrive-Personnel/","THESE/CHAPITRE 2/git/chapter-2/HMM/","HMM on ACC-classified behaviors/Models-fitted")
# 
# dir.create(step3_output_directory,recursive = TRUE,showWarnings = FALSE)

# saveRDS(model_bundle_all_20,file = file.path(step3_output_directory,"candidate_multinomial_models_20min.rds"),ompress = FALSE)

# models <- readRDS("/Users/louisefaure/Library/CloudStorage/OneDrive-Personnel/THESE/CHAPITRE 2/git/chapter-2/HMM/HMM on ACC-classified behaviors/Models-fitted/candidate_multinomial_models_20min.rds")

# Inspect best models details 
candidate_multinomial_models_20min <- readRDS("/Users/louisefaure/Library/CloudStorage/OneDrive-Personnel/THESE/CHAPITRE 2/git/chapter-2/HMM/HMM on ACC-classified behaviors/Models-fitted/candidate_multinomial_models_20min.rds")

best_model_name_20 <- model_comparison_all_20$model[[1]]
best_model_20 <- model_objects_20[[best_model_name_20]]

cat(
  "Best model:",
  best_model_name_20,
  "\n\nModel specification:\n")

print(
  model_specs_20 |>
    dplyr::filter(
      model == best_model_name_20))

cat("\nModel summary:\n")
print(summary(best_model_20))
# Family: multinom 
# Link function: 
#   
#   Formula:
#   arrival_state_id ~ behavior_grouped * (age_days_z + distance_to_nest_km_z + 
#                                            distance_ridgeline_z + ruggedness_z + cos_Diel + sin_Time) + 
#   s(hfi_mean_1000m_z, by = behavior_grouped, bs = "cr", k = 6) + 
#   s(individual.local.identifier, bs = "re")
# <environment: 0x75c27ac10>
#   ~behavior_grouped * (age_days_z + distance_to_nest_km_z + distance_ridgeline_z + 
#                          ruggedness_z + cos_Diel + sin_Time) + s(hfi_mean_1000m_z, 
#                                                                  by = behavior_grouped, bs = "cr", k = 6) + s(individual.local.identifier, 
#                                                                                                               bs = "re")
# <environment: 0x75c27ac10>
#   
#   Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)                                      0.87971    0.06553  13.425  < 2e-16 ***
#   behavior_groupedresting                          2.18237    0.06621  32.960  < 2e-16 ***
#   behavior_groupedfeeding                          1.21898    0.15434   7.898 2.83e-15 ***
#   age_days_z                                       0.13882    0.01573   8.827  < 2e-16 ***
#   distance_to_nest_km_z                           -0.10099    0.01687  -5.987 2.14e-09 ***
#   distance_ridgeline_z                             0.01987    0.01484   1.339 0.180623    
# ruggedness_z                                     0.06743    0.01550   4.351 1.36e-05 ***
#   cos_Diel                                         0.92036    0.06485  14.192  < 2e-16 ***
#   sin_Time                                         0.38794    0.03424  11.329  < 2e-16 ***
#   behavior_groupedresting:age_days_z              -0.09961    0.01950  -5.107 3.27e-07 ***
#   behavior_groupedfeeding:age_days_z              -0.14396    0.04867  -2.958 0.003099 ** 
#   behavior_groupedresting:distance_to_nest_km_z    0.12455    0.01878   6.634 3.28e-11 ***
#   behavior_groupedfeeding:distance_to_nest_km_z    0.20103    0.04895   4.107 4.01e-05 ***
#   behavior_groupedresting:distance_ridgeline_z     0.13750    0.02017   6.815 9.40e-12 ***
#   behavior_groupedfeeding:distance_ridgeline_z     0.07517    0.05154   1.459 0.144684    
# behavior_groupedresting:ruggedness_z            -0.03705    0.01988  -1.864 0.062291 .  
# behavior_groupedfeeding:ruggedness_z             0.02501    0.04947   0.506 0.613056    
# behavior_groupedresting:cos_Diel                 0.74675    0.07792   9.584  < 2e-16 ***
#   behavior_groupedfeeding:cos_Diel                -0.10326    0.18188  -0.568 0.570200    
# behavior_groupedresting:sin_Time                 0.03729    0.04022   0.927 0.353871    
# behavior_groupedfeeding:sin_Time                -0.10616    0.09853  -1.077 0.281309    
# (Intercept).1                                   -2.01288    0.14526 -13.857  < 2e-16 ***
#   behavior_groupedresting.1                        1.40980    0.15200   9.275  < 2e-16 ***
#   behavior_groupedfeeding.1                        4.18941    0.19915  21.037  < 2e-16 ***
#   age_days_z.1                                     0.09906    0.03946   2.510 0.012062 *  
#   distance_to_nest_km_z.1                         -0.14798    0.04201  -3.523 0.000427 ***
#   distance_ridgeline_z.1                           0.08107    0.03623   2.238 0.025227 *  
#   ruggedness_z.1                                   0.04608    0.03984   1.157 0.247412    
# cos_Diel.1                                       0.59156    0.16067   3.682 0.000231 ***
#   sin_Time.1                                       0.17930    0.08651   2.073 0.038216 *  
#   behavior_groupedresting:age_days_z.1            -0.14721    0.04546  -3.238 0.001202 ** 
#   behavior_groupedfeeding:age_days_z.1            -0.21625    0.06008  -3.599 0.000319 ***
#   behavior_groupedresting:distance_to_nest_km_z.1  0.16310    0.04505   3.621 0.000294 ***
#   behavior_groupedfeeding:distance_to_nest_km_z.1  0.21666    0.06058   3.576 0.000349 ***
#   behavior_groupedresting:distance_ridgeline_z.1   0.27997    0.04280   6.541 6.12e-11 ***
#   behavior_groupedfeeding:distance_ridgeline_z.1   0.10580    0.06070   1.743 0.081348 .  
# behavior_groupedresting:ruggedness_z.1          -0.01307    0.04673  -0.280 0.779725    
# behavior_groupedfeeding:ruggedness_z.1           0.09345    0.06124   1.526 0.127039    
# behavior_groupedresting:cos_Diel.1               0.34455    0.17851   1.930 0.053592 .  
# behavior_groupedfeeding:cos_Diel.1               0.35268    0.23266   1.516 0.129554    
# behavior_groupedresting:sin_Time.1              -0.13706    0.09510  -1.441 0.149510    
# behavior_groupedfeeding:sin_Time.1               0.00710    0.12614   0.056 0.955114    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df   Chi.sq  p-value    
# s(hfi_mean_1000m_z):behavior_groupedflight     1.487  1.794   78.999  < 2e-16 ***
#   s(hfi_mean_1000m_z):behavior_groupedresting    3.516  3.949  124.175  < 2e-16 ***
#   s(hfi_mean_1000m_z):behavior_groupedfeeding    1.038  1.074    4.129   0.0443 *  
#   s(individual.local.identifier)                49.074 55.000 1532.007  < 2e-16 ***
#   s.1(hfi_mean_1000m_z):behavior_groupedflight   1.014  1.029   26.983 1.17e-06 ***
#   s.1(hfi_mean_1000m_z):behavior_groupedresting  2.439  2.902    7.965   0.0317 *  
#   s.1(hfi_mean_1000m_z):behavior_groupedfeeding  1.018  1.035    0.002   0.9952    
# s.1(individual.local.identifier)              43.081 55.000  854.663  < 2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Deviance explained = 17.1%
# -REML =  38668  Scale est. = 1         n = 93699





#'------------------------------------------------------------------------------
##### Step 3: predictions and plots for the best AIC model ----
#'
#' The model with the lowest AIC is selected. Predictions are marginalized over 
#' the observed values of age, distance to nest, distance to ridgeline, ruggedness, 
#' and time of day.
#'
#' Predictions are averaged first within individuals and then equally across
#' individuals.



state_levels_3 <- c(
  "flight",
  "resting",
  "feeding"
)

percentile <- 0.05


# 3.1 Select the three models with the lowest AIC ----

top_model_info_20 <- model_comparison_all_20 |>
  dplyr::slice_head(n = 3) |>
  dplyr::mutate(
    AIC_rank = dplyr::row_number()
  ) |>
  dplyr::left_join(
    model_specs_20 |>
      dplyr::select(
        model,
        hfi_variable
      ),
    by = "model"
  ) |>
  dplyr::mutate(
    model_label = paste0(
      "Rank ",
      AIC_rank,
      ": ",
      hfi_statistic,
      " | ",
      hfi_scale,
      " | ",
      functional_form,
      "\nAIC = ",
      sprintf("%.1f", AIC),
      "; ΔAIC = ",
      sprintf("%.1f", delta_AIC)
    )
  )

print(top_model_info_20)


# 3.2 Identify individual random-effect terms ----

get_individual_random_terms_20 <- function(model) {
  
  random_terms <- vapply(
    model$smooth,
    function(smooth_object) {
      smooth_object$label
    },
    character(1)
  )
  
  random_terms <- unique(
    random_terms[
      grepl(
        "individual.local.identifier",
        random_terms,
        fixed = TRUE
      )
    ]
  )
  
  if (length(random_terms) == 0L) {
    random_terms <- NULL
  }
  
  random_terms
}


# 3.3 Marginalized predictions at one HFI value ----

predict_marginalized_hfi_20 <- function(
    model,
    model_data,
    hfi_variable,
    hfi_value,
    hfi_level,
    states,
    excluded_terms = NULL
) {
  
  departure_results <- lapply(
    states,
    function(from_state_i) {
      
      observed_rows <- model_data |>
        dplyr::filter(
          behavior_grouped == from_state_i
        )
      
      prediction_data <- observed_rows
      
      # Change only the HFI variable used by the selected model
      prediction_data[[hfi_variable]] <- hfi_value
      
      predicted_probabilities <- stats::predict(
        model,
        newdata = prediction_data,
        type = "response",
        exclude = excluded_terms
      ) |>
        as.data.frame()
      
      names(predicted_probabilities) <- states
      
      dplyr::bind_cols(
        prediction_data |>
          dplyr::select(
            individual.local.identifier
          ),
        predicted_probabilities
      ) |>
        tidyr::pivot_longer(
          cols = dplyr::all_of(states),
          names_to = "to_state",
          values_to = "transition_probability"
        ) |>
        
        # Average over observed covariates within each individual
        dplyr::group_by(
          individual.local.identifier,
          to_state
        ) |>
        dplyr::summarise(
          transition_probability = mean(
            transition_probability,
            na.rm = TRUE
          ),
          .groups = "drop"
        ) |>
        
        # Give each individual equal weight
        dplyr::group_by(to_state) |>
        dplyr::summarise(
          transition_probability = mean(
            transition_probability,
            na.rm = TRUE
          ),
          n_individuals = dplyr::n(),
          .groups = "drop"
        ) |>
        dplyr::mutate(
          from_state = from_state_i,
          hfi_level = hfi_level,
          hfi_value = hfi_value
        )
    }
  )
  
  dplyr::bind_rows(departure_results)
}


# 3.4 Calculate low-versus-high HFI effects for one model ----

calculate_model_hfi_effect_20 <- function(
    model_name,
    hfi_variable,
    model_label,
    AIC_rank,
    model_data = transition_model_all_20
) {
  
  model <- model_objects_20[[model_name]]
  
  hfi_values <- stats::quantile(
    model_data[[hfi_variable]],
    probs = c(
      percentile,
      1 - percentile
    ),
    na.rm = TRUE,
    names = FALSE
  )
  
  names(hfi_values) <- c(
    "low",
    "high"
  )
  
  individual_random_terms <- get_individual_random_terms_20(
    model
  )
  
  marginal_probabilities <- dplyr::bind_rows(
    predict_marginalized_hfi_20(
      model = model,
      model_data = model_data,
      hfi_variable = hfi_variable,
      hfi_value = hfi_values[["low"]],
      hfi_level = "low",
      states = state_levels_3,
      excluded_terms = individual_random_terms
    ),
    
    predict_marginalized_hfi_20(
      model = model,
      model_data = model_data,
      hfi_variable = hfi_variable,
      hfi_value = hfi_values[["high"]],
      hfi_level = "high",
      states = state_levels_3,
      excluded_terms = individual_random_terms
    )
  )
  
  marginal_probabilities |>
    dplyr::select(
      from_state,
      to_state,
      hfi_level,
      hfi_value,
      transition_probability,
      n_individuals
    ) |>
    tidyr::pivot_wider(
      names_from = hfi_level,
      values_from = c(
        hfi_value,
        transition_probability,
        n_individuals
      ),
      names_glue = "{hfi_level}_{.value}"
    ) |>
    dplyr::rename(
      low_HFI = low_hfi_value,
      high_HFI = high_hfi_value,
      low_probability = low_transition_probability,
      high_probability = high_transition_probability
    ) |>
    dplyr::mutate(
      model_name = model_name,
      model_label = model_label,
      AIC_rank = AIC_rank,
      
      absolute_change =
        high_probability - low_probability,
      
      relative_change =
        dplyr::if_else(
          low_probability > 1e-8,
          absolute_change / low_probability,
          NA_real_
        ),
      
      relative_change_percent =
        100 * relative_change,
      
      from_state = factor(
        from_state,
        levels = state_levels_3
      ),
      
      to_state = factor(
        to_state,
        levels = state_levels_3
      ),
      
      absolute_label = sprintf(
        "%+.3f",
        absolute_change
      ),
      
      relative_label = dplyr::if_else(
        is.finite(relative_change_percent),
        sprintf(
          "%+.1f%%",
          relative_change_percent
        ),
        "NA"
      )
    )
}


# 3.5 Calculate effects for the three best models ----
transition_effects_top3_20 <- lapply(
  seq_len(nrow(top_model_info_20)),
  function(i) {
    
    calculate_model_hfi_effect_20(
      model_name =
        top_model_info_20$model[[i]],
      
      hfi_variable =
        top_model_info_20$hfi_variable[[i]],
      
      model_label =
        top_model_info_20$model_label[[i]],
      
      AIC_rank =
        top_model_info_20$AIC_rank[[i]],
      
      model_data =
        transition_model_all_20
    )
  }
) |>
  dplyr::bind_rows() |>
  dplyr::mutate(
    model_label = factor(
      model_label,
      levels = top_model_info_20$model_label
    )
  ) |>
  dplyr::arrange(
    AIC_rank,
    from_state,
    to_state
  )

print(
  transition_effects_top3_20,
  n = Inf
)


# 3.6 Plot absolute changes in transition probabilities ----
absolute_limit_top3_20 <- max(
  abs(
    transition_effects_top3_20$absolute_change
  ),
  na.rm = TRUE
)

if (
  !is.finite(absolute_limit_top3_20) ||
  absolute_limit_top3_20 == 0
) {
  absolute_limit_top3_20 <- 1e-6
}

p_absolute_top3_20 <- ggplot2::ggplot(
  transition_effects_top3_20,
  ggplot2::aes(
    x = from_state,
    y = to_state,
    fill = absolute_change
  )
) +
  ggplot2::geom_tile(
    color = "grey70",
    linewidth = 0.5
  ) +
  ggplot2::geom_text(
    ggplot2::aes(
      label = absolute_label
    ),
    size = 3.5
  ) +
  ggplot2::facet_wrap(
    ~ model_label,
    nrow = 1
  ) +
  ggplot2::scale_fill_gradient2(
    low = "#3B5BDB",
    mid = "grey95",
    high = "#FF6B4A",
    midpoint = 0,
    limits = c(
      -absolute_limit_top3_20,
      absolute_limit_top3_20
    ),
    name = expression(
      Delta * " probability"
    )
  ) +
  ggplot2::coord_equal() +
  ggplot2::labs(
    x = "Departure state",
    y = "Arrival state",
    
    title =
      "Absolute change in predicted transition probabilities",
    
    subtitle = paste0(
      "Three models with the lowest AIC; P",
      percentile * 100,
      " to P",
      (1 - percentile) * 100,
      "; P(high HFI) − P(low HFI)"
    )
  ) +
  ggplot2::theme_bw(
    base_size = 12
  ) +
  ggplot2::theme(
    panel.grid = ggplot2::element_blank(),
    
    strip.background = ggplot2::element_rect(
      fill = "white",
      color = "grey70"
    ),
    
    strip.text = ggplot2::element_text(
      face = "bold",
      size = 10
    ),
    
    axis.text.x = ggplot2::element_text(
      angle = 45,
      hjust = 1
    )
  )


# 3.7 Plot relative changes in transition probabilities ----
relative_limit_top3_20 <- max(
  abs(
    transition_effects_top3_20$
      relative_change_percent
  ),
  na.rm = TRUE
)

if (
  !is.finite(relative_limit_top3_20) ||
  relative_limit_top3_20 == 0
) {
  relative_limit_top3_20 <- 1e-6
}

p_relative_top3_20 <- ggplot2::ggplot(
  transition_effects_top3_20,
  ggplot2::aes(
    x = from_state,
    y = to_state,
    fill = relative_change_percent
  )
) +
  ggplot2::geom_tile(
    color = "grey70",
    linewidth = 0.5
  ) +
  ggplot2::geom_text(
    ggplot2::aes(
      label = relative_label
    ),
    size = 3.5
  ) +
  ggplot2::facet_wrap(
    ~ model_label,
    nrow = 1
  ) +
  ggplot2::scale_fill_gradient2(
    low = "#3B5BDB",
    mid = "grey95",
    high = "#FF6B4A",
    midpoint = 0,
    limits = c(
      -relative_limit_top3_20,
      relative_limit_top3_20
    ),
    name = "Relative change (%)"
  ) +
  ggplot2::coord_equal() +
  ggplot2::labs(
    x = "Departure state",
    y = "Arrival state",
    
    title =
      "Relative change in predicted transition probabilities",
    
    subtitle = paste0(
      "Three models with the lowest AIC; P",
      percentile * 100,
      " to P",
      (1 - percentile) * 100,
      "; 100 × [P(high HFI) − P(low HFI)] / P(low HFI)"
    )
  ) +
  ggplot2::theme_bw(
    base_size = 12
  ) +
  ggplot2::theme(
    panel.grid = ggplot2::element_blank(),
    
    strip.background = ggplot2::element_rect(
      fill = "white",
      color = "grey70"
    ),
    
    strip.text = ggplot2::element_text(
      face = "bold",
      size = 10
    ),
    
    axis.text.x = ggplot2::element_text(
      angle = 45,
      hjust = 1
    )
  )


# 3.8 Display plots ----
print(p_absolute_top3_20)
print(p_relative_top3_20)

ggplot2::ggsave(
  filename = "/Users/louisefaure/Library/CloudStorage/OneDrive-Personnel/THESE/CHAPITRE 2/git/chapter-2/HMM/HMM on ACC-classified behaviors/Models-fitted/absolute_transition_changes_top3_models_20min.png",
  plot = p_absolute_top3_20,
  width = 12,
  height = 5,
  units = "in",
  dpi = 300
)

ggplot2::ggsave(
  filename = "/Users/louisefaure/Library/CloudStorage/OneDrive-Personnel/THESE/CHAPITRE 2/git/chapter-2/HMM/HMM on ACC-classified behaviors/Models-fitted/relative_transition_changes_top3_models_20min.png",
  plot = p_relative_top3_20,
  width = 12,
  height = 5,
  units = "in",
  dpi = 300
)


# Identify individuals with no observed transitions at or above global Q95 HFI
percentile <- 0.05
best_model_name_20 <- model_comparison_all_20$model[[1]]


best_hfi_variable_20 <- model_specs_20 |>
  dplyr::filter(
    model == best_model_name_20
  ) |>
  dplyr::pull(
    hfi_variable
  ) |>
  .[[1]]

high_hfi_threshold_20 <- stats::quantile(
  transition_model_all_20[[best_hfi_variable_20]],
  probs = 1 - percentile,
  na.rm = TRUE,
  names = FALSE
)

individuals_without_high_hfi_20 <- transition_model_all_20 |>
  dplyr::group_by(
    individual.local.identifier
  ) |>
  dplyr::summarise(
    n_high_HFI = sum(
      .data[[best_hfi_variable_20]] >= high_hfi_threshold_20,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) |>
  dplyr::filter(
    n_high_HFI == 0
  ) |>
  dplyr::pull(
    individual.local.identifier
  )

cat(
  "Global Q95 HFI threshold:",
  high_hfi_threshold_20,
  "\n\nIndividuals with no observations at or above Q95:\n")

print(individuals_without_high_hfi_20)






