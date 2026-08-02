#' ----------------------------------------------------------------------------- 
# Title: Data Processing and Annotation with Static covariates ----
#' Authors : Louise Faure
#' Date : 28.07.26
#' Info : this script follow the gps_data_preparation.R script where data are thinned
#' at a 60 minutes interval. 
#' **Purpose :** 
#' (1) identify turning angle and step lenght of landing steps 
#' (2) generate landing random steps
#' (3) annotates both random and observed data point with environmental covariates
#' -----------------------------------------------------------------------------

# libararies
library(tidyverse)
library(lubridate)
library(move)
library(sf)
library(mapview)
library(parallel)
library(CircStats)
library(circular)
library(fitdistrplus)
library(terra)

# Golden eagle gps data thinned at a 60 minute resolution
ge_thinned <- readRDS("/Users/louisefaure/Desktop/dossier sans titre/donnees filtree/GE_60_min_thinned.rds")


#------------------------------------------------------------------------------- STEP 1: step selection prep- generate alternative steps ----
#' **Steps**: 
#' (i) calculate step lenght and turning angle for both landing 
#' (ii) estimate step length and turning angle distributions
#' (iii) produce alternative steps
#' (iv) annotate with life stages.  



# Parameters ----
# move() expects an sp::CRS object
wgs <- sp::CRS("+proj=longlat +datum=WGS84 +no_defs")

# Used later when alternative endpoints are generated
meters_proj <- sp::CRS("+proj=utm +zone=32 +datum=WGS84 +units=m +no_defs")
n_alt <- 50


# 1.0 Prepare the already thinned 60-min dataset ----
# Remove the sf geometry column if ge_thinned is an sf object.
# Coordinates are already available in location.long and location.lat.
ge_thinned_prepared_60 <- if (inherits(ge_thinned, "sf")) {sf::st_drop_geometry(ge_thinned)} else {as.data.frame(ge_thinned)}

required_columns_60 <- c("location.long","location.lat","timestamp","individual.local.identifier","behavior_binary","behavior_cluster","new_burst",
  "burst_n","burst_id","n_points_burst","row_in_burst", "cos_diel", "sin_time")

# Order observations and independently recalculate burst indices
ge_thinned_prepared_60 <- ge_thinned_prepared_60 %>%
  dplyr::mutate( behavior_binary =stringr::str_to_lower(as.character(behavior_binary))) %>%
  dplyr::filter(
    is.finite(lon),
    is.finite(lat),
    !is.na(timestamp),
    !is.na(individual.local.identifier),
    !is.na(burst_id)) %>%
  dplyr::arrange(individual.local.identifier,timestamp) %>%
  dplyr::group_by(individual.local.identifier,burst_id) %>%
  dplyr::mutate(
    n_points_burst_calculated = dplyr::n(),
    row_in_burst_calculated   = dplyr::row_number(),
    new_burst_calculated      = dplyr::row_number() == 1) %>%
  dplyr::ungroup()

# 1.1 Create one Move object per individual ----
ge_thinned_move_60 <- ge_thinned_prepared_60 %>%
  dplyr::arrange(individual.local.identifier,timestamp) %>%
  dplyr::mutate(individual.local.identifier = as.character(individual.local.identifier)) %>%
  as.data.frame()

rownames(ge_thinned_move_60) <- NULL

# control: summarize sampling effort by individual
control_individual_sampling_60 <- ge_thinned_move_60 %>%
  dplyr::count(
    individual.local.identifier,
    burst_id,
    name = "n_points_in_burst") %>%
  dplyr::group_by(
    individual.local.identifier) %>%
  dplyr::summarise(
    n_points =
      sum(n_points_in_burst),
    n_bursts =
      dplyr::n(),
    n_bursts_with_3_points =
      sum(n_points_in_burst >= 3),
    median_points_per_burst =
      stats::median(n_points_in_burst),
    max_points_per_burst =
      max(n_points_in_burst),
    .groups = "drop") %>%
  dplyr::arrange(n_points,n_bursts)

print(control_individual_sampling_60,n = Inf)

# we exclude the individual Langgries which has only 4 gps points
ge_thinned_move_60 <- ge_thinned_move_60 %>%
  dplyr::filter(individual.local.identifier !="Langgries21 (eobs 7586)")

# Split the filtered data frame before constructing Move objects 
track_data_list_60 <- split(
  ge_thinned_move_60,
  ge_thinned_move_60$individual.local.identifier,
  drop = TRUE)


# Create one Move object for each individual ----
track_list_60 <- lapply(
  track_data_list_60,
  function(track_data) {
    
    individual_id <- unique(
      track_data$individual.local.identifier
    )
    
    if (length(individual_id) != 1) {
      stop(
        "More than one individual found in a track: ",
        paste(individual_id, collapse = ", ")
      )
    }
    
    move::move(
      x = track_data$lon,
      y = track_data$lat,
      time = track_data$timestamp,
      data = track_data,
      proj = wgs,
      animal = individual_id,
      sensor = "GPS"
    )
  }
)


# Assign explicit names to the list ----
names(track_list_60) <- names(
  track_data_list_60
)



# 1.2 Calculate step lengths and turning angles by burst ----
mycl <- parallel::makeCluster(10)
parallel::clusterEvalQ(
  mycl,
  {
    library(tidyverse)
    library(move)
  }
)

calculation_start_60 <- Sys.time()


sp_obj_ls <- parallel::parLapply(
  mycl,
  track_list_60,
  function(track) {
    
    # Split row indices according to the existing burst_id.
    # No thinning and no reconstruction of bursts are performed here.
    burst_indices <- split(
      seq_len(nrow(track@data)),
      as.character(track@data$burst_id)
    )
    
    # A turning angle requires three successive positions.
    burst_indices <- burst_indices[
      lengths(burst_indices) >= 3
    ]
    
    if (length(burst_indices) == 0) {
      return(NULL)
    }
    
    
    burst_list <- lapply(
      burst_indices,
      function(indices) {
        
        burst <- track[indices, ]
        
        # Distance of the step ending at the current row:
        # row i contains the distance from row i - 1 to row i.
        burst$step_length_m <- c(
          NA_real_,
          move::distance(burst)
        )
        
        
        turning_angle_raw <- move::turnAngleGc(
          burst
        )
        
        
        # Angle located at the turning vertex.
        # This reproduces the formal output alignment recommended by move:
        # point i contains the angle between segments i-1 -> i and i -> i+1.
        burst$turning_angle_vertex_deg <- c(
          NA_real_,
          turning_angle_raw,
          NA_real_
        )
        
        
        # Step-aligned turning angle.
        #
        # Row i contains:
        # - the step length from i-1 to i;
        # - the turning angle of that same step relative to the previous step.
        #
        # This alignment is more convenient for subsequent SSF generation,
        # because the observed endpoint carries both movement variables.
        burst$turning_angle_deg <- c(
          NA_real_,
          NA_real_,
          turning_angle_raw
        )
        
        
        # Define the behavior at the beginning and end of each step.
        behavior_end <- as.character(
          burst$behavior_binary
        )
        
        behavior_start <- dplyr::lag(
          behavior_end
        )
        
        burst$behavior_start <- behavior_start
        burst$behavior_end   <- behavior_end
        
        
        # Classify the transition represented by each endpoint row.
        burst$step_type <- dplyr::case_when(
          behavior_start == "aerial" &
            behavior_end == "terrestrial" ~ "landing",
          
          behavior_start == "aerial" &
            behavior_end == "aerial" ~ "aerial_persistence",
          
          behavior_start == "terrestrial" &
            behavior_end == "aerial" ~ "takeoff",
          
          behavior_start == "terrestrial" &
            behavior_end == "terrestrial" ~ "terrestrial_persistence",
          
          TRUE ~ NA_character_
        )
        
        burst
      }
    )
    
    
    # Recombine bursts belonging to the same individual
    do.call(
      rbind,
      burst_list
    )
  }
)


parallel::stopCluster(mycl)

Sys.time() - calculation_start_60


# Remove individuals with no burst containing at least three points
sp_obj_ls <- Filter(Negate(is.null),sp_obj_ls)

# 2.0 Combine all individuals in one data frame ----
stopifnot(!is.null(names(sp_obj_ls)), all(names(sp_obj_ls) != ""))

bursted_df_list_60 <- lapply(
  seq_along(sp_obj_ls),
  function(i) {
    as.data.frame(sp_obj_ls[[i]]) %>%
      dplyr::mutate(
        individual.local.identifier = names(sp_obj_ls)[i]
      )
  }
)

bursted_df <- dplyr::bind_rows(bursted_df_list_60) %>%
  dplyr::select(
    -dplyr::any_of(c("coords.x1", "coords.x2"))
  )

# Retain the movement steps originating from an aerial location.
aerial_origin_steps_60 <- bursted_df %>%
  dplyr::filter(
    step_type %in%
      c("landing","aerial_persistence"))

landing_steps_60 <- aerial_origin_steps_60 %>% dplyr::filter(step_type == "landing")

aerial_persistence_steps_60 <- aerial_origin_steps_60 %>% dplyr::filter(step_type == "aerial_persistence")

# control: number of steps by transition type
control_step_type_counts_60 <- bursted_df %>%
  dplyr::count(
    step_type,
    name = "n_steps") %>% dplyr::arrange(dplyr::desc(n_steps))

print(control_step_type_counts_60,)

# control: number of usable landing and persistence steps 
control_aerial_step_summary_60 <- aerial_origin_steps_60 %>%
  dplyr::group_by(step_type) %>%
  dplyr::summarise(
    n_steps = dplyr::n(),
    n_with_step_length = sum( !is.na(step_length_m) & step_length_m > 0),
    n_with_turning_angle = sum(!is.na(turning_angle_deg)),
    n_complete_movement_steps =
      sum(
        !is.na(step_length_m) &
          step_length_m > 0 &
          !is.na(turning_angle_deg)),
    .groups = "drop")

print(control_aerial_step_summary_60,n = Inf)

# 2.1 Select landing steps used to estimate movement distributions ----
# Estimate the movement kernel only from aerial-to-terrestrial transitions:
# conditional on landing, where does the eagle land?
landing_distribution_steps_60 <- landing_steps_60 %>%
  dplyr::filter(!is.na(step_length_m), step_length_m > 0, !is.na(turning_angle_deg))

# control: number of complete landing steps ----
control_landing_distribution_60 <- landing_distribution_steps_60 %>%
  dplyr::summarise(n_landing_steps = dplyr::n(),
                   min_step_length_km = min(step_length_m) / 1000,
                   median_step_length_km = stats::median(step_length_m) / 1000,
                   mean_step_length_km = mean(step_length_m) / 1000,
                   max_step_length_km = max(step_length_m) / 1000)
print(control_landing_distribution_60)

# 2.2 Prepare landing movement distributions ----
# Strictly positive step lengths
landing_step_length_km_60 <- landing_distribution_steps_60$step_length_m / 1000
landing_step_length_km_60 <- landing_step_length_km_60[ is.finite(landing_step_length_km_60) & landing_step_length_km_60 > 0]

# Convert degrees to radians and wrap angles to [-pi, pi]
landing_turning_angles_rad_60 <- landing_distribution_steps_60$turning_angle_deg * pi / 180
landing_turning_angles_rad_60 <- atan2(sin(landing_turning_angles_rad_60),cos(landing_turning_angles_rad_60))
landing_turning_angles_rad_60 <-landing_turning_angles_rad_60[is.finite(landing_turning_angles_rad_60)]

# 2.3 Compare gamma and exponential step-length distributions ----
landing_gamma_fit_60 <- fitdistrplus::fitdist(landing_step_length_km_60,distr = "gamma",method = "mle")
landing_exponential_fit_60 <-fitdistrplus::fitdist(landing_step_length_km_60,distr = "exp",method = "mle")

landing_length_fits_60 <- list(Gamma = landing_gamma_fit_60,Exponential = landing_exponential_fit_60)

# Model comparison
landing_length_model_comparison_60 <- tibble::tibble(
  distribution = c("Gamma", "Exponential"),
  n_parameters = c(2L, 1L),
  log_likelihood = c(
    as.numeric(stats::logLik(landing_gamma_fit_60)),
    as.numeric(stats::logLik(landing_exponential_fit_60))),
  AIC = c(
    stats::AIC(landing_gamma_fit_60),
    stats::AIC(landing_exponential_fit_60)),
  BIC = c(
    stats::BIC(landing_gamma_fit_60),
    stats::BIC(landing_exponential_fit_60))) %>%
  dplyr::mutate(
    delta_AIC = AIC - min(AIC),
    delta_BIC = BIC - min(BIC)) %>%
  dplyr::arrange(AIC)

# distribution n_parameters log_likelihood    AIC    BIC delta_AIC delta_BIC
# Gamma                   2        -14179. 28361. 28374.        0         0 
# Exponential             1        -14424. 28851. 28858.      490.      483.

# Additional goodness-of-fit statistics
landing_length_gof_60 <-fitdistrplus::gofstat(landing_length_fits_60,fitnames = c("Gamma", "Exponential"))
# Goodness-of-fit statistics
# Gamma  Exponential
# Kolmogorov-Smirnov statistic 0.02906608   0.09160358
# Cramer-von Mises statistic   1.18109643  18.90251197
# Anderson-Darling statistic   6.22622608 121.34001566
# 
# Goodness-of-fit criteria
# Gamma Exponential
# Akaike's Information Criterion 28361.11    28850.88
# Bayesian Information Criterion 28374.43    28857.54

# Gamma has a smaller BIC but also outperform the other distribution, we choose gamma distribution for step lenght. 


# 2.4 Compare von Mises and uniform turning-angle distributions ----
landing_turning_angles_circular_60 <-circular::circular(landing_turning_angles_rad_60,units = "radians")

# Von Mises maximum-likelihood fit
landing_vonmises_fit_60 <-circular::mle.vonmises(landing_turning_angles_circular_60)
landing_mu_60 <-landing_vonmises_fit_60$mu
landing_kappa_60 <-landing_vonmises_fit_60$kappa

# Von Mises log-likelihood
landing_angle_loglik_vonmises_60 <-
  sum(
    log(
      pmax(
        circular::dvonmises(
          landing_turning_angles_rad_60,
          mu = landing_mu_60,
          kappa = landing_kappa_60
        ),
        .Machine$double.xmin)))

# Uniform log-likelihood on [-pi, pi]
landing_angle_loglik_uniform_60 <- length(landing_turning_angles_rad_60) *log(1 / (2 * pi))

n_landing_angles_60 <-length(landing_turning_angles_rad_60)

# Uniform has no estimated parameter when support is fixed to [-pi, pi].
# Von Mises estimates mu and kappa.
landing_angle_model_comparison_60 <- tibble::tibble(
  distribution = c("Von Mises", "Uniform"),
  n_parameters = c(2L, 0L),
  log_likelihood = c(
    landing_angle_loglik_vonmises_60,
    landing_angle_loglik_uniform_60),
  AIC = c(
    -2 * landing_angle_loglik_vonmises_60 + 2 * 2,
    -2 * landing_angle_loglik_uniform_60),
  BIC = c(
    -2 * landing_angle_loglik_vonmises_60 +
      log(n_landing_angles_60) * 2,
    -2 * landing_angle_loglik_uniform_60
  )) %>%
  dplyr::mutate(
    delta_AIC = AIC - min(AIC),
    delta_BIC = BIC - min(BIC)
  ) %>%
  dplyr::arrange(AIC)


# 2.5 Display fitted parameters and comparison tables ----
landing_gamma_fit_60$estimate
landing_exponential_fit_60$estimate

landing_mu_60
landing_kappa_60
print(landing_angle_model_comparison_60)
# Circular Data: 
# Type = angles 
# Units = radians 
# Template = none 
# Modulo = asis 
# Zero = 0 
# Rotation = counter 
# [1] 0.9574975
# landing_kappa_60 : 0.007229684

# distribution n_parameters log_likelihood    AIC    BIC delta_AIC delta_BIC
#   1 Uniform                 0        -10617. 21235. 21235.      0          0  
# 2 Von Mises               2        -10617. 21239. 21252.      3.85      17.2
# we choose a uniform distribution of turning angle. 


# 2.6 Plot the selected landing movement distributions ----
# The 99.5% quantile is used only to improve visualization.
# All observations were retained when fitting the distributions.
landing_step_length_plot_max_60 <- stats::quantile(
  landing_step_length_km_60,
  probs = 0.995,
  na.rm = TRUE
)

graphics::par(
  mfrow = c(1, 2),
  mar = c(4.5, 4.5, 1, 1)
)

# Gamma distribution retained for landing step lengths
graphics::hist(
  landing_step_length_km_60,
  probability = TRUE,
  breaks = 50,
  xlim = c(0, landing_step_length_plot_max_60),
  main = "",
  xlab = "Landing step length (km)",
  ylab = "Density"
)

graphics::curve(
  stats::dgamma(
    x,
    shape = landing_gamma_fit_60$estimate[["shape"]],
    rate = landing_gamma_fit_60$estimate[["rate"]]
  ),
  add = TRUE,
  from = 0,
  to = landing_step_length_plot_max_60,
  col = "blue",
  lwd = 2
)

graphics::legend(
  "topright",
  legend = "Gamma distribution",
  col = "blue",
  lwd = 2,
  bty = "n"
)

# Uniform distribution retained for landing turning angles
graphics::hist(
  landing_turning_angles_rad_60,
  probability = TRUE,
  breaks = seq(-pi, pi, length.out = 37),
  xlim = c(-pi, pi),
  main = "",
  xlab = "Landing turning angle (radians)",
  ylab = "Density"
)

graphics::abline(
  h = 1 / (2 * pi),
  col = "red",
  lwd = 2
)

graphics::legend(
  "topright",
  legend = "Uniform distribution",
  col = "red",
  lwd = 2,
  bty = "n"
)

# Reset graphical layout
graphics::par(mfrow = c(1, 1))

#------------------------------------------------------------------------------ STEP 3: Generate observed and available landing destinations ----

# 3.0 Parameters ----
n_alt_60 <- 50L
landing_gamma_shape_60 <- unname(landing_gamma_fit_60$estimate[["shape"]])
landing_gamma_rate_60 <- unname(landing_gamma_fit_60$estimate[["rate"]])

intermediate_dataset_directory <- paste0(
  "/Users/louisefaure/Library/CloudStorage/OneDrive-Personnel/THESE/",
  "CHAPITRE 2/git/chapter-2/ACC&GPS_HMM/Results/Intermediate_dataset")

dir.create(intermediate_dataset_directory, recursive = TRUE, showWarnings = FALSE)
set.seed(20260728)

# 3.1 Prepare observed landing steps ----
# Each landing row is the terrestrial endpoint.
# The preceding row is the aerial origin and the row before that defines
# the direction of movement preceding the landing step.
# This step calculate the angle in UTM. 
landing_generation_input_60 <- bursted_df %>%
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
    timestamp_prior = dplyr::lag(timestamp,2),
    timestamp_origin = dplyr::lag(timestamp,1),
    timestamp_destination = timestamp,
    
    prior_lon = dplyr::lag(lon,2),
    prior_lat = dplyr::lag(lat,2),
    
    origin_lon = dplyr::lag(lon,1),
    origin_lat = dplyr::lag(lat,1),
    
    observed_destination_lon = lon,
    observed_destination_lat = lat
  ) %>%
  dplyr::ungroup() %>%
  dplyr::filter(
    step_type == "landing",
    !is.na(timestamp_prior),
    !is.na(timestamp_origin),
    dplyr::if_all(
      dplyr::all_of(c(
        "prior_lon","prior_lat",
        "origin_lon","origin_lat",
        "observed_destination_lon",
        "observed_destination_lat"
      )),
      is.finite
    )
  )

# Project coordinates to the metric CRS defined for step generation ----
project_xy_60 <- function(longitude,latitude) {
  points <- sf::st_as_sf(
    data.frame(longitude = longitude,latitude = latitude),
    coords = c("longitude","latitude"),
    crs = 4326
  )
  
  sf::st_coordinates(
    sf::st_transform(points,crs = 3035)
  )
}

prior_xy_60 <- project_xy_60(
  landing_generation_input_60$prior_lon,
  landing_generation_input_60$prior_lat
)

origin_xy_60 <- project_xy_60(
  landing_generation_input_60$origin_lon,
  landing_generation_input_60$origin_lat
)

destination_xy_60 <- project_xy_60(
  landing_generation_input_60$observed_destination_lon,
  landing_generation_input_60$observed_destination_lat
)

# Calculate observed movement geometry in metres and radians ----
landing_generation_input_60 <- landing_generation_input_60 %>%
  dplyr::mutate(
    prior_x_m = prior_xy_60[,1],
    prior_y_m = prior_xy_60[,2],
    
    origin_x_m = origin_xy_60[,1],
    origin_y_m = origin_xy_60[,2],
    
    observed_destination_x_m = destination_xy_60[,1],
    observed_destination_y_m = destination_xy_60[,2],
    
    previous_heading_rad = atan2(
      origin_y_m - prior_y_m,
      origin_x_m - prior_x_m
    ),
    
    observed_heading_rad = atan2(
      observed_destination_y_m - origin_y_m,
      observed_destination_x_m - origin_x_m
    ),
    
    observed_step_length_m = sqrt(
      (observed_destination_x_m - origin_x_m)^2 +
        (observed_destination_y_m - origin_y_m)^2
    ),
    
    observed_turning_angle_rad = atan2(
      sin(observed_heading_rad - previous_heading_rad),
      cos(observed_heading_rad - previous_heading_rad)
    )
  ) %>%
  dplyr::filter(
    is.finite(observed_step_length_m),
    observed_step_length_m > 0,
    is.finite(observed_turning_angle_rad)
  ) %>%
  dplyr::mutate(
    landing_id = dplyr::row_number(),
    stratum = paste0("landing_",landing_id)
  )

# Columns that remain constant within each choice set
landing_metadata_columns_60 <- c(
  "individual.local.identifier", "burst_id", "row_in_burst",
  "row_in_burst_calculated", "landing_step_id", "stratum",
  "timestamp_origin", "timestamp_destination", "weeks_since_emig",
  "behavior_start", "behavior_end", "step_type",
  "origin_lon", "origin_lat", "origin_x_m", "origin_y_m",
  "previous_heading_rad"
)

# 3.3 Format observed landing destinations ----
observed_landing_locations_60 <- landing_generation_input_60 %>%
  dplyr::select(
    dplyr::any_of(landing_metadata_columns_60),
    
    destination_x_m = observed_destination_x_m,
    destination_y_m = observed_destination_y_m,
    
    step_length_m = observed_step_length_m,
    turning_angle_rad = observed_turning_angle_rad,
    heading_rad = observed_heading_rad
  ) %>%
  dplyr::mutate(
    used = 1L,
    alternative_id = 0L,
    location_type = "observed"
  )

# 3.4 Generate available landing destinations ----
n_landing_steps_60 <- nrow(landing_generation_input_60)
n_available_steps_60 <- n_landing_steps_60 * n_alt_60

available_landing_locations_60 <-
  landing_generation_input_60[
    rep(seq_len(n_landing_steps_60), each = n_alt_60),
  ] %>%
  dplyr::select(dplyr::any_of(landing_metadata_columns_60)) %>%
  dplyr::mutate(
    alternative_id = rep(seq_len(n_alt_60), times = n_landing_steps_60),
    step_length_m = stats::rgamma(
      n_available_steps_60,
      shape = landing_gamma_shape_60,
      rate = landing_gamma_rate_60
    ) * 1000,
    turning_angle_rad = stats::runif(
      n_available_steps_60,
      min = -pi,
      max = pi
    ),
    heading_rad = atan2(
      sin(previous_heading_rad + turning_angle_rad),
      cos(previous_heading_rad + turning_angle_rad)
    ),
    destination_x_m = origin_x_m + step_length_m * cos(heading_rad),
    destination_y_m = origin_y_m + step_length_m * sin(heading_rad),
    used = 0L,
    location_type = "available"
  )

# 3.5 Combine observed and available destinations ----
issf_generated_projected_60 <- dplyr::bind_rows(
  observed_landing_locations_60,
  available_landing_locations_60
)

# Convert destination coordinates back to longitude and latitude
issf_generated_sf_60 <- issf_generated_projected_60 %>%
  sf::st_as_sf(
    coords = c("destination_x_m","destination_y_m"),
    crs = 3035,
    remove = FALSE
  ) %>%
  sf::st_transform(4326)

destination_coordinates_60 <- sf::st_coordinates(issf_generated_sf_60)

issf_generated_observed_location <- issf_generated_sf_60 %>%
  sf::st_drop_geometry() %>%
  dplyr::mutate(
    location.long = destination_coordinates_60[, 1],
    location.lat = destination_coordinates_60[, 2],
    step_length_km = step_length_m / 1000,
    turning_angle_deg = turning_angle_rad * 180 / pi,
    heading_deg = (90 - heading_rad * 180 / pi) %% 360
  ) %>%
  dplyr::arrange(
    individual.local.identifier, timestamp_origin,
    dplyr::desc(used), alternative_id
  )

print(
  issf_generated_observed_location %>%
    dplyr::summarise(
      n_landing_steps = dplyr::n_distinct(stratum),
      n_observed = sum(used == 1L),
      n_available = sum(used == 0L),
      n_individuals = dplyr::n_distinct(individual.local.identifier)
    )
)



# 3.7 Save generated iSSF dataset ----
saveRDS(
  issf_generated_observed_location,
  file = file.path(
    intermediate_dataset_directory,
    "issf_generated_observed_location.rds"
  )
)

#------------------------------------------------------------------------------
# STEP 4: Summary statistics ----
# 4.1 Sampling effort per individual ----
issf_generated_observed_location %>%
  dplyr::mutate(
    yr_mn = paste(lubridate::year(timestamp_destination),
                  lubridate::month(timestamp_destination), sep = "_"),
    yr_day = paste(lubridate::year(timestamp_destination),
                   lubridate::yday(timestamp_destination), sep = "_"),
    yr = lubridate::year(timestamp_destination)
  ) %>%
  dplyr::group_by(individual.local.identifier) %>%
  dplyr::summarise(
    n_d = dplyr::n_distinct(yr_day),
    n_m = dplyr::n_distinct(yr_mn),
    n_yr = dplyr::n_distinct(yr),
    min_yr = min(yr),
    .groups = "drop"
  ) %>%
  dplyr::arrange(n_d) %>%
  as.data.frame()

# 4.2 Data quantity per month ----
mn_summary_60 <- issf_generated_observed_location %>%
  dplyr::mutate(
    mn = lubridate::month(timestamp_destination)
  ) %>%
  dplyr::group_by(individual.local.identifier, mn) %>%
  dplyr::summarise(
    data = dplyr::n(),
    .groups = "drop"
  )

graphics::barplot(
  names.arg = mn_summary_60$mn,
  height = mn_summary_60$data,
  col = as.factor(mn_summary_60$individual.local.identifier),
  beside = FALSE
)


# STEP 5: annotation: life stages ---------------------------------------------

emig_dates <- readRDS("/Users/louisefaure/Library/CloudStorage/OneDrive-Personnel/THESE/CHAPITRE 2/git/chapter-2/DONNEES AIGLES/emigration dates/emigration_dates_20250417.rds")

issf_generated_observed_location_with_age <- issf_generated_observed_location %>%
  dplyr::mutate(timestamp = as.POSIXct(timestamp_destination, tz = "UTC")) %>%
  dplyr::left_join(
    emig_dates %>% dplyr::select(individual.local.identifier, dispersal_date),
    by = "individual.local.identifier"
  ) %>%
  dplyr::mutate(
    days_since_emig = ceiling(as.numeric(difftime(timestamp, dispersal_date, units = "days"))),
    weeks_since_emig = ceiling(as.numeric(difftime(timestamp, dispersal_date, units = "weeks")))
  )


# STEP 6: annotation: static --------------------------------------------------

# all raster layers have the same CRS and resolution (100 m)
settlement_density <- terra::rast("/Users/louisefaure/Desktop/dossier sans titre/Rasters/settlement_density_1km2_100m.tif")
population_density <- terra::rast("/Users/louisefaure/Desktop/dossier sans titre/Rasters/population_density_1km2_100m.tif")
ruggedness <- terra::rast("/Users/louisefaure/Desktop/dossier sans titre/Rasters/ruggedness_100m.tif")
dist_ridgeline <- terra::rast("/Users/louisefaure/Desktop/dossier sans titre/Rasters/distance_to_ridgeline_100m.tif")
landcover <- terra::rast("/Users/louisefaure/Desktop/dossier sans titre/Rasters/landcover_100m.tif")
elevation <- terra::rast("/Users/louisefaure/Desktop/dossier sans titre/Rasters/elevation_100m.tif")

terra::terraOptions(threads = 5)

# 6.1 Function: proportions in the central cell and four rook neighbours ----
extract_landcover_5cells <- function(landcover_raster, points_raster) {
  central_cells <- terra::cellFromXY(
    landcover_raster,
    terra::crds(points_raster)
  )
  
  cells_5 <- terra::adjacent(
    landcover_raster,
    cells = central_cells,
    directions = "rook",
    pairs = FALSE,
    include = TRUE
  )
  
  cells_vector <- as.vector(t(cells_5))
  values_vector <- rep(NA_real_, length(cells_vector))
  valid_cells <- is.finite(cells_vector)
  
  values_vector[valid_cells] <- landcover_raster[
    cells_vector[valid_cells]
  ][, 1]
  
  values_5 <- matrix(
    values_vector,
    nrow = nrow(cells_5),
    ncol = ncol(cells_5),
    byrow = TRUE
  )
  
  n_forest <- rowSums(
    values_5 == 2 | values_5 == 3 | values_5 == 4,
    na.rm = TRUE
  )
  
  n_low_vegetation <- rowSums(
    values_5 == 5 | values_5 == 6 | values_5 == 7,
    na.rm = TRUE
  )
  
  n_rocky_terrain <- rowSums(
    values_5 == 8 | values_5 == 9 | values_5 == 11,
    na.rm = TRUE
  )
  
  complete_neighbourhood <- rowSums(!is.na(values_5)) == 5
  n_other <- 5 - n_forest - n_low_vegetation - n_rocky_terrain
  
  data.frame(
    prop_forest_5cells = ifelse(
      complete_neighbourhood, n_forest / 5, NA_real_
    ),
    prop_low_vegetation_5cells = ifelse(
      complete_neighbourhood, n_low_vegetation / 5, NA_real_
    ),
    prop_rocky_terrain_5cells = ifelse(
      complete_neighbourhood, n_rocky_terrain / 5, NA_real_
    ),
    prop_other_5cells = ifelse(
      complete_neighbourhood, n_other / 5, NA_real_
    )
  )
}
# 6.2 Split by individual to limit memory use ----
issf_annotation_list_60 <- split(
  issf_generated_observed_location_with_age,
  as.character(issf_generated_observed_location_with_age$individual.local.identifier)
)

# 6.3 Extract point-level covariates and mean HFI within 1,000 m ----
issf_annotated_list_60 <- lapply(seq_along(issf_annotation_list_60), function(i) {
  individual_name <- names(issf_annotation_list_60)[i]
  cat("Processing individual", i, "of", length(issf_annotation_list_60),
      ":", individual_name, "\n")
  flush.console()
  
  df <- issf_annotation_list_60[[i]]
  if (nrow(df) == 0) return(df)
  
  pts_raster <- terra::vect(df,geom = c("destination_x_m","destination_y_m"),crs = "EPSG:3035")
  
  # Values under each observed or available destination
  df$elevation_100m <- terra::extract(elevation, pts_raster,
                                      method = "simple", ID = FALSE)[, 1]
  df$ruggedness_100m <- terra::extract(ruggedness, pts_raster,
                                       method = "simple", ID = FALSE)[, 1]
  df$distance_to_ridgeline_100m <- terra::extract(
    dist_ridgeline, pts_raster, method = "simple", ID = FALSE
  )[, 1]
  
  df$settlement_density <- terra::extract(
    settlement_density,
    pts_raster,
    method = "simple",
    ID = FALSE
  )[,1]
  
  df$population_density <- terra::extract(
    population_density,
    pts_raster,
    method = "simple",
    ID = FALSE
  )[,1]
  
  # Land-cover proportions over the central cell and four neighbours
  landcover_5cells <- extract_landcover_5cells(landcover, pts_raster)
  df <- dplyr::bind_cols(df, landcover_5cells)
  
  rm(pts, pts_raster, landcover_5cells)
  gc()
  df
})

names(issf_annotated_list_60) <- names(issf_annotation_list_60)

# 6.4 Recombine annotated individuals ----
issf_generated_observed_location_annotated <- dplyr::bind_rows(
  issf_annotated_list_60
) %>%
  dplyr::arrange(
    individual.local.identifier, timestamp_origin,
    dplyr::desc(used), alternative_id
  ) %>%
  as.data.frame()

# 6.5 Controls ----
issf_annotation_control_60 <- issf_generated_observed_location_annotated %>%
  dplyr::summarise(
    n_rows = dplyr::n(),
    n_strata = dplyr::n_distinct(stratum),
    n_individuals = dplyr::n_distinct(individual.local.identifier),
    missing_hfi = sum(is.na(settlement_density)),
    missing_elevation = sum(is.na(elevation_100m)),
    missing_ruggedness = sum(is.na(ruggedness_100m)),
    missing_ridgeline = sum(is.na(distance_to_ridgeline_100m)),
    missing_landcover = sum(is.na(prop_low_vegetation_5cells))
  )

print(issf_annotation_control_60)

# Check that observed and available destinations have environmental variation
issf_within_stratum_control_60 <- issf_generated_observed_location_annotated %>%
  dplyr::group_by(stratum) %>%
  dplyr::summarise(
    hfi_range = diff(range(settlement_density, na.rm = TRUE)),
    elevation_range = diff(range(settlement_density, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  dplyr::summarise(
    n_strata = dplyr::n(),
    n_strata_with_hfi_variation = sum(is.finite(hfi_range) & hfi_range > 0),
    n_strata_with_elevation_variation =
      sum(is.finite(elevation_range) & elevation_range > 0)
  )

print(issf_within_stratum_control_60)

# 6.6 Save annotated iSSF dataset ----
saveRDS(
  issf_generated_observed_location_annotated,
  file = "/Users/louisefaure/Library/CloudStorage/OneDrive-Personnel/THESE/CHAPITRE 2/git/chapter-2/ACC&GPS_HMM/Results/Intermediate_dataset/issf_generated_observed_location_annotated(2).rds",
  compress = "gzip"
)

