#' ####################################################################### 
#' Title: Preparation of GPS data for HHM
#' Authors : Louise Faure
#' Date : 25.06.26
#' Purpose : 
#' (1) filter location to the first fifteen weeks of dispersal, 
#' (2) classify behaviours,
#' (3) thin the data at two temporal resolution (20min and 60min) and split into
#' burst,
#' (4) extract environmental covariates at the two temporal scale and within 
#' three different buffers
#' #############################################################################



# 0. Setup ----
library(move2)
library(sf)
library(dplyr)
library(tidyverse)
library(lubridate)
library(ggplot2)
library(ggspatial)
library(rnaturalearth)
library(atlastools)
library(data.table)
library(EMbC)
library(suncalc)

# Load GPS data 
output_dir <- "C:/Users/lfaure7/Documents/git/chapter-2/HMM/preparation HMM/donnees intermediaire"

# Géoïde and digital elevation model
geo <- terra::rast("C:/Users/lfaure7/OneDrive/MEMOIRE M2/eagle_projet/data/geoide/us_nga_egm96_15.tif")
dem <- terra::rast("C:/Users/lfaure7/OneDrive/MEMOIRE M2/eagle_projet/data/pretraitements/Region-Alpes-Dem/region-alpes-dem.tif")

# raster layers
human_footprint <- terra::rast("C:/Users/lfaure7/Desktop/COUCHES QGIS/settlements/Overture/human_footprint_index_building_pop_builtprop_100m.tif")

# emigration dates and golden eagle data
emig_dates <- readRDS("C:/Users/lfaure7/Documents/git/chapter-2/preparation donnee aigle ssf/donnees/emigration dates/emigration_dates_20250417.rds")

rds_dir <- "C:/Users/lfaure7/Documents/git/chapter-2/preparation donnee aigle ssf/donnees/no_burst_GE/no_burst_GE"
rds_files <- tools::file_path_sans_ext(list.files(rds_dir, pattern = "\\.rds$", ignore.case = TRUE))
rds_ids   <- as.numeric(gsub("_gpsNoDup_moveObj", "", rds_files))

# Filters individuals to retain those that have an emigration date
emig_dates_filtered <- emig_dates[emig_dates$did_disperse == TRUE &
                                    emig_dates$individual.id %in% rds_ids, ]
rds_files_filtered  <- rds_files[rds_ids %in% emig_dates_filtered$individual.id]


################################################################################
#' ### STEP 1 : spatio temporal filtering 

#' (i) keep only GPS locations within the first 15 weeks after emigration date
#' (ii) compute height above ground (geoid + DEM correction)
#' (iii) compute sunrise and sunset
################################################################################

tz_loc <- "Europe/Zurich"

dispersal_data <- lapply(rds_files_filtered, function(f) {
  
  id  <- as.numeric(gsub("_gpsNoDup_moveObj", "", f))
  obj <- readRDS(file.path(rds_dir, paste0(f, ".rds")))
  df  <- as.data.frame(obj)
  
  # filter column
  df <- df[, c("individual.local.identifier", "timestamp",
               "location.long", "location.lat",
               "height.above.ellipsoid", "ground.speed",
               "eobs.speed.accuracy.estimate", "eobs.horizontal.accuracy.estimate",
               "gps.dop", "vert.speed",
               "turn.angle", "step.length", "gr.speed", "eobs.type.of.fix")]
  
  df$timestamp <- as.POSIXct(df$timestamp, tz = "UTC")
  emig_dt <- emig_dates_filtered$dispersal_date[emig_dates_filtered$individual.id == id]
  
  # retain only the first fifteen weeks of dispersal
  df <- df[df$timestamp >= emig_dt &
             df$timestamp <= emig_dt + 105 * 24 * 3600, ]
  
  df$individual.id   <- id
  df$dispersal_date  <- emig_dt
  df$days_since_emig <- as.numeric(difftime(df$timestamp, emig_dt, units = "days"))
  
  # calculate altitude (height above mean sea level) using the geoide
  xy_ll <- as.matrix(df[, c("location.long", "location.lat")])
  N     <- terra::extract(geo, xy_ll)[, 1]
  
  df$height_msl <- df$height.above.ellipsoid - N
  
  # applied the dem crs
  pts     <- terra::vect(df[, c("location.long", "location.lat")],
                         geom = c("location.long", "location.lat"),
                         crs  = "EPSG:4326")
  pts_dem <- terra::project(pts, crs(dem))
  df$dem_elevation <- terra::extract(dem, pts_dem)[, 2]  
  
  # calculate flight height
  df$height_above_ground <- df$height_msl - df$dem_elevation
  
  # filter the locations to remove all the point below -200m and above 3000m
  df <- df[!is.na(df$height_above_ground) &
             df$height_above_ground > -200 &
             df$height_above_ground < 3000, ]
  # If no data remain after height filtering, return an empty data frame
  if (nrow(df) == 0) {
    return(df)
  }
  
  # Convert timestamp to local time.
  # The date used for sunrise/sunset must be the local date, not the UTC date.
  df$time_local <- lubridate::with_tz(df$timestamp, tz = tz_loc)
  df$date_local <- as.Date(df$time_local)
  
  # Compute sunrise and sunset per row, using local date and GPS location
  sun_times <- suncalc::getSunlightTimes(
    data = data.frame(
      date = df$date_local,
      lat  = df$location.lat,
      lon  = df$location.long
    ),
    keep = c("sunrise", "sunset"),
    tz = tz_loc
  )
  
  # Attach sunrise and sunset to the dataset
  df$sunrise <- sun_times$sunrise
  df$sunset  <- sun_times$sunset
  
  # Define daylight / night for each GPS point
  df$is_daylight <- df$time_local >= df$sunrise & df$time_local <= df$sunset
  df$is_night    <- !df$is_daylight
  
  df
  
})

dispersal_data <- bind_rows(dispersal_data)
rownames(dispersal_data) <- NULL


################################################################################
#' STEP 2 : EMbC behavioural classification on the full filtered dataset
#'
#' **Ref.** Garriga et al., 2016; Nourani et al. workflow
#'
#' (1) Classification through embc function
#' (2) identification of overnight resting versus day resting sites
################################################################################

dispersal_data <- dispersal_data %>%
  arrange(individual.local.identifier, timestamp)

# Bivariate matrix: exactly two variables
behavioural_classification <- data.frame(
  ground.speed = as.numeric(dispersal_data$ground.speed),
  height_above_ground = as.numeric(dispersal_data$height_above_ground)
)

complete_idx <- which(
  complete.cases(behavioural_classification) &
    is.finite(behavioural_classification$ground.speed) &
    is.finite(behavioural_classification$height_above_ground)
)

behavioural_classification <- data.matrix(
  behavioural_classification[complete_idx, c("ground.speed", "height_above_ground")]
)

# Call EMbC
embc <- EMbC::embc(behavioural_classification)

# Smooth EMbC classification, following the published workflow
embc_smoothed <- EMbC::smth(embc, dlta = 0.7)

# Attach raw and smoothed cluster labels to the full filtered dataset
dispersal_data$behavior_cluster_raw <- NA_integer_
dispersal_data$behavior_cluster <- NA_integer_

dispersal_data$behavior_cluster_raw[complete_idx] <- embc@A
dispersal_data$behavior_cluster[complete_idx] <- embc_smoothed@A

# Check the four clusters
cluster_summary_full <- dispersal_data %>%
  filter(!is.na(behavior_cluster)) %>%
  group_by(behavior_cluster) %>%
  summarise(
    speed_med  = median(ground.speed, na.rm = TRUE),
    height_med = median(height_above_ground, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  ) %>%
  arrange(behavior_cluster)

print(cluster_summary_full)

#' We obtain and decide to consider cluster 1 and 2 as terrestrial
#' behavior_cluster speed_med height_med      n
#'  1               0.12       13.5        166905
#'  2               0.81       30.7        20355
#'  3              10.1        34.1        33613
#'  4              13.1       361.0         17154


behavior_cluster_labels <- c(
  "1" = "low_speed_low_height",
  "2" = "low_speed_moderate_height",
  "3" = "high_speed_low_height",
  "4" = "high_speed_high_height"
)

terrestrial_clusters <- c(1, 2)

dispersal_data <- dispersal_data %>%
  dplyr::mutate(
    behavior_cluster_name = dplyr::recode(
      as.character(behavior_cluster),
      !!!behavior_cluster_labels,
      .default = NA_character_
    ),
    
    behavior_broad = dplyr::case_when(
      behavior_cluster %in% terrestrial_clusters ~ "terrestrial",
      behavior_cluster %in% c(3, 4) ~ "flight",
      TRUE ~ NA_character_
    ),
    
    behavior_final = dplyr::case_when(
      behavior_cluster %in% terrestrial_clusters & is_night ~ "overnight_resting",
      behavior_cluster %in% terrestrial_clusters & is_daylight ~ "daily_resting",
      behavior_cluster == 3 ~ "high_speed_low_height",
      behavior_cluster == 4 ~ "high_speed_high_height",
      TRUE ~ NA_character_
    )
  )

# check final classification
behavior_summary_full <- dispersal_data %>%
  dplyr::count(
    behavior_cluster,
    behavior_cluster_name,
    behavior_broad,
    behavior_final
  ) %>%
  dplyr::group_by(behavior_broad) %>%
  dplyr::mutate(prop_within_broad = n / sum(n)) %>%
  dplyr::ungroup()

print(behavior_summary_full)

cluster_summary_named <- dispersal_data %>%
  dplyr::filter(!is.na(behavior_cluster)) %>%
  dplyr::group_by(behavior_cluster, behavior_cluster_name, behavior_broad) %>%
  dplyr::summarise(
    speed_med = median(ground.speed, na.rm = TRUE),
    speed_mean = mean(ground.speed, na.rm = TRUE),
    height_med = median(height_above_ground, na.rm = TRUE),
    height_mean = mean(height_above_ground, na.rm = TRUE),
    n = dplyr::n(),
    .groups = "drop"
  ) %>%
  dplyr::arrange(behavior_cluster)

print(cluster_summary_named)


################################################################################
#' Step 3. Data cleaning
#' 
#' (1) inspect raw GPS sampling intervals;
#' (2) create one high-resolution dataset at 20 min;
#' (3) create one intermediate-resolution dataset at 60 min;
#' (4) split tracks into bursts
################################################################################


# 3.1 Inspect raw GPS data ----
# Create move2 object from dispersal_data
locs <- dispersal_data %>%
  mutate(
    individual.id = as.character(individual.id),
    timestamp = as.POSIXct(timestamp, tz = "UTC"),
    lon = location.long,
    lat = location.lat
  ) %>%
  filter(
    !is.na(individual.id),
    !is.na(timestamp),
    !is.na(lon),
    !is.na(lat)
  ) %>%
  arrange(individual.id, timestamp) %>%
  move2::mt_as_move2(
    coords = c("lon", "lat"),
    time_column = "timestamp",
    track_id_column = "individual.id",
    crs = 4326,
    remove = FALSE
  )

# Remove duplicate records with same individual and timestamp
locs <- locs %>%
  move2::mt_filter_unique(criterion = "first") %>%
  arrange(individual.id, timestamp)

# Project to EPSG:3035 for metric distances
locs_3035 <- locs %>%
  sf::st_transform(3035)

xy_3035 <- sf::st_coordinates(locs_3035)

locs_tbl <- locs_3035 %>%
  mutate(
    x_3035 = xy_3035[, 1],
    y_3035 = xy_3035[, 2]
  ) %>%
  sf::st_drop_geometry() %>%
  arrange(individual.id, timestamp)


# time lags
time_lags <- locs_tbl %>%
  group_by(individual.id) %>%
  arrange(timestamp, .by_group = TRUE) %>%
  mutate(
    dt_min = as.numeric(difftime(timestamp, lag(timestamp), units = "mins"))
  ) %>%
  ungroup() %>%
  filter(!is.na(dt_min), dt_min > 0)

summary(time_lags$dt_min) 

quantile(
  time_lags$dt_min,
  probs = c(0.01, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 0.99),
  na.rm = TRUE
)

#' We found 
#' Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#' 0.10    10.53    15.45    38.62    20.03 24659.60
#' We choose a fine scale resolution of 20 minutes because it is close to the 
#' upper quartile of the main fine-scale sampling distribution and 60 minutes 
#' because the 90% and 95% quantiles are around 60 minutes.


# 3.2 Select observed gps points at 20-min and 60-min interavals using move2
high_res_min <- 20
intermediate_res_min <- 60

high_res_tolerance_min <- 5 # tolerance of 5min for the 20min interval
intermediate_tolerance_min <- 10 # tolerance of 10min for the 20min interval

min_points_per_burst <- 3

format_interval_name <- function(interval_min) {
  gsub("\\.", "p", paste0(interval_min, "min"))
}

# Thin gps data
thin_move2_interval <- function(locs,
                                interval_unit,
                                resolution_min,
                                resolution_label) {
  
  locs_thin <- locs %>%
    arrange(mt_track_id(.), mt_time(.)) %>%
    move2::mt_filter_per_interval(
      unit = interval_unit,
      criterion = "first"
    )
  
  xy <- sf::st_coordinates(locs_thin)
  
  locs_thin_tbl <- locs_thin %>%
    mutate(
      individual.id = as.character(move2::mt_track_id(locs_thin)),
      timestamp = as.POSIXct(move2::mt_time(locs_thin), tz = "UTC"),
      x_3035 = xy[, 1],
      y_3035 = xy[, 2],
      resolution_min = resolution_min,
      resolution_class = resolution_label
    ) %>%
    sf::st_drop_geometry() %>%
    arrange(individual.id, timestamp)
  
  locs_thin_tbl
}

# 3.3 Split retained tracks into bursts ----
split_into_regular_bursts <- function(df,
                                      interval_min,
                                      tolerance_min,
                                      min_points_per_burst = 3) {
  
  df_burst <- df %>%
    arrange(individual.id, timestamp) %>%
    group_by(individual.id) %>%
    mutate(
      
      # Time lag between retained observed GPS fixes
      dt_prev_min = as.numeric(
        difftime(timestamp, lag(timestamp), units = "mins")
      ),
      
      # A new burst starts when: a) its is the first point of the ind, b)
      # the previous retained point is not close enough to the target interval.
      new_burst = is.na(dt_prev_min) |
        abs(dt_prev_min - interval_min) > tolerance_min,
      
      burst_n = cumsum(new_burst),
      
      burst_id = paste(
        individual.id,
        format_interval_name(interval_min),
        sprintf("%04d", burst_n),
        sep = "_"
      )
    ) %>%
    ungroup()
  
  # Remove very short bursts
  df_burst <- df_burst %>%
    group_by(individual.id, burst_id) %>%
    mutate(
      n_points_burst = n()
    ) %>%
    ungroup() %>%
    filter(n_points_burst >= min_points_per_burst)
  
  # Recompute next-step timing after removing short bursts
  df_burst <- df_burst %>%
    arrange(individual.id, burst_id, timestamp) %>%
    group_by(individual.id, burst_id) %>%
    mutate(
      row_in_burst = row_number(),
      
      dt_prev_burst_min = as.numeric(
        difftime(timestamp, lag(timestamp), units = "mins")
      ),
      
      dt_next_burst_min = as.numeric(
        difftime(lead(timestamp), timestamp, units = "mins")
      ),
      
      has_next_regular = !is.na(dt_next_burst_min) &
        abs(dt_next_burst_min - interval_min) <= tolerance_min,
      
      lag_deviation_next_min = dt_next_burst_min - interval_min
    ) %>%
    ungroup()
  
  df_burst
}

# 3.4 create 20-min and 60-min datasets ----

thin_20_tbl <- thin_move2_interval(
  locs = locs_3035,
  interval_unit = "20 minutes",
  resolution_min = high_res_min,
  resolution_label = "20min_high_resolution"
)

thin_60_tbl <- thin_move2_interval(
  locs = locs_3035,
  interval_unit = "1 hour",
  resolution_min = intermediate_res_min,
  resolution_label = "60min_intermediate_resolution"
)

regular_20_tbl <- split_into_regular_bursts(
  df = thin_20_tbl,
  interval_min = high_res_min,
  tolerance_min = high_res_tolerance_min,
  min_points_per_burst = min_points_per_burst
)

regular_60_tbl <- split_into_regular_bursts(
  df = thin_60_tbl,
  interval_min = intermediate_res_min,
  tolerance_min = intermediate_tolerance_min,
  min_points_per_burst = min_points_per_burst
)

# 3.5 Control the dataset filtered ----
check_timing <- function(df) {
  
  df %>%
    filter(has_next_regular == TRUE) %>%
    summarise(
      resolution_min = first(resolution_min),
      n_individuals = n_distinct(individual.local.identifier),
      n_points = n(),
      n_bursts = n_distinct(burst_id),
      n_valid_transitions = sum(has_next_regular, na.rm = TRUE),
      median_dt_next_min = median(dt_next_burst_min, na.rm = TRUE),
      p05_dt_next_min = quantile(dt_next_burst_min, 0.05, na.rm = TRUE),
      p95_dt_next_min = quantile(dt_next_burst_min, 0.95, na.rm = TRUE),
      min_dt_next_min = min(dt_next_burst_min, na.rm = TRUE),
      max_dt_next_min = max(dt_next_burst_min, na.rm = TRUE),
      median_abs_deviation_min = median(abs(lag_deviation_next_min), na.rm = TRUE),
      p95_abs_deviation_min = quantile(abs(lag_deviation_next_min), 0.95, na.rm = TRUE)
    )
}

timing_20 <- check_timing(regular_20_tbl)
timing_60 <- check_timing(regular_60_tbl)

timing_summary <- bind_rows(
  timing_20 %>% mutate(dataset = "20min_high_resolution"),
  timing_60 %>% mutate(dataset = "60min_intermediate_resolution")
)

print(timing_summary)

# inspect individual differences
summarise_by_individual <- function(df) {
  
  df %>%
    group_by(resolution_class, individual.local.identifier) %>%
    summarise(
      n_points_retained = n(),
      n_bursts = n_distinct(burst_id),
      n_valid_transitions = sum(has_next_regular, na.rm = TRUE),
      median_points_per_burst = median(n_points_burst, na.rm = TRUE),
      max_points_per_burst = max(n_points_burst, na.rm = TRUE),
      .groups = "drop"
    )
}

summary_by_id <- bind_rows(
  summarise_by_individual(regular_20_tbl),
  summarise_by_individual(regular_60_tbl)
)

print(summary_by_id, n = 128)

low_data_individuals <- summary_by_id %>%
  filter(n_valid_transitions < 200) %>%
  arrange(resolution_class, n_valid_transitions)

print(low_data_individuals)

# 3.6 Export datasets ----
regular_20_sf <- sf::st_as_sf(
  regular_20_tbl,
  coords = c("x_3035", "y_3035"),
  crs = 3035,
  remove = FALSE
)

regular_60_sf <- sf::st_as_sf(
  regular_60_tbl,
  coords = c("x_3035", "y_3035"),
  crs = 3035,
  remove = FALSE
)

saveRDS(
  regular_20_sf,
  file.path("C:/Users/lfaure7/Documents/git/chapter-2/HMM/preparation HMM/donnees intermediaire/GE_20_min_thinned.rds"))

saveRDS(
  regular_60_sf,
  file.path("C:/Users/lfaure7/Documents/git/chapter-2/HMM/preparation HMM/donnees intermediaire/GE_60_min_thinned.rds"))