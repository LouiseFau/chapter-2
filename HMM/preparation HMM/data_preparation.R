#' ####################################################################### 
#' Title: Preparation of GPS data for HHM
#' Authors : Louise Faure
#' Date : 25.06.26
#' Purpose : (1) create several dataset of various
#' temporal scale, (2) classify behavior within these dataset, (3) extract 
#' covariate for each temporal resolution at different scales
#' #######################################################################



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
################################################################################

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
  
  df
})

dispersal_data <- bind_rows(dispersal_data)
rownames(dispersal_data) <- NULL



################################################################################
#' Step 2. Data cleaning
#' 
#' Goal:
#   (i) inspect raw GPS sampling intervals;
#   (ii) create one high-resolution dataset at 20 min;
#   (iii) create one intermediate-resolution dataset at 60 min;
#   (iv) split tracks into bursts
################################################################################


# 2.1 Inspect raw GPS data ----
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


# 2.2 Select observed gps points at 20-min and 60-min interavals using move2
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

# 2.3 Split retained tracks into bursts ----
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

# 2.4 create 20-min and 60-min datasets ----

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

# 2.5 Control the dataset filtered ----
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

# 2.6 Export datasets ----
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


#'##############################################################################
#'Step 3 : behavioral classification 
#' 
#' ** Strategy**: (1) classify behaviors using embc function (2) sub classify 
#' terrestrial behavior to distinguish between overnight encroachment and short 
#' term resting 
################################################################################

# 3.1: EMbC classification ----
input_vars   <- c("ground.speed", "height_above_ground")
complete_idx <- which(complete.cases(dispersal_data[, input_vars]))
bc_matrix    <- data.matrix(dispersal_data[complete_idx, input_vars])

embc_fit      <- EMbC::embc(bc_matrix)
embc_smoothed <- EMbC::smth(embc_fit, dlta = 0.7)

# complete embc cluster with turning angles
cluster_summary_extended <- dispersal_data[complete_idx, ] %>%
  dplyr::mutate(behavior_cluster = embc_smoothed@A) %>%
  dplyr::group_by(behavior_cluster) %>%
  dplyr::summarise(
    speed_med      = median(ground.speed,        na.rm = TRUE),
    height_med     = median(height_above_ground, na.rm = TRUE),
    turn_angle_med = median(abs(turn.angle),     na.rm = TRUE),
    vert_speed_med = median(abs(vert.speed),     na.rm = TRUE),
    n = n(), .groups = "drop")
print(cluster_summary_extended)

dispersal_data$behavior_cluster <- NA_integer_
dispersal_data$behavior_cluster[complete_idx] <- embc_smoothed@A

# adjust cluster ids below if cluster_summary shows a different ordering
behavior_labels <- c("1" = "terrestrial",
                     "2" = "low soaring",
                     "3" = "fast flight at low elevation",
                     "4" = "fast commuting flight at high elevation")

dispersal_data <- dispersal_data %>%
  dplyr::mutate(behavior = dplyr::recode(as.character(behavior_cluster),
                                         !!!behavior_labels))

# 3.2 : sub-classification of terrestrial position ----
Before setting the time threshold at which we can defined a terrestrial position as
an overnight roosting, explore the time lags between terrestrial position and justify your choice 
timethreshold. 





dispersal_data <- dispersal_data %>%
  dplyr::arrange(individual.local.identifier, timestamp) %>%
  dplyr::group_by(individual.local.identifier) %>%
  dplyr::mutate(
    is_terrestrial = !is.na(behavior) & behavior == "terrestrial",
    bout_change    = is_terrestrial != dplyr::lag(is_terrestrial,
                                                  default = first(is_terrestrial)),
    bout_id        = cumsum(bout_change)) %>%
  dplyr::ungroup()

bout_durations <- dispersal_data %>%
  dplyr::filter(is_terrestrial) %>%
  dplyr::group_by(individual.local.identifier, bout_id) %>%
  dplyr::summarise(
    bout_duration_h = as.numeric(difftime(max(timestamp), min(timestamp),
                                          units = "hours")),
    n_fixes_bout    = n(), .groups = "drop")

dispersal_data <- dispersal_data %>%
  dplyr::left_join(
    bout_durations %>%
      dplyr::select(individual.local.identifier, bout_id, bout_duration_h),
    by = c("individual.local.identifier", "bout_id")) %>%
  dplyr::mutate(
    behavior_refined = dplyr::case_when(
      behavior == "terrestrial" & bout_duration_h >= 10 ~ "overnight roosting",
      behavior == "terrestrial" & bout_duration_h <  10 ~ "short resting",
      TRUE ~ behavior))

dispersal_data %>%
  dplyr::filter(is_terrestrial) %>%
  dplyr::count(behavior_refined) %>%
  print()
