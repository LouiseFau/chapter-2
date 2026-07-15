#' ----------------------------------------------------------------------------- 
# Title: Preparation of GPS based behavioral classification data for HHM ----
#' Authors : Louise Faure
#' Date : 25.06.26
#' Purpose : 
#' (1) filter location to the first fifteen weeks of dispersal, 
#' (2) classify behaviours based on gps position,
#' (3) thin the data at two temporal resolution (20min and 60min) and split into
#' burst,
#' (4) extract environmental covariates at the two temporal scale and within 
#' three different buffers (only for HFI)
#' -----------------------------------------------------------------------------




# 0. Setup ----
library(move2)
library(sf)
library(terra)
library(dplyr)
library(tidyverse)
library(lubridate)
library(data.table)
library(EMbC)

# Define out put directory
output_dir <- "/Users/louisefaure/Library/CloudStorage/OneDrive-Personnel/THESE/CHAPITRE 2/git/chapter-2/HMM/preparation HMM/donnees intermediaire"

# Géoïde, digital elevation model and other ground based informations
geo <- terra::rast("/Users/louisefaure/Library/CloudStorage/OneDrive-Personnel/MEMOIRE M2/eagle_projet/data/geoide/us_nga_egm96_15.tif")
dem <- terra::rast("/Users/louisefaure/Library/CloudStorage/OneDrive-Personnel/MEMOIRE M2/eagle_projet/data/pretraitements/Region-Alpes-Dem/region-alpes-dem.tif")
human_footprint <- terra::rast("/Users/louisefaure/Library/CloudStorage/OneDrive-Personnel/THESE/COUCHES QGIS/COUCHES QGIS/settlements/Overture/human_footprint_index_building_pop_builtprop_100m.tif")


ruggedness < terra::rast("/Users/louisefaure/Library/CloudStorage/OneDrive-Personnel/MEMOIRE M2/eagle_projet/data/pretraitements/DEM_25m/TRI/TRI.tif")
slope <- terra::rast("/Users/louisefaure/Library/CloudStorage/OneDrive-Personnel/MEMOIRE M2/eagle_projet/data/pretraitements/DEM_25m/slope/slope_25.tif")
dist_ridgeline <- terra::rast("Users/louisefaure/Library/CloudStorage/OneDrive-Personnel/MEMOIRE M2/donnees/raster/topography/distance_to_ridge_line_complete_version.tif")
landcover <- rast("/Users/louisefaure/Library/CloudStorage/OneDrive-Personnel/MEMOIRE M2/eagle_projet/data_landuse/Data_CLC/CLC_Alps/CLC_longlat/CLC_longlat.tif")

# Emigration dates and golden eagle data
emig_dates <- readRDS("/Users/louisefaure/Library/CloudStorage/OneDrive-Personnel/THESE/CHAPITRE 2/git/chapter-2/DONNEES AIGLES/emigration dates/emigration_dates_20250417.rds")
nest_site <- readRDS("/Users/louisefaure/Library/CloudStorage/OneDrive-Personnel/THESE/CHAPITRE 2/git/chapter-2/DONNEES AIGLES/nest site location/nest_site_location/nest_site_location.rds")

rds_dir <- "/Users/louisefaure/Library/CloudStorage/OneDrive-Personnel/THESE/CHAPITRE 2/git/chapter-2/DONNEES AIGLES/no_burst_GE/no_burst_GE"
rds_files <- tools::file_path_sans_ext(list.files(rds_dir, pattern = "\\.rds$", ignore.case = TRUE))
rds_ids   <- as.numeric(gsub("_gpsNoDup_moveObj", "", rds_files))

# Filters individuals to retain those that have an emigration date
emig_dates_filtered <- emig_dates[emig_dates$did_disperse == TRUE &
                                    emig_dates$individual.id %in% rds_ids, ]
rds_files_filtered  <- rds_files[rds_ids %in% emig_dates_filtered$individual.id]

# Parameters
gps_dop_max <- 10
terra::terraOptions(threads = 5) 




#' -----------------------------------------------------------------------------
# STEP 1 : spatio-temporal filtering ----
#'
#' **Steps:**
#' (i) keep only GPS locations within the first 15 weeks after emigration date
#' (ii) compute height above ground using geoid and DEM corrections
#' (iii) remove GPS locations with poor positional accuracy
#' (iv) remove location at weird altitude below or above the ground


# 1.1 Loop through the files ----
dispersal_data <- lapply(rds_files_filtered, function(f) {
  
  id  <- as.numeric(gsub("_gpsNoDup_moveObj", "", f))
  obj <- readRDS(file.path(rds_dir, paste0(f, ".rds")))
  df  <- as.data.frame(obj)
  
  df <- df[, c(
    "individual.local.identifier",
    "timestamp",
    "location.long",
    "location.lat",
    "height.above.ellipsoid",
    "ground.speed",
    "eobs.speed.accuracy.estimate",
    "eobs.horizontal.accuracy.estimate",
    "gps.dop",
    "vert.speed",
    "turn.angle",
    "step.length",
    "gr.speed",
    "eobs.type.of.fix")]
  
  df$timestamp <- as.POSIXct(df$timestamp, tz = "UTC")
  
  emig_dt <- emig_dates_filtered$dispersal_date[emig_dates_filtered$individual.id == id]
  
  # 1.1.1 Retain only the first fifteen weeks of dispersal ----
  df <- df[
    df$timestamp >= emig_dt &
      df$timestamp <= emig_dt + 105 * 24 * 3600,]
  
  # Stop processing this individual if no point remains
  if (nrow(df) == 0) {return(df)}
  
  df$individual.id   <- id
  df$dispersal_date  <- emig_dt
  df$days_since_emig <- as.numeric(difftime(df$timestamp, emig_dt, units = "days"))
  
  
  # 1.1.2 Calculate height above ground ----
  xy_ll <- as.matrix(df[, c("location.long", "location.lat")])
  
  N <- terra::extract(geo, xy_ll)[, 1]
  
  # Convert ellipsoid height to height above mean sea level
  df$height_msl <- df$height.above.ellipsoid - N
  
  # Project GPS locations to the DEM coordinate system
  pts <- terra::vect(
    df[, c("location.long", "location.lat")],
    geom = c("location.long", "location.lat"),
    crs = "EPSG:4326")
  
  pts_dem <- terra::project(pts, terra::crs(dem))
  df$dem_elevation <- terra::extract(dem, pts_dem)[, 2]
  
  # Calculate height above ground
  df$height_above_ground <- df$height_msl - df$dem_elevation
  
  # 1.1.3 Remove locations with poor GPS accuracy ----
  df <- df %>%
    mutate(
      poor_gps_dop =
        !is.na(gps.dop) &
        gps.dop > gps_dop_max
    ) %>%
    filter(!poor_gps_dop)
  
  
  # Stop processing this individual if no point remains
  if (nrow(df) == 0) {return(df)}
  
  # 1.1.4 Retain only locations with available height information ----
  df <- df %>%
    filter(
      !is.na(height_above_ground),
      is.finite(height_above_ground))
  
  # Stop processing this individual if no point remains
  if (nrow(df) == 0) {return(df)}
  
  # 1.1.5 Calculate age since emigration ----
  df$age_since_emig_days <- as.numeric(
    difftime(
      df$timestamp,
      emig_dt,
      units = "days"))
  
  df$age_since_emig_weeks <- df$age_since_emig_days / 7
  
  
  # 1.1.6 Compute distance to nest site ----
  pts_3035 <- sf::st_as_sf(
    df,
    coords = c("location.long", "location.lat"),
    crs = 4326,
    remove = FALSE
  ) %>%
    sf::st_transform(3035)
  
  # Keep projected GPS coordinates
  xy_3035 <- sf::st_coordinates(pts_3035)
  
  df$x_3035 <- xy_3035[, 1]
  df$y_3035 <- xy_3035[, 2]
  
  # Select nest site for the current individual
  nest_i <- nest_site %>%
    dplyr::filter(
      as.character(individual.local.identifier) ==
        as.character(df$individual.local.identifier[1])
    )
  
  # Distance from each GPS point to the nest, in kilometres
  df$distance_to_nest_km <- as.numeric(
    sf::st_distance(
      pts_3035,
      nest_i)) / 1000
  
  df
  
})


# Combine all individuals
dispersal_data <- dplyr::bind_rows(dispersal_data)
rownames(dispersal_data) <- NULL


# 1.2 Inspect and filter extreme height-above-ground values ----
# Compile the lower 1% quantile and the observed maximum before filtering
height_limits_before_filter <- dispersal_data %>%
  summarise(
    height_q01 = quantile(
      height_above_ground,
      probs = 0.01,
      na.rm = TRUE
    ),
    height_max = max(
      height_above_ground,
      na.rm = TRUE))

print(height_limits_before_filter)


# Extract the q01 value as a numeric threshold
height_q01 <- height_limits_before_filter$height_q01


# Remove locations below q01 and above 3000 m
dispersal_data <- dispersal_data %>%
  filter(
    !is.na(height_above_ground),
    is.finite(height_above_ground),
    height_above_ground >= height_q01,
    height_above_ground <= 3000)






#' -----------------------------------------------------------------------------
# STEP 2 : EMbC behavioural classification on the full filtered dataset ----
#'
#' **Ref.** Garriga et al., 2016; Nourani et al. workflow
#'
#' (1) Classification through embc function
#' (2) Reclassify cluster 2



# 2.0 Classify behaviors ----
dispersal_data <- dispersal_data %>% arrange(individual.local.identifier, timestamp)

# Bivariate matrix: exactly two variables
behavioural_classification <- data.frame(
  ground.speed = as.numeric(dispersal_data$ground.speed),
  height_above_ground = as.numeric(dispersal_data$height_above_ground))

complete_idx <- which(
  complete.cases(behavioural_classification) &
    is.finite(behavioural_classification$ground.speed) &
    is.finite(behavioural_classification$height_above_ground))

behavioural_classification <- data.matrix(behavioural_classification[complete_idx, c("ground.speed", "height_above_ground")])

# Call EMbC
embc <- EMbC::embc(behavioural_classification)

# Smooth EMbC classification, following the published workflow
embc_smoothed <- EMbC::smth(embc, dlta = 0.7)

# Attach raw and smoothed cluster labels to the full filtered dataset
dispersal_data$behavior_cluster_raw <- NA_integer_
dispersal_data$behavior_cluster <- NA_integer_

dispersal_data$behavior_cluster_raw[complete_idx] <- embc@A
dispersal_data$behavior_cluster[complete_idx] <- embc_smoothed@A

# Check the four clusters after height filtering
cluster_summary_full <- dispersal_data %>%
  filter(!is.na(behavior_cluster)) %>%
  group_by(behavior_cluster) %>%
  summarise(
    speed_med = median(
      ground.speed,
      na.rm = TRUE
    ),
    height_med = median(
      height_above_ground,
      na.rm = TRUE
    ),
    n = n(),
    .groups = "drop"
  ) %>%
  arrange(behavior_cluster)

print(cluster_summary_full)
# behavior_cluster speed_med height_med      n
# 1                  0.11       13.8        149941
# 2                  0.78       31.5        17195
# 3                 10.4        33.4        28146
# 4                 12.7       297.         15589

#' Cluster 2 is difficult to interpret. 

# 2.1 Reclassification of cluster 2 ----
# Extract cluster 2
cluster_2 <- dispersal_data %>%
  filter(behavior_cluster == 2)

# Compile Q25, median and Q90 for speed and height
cluster_2_distribution <- cluster_2 %>%
  summarise(
    across(
      c(ground.speed, height_above_ground),
      list(
        q25 = ~ quantile(.x, 0.25, na.rm = TRUE),
        median = ~ median(.x, na.rm = TRUE),
        q90 = ~ quantile(.x, 0.90, na.rm = TRUE)),
      .names = "{.col}_{.fn}"))

print(cluster_2_distribution) 
# we can see that elevation have weird values while speed is always lower than 5 m/s.

  
# Estimate an empirical vertical-error threshold from presumed terrestrial points
# Estimate the empirical vertical-accuracy threshold from cluster 1
vertical_accuracy_max <- dispersal_data %>%
  filter(
    behavior_cluster == 1,
    is.finite(height_above_ground)
  ) %>%
  summarise(
    threshold = quantile(
      abs(height_above_ground),
      probs = 0.90,
      na.rm = TRUE
    )
  ) %>%
  pull(threshold)

# Add a 20% tolerance margin
vertical_accuracy_max <- vertical_accuracy_max * 1.20
print(vertical_accuracy_max)

# Count observations before filtering
n_before <- nrow(dispersal_data)

# Remove cluster-2 points above the empirical vertical-accuracy threshold
dispersal_data <- dispersal_data %>%
  filter(
    behavior_cluster != 2 |
      height_above_ground <= vertical_accuracy_max)

# Count and print removed GPS locations
n_removed <- n_before - nrow(dispersal_data)
cat("Number of GPS locations removed:", n_removed, "\n")

# 2.2 Attribute a name to the cluster ----
dispersal_data <- dispersal_data %>%
  mutate(
    behavior_binary = case_when(
      behavior_cluster %in% c(1, 2) ~ "terrestrial",
      behavior_cluster %in% c(3, 4) ~ "aerial",
      TRUE                          ~ NA_character_
    ),
    behavior_binary = factor(
      behavior_binary,
      levels = c("terrestrial", "aerial")))

# Check the number of locations in each behavioral category
print(table(dispersal_data$behavior_binary, useNA = "ifany"))





#'------------------------------------------------------------------------------
# Step 3. Data thining ----
#' 
#' (1) inspect raw GPS sampling intervals;
#' (2) create one high-resolution dataset at 20 min;
#' (3) create one intermediate-resolution dataset at 60 min;
#' (4) split tracks into bursts



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
  na.rm = TRUE)

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
  timing_60 %>% mutate(dataset = "60min_intermediate_resolution"))

print(timing_summary)
#' we obtain 
#' resolution_min n_individuals n_points n_bursts n_valid_transitions median_dt_next_min p05_dt_next_min p95_dt_next_min min_dt_next_min max_dt_next_min
#' 20            62              76613    12914               76613                 20            15.0            20.3              15              25                   0.0833
#' 60            66              40538     8031               40538                 60            59.4            60.5              50              70                   0.150 

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
  summarise_by_individual(regular_60_tbl))

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
  remove = FALSE)

regular_60_sf <- sf::st_as_sf(
  regular_60_tbl,
  coords = c("x_3035", "y_3035"),
  crs = 3035,
  remove = FALSE)

saveRDS(regular_20_sf, file.path("/Users/louisefaure/Library/CloudStorage/OneDrive-Personnel/THESE/CHAPITRE 2/git/chapter-2/HMM/preparation HMM/donnees intermediaire/GE_20_min_thinned.rds"))

saveRDS(regular_60_sf, file.path("/Users/louisefaure/Library/CloudStorage/OneDrive-Personnel/THESE/CHAPITRE 2/git/chapter-2/HMM/preparation HMM/donnees intermediaire/GE_60_min_thinned.rds"))





#'------------------------------------------------------------------------------
# Step 4 : extract human footprint index and control covariates ----
#'
#'** Steps:**
#' (1) evalute the horizontal accuracy to define the buffer around each GPS point
#' that should be used to extract the environmental covariates
#' (2) extract all numerical control covariates and HFI values in a 30m buffer
#' (3) for categorical covariates, such as landcover, create three categorical 
#' variables that represent the percentage of landcover in a 30m buffer for forest, 
#' grassland and bare rocks. 
#' (4) extract hfi values at two other buffer scale :
#' - hfi_mean_500m  = mean HFI of pixels within 500 m of that pixel
#' - hfi_mean_1000m = mean HFI of pixels within 1000 m of that pixel
#' (5) calculate continuous temporal covariates for the diurnal rythm
#' - the variable cos(Diel) represents diurnal (negative values) and nocturnal (positive values) periods,
#' - sin(Time) represents midnight until 11:59 am (positive values) and noon until the following 11:59 pm (negative values).
#' (6) export the RDS file


# 4.1 Evaluate the horizontal error ----

# 4.2 calculate cos(Diel) and sin(Time) ---

# 4.3 Extract environmental variables ----


# Function for covariate extraction
hfi_crs <- terra::crs(human_footprint)

extract_hfi <- function(df_sf,
                        raster    = human_footprint,
                        buffers_m = c(500, 1000)) {
  
  # Remove geometry, create a matrix for coordinates
  df  <- as.data.frame(sf::st_drop_geometry(df_sf))
  pts <- terra::vect(as.matrix(df[, c("x_3035", "y_3035")]),
                     type = "points",
                     crs  = "EPSG:3035")
  
  # Extraction of the HFI value below GPS points
  pts_r        <- terra::project(pts, hfi_crs)         # alignement sur le CRS raster
  df$hfi_point <- terra::extract(raster, pts_r, ID = FALSE)[, 1]
  
  # Buffer at 500m and 1000m, mean of the values
  for (r in buffers_m) {
    buf   <- terra::buffer(pts, width = r)             
    buf_r <- terra::project(buf, hfi_crs)
    val   <- terra::extract(raster, buf_r,
                            fun = mean, na.rm = TRUE, ID = FALSE)[, 1]
    df[[paste0("hfi_mean_", r, "m")]] <- val
  }
  
  df
}

regular_20_hfi <- extract_hfi(regular_20_sf)
regular_60_hfi <- extract_hfi(regular_60_sf)

# save and export
saveRDS(regular_20_hfi, file.path(output_dir, "GE_20_min_thinned_hfi.rds"))
saveRDS(regular_60_hfi, file.path(output_dir, "GE_60_min_thinned_hfi.rds"))
