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

# Load GPS data 
output_dir <- "C:/Users/lfaure7/Documents/git/chapter-2/preparation donnee aigle ssf/resultats-intermediaire"

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
#   (iv) split tracks so transitions do not cross long gaps.
################################################################################

library(data.table)

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


# A. Inspect raw time lags ----
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

high_res_min <- 20
intermediate_res_min <- 60
high_res_tolerance_min <- 5
intermediate_tolerance_min <- 10

# Long-gap rule:
# A raw temporal gap larger than 3 times the target interval starts a new segment.
high_res_max_gap_min <- high_res_min * 3          # 60 min
intermediate_max_gap_min <- intermediate_res_min * 3  # 180 min

min_points_per_segment <- 3

# B. Function: strict temporal regularisation ----

format_interval_name <- function(interval_min) {
  gsub("\\.", "p", paste0(interval_min, "min"))
}

make_regular_dataset <- function(data,
                                 interval_min,
                                 tolerance_min,
                                 max_gap_min,
                                 min_points_per_segment = 3) {
  
  required_cols <- c(
    "individual.id", "timestamp", "x_3035", "y_3035"
  )
  
  stopifnot(all(required_cols %in% names(data)))
  
  dt <- as.data.table(data)
  
  dt <- dt[
    !is.na(individual.id) &
      !is.na(timestamp) &
      !is.na(x_3035) &
      !is.na(y_3035)
  ]
  
  dt[, individual.id := as.character(individual.id)]
  setorder(dt, individual.id, timestamp)
  
  # Keep original row identifier
  dt[, row_id_original := .I]
  
  # Calculate raw temporal gaps
  dt[, raw_dt_min := as.numeric(
    difftime(timestamp, shift(timestamp), units = "mins")
  ), by = individual.id]
  
  # Create technical segments at long gaps.
  dt[, new_raw_segment := is.na(raw_dt_min) | raw_dt_min > max_gap_min]
  dt[, raw_segment_n := cumsum(new_raw_segment), by = individual.id]
  dt[, raw_segment_id := paste(
    individual.id,
    sprintf("%04d", raw_segment_n),
    sep = "_"
  )]
  
  # Build one regular timestamp grid per raw segment.
  # The grid is anchored at the first timestamp of the segment.
  # This avoids problems with arbitrary clock alignment.
  step_sec <- interval_min * 60
  
  segment_ranges <- dt[, .(
    start_time = min(timestamp),
    end_time   = max(timestamp),
    n_raw      = .N
  ), by = .(individual.id, raw_segment_id)]
  
  grid <- segment_ranges[, {
    grid_times <- seq(
      from = start_time[1],
      to   = end_time[1],
      by   = paste(step_sec, "secs")
    )
    
    .(timestamp_regular = grid_times)
    
  }, by = .(individual.id, raw_segment_id)]
  
  if (nrow(grid) == 0) {
    diagnostics <- tibble(
      interval_min = interval_min,
      tolerance_min = tolerance_min,
      max_gap_min = max_gap_min,
      n_raw = nrow(dt),
      n_grid = 0L,
      n_matched = 0L,
      grid_coverage = NA_real_,
      n_ids_raw = dplyr::n_distinct(dt$individual.id),
      n_ids_matched = 0L,
      n_valid_transitions = 0L,
      n_regular_segments = 0L,
      median_time_error_min = NA_real_,
      p95_time_error_min = NA_real_
    )
    
    return(list(
      data = tibble(),
      diagnostics = diagnostics,
      diagnostics_by_id = tibble()
    ))
  }
  
  grid[, grid_row_id := .I]
  grid[, time_join := timestamp_regular]
  
  # Prepare observed fixes for nearest-time matching
  dt[, timestamp_obs := timestamp]
  dt[, time_join := timestamp]
  
  setkey(dt, individual.id, raw_segment_id, time_join)
  setkey(grid, individual.id, raw_segment_id, time_join)
  
  # Match each regular timestamp to the nearest observed GPS fix
  matched <- dt[
    grid,
    on = .(individual.id, raw_segment_id, time_join),
    roll = "nearest"
  ]
  
  matched[, time_error_min := abs(as.numeric(
    difftime(timestamp_obs, timestamp_regular, units = "mins")
  ))]
  
  # Keep only observed fixes close enough to the target regular timestamp
  matched <- matched[
    !is.na(timestamp_obs) &
      time_error_min <= tolerance_min
  ]
  
  if (nrow(matched) == 0) {
    diagnostics <- tibble(
      interval_min = interval_min,
      tolerance_min = tolerance_min,
      max_gap_min = max_gap_min,
      n_raw = nrow(dt),
      n_grid = nrow(grid),
      n_matched = 0L,
      grid_coverage = 0,
      n_ids_raw = dplyr::n_distinct(dt$individual.id),
      n_ids_matched = 0L,
      n_valid_transitions = 0L,
      n_regular_segments = 0L,
      median_time_error_min = NA_real_,
      p95_time_error_min = NA_real_
    )
    
    return(list(
      data = tibble(),
      diagnostics = diagnostics,
      diagnostics_by_id = tibble()
    ))
  }
  
  # Enforce one GPS fix per regular timestamp
  setorder(matched, individual.id, raw_segment_id,
           timestamp_regular, time_error_min)
  
  matched <- unique(
    matched,
    by = c("individual.id", "raw_segment_id", "timestamp_regular")
  )
  
  # Enforce one regular timestamp per GPS fix
  setorder(matched, individual.id, raw_segment_id,
           row_id_original, time_error_min)
  
  matched <- unique(
    matched,
    by = c("individual.id", "raw_segment_id", "row_id_original")
  )
  
  # Split again when regular grid cells are missing.
  # This guarantees that transitions do not cross missing regular fixes.
  setorder(matched, individual.id, raw_segment_id, timestamp_regular)
  
  matched[, dt_regular_prev_min := as.numeric(
    difftime(timestamp_regular, shift(timestamp_regular), units = "mins")
  ), by = .(individual.id, raw_segment_id)]
  
  matched[, new_regular_segment :=
            is.na(dt_regular_prev_min) |
            abs(dt_regular_prev_min - interval_min) > 1e-6]
  
  matched[, regular_segment_n := cumsum(new_regular_segment),
          by = .(individual.id, raw_segment_id)]
  
  matched[, regular_segment_id := paste(
    raw_segment_id,
    sprintf("%04d", regular_segment_n),
    sep = "_"
  )]
  
  matched[, n_regular_segment := .N,
          by = .(individual.id, regular_segment_id)]
  
  matched <- matched[n_regular_segment >= min_points_per_segment]
  
  if (nrow(matched) == 0) {
    diagnostics <- tibble(
      interval_min = interval_min,
      tolerance_min = tolerance_min,
      max_gap_min = max_gap_min,
      n_raw = nrow(dt),
      n_grid = nrow(grid),
      n_matched = 0L,
      grid_coverage = 0,
      n_ids_raw = dplyr::n_distinct(dt$individual.id),
      n_ids_matched = 0L,
      n_valid_transitions = 0L,
      n_regular_segments = 0L,
      median_time_error_min = NA_real_,
      p95_time_error_min = NA_real_
    )
    
    return(list(
      data = tibble(),
      diagnostics = diagnostics,
      diagnostics_by_id = tibble()
    ))
  }
  
  # Flag valid next-step transitions
  setorder(matched, individual.id, regular_segment_id, timestamp_regular)
  
  matched[, dt_next_regular_min := as.numeric(
    difftime(
      shift(timestamp_regular, type = "lead"),
      timestamp_regular,
      units = "mins"
    )
  ), by = .(individual.id, regular_segment_id)]
  
  matched[, has_next_regular :=
            !is.na(dt_next_regular_min) &
            abs(dt_next_regular_min - interval_min) < 1e-6]
  
  matched[, resolution_min := interval_min]
  matched[, tolerance_min := tolerance_min]
  matched[, max_gap_min := max_gap_min]
  
  # Diagnostics by individual
  grid_id <- grid[, .(
    n_grid = .N
  ), by = individual.id]
  
  matched_id <- matched[, .(
    n_matched = .N,
    n_valid_transitions = sum(has_next_regular, na.rm = TRUE),
    n_regular_segments = uniqueN(regular_segment_id),
    median_time_error_min = median(time_error_min, na.rm = TRUE),
    p95_time_error_min = as.numeric(
      quantile(time_error_min, probs = 0.95, na.rm = TRUE)
    )
  ), by = individual.id]
  
  diagnostics_by_id <- merge(
    grid_id,
    matched_id,
    by = "individual.id",
    all.x = TRUE
  )
  
  diagnostics_by_id[is.na(n_matched), n_matched := 0L]
  diagnostics_by_id[is.na(n_valid_transitions), n_valid_transitions := 0L]
  diagnostics_by_id[is.na(n_regular_segments), n_regular_segments := 0L]
  
  diagnostics_by_id[, grid_coverage := n_matched / n_grid]
  
  diagnostics <- tibble(
    interval_min = interval_min,
    tolerance_min = tolerance_min,
    max_gap_min = max_gap_min,
    n_raw = nrow(dt),
    n_grid = nrow(grid),
    n_matched = nrow(matched),
    grid_coverage = nrow(matched) / nrow(grid),
    n_ids_raw = dplyr::n_distinct(dt$individual.id),
    n_ids_matched = dplyr::n_distinct(matched$individual.id),
    n_valid_transitions = matched[, sum(has_next_regular, na.rm = TRUE)],
    n_regular_segments = matched[, uniqueN(regular_segment_id)],
    median_time_error_min = median(matched$time_error_min, na.rm = TRUE),
    p95_time_error_min = as.numeric(
      quantile(matched$time_error_min, probs = 0.95, na.rm = TRUE)
    )
  )
  
  list(
    data = as_tibble(matched),
    diagnostics = diagnostics,
    diagnostics_by_id = as_tibble(diagnostics_by_id)
  )
}

################################################################################
# 2.4 Create 20-min and 60-min regular datasets
################################################################################

regular_20 <- make_regular_dataset(
  data = locs_tbl,
  interval_min = high_res_min,
  tolerance_min = high_res_tolerance_min,
  max_gap_min = high_res_max_gap_min,
  min_points_per_segment = min_points_per_segment
)

regular_60 <- make_regular_dataset(
  data = locs_tbl,
  interval_min = intermediate_res_min,
  tolerance_min = intermediate_tolerance_min,
  max_gap_min = intermediate_max_gap_min,
  min_points_per_segment = min_points_per_segment
)

regular_20_tbl <- regular_20$data %>%
  mutate(resolution_class = "20min_high_resolution")

regular_60_tbl <- regular_60$data %>%
  mutate(resolution_class = "60min_intermediate_resolution")

regularisation_diagnostics <- bind_rows(
  regular_20$diagnostics %>%
    mutate(dataset = "20min_high_resolution"),
  regular_60$diagnostics %>%
    mutate(dataset = "60min_intermediate_resolution")
)

print(regularisation_diagnostics)

readr::write_csv(
  regularisation_diagnostics,
  file.path(output_dir, "regularisation_diagnostics_20min_60min.csv")
)

readr::write_csv(
  regular_20$diagnostics_by_id,
  file.path(output_dir, "regularisation_diagnostics_by_id_20min.csv")
)

readr::write_csv(
  regular_60$diagnostics_by_id,
  file.path(output_dir, "regularisation_diagnostics_by_id_60min.csv")
)

################################################################################
# 2.5 Recompute movement variables after temporal regularisation
################################################################################

add_regular_movement_metrics <- function(df) {
  
  df %>%
    arrange(individual.id, regular_segment_id, timestamp_regular) %>%
    group_by(individual.id, regular_segment_id) %>%
    mutate(
      x_prev = lag(x_3035),
      y_prev = lag(y_3035),
      x_next = lead(x_3035),
      y_next = lead(y_3035),
      
      step_length_m_reg = sqrt(
        (x_3035 - x_prev)^2 +
          (y_3035 - y_prev)^2
      ),
      
      ground_speed_m_s_reg = step_length_m_reg / (resolution_min * 60),
      
      bearing_in = atan2(y_3035 - y_prev, x_3035 - x_prev),
      bearing_out = atan2(y_next - y_3035, x_next - x_3035),
      
      turn_angle_reg = atan2(
        sin(bearing_out - bearing_in),
        cos(bearing_out - bearing_in)
      ),
      
      vertical_speed_m_s_reg = (
        height_above_ground - lag(height_above_ground)
      ) / (resolution_min * 60)
    ) %>%
    ungroup()
}

regular_20_tbl <- add_regular_movement_metrics(regular_20_tbl)
regular_60_tbl <- add_regular_movement_metrics(regular_60_tbl)


################################################################################
# 2.6 Convert back to sf objects and save
################################################################################

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
  file.path(
    output_dir,
    paste0("golden_eagle_regular_", format_interval_name(high_res_min), ".rds")
  )
)

saveRDS(
  regular_60_sf,
  file.path(
    output_dir,
    paste0("golden_eagle_regular_", format_interval_name(intermediate_res_min), ".rds")
  )
)

################################################################################
# 2.7 Final checks
################################################################################

check_regular_lags <- function(df) {
  
  df %>%
    arrange(individual.id, regular_segment_id, timestamp_regular) %>%
    group_by(individual.id, regular_segment_id) %>%
    mutate(
      dt_check_min = as.numeric(
        difftime(timestamp_regular, lag(timestamp_regular), units = "mins")
      )
    ) %>%
    ungroup() %>%
    filter(!is.na(dt_check_min)) %>%
    summarise(
      resolution_min = first(resolution_min),
      n_steps = n(),
      min_dt = min(dt_check_min, na.rm = TRUE),
      median_dt = median(dt_check_min, na.rm = TRUE),
      max_dt = max(dt_check_min, na.rm = TRUE),
      n_non_regular = sum(abs(dt_check_min - resolution_min) > 1e-6),
      n_valid_transitions = sum(has_next_regular, na.rm = TRUE),
      n_segments = n_distinct(regular_segment_id),
      n_individuals = n_distinct(individual.id)
    )
}

check_20 <- check_regular_lags(regular_20_tbl)
check_60 <- check_regular_lags(regular_60_tbl)

print(check_20)
print(check_60)

# Distribution of time-matching error
regular_20_tbl %>%
  summarise(
    median_time_error_min = median(time_error_min, na.rm = TRUE),
    p95_time_error_min = quantile(time_error_min, 0.95, na.rm = TRUE),
    max_time_error_min = max(time_error_min, na.rm = TRUE)
  )

regular_60_tbl %>%
  summarise(
    median_time_error_min = median(time_error_min, na.rm = TRUE),
    p95_time_error_min = quantile(time_error_min, 0.95, na.rm = TRUE),
    max_time_error_min = max(time_error_min, na.rm = TRUE)
  )
