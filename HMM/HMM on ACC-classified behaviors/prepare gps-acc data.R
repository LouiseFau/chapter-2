#' ----------------------------------------------------------------------------- 
#' Title: Preparation of acc-classified data 
#' Authors : Louise Faure
#' Date : 02.07.26
#' Purpose : 
#' (1) filter location to the first fifteen weeks of dispersal, 
#' (2) associate behaviors to one location
#' (3) thin the data at two temporal resolution (20min and 60min) and split into
#' burst,
#' (4) extract environmental covariates at the two temporal scale and within 
#' three different buffers
#' -----------------------------------------------------------------------------




library(terra)           # for raster
library(data.table)      # to read the acc-classified data
library(move)
library(dplyr)
library(sf)
library(move2)

#' Step 0 : load the data ----
output_dir <- "C:/Users/lfaure7/OneDrive/THESE/CHAPITRE 2/git/chapter-2/HMM/HMM on ACC-classified behaviors/donnes filtree intermediaire"

# Géoïde, digital elevation model, human footprint index
geo <- terra::rast("C:/Users/lfaure7/OneDrive/MEMOIRE M2/eagle_projet/data/geoide/us_nga_egm96_15.tif")
dem <- terra::rast("C:/Users/lfaure7/OneDrive/MEMOIRE M2/eagle_projet/data/pretraitements/Region-Alpes-Dem/region-alpes-dem.tif")
human_footprint <- terra::rast("C:/Users/lfaure7/OneDrive/THESE/COUCHES QGIS/COUCHES QGIS/settlements/Overture/human_footprint_index_building_pop_builtprop_100m.tif")

# emigration dates and golden eagle data
emig_dates <- readRDS("C:/Users/lfaure7/Documents/git/chapter-2/preparation donnee aigle ssf/donnees/emigration dates/emigration_dates_20250417.rds")
classified_path <- "C:/Users/lfaure7/OneDrive/THESE/CHAPITRE 2/donnees/acc-data/Louise_PhD/classified_acc_data/2024_01_24_alldata_allbirds_merged_rf_raw.csv"
new_classified_dir <- "C:/Users/lfaure7/Desktop/aigle non classifié/rf_assigned"





#'------------------------------------------------------------------------------
#'Step 1: prepare dataset and select first 15th week of the dispersal 
#'
#'**Steps:**
#'(1) merge recently classified individuals with already classified ones
#'(2) prepare dispersal date
#'(3) create age since dispersal day and retain first 15 weeks


# 1.1 load and merge the two behavioral classification ----
acc_old <- data.table::fread(
  classified_path,
  showProgress = TRUE)

acc_old[, source_dataset := "previous_classification"]
acc_old[, source_file := basename(classified_path)]

new_files <- list.files(
  path = new_classified_dir,
  pattern = "^rf_ss_.*\\.csv$",
  full.names = TRUE)


acc_new_list <- lapply(new_files, function(f) {
  
  x <- data.table::fread(f, showProgress = FALSE)
  
  x[, source_dataset := "new_classification"]
  x[, source_file := basename(f)]
  
  x
})

acc_new <- data.table::rbindlist(
  acc_new_list,
  use.names = TRUE,
  fill = TRUE
)


# Merge
acc_classified_all <- data.table::rbindlist(
  list(acc_old, acc_new),
  use.names = TRUE,
  fill = TRUE
)

acc_classified_all[, individualID := as.character(individualID)]
acc_classified_all[, rf8fitted := as.character(rf8fitted)]
acc_classified_all[, timestamp := as.POSIXct(timestamp, tz = "UTC")]

acc_classified_all <- acc_classified_all[!is.na(timestamp)]

# 1.2 Prepare dispersal dates ----
emig_dt <- data.table::as.data.table(emig_dates)

emig_dt[, individual.local.identifier := as.character(individual.local.identifier)]
emig_dt[, dispersal_date := as.POSIXct(dispersal_date, tz = "UTC")]

emig_filtered <- emig_dt[
  did_disperse == TRUE & !is.na(dispersal_date),
  .(
    individual.local.identifier,
    dispersal_date
  )
]

emig_filtered <- unique(
  emig_filtered,
  by = "individual.local.identifier"
)


# 1.3 Join dispersal dates and retain first 15 weeks ----
acc_classified_all[
  emig_filtered,
  dispersal_date := i.dispersal_date,
  on = c("individualID" = "individual.local.identifier")
]

acc_classified_all <- acc_classified_all[!is.na(dispersal_date)]

acc_classified_all[
  ,
  age_days := as.numeric(
    difftime(timestamp, dispersal_date, units = "days")
  )]

acc_classified_all[, age_weeks := age_days / 7]

acc_15w <- acc_classified_all[
  age_days >= 0 &
    age_days < 105] # first fifteen weeks of dispersal

acc_15w[, timespan := data.table::as.IDate(timestamp)]

saveRDS(
  acc_15w,
  file.path(output_dir, "acc_classified_first_15_weeks_all_individuals_merged_old_new.rds"),
  compress = "gzip"
)



#' -----------------------------------------------------------------------------
#' Step 2: assign ACC-classified bursts to the nearest GPS location ----
#'
#' **Purpose:** Associate each ACC-classified burst with the temporally nearest 
#' GPS location from the same individual.
#'
#' Rule:
#' (1) for each ACC burst, find the nearest GPS timestamp from the same individual;
#' (2) compute the absolute temporal difference between ACC burst and GPS location;
#' (3) keep the assignment only if the GPS location is within 60 min;
#' (4) ACC bursts with no GPS location within 60 min are kept as diagnostics but
#'     excluded from spatial analyses;
#' (5) when several ACC bursts are assigned to the same GPS location, retain the
#'     closest ACC burst as the behavior attached to that GPS point.



#' Step 2.0: load GPS files without burst segmentation ----
gps_dir <- "C:/Users/lfaure7/OneDrive/THESE/CHAPITRE 2/git/chapter-2/preparation donnee aigle ssf/donnees/no_burst_GE/no_burst_GE"

gps_files <- list.files(
  path = gps_dir,
  pattern = "\\.rds$",
  full.names = TRUE,
  recursive = TRUE,
  ignore.case = TRUE
)


read_one_gps_move <- function(f) {
  
  move_obj <- readRDS(f)
  
  # Extract GPS relocation table
  x <- data.table::as.data.table(move_obj)
  
  # Extract individual metadata stored in the Move object's idData slot
  id_data <- data.table::as.data.table(
    methods::slot(move_obj, "idData")
  )
  
  # Add individual identity to every GPS row
  x[, individual.local.identifier := as.character(id_data$individual.local.identifier[1])]
  x[, tag.local.identifier := as.integer(id_data$tag.local.identifier[1])]
  x[, source_gps_file := basename(f)]
  
  # Ensure timestamp is POSIXct
  x[, timestamp := as.POSIXct(timestamp, tz = "UTC")]
  
  x
}


gps_by_file <- lapply(gps_files, read_one_gps_move)

names(gps_by_file) <- make.unique(
  vapply(
    gps_by_file,
    function(x) x$individual.local.identifier[1],
    character(1)
  )
)



#' Step 2.1: merge GPS files and retain the first 15 weeks of dispersal ----
gps_raw_all <- data.table::rbindlist(
  gps_by_file,
  use.names = TRUE,
  fill = TRUE
)

gps_raw_all[
  emig_filtered,
  dispersal_date := i.dispersal_date,
  on = "individual.local.identifier"
]

gps_raw_all <- gps_raw_all[!is.na(dispersal_date)]

gps_raw_all[
  ,
  age_days := as.numeric(
    difftime(timestamp, dispersal_date, units = "days")
  )
]

gps_raw_all[, age_weeks := age_days / 7]

gps_15w <- gps_raw_all[
  age_days >= 0 &
    age_days < 105 &
    !is.na(location.long) &
    !is.na(location.lat)
]

gps_15w[, timespan := data.table::as.IDate(timestamp)]

# Keep only GPS individuals that also have ACC-classified behaviors
gps_15w <- gps_15w[
  individual.local.identifier %in% unique(acc_15w$individualID)
]

data.table::setorder(
  gps_15w,
  individual.local.identifier,
  timestamp
)

# Remove exact duplicated GPS records.
# This removes repeated records only when individual, timestamp and event.id are identical.
gps_15w <- unique(
  gps_15w,
  by = c("individual.local.identifier", "timestamp", "event.id")
)

# Create one unique GPS row identifier after all GPS filtering
gps_15w[, gps_row_id := .I]



#' Step 2.2: prepare GPS points for nearest-neighbour temporal join ----
gps_points <- gps_15w[
  ,
  .(
    individual.local.identifier,
    tag.local.identifier,
    source_gps_file,
    gps_row_id,
    gps_event_id = event.id,
    gps_timestamp = timestamp,
    timespan,
    location.long,
    location.lat,
    height.above.ellipsoid,
    ground.speed,
    heading,
    dispersal_date,
    age_days,
    age_weeks
  )
]

gps_points[, join_time := gps_timestamp]

data.table::setorder(
  gps_points,
  individual.local.identifier,
  join_time
)


#' Step 2.3: prepare ACC-classified bursts ----
acc_points <- acc_15w[
  ,
  .(
    individual.local.identifier = individualID,
    acc_event_id = event.id,
    acc_burstID = burstID,
    acc_timestamp = timestamp,
    rf8fitted,
    pro_rf8fitted,
    odbaAvg,
    rollanimaltrack,
    pitchanimaltrack,
    burstmeanx,
    burstmeany,
    burstmeanz
  )
]

acc_points <- acc_points[
  !is.na(individual.local.identifier) &
    !is.na(acc_timestamp) &
    !is.na(rf8fitted)
]

acc_points[, join_time := acc_timestamp]

data.table::setorder(
  acc_points,
  individual.local.identifier,
  join_time
)


#' Step 2.4: assign each ACC burst to the nearest GPS location ----
max_assignment_gap_min <- 60

data.table::setkey(
  gps_points,
  individual.local.identifier,
  join_time
)

data.table::setkey(
  acc_points,
  individual.local.identifier,
  join_time
)

# For each ACC burst, find the nearest GPS point from the same individual.
# The resulting table has one row per ACC burst.
acc_nearest_gps <- gps_points[
  acc_points,
  on = .(individual.local.identifier, join_time),
  roll = "nearest"
]

acc_nearest_gps[
  ,
  abs_time_diff_min := abs(
    as.numeric(
      difftime(acc_timestamp, gps_timestamp, units = "mins")
    )
  )
]

acc_nearest_gps[
  ,
  gps_assigned_60min := !is.na(gps_row_id) &
    abs_time_diff_min <= max_assignment_gap_min
]

acc_with_gps_60min <- acc_nearest_gps[
  gps_assigned_60min == TRUE
]

acc_without_gps_60min <- acc_nearest_gps[
  is.na(gps_assigned_60min) |
    gps_assigned_60min == FALSE
]



#' Step 2.5: reduce to one behavior per GPS location ----
#'
#' Several ACC bursts can be assigned to the same GPS point.
#' For a GPS-level trajectory, retain the ACC burst closest in time to the GPS fix.
#' If two ACC bursts are equally close, retain the one with the highest RF probability.

acc_with_gps_60min[
  ,
  pro_rf8fitted_order := data.table::fifelse(
    is.na(pro_rf8fitted),
    -Inf,
    pro_rf8fitted
  )
]

data.table::setorder(
  acc_with_gps_60min,
  gps_row_id,
  abs_time_diff_min,
  -pro_rf8fitted_order
)

gps_behavior_from_acc <- acc_with_gps_60min[
  ,
  .SD[1],
  by = gps_row_id
]

gps_behavior_from_acc <- gps_behavior_from_acc[
  ,
  .(
    gps_row_id,
    rf8fitted,
    pro_rf8fitted,
    acc_event_id,
    acc_burstID,
    acc_timestamp,
    abs_time_diff_min,
    odbaAvg,
    rollanimaltrack,
    pitchanimaltrack,
    burstmeanx,
    burstmeany,
    burstmeanz
  )
]


#' Step 2.7: merge assigned behavior back to GPS locations ----
gps_beh_15w <- merge(
  gps_points[
    ,
    .(
      individual.local.identifier,
      tag.local.identifier,
      source_gps_file,
      gps_row_id,
      gps_event_id,
      gps_timestamp,
      timespan,
      location.long,
      location.lat,
      height.above.ellipsoid,
      ground.speed,
      heading,
      dispersal_date,
      age_days,
      age_weeks
    )
  ],
  gps_behavior_from_acc,
  by = "gps_row_id",
  all.x = TRUE
)

gps_beh_15w[
  ,
  behavior_assigned := !is.na(rf8fitted)
]

data.table::setorder(
  gps_beh_15w,
  individual.local.identifier,
  gps_timestamp
)


#' Step 2.8: final GPS-level diagnostics ----
gps_behavior_assignment_summary <- gps_beh_15w[
  ,
  .(
    n_gps = .N,
    n_gps_with_behavior = sum(behavior_assigned),
    prop_gps_with_behavior = mean(behavior_assigned),
    median_abs_time_diff_min = median(abs_time_diff_min, na.rm = TRUE)
  ),
  by = individual.local.identifier
][order(individual.local.identifier)]

gps_behavior_assignment_summary


#' Step 2.9: save outputs ----
saveRDS(
  gps_beh_15w,
  file.path(output_dir, "gps_with_nearest_acc_behavior_within_60min.rds"),
  compress = "gzip"
)



#' -----------------------------------------------------------------------------
#' Step 3 : creation of two temporally regular dataset
#' 
#' **Steps:**
#' (1) inspect raw GPS sampling intervals;
#' (2) create one high-resolution dataset at 20 min;
#' (3) create one intermediate-resolution dataset at 60 min;
#' (4) split tracks into bursts


# 3.1 Inspect raw GPS data ----
# Create move2 object from dispersal_data
step3_input <- gps_beh_15w[behavior_assigned == TRUE]

locs <- step3_input %>%
  mutate(
    individual.id = as.character(individual.local.identifier),
    timestamp = as.POSIXct(gps_timestamp, tz = "UTC"),
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
#' 0.133    14.9    19.76     43.9    20.05 24839.40
#' We choose a fine scale resolution of 20 minutes because it is close to the 
#' upper quartile of the main fine-scale sampling distribution and 60 minutes 
#' because 95% quantiles are around 60 minutes.


#' Step 3.2 Thin GPS data at 20-min and 60-min intervals ----
high_res_min <- 20
intermediate_res_min <- 60

high_res_tolerance_min <- 5
intermediate_tolerance_min <- 10

min_points_per_burst <- 3

format_interval_name <- function(interval_min) {
  gsub("\\.", "p", paste0(interval_min, "min"))
}


#' Thin GPS data
thin_move2_interval <- function(locs,
                                interval_unit,
                                resolution_min,
                                resolution_label) {
  
  locs_thin <- locs %>%
    dplyr::arrange(individual.id, timestamp) %>%
    move2::mt_filter_per_interval(
      unit = interval_unit,
      criterion = "first"
    )
  
  xy <- sf::st_coordinates(locs_thin)
  
  locs_thin_tbl <- locs_thin %>%
    dplyr::mutate(
      individual.id = as.character(individual.id),
      individual.local.identifier = as.character(individual.id),
      timestamp = as.POSIXct(timestamp, tz = "UTC"),
      x_3035 = as.numeric(xy[, 1]),
      y_3035 = as.numeric(xy[, 2]),
      resolution_min = resolution_min,
      resolution_class = resolution_label
    ) %>%
    sf::st_drop_geometry() %>%
    as.data.frame() %>%
    dplyr::arrange(individual.id, timestamp)
  
  locs_thin_tbl
}


#' Step 3.3 Split retained tracks into regular bursts ----
split_into_regular_bursts <- function(df,
                                      interval_min,
                                      tolerance_min,
                                      min_points_per_burst = 3) {
  
  df_burst <- df %>%
    dplyr::arrange(individual.id, timestamp) %>%
    dplyr::group_by(individual.id) %>%
    dplyr::mutate(
      dt_prev_min = as.numeric(
        difftime(timestamp, dplyr::lag(timestamp), units = "mins")
      ),
      
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
    dplyr::ungroup()
  
  df_burst <- df_burst %>%
    dplyr::group_by(individual.id, burst_id) %>%
    dplyr::mutate(
      n_points_burst = dplyr::n()
    ) %>%
    dplyr::ungroup() %>%
    dplyr::filter(n_points_burst >= min_points_per_burst)
  
  df_burst <- df_burst %>%
    dplyr::arrange(individual.id, burst_id, timestamp) %>%
    dplyr::group_by(individual.id, burst_id) %>%
    dplyr::mutate(
      row_in_burst = dplyr::row_number(),
      
      dt_prev_burst_min = as.numeric(
        difftime(timestamp, dplyr::lag(timestamp), units = "mins")
      ),
      
      dt_next_burst_min = as.numeric(
        difftime(dplyr::lead(timestamp), timestamp, units = "mins")
      ),
      
      has_next_regular = !is.na(dt_next_burst_min) &
        abs(dt_next_burst_min - interval_min) <= tolerance_min,
      
      lag_deviation_next_min = dt_next_burst_min - interval_min
    ) %>%
    dplyr::ungroup()
  
  df_burst
}


#' Step 3.4 Create 20-min and 60-min datasets ----
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



#' Step 3.5 Control the filtered datasets ----
check_timing <- function(df) {
  
  df %>%
    dplyr::filter(has_next_regular == TRUE) %>%
    dplyr::summarise(
      resolution_min = dplyr::first(resolution_min),
      n_individuals = dplyr::n_distinct(individual.local.identifier),
      n_points = dplyr::n(),
      n_bursts = dplyr::n_distinct(burst_id),
      n_valid_transitions = sum(has_next_regular, na.rm = TRUE),
      median_dt_next_min = median(dt_next_burst_min, na.rm = TRUE),
      p05_dt_next_min = stats::quantile(dt_next_burst_min, 0.05, na.rm = TRUE),
      p95_dt_next_min = stats::quantile(dt_next_burst_min, 0.95, na.rm = TRUE),
      min_dt_next_min = min(dt_next_burst_min, na.rm = TRUE),
      max_dt_next_min = max(dt_next_burst_min, na.rm = TRUE),
      median_abs_deviation_min = median(abs(lag_deviation_next_min), na.rm = TRUE),
      p95_abs_deviation_min = stats::quantile(abs(lag_deviation_next_min), 0.95, na.rm = TRUE)
    )
}

timing_20 <- check_timing(regular_20_tbl)
timing_60 <- check_timing(regular_60_tbl)

timing_summary <- dplyr::bind_rows(
  timing_20 %>% dplyr::mutate(dataset = "20min_high_resolution"),
  timing_60 %>% dplyr::mutate(dataset = "60min_intermediate_resolution")
)

print(timing_summary)


#' Step 3.5b Inspect individual differences ----
summarise_by_individual <- function(df) {
  
  df %>%
    dplyr::group_by(resolution_class, individual.local.identifier) %>%
    dplyr::summarise(
      n_points_retained = dplyr::n(),
      n_bursts = dplyr::n_distinct(burst_id),
      n_valid_transitions = sum(has_next_regular, na.rm = TRUE),
      median_points_per_burst = median(n_points_burst, na.rm = TRUE),
      max_points_per_burst = max(n_points_burst, na.rm = TRUE),
      .groups = "drop"
    )
}

summary_by_id <- dplyr::bind_rows(
  summarise_by_individual(regular_20_tbl),
  summarise_by_individual(regular_60_tbl)
)

print(summary_by_id, n = 128)

low_data_individuals <- summary_by_id %>%
  dplyr::filter(n_valid_transitions < 200) %>%
  dplyr::arrange(resolution_class, n_valid_transitions)

print(low_data_individuals)


#' Step 3.6 Export datasets ----
regular_20_tbl <- regular_20_tbl %>%
  dplyr::mutate(
    x_3035 = as.numeric(x_3035),
    y_3035 = as.numeric(y_3035)
  ) %>%
  dplyr::filter(
    !is.na(x_3035),
    !is.na(y_3035),
    is.finite(x_3035),
    is.finite(y_3035)
  )

regular_60_tbl <- regular_60_tbl %>%
  dplyr::mutate(
    x_3035 = as.numeric(x_3035),
    y_3035 = as.numeric(y_3035)
  ) %>%
  dplyr::filter(
    !is.na(x_3035),
    !is.na(y_3035),
    is.finite(x_3035),
    is.finite(y_3035)
  )

regular_20_tbl %>%
  dplyr::summarise(
    n_rows = dplyr::n(),
    n_na_x = sum(is.na(x_3035)),
    n_na_y = sum(is.na(y_3035)),
    n_non_finite_x = sum(!is.finite(x_3035)),
    n_non_finite_y = sum(!is.finite(y_3035))
  )

regular_60_tbl %>%
  dplyr::summarise(
    n_rows = dplyr::n(),
    n_na_x = sum(is.na(x_3035)),
    n_na_y = sum(is.na(y_3035)),
    n_non_finite_x = sum(!is.finite(x_3035)),
    n_non_finite_y = sum(!is.finite(y_3035))
  )


make_sf_from_xy <- function(df, crs_value = 3035) {
  
  df <- as.data.frame(df)
  
  x <- as.numeric(df$x_3035)
  y <- as.numeric(df$y_3035)
  
  keep <- !is.na(x) & !is.na(y) & is.finite(x) & is.finite(y)
  
  df <- df[keep, , drop = FALSE]
  x <- x[keep]
  y <- y[keep]
  
  geometry <- sf::st_sfc(
    lapply(
      seq_along(x),
      function(i) sf::st_point(c(x[i], y[i]))
    ),
    crs = crs_value
  )
  
  sf::st_sf(
    df,
    geometry = geometry
  )
}

regular_20_sf <- make_sf_from_xy(regular_20_tbl, crs_value = 3035)
regular_60_sf <- make_sf_from_xy(regular_60_tbl, crs_value = 3035)

sf::st_crs(regular_20_sf)
sf::st_crs(regular_60_sf)

sum(sf::st_is_empty(regular_20_sf))
sum(sf::st_is_empty(regular_60_sf))


#' Step 3.7 Save outputs ----
saveRDS(regular_20_sf,
  file.path(output_dir, "GE_20_min_thinned_behavior_assigned.rds"),
  compress = "gzip")

saveRDS(regular_60_sf,
  file.path(output_dir, "GE_60_min_thinned_behavior_assigned.rds"),
  compress = "gzip")



#'------------------------------------------------------------------------------
#' Step 4 : extract human footprint index around GPS locations
#'
#' (1) extract covariates values at three buffer size : 
#' - hfi_mean_100m  = HFI value of the pixel containing the GPS point
#' - hfi_mean_500m  = mean HFI of pixels within 500 m of that pixel
#' - hfi_mean_1000m = mean HFI of pixels within 1000 m of that pixel
#' (2) export the RDS file


# Terra parameters
terra::terraOptions(threads = 5)

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
