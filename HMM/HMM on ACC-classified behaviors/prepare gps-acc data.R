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
    median_abs_time_diff_min = median(abs_time_diff_min, na.rm = TRUE),
    q95_abs_time_diff_min = safe_q95(abs_time_diff_min)
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


