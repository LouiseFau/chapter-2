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
#' Step 2: associate acc-classified behaviors to gps location
#' **Steps:**
#' (1) estimate the typical GPS interval for each individual;
#' (2) create allocation windows around each GPS location;
#' (3) assign ACC-classified bursts to GPS locations according to these windows;
#' (4) retain the majority behavior per GPS location.


gps_dir <- "C:/Users/lfaure7/OneDrive/THESE/CHAPITRE 2/git/chapter-2/preparation donnee aigle ssf/donnees/no_burst_GE/no_burst_GE"



# Path to GPS raw files.
# This folder must contain the GPS files for all individuals, old and new.
gps_dir <- "C:/Users/lfaure7/Desktop/dossier_gps_a_completer"

gps_files <- list.files(
  path = gps_dir,
  pattern = "\\.csv$",
  full.names = TRUE
)

# Load GPS data
gps_raw <- data.table::rbindlist(
  lapply(gps_files, data.table::fread, showProgress = FALSE),
  use.names = TRUE,
  fill = TRUE
)

# Standardise GPS columns
gps_dt <- gps_raw[
  ,
  .(
    individualID = as.character(individual.local.identifier),
    gps_timestamp = as.POSIXct(timestamp, tz = "UTC"),
    gps_event_id = event.id,
    longitude = as.numeric(location.long),
    latitude = as.numeric(location.lat)
  )
]

# Retain valid GPS locations
gps_dt <- gps_dt[
  !is.na(individualID) &
    !is.na(gps_timestamp) &
    !is.na(longitude) &
    !is.na(latitude)
]

# Remove duplicated GPS points if present
gps_dt <- unique(
  gps_dt,
  by = c("individualID", "gps_timestamp", "longitude", "latitude")
)

# Add dispersal date to GPS locations
gps_dt[
  emig_filtered,
  dispersal_date := i.dispersal_date,
  on = c("individualID" = "individual.local.identifier")
]

gps_dt <- gps_dt[!is.na(dispersal_date)]

# Create GPS age since dispersal
gps_dt[
  ,
  age_days := as.numeric(
    difftime(gps_timestamp, dispersal_date, units = "days")
  )
]

gps_dt[, age_weeks := age_days / 7]

# Retain GPS locations in the first 15 weeks of dispersal
gps_15w <- gps_dt[
  age_days >= 0 &
    age_days < 105
]

gps_15w[, timespan := data.table::as.IDate(gps_timestamp)]

# Keep only individuals that also have ACC-classified behaviors
gps_15w <- gps_15w[
  individualID %in% unique(acc_15w$individualID)
]

# Sort GPS locations
data.table::setorder(gps_15w, individualID, gps_timestamp)

# Estimate GPS interval per individual
gps_interval_summary <- gps_15w[
  ,
  {
    dt_min <- as.numeric(diff(gps_timestamp), units = "mins")
    
    .(
      n_gps = .N,
      median_gps_interval_min = median(dt_min[dt_min > 0], na.rm = TRUE),
      p90_gps_interval_min = as.numeric(
        quantile(dt_min[dt_min > 0], probs = 0.90, na.rm = TRUE)
      )
    )
  },
  by = individualID
]

# Add GPS interval to each GPS location
gps_15w[
  gps_interval_summary,
  median_gps_interval_min := i.median_gps_interval_min,
  on = "individualID"
]

# Maximum tolerated allocation window.
# If GPS fixes are missing, we do not allocate ACC bursts farther than:
# min(individual median GPS interval, 60 minutes).
gps_15w[
  ,
  allocation_tolerance_min := pmin(median_gps_interval_min, 60)
]

# Previous and next GPS timestamps
gps_15w[
  ,
  previous_gps_timestamp := data.table::shift(gps_timestamp, type = "lag"),
  by = individualID
]

gps_15w[
  ,
  next_gps_timestamp := data.table::shift(gps_timestamp, type = "lead"),
  by = individualID
]

gps_15w[
  ,
  dt_previous_min := as.numeric(
    difftime(gps_timestamp, previous_gps_timestamp, units = "mins")
  )
]

gps_15w[
  ,
  dt_next_min := as.numeric(
    difftime(next_gps_timestamp, gps_timestamp, units = "mins")
  )
]

# Allocation window before and after each GPS point.
# For regular GPS intervals, this cuts at the temporal midpoint.
# For long intervals caused by missing GPS fixes, the window is capped.
gps_15w[
  ,
  allocation_before_min := data.table::fifelse(
    is.na(dt_previous_min),
    allocation_tolerance_min,
    pmin(dt_previous_min / 2, allocation_tolerance_min)
  )
]

gps_15w[
  ,
  allocation_after_min := data.table::fifelse(
    is.na(dt_next_min),
    allocation_tolerance_min,
    pmin(dt_next_min / 2, allocation_tolerance_min)
  )
]

gps_15w[
  ,
  allocation_start := gps_timestamp - allocation_before_min * 60
]

gps_15w[
  ,
  allocation_end := gps_timestamp + allocation_after_min * 60
]

gps_15w[, gps_row_id := .I]

# Prepare GPS windows for temporal overlap
gps_windows <- gps_15w[
  !is.na(allocation_start) &
    !is.na(allocation_end) &
    allocation_start <= allocation_end
]

# Prepare ACC-classified bursts as temporal points
acc_points <- acc_15w[
  ,
  .(
    individualID,
    acc_timestamp = timestamp,
    acc_start = timestamp,
    acc_end = timestamp,
    acc_event_id = event.id,
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
  !is.na(individualID) &
    !is.na(acc_timestamp) &
    !is.na(rf8fitted)
]

# Temporal overlap: assign each ACC burst to a GPS allocation window
data.table::setkey(
  gps_windows,
  individualID,
  allocation_start,
  allocation_end
)

data.table::setkey(
  acc_points,
  individualID,
  acc_start,
  acc_end
)

acc_to_gps <- data.table::foverlaps(
  x = acc_points,
  y = gps_windows,
  by.x = c("individualID", "acc_start", "acc_end"),
  by.y = c("individualID", "allocation_start", "allocation_end"),
  type = "within",
  nomatch = 0L
)

# Count behaviors assigned to each GPS location
beh_by_gps <- acc_to_gps[
  ,
  .(
    n_acc_bursts_behavior = .N,
    mean_pro_rf8fitted = mean(pro_rf8fitted, na.rm = TRUE),
    median_pro_rf8fitted = median(pro_rf8fitted, na.rm = TRUE),
    first_acc_timestamp = min(acc_timestamp, na.rm = TRUE),
    last_acc_timestamp = max(acc_timestamp, na.rm = TRUE)
  ),
  by = .(
    gps_row_id,
    individualID,
    gps_timestamp,
    rf8fitted
  )
]

# Select majority behavior per GPS location.
# Tie-breaking rule: highest mean Random Forest probability.
data.table::setorder(
  beh_by_gps,
  gps_row_id,
  -n_acc_bursts_behavior,
  -mean_pro_rf8fitted
)

gps_behavior_majority <- beh_by_gps[
  ,
  .SD[1],
  by = gps_row_id
]

data.table::setnames(
  gps_behavior_majority,
  old = "rf8fitted",
  new = "rf8fitted_majority"
)

# Merge majority behavior back to GPS locations
gps_beh_15w <- merge(
  gps_windows,
  gps_behavior_majority[
    ,
    .(
      gps_row_id,
      rf8fitted_majority,
      n_acc_bursts_behavior,
      mean_pro_rf8fitted,
      median_pro_rf8fitted,
      first_acc_timestamp,
      last_acc_timestamp
    )
  ],
  by = "gps_row_id",
  all.x = TRUE
)

gps_beh_15w[
  ,
  behavior_assigned := !is.na(rf8fitted_majority)
]

# Rename for clarity
data.table::setnames(
  gps_beh_15w,
  old = "rf8fitted_majority",
  new = "rf8fitted"
)

# Keep useful columns for the next steps
gps_beh_15w <- gps_beh_15w[
  ,
  .(
    individualID,
    gps_row_id,
    gps_event_id,
    gps_timestamp,
    timespan,
    longitude,
    latitude,
    dispersal_date,
    age_days,
    age_weeks,
    median_gps_interval_min,
    allocation_tolerance_min,
    allocation_start,
    allocation_end,
    behavior_assigned,
    rf8fitted,
    n_acc_bursts_behavior,
    mean_pro_rf8fitted,
    median_pro_rf8fitted,
    first_acc_timestamp,
    last_acc_timestamp
  )
]

data.table::setorder(
  gps_beh_15w,
  individualID,
  gps_timestamp
)