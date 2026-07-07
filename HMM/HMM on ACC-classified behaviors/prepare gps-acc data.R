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
#' Step 2: associate ACC-classified behaviors to GPS locations ----
#'
#' **Steps:**
#' (1) load GPS files without burst segmentation;
#' (2) retain GPS locations within the first 15 weeks of dispersal;
#' (3) inspect GPS and ACC temporal gaps;
#' (4) define two empirical temporal windows (10 min and 20 min);
#' (5) assign ACC-classified bursts to GPS locations within these windows;
#' (6) retain the majority behavior in each window;
#' (7) keep a nearest-ACC-within-60-min fallback as conservative backup.



#' Step 2.0: load GPS files without burst segmentation ----
gps_dir <- "C:/Users/lfaure7/OneDrive/THESE/CHAPITRE 2/git/chapter-2/preparation donnee aigle ssf/donnees/no_burst_GE/no_burst_GE"

gps_files <- list.files(
  path = gps_dir,
  pattern = "\\.rds$",
  full.names = TRUE,
  recursive = TRUE,
  ignore.case = TRUE
)

# Read one Move object and add individual metadata stored in idData
read_one_gps_move <- function(f) {
  
  move_obj <- readRDS(f)
  
  # Extract GPS relocation table
  x <- data.table::as.data.table(move_obj)
  
  # Extract individual metadata from the Move object's idData slot
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

gps_by_individual <- lapply(gps_files, read_one_gps_move)

names(gps_by_individual) <- make.unique(
  vapply(
    gps_by_individual,
    function(x) x$individual.local.identifier[1],
    character(1)
  )
)


#' Step 2.1: retain GPS locations within the first 15 weeks of dispersal ----
gps_by_individual_15w <- lapply(gps_by_individual, function(x) {
  
  x[
    emig_filtered,
    dispersal_date := i.dispersal_date,
    on = "individual.local.identifier"
  ]
  
  x <- x[!is.na(dispersal_date)]
  
  x[
    ,
    age_days := as.numeric(
      difftime(timestamp, dispersal_date, units = "days")
    )
  ]
  
  x[, age_weeks := age_days / 7]
  
  x <- x[
    age_days >= 0 &
      age_days < 105 &
      !is.na(location.long) &
      !is.na(location.lat)
  ]
  
  x[, timespan := data.table::as.IDate(timestamp)]
  
  data.table::setorder(x, timestamp)
  
  x
})

# Remove individuals with no GPS location in the 15-week window
gps_by_individual_15w <- gps_by_individual_15w[
  vapply(gps_by_individual_15w, nrow, integer(1)) > 0
]

gps_15w <- data.table::rbindlist(
  gps_by_individual_15w,
  use.names = TRUE,
  fill = TRUE
)

# Keep only GPS individuals that also have ACC-classified behaviors
gps_15w <- gps_15w[
  individual.local.identifier %in% unique(acc_15w$individualID)
]

data.table::setorder(
  gps_15w,
  individual.local.identifier,
  timestamp
)

gps_15w <- unique(
  gps_15w,
  by = c("individual.local.identifier", "timestamp", "event.id")
)

gps_15w[, gps_row_id := .I]


#' Step 2.2: inspect GPS and ACC temporal gaps ----
gps_gap_dt <- gps_15w[
  ,
  .(
    timestamp,
    next_timestamp = data.table::shift(timestamp, type = "lead")
  ),
  by = individual.local.identifier
]

gps_gap_dt[
  ,
  gps_gap_min := as.numeric(
    difftime(next_timestamp, timestamp, units = "mins")
  )
]

gps_gap_dt <- gps_gap_dt[
  !is.na(gps_gap_min) &
    gps_gap_min > 0
]

acc_ts <- unique(
  acc_15w[
    !is.na(individualID) &
      !is.na(timestamp),
    .(
      individualID,
      timestamp
    )
  ]
)

data.table::setorder(acc_ts, individualID, timestamp)

acc_gap_dt <- acc_ts[
  ,
  .(
    timestamp,
    next_timestamp = data.table::shift(timestamp, type = "lead")
  ),
  by = individualID
]

acc_gap_dt[
  ,
  acc_gap_min := as.numeric(
    difftime(next_timestamp, timestamp, units = "mins")
  )
]

acc_gap_dt <- acc_gap_dt[
  !is.na(acc_gap_min) &
    acc_gap_min > 0
]

# table 
summarise_time_gaps <- function(dt, gap_col, data_type) {
  
  x <- dt[[gap_col]]
  x <- x[!is.na(x) & x > 0]
  
  data.table::data.table(
    data_type = data_type,
    n_gaps = length(x),
    q25_gap_min = as.numeric(quantile(x, 0.25)),
    median_gap_min = median(x),
    q75_gap_min = as.numeric(quantile(x, 0.75)),
    q95_gap_min = as.numeric(quantile(x, 0.95)),
    prop_gap_le_10 = mean(x <= 10),
    prop_gap_le_20 = mean(x <= 20),
    prop_gap_gt_60 = mean(x > 60)
  )
}

gap_summary_compact <- data.table::rbindlist(
  list(
    summarise_time_gaps(gps_gap_dt, "gps_gap_min", "GPS"),
    summarise_time_gaps(acc_gap_dt, "acc_gap_min", "ACC")
  )
)

gap_summary_compact

#' we obtain 
#' data_type n_gaps q25_gap_min median_gap_min q75_gap_min q95_gap_min prop_gap_le_10 prop_gap_le_20 prop_gap_gt_60
#' GPS 231706   10.333333       15.30002   20.033333    60.46667      0.2431832      0.7126056       0.070982193
#' ACC 667189    4.983333        5.00000    5.183333    15.00000      0.8986044      0.9884321       0.008965975


#' Step 2.4: define empirical assignment thresholds ----
# fine_window_min = strict behavior assignment window around dense GPS fixes.
# regular_window_min = main assignment window matching the dominant GPS interval.
# max_assignment_gap_min = maximum accepted temporal distance in case of sparse GPS.

fine_window_min <- 10
regular_window_min <- 20
max_assignment_gap_min <- 60

#' Step 2.5: create GPS temporal structure ----
data.table::setorder(
  gps_15w,
  individual.local.identifier,
  timestamp
)

gps_15w[
  ,
  previous_gps_timestamp := data.table::shift(timestamp, type = "lag"),
  by = individual.local.identifier
]

gps_15w[
  ,
  next_gps_timestamp := data.table::shift(timestamp, type = "lead"),
  by = individual.local.identifier
]

gps_15w[
  ,
  dt_previous_min := as.numeric(
    difftime(timestamp, previous_gps_timestamp, units = "mins")
  )
]

gps_15w[
  ,
  dt_next_min := as.numeric(
    difftime(next_gps_timestamp, timestamp, units = "mins")
  )
]

# GPS burst identifiers at both temporal scales
gps_15w[
  ,
  new_gps_burst_10min := is.na(dt_previous_min) |
    dt_previous_min > fine_window_min,
  by = individual.local.identifier
]

gps_15w[
  ,
  gps_burst_id_10min := cumsum(new_gps_burst_10min),
  by = individual.local.identifier
]

gps_15w[
  ,
  new_gps_burst_20min := is.na(dt_previous_min) |
    dt_previous_min > regular_window_min,
  by = individual.local.identifier
]

gps_15w[
  ,
  gps_burst_id_20min := cumsum(new_gps_burst_20min),
  by = individual.local.identifier
]


#' Step 2.6: build allocation windows for a given temporal scale ----

build_gps_windows <- function(gps_dt, window_min, scale_name, max_gap_min = 60) {
  
  w <- data.table::copy(gps_dt)
  
  half_window_min <- window_min / 2
  
  w[
    ,
    before_midpoint_limit_min := data.table::fifelse(
      is.na(dt_previous_min),
      max_gap_min,
      dt_previous_min / 2
    )
  ]
  
  w[
    ,
    after_midpoint_limit_min := data.table::fifelse(
      is.na(dt_next_min),
      max_gap_min,
      dt_next_min / 2
    )
  ]
  
  # The window is centred on the GPS timestamp,
  # but it cannot cross the midpoint with neighbouring GPS fixes.
  w[
    ,
    allocation_before_min := pmin(
      half_window_min,
      before_midpoint_limit_min,
      max_gap_min,
      na.rm = TRUE
    )
  ]
  
  w[
    ,
    allocation_after_min := pmin(
      half_window_min,
      after_midpoint_limit_min,
      max_gap_min,
      na.rm = TRUE
    )
  ]
  
  w[
    ,
    dispersal_end := dispersal_date + 105 * 24 * 60 * 60
  ]
  
  w[
    ,
    allocation_start := pmax(
      timestamp - allocation_before_min * 60,
      dispersal_date
    )
  ]
  
  w[
    ,
    allocation_end := pmin(
      timestamp + allocation_after_min * 60,
      dispersal_end
    )
  ]
  
  w[
    ,
    allocation_window_duration_min := as.numeric(
      difftime(allocation_end, allocation_start, units = "mins")
    )
  ]
  
  w[, scale_name := scale_name]
  w[, window_min := window_min]
  
  w[
    !is.na(allocation_start) &
      !is.na(allocation_end) &
      allocation_start <= allocation_end
  ]
}

#' Step 2.7: create 10-min, 20-min, and 60-min windows ----

gps_windows_10 <- build_gps_windows(
  gps_dt = gps_15w,
  window_min = fine_window_min,
  scale_name = "fine_10min",
  max_gap_min = max_assignment_gap_min
)

gps_windows_20 <- build_gps_windows(
  gps_dt = gps_15w,
  window_min = regular_window_min,
  scale_name = "regular_20min",
  max_gap_min = max_assignment_gap_min
)

# 60-min nearest window:
# here window_min = 120 gives a maximum of 60 min before and 60 min after,
# while still preventing the window from crossing midpoints with neighbouring GPS fixes.
gps_windows_60 <- build_gps_windows(
  gps_dt = gps_15w,
  window_min = 2 * max_assignment_gap_min,
  scale_name = "nearest_60min",
  max_gap_min = max_assignment_gap_min
)

gps_windows_multiscale <- data.table::rbindlist(
  list(gps_windows_10, gps_windows_20),
  use.names = TRUE,
  fill = TRUE
)

#' Step 2.8: prepare ACC-classified bursts as temporal points ----
acc_points <- acc_15w[
  ,
  .(
    individual.local.identifier = individualID,
    acc_timestamp = timestamp,
    acc_start = timestamp,
    acc_end = timestamp,
    acc_event_id = event.id,
    rf8fitted,
    pro_rf8fitted
  )
]

acc_points <- acc_points[
  !is.na(individual.local.identifier) &
    !is.na(acc_timestamp) &
    !is.na(rf8fitted)
]


#' Step 2.9: assign ACC bursts to 10-min and 20-min GPS windows ----
data.table::setkey(
  gps_windows_multiscale,
  individual.local.identifier,
  allocation_start,
  allocation_end
)

data.table::setkey(
  acc_points,
  individual.local.identifier,
  acc_start,
  acc_end
)

acc_to_gps_multiscale <- data.table::foverlaps(
  x = acc_points,
  y = gps_windows_multiscale,
  by.x = c("individual.local.identifier", "acc_start", "acc_end"),
  by.y = c("individual.local.identifier", "allocation_start", "allocation_end"),
  type = "within",
  nomatch = 0L
)

acc_to_gps_multiscale[
  ,
  abs_time_diff_min := abs(
    as.numeric(
      difftime(acc_timestamp, timestamp, units = "mins")
    )
  )
]

#' Step 2.10: compute majority behavior within each temporal window ----
#'
#' Because behavior is categorical, we do not compute a numeric mean of labels.
#' We compute behavior proportions and retain the majority behavior.

behavior_scale_detail <- acc_to_gps_multiscale[
  ,
  .(
    n_acc_bursts_behavior = .N,
    mean_pro_rf8fitted = mean(pro_rf8fitted, na.rm = TRUE),
    median_pro_rf8fitted = median(pro_rf8fitted, na.rm = TRUE),
    mean_abs_time_diff_min = mean(abs_time_diff_min, na.rm = TRUE),
    min_abs_time_diff_min = min(abs_time_diff_min, na.rm = TRUE)
  ),
  by = .(
    gps_row_id,
    scale_name,
    window_min,
    rf8fitted
  )
]

behavior_scale_detail[
  ,
  n_acc_bursts_scale := sum(n_acc_bursts_behavior),
  by = .(
    gps_row_id,
    scale_name
  )
]

behavior_scale_detail[
  ,
  prop_behavior := n_acc_bursts_behavior / n_acc_bursts_scale
]

data.table::setorder(
  behavior_scale_detail,
  gps_row_id,
  scale_name,
  -n_acc_bursts_behavior,
  -prop_behavior,
  -mean_pro_rf8fitted,
  min_abs_time_diff_min
)

behavior_scale_majority <- behavior_scale_detail[
  ,
  .SD[1],
  by = .(
    gps_row_id,
    scale_name
  )
]

behavior_scale_majority_wide <- data.table::dcast(
  behavior_scale_majority,
  gps_row_id ~ scale_name,
  value.var = c(
    "rf8fitted",
    "prop_behavior",
    "n_acc_bursts_scale",
    "mean_pro_rf8fitted",
    "median_pro_rf8fitted",
    "mean_abs_time_diff_min",
    "min_abs_time_diff_min"
  )
)

#' Step 2.11: compute nearest ACC behavior within 60 min ----
#'
#' This is not used to average behavior.
#' It is a fallback / diagnostic for GPS locations without ACC behavior
#' in the 10-min or 20-min windows.

data.table::setkey(
  gps_windows_60,
  individual.local.identifier,
  allocation_start,
  allocation_end
)

acc_to_gps_60 <- data.table::foverlaps(
  x = acc_points,
  y = gps_windows_60,
  by.x = c("individual.local.identifier", "acc_start", "acc_end"),
  by.y = c("individual.local.identifier", "allocation_start", "allocation_end"),
  type = "within",
  nomatch = 0L
)

acc_to_gps_60[
  ,
  abs_time_diff_min := abs(
    as.numeric(
      difftime(acc_timestamp, timestamp, units = "mins")
    )
  )
]

data.table::setorder(
  acc_to_gps_60,
  gps_row_id,
  abs_time_diff_min,
  -pro_rf8fitted
)

gps_behavior_nearest_60 <- acc_to_gps_60[
  ,
  .SD[1],
  by = gps_row_id
]

gps_behavior_nearest_60 <- gps_behavior_nearest_60[
  ,
  .(
    gps_row_id,
    rf8fitted_nearest_60min = rf8fitted,
    pro_rf8fitted_nearest_60min = pro_rf8fitted,
    nearest_acc_timestamp_60min = acc_timestamp,
    nearest_acc_event_id_60min = acc_event_id,
    abs_time_diff_nearest_60min = abs_time_diff_min
  )
]

#' Step 2.12: merge behavioral assignments back to GPS locations ----

gps_beh_15w <- merge(
  gps_15w,
  behavior_scale_majority_wide,
  by = "gps_row_id",
  all.x = TRUE
)

gps_beh_15w <- merge(
  gps_beh_15w,
  gps_behavior_nearest_60,
  by = "gps_row_id",
  all.x = TRUE
)

# Recommended final behavior:
# use the 20-min majority behavior when available;
# if not available, use nearest ACC behavior within 60 min.
gps_beh_15w[
  ,
  rf8fitted_final := data.table::fcoalesce(
    rf8fitted_regular_20min,
    rf8fitted_nearest_60min
  )
]

gps_beh_15w[
  ,
  behavior_assigned := !is.na(rf8fitted_final)
]

gps_beh_15w[
  ,
  assignment_rule := data.table::fcase(
    !is.na(rf8fitted_regular_20min),
    "regular_20min_majority",
    
    is.na(rf8fitted_regular_20min) & !is.na(rf8fitted_nearest_60min),
    "nearest_60min_fallback",
    
    default = NA_character_
  )
]


#' Step 2.13: keep useful columns for next steps ----

gps_beh_15w <- gps_beh_15w[
  ,
  .(
    individual.local.identifier,
    tag.local.identifier,
    source_gps_file,
    gps_row_id,
    event.id,
    timestamp,
    timespan,
    location.long,
    location.lat,
    height.above.ellipsoid,
    ground.speed,
    heading,
    dispersal_date,
    age_days,
    age_weeks,
    gps_burst_id_10min,
    gps_burst_id_20min,
    dt_previous_min,
    dt_next_min,
    behavior_assigned,
    assignment_rule,
    rf8fitted_final,
    rf8fitted_fine_10min,
    rf8fitted_regular_20min,
    rf8fitted_nearest_60min,
    prop_behavior_fine_10min,
    prop_behavior_regular_20min,
    n_acc_bursts_scale_fine_10min,
    n_acc_bursts_scale_regular_20min,
    mean_pro_rf8fitted_fine_10min,
    mean_pro_rf8fitted_regular_20min,
    median_pro_rf8fitted_fine_10min,
    median_pro_rf8fitted_regular_20min,
    mean_abs_time_diff_min_fine_10min,
    mean_abs_time_diff_min_regular_20min,
    min_abs_time_diff_min_fine_10min,
    min_abs_time_diff_min_regular_20min,
    pro_rf8fitted_nearest_60min,
    abs_time_diff_nearest_60min,
    nearest_acc_timestamp_60min,
    nearest_acc_event_id_60min
  )
]

data.table::setorder(
  gps_beh_15w,
  individual.local.identifier,
  timestamp
)




