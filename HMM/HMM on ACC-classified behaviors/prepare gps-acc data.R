#' ----------------------------------------------------------------------------- 
#' Title: Preparation of acc-classified data 
#' Authors : Louise Faure
#' Date : 02.07.26
#' Purpose : 
#' (1) filter location to the first fifteen weeks of dispersal, 
#' (2) associate acc-burst (with behaviors) to one location
#' (3) reclassify behaviors before thining
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
output_dir <- "/Users/louisefaure/Library/CloudStorage/OneDrive-Personnel/THESE/CHAPITRE 2/git/chapter-2/HMM/HMM on ACC-classified behaviors/donnees intermediaire (2)"

# emigration dates, golden eagle data, human footprint index
human_footprint <- terra::rast("/Users/louisefaure/Library/CloudStorage/OneDrive-Personnel/THESE/COUCHES QGIS/COUCHES QGIS/settlements/Overture/human_footprint_index_building_pop_builtprop_100m.tif")
emig_dates <- readRDS("/Users/louisefaure/Library/CloudStorage/OneDrive-Personnel/THESE/CHAPITRE 2/git/chapter-2/DONNEES AIGLES/emigration dates/emigration_dates_20250417.rds")
classified_path <- "/Users/louisefaure/Library/CloudStorage/OneDrive-Personnel/THESE/CHAPITRE 2/git/chapter-2/DONNEES AIGLES/classified_acc_data/2024_01_24_alldata_allbirds_merged_rf_raw.csv"
new_classified_dir <- "/Users/louisefaure/Library/CloudStorage/OneDrive-Personnel/THESE/CHAPITRE 2/git/chapter-2/DONNEES AIGLES/Individus non classifies/rf_assigned"





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
  
  x})

acc_new <- data.table::rbindlist(
  acc_new_list,
  use.names = TRUE,
  fill = TRUE)


# Merge
acc_classified_all <- data.table::rbindlist(
  list(acc_old, acc_new),
  use.names = TRUE,
  fill = TRUE)

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
  )]

emig_filtered <- unique(
  emig_filtered,
  by = "individual.local.identifier")


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
    difftime(timestamp, dispersal_date, units = "days"))]

acc_classified_all[, age_weeks := age_days / 7]

acc_15w <- acc_classified_all[
  age_days >= 0 &
    age_days < 105] # first fifteen weeks of dispersal

acc_15w[, timespan := data.table::as.IDate(timestamp)]

saveRDS(
  acc_15w,
  file.path(output_dir, "acc_classified_first_15_weeks_all_individuals_merged.rds"),
  compress = "gzip"
) 
#' Droslöng17 (eobs 5704) and Viluoch17 (eobs 4570) does not have ACC data for
#' the first fifteen weeks of the dispersal period. N. of individuals : 64




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


#'------------------------------------------------------------------------------
# STEP 2: assign ACC-classified behaviour to GPS locations ----
#'------------------------------------------------------------------------------

library(data.table)
library(dplyr)
library(sf)


#'------------------------------------------------------------------------------
# 2.1 Load the two datasets ----
#'------------------------------------------------------------------------------

gps <- readRDS(
  "/Users/louisefaure/Desktop/dossier sans titre/donnees aigles gps burst/gps_clean_by_individual.rds"
)

acc <- readRDS(
  paste0(
    "/Users/louisefaure/Library/CloudStorage/OneDrive-Personnel/",
    "THESE/CHAPITRE 2/git/chapter-2/HMM/",
    "HMM on ACC-classified behaviors/donnees intermediaire (2)/",
    "acc_classified_first_15_weeks_all_individuals_merged.rds"
  )
)


# Maximum accepted temporal difference between an ACC record
# and a GPS location.

max_assignment_gap_min <- 60


#'------------------------------------------------------------------------------
# 2.2 Combine the list of GPS objects while preserving geometry ----
#'------------------------------------------------------------------------------

# The code also works if gps is accidentally a single sf or move2 object.

if (
  inherits(gps, "sf") ||
  inherits(gps, "move2")
) {
  gps_list <- list(gps)
} else {
  gps_list <- gps
}


if (
  length(gps_list) == 0L
) {
  stop("The GPS list is empty.")
}


# Check that all GPS objects use the same CRS.

reference_crs <- sf::st_crs(
  gps_list[[1]]
)

same_crs <- vapply(
  gps_list,
  function(x) {
    isTRUE(
      sf::st_crs(x) ==
        reference_crs
    )
  },
  logical(1)
)

if (
  !all(same_crs)
) {
  stop(
    paste0(
      "The GPS objects do not all use the same CRS. ",
      "Reproject them before combining the list."
    )
  )
}


# Create names for unnamed list elements.

gps_list_names <- names(
  gps_list
)

if (
  is.null(gps_list_names)
) {
  gps_list_names <- rep(
    "",
    length(
      gps_list
    )
  )
}

missing_list_names <-
  is.na(gps_list_names) |
  gps_list_names == ""

gps_list_names[
  missing_list_names
] <- paste0(
  "gps_object_",
  which(
    missing_list_names
  )
)


# Convert every object to sf and retain the original list-element name.

gps_sf_list <- lapply(
  seq_along(
    gps_list
  ),
  function(i) {
    
    x <- sf::st_as_sf(
      gps_list[[i]]
    )
    
    x$source_gps_object <-
      gps_list_names[[i]]
    
    x
  }
)


gps_sf <- dplyr::bind_rows(
  gps_sf_list
)

gps_sf <- sf::st_as_sf(
  gps_sf
)


#'------------------------------------------------------------------------------
# 2.3 Harmonise GPS column names and timestamps ----
#'------------------------------------------------------------------------------

# Support both Movebank naming conventions:
# individual_local_identifier and individual.local.identifier.

if (
  !"individual_local_identifier" %in%
  names(gps_sf) &&
  "individual.local.identifier" %in%
  names(gps_sf)
) {
  gps_sf <- gps_sf %>%
    dplyr::rename(
      individual_local_identifier =
        individual.local.identifier
    )
}


if (
  !"individual_local_identifier" %in%
  names(gps_sf)
) {
  stop(
    paste0(
      "The GPS data do not contain ",
      "'individual_local_identifier'."
    )
  )
}


if (
  !"timestamp" %in%
  names(gps_sf)
) {
  stop(
    "The GPS data do not contain a timestamp column."
  )
}


gps_sf <- gps_sf %>%
  dplyr::mutate(
    
    individual_local_identifier =
      trimws(
        as.character(
          individual_local_identifier
        )
      ),
    
    # This assumes that the timestamps represent UTC.
    # Change tz if the original timestamps were recorded differently.
    gps_timestamp =
      as.POSIXct(
        timestamp,
        tz = "UTC"
      ),
    
    # Stable identifier used to merge the behaviour back
    # without disturbing the geometry.
    gps_row_id =
      dplyr::row_number()
  )


if (
  anyNA(
    gps_sf$gps_timestamp
  )
) {
  warning(
    "Some GPS timestamps could not be converted to POSIXct."
  )
}


#'------------------------------------------------------------------------------
# 2.4 Prepare the ACC-classified behaviour data ----
#'------------------------------------------------------------------------------

acc_dt <- data.table::as.data.table(
  data.table::copy(
    acc
  )
)


required_acc_columns <- c(
  "individualID",
  "timestamp",
  "rf8fitted"
)

missing_acc_columns <- setdiff(
  required_acc_columns,
  names(
    acc_dt
  )
)

if (
  length(
    missing_acc_columns
  ) > 0L
) {
  stop(
    paste0(
      "The following required ACC columns are missing: ",
      paste(
        missing_acc_columns,
        collapse = ", "
      )
    )
  )
}


acc_dt[
  ,
  individualID :=
    trimws(
      as.character(
        individualID
      )
    )
]

acc_dt[
  ,
  acc_timestamp :=
    as.POSIXct(
      timestamp,
      tz = "UTC"
    )
]

acc_dt[
  ,
  rf8fitted :=
    as.character(
      rf8fitted
    )
]

acc_dt[
  ,
  acc_row_id :=
    .I
]


# Keep the required variables and any optional ACC variables
# that are present in the dataset.

optional_acc_columns <- intersect(
  c(
    "event.id",
    "burstID",
    "pro_rf8fitted",
    "odbaAvg",
    "rollanimaltrack",
    "pitchanimaltrack",
    "burstmeanx",
    "burstmeany",
    "burstmeanz"
  ),
  names(
    acc_dt
  )
)

acc_columns_to_keep <- c(
  "acc_row_id",
  "individualID",
  "acc_timestamp",
  "rf8fitted",
  optional_acc_columns
)

acc_points <- acc_dt[
  ,
  ..acc_columns_to_keep
]


# Rename the individual identifier to match the GPS data.

data.table::setnames(
  acc_points,
  old =
    "individualID",
  new =
    "individual_local_identifier"
)


# Rename optional identifying columns when available.

if (
  "event.id" %in%
  names(acc_points)
) {
  data.table::setnames(
    acc_points,
    old =
      "event.id",
    new =
      "acc_event_id"
  )
}

if (
  "burstID" %in%
  names(acc_points)
) {
  data.table::setnames(
    acc_points,
    old =
      "burstID",
    new =
      "acc_burstID"
  )
}


# Remove ACC rows that cannot be used for matching.

acc_points <- acc_points[
  !is.na(
    individual_local_identifier
  ) &
    individual_local_identifier != "" &
    !is.na(
      acc_timestamp
    ) &
    !is.na(
      rf8fitted
    )
]


# Create a numeric time variable for the rolling join.
# Numeric Unix time avoids inconsistencies caused only by displayed time zones.

acc_points[
  ,
  join_time :=
    as.numeric(
      acc_timestamp
    )
]


#'------------------------------------------------------------------------------
# 2.5 Compare individual names between datasets ----
#'------------------------------------------------------------------------------

gps_individuals <- sort(
  unique(
    gps_sf$
      individual_local_identifier
  )
)

acc_individuals <- sort(
  unique(
    acc_points$
      individual_local_identifier
  )
)


acc_individuals_without_gps <- setdiff(
  acc_individuals,
  gps_individuals
)

gps_individuals_without_acc <- setdiff(
  gps_individuals,
  acc_individuals
)


individual_name_diagnostics <- list(
  
  n_gps_individuals =
    length(
      gps_individuals
    ),
  
  n_acc_individuals =
    length(
      acc_individuals
    ),
  
  acc_individuals_without_gps =
    acc_individuals_without_gps,
  
  gps_individuals_without_acc =
    gps_individuals_without_acc
)


print(
  individual_name_diagnostics
)


# According to your description, this object should contain:
# "Droslöng17 (eobs 5704)"
# "Viluoch17 (eobs 4570)"

print(
  acc_individuals_without_gps
)


#'------------------------------------------------------------------------------
# 2.6 Prepare a non-spatial GPS table for the temporal join ----
#'------------------------------------------------------------------------------

gps_points <- data.table::data.table(
  
  gps_row_id =
    gps_sf$gps_row_id,
  
  individual_local_identifier =
    gps_sf$
    individual_local_identifier,
  
  gps_timestamp =
    gps_sf$gps_timestamp
)


gps_points[
  ,
  join_time :=
    as.numeric(
      gps_timestamp
    )
]


gps_points <- gps_points[
  !is.na(
    individual_local_identifier
  ) &
    individual_local_identifier != "" &
    !is.na(
      gps_timestamp
    )
]


data.table::setkey(
  acc_points,
  individual_local_identifier,
  join_time
)

data.table::setkey(
  gps_points,
  individual_local_identifier,
  join_time
)


#------------------------------------------------------------------------------
# Assign to each GPS point the closest ACC behaviour in time
#------------------------------------------------------------------------------

max_assignment_gap_min <- 60


# Ensure that timestamps are POSIXct and use the same time zone.

gps_points[
  ,
  gps_timestamp :=
    as.POSIXct(
      gps_timestamp,
      tz = "UTC"
    )
]

acc_points[
  ,
  acc_timestamp :=
    as.POSIXct(
      acc_timestamp,
      tz = "UTC"
    )
]


# Use numeric timestamps for the rolling joins.

gps_points[
  ,
  join_time :=
    as.numeric(
      gps_timestamp
    )
]

acc_points[
  ,
  join_time :=
    as.numeric(
      acc_timestamp
    )
]


# Each GPS point must have a unique identifier.

if (
  anyDuplicated(
    gps_points$gps_row_id
  ) > 0L
) {
  stop(
    "gps_row_id is not unique in gps_points."
  )
}


data.table::setkey(
  acc_points,
  individual_local_identifier,
  join_time
)

data.table::setkey(
  gps_points,
  individual_local_identifier,
  join_time
)


#------------------------------------------------------------------------------
# 1. Closest ACC record occurring before the GPS timestamp
#------------------------------------------------------------------------------

acc_previous <- acc_points[
  gps_points,
  
  on = .(
    individual_local_identifier,
    join_time
  ),
  
  roll = Inf,
  
  # Prevent duplicate ACC timestamps from producing several rows.
  mult = "last",
  
  .(
    gps_row_id =
      i.gps_row_id,
    
    individual_local_identifier =
      i.individual_local_identifier,
    
    gps_timestamp =
      i.gps_timestamp,
    
    acc_row_id,
    
    acc_timestamp,
    
    rf8fitted,
    
    pro_rf8fitted,
    
    acc_event_id,
    
    acc_burstID,
    
    odbaAvg,
    
    rollanimaltrack,
    
    pitchanimaltrack,
    
    burstmeanx,
    
    burstmeany,
    
    burstmeanz
  )
]


acc_previous[
  ,
  temporal_direction :=
    "previous"
]


#------------------------------------------------------------------------------
# 2. Closest ACC record occurring after the GPS timestamp
#------------------------------------------------------------------------------

acc_next <- acc_points[
  gps_points,
  
  on = .(
    individual_local_identifier,
    join_time
  ),
  
  roll = -Inf,
  
  # Prevent duplicate ACC timestamps from producing several rows.
  mult = "first",
  
  .(
    gps_row_id =
      i.gps_row_id,
    
    individual_local_identifier =
      i.individual_local_identifier,
    
    gps_timestamp =
      i.gps_timestamp,
    
    acc_row_id,
    
    acc_timestamp,
    
    rf8fitted,
    
    pro_rf8fitted,
    
    acc_event_id,
    
    acc_burstID,
    
    odbaAvg,
    
    rollanimaltrack,
    
    pitchanimaltrack,
    
    burstmeanx,
    
    burstmeany,
    
    burstmeanz
  )
]


acc_next[
  ,
  temporal_direction :=
    "next"
]


#------------------------------------------------------------------------------
# 3. Combine previous and next candidates
#------------------------------------------------------------------------------

gps_acc_candidates <- data.table::rbindlist(
  list(
    acc_previous,
    acc_next
  ),
  
  use.names =
    TRUE,
  
  fill =
    TRUE
)


# Calculate the absolute temporal distance.

gps_acc_candidates[
  ,
  abs_time_diff_min :=
    abs(
      as.numeric(
        difftime(
          gps_timestamp,
          acc_timestamp,
          units = "mins"
        )
      )
    )
]


# Missing ACC candidates must be ranked after valid candidates.

gps_acc_candidates[
  ,
  temporal_distance_order :=
    data.table::fifelse(
      is.na(
        abs_time_diff_min
      ),
      Inf,
      abs_time_diff_min
    )
]


# Use RF probability as a tie-breaking criterion.

gps_acc_candidates[
  ,
  probability_order :=
    data.table::fifelse(
      is.na(
        pro_rf8fitted
      ),
      -Inf,
      pro_rf8fitted
    )
]


# Deterministic final tie-breaking rule:
# when distance and probability are equal, retain the previous ACC record.

gps_acc_candidates[
  ,
  direction_order :=
    data.table::fifelse(
      temporal_direction ==
        "previous",
      1L,
      2L
    )
]


data.table::setorder(
  gps_acc_candidates,
  
  gps_row_id,
  
  temporal_distance_order,
  
  -probability_order,
  
  direction_order
)


#------------------------------------------------------------------------------
# 4. Retain exactly one nearest ACC record for each GPS point
#------------------------------------------------------------------------------

gps_nearest_acc <- gps_acc_candidates[
  ,
  .SD[1L],
  by =
    gps_row_id
]


# The result must now contain exactly one row per GPS point.

stopifnot(
  nrow(
    gps_nearest_acc
  ) ==
    nrow(
      gps_points
    )
)

stopifnot(
  anyDuplicated(
    gps_nearest_acc$gps_row_id
  ) == 0L
)


#------------------------------------------------------------------------------
# Step 2.5: Apply the temporal assignment rule ----
#------------------------------------------------------------------------------
# Rule:
# For each GPS point:
# - retain the closest ACC behaviour from the same individual
# - keep the assignment only if the temporal distance is <= 60 min
# - otherwise leave behaviour as NA
#------------------------------------------------------------------------------


max_assignment_gap_min <- 60


# Identify valid behaviour assignments

gps_nearest_acc[
  ,
  behavior_assigned :=
    !is.na(acc_timestamp) &
    !is.na(rf8fitted) &
    !is.na(abs_time_diff_min) &
    abs_time_diff_min <= max_assignment_gap_min
]


# Record how the behaviour was assigned

gps_nearest_acc[
  ,
  behavior_assignment_method :=
    data.table::fifelse(
      
      behavior_assigned,
      
      paste0(
        "nearest_ACC_within_",
        max_assignment_gap_min,
        "_min"
      ),
      
      "no_ACC_within_maximum_gap"
    )
]


# Remove behaviour labels that are outside the accepted temporal window.
# Keep the temporal information for diagnostics.

gps_nearest_acc[
  behavior_assigned == FALSE,
  `:=`(
    
    rf8fitted =
      NA_character_,
    
    pro_rf8fitted =
      NA_real_
  )
]


# Order final GPS-level dataset

data.table::setorder(
  gps_nearest_acc,
  individual_local_identifier,
  gps_timestamp
)


#------------------------------------------------------------------------------
# Step 2.6: Diagnostic of behaviour assignment ----
#------------------------------------------------------------------------------

gps_behavior_assignment_overall <-
  
  gps_nearest_acc[
    ,
    .(
      
      n_gps =
        .N,
      
      n_gps_with_behavior =
        sum(
          !is.na(rf8fitted)
        ),
      
      n_gps_without_behavior =
        sum(
          is.na(rf8fitted)
        ),
      
      prop_gps_with_behavior =
        mean(
          !is.na(rf8fitted)
        ),
      
      median_abs_time_diff_min =
        median(
          abs_time_diff_min,
          na.rm = TRUE
        ),
      
      maximum_abs_time_diff_min =
        max(
          abs_time_diff_min,
          na.rm = TRUE
        )
    )
  ]


print(
  gps_behavior_assignment_overall
)



#------------------------------------------------------------------------------
# Step 2.7: Create the final GPS dataset with behaviour + geometry ----
#------------------------------------------------------------------------------

# Keep only the variables that must be transferred back to GPS.

gps_behavior_assignment <-
  
  gps_nearest_acc[
    ,
    .(
      
      gps_row_id,
      
      rf8fitted,
      
      pro_rf8fitted,
      
      acc_row_id,
      
      acc_timestamp,
      
      abs_time_diff_min,
      
      behavior_assigned,
      
      behavior_assignment_method
    )
  ]


# Remove possible old behaviour columns before joining.
# This avoids creating rf8fitted.x / rf8fitted.y.

gps_with_behavior <-
  
  gps_sf %>%
  
  dplyr::select(
    
    -dplyr::any_of(
      c(
        "rf8fitted",
        "pro_rf8fitted",
        "behavior_assigned",
        "acc_timestamp",
        "abs_time_diff_min",
        "behavior_assignment_method"
      )
    )
  ) %>%
  
  dplyr::left_join(
    
    gps_behavior_assignment,
    
    by =
      "gps_row_id"
  )



# Verify that the spatial dataset contains the expected number
# of assigned behaviours.

print(
  gps_with_behavior %>%
    sf::st_drop_geometry() %>%
    dplyr::summarise(
      
      n_gps =
        dplyr::n(),
      
      n_with_behavior =
        sum(
          !is.na(rf8fitted)
        ),
      
      n_without_behavior =
        sum(
          is.na(rf8fitted)
        ),
      
      proportion_with_behavior =
        mean(
          !is.na(rf8fitted)
        )
    )
)



#------------------------------------------------------------------------------
# Step 2.8: Create one spatial object per individual ----
#------------------------------------------------------------------------------

row_indices_by_individual <-
  
  split(
    
    seq_len(
      nrow(
        gps_with_behavior
      )
    ),
    
    gps_with_behavior$
      individual_local_identifier
  )


gps_behavior_by_individual <-
  
  lapply(
    
    row_indices_by_individual,
    
    function(row_indices){
      
      gps_with_behavior[
        row_indices,
      ]
    }
  )


print(
  length(
    gps_behavior_by_individual
  )
)



#------------------------------------------------------------------------------
# Step 2.9: Prepare dataset for behaviour reclassification ----
#------------------------------------------------------------------------------

gps_beh_15w <-
  
  gps_with_behavior %>%
  
  sf::st_drop_geometry() %>%
  
  dplyr::rename(
    
    individual.local.identifier =
      individual_local_identifier
  )


data.table::setDT(
  gps_beh_15w
)


# Remove GPS points without assigned behaviour.
# They cannot enter a Markov chain because their state is unknown.

gps_beh_15w <-
  
  gps_beh_15w[
    !is.na(rf8fitted)
  ]


# Harmonise data types for Step 3

gps_beh_15w[
  ,
  individual.local.identifier :=
    as.character(
      individual.local.identifier
    )
]


gps_beh_15w[
  ,
  gps_timestamp :=
    as.POSIXct(
      gps_timestamp,
      tz = "UTC"
    )
]


gps_beh_15w[
  ,
  rf8fitted :=
    as.character(
      rf8fitted
    )
]


data.table::setorder(
  gps_beh_15w,
  individual.local.identifier,
  gps_timestamp
)



# Diagnostic before entering Step 3

print(
  gps_beh_15w[
    ,
    .(
      n_rows =
        .N,
      
      n_individuals =
        uniqueN(
          individual.local.identifier
        ),
      
      n_behaviour_NA =
        sum(
          is.na(rf8fitted)
        )
    )
  ]
)
    



#' -----------------------------------------------------------------------------
#' Step 3 : reclassify behaviors
#' 
#' **Philosophy**: we group certain behavioural categories together in order 
#' to ensure significantly large classes of behaviors to calculate transition 
#' probability later. We focus on three groups : feeding, resting, flight. 
#' 
#' **Steps**:
#' (1) create a new column called behaviors reclassified 
#' (2) reclassify the behaviors : 
#'       - flight correspond to the category 'undulating' and 'active'
#'       - resting correspond to the category 'bodycare', 'standing' and 'passive'
#'       - feeding correspond to 'feeding'
#' (3) for 'walking' : (a) identify consecutive walking events within individual
#' temporal sequences; (b) find the nearest non-walking terrestrial states before
#' and after the block;(c) if these behavior are the same, assign the whole 
#' walking block to that state and if not, use the majority among the 2 nearest
#' terrestrial events before and the 2 nearest terrestrial events after; (d) if
#' equality, use the temporally closest terrestrial event; (e) if unresolved, 
#' keep NA and exclude from the Markov model later.


data.table::setDT(gps_beh_15w)

required_cols <- c(
  "individual.local.identifier",
  "gps_timestamp",
  "rf8fitted"
)

missing_cols <- setdiff(required_cols, names(gps_beh_15w))

if (length(missing_cols) > 0) {
  stop(
    "Missing required columns in gps_beh_15w: ",
    paste(missing_cols, collapse = ", ")
  )
}

gps_beh_15w[, individual.local.identifier := as.character(individual.local.identifier)]
gps_beh_15w[, gps_timestamp := as.POSIXct(gps_timestamp, tz = "UTC")]
gps_beh_15w[, rf8fitted := as.character(rf8fitted)]

data.table::setorder(
  gps_beh_15w,
  individual.local.identifier,
  gps_timestamp
)

# Unique row identifier after ordering
gps_beh_15w[, row_index := .I]

# Maximum gap allowed when using neighboring behavior to reclassify walking.
# This prevents classifying a walking event using a distant context.
if (!exists("max_assignment_gap_min")) {
  max_assignment_gap_min <- 60
}

max_context_gap_min <- 2 * max_assignment_gap_min

# Number of terrestrial non-walking events used on each side of a walking block
# when immediate flanks disagree.
k_context <- 2


#' Step 3.1: direct reclassification of non-walking behaviors ----
gps_beh_15w[,
  behavior_base := data.table::fcase(
    rf8fitted %in% c("Active", "Undulating", "Passive"),
    "flight",
    
    rf8fitted %in% c("Bodycare", "Standing"),
    "resting",
    
    rf8fitted == "Feeding",
    "feeding",
    
    rf8fitted == "Walking",
    "walking",
    
    default = NA_character_
  )
]

gps_beh_15w[
  ,
  behavior_reclassified := behavior_base
]

gps_beh_15w[
  ,
  behavior_reclassification_rule := data.table::fcase(
    behavior_base %in% c("flight", "resting", "feeding"),
    "direct_reclassification",
    
    behavior_base == "walking",
    "walking_to_be_reclassified",
    
    is.na(behavior_base),
    "no_behavior_assigned",
    
    default = NA_character_
  )
]

# Walking is not a final Markov state in this three-state version.
# It will be reassigned block by block below.
gps_beh_15w[
  behavior_base == "walking",
  behavior_reclassified := NA_character_
]


#' Step 3.2: define temporal sequences within individuals ----
# A new temporal sequence starts when there is a large gap between consecutive
# GPS points. This avoids using an old behavior to classify a later walking event.

gps_beh_15w[
  ,
  time_gap_from_previous_min := as.numeric(
    difftime(
      gps_timestamp,
      data.table::shift(gps_timestamp),
      units = "mins"
    )
  ),
  by = individual.local.identifier
]

gps_beh_15w[
  ,
  new_context_sequence := is.na(time_gap_from_previous_min) |
    time_gap_from_previous_min > max_context_gap_min,
  by = individual.local.identifier
]

gps_beh_15w[
  ,
  context_sequence_id := cumsum(new_context_sequence),
  by = individual.local.identifier
]

gps_beh_15w[
  ,
  row_in_context_sequence := seq_len(.N),
  by = .(individual.local.identifier, context_sequence_id)
]


#' Step 3.3: identify consecutive walking blocks ----
gps_beh_15w[
  ,
  is_walking := !is.na(behavior_base) & behavior_base == "walking"
]

gps_beh_15w[
  ,
  walking_block_id := data.table::fifelse(
    is_walking,
    cumsum(
      is_walking &
        !data.table::shift(is_walking, fill = FALSE)
    ),
    NA_integer_
  ),
  by = .(individual.local.identifier, context_sequence_id)
]

walking_blocks <- gps_beh_15w[
  is_walking == TRUE,
  .(
    block_start_row = min(row_index),
    block_end_row = max(row_index),
    block_start_time = min(gps_timestamp),
    block_end_time = max(gps_timestamp),
    first_row_in_context = min(row_in_context_sequence),
    last_row_in_context = max(row_in_context_sequence),
    n_walking_in_block = .N
  ),
  by = .(
    individual.local.identifier,
    context_sequence_id,
    walking_block_id
  )
]


#' Step 3.4: reclassify walking blocks from local terrestrial context ----

terrestrial_states <- c("resting", "feeding")

for (b in seq_len(nrow(walking_blocks))) {
  
  this_ind <- walking_blocks$individual.local.identifier[b]
  this_context <- walking_blocks$context_sequence_id[b]
  this_block <- walking_blocks$walking_block_id[b]
  
  this_first_pos <- walking_blocks$first_row_in_context[b]
  this_last_pos <- walking_blocks$last_row_in_context[b]
  this_start_time <- walking_blocks$block_start_time[b]
  this_end_time <- walking_blocks$block_end_time[b]
  
  block_rows <- gps_beh_15w[
    individual.local.identifier == this_ind &
      context_sequence_id == this_context &
      walking_block_id == this_block,
    row_index
  ]
  
  context_dt <- gps_beh_15w[
    individual.local.identifier == this_ind &
      context_sequence_id == this_context
  ]
  
  before_candidates <- context_dt[
    row_in_context_sequence < this_first_pos &
      behavior_base %in% terrestrial_states
  ][
    order(-row_in_context_sequence)
  ]
  
  after_candidates <- context_dt[
    row_in_context_sequence > this_last_pos &
      behavior_base %in% terrestrial_states
  ][
    order(row_in_context_sequence)
  ]
  
  before_candidates <- before_candidates[
    seq_len(min(nrow(before_candidates), k_context))
  ]
  
  after_candidates <- after_candidates[
    seq_len(min(nrow(after_candidates), k_context))
  ]
  
  left_state <- if (nrow(before_candidates) >= 1) {
    before_candidates$behavior_base[1]
  } else {
    NA_character_
  }
  
  right_state <- if (nrow(after_candidates) >= 1) {
    after_candidates$behavior_base[1]
  } else {
    NA_character_
  }
  
  assigned_state <- NA_character_
  assigned_rule <- "walking_unresolved_no_terrestrial_context"
  
  # Rule 1: nearest terrestrial flanks agree
  if (
    !is.na(left_state) &&
    !is.na(right_state) &&
    left_state == right_state
  ) {
    
    assigned_state <- left_state
    assigned_rule <- "walking_flanked_by_same_terrestrial_state"
    
  } else {
    
    local_candidates <- data.table::rbindlist(
      list(before_candidates, after_candidates),
      use.names = TRUE,
      fill = TRUE
    )
    
    if (nrow(local_candidates) > 0) {
      
      n_resting <- sum(local_candidates$behavior_base == "resting")
      n_feeding <- sum(local_candidates$behavior_base == "feeding")
      
      # Rule 2: majority among local terrestrial context
      if (n_resting > n_feeding) {
        
        assigned_state <- "resting"
        assigned_rule <- "walking_reclassified_by_local_majority"
        
      } else if (n_feeding > n_resting) {
        
        assigned_state <- "feeding"
        assigned_rule <- "walking_reclassified_by_local_majority"
        
      } else {
        
        # Rule 3: tie resolved by nearest terrestrial event in time
        local_candidates[
          row_in_context_sequence < this_first_pos,
          distance_to_walking_block_min := as.numeric(
            difftime(this_start_time, gps_timestamp, units = "mins")
          )
        ]
        
        local_candidates[
          row_in_context_sequence > this_last_pos,
          distance_to_walking_block_min := as.numeric(
            difftime(gps_timestamp, this_end_time, units = "mins")
          )
        ]
        
        data.table::setorder(
          local_candidates,
          distance_to_walking_block_min
        )
        
        min_distance <- local_candidates$distance_to_walking_block_min[1]
        
        closest_candidates <- local_candidates[
          distance_to_walking_block_min == min_distance
        ]
        
        closest_states <- unique(closest_candidates$behavior_base)
        
        if (length(closest_states) == 1) {
          
          assigned_state <- closest_states[1]
          assigned_rule <- "walking_tie_resolved_by_nearest_time"
          
        } else {
          
          assigned_state <- NA_character_
          assigned_rule <- "walking_unresolved_equal_context"
        }
      }
    }
  }
  
  gps_beh_15w[
    row_index %in% block_rows,
    `:=`(
      behavior_reclassified = assigned_state,
      behavior_reclassification_rule = assigned_rule,
      walking_block_size = length(block_rows)
    )
  ]
}


#' Step 3.5: final formatting ----

state_levels_3 <- c(
  "flight",
  "resting",
  "feeding"
)

gps_beh_15w[
  ,
  behavior_reclassified := factor(
    behavior_reclassified,
    levels = state_levels_3
  )
]


#' Step 3.6: diagnostics ----
behavior_counts_raw <- gps_beh_15w[
  ,
  .N,
  by = rf8fitted
][
  order(-N)
]

behavior_counts_base <- gps_beh_15w[
  ,
  .N,
  by = behavior_base
][
  order(-N)
]

behavior_counts_reclassified <- gps_beh_15w[
  ,
  .N,
  by = behavior_reclassified
][
  order(-N)
]

print(behavior_counts_raw)
print(behavior_counts_base)
print(behavior_counts_reclassified)

#' Step 3.7: object for later thinning and Markov modelling ----
# Rows with NA behavior_reclassified are kept in gps_beh_15w for diagnostics,
# but should be removed before constructing Markov transitions.
gps_beh_15w_markov_ready <- gps_beh_15w[
  !is.na(behavior_reclassified)]

saveRDS(
  gps_beh_15w_markov_ready,
  file.path(
    output_dir,
    "acc-joint-gps-15weeks-raw.rds"
  ),
  compress = "gzip")




#' -----------------------------------------------------------------------------
#' Step 4 : creation of two temporally regular dataset
#' 
#' **Steps:**
#' (1) inspect raw GPS sampling intervals;
#' (2) create one high-resolution dataset at 20 min;
#' (3) create one intermediate-resolution dataset at 60 min;
#' (4) split tracks into bursts

gps_beh_15w_markov_ready <- readRDS("/Users/louisefaure/Library/CloudStorage/OneDrive-Personnel/THESE/CHAPITRE 2/git/chapter-2/HMM/HMM on ACC-classified behaviors/donnees intermediaire (2)/acc-joint-gps-15weeks-raw.rds")


library(data.table)
library(dplyr)
library(sf)


# 4.1 Prepare the non-spatial input table ----
step4_input <- data.table::as.data.table(
  data.table::copy(
    gps_beh_15w_markov_ready
  )
)

step4_input <- step4_input[
  behavior_assigned %in% TRUE &
    !is.na(behavior_reclassified) &
    !is.na(individual.local.identifier) &
    !is.na(gps_timestamp) &
    !is.na(lon) &
    !is.na(lat)
]


step4_input[
  ,
  individual.id :=
    as.character(
      individual.local.identifier
    )
]

step4_input[
  ,
  timestamp :=
    as.POSIXct(
      gps_timestamp,
      tz = "UTC"
    )
]

step4_input[
  ,
  lon :=
    as.numeric(
      lon
    )
]

step4_input[
  ,
  lat :=
    as.numeric(
      lat
    )
]


# Remove exact duplicate individual-timestamp combinations.
data.table::setorder(
  step4_input,
  individual.id,
  timestamp
)

step4_input <- unique(
  step4_input,
  by = c(
    "individual.id",
    "timestamp"
  )
)


# Check the number of rows before thinning.
cat(
  "Rows before 60-min thinning:",
  nrow(step4_input),
  "\n"
)


# 4.2 Select one GPS point per 60-min interval ----
resolution_min_60 <- 60
tolerance_min_60 <- 10
min_points_per_burst <- 3L


# Unix time expressed in seconds.

step4_input[
  ,
  timestamp_numeric :=
    as.numeric(
      timestamp
    )
]


# Create fixed one-hour calendar bins in UTC.
#
# floor(timestamp / 3600) assigns all points between, for example,
# 10:00:00 and 10:59:59 to the same temporal interval.

step4_input[
  ,
  interval_60_id :=
    floor(
      timestamp_numeric /
        (resolution_min_60 * 60)
    )
]


# Retain the first GPS observation in every one-hour interval
# for each individual.

thin_60_dt <- step4_input[
  order(
    individual.id,
    timestamp
  ),
  .SD[1L],
  by = .(
    individual.id,
    interval_60_id
  )
]


data.table::setorder(
  thin_60_dt,
  individual.id,
  timestamp
)


thin_60_dt[
  ,
  `:=`(
    individual.local.identifier =
      individual.id,
    
    resolution_min =
      resolution_min_60,
    
    resolution_class =
      "60min_intermediate_resolution"
  )
]


cat(
  "Rows after 60-min thinning:",
  nrow(thin_60_dt),
  "\n"
)

cat(
  "Proportion retained:",
  nrow(thin_60_dt) /
    nrow(step4_input),
  "\n"
)


#'------------------------------------------------------------------------------
# 4.3 Split the thinned data into regular temporal bursts ----
#'------------------------------------------------------------------------------

thin_60_dt[
  ,
  dt_prev_min :=
    as.numeric(
      difftime(
        timestamp,
        data.table::shift(
          timestamp
        ),
        units = "mins"
      )
    ),
  by =
    individual.id
]


# A new burst begins when:
# - the row is the first observation of the individual;
# - or the interval differs from 60 min by more than 10 min.

thin_60_dt[
  ,
  new_burst :=
    is.na(
      dt_prev_min
    ) |
    abs(
      dt_prev_min -
        resolution_min_60
    ) >
    tolerance_min_60
]


thin_60_dt[
  ,
  burst_n :=
    cumsum(
      new_burst
    ),
  by =
    individual.id
]


thin_60_dt[
  ,
  burst_id :=
    paste(
      individual.id,
      "60min",
      sprintf(
        "%04d",
        burst_n
      ),
      sep = "_"
    )
]


# Count observations per burst.

thin_60_dt[
  ,
  n_points_burst :=
    .N,
  by = .(
    individual.id,
    burst_id
  )
]


# Retain bursts containing at least three points.

regular_60_dt <- thin_60_dt[
  n_points_burst >=
    min_points_per_burst
]


data.table::setorder(
  regular_60_dt,
  individual.id,
  burst_id,
  timestamp
)


#'------------------------------------------------------------------------------
# 4.4 Add within-burst temporal diagnostics ----
#'------------------------------------------------------------------------------

regular_60_dt[
  ,
  row_in_burst :=
    seq_len(
      .N
    ),
  by = .(
    individual.id,
    burst_id
  )
]


regular_60_dt[
  ,
  dt_prev_burst_min :=
    as.numeric(
      difftime(
        timestamp,
        data.table::shift(
          timestamp
        ),
        units = "mins"
      )
    ),
  by = .(
    individual.id,
    burst_id
  )
]


regular_60_dt[
  ,
  dt_next_burst_min :=
    as.numeric(
      difftime(
        data.table::shift(
          timestamp,
          type = "lead"
        ),
        timestamp,
        units = "mins"
      )
    ),
  by = .(
    individual.id,
    burst_id
  )
]


regular_60_dt[
  ,
  has_next_regular :=
    !is.na(
      dt_next_burst_min
    ) &
    abs(
      dt_next_burst_min -
        resolution_min_60
    ) <=
    tolerance_min_60
]


regular_60_dt[
  ,
  lag_deviation_next_min :=
    dt_next_burst_min -
    resolution_min_60
]


#'------------------------------------------------------------------------------
# 4.5 Control the resulting dataset ----
#'------------------------------------------------------------------------------

timing_60 <- regular_60_dt[
  ,
  .(
    resolution_min =
      data.table::first(
        resolution_min
      ),
    
    n_individuals =
      data.table::uniqueN(
        individual.local.identifier
      ),
    
    n_points_retained =
      .N,
    
    n_bursts =
      data.table::uniqueN(
        burst_id
      ),
    
    n_valid_transitions =
      sum(
        has_next_regular,
        na.rm = TRUE
      ),
    
    n_terminal_points =
      sum(
        !has_next_regular |
          is.na(
            has_next_regular
          )
      )
  )
]


timing_60[
  ,
  expected_valid_transitions :=
    n_points_retained -
    n_bursts
]


print(
  timing_60
)


# Individual-level diagnostic.

summary_by_individual_60 <- regular_60_dt[
  ,
  .(
    n_points_retained =
      .N,
    
    n_bursts =
      data.table::uniqueN(
        burst_id
      ),
    
    n_valid_transitions =
      sum(
        has_next_regular,
        na.rm = TRUE
      ),
    
    median_points_per_burst =
      as.numeric(
        stats::median(
          n_points_burst,
          na.rm = TRUE
        )
      ),
    
    max_points_per_burst =
      max(
        n_points_burst,
        na.rm = TRUE
      )
  ),
  by = .(
    resolution_class,
    individual.local.identifier
  )
]


low_data_individuals_60 <-
  summary_by_individual_60[
    n_valid_transitions < 200
  ][
    order(
      n_valid_transitions
    )
  ]


print(
  low_data_individuals_60,
  nrows = Inf
)


#'------------------------------------------------------------------------------
# 4.6 Create the spatial dataset only after thinning ----
#'------------------------------------------------------------------------------

# Remove temporary variables that are not required downstream.

regular_60_tbl <- data.table::copy(
  regular_60_dt
)

regular_60_tbl[
  ,
  c(
    "timestamp_numeric",
    "interval_60_id",
    "new_burst"
  ) :=
    NULL
]


# Create points initially in geographic coordinates.

regular_60_sf <- sf::st_as_sf(
  as.data.frame(
    regular_60_tbl
  ),
  coords = c(
    "lon",
    "lat"
  ),
  crs = 4326,
  remove = FALSE
)


# Project only the retained points to EPSG:3035.

regular_60_sf <- sf::st_transform(
  regular_60_sf,
  3035
)


# Add projected coordinates required by later analyses.

xy_3035 <- sf::st_coordinates(
  regular_60_sf
)

regular_60_sf$x_3035 <-
  as.numeric(
    xy_3035[, 1]
  )

regular_60_sf$y_3035 <-
  as.numeric(
    xy_3035[, 2]
  )


# Final controls.

stopifnot(
  all(
    is.finite(
      regular_60_sf$x_3035
    )
  ),
  
  all(
    is.finite(
      regular_60_sf$y_3035
    )
  ),
  
  sf::st_crs(
    regular_60_sf
  )$epsg ==
    3035
)


print(
  regular_60_sf
)


#'------------------------------------------------------------------------------
# 4.7 Save the 60-min dataset ----
#'------------------------------------------------------------------------------

saveRDS(
  regular_60_sf,
  file.path(
    output_dir,
    "GE_60_min_thinned_behavior_assigned2.rds"
  ),
  compress = "gzip"
)




#'------------------------------------------------------------------------------
#' Step 5 : extract human footprint index around GPS locations
#' 
#' **Extract:**
#' - HFI value of the pixel containing the GPS point
#' - mean, maximum, and 75th percentile HFI within 500 m
#' - mean, maximum, and 75th percentile HFI within 1000 m



regular_20_sf <- readRDS("/Users/louisefaure/Library/CloudStorage/OneDrive-Personnel/THESE/CHAPITRE 2/git/chapter-2/HMM/HMM on ACC-classified behaviors/donnees intermediaire (2)/GE_20_min_thinned_behavior_assigned.rds")
terra::terraOptions(threads = 5)

# Function for HFI extraction
terra::terraOptions(threads = 5)

# Function for HFI extraction
extract_hfi <- function(df_sf,
                        raster = human_footprint,
                        buffers_m = c(500, 1000)) {
  
  # Keep attribute table without geometry
  df <- as.data.frame(sf::st_drop_geometry(df_sf))
  
  # Create terra points directly from sf geometry
  pts <- terra::vect(df_sf)
  
  # Project points to raster CRS if needed
  pts_r <- terra::project(pts, terra::crs(raster))
  
  # Extract HFI value at GPS point
  df$hfi_point <- terra::extract(
    raster,
    pts_r,
    ID = FALSE
  )[[1]]
  
  # Extract HFI statistics in buffers
  for (r in buffers_m) {
    
    # Create buffers in the metric CRS of the GPS data, EPSG:3035
    buf <- terra::buffer(pts, width = r)
    
    # Project buffers to raster CRS
    buf_r <- terra::project(buf, terra::crs(raster))
    
    # Mean HFI inside buffer
    df[[paste0("hfi_mean_", r, "m")]] <- terra::extract(
      raster,
      buf_r,
      fun = mean,
      na.rm = TRUE,
      ID = FALSE
    )[[1]]
    
    # Third quartile of HFI inside buffer
    df[[paste0("hfi_q75_", r, "m")]] <- terra::extract(
      raster,
      buf_r,
      fun = function(x, ...) {
        stats::quantile(
          x,
          probs = 0.90,
          na.rm = TRUE,
          names = FALSE
        )
      },
      ID = FALSE
    )[[1]]
    
    # Third quartile of HFI inside buffer
    df[[paste0("hfi_q75_", r, "m")]] <- terra::extract(
      raster,
      buf_r,
      fun = function(x, ...) {
        stats::quantile(
          x,
          probs = 0.75,
          na.rm = TRUE,
          names = FALSE
        )
      },
      ID = FALSE
    )[[1]]
  }
  
  df
}

regular_20_hfi <- extract_hfi(regular_20_sf)
regular_60_hfi <- extract_hfi(regular_60_sf)

# Controls
hfi_columns <- c(
  "hfi_point",
  "hfi_mean_500m",
  "hfi_q75_500m",
  "hfi_mean_1000m",
  "hfi_q75_1000m", 
  "hfi_q90_1000m",
  "hfi_q90_500m"
)

regular_20_hfi[, hfi_columns] |>
  summary()

regular_60_hfi[, hfi_columns] |>
  summary()

# clean 


regular_60_hfi_clean <- regular_60_hfi


# Harmonise individual identifier.

if (
  !"individual.local.identifier" %in%
  names(regular_60_hfi_clean) &&
  "individual_local_identifier" %in%
  names(regular_60_hfi_clean)
) {
  regular_60_hfi_clean <-
    regular_60_hfi_clean %>%
    dplyr::rename(
      individual.local.identifier =
        individual_local_identifier
    )
}


# Harmonise ground-speed name.

if (
  !"ground_speed" %in%
  names(regular_60_hfi_clean) &&
  "ground.speed" %in%
  names(regular_60_hfi_clean)
) {
  regular_60_hfi_clean <-
    regular_60_hfi_clean %>%
    dplyr::rename(
      ground_speed =
        ground.speed
    )
}


# Harmonise height name.

if (
  !"height_above_ellipsoid" %in%
  names(regular_60_hfi_clean) &&
  "height.above.ellipsoid" %in%
  names(regular_60_hfi_clean)
) {
  regular_60_hfi_clean <-
    regular_60_hfi_clean %>%
    dplyr::rename(
      height_above_ellipsoid =
        height.above.ellipsoid
    )
}


# Harmonise timestamp.

if (
  !"timestamp" %in%
  names(regular_60_hfi_clean) &&
  "gps_timestamp" %in%
  names(regular_60_hfi_clean)
) {
  regular_60_hfi_clean <-
    regular_60_hfi_clean %>%
    dplyr::rename(
      timestamp =
        gps_timestamp
    )
}

columns_to_keep_60 <- c(
  "sensor_type_id",
  "individual.local.identifier",
  "eobs_horizontal_accuracy_estimate",
  "eobs_speed_accuracy_estimate",
  "eobs_type_of_fix",
  "gps_dop",
  "gps_satellite_count",
  "ground_speed",
  "height_above_ellipsoid",
  "timestamp",
  "lon",
  "lat",
  "dist.traveled",
  "rf8fitted",
  "behavior_assignment_method",
  "behavior_base",
  "behavior_reclassified",
  "behavior_reclassification_rule",
  "burst_n",
  "row_in_burst",
  "hfi_point",
  "hfi_mean_500m",
  "hfi_q75_500m",
  "hfi_mean_1000m",
  "hfi_q75_1000m"
)

missing_final_columns_60 <- setdiff(
  columns_to_keep_60,
  names(regular_60_hfi_clean)
)

print(
  missing_final_columns_60
)

if (
  length(
    missing_final_columns_60
  ) > 0L
) {
  stop(
    paste0(
      "Missing required output columns: ",
      paste(
        missing_final_columns_60,
        collapse = ", "
      )
    )
  )
}


regular_60_hfi_clean <-
  regular_60_hfi_clean %>%
  dplyr::select(
    dplyr::all_of(
      columns_to_keep_60
    )
  )


hfi_columns <- c(
  "hfi_point",
  "hfi_mean_500m",
  "hfi_q75_500m",
  "hfi_mean_1000m",
  "hfi_q75_1000m"
)


summary(
  regular_60_hfi_clean[
    ,
    hfi_columns,
    drop = FALSE
  ]
)



# Save outputs
saveRDS(
  regular_20_hfi,
  file.path(output_dir, "GE_20_min_thinned_behavior_assigned_hfi.rds"),
  compress = "gzip"
)

saveRDS(
  regular_60_hfi_clean,
  file.path(output_dir, "GE_60_min_thinned_behavior_assigned_hfi.rds"),
  compress = "gzip"
)
