#' ############################################################################# 
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
#' #############################################################################



library(terra)           # for raster
library(data.table)      # to read the acc-classified data



#' ### Step 0 : load the data ----
output_dir <- "C:/Users/lfaure7/OneDrive/THESE/CHAPITRE 2/git/chapter-2/HMM/HMM on ACC-classified behaviors/donnes filtree intermediaire"

# Géoïde, digital elevation model, human footprint index
geo <- terra::rast("C:/Users/lfaure7/OneDrive/MEMOIRE M2/eagle_projet/data/geoide/us_nga_egm96_15.tif")
dem <- terra::rast("C:/Users/lfaure7/OneDrive/MEMOIRE M2/eagle_projet/data/pretraitements/Region-Alpes-Dem/region-alpes-dem.tif")
human_footprint <- terra::rast("C:/Users/lfaure7/OneDrive/THESE/COUCHES QGIS/COUCHES QGIS/settlements/Overture/human_footprint_index_building_pop_builtprop_100m.tif")

# emigration dates and golden eagle data
emig_dates <- readRDS("C:/Users/lfaure7/Documents/git/chapter-2/preparation donnee aigle ssf/donnees/emigration dates/emigration_dates_20250417.rds")
classified_path <- "C:/Users/lfaure7/OneDrive/THESE/CHAPITRE 2/donnees/acc-data/Louise_PhD/classified_acc_data/2024_01_24_alldata_allbirds_merged_rf_raw.csv"


#'### Step 1. Retain the first fifteen weeks of the dispersal phase ----
#'**Steps:** 
#' (1) define column names used in the ACC and emigration files,
#' (2) load the full ACC-classified dataset,
#' (3) prepare emigration dates,
#' (4) retain only the first 15 weeks after dispersal.

# 1.0 Define colums names in the ACC classification dataset ----
acc_id_col       <- "individualID"     # individuals
acc_time_col     <- "timestamp"        # time
acc_behavior_col <- "rf8fitted"        # behavior
acc_prob_col     <- "pro_rf8fitted"    # classification probability

# Load the full classified ACC file
acc_classified <- data.table::fread(
  classified_path,
  showProgress = TRUE,
  tz = "UTC"
)

# Ensure data.table format
data.table::setDT(acc_classified)


# 1.1 Prepare emigration dates ----
# Convert emig_dates to data.table
emig_dt <- data.table::as.data.table(emig_dates)

# Column names in emig_dates
emig_id_col        <- "individual.local.identifier"
did_disperse_col   <- "did_disperse"
dispersal_date_col <- "dispersal_date"

#'==============================================================================
#' 1.3 Check and convert dispersal date format if needed ----
#'==============================================================================

# CONTROL: check dispersal_date temporal format
# Expected result:
# - TRUE for dispersal_date_is_posix if dispersal_date is already POSIXct.
# - If TRUE, skip conversion and copy the column into dispersal_start_utc.
# - If FALSE, conversion is necessary.
dispersal_date_is_posix <- inherits(emig_dt[[dispersal_date_col]], "POSIXct")
dispersal_date_is_date  <- inherits(emig_dt[[dispersal_date_col]], "Date")

cat("dispersal_date already POSIXct:", dispersal_date_is_posix, "\n")
cat("dispersal_date is Date:", dispersal_date_is_date, "\n")
cat("dispersal_date class:", paste(class(emig_dt[[dispersal_date_col]]), collapse = " "), "\n")
print(head(emig_dt[[dispersal_date_col]]))

# If dispersal_date is already POSIXct, no conversion is needed.
# If it is Date, convert it to POSIXct at midnight UTC.
# Otherwise, use parse_datetime_utc().
if (dispersal_date_is_posix) {
  
  # No temporal conversion needed.
  emig_dt[, dispersal_start_utc := get(dispersal_date_col)]
  
} else if (dispersal_date_is_date) {
  
  # Date has no hour/min/sec; convert to midnight UTC.
  emig_dt[, dispersal_start_utc := as.POSIXct(as.Date(get(dispersal_date_col)), tz = "UTC")]
  
} else {
  
  # Conversion needed.
  emig_dt[, dispersal_start_utc := parse_datetime_utc(get(dispersal_date_col))]
}

# CONTROL: check dispersal_start_utc after conversion/copy
# Expected result:
# - no or very few missing values.
# - dates should correspond to biologically plausible dispersal dates.
cat("Missing dispersal_start_utc:", sum(is.na(emig_dt$dispersal_start_utc)), "\n")
print(head(emig_dt[, .(individual.local.identifier, dispersal_date, dispersal_start_utc)]))




# 1.2 Keep only individuals with a confirmed dispersal date ----
# Convert did_disperse to logical without modifying the original column.
# This is robust if did_disperse is TRUE/FALSE, 1/0, or character.
did_raw <- emig_dt[[did_disperse_col]]

if (is.logical(did_raw)) {
  did_logical <- did_raw
} else if (is.numeric(did_raw)) {
  did_logical <- did_raw == 1
} else {
  did_logical <- tolower(trimws(as.character(did_raw))) %in%
    c("true", "t", "yes", "y", "1")
}

# Keep only individuals that dispersed and have a dispersal date
emig_filtered <- emig_dt[
  did_logical == TRUE &
    !is.na(get(dispersal_date_col))
]



# 1.3 Convert dispersal date to POSIXct if needed ----
# If parse_datetime_utc() was defined earlier, keep using it.
# If dispersal_date is already POSIXct, this will keep it in the correct format.
emig_filtered[
  ,
  dispersal_start_utc := parse_datetime_utc(get(dispersal_date_col))
]

# Define end of the first 15 weeks
# 15 weeks = 105 days
emig_filtered[
  ,
  dispersal_end_utc := dispersal_start_utc + 15 * 7 * 24 * 60 * 60
]

# Keep one row per individual if duplicate rows have the same dispersal date
emig_filtered <- unique(emig_filtered, by = emig_id_col)


#'==============================================================================
#' 1.4 Retain only individuals present in both files
#'==============================================================================

common_ids <- intersect(
  unique(acc_classified[[acc_id_col]]),
  unique(emig_filtered[[emig_id_col]])
)

# Keep only ACC rows from individuals with a known emigration date
acc_classified <- acc_classified[get(acc_id_col) %chin% common_ids]
emig_filtered  <- emig_filtered[get(emig_id_col) %chin% common_ids]


#'==============================================================================
#' 1.5 Add dispersal start/end to ACC data
#'==============================================================================

# Create lookup vectors indexed by individual ID
dispersal_start_lookup <- emig_filtered$dispersal_start_utc
names(dispersal_start_lookup) <- emig_filtered[[emig_id_col]]

dispersal_end_lookup <- emig_filtered$dispersal_end_utc
names(dispersal_end_lookup) <- emig_filtered[[emig_id_col]]

# Add dispersal period to each ACC row
acc_classified[
  ,
  dispersal_start_utc := dispersal_start_lookup[get(acc_id_col)]
]

acc_classified[
  ,
  dispersal_end_utc := dispersal_end_lookup[get(acc_id_col)]
]


#'==============================================================================
#' 1.6 Retain only the first 15 weeks after dispersal
#'==============================================================================

# Keep records in the interval [dispersal_start_utc, dispersal_end_utc)
acc_15w <- acc_classified[
  !is.na(get(acc_time_col)) &
    get(acc_time_col) >= dispersal_start_utc &
    get(acc_time_col) < dispersal_end_utc
]

# Add age since dispersal in days
acc_15w[
  ,
  days_since_dispersal := as.numeric(
    difftime(get(acc_time_col), dispersal_start_utc, units = "days")
  )
]



#'==============================================================================
#' 1.8 Save filtered dataset
#'==============================================================================

saveRDS(
  acc_15w,
  file.path(output_dir, "acc_classified_first_15_weeks_all_individuals.rds"),
  compress = "gzip"
)

cat("Filtered ACC dataset saved.\n")







# Robust conversion of did_disperse to logical, without changing original column
did_raw <- emig_dt[[did_disperse_col]]

if (is.logical(did_raw)) {
  did_logical <- did_raw
} else if (is.numeric(did_raw)) {
  did_logical <- did_raw == 1
} else {
  did_logical <- tolower(trimws(as.character(did_raw))) %in%
    c("true", "t", "yes", "y", "1")
}

#' CONTROL 12: check did_disperse interpretation
#' Expected result:
#' - TRUE and FALSE values should be present or at least TRUE values.
#' - Number of TRUE should correspond to dispersing individuals.
cat("Number of dispersing individuals/rows:", sum(did_logical, na.rm = TRUE), "\n")
print(table(did_logical, useNA = "ifany"))

#' CONTROL 13: check dispersal_date format
#' Expected result: TRUE if already POSIXct.
#' If TRUE, conversion is skipped.
#' If FALSE, conversion to POSIXct UTC is necessary and will be applied.
dispersal_date_already_posix <- inherits(emig_dt[[dispersal_date_col]], "POSIXct")
cat("dispersal_date already POSIXct:", dispersal_date_already_posix, "\n")
cat("dispersal_date class:", class(emig_dt[[dispersal_date_col]]), "\n")

if (dispersal_date_already_posix) {
  emig_dt[, dispersal_start_utc := dispersal_date]
} else {
  emig_dt[, dispersal_start_utc := parse_datetime_utc(dispersal_date)]
}

# Define end of the first 15 weeks
# 15 weeks = 105 days
emig_dt[, dispersal_end_utc := dispersal_start_utc + 15 * 7 * 24 * 60 * 60]

# Keep only individuals that dispersed and have a valid dispersal date
emig_filtered <- emig_dt[
  did_logical == TRUE &
    !is.na(dispersal_start_utc)
]

#' CONTROL 14: check dispersal-date conversion
#' Expected result:
#' - no missing dispersal_start_utc among retained individuals
#' - dispersal_end_utc should be 105 days after dispersal_start_utc
cat("Number of retained emigration rows:", nrow(emig_filtered), "\n")
cat("Number of missing dispersal_start_utc:", sum(is.na(emig_filtered$dispersal_start_utc)), "\n")

print(
  emig_filtered[
    ,
    .(
      individual.id,
      dispersal_date,
      dispersal_start_utc,
      dispersal_end_utc
    )
  ][1:min(.N, 10)]
)

if (nrow(emig_filtered) == 0) {
  stop("No dispersing individual with valid dispersal date was retained.", call. = FALSE)
}


#'==============================================================================
#' 1.8 Check duplicate emigration dates per individual
#'==============================================================================

#' CONTROL 15: check duplicate dispersal dates
#' Expected result: empty table.
#' If empty, continue.
#' If not empty with different dates, inspect emig_dates before proceeding.
dup_date_control <- emig_filtered[
  ,
  .(
    n_rows = .N,
    n_unique_dispersal_dates = data.table::uniqueN(dispersal_start_utc)
  ),
  by = individual.id
][n_unique_dispersal_dates > 1]

cat("Individuals with multiple different dispersal dates:\n")
print(dup_date_control)

if (nrow(dup_date_control) > 0) {
  stop(
    "Some individuals have multiple different dispersal dates. Inspect emig_dates.",
    call. = FALSE
  )
}

# If duplicated rows have the same dispersal date, keep one row per individual.
emig_filtered <- unique(emig_filtered, by = "individual.id")


#'==============================================================================
#' 1.9 Retain only individuals present in both files
#'==============================================================================

acc_ids <- unique(acc_classified$IndividualID)
emig_ids <- unique(emig_filtered$individual.id)

common_ids <- intersect(acc_ids, emig_ids)

#' CONTROL 16: check overlap between ACC file and emigration file
#' Expected result: number > 0.
#' Ideally, this number should be close to the number of dispersing individuals
#' for which you have classified ACC data.
cat("Number of individuals in ACC file:", length(acc_ids), "\n")
cat("Number of dispersing individuals in emig_dates:", length(emig_ids), "\n")
cat("Number of common individuals:", length(common_ids), "\n")
print(common_ids)

if (length(common_ids) == 0) {
  stop(
    "No common IndividualID / individual.id between ACC file and emig_dates.",
    call. = FALSE
  )
}

# Keep only ACC rows from individuals with a known emigration date
acc_classified <- acc_classified[IndividualID %chin% common_ids]


#'==============================================================================
#' 1.10 Add dispersal start/end to ACC data without renaming ACC columns
#'==============================================================================

# Named lookup vectors
dispersal_start_lookup <- emig_filtered$dispersal_start_utc
names(dispersal_start_lookup) <- emig_filtered$individual.id

dispersal_end_lookup <- emig_filtered$dispersal_end_utc
names(dispersal_end_lookup) <- emig_filtered$individual.id

# Add dispersal dates to ACC rows
acc_classified[
  ,
  dispersal_start_utc := dispersal_start_lookup[IndividualID]
]

acc_classified[
  ,
  dispersal_end_utc := dispersal_end_lookup[IndividualID]
]

#' CONTROL 17: check that dispersal dates were assigned to ACC rows
#' Expected result:
#' - n_missing_dispersal_start_utc = 0
#' - n_missing_dispersal_end_utc = 0
assignment_control <- acc_classified[
  ,
  .(
    n_rows = .N,
    n_missing_dispersal_start_utc = sum(is.na(dispersal_start_utc)),
    n_missing_dispersal_end_utc = sum(is.na(dispersal_end_utc))
  )
]

print(assignment_control)

if (
  assignment_control$n_missing_dispersal_start_utc > 0 |
  assignment_control$n_missing_dispersal_end_utc > 0
) {
  stop("Some ACC rows did not receive dispersal dates. Check ID matching.", call. = FALSE)
}


#'==============================================================================
#' 1.11 Retain only the first 15 weeks after dispersal
#'==============================================================================

# Keep records in the interval:
# [dispersal_start_utc, dispersal_end_utc)
acc_15w <- acc_classified[
  !is.na(timestamp) &
    timestamp >= dispersal_start_utc &
    timestamp < dispersal_end_utc
]

# Stop if no data remain
if (nrow(acc_15w) == 0) {
  stop("No ACC rows were retained within the first 15 weeks of dispersal.", call. = FALSE)
}

# Add age since dispersal, in days
acc_15w[
  ,
  days_since_dispersal := as.numeric(
    difftime(timestamp, dispersal_start_utc, units = "days")
  )
]

# Free memory
rm(acc_classified)
gc()


#'==============================================================================
#' 1.12 Controls after filtering to first 15 weeks
#'==============================================================================

#' CONTROL 18: global check of the 15-week filter
#' Expected result:
#' - min_days_since_dispersal >= 0
#' - max_days_since_dispersal < 105
#' - n_individuals > 0
global_control <- acc_15w[
  ,
  .(
    n_individuals = data.table::uniqueN(IndividualID),
    n_rows = .N,
    min_days_since_dispersal = min(days_since_dispersal, na.rm = TRUE),
    max_days_since_dispersal = max(days_since_dispersal, na.rm = TRUE),
    n_missing_timestamp = sum(is.na(timestamp))
  )
]

print(global_control)

if (global_control$min_days_since_dispersal < 0) {
  stop("Some rows occur before dispersal_start_utc. Filtering failed.", call. = FALSE)
}

if (global_control$max_days_since_dispersal >= 105) {
  stop("Some rows occur after the first 15 weeks. Filtering failed.", call. = FALSE)
}


#' CONTROL 19: individual-level check
#' Expected result:
#' - all min_days_since_dispersal >= 0
#' - all max_days_since_dispersal < 105
#' - n_rows > 0 for retained individuals
individual_control <- acc_15w[
  ,
  .(
    n_rows = .N,
    first_timestamp = min(timestamp, na.rm = TRUE),
    last_timestamp = max(timestamp, na.rm = TRUE),
    min_days_since_dispersal = min(days_since_dispersal, na.rm = TRUE),
    max_days_since_dispersal = max(days_since_dispersal, na.rm = TRUE)
  ),
  by = IndividualID
][order(IndividualID)]

print(individual_control)

# Save individual-level control table
data.table::fwrite(
  individual_control,
  file.path(output_dir, "control_first_15_weeks_by_IndividualID.csv")
)


#' CONTROL 20: check behavioural classes after filtering
#' Expected result:
#' - rf8fitted should contain expected behavioural classes
#' - no unexpected empty or NA-dominated classification
behavior_control <- acc_15w[
  ,
  .N,
  by = rf8fitted
][order(-N)]

print(behavior_control)

# Save behavioural control table
data.table::fwrite(
  behavior_control,
  file.path(output_dir, "control_rf8fitted_first_15_weeks.csv")
)


#' CONTROL 21: check RF confidence after filtering
#' Expected result:
#' - pro_rf8fitted values should be between 0 and 1
#' - missing values should be rare or absent
prob_control_15w <- acc_15w[
  ,
  .(
    min_pro_rf8fitted = min(pro_rf8fitted, na.rm = TRUE),
    mean_pro_rf8fitted = mean(pro_rf8fitted, na.rm = TRUE),
    median_pro_rf8fitted = median(pro_rf8fitted, na.rm = TRUE),
    max_pro_rf8fitted = max(pro_rf8fitted, na.rm = TRUE),
    n_missing_pro_rf8fitted = sum(is.na(pro_rf8fitted))
  )
]

print(prob_control_15w)

data.table::fwrite(
  prob_control_15w,
  file.path(output_dir, "control_pro_rf8fitted_first_15_weeks.csv")
)











#'==============================================================================
#' 1.13 Save filtered classified ACC dataset
#'==============================================================================

# Save as RDS to preserve column types, especially POSIXct timestamps.
out_rds <- file.path(
  output_dir,
  "acc_classified_first_15_weeks_all_individuals.rds"
)

saveRDS(acc_15w, out_rds, compress = "gzip")

cat("Filtered ACC dataset saved here:\n")
cat(out_rds, "\n")