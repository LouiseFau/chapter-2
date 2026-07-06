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
missing_eagle_classified <- 

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


# 1.1 Prepare emigration dates ----
# Convert emig_dates to data.table
emig_dt <- data.table::as.data.table(emig_dates)

# Column names in emig_dates
emig_id_col        <- "individual.local.identifier"
did_disperse_col   <- "did_disperse"
dispersal_date_col <- "dispersal_date"




# 1.2 Keep only individuals with a confirmed dispersal date ----
did_disperse == TRUE

# Keep only individuals that dispersed and have a dispersal date
emig_filtered <- emig_dt[
  did_logical == TRUE &
    !is.na(get(dispersal_date_col))
]



#' 1.3 Prepare emigration dates and match them with ACC data ----
# Convert individual identifiers to character in both files.
acc_classified[, individualID := as.character(individualID)]
emig_dt[, individual.local.identifier := as.character(individual.local.identifier)]

# Keep only individuals that dispersed and have a valid dispersal date.
emig_filtered <- emig_dt[
  did_disperse == TRUE & !is.na(dispersal_date),
  list(
    individual.local.identifier = individual.local.identifier,
    dispersal_date = dispersal_date
  )
]



#' 1.4 Add dispersal date to ACC data ----
# match acc_classified$individualID with emig_filtered$individual.local.identifier.
acc_classified[
  emig_filtered,
  dispersal_date := i.dispersal_date,
  on = c("individualID" = "individual.local.identifier")
]

# Keep only ACC rows belonging to individuals with a valid dispersal date.
acc_classified <- acc_classified[!is.na(dispersal_date)]


#' 15 Create age since emigration and retain the first 15 weeks ----
acc_classified[
  ,
  days_since_emig := as.numeric(
    difftime(timestamp, dispersal_date, units = "days")
  )
]

# Retain only rows between departure and 105 days after departure.
# 15 weeks = 105 days.
acc_15w <- acc_classified[
  days_since_emig >= 0 &
    days_since_emig < 105
]


#' 1.6 Save filtered dataset ----
saveRDS(
  acc_15w,
  file.path(output_dir, "acc_classified_first_15_weeks_all_individuals.rds"),
  compress = "gzip"
)