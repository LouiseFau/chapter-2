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


# 1. Tracks exploration ----
# Ensure timestamps are ordered within tracks
locs <- locs %>%
  arrange(mt_track_id(.), mt_time(.))

# Ensure geographic CRS 
# for calculating spatial metrics
locs <- locs %>%
  st_transform(3035)

## 1.1 Time lags ----
# move2
time_lag <- mt_time_lags(locs, units = "min")  # minutes
summary(time_lag)

# visualise but cutting above 30min 
# for graphical ease
ggplot(tibble(track = mt_track_id(locs), time_lag = as.numeric(time_lag))) +
  geom_histogram(aes(time_lag), bins = 200) +
  labs(x = "Time lag (min)", y = "Count", title = "Distribution of time lags") +
  xlim(0, 30) +
  theme_minimal()



################################################################################
#' Step 2. Data cleaning
#' 
#' Identify two scale of temporal resolution and handles gap in the dataset
#' #############################################################################

locs <- mt_filter_unique(locs)


# Handle gaps: split bursts so steps do not cross long missing periods
# Decide a maximum acceptable gap (species + study design).

# Regularize to a fixed timestep (common prerequisite for step-based analyses) 
# This does NOT "invent" positions; it standardizes timestamps and flags missingness.
# Choose a timestep consistent with your data density (e.g., 30 min).
# using lubridrate


