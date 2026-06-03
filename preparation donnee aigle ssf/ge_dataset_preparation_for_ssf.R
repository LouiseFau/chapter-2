#' ---
#' title: "Preparing Golden Eagle dataset for ssf"
#' author: "Louise Faure"
#' date: 02.06.2026
#' details: (i) select all the golden eagles gps values for the first 15 weeks of the dispersal phase, 
#' (ii) define a threshold value to keep the location close to the group
#' (iii) identify the time at which each location become independent    
#' ---   


#' #library
library(dplyr)
library(terra)
library(move)

#' # STEP 0 : load the data ----------------------------------------------------
output_dir <- "C:/Users/lfaure7/Documents/git/chapter-2/preparation donnee aigle ssf/resultats-intermediaire"

# Géoïde and digital elevation model
geo <- terra::rast("C:/Users/lfaure7/OneDrive/MEMOIRE M2/eagle_projet/data/geoide/us_nga_egm96_15.tif")
dem <- terra::rast("C:/Users/lfaure7/OneDrive/MEMOIRE M2/eagle_projet/data/pretraitements/Region-Alpes-Dem/region-alpes-dem.tif")

# emigration dates and golden eagle data
emig_dates <- readRDS("C:/Users/lfaure7/Documents/git/chapter-2/preparation donnee aigle ssf/donnees/emigration dates/emigration_dates_20250417.rds")

rds_dir <- "C:/Users/lfaure7/Documents/git/chapter-2/preparation donnee aigle ssf/donnees/no_burst_GE/no_burst_GE"
rds_files <- tools::file_path_sans_ext(list.files(rds_dir, pattern = "\\.rds$", ignore.case = TRUE))
rds_ids   <- as.numeric(gsub("_gpsNoDup_moveObj", "", rds_files))

# Filters individuals to retain those that have an emigration date
emig_dates_filtered <- emig_dates[emig_dates$did_disperse == TRUE &
                                    emig_dates$individual.id %in% rds_ids, ]
rds_files_filtered  <- rds_files[rds_ids %in% emig_dates_filtered$individual.id]


#' # STEP 1 : filtration and calculation of height values
#' (i) keep only GPS locations within the first 15 weeks after emigration date
#' (ii) compute height above ground (geoid + DEM correction)
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

#'STEP 2 : sensibility test to identify location points associated with terrestrial behaviors (resting, walking, running)
#'(i) identification of a speed threshold 
#'(ii) identification of an height above the ground threshold
grid <- expand.grid(v = c(1.5, 2, 3), h = c(30, 50, 100))
grid$n <- mapply(function(v, h)
  sum(dispersal_data$ground.speed < v & dispersal_data$height_above_ground < h, na.rm = TRUE),
  grid$v, grid$h)
grid 

# We find that at constant height, the nbr of gps point does not change substantially, while at constant speed, the nbr of gps points varies.
# We conclude that the threshold of speed has little influenced on the selection of location points. 
# To define a height threshold, we explore the dispersal of height values for location that are certainly stationary (e.i., < 0.5 m/s)

dispersal_data %>%
  filter(ground.speed < 0.5) %>%  
  pull(height_above_ground) %>%
  quantile(c(0.5, 0.9, 0.95, 0.99), na.rm = TRUE) # 90% of height values are below 59m, we can therefore keep the height threshold at 50m

dispersal_terrestrial <- dispersal_data %>%
  filter(!is.na(ground.speed),
         ground.speed        < 2,
         height_above_ground < 50)

#'STEP 3: identify the time at which each location point become independant 
#'
