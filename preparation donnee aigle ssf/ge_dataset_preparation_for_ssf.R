#' ---
#' title: "Preparing Golden Eagle dataset for ssf"
#' author: "Louise Faure"
#' date: 02.06.2026
#' details: (i) select all the golden eagles gps values for the first 15 weeks of the dispersal phase, 
#' (ii) define a threshold value to keep the location close to the group
#' (iii) identify the time at which each location become independent    
#' ---   


#' #library
library(dplyr)          # manipulation de data.frames (filter, mutate, etc.)
library(tidyr)          # nest / unnest
library(purrr)          # map / map_dfr
library(lubridate)      # floor_date, hours(), minutes()
library(sf)             # st_as_sf, st_buffer, st_transform
library(terra)          # rast, vect, project
library(amt)            # make_track, track_resample, steps_by_burst, random_steps, fit_distr
library(EMbC)           # embc, smth
#library(CircStats)      # von Mises (utilisé par amt en interne)
#library(circular)       # statistiques circulaires (idem)
#library(fitdistrplus)   # ajustement de distributions (idem)
#library(units)          # gestion des unités (idem)
library(exactextractr)  # NOUVEAU : extraction fractionnelle dans un buffer


#' ### STEP 0 : load the data --------------------------------------------------

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



#' ### STEP 1 : filtration and calculation of height values --------------------

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


#' ### STEP 2 : behavioral classification --------------------------------------

#' Ref. Garriga et al., 2016, PLOS One

# bivariate matrix
behavioural_classification <- data.matrix(dispersal_data %>%
  select(ground.speed, height_above_ground))

# call embc
embc <- embc(behavioural_classification)

# investigate the bc
X11(); sctr (bc)

embc_smoothed <- smth(embc, dlta = 0.7)
X11();sctr(bc_smth)

# link the cluster labels (1: low speed, low height, 2: low speed, height height, 3: heigh speed, low height, 4: height speed, heigh height) to the original data, e.i, dispersal_data
dispersal_data$behavior_cluster <- embc_smoothed@A

# check that cluster 1 correspond to terrestrial behaviors
dispersal_data %>%
  group_by(behavior_cluster) %>%
  summarise(speed_med  = median(ground.speed,        na.rm = TRUE),
            height_med = median(height_above_ground, na.rm = TRUE),
            n          = n())
# terrestrial behaviour should have a speed ~0, and a height ~0 



#' ### STEP 3 : step selection preparation - generate alternative steps --------

#' **Philosophy:** we retain all gps location associated to terrestrial and non-terrestrial 
#' behaviours. It is only later that we will consider in the definition of the landscape
#' caracteristics used that we will consider the steps that end with a terrestrial behaviour

# reproject the data to a metric crs
crs_metric <- "+proj=utm +zone=32 +datum=WGS84 +units=m +no_defs"
dd_sf <- st_as_sf(dispersal_data,
                  coords = c("location.long", "location.lat"),
                  crs    = 4326) %>%
  st_transform(crs_metric)
dispersal_data$x_m <- st_coordinates(dd_sf)[, 1]
dispersal_data$y_m <- st_coordinates(dd_sf)[, 2]

# create a track
gps_track <- make_track(dispersal_data,
                        x_m, y_m, timestamp,
                        id               = individual.local.identifier,
                        behavior_cluster = behavior_cluster,   # on conserve l'étiquette
                        crs              = 3035)

# hourly re sampling of burst 
steps_by_individual <- gps_track %>%
  nest(data = -id) %>%
  mutate(
    data_hourly = map(data, ~ track_resample(.x,
                                             rate      = hours(1),
                                             tolerance = minutes(10)) %>% # 1 hour + 10 minutes tolerance
                        filter_min_n_burst(min_n = 3)), # remove burst with less than 3 points
    steps = map(data_hourly, steps_by_burst, keep_cols = "end")  # calculate turning angle and step lenghts
  ) %>%
  select(id, steps) %>%
  unnest(cols = steps)

# diagnostic
nrow(steps_by_individual)                            # nombre total de pas
table(steps_by_individual$behavior_cluster)          # distribution des clusters d'ARRIVÉE


#' ### STEP 4 : estimate turning angles and step length distributions -----------

# gamma sur les longueurs de pas
sl_fit <- fit_distr(steps_by_individual$sl_, "gamma")

# von Mises sur les angles de virage
ta_fit <- fit_distr(steps_by_individual$ta_, "vonmises")

# inspection visuelle
par(mfrow = c(1, 2))
hist(steps_by_individual$sl_, breaks = 50, freq = FALSE,
     xlab = "Step length (m)", main = "")
curve(dgamma(x, shape = sl_fit$params$shape, scale = sl_fit$params$scale),
      add = TRUE, col = "red")
hist(steps_by_individual$ta_, breaks = 30, freq = FALSE,
     xlab = "Turning angle (rad)", main = "")
curve(dvonmises(x, mu = ta_fit$params$mu, kappa = ta_fit$params$kappa),
      add = TRUE, col = "red")



#' ### STEP 5 : generate alternative steps -------------------------------------

#' **Philosophy:** us is defined when the end of a step correspond to a terrestrial 
#' behavior

# filter steps by behaviors
used_steps <- steps_by_individual %>%
  filter(behavior_cluster == 1)         

# amt::random_steps
ssf_data <- used_steps %>%
  random_steps(n_control = 50,          # 50 random steps per observed steps
               sl_distr  = sl_fit,
               ta_distr  = ta_fit) %>%
  # stratum id (used + 50 random non-used = 1 stratum)
  mutate(stratum_id = paste(id, burst_, step_id_, sep = "_"))

# vérification : chaque strate doit contenir 1 used (case_ = TRUE) et 50 available
ssf_data %>% group_by(stratum_id) %>%
  summarise(n_used = sum(case_), n_avail = sum(!case_)) %>%
  count(n_used, n_avail)                # attendu : n_used=1, n_avail=50 partout



#' ### STEP 6: annotation with topographic and human info ----------------------

#' For each end point (used of available), we extract the raster values from the 
#' rasters

# open raster layers
covar_1_elevation <- terra::rast("CHEMIN/VERS/DIST_TO_BUILDINGS.tif")
covar_2_distance_ridgeline <- terra::rast("C:/Users/lfaure7/OneDrive/MEMOIRE M2/donnees/raster/topography/distance_to_ridge_line_complete_version.tif")
covar_3_TRI <- terra::rast("C:/Users/lfaure7/OneDrive/MEMOIRE M2/eagle_projet/data/pretraitements/DEM_25m/TRI/TRI.tif")
covar_4_slope <- terra::rast("C:/Users/lfaure7/OneDrive/MEMOIRE M2/eagle_projet/data/pretraitements/DEM_25m/slope/slope_25.tif")
covar_5_TPI <- terra::rast("C:/Users/lfaure7/OneDrive/MEMOIRE M2/eagle_projet/data/pretraitements/DEM_25m/TPI/Topographic Position Index.tif")
covar_6_distance_settlement <- terra::rast("C:/Users/lfaure7/Desktop/COUCHES QGIS/settlements/distance to settlements 30 m cropped/dist_to_settlements_30m_croped.tif")
covar_7_density_settlement <- terra::rast("C:/Users/lfaure7/Desktop/COUCHES QGIS/settlements/density_rad_495m.tif")
covar_8_landcover <- terra::rast("C:/Users/lfaure7/OneDrive/THESE/Conférences/Sempach workshop/Donnees/Landscape layers/CLC_longlat_10m.tif")

#reproject tracking data to match raster crs, extract values, convert back to wgs and save as a dataframe
ssf_annotated <- ssf_data %>%
  mutate(
    dist_buildings = terra::extract(
      covar_raster_1,
      terra::vect(cbind(x2_, y2_), crs = "EPSG:3035") |>      # point d'arrivée
        terra::project(crs(covar_raster_1))
    )[, 2],
    TRI_100 = terra::extract(
      covar_raster_2,
      terra::vect(cbind(x2_, y2_), crs = "EPSG:3035") |>
        terra::project(crs(covar_raster_2))
    )[, 2]
  )

# export as a csv. 