#' ---
#' title: "Behaviors expressed in human-dominated landscape"
#' author: "Louise Faure"
#' date of creation: 16.06.2026
#' details: (i) select all the golden eagles gps values for the first 15 weeks of the dispersal phase, 
#' (ii) extract location points in the vicinity of human infrastructures
#' (iii) identify movement characteristics (flight height, turning angle, speed) and attribute a behavior  
#' ---   

#' #library
library(dplyr)       # data manipulation
library(tidyr)       # pivot_longer for visualisation
library(sf)          # spatial operations (st_as_sf, st_buffer, st_transform)
library(terra)       # raster operations (rast, extract, project)
library(ggplot2)     # visualisation
library(move)

################################################################################
#' ### STEP 0 : load the data --------------------------------------------------
################################################################################


output_dir <- "C:/Users/lfaure7/Documents/git/chapter-2/preparation donnee aigle ssf/resultats-intermediaire"

# Géoïde and digital elevation model
geo <- terra::rast("C:/Users/lfaure7/OneDrive/MEMOIRE M2/eagle_projet/data/geoide/us_nga_egm96_15.tif")
dem <- terra::rast("C:/Users/lfaure7/OneDrive/MEMOIRE M2/eagle_projet/data/pretraitements/Region-Alpes-Dem/region-alpes-dem.tif")

# raster layers
distance_settlement <- terra::rast("C:/Users/lfaure7/Desktop/COUCHES QGIS/settlements/distance to settlements 30 m cropped/dist_to_settlements_30m_croped.tif")
density_settlement <- terra::rast("C:/Users/lfaure7/Desktop/COUCHES QGIS/density-resident-pop-km2/density_settlement_km2_m/density_settlement_km2_meters_crs.tif")
density_pop_km2 <- terra::rast("C:/Users/lfaure7/Desktop/COUCHES QGIS/density-resident-pop-km2/density_pop_km2_m/density_pop_km2_meters_crs.tif")
landcover <- terra::rast("C:/Users/lfaure7/OneDrive/THESE/Conférences/Sempach workshop/Donnees/Landscape layers/CLC_longlat_10m.tif")


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
#' ### STEP 1 : filtration and calculation of height values --------------------

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

################################################################################
#' ### Step 2 : extract location points in the vicinity of human infrastructures
#' 
#' **Philosophy:** preliminary analysis revealed that distance to settlement, 
#' human density and settlement density were highly correlated (see annex 1).
#' We can therefore aggregate all these raster layer in one. Here are the 
#' correlation coefficient obtained. 
#' 
#'              dist_settlement dens_settlement dens_pop_km2
#'dist_settlement   1.0000000      -0.8647500   -0.8249338
#'dens_settlement  -0.8647500       1.0000000    0.8762607
#'dens_pop_km2     -0.8249338       0.8762607    1.0000000
#' 
#' (1) import the human footprint raster which has been created following the
#' steps in annex 2
#' (2) extract location points in the vicinity of human footprint using at
#' least to different scale
#' (3) inspect the number of location retained per individuals
################################################################################

hfi_masked <- terra::rast("C:/Users/lfaure7/Desktop/COUCHES QGIS/Human footprint index/human_footprint_index.tif")
crs_hfi <- terra::crs(hfi_masked) # project point in the crs of hfi

pts_proj <- dispersal_data %>%
  st_as_sf(coords = c("location.long", "location.lat"),
           crs = 4326, remove = FALSE) %>%
  st_transform(crs_hfi) %>%
  terra::vect()

# extract covariate
dispersal_data$hfi              <- terra::extract(hfi_masked,          pts_proj)[, 2]
dispersal_data$dist_settlement  <- terra::extract(distance_settlement,  pts_proj)[, 2]
dispersal_data$dens_settlement  <- terra::extract(density_settlement, pts_proj)[, 2]
dispersal_data$dens_pop_km2     <- terra::extract(density_pop_km2,    pts_proj)[, 2]

# define threshold in meters
dispersal_data <- dispersal_data %>%
  dplyr::mutate(
    near_human_100m = !is.na(dist_settlement) & dist_settlement <= 100,
    near_human_500m = !is.na(dist_settlement) & dist_settlement <= 500,
    near_human_1km  = !is.na(dist_settlement) & dist_settlement <= 1000
  )

# inspect the nbr of point per meter threshold
dispersal_data %>%
  dplyr::summarise(
    n_total     = n(),
    n_100m      = sum(near_human_100m, na.rm = TRUE),
    n_500m      = sum(near_human_500m, na.rm = TRUE),
    n_1km       = sum(near_human_1km,  na.rm = TRUE),
    pct_100m    = round(100 * n_100m / n_total, 1),
    pct_500m    = round(100 * n_500m / n_total, 1),
    pct_1km     = round(100 * n_1km  / n_total, 1)
  ) %>% print()

# inspect the nbr of point per individual and threshold

sample_summary <- dispersal_data %>%
  dplyr::group_by(individual.local.identifier) %>%
  dplyr::summarise(
    n_total  = n(),
    n_100m   = sum(near_human_100m, na.rm = TRUE),
    n_500m   = sum(near_human_500m, na.rm = TRUE),
    n_1km    = sum(near_human_1km,  na.rm = TRUE),
    pct_100m = round(100 * n_100m / n_total, 1),
    pct_500m = round(100 * n_500m / n_total, 1),
    pct_1km  = round(100 * n_1km  / n_total, 1),
    .groups  = "drop"
  )

print(sample_summary, n = 66)

# Plot the results
plot_data <- sample_summary %>%
  dplyr::select(individual.local.identifier, n_100m, n_500m, n_1km) %>%
  tidyr::pivot_longer(-individual.local.identifier,
                      names_to  = "threshold",
                      values_to = "n_fixes") %>%
  dplyr::mutate(threshold = factor(threshold,
                                   levels = c("n_100m", "n_500m", "n_1km"),
                                   labels = c("≤ 100 m", "≤ 500 m", "≤ 1 km")))

ggplot(plot_data,
       aes(x     = reorder(individual.local.identifier, -n_fixes),
           y     = n_fixes,
           fill  = threshold)) +
  geom_col(position = "dodge") +
  facet_wrap(~ threshold, ncol = 1, scales = "free_y") +
  scale_fill_brewer(palette = "Blues", guide = "none") +
  labs(x     = "Individu",
       y     = "Nombre de fixes GPS",
       title = "Points GPS par proximité aux settlements") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

################################################################################
#' ### STEP 3 : movement caracteristics in the vicinity of humans 
################################################################################









################################################################################
#' ### ANNEX 1 : correlation coefficent for anthropogenic features
################################################################################

# raster de référence : définit la grille cible
ref <- distance_settlement

# reprojeter et rééchantillonner sur la grille de référence
# method = "bilinear" pour les valeurs continues
density_settlement_r <- terra::project(density_settlement, ref, method = "bilinear")
density_pop_km2_r    <- terra::project(density_pop_km2,    ref, method = "bilinear")

# vérification
terra::ext(ref)
terra::ext(density_settlement_r)
terra::ext(density_pop_km2_r)

# maintenant tu peux les combiner
samp <- terra::spatSample(c(ref,
                            density_settlement_r,
                            density_pop_km2_r),
                          size   = 10000,
                          method = "random",
                          na.rm  = TRUE)

names(distance_settlement)   <- "dist_settlement"
names(density_settlement_r)  <- "dens_settlement"
names(density_pop_km2_r)     <- "dens_pop_km2"

cor(samp, method = "spearman")




################################################################################
#' ### ANNEX 2 : prepare Human footprint index 
################################################################################

# transform distance to proximity and reduce extrem value by using log
prox_settlement <- 1 / (distance_settlement + 1)
dens_s_log <- log1p(density_settlement)
dens_p_log <- log1p(density_pop_km2)

# normalisation
normalize <- function(r) {
  mn <- terra::global(r, "min", na.rm = TRUE)[[1]]
  mx <- terra::global(r, "max", na.rm = TRUE)[[1]]
  (r - mn) / (mx - mn)
}

prox_norm  <- normalize(prox_settlement)
dens_s_norm <- normalize(dens_s_log)
dens_p_norm <- normalize(dens_p_log)

# combine
hfi <- (prox_norm + dens_s_norm + dens_p_norm) / 3

names(hfi) <- "hfi"

# crop
alpes_v <- terra::vect(
  sf::st_read("C:/Users/lfaure7/Desktop/COUCHES QGIS/alpine_area/alpine_area_corrected_geometry.shp")
) %>% terra::project(terra::crs(distance_settlement))

hfi_masked <- terra::crop(hfi, alpes_v) %>%
  terra::mask(alpes_v)
