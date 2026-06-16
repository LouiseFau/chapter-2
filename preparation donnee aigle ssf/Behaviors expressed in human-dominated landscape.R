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


################################################################################
#' ### STEP 0 : load the data --------------------------------------------------
################################################################################


output_dir <- "C:/Users/lfaure7/Documents/git/chapter-2/preparation donnee aigle ssf/resultats-intermediaire"

# Géoïde and digital elevation model
geo <- terra::rast("C:/Users/lfaure7/OneDrive/MEMOIRE M2/eagle_projet/data/geoide/us_nga_egm96_15.tif")
dem <- terra::rast("C:/Users/lfaure7/OneDrive/MEMOIRE M2/eagle_projet/data/pretraitements/Region-Alpes-Dem/region-alpes-dem.tif")

# raster layers
distance_settlement <- terra::rast("C:/Users/lfaure7/Desktop/COUCHES QGIS/settlements/distance to settlements 30 m cropped/dist_to_settlements_30m_croped.tif")
density_settlement <- terra::rast("C:/Users/lfaure7/Desktop/COUCHES QGIS/density-resident-pop-km2/density_pop_km2_meters/density_settlement_km2_meters_crs.tif")
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
#' We can therefore aggregate all these raster layer in one. Here are the correlation 
#' coefficient obtained. 
#' 
#'              dist_settlement dens_settlement dens_pop_km2
#'dist_settlement   1.0000000      -0.8647500   -0.8249338
#'dens_settlement  -0.8647500       1.0000000    0.8762607
#'dens_pop_km2     -0.8249338       0.8762607    1.0000000
#' 
#' (1) create human footprint raster
#' (2) extract location points in the vicinity of human infrastructure using at
#' least to different scale
#' (3) inspect the number of location retained per individuals
################################################################################

# Transform

# distance_settlement : inverser (proximité = 1/distance)
# +1 pour éviter division par zéro quand dist = 0
prox_settlement <- 1 / (distance_settlement + 1)

# log-transformer les densités pour réduire l'effet des extrêmes urbains
# (Genève ou Turin ne doivent pas écraser toute la variabilité alpine)
dens_s_log <- log1p(density_settlement_r)
dens_p_log <- log1p(density_pop_km2_r)

# --- A3 : normaliser chaque couche entre 0 et 1 ---
# (valeur - min) / (max - min)
# nécessaire pour rendre les trois couches comparables avant combinaison

normalize <- function(r) {
  mn <- terra::global(r, "min", na.rm = TRUE)[[1]]
  mx <- terra::global(r, "max", na.rm = TRUE)[[1]]
  (r - mn) / (mx - mn)
}

prox_norm  <- normalize(prox_settlement)
dens_s_norm <- normalize(dens_s_log)
dens_p_norm <- normalize(dens_p_log)

# --- A4 : combiner en HFI (moyenne des trois composantes) ---
# poids égaux : hypothèse conservative, à mentionner dans la méthode
hfi <- (prox_norm + dens_s_norm + dens_p_norm) / 3

names(hfi) <- "hfi"

# --- A5 : masquer sur les Alpes ---
alpes_v <- terra::vect(
  sf::st_read("C:/Users/lfaure7/Desktop/COUCHES QGIS/alpine_area/alpine_area_corrected_geometry.shp")
) %>% terra::project(terra::crs(distance_settlement))

hfi_masked <- terra::crop(hfi, alpes_v) %>%
  terra::mask(alpes_v)

terra::writeRaster(hfi_masked,
                   file.path(output_dir, "human_footprint_index_alpes.tif"),
                   datatype  = "FLT4S",
                   overwrite = TRUE)

################################################################################
# STEP B : extraire le HFI à chaque point GPS et définir les seuils
################################################################################

# --- B1 : projeter les points GPS dans le CRS du HFI ---
crs_hfi <- terra::crs(hfi_masked)

pts_proj <- dispersal_data %>%
  st_as_sf(coords = c("location.long", "location.lat"),
           crs = 4326, remove = FALSE) %>%
  st_transform(crs_hfi) %>%
  terra::vect()

# --- B2 : extraire la valeur HFI et les covariables brutes à chaque point ---
# on garde aussi les covariables brutes pour pouvoir les utiliser
# comme descripteurs des points retenus (ta deuxième approche)

dispersal_data$hfi              <- terra::extract(hfi_masked,          pts_proj)[, 2]
dispersal_data$dist_settlement  <- terra::extract(distance_settlement,  pts_proj)[, 2]
dispersal_data$dens_settlement  <- terra::extract(density_settlement_r, pts_proj)[, 2]
dispersal_data$dens_pop_km2     <- terra::extract(density_pop_km2_r,    pts_proj)[, 2]

# --- B3 : définir les seuils de proximité sur le HFI ---
# le HFI est entre 0 (aucune empreinte humaine) et 1 (maximum)
# les seuils sont des quantiles de la distribution du HFI dans les Alpes
# ce qui est plus défendable que des valeurs arbitraires absolues :
# "proche des humains" = dans le quartile supérieur du HFI alpin

q_hfi <- terra::global(hfi_masked, fun = quantile,
                       probs = c(0.50, 0.75, 0.90),
                       na.rm = TRUE)
print(q_hfi)  # vérifie les valeurs avant de fixer les seuils

# trois niveaux de "proximité humaine" basés sur les quantiles du HFI
dispersal_data <- dispersal_data %>%
  dplyr::mutate(
    near_human_q50 = hfi >= q_hfi[1, 1],   # > médiane du HFI alpin
    near_human_q75 = hfi >= q_hfi[2, 1],   # > 3ème quartile
    near_human_q90 = hfi >= q_hfi[3, 1]    # > 90ème percentile
  )

################################################################################
# STEP C : inspecter le nombre de points retenus par individu et par seuil
################################################################################

sample_summary <- dispersal_data %>%
  dplyr::group_by(individual.id) %>%
  dplyr::summarise(
    n_total    = n(),
    n_q50      = sum(near_human_q50, na.rm = TRUE),
    n_q75      = sum(near_human_q75, na.rm = TRUE),
    n_q90      = sum(near_human_q90, na.rm = TRUE),
    pct_q50    = round(100 * n_q50 / n_total, 1),
    pct_q75    = round(100 * n_q75 / n_total, 1),
    pct_q90    = round(100 * n_q90 / n_total, 1),
    .groups    = "drop"
  )

print(sample_summary)

# visualisation : nombre de points par individu et par seuil
plot_data <- sample_summary %>%
  dplyr::select(individual.id, n_q50, n_q75, n_q90) %>%
  tidyr::pivot_longer(-individual.id,
                      names_to  = "threshold",
                      values_to = "n_fixes") %>%
  dplyr::mutate(threshold = factor(threshold,
                                   levels = c("n_q50", "n_q75", "n_q90"),
                                   labels = c("HFI > Q50", "HFI > Q75", "HFI > Q90")))

ggplot(plot_data,
       aes(x = reorder(factor(individual.id), -n_fixes), y = n_fixes)) +
  geom_col(fill = "steelblue") +
  facet_wrap(~ threshold, ncol = 1, scales = "free_y") +
  labs(x     = "Individu",
       y     = "Nombre de fixes GPS",
       title = "Points GPS en zone à forte empreinte humaine") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

saveRDS(dispersal_data,
        file.path(output_dir, "dispersal_data_hfi.rds"))









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
