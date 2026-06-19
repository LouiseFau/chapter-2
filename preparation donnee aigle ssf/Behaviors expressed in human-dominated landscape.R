#' ---
#' title: "Behaviors expressed in human-dominated landscape"
#' author: "Louise Faure"
#' date of creation: 16.06.2026
#' details: (i) select all the golden eagles gps values for the first 15 weeks of the dispersal phase, 
#' (ii) classify behavior and applied a spatial thinning to the gps points
#' (iii) identify behavior characteristics in human-dominated landscape  
#' ---   

#' #library
library(dplyr)       # data manipulation
library(tidyr)       # pivot_longer for visualisation
library(sf)          # spatial operations (st_as_sf, st_buffer, st_transform)
library(terra)       # raster operations (rast, extract, project)
library(ggplot2)     # visualisation
library(move)

library(EMbC)
library(lubridate)
library(exactextractr)

################################################################################
#' ### STEP 0 : load the data --------------------------------------------------
################################################################################


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
#' ### Step 2 : behavioral classification
#' 
#' ** Strategy**: (1) classify behaviors using embc function (2) sub classify 
#' terrestrial behavior to distinguish between overnight encroachment and short 
#' term resting 
################################################################################

# --- 2a : EMbC classification -------------------------------------------------
input_vars   <- c("ground.speed", "height_above_ground")
complete_idx <- which(complete.cases(dispersal_data[, input_vars]))
bc_matrix    <- data.matrix(dispersal_data[complete_idx, input_vars])

embc_fit      <- EMbC::embc(bc_matrix)
embc_smoothed <- EMbC::smth(embc_fit, dlta = 0.7)

# complete embc cluster with turning angles
cluster_summary_extended <- dispersal_data[complete_idx, ] %>%
  dplyr::mutate(behavior_cluster = embc_smoothed@A) %>%
  dplyr::group_by(behavior_cluster) %>%
  dplyr::summarise(
    speed_med      = median(ground.speed,        na.rm = TRUE),
    height_med     = median(height_above_ground, na.rm = TRUE),
    turn_angle_med = median(abs(turn.angle),     na.rm = TRUE),
    vert_speed_med = median(abs(vert.speed),     na.rm = TRUE),
    n = n(), .groups = "drop")
print(cluster_summary_extended)

dispersal_data$behavior_cluster <- NA_integer_
dispersal_data$behavior_cluster[complete_idx] <- embc_smoothed@A

# adjust cluster ids below if cluster_summary shows a different ordering
behavior_labels <- c("1" = "terrestrial",
                     "2" = "low soaring",
                     "3" = "fast flight at low elevation",
                     "4" = "fast commuting flight at high elevation")

dispersal_data <- dispersal_data %>%
  dplyr::mutate(behavior = dplyr::recode(as.character(behavior_cluster),
                                         !!!behavior_labels))

# --- 2b : sub-classification of terrestrial bouts ----------------------------
dispersal_data <- dispersal_data %>%
  dplyr::arrange(individual.local.identifier, timestamp) %>%
  dplyr::group_by(individual.local.identifier) %>%
  dplyr::mutate(
    is_terrestrial = !is.na(behavior) & behavior == "terrestrial",
    bout_change    = is_terrestrial != dplyr::lag(is_terrestrial,
                                                  default = first(is_terrestrial)),
    bout_id        = cumsum(bout_change)) %>%
  dplyr::ungroup()

bout_durations <- dispersal_data %>%
  dplyr::filter(is_terrestrial) %>%
  dplyr::group_by(individual.local.identifier, bout_id) %>%
  dplyr::summarise(
    bout_duration_h = as.numeric(difftime(max(timestamp), min(timestamp),
                                          units = "hours")),
    n_fixes_bout    = n(), .groups = "drop")

dispersal_data <- dispersal_data %>%
  dplyr::left_join(
    bout_durations %>%
      dplyr::select(individual.local.identifier, bout_id, bout_duration_h),
    by = c("individual.local.identifier", "bout_id")) %>%
  dplyr::mutate(
    behavior_refined = dplyr::case_when(
      behavior == "terrestrial" & bout_duration_h >= 6 ~ "overnight roosting",
      behavior == "terrestrial" & bout_duration_h <  6 ~ "short resting",
      TRUE ~ behavior))

dispersal_data %>%
  dplyr::filter(is_terrestrial) %>%
  dplyr::count(behavior_refined) %>%
  print()

################################################################################
#' ### STEP 3 : covariate extraction and proximity to human infrastructure
#'
#' **Proximity definition** : A GPS point is "near humans" if the maximum HFI 
#' within a buffer of radius r exceeds the 75th percentile of the HFI alpine 
#' landscape. This avoids flagging points near isolated chalets (low surrounding
#' HFI) while correctly flagging points near village edges (high surrounding HFI).
################################################################################


crs_hfi    <- terra::crs(human_footprint)

# --- 3a : project GPS points --------------------------------------------------
pts_sf <- dispersal_data %>%
  sf::st_as_sf(coords = c("location.long", "location.lat"),
               crs = 4326, remove = FALSE) %>%
  sf::st_transform(crs_hfi)

pts_vect <- terra::vect(pts_sf)

# --- 3b : point-level covariates (no buffer, fast) ---------------------------
dispersal_data$hfi             <- terra::extract(human_footprint,         pts_vect)[, 2]
dispersal_data$dist_settlement <- terra::extract(distance_settlement, pts_vect)[, 2]
dispersal_data$dens_settlement <- terra::extract(density_settlement,  pts_vect)[, 2]
dispersal_data$dens_pop_km2    <- terra::extract(density_pop_km2,     pts_vect)[, 2]

# --- 3c : maximum HFI within buffers -----------------------------------------
#
# exact_extract with fun = "max" returns the maximum pixel value within each
# buffer polygon. A high max_HFI means there is at least one anthropised
# patch within radius r of the GPS point — regardless of how many isolated
# buildings are nearby.

buffers_100m <- sf::st_buffer(pts_sf, dist = 100)
buffers_500m <- sf::st_buffer(pts_sf, dist = 500)
buffers_1km  <- sf::st_buffer(pts_sf, dist = 1000)

message("Extracting max HFI in 100m buffers...")
dispersal_data$hfi_max_100m <- exactextractr::exact_extract(
  human_footprint, buffers_100m, fun = "max", progress = FALSE)

message("Extracting max HFI in 500m buffers...")
dispersal_data$hfi_max_500m <- exactextractr::exact_extract(
  human_footprint, buffers_500m, fun = "max", progress = FALSE)

message("Extracting max HFI in 1km buffers...")
dispersal_data$hfi_max_1km <- exactextractr::exact_extract(
  human_footprint, buffers_1km, fun = "max", progress = FALSE)

# --- 3d : define proximity threshold ------------------------------------------
#
# Reference : 75th percentile of HFI over the entire Alpine raster.
# Using the landscape distribution (not GPS points) avoids circularity —
# eagles already avoid humans, so their HFI distribution is not representative
# of what "high HFI" means in the Alps.

hfi_q75 <- terra::global(human_footprint, fun = quantile,
                         probs = 0.75, na.rm = TRUE)[[1]]
message(sprintf("HFI landscape Q75 = %.4f (used as near-human threshold)", hfi_q75))

dispersal_data <- dispersal_data %>%
  dplyr::mutate(
    near_human_100m = !is.na(hfi_max_100m) & hfi_max_100m >= hfi_q75,
    near_human_500m = !is.na(hfi_max_500m) & hfi_max_500m >= hfi_q75,
    near_human_1km  = !is.na(hfi_max_1km)  & hfi_max_1km  >= hfi_q75
  )

# --- 3e : inspect retention rates ---------------------------------------------
retention <- dispersal_data %>%
  dplyr::summarise(
    n_total  = n(),
    n_100m   = sum(near_human_100m, na.rm = TRUE),
    n_500m   = sum(near_human_500m, na.rm = TRUE),
    n_1km    = sum(near_human_1km,  na.rm = TRUE),
    pct_100m = round(100 * n_100m / n_total, 1),
    pct_500m = round(100 * n_500m / n_total, 1),
    pct_1km  = round(100 * n_1km  / n_total, 1))
print(retention)
# expected : pct_100m very low (<1%), pct_500m 2-5%, pct_1km 5-15%
# if all three are 0 : lower hfi_q75 to hfi_q50 and re-run

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
    .groups  = "drop")
print(sample_summary, n = 66)

# --- 3f : spatial thinning AFTER extraction -----------------------------------
#
# dispersal_data now has all covariate and near_human columns.
# Thinning filters rows -> dispersal_thinned inherits everything.

greedy_thin_idx <- function(xy, d) {
  keep <- logical(nrow(xy))
  kept <- matrix(nrow = 0, ncol = 2)
  for (i in seq_len(nrow(xy))) {
    if (nrow(kept) == 0 ||
        all(sqrt((kept[,1]-xy[i,1])^2 + (kept[,2]-xy[i,2])^2) >= d)) {
      keep[i] <- TRUE
      kept    <- rbind(kept, xy[i,])
    }
  }
  keep
}

crs_proj <- terra::crs(distance_settlement)

dispersal_sf_full <- dispersal_data %>%
  dplyr::filter(!is.na(location.long), !is.na(location.lat)) %>%
  sf::st_as_sf(coords = c("location.long", "location.lat"),
               crs = 4326, remove = FALSE) %>%
  sf::st_transform(crs_proj)

dispersal_sf_full$keep_thinned <- FALSE
for (id in unique(dispersal_sf_full$individual.local.identifier)) {
  ix <- which(dispersal_sf_full$individual.local.identifier == id)
  ix <- ix[order(dispersal_sf_full$timestamp[ix])]
  xy <- sf::st_coordinates(dispersal_sf_full)[ix, , drop = FALSE]
  dispersal_sf_full$keep_thinned[ix] <- greedy_thin_idx(xy, 100)
}

dispersal_thinned <- dispersal_sf_full %>%
  dplyr::filter(keep_thinned) %>%
  sf::st_drop_geometry()

message(sprintf("Thinning: %d -> %d points (%.1f%% retained)",
                nrow(dispersal_sf_full), nrow(dispersal_thinned),
                100 * nrow(dispersal_thinned) / nrow(dispersal_sf_full)))

saveRDS(dispersal_data,    "C:/Users/lfaure7/Documents/git/chapter-2/preparation donnee aigle ssf/resultats-intermediaire/behavioral_classification_near_humans/dispersal_data_full_Ov.rds")
saveRDS(dispersal_thinned, "C:/Users/lfaure7/Documents/git/chapter-2/preparation donnee aigle ssf/resultats-intermediaire/behavioral_classification_near_humans/dispersal_thinned_Ov.rds")


################################################################################
#' ### STEP 4 : behavioral composition near human infrastructure
################################################################################

# behavior_colors must match behavior_labels exactly
behavior_colors <- c(
  "overnight roosting"                      = "#4e6b8c",
  "short resting"                           = "#91b4d4",
  "low soaring"                             = "#f0c04a",
  "fast flight at low elevation"            = "#e07b39",
  "fast commuting flight at high elevation" = "#c0392b"
)

# Remove final parenthetical suffix from individual names in plot labels
clean_individual_name <- function(x) {
  gsub("\\s*\\([^)]*\\)\\s*$", "", x)
}

behavior_near <- dispersal_thinned %>%
  dplyr::filter(!is.na(behavior_refined)) %>%
  tidyr::pivot_longer(cols      = c(near_human_100m, near_human_500m, near_human_1km),
                      names_to  = "threshold",
                      values_to = "near_human") %>%
  dplyr::filter(near_human) %>%
  dplyr::mutate(threshold = factor(threshold,
                                   levels = c("near_human_100m",
                                              "near_human_500m",
                                              "near_human_1km"),
                                   labels = c("≤ 100 m", "≤ 500 m", "≤ 1 km")))

behavior_composition <- behavior_near %>%
  dplyr::group_by(threshold, behavior_refined) %>%
  dplyr::summarise(n = n(), .groups = "drop") %>%
  dplyr::group_by(threshold) %>%
  dplyr::mutate(pct = 100 * n / sum(n))

p_composition <- ggplot(behavior_composition,
                        aes(x = threshold, y = pct, fill = behavior_refined)) +
  geom_col(position = "stack") +
  scale_fill_manual(values = behavior_colors) +
  labs(x     = "Buffer radius (max HFI threshold)",
       y     = "% of GPS fixes",
       fill  = "Behavior",
       title = "Behavioral composition near anthropised areas") +
  theme_minimal(base_size = 11)
print(p_composition)


# compile les résultats par individus
behavior_individual <- behavior_near %>%
  dplyr::group_by(individual.local.identifier, threshold, behavior_refined) %>%
  dplyr::summarise(n = n(), .groups = "drop") %>%
  dplyr::group_by(individual.local.identifier, threshold) %>%
  dplyr::mutate(pct = 100 * n / sum(n))

# Keep only GPS fixes near anthropised areas at ≤ 100 m
behavior_near_100m <- dispersal_thinned %>%
  dplyr::filter(
    !is.na(behavior_refined),
    near_human_100m %in% TRUE
  )

# Compute behavioral composition per individual
behavior_individual_100m <- behavior_near_100m %>%
  dplyr::group_by(individual.local.identifier, behavior_refined) %>%
  dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
  dplyr::group_by(individual.local.identifier) %>%
  dplyr::mutate(pct = 100 * n / sum(n)) %>%
  dplyr::ungroup()

# Keep only individuals with at least 30 retained GPS fixes at ≤ 100 m
individuals_min30_100m <- behavior_near_100m %>%
  dplyr::count(individual.local.identifier, name = "n_retained") %>%
  dplyr::filter(n_retained >= 30) %>%
  dplyr::pull(individual.local.identifier)

behavior_individual_100m <- behavior_individual_100m %>%
  dplyr::filter(individual.local.identifier %in% individuals_min30_100m)

# Order individuals by proportion of overnight roosting
individual_order <- behavior_individual_100m %>%
  dplyr::group_by(individual.local.identifier) %>%
  dplyr::summarise(
    pct_overnight_roosting = sum(
      pct[behavior_refined == "overnight roosting"],
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  dplyr::arrange(pct_overnight_roosting) %>%
  dplyr::pull(individual.local.identifier)

behavior_individual_100m <- behavior_individual_100m %>%
  dplyr::mutate(
    individual.local.identifier = factor(
      individual.local.identifier,
      levels = individual_order
    )
  )

# Plot for 100m 
p_individual_100m <- ggplot(
  behavior_individual_100m,
  aes(
    x = individual.local.identifier,
    y = pct,
    fill = behavior_refined
  )
) +
  geom_col(position = "stack") +
  coord_flip() +
  scale_x_discrete(labels = clean_individual_name) +
  scale_fill_manual(
    values = behavior_colors,
    drop = FALSE
  ) +
  labs(
    x = "Individual",
    y = "% of GPS fixes near anthropised areas",
    fill = "Behavior",
    title = "Individual variation in behavioral traits expressed near human-dominated areas",
    subtitle = "GPS fixes within ≤ 100 m of artificial areas; individuals ordered by overnight roosting"
  ) +
  theme_minimal(base_size = 11)

print(p_individual_100m)

# export near-human points (≤ 100m) as shapefile for QGIS inspection
dispersal_thinned %>%
  dplyr::filter(near_human_100m) %>%
  sf::st_as_sf(coords = c("location.long", "location.lat"),
               crs    = 4326,
               remove = FALSE) %>%
  sf::st_transform(crs_hfi) %>%
  dplyr::select(individual.local.identifier, timestamp,
                behavior_refined, hfi, hfi_max_100m,
                dist_settlement, dens_pop_km2,
                location.long, location.lat) %>%
  sf::st_write("C:/Users/lfaure7/Documents/git/chapter-2/preparation donnee aigle ssf/resultats-intermediaire/behavioral_classification_near_humans/comportement proche humains shape/behavior_close_humans.shp",
               delete_dsn = TRUE)


# Plot for 500m 
# Keep only GPS fixes near anthropised areas at ≤ 500 m
behavior_near_500m <- dispersal_thinned %>%
  dplyr::filter(
    !is.na(behavior_refined),
    near_human_500m %in% TRUE
  )

# Compute behavioral composition per individual
behavior_individual_500m <- behavior_near_500m %>%
  dplyr::group_by(individual.local.identifier, behavior_refined) %>%
  dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
  dplyr::group_by(individual.local.identifier) %>%
  dplyr::mutate(pct = 100 * n / sum(n)) %>%
  dplyr::ungroup()

# Keep only individuals with at least 30 retained GPS fixes at ≤ 500 m
individuals_min30_500m <- behavior_near_500m %>%
  dplyr::count(individual.local.identifier, name = "n_retained") %>%
  dplyr::filter(n_retained >= 30) %>%
  dplyr::pull(individual.local.identifier)

behavior_individual_500m <- behavior_individual_500m %>%
  dplyr::filter(individual.local.identifier %in% individuals_min30_500m)

# Order individuals by proportion of overnight roosting
individual_order_500m <- behavior_individual_500m %>%
  dplyr::group_by(individual.local.identifier) %>%
  dplyr::summarise(
    pct_overnight_roosting = sum(
      pct[behavior_refined == "overnight roosting"],
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  dplyr::arrange(pct_overnight_roosting) %>%
  dplyr::pull(individual.local.identifier)

behavior_individual_500m <- behavior_individual_500m %>%
  dplyr::mutate(
    individual.local.identifier = factor(
      individual.local.identifier,
      levels = individual_order_500m
    )
  )

# Plot for 500m
p_individual_500m <- ggplot(
  behavior_individual_500m,
  aes(
    x = individual.local.identifier,
    y = pct,
    fill = behavior_refined
  )
) +
  geom_col(position = "stack") +
  coord_flip() +
  scale_x_discrete(labels = clean_individual_name) +
  scale_fill_manual(
    values = behavior_colors,
    drop = FALSE
  ) +
  labs(
    x = "Individual",
    y = "% of GPS fixes near anthropised areas",
    fill = "Behavior",
    title = "Individual variation in behavioral traits expressed near human dominated areas",
    subtitle = "GPS fixes within ≤ 500 m of anthropised areas; individuals ordered by overnight roosting"
  ) +
  theme_minimal(base_size = 11)

print(p_individual_500m)


################################################################################
#' ### ANNEX 1 : correlation coefficent for anthropogenic features
################################################################################

# load raster data
distance_settlement <- rast("C:/Users/lfaure7/Desktop/COUCHES QGIS/settlements/Overture/distance_to_building_m_3035.tif")
density_settlement <- rast("C:/Users/lfaure7/Desktop/COUCHES QGIS/settlements/Overture/building_density_count_per_km2_100m_grass.tif")
proportion_building <- rast("C:/Users/lfaure7/Desktop/COUCHES QGIS/settlements/Overture/built_proportion_1km_window_100m_grass.tif")
density_pop_km2 <- rast("C:/Users/lfaure7/Desktop/COUCHES QGIS/density-resident-pop-km2/density_pop_km2_m/density_pop_km2_meters_crs.tif")

# Remove negative distances
distance_settlement <- ifel(distance_settlement < 0, NA, distance_settlement)

# Selection of reference grid at 100m resolution and reprojection of distance to
# settlement raster
ref <- proportion_building

# aggregation of distance to settlement and population density at 100m resolution
density_pop_km2_r <- project(
  density_pop_km2,
  ref,
  method = "bilinear"
)

distance_settlement_r <- project(
  distance_settlement,
  ref,
  method = "bilinear"
)

proportion_building_r <- ref

names(distance_settlement_r)   <- "dist_settlement_m"
names(density_settlement_r)    <- "dens_building_km2"
names(density_pop_km2_r)       <- "dens_pop_km2"
names(proportion_building_r)   <- "built_prop_1km"

# Correlation matrix of anthropogenic features of the landscape
corr_stack <- c(
  distance_settlement_r,
  density_settlement_r,
  density_pop_km2_r,
  proportion_building_r
)

names(corr_stack) <- c(
  "dist_settlement_m",
  "dens_building_km2",
  "dens_pop_km2",
  "built_prop_1km"
)

set.seed(123)

samp_corr <- spatSample(
  corr_stack,
  size = 100000,
  method = "random",
  na.rm = TRUE,
  as.df = TRUE
)

cor_spearman <- cor(
  samp_corr,
  method = "spearman",
  use = "complete.obs"
)

print(cor_spearman)

# Results ######################################################################
#'                        dist_settlement_m dens_building_km2 dens_pop_km2 built_prop_1km
#' dist_settlement_m         1.0000000        -0.8720563   -0.7197879     -0.8841070
#' dens_building_km2        -0.8720563         1.0000000    0.8444852      0.9894498
#' dens_pop_km2             -0.7197879         0.8444852    1.0000000      0.8470960
#' built_prop_1km           -0.8841070         0.9894498    0.8470960      1.0000000
#' #############################################################################



################################################################################
#' ### ANNEX 2 : prepare the human footprint index
################################################################################

# log
dens_building_log <- log1p(density_settlement_r)
dens_pop_log      <- log1p(density_pop_km2_r)
built_prop_trans  <- proportion_building_r

names(dens_building_log) <- "log1p_dens_building_km2"
names(dens_pop_log)      <- "log1p_dens_pop_km2"
names(built_prop_trans)  <- "built_prop_1km"

# normalisation
normalize_robust <- function(r, name_out, p_low = 0.01, p_high = 0.99) {
  
  set.seed(123)
  
  vals <- spatSample(
    r,
    size = min(1000000, ncell(r)),
    method = "random",
    na.rm = TRUE,
    as.df = TRUE
  )[[1]]
  
  q_low  <- as.numeric(quantile(vals, probs = p_low,  na.rm = TRUE))
  q_high <- as.numeric(quantile(vals, probs = p_high, na.rm = TRUE))
  
  message(name_out, " : q", p_low, " = ", q_low,
          " ; q", p_high, " = ", q_high)
  
  if (!is.finite(q_low) || !is.finite(q_high) || q_high <= q_low) {
    stop(paste("Normalisation impossible pour", name_out))
  }
  
  out <- (r - q_low) / (q_high - q_low)
  out <- ifel(out < 0, 0, ifel(out > 1, 1, out))
  
  names(out) <- name_out
  
  return(out)
}

dens_building_norm <- normalize_robust(
  dens_building_log,
  name_out = "dens_building_norm"
)

dens_pop_norm <- normalize_robust(
  dens_pop_log,
  name_out = "dens_pop_norm"
)

built_prop_norm <- normalize_robust(
  built_prop_trans,
  name_out = "built_prop_norm"
)

# combine anthropogenic variables 
hfi <- (dens_building_norm + dens_pop_norm + built_prop_norm) / 3
names(hfi) <- "hfi"

# crop to the Alps extend 
alpes_v <- vect(
  st_read(
    "C:/Users/lfaure7/Desktop/COUCHES QGIS/alpine_area/alpine_area_corrected_geometry.shp",
    quiet = TRUE
  )
)

alpes_v <- project(alpes_v, crs(hfi))

hfi_masked <- mask(hfi, alpes_v)
names(hfi_masked) <- "hfi"
plot(hfi_masked)
