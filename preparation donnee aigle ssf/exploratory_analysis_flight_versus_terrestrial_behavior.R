#' ---
#' title: "Preparing Golden Eagle dataset for ssf"
#' author: "Louise Faure"
#' date: 02.06.2026
#' details: (i) select all the golden eagles gps values for the first 15 weeks of the dispersal phase, 
#' (ii) classify behavior using the embc package
#' (iii) identify the landscape features associated to the expression of each behavior
#' (iv) generate random point within a 10 000m radius for both terrestrial and flight behavior and compared to the "used" terrestrial and 
#' flight point for one individual (Adamello)    
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
library(circular)       # statistiques circulaires (idem)
#library(fitdistrplus)   # ajustement de distributions (idem)
#library(units)          # gestion des unités (idem)
library(exactextractr)  # extraction fractionnelle dans un buffer

library(ggplot2)
library(ggridges)     # ridgeline plots (distribution shape per individual)
library(ggh4x)


#' ### STEP 0 : load the data --------------------------------------------------

output_dir <- "C:/Users/lfaure7/Documents/git/chapter-2/preparation donnee aigle ssf/resultats-intermediaire"

# Géoïde and digital elevation model
geo <- terra::rast("C:/Users/lfaure7/OneDrive/MEMOIRE M2/eagle_projet/data/geoide/us_nga_egm96_15.tif")
dem <- terra::rast("C:/Users/lfaure7/OneDrive/MEMOIRE M2/eagle_projet/data/pretraitements/Region-Alpes-Dem/region-alpes-dem.tif")

# raster layers
covar_1_elevation <- terra::rast("C:/Users/lfaure7/OneDrive/MEMOIRE M2/eagle_projet/data/pretraitements/Region-Alpes-Dem/region-alpes-dem.tif")
covar_2_distance_ridgeline <- terra::rast("C:/Users/lfaure7/OneDrive/MEMOIRE M2/donnees/raster/topography/distance_to_ridge_line_complete_version.tif")
covar_3_TRI <- terra::rast("C:/Users/lfaure7/OneDrive/MEMOIRE M2/eagle_projet/data/pretraitements/DEM_25m/TRI/TRI.tif")
covar_4_slope <- terra::rast("C:/Users/lfaure7/OneDrive/MEMOIRE M2/eagle_projet/data/pretraitements/DEM_25m/slope/slope_25.tif")
covar_5_TPI <- terra::rast("C:/Users/lfaure7/OneDrive/MEMOIRE M2/eagle_projet/data/pretraitements/DEM_25m/TPI/Topographic Position Index.tif")
covar_6_distance_settlement <- terra::rast("C:/Users/lfaure7/Desktop/COUCHES QGIS/settlements/distance to settlements 30 m cropped/dist_to_settlements_30m_croped.tif")
covar_7_density_settlement <- terra::rast("C:/Users/lfaure7/Desktop/COUCHES QGIS/settlements/density_rad_495m.tif")
covar_8_density_pop_hectar <- terra::rast("C:/Users/lfaure7/Desktop/COUCHES QGIS/density-resident-pop-hectar/density_pop_hectar.tif")
covar_9_density_pop_km2 <- terra::rast("C:/Users/lfaure7/Desktop/COUCHES QGIS/density-resident-pop-km2/density_pop_km2.tif")
covar_10_landcover <- terra::rast("C:/Users/lfaure7/OneDrive/THESE/Conférences/Sempach workshop/Donnees/Landscape layers/CLC_longlat_10m.tif")


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
behavioural_classification <- data.matrix(dispersal_data[, c("ground.speed", "height_above_ground")])

# call embc
embc<- embc(behavioural_classification)

# investigate the bc
X11(); sctr (embc)

embc_smoothed <- smth(embc, dlta = 0.7)
X11();sctr(embc_smoothed)

# link the cluster labels (1: low speed, low height, 2: low speed, height height, 3: heigh speed, low height, 4: height speed, heigh height) to the original data, e.i, dispersal_data
dispersal_data$behavior_cluster <- embc_smoothed@A

# check that cluster 1 correspond to terrestrial behaviors
dispersal_data %>%
  group_by(behavior_cluster) %>%
  summarise(speed_med  = median(ground.speed,        na.rm = TRUE),
            height_med = median(height_above_ground, na.rm = TRUE),
            n          = n()) # terrestrial behaviour should have a speed ~0, and a height ~0 





#' ### STEP 3 : Landscape characteristics for terrestrial and flying behaviours 


# Projected CRS in metres for distance and buffer operations
crs_proj <- terra::crs(dem)


#' ---------------------------------------------------------------------------
#' Helper functions
#' ---------------------------------------------------------------------------

make_points_sf <- function(df) {
  df %>%
    dplyr::filter(!is.na(location.long), !is.na(location.lat)) %>%
    sf::st_as_sf(
      coords = c("location.long", "location.lat"),
      crs = 4326,
      remove = FALSE
    ) %>%
    sf::st_transform(crs_proj)
}

consec_dist <- function(geom) {
  n <- length(geom)
  if (n < 2) return(rep(NA_real_, n))
  c(NA_real_, as.numeric(sf::st_distance(geom[-n], geom[-1], by_element = TRUE)))
}

add_step_distance <- function(points_sf) {
  points_sf$step_m <- NA_real_
  
  for (id in unique(points_sf$individual.local.identifier)) {
    ix <- which(points_sf$individual.local.identifier == id)
    ix <- ix[order(points_sf$timestamp[ix])]
    points_sf$step_m[ix] <- consec_dist(sf::st_geometry(points_sf)[ix])
  }
  
  points_sf
}

greedy_thin_idx <- function(xy, d) {
  keep <- logical(nrow(xy))
  kept <- matrix(nrow = 0, ncol = 2)
  
  for (i in seq_len(nrow(xy))) {
    if (
      nrow(kept) == 0 ||
      all(sqrt((kept[, 1] - xy[i, 1])^2 + (kept[, 2] - xy[i, 2])^2) >= d)
    ) {
      keep[i] <- TRUE
      kept <- rbind(kept, xy[i, ])
    }
  }
  
  keep
}

thin_by_distance <- function(points_sf, d) {
  points_sf$keep <- FALSE
  
  for (id in unique(points_sf$individual.local.identifier)) {
    ix <- which(points_sf$individual.local.identifier == id)
    ix <- ix[order(points_sf$timestamp[ix])]
    xy <- sf::st_coordinates(points_sf)[ix, , drop = FALSE]
    points_sf$keep[ix] <- greedy_thin_idx(xy, d)
  }
  
  points_sf %>%
    dplyr::filter(keep)
}

extract_into <- function(buffers_proj, r, categorical = FALSE, fun = "mean") {
  b <- sf::st_transform(buffers_proj, terra::crs(r))
  
  if (categorical) {
    exactextractr::exact_extract(r, b, "frac", progress = FALSE)
  } else {
    exactextractr::exact_extract(r, b, fun, progress = FALSE)
  }
}

extract_landscape <- function(points_sf, buf_r = 100) {
  if (nrow(points_sf) == 0) {
    stop("No points supplied to extract_landscape().")
  }
  
  buffers <- sf::st_buffer(points_sf, dist = buf_r)
  
  points_sf$elevation       <- extract_into(buffers, covar_1_elevation)
  points_sf$dist_ridgeline  <- extract_into(buffers, covar_2_distance_ridgeline)
  points_sf$tri             <- extract_into(buffers, covar_3_TRI)
  points_sf$slope           <- extract_into(buffers, covar_4_slope)
  points_sf$tpi             <- extract_into(buffers, covar_5_TPI)
  points_sf$dist_settlement <- extract_into(buffers, covar_6_distance_settlement)
  points_sf$dens_settlement <- extract_into(buffers, covar_7_density_settlement)
  points_sf$pop_dens_hectar        <- extract_into(buffers, covar_8_density_pop_hectar)
  points_sf$pop_dens_km2        <- extract_into(buffers, covar_9_density_pop_km2)
  
  lc_frac <- extract_into(
    buffers,
    covar_10_landcover,
    categorical = TRUE
  )
  
  points_sf %>%
    sf::st_drop_geometry() %>%
    dplyr::bind_cols(as.data.frame(lc_frac))
}

#' ---------------------------------------------------------------------------
#' STEP 3a : terrestrial points, thinned spatially
#' ---------------------------------------------------------------------------

terr <- dispersal_data %>%
  dplyr::filter(behavior_cluster == 1) %>%
  dplyr::mutate(
    behavior_type = "terrestrial",
    behavior_subtype = "low_speed_low_height"
  )

terr_sf <- make_points_sf(terr) %>%
  add_step_distance()

med_step <- terr_sf %>%
  sf::st_drop_geometry() %>%
  dplyr::group_by(individual.local.identifier) %>%
  dplyr::summarise(
    median_step_m = median(step_m, na.rm = TRUE),
    q25 = quantile(step_m, 0.25, na.rm = TRUE),
    q75 = quantile(step_m, 0.75, na.rm = TRUE),
    n = dplyr::n(),
    .groups = "drop"
  )

print(med_step)

d_thin <- median(med_step$median_step_m, na.rm = TRUE)

if (is.na(d_thin) || !is.finite(d_thin)) {
  stop("d_thin is NA or infinite. Check terrestrial step distances.")
}

terr_thin <- thin_by_distance(terr_sf, d_thin)

message(sprintf(
  "Terrestrial thinning at d = %.0f m: %d -> %d points",
  d_thin, nrow(terr_sf), nrow(terr_thin)
))

#' ---------------------------------------------------------------------------
#' STEP 3b : flying points, clusters 3 and 4, not thinned
#' ---------------------------------------------------------------------------

fly <- dispersal_data %>%
  dplyr::filter(behavior_cluster %in% c(3, 4)) %>%
  dplyr::mutate(
    behavior_type = "flying",
    behavior_subtype = dplyr::case_when(
      behavior_cluster == 3 ~ "high_speed_low_height",
      behavior_cluster == 4 ~ "high_speed_high_height",
      TRUE ~ NA_character_
    )
  )

fly_sf <- make_points_sf(fly)


# Add consecutive distances, only for inspection.
# The thinning itself uses the same greedy distance rule as terrestrial points.
fly_sf <- add_step_distance(fly_sf)

# Number of flying points per individual before thinning
fly_n_before <- fly_sf %>%
  sf::st_drop_geometry() %>%
  dplyr::count(
    individual.local.identifier,
    name = "n_flying_before_thinning"
  )

print(fly_n_before)

# Apply the same thinning distance as for terrestrial points
fly_thin <- thin_by_distance(fly_sf, d_thin)

# Number of flying points per individual after thinning
fly_n_after <- fly_thin %>%
  sf::st_drop_geometry() %>%
  dplyr::count(
    individual.local.identifier,
    name = "n_flying_after_thinning"
  )

# Summary table: before vs after
fly_n_by_ind <- fly_n_before %>%
  dplyr::full_join(
    fly_n_after,
    by = "individual.local.identifier"
  ) %>%
  dplyr::mutate(
    n_flying_before_thinning = tidyr::replace_na(n_flying_before_thinning, 0L),
    n_flying_after_thinning  = tidyr::replace_na(n_flying_after_thinning, 0L),
    n_removed = n_flying_before_thinning - n_flying_after_thinning,
    pct_retained = 100 * n_flying_after_thinning / n_flying_before_thinning
  ) %>%
  dplyr::arrange(individual.local.identifier)

print(fly_n_by_ind)

message(sprintf(
  "Flying thinning at d = %.0f m: %d -> %d points",
  d_thin, nrow(fly_sf), nrow(fly_thin)
))


#' ---------------------------------------------------------------------------
#' STEP 3c : extract landscape characteristics within buffers
#' ---------------------------------------------------------------------------

buf_r <- 100

terr_covars <- extract_landscape(terr_thin, buf_r = buf_r)
fly_covars  <- extract_landscape(fly_thin,    buf_r = buf_r)

# Combined table for all analyses
covar_all <- dplyr::bind_rows(terr_covars, fly_covars)

# Landcover fraction columns that are absent in one behaviour after bind_rows()
# should be treated as 0 fractional cover, not as missing values.
frac_vars <- names(covar_all)[startsWith(names(covar_all), "frac_")]

covar_all <- covar_all %>%
  dplyr::mutate(
    dplyr::across(dplyr::all_of(frac_vars), ~ tidyr::replace_na(.x, 0))
  )

# Keep separate objects for compatibility with your previous code
covar_tbl     <- covar_all %>% dplyr::filter(behavior_type == "terrestrial")
fly_covar_tbl <- covar_all %>% dplyr::filter(behavior_type == "flying")

saveRDS(covar_tbl,     file.path(output_dir, "terrestrial_covariates.rds"))
saveRDS(fly_covar_tbl, file.path(output_dir, "flying_covariates_clusters_3_4.rds"))
saveRDS(covar_all,     file.path(output_dir, "terrestrial_and_flying_covariates.rds"))

message(sprintf(
  "Landscape extraction complete: %d terrestrial + %d flying = %d retained points",
  nrow(covar_tbl), nrow(fly_covar_tbl), nrow(covar_all)
))

#' #############################################################################
#' STEP 4 : Human-related variables, terrestrial vs flying
#' #############################################################################

human_vars <- c(
  "dist_settlement",
  "dens_settlement",
  "pop_dens_hectar",
  "pop_dens_km2"
)

missing_human_vars <- setdiff(human_vars, names(covar_all))

if (length(missing_human_vars) > 0) {
  stop(
    "Missing human variables in covar_all: ",
    paste(missing_human_vars, collapse = ", ")
  )
}

human_var_labels <- c(
  dist_settlement = "Distance to settlement",
  dens_settlement = "Settlement density",
  pop_dens_hectar = "Micro local population density",
  pop_dens_km2 = "Local population density"
)

# Function to clean individual names:
# "Almen18 (tag)" becomes "Almen18"
clean_individual_name <- function(x) {
  x %>%
    as.character() %>%
    stringr::str_replace("\\s*\\([^\\)]*\\)\\s*$", "") %>%
    stringr::str_trim()
}

long_human <- covar_all %>%
  dplyr::mutate(
    individual.clean = clean_individual_name(individual.local.identifier)
  ) %>%
  dplyr::filter(
    !individual.clean %in% c("Almen18", "Langgries21")
  ) %>%
  dplyr::select(
    behavior_type,
    individual.local.identifier,
    individual.clean,
    dplyr::all_of(human_vars)
  ) %>%
  tidyr::pivot_longer(
    cols = dplyr::all_of(human_vars),
    names_to = "variable",
    values_to = "value"
  ) %>%
  dplyr::filter(!is.na(value)) %>%
  dplyr::mutate(
    behavior_type = factor(
      behavior_type,
      levels = c("terrestrial", "flying")
    ),
    variable_label = factor(
      unname(human_var_labels[variable]),
      levels = unname(human_var_labels)
    ),
    individual.clean = factor(individual.clean)
  )

make_one_human_ridge_plot <- function(data, var_name) {
  
  if (!var_name %in% unique(data$variable)) {
    stop("Variable not found in long_human: ", var_name)
  }
  
  plot_data <- data %>%
    dplyr::filter(variable == var_name) %>%
    dplyr::mutate(
      individual.clean = factor(individual.clean)
    )
  
  plot_title <- human_var_labels[[var_name]]
  
  p <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = value,
      y = individual.clean,
      fill = behavior_type
    )
  ) +
    ggridges::geom_density_ridges(
      scale = 1.2,
      alpha = 0.6,
      rel_min_height = 0.01,
      quantile_lines = TRUE,
      quantiles = 2
    ) +
    ggplot2::guides(fill = "none") +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::labs(
      y = "individual",
      x = "covariate value within 100 m buffer",
      title = plot_title
    )
  
  if (requireNamespace("ggh4x", quietly = TRUE)) {
    p <- p +
      ggh4x::facet_grid2(
        rows = ggplot2::vars(variable_label),
        cols = ggplot2::vars(behavior_type),
        scales = "free_x",
        independent = "x",
        labeller = ggplot2::labeller(.default = ggplot2::label_value)
      )
  } else {
    warning(
      "Package 'ggh4x' is not installed. Using facet_wrap fallback. ",
      "Install ggh4x for exact column layout with independent x-scales."
    )
    
    p <- p +
      ggplot2::facet_wrap(
        ~ behavior_type,
        scales = "free_x",
        ncol = 2,
        labeller = ggplot2::labeller(.default = ggplot2::label_value)
      )
  }
  
  p
}

p_dist_settlement <- make_one_human_ridge_plot(
  data = long_human,
  var_name = "dist_settlement"
)

p_dens_settlement <- make_one_human_ridge_plot(
  data = long_human,
  var_name = "dens_settlement"
)

p_pop_dens_hectar <- make_one_human_ridge_plot(
  data = long_human,
  var_name = "pop_dens_hectar"
)

p_pop_dens_km2 <- make_one_human_ridge_plot(
  data = long_human,
  var_name = "pop_dens_km2"
)

p_dist_settlement
p_dens_settlement
p_pop_dens_hectar
p_pop_dens_km2


#' ---------------------------------------------------------------------------
#' Save STEP 4 plots
#' ---------------------------------------------------------------------------


plot_dir <- "C:/Users/lfaure7/Documents/git/chapter-2/preparation donnee aigle ssf/dossier de plots/resultats preliminaire"

if (!dir.exists(plot_dir)) {
  dir.create(plot_dir, recursive = TRUE)
}

ggplot2::ggsave(
  filename = file.path(plot_dir, "human_dist_settlement_terrestrial_vs_flying.png"),
  plot = p_dist_settlement,
  width = 10,
  height = 8,
  dpi = 300
)

ggplot2::ggsave(
  filename = file.path(plot_dir, "human_dens_settlement_terrestrial_vs_flying.png"),
  plot = p_dens_settlement,
  width = 10,
  height = 8,
  dpi = 300
)

ggplot2::ggsave(
  filename = file.path(plot_dir, "human_pop_dens_hectar_terrestrial_vs_flying.png"),
  plot = p_pop_dens_hectar,
  width = 10,
  height = 8,
  dpi = 300
)

ggplot2::ggsave(
  filename = file.path(plot_dir, "human_pop_dens_km2_terrestrial_vs_flying.png"),
  plot = p_pop_dens_km2,
  width = 10,
  height = 8,
  dpi = 300
)


#' #############################################################################
#' STEP 5 : Correlation matrix for landscape variables
#' variables 

covar_all <- readRDS("resultats-intermediaire/terrestrial_and_flying_covariates.rds")

#' #############################################################################

# These are the continuous landscape variables.
# height_above_ground is deliberately excluded because it describes the eagle's
# position above the ground, not the landscape itself.
landscape_cont_vars <- c(
  "elevation",
  "dist_ridgeline",
  "tri",
  "slope",
  "tpi",
  "dist_settlement",
  "dens_settlement",
  "pop_dens_hectar", 
  "pop_dens_km2"
)

missing_landscape_vars <- setdiff(landscape_cont_vars, names(covar_all))
if (length(missing_landscape_vars) > 0) {
  warning(
    "Missing landscape variables: ",
    paste(missing_landscape_vars, collapse = ", ")
  )
}

landscape_cont_vars <- intersect(landscape_cont_vars, names(covar_all))

corr_dat <- covar_all %>%
  dplyr::select(dplyr::all_of(landscape_cont_vars)) %>%
  dplyr::mutate(dplyr::across(dplyr::everything(), as.numeric))

M_spearman <- cor(
  corr_dat,
  use = "pairwise.complete.obs",
  method = "spearman"
)

print(round(M_spearman, 2))

corr_long <- as.data.frame(as.table(M_spearman)) %>%
  dplyr::as_tibble() %>%
  dplyr::rename(
    var_x = Var1,
    var_y = Var2,
    rho = Freq
  )

p_corr <- ggplot2::ggplot(
  corr_long,
  ggplot2::aes(x = var_x, y = var_y, fill = rho)
) +
  ggplot2::geom_tile() +
  ggplot2::geom_text(
    ggplot2::aes(label = round(rho, 2)),
    size = 3
  ) +
  ggplot2::scale_fill_gradient2(
    limits = c(-1, 1),
    midpoint = 0,
    na.value = "grey90"
  ) +
  ggplot2::coord_equal() +
  ggplot2::theme_minimal(base_size = 11) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
  ) +
  ggplot2::labs(
    x = NULL,
    y = NULL,
    fill = "Spearman rho",
    title = "Correlation matrix of continuous landscape variables"
  )

p_corr


ggplot2::ggsave(
  filename = file.path(plot_dir, "correlation_matrix_continuous_landscape_spearman.png"),
  plot = p_corr,
  width = 10,
  height = 8,
  dpi = 300
)


#' #############################################################################
#' STEP 6 : Landcover composition, terrestrial vs flying
#' #############################################################################

if (length(frac_vars) == 0) {
  stop("No landcover fraction columns found. Expected columns starting with 'frac_'.")
}

# CORINE Land Cover lookup.
# This assumes your raster values are coded as 1:44 in standard CLC order.
# The lookup also supports standard CLC codes such as 111, 112, 121, etc.
# CORINE / CLCplus Backbone lookup.
# Raster codes according to CLCplus Backbone raster product:
# 1-11 = land-cover classes
# 253 = coastal seawater buffer
# 254 = outside area
# 255 = no data

clc_lookup <- tibble::tibble(
  class_code = c(
    1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 253, 254, 255
  ),
  landcover_name = c(
    "Sealed",
    "Woody - needle leaved trees",
    "Woody - broadleaved deciduous trees",
    "Woody - broadleaved evergreen trees",
    "Low-growing woody plants",
    "Permanent herbaceous",
    "Periodically herbaceous",
    "Lichens and mosses",
    "Non- and sparsely-vegetated",
    "Water",
    "Snow and ice",
    "Coastal seawater buffer",
    "Outside area",
    "No data"
  ),
  landcover_level1 = c(
    "Artificial surfaces",             # 1 Sealed
    "Forest",                          # 2 Needle-leaved trees
    "Forest",                          # 3 Broadleaved deciduous trees
    "Forest",                          # 4 Broadleaved evergreen trees
    "Shrubland / low woody vegetation",# 5 Low-growing woody plants
    "Herbaceous vegetation",           # 6 Permanent herbaceous
    "Herbaceous vegetation",           # 7 Periodically herbaceous
    "Sparse / non-vascular vegetation",# 8 Lichens and mosses
    "Bare or sparsely vegetated",      # 9 Non- and sparsely-vegetated
    "Water bodies",                    # 10 Water
    "Snow and ice",                    # 11 Snow and ice
    "Water bodies",                    # 253 Coastal seawater buffer
    "Outside area",                    # 254 Outside area
    "No data"                          # 255 No data
  )
)

lc_long <- covar_all %>%
  dplyr::select(
    behavior_type,
    individual.local.identifier,
    dplyr::all_of(frac_vars)
  ) %>%
  dplyr::group_by(behavior_type, individual.local.identifier) %>%
  dplyr::summarise(
    dplyr::across(
      dplyr::all_of(frac_vars),
      ~ mean(.x, na.rm = TRUE)
    ),
    .groups = "drop"
  ) %>%
  tidyr::pivot_longer(
    cols = dplyr::all_of(frac_vars),
    names_to = "class",
    values_to = "frac"
  ) %>%
  dplyr::mutate(
    class_code = as.integer(sub("^frac_", "", class))
  ) %>%
  dplyr::left_join(clc_lookup, by = "class_code") %>%
  dplyr::mutate(
    landcover_name = dplyr::if_else(
      is.na(landcover_name),
      class,
      landcover_name
    ),
    landcover_level1 = dplyr::if_else(
      is.na(landcover_level1),
      "Unmapped",
      landcover_level1
    ),
    behavior_type = factor(
      behavior_type,
      levels = c("terrestrial", "flying")
    )
  ) %>%
  dplyr::filter(!is.na(frac), frac > 0)



# Optional coarser plot: CLC level-1 classes, easier to read
lc_level1 <- lc_long %>%
  dplyr::group_by(
    behavior_type,
    individual.local.identifier,
    landcover_level1
  ) %>%
  dplyr::summarise(
    frac = sum(frac, na.rm = TRUE),
    .groups = "drop"
  )


# color code : 
# Color code for CLC level-1 classes
lc_level1_colors <- c(
  "Artificial surfaces" = "snow4",
  "Bare or sparsely vegetated" = "bisque2",
  "Forest" = "seagreen",
  "Herbaceous vegetation" = "darkolivegreen2",
  "Shrubland / low woody vegetation" = "lightgoldenrod1",
  "Snow and ice" = "lightskyblue1",
  "Water bodies" = "skyblue4"
)


p_lc_level1 <- ggplot2::ggplot(
  lc_level1,
  ggplot2::aes(
    x = factor(individual.local.identifier),
    y = frac,
    fill = landcover_level1
  )
) +
  ggplot2::geom_col() +
  ggplot2::coord_flip() +
  ggplot2::facet_grid(. ~ behavior_type) +
  ggplot2::scale_fill_manual(
    values = lc_level1_colors,
    na.value = "grey80"
  ) +
  ggplot2::theme_minimal(base_size = 11) +
  ggplot2::labs(
    x = "individual",
    y = "mean fractional cover in 100 m buffer",
    fill = "CLC level 1",
    title = "Landcover composition by broad CLC class"
  )

p_lc_level1

ggplot2::ggsave(
  filename = file.path(plot_dir, "landcover_per_behavior_type_per_ind.png"),
  plot = p_lc_level1,
  width = 10,
  height = 8,
  dpi = 300
)

#' test for one individual

#' #############################################################################
#' STEP 7 : Used vs available analysis for one individual
#' #############################################################################

# Individual to test
target_ind <- "Adamello20"

# Large availability buffer around terrestrial points, in metres
available_buffer_m <- 10000

# Number of random available points
n_available <- 5000

# Buffer used for covariate extraction around random available points
buf_r <- 100

# Clean individual names if tags are written in parentheses
clean_individual_name <- function(x) {
  trimws(gsub("\\s*\\([^\\)]*\\)\\s*$", "", as.character(x)))
}

# Check available individual names
individual_names <- dispersal_data %>%
  dplyr::mutate(individual.clean = clean_individual_name(individual.local.identifier)) %>%
  dplyr::distinct(individual.local.identifier, individual.clean) %>%
  dplyr::arrange(individual.clean)

print(individual_names)

# Select Adamello
adamello_ids <- individual_names %>%
  dplyr::filter(grepl(target_ind, individual.clean, ignore.case = TRUE))

print(adamello_ids)

if (nrow(adamello_ids) == 0) {
  stop("No individual matching '", target_ind, "' was found.")
}

if (nrow(adamello_ids) > 1) {
  stop(
    "Several individuals match '", target_ind, "'. ",
    "Please choose one exact individual.local.identifier."
  )
}

target_id_raw <- adamello_ids$individual.local.identifier[1]
target_id_clean <- adamello_ids$individual.clean[1]

message("Selected individual: ", target_id_clean)

#' ---------------------------------------------------------------------------
#' Separate availability areas for terrestrial and flying behaviour
#' ---------------------------------------------------------------------------

available_buffer_terr_m <- 10000
available_buffer_fly_m  <- 10000

n_available_terr <- 5000
n_available_fly  <- 5000

set.seed(123)

# Terrestrial used points for Adamello
adamello_terr_sf <- terr_thin %>%
  dplyr::mutate(
    individual.clean = clean_individual_name(individual.local.identifier)
  ) %>%
  dplyr::filter(individual.clean == target_id_clean)

# Flying used points for Adamello
adamello_fly_sf <- fly_thin %>%
  dplyr::mutate(
    individual.clean = clean_individual_name(individual.local.identifier)
  ) %>%
  dplyr::filter(individual.clean == target_id_clean)

if (nrow(adamello_terr_sf) == 0) {
  stop("No terrestrial points found for ", target_id_clean)
}

if (nrow(adamello_fly_sf) == 0) {
  stop("No flying points found for ", target_id_clean)
}

# Availability area around terrestrial points
adamello_available_terr_area <- adamello_terr_sf %>%
  sf::st_union() %>%
  sf::st_buffer(dist = available_buffer_terr_m)

# Availability area around flying points
adamello_available_fly_area <- adamello_fly_sf %>%
  sf::st_union() %>%
  sf::st_buffer(dist = available_buffer_fly_m)

# Random available terrestrial points
available_terr_geom <- sf::st_sample(
  adamello_available_terr_area,
  size = n_available_terr,
  type = "random"
)

adamello_available_terr_sf <- sf::st_sf(
  individual.local.identifier = target_id_raw,
  individual.clean = target_id_clean,
  behavior_type = "available",
  behavior_subtype = "available_terrestrial_buffer",
  availability_class = "available_terrestrial",
  geometry = available_terr_geom,
  crs = sf::st_crs(adamello_terr_sf)
)

# Random available flying points
available_fly_geom <- sf::st_sample(
  adamello_available_fly_area,
  size = n_available_fly,
  type = "random"
)

adamello_available_fly_sf <- sf::st_sf(
  individual.local.identifier = target_id_raw,
  individual.clean = target_id_clean,
  behavior_type = "available",
  behavior_subtype = "available_flying_buffer",
  availability_class = "available_flying",
  geometry = available_fly_geom,
  crs = sf::st_crs(adamello_fly_sf)
)


adamello_available_terr_covars <- extract_landscape(
  adamello_available_terr_sf,
  buf_r = buf_r
) %>%
  dplyr::mutate(
    individual.clean = target_id_clean,
    availability_class = "available_terrestrial"
  )

adamello_available_fly_covars <- extract_landscape(
  adamello_available_fly_sf,
  buf_r = buf_r
) %>%
  dplyr::mutate(
    individual.clean = target_id_clean,
    availability_class = "available_flying"
  )


adamello_used_terr <- covar_all %>%
  dplyr::mutate(
    individual.clean = clean_individual_name(individual.local.identifier)
  ) %>%
  dplyr::filter(
    individual.clean == target_id_clean,
    behavior_type == "terrestrial"
  ) %>%
  dplyr::mutate(
    availability_class = "used_terrestrial"
  )

adamello_used_fly <- covar_all %>%
  dplyr::mutate(
    individual.clean = clean_individual_name(individual.local.identifier)
  ) %>%
  dplyr::filter(
    individual.clean == target_id_clean,
    behavior_type == "flying"
  ) %>%
  dplyr::mutate(
    availability_class = "used_flying"
  )

# combine 
adamello_used_available_2scale <- dplyr::bind_rows(
  adamello_available_terr_covars,
  adamello_used_terr,
  adamello_available_fly_covars,
  adamello_used_fly
)

frac_vars_adamello_2scale <- names(adamello_used_available_2scale)[
  startsWith(names(adamello_used_available_2scale), "frac_")
]

adamello_used_available_2scale <- adamello_used_available_2scale %>%
  dplyr::mutate(
    dplyr::across(
      dplyr::all_of(frac_vars_adamello_2scale),
      ~ tidyr::replace_na(.x, 0)
    ),
    availability_class = factor(
      availability_class,
      levels = c(
        "available_terrestrial",
        "used_terrestrial",
        "available_flying",
        "used_flying"
      )
    )
  )

table(adamello_used_available_2scale$availability_class)

adamello_lc_long_2scale <- adamello_used_available_2scale %>%
  dplyr::select(
    availability_class,
    dplyr::all_of(frac_vars_adamello_2scale)
  ) %>%
  dplyr::group_by(availability_class) %>%
  dplyr::summarise(
    dplyr::across(
      dplyr::all_of(frac_vars_adamello_2scale),
      ~ mean(.x, na.rm = TRUE)
    ),
    .groups = "drop"
  ) %>%
  tidyr::pivot_longer(
    cols = dplyr::all_of(frac_vars_adamello_2scale),
    names_to = "class",
    values_to = "frac"
  ) %>%
  dplyr::mutate(
    class_code = as.integer(sub("^frac_", "", class))
  ) %>%
  dplyr::left_join(clc_lookup, by = "class_code") %>%
  dplyr::mutate(
    landcover_name = dplyr::if_else(
      is.na(landcover_name),
      class,
      landcover_name
    ),
    landcover_level1 = dplyr::if_else(
      is.na(landcover_level1),
      "Unmapped",
      landcover_level1
    )
  ) %>%
  dplyr::filter(
    !is.na(frac),
    frac > 0,
    !landcover_level1 %in% c("Outside area", "No data")
  )

adamello_lc_level1_2scale <- adamello_lc_long_2scale %>%
  dplyr::group_by(
    availability_class,
    landcover_level1
  ) %>%
  dplyr::summarise(
    frac = sum(frac, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::group_by(availability_class) %>%
  dplyr::mutate(
    prop = frac / sum(frac, na.rm = TRUE)
  ) %>%
  dplyr::ungroup()

p_adamello_lc_2scale <- ggplot2::ggplot(
  adamello_lc_level1_2scale,
  ggplot2::aes(
    x = availability_class,
    y = prop,
    fill = landcover_level1
  )
) +
  ggplot2::geom_col() +
  ggplot2::scale_fill_manual(
    values = lc_level1_colors,
    na.value = "grey80"
  ) +
  ggplot2::theme_minimal(base_size = 11) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 30, hjust = 1)
  ) +
  ggplot2::labs(
    x = NULL,
    y = "proportion of mean fractional cover",
    fill = "landcover level 1",
    title = paste0(
      "Used vs available landcover composition for ",
      target_id_clean
    ),
    subtitle = "Availability defined separately for terrestrial and flying behaviours"
  )

p_adamello_lc_2scale

ggplot2::ggsave(
  filename = file.path(plot_dir, "adamello.png"),
  plot = p_adamello_lc_2scale,
  width = 10,
  height = 8,
  dpi = 300
)