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


saveRDS(
  M_spearman,
  file.path(output_dir, "correlation_matrix_continuous_landscape_spearman.rds")
)

write.csv(
  round(M_spearman, 3),
  file.path(output_dir, "correlation_matrix_continuous_landscape_spearman.csv")
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
clc_lookup_index <- tibble::tibble(
  class_code = 1:44,
  clc_code = c(
    111, 112, 121, 122, 123, 124, 131, 132, 133, 141, 142,
    211, 212, 213, 221, 222, 223, 231, 241, 242, 243, 244,
    311, 312, 313, 321, 322, 323, 324, 331, 332, 333, 334, 335,
    411, 412, 421, 422, 423,
    511, 512, 521, 522, 523
  ),
  landcover_name = c(
    "Continuous urban fabric",
    "Discontinuous urban fabric",
    "Industrial or commercial units",
    "Road and rail networks and associated land",
    "Port areas",
    "Airports",
    "Mineral extraction sites",
    "Dump sites",
    "Construction sites",
    "Green urban areas",
    "Sport and leisure facilities",
    "Non-irrigated arable land",
    "Permanently irrigated land",
    "Rice fields",
    "Vineyards",
    "Fruit trees and berry plantations",
    "Olive groves",
    "Pastures",
    "Annual crops associated with permanent crops",
    "Complex cultivation patterns",
    "Agriculture with natural vegetation",
    "Agro-forestry areas",
    "Broad-leaved forest",
    "Coniferous forest",
    "Mixed forest",
    "Natural grasslands",
    "Moors and heathland",
    "Sclerophyllous vegetation",
    "Transitional woodland-shrub",
    "Beaches, dunes, sands",
    "Bare rocks",
    "Sparsely vegetated areas",
    "Burnt areas",
    "Glaciers and perpetual snow",
    "Inland marshes",
    "Peat bogs",
    "Salt marshes",
    "Salines",
    "Intertidal flats",
    "Water courses",
    "Water bodies",
    "Coastal lagoons",
    "Estuaries",
    "Sea and ocean"
  )
) %>%
  dplyr::mutate(
    landcover_level1 = dplyr::case_when(
      clc_code >= 100 & clc_code < 200 ~ "Artificial surfaces",
      clc_code >= 200 & clc_code < 300 ~ "Agricultural areas",
      clc_code >= 300 & clc_code < 400 ~ "Forest and semi-natural areas",
      clc_code >= 400 & clc_code < 500 ~ "Wetlands",
      clc_code >= 500 & clc_code < 600 ~ "Water bodies",
      TRUE ~ "Unclassified"
    )
  )

clc_lookup <- dplyr::bind_rows(
  clc_lookup_index,
  clc_lookup_index %>%
    dplyr::transmute(
      class_code = clc_code,
      clc_code = clc_code,
      landcover_name = landcover_name,
      landcover_level1 = landcover_level1
    )
) %>%
  dplyr::distinct(class_code, .keep_all = TRUE)

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
    class_code = as.integer(as.numeric(sub("^frac_", "", class)))
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

# Detailed CLC class composition per individual
p_lc_individual <- ggplot2::ggplot(
  lc_long,
  ggplot2::aes(
    x = factor(individual.local.identifier),
    y = frac,
    fill = landcover_name
  )
) +
  ggplot2::geom_col() +
  ggplot2::coord_flip() +
  ggplot2::facet_grid(. ~ behavior_type) +
  ggplot2::theme_minimal(base_size = 11) +
  ggplot2::labs(
    x = "individual",
    y = "mean fractional cover in 100 m buffer",
    fill = "landcover class",
    title = "Landcover composition by individual and behaviour"
  )

p_lc_individual

# Mean CLC class composition by behaviour
lc_behavior_mean <- lc_long %>%
  dplyr::group_by(behavior_type, landcover_name, landcover_level1) %>%
  dplyr::summarise(
    mean_frac = mean(frac, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::group_by(behavior_type) %>%
  dplyr::mutate(
    mean_frac = mean_frac / sum(mean_frac, na.rm = TRUE)
  ) %>%
  dplyr::ungroup()

p_lc_behavior_mean <- ggplot2::ggplot(
  lc_behavior_mean,
  ggplot2::aes(
    x = behavior_type,
    y = mean_frac,
    fill = landcover_name
  )
) +
  ggplot2::geom_col() +
  ggplot2::theme_minimal(base_size = 11) +
  ggplot2::labs(
    x = "behaviour",
    y = "mean fractional cover",
    fill = "landcover class",
    title = "Mean landcover composition: terrestrial vs flying"
  )

p_lc_behavior_mean

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
  ggplot2::theme_minimal(base_size = 11) +
  ggplot2::labs(
    x = "individual",
    y = "mean fractional cover in 100 m buffer",
    fill = "CLC level 1",
    title = "Landcover composition by broad CLC class"
  )

p_lc_level1



#' end alternative step 3






























































#' STEP 3 : Landscape characteristics under terrestrial behaviour
#'   3a  spatial thinning of terrestrial locations (per individual)
#'   3b  buffer extraction of covariates (area-weighted mean / fractional)
#'   3c  exploratory visualisation (population & individual level)
#'   3d  OPTIONAL available/background sample (use vs. availability context)
#' #############################################################################


#' Projected CRS in METRES for all distance / buffer operations.
crs_proj <- terra::crs(dem)

#' ### STEP 3a : spatial thinning ----------------------------------------------
# filter dispersal data to retain only the terrestrial points
terr <- dispersal_data %>%
  filter(behavior_cluster == 1,
         !is.na(location.long), !is.na(location.lat))

# transform to a spatial object with a crs
terr_sf <- st_as_sf(terr,
                    coords = c("location.long", "location.lat"),
                    crs = 4326, remove = FALSE) %>%
  st_transform(crs_proj)

## calculate distance between time-ordered fixes per individuals
consec_dist <- function(geom) {
  n <- length(geom)
  if (n < 2) return(rep(NA_real_, n))
  c(NA, as.numeric(st_distance(geom[-n], geom[-1], by_element = TRUE)))
}

## fill step_m per individual, ordered by time
terr_sf$step_m <- NA_real_
for (id in unique(terr_sf$individual.local.identifier)) {
  ix <- which(terr_sf$individual.local.identifier == id)
  ix <- ix[order(terr_sf$timestamp[ix])]
  terr_sf$step_m[ix] <- consec_dist(st_geometry(terr_sf)[ix])
}

## inspect the per-individual distance structure BEFORE fixing a threshold.
med_step <- terr_sf %>%
  st_drop_geometry() %>%
  group_by(individual.local.identifier) %>%
  summarise(median_step_m = median(step_m, na.rm = TRUE),
            q25 = quantile(step_m, .25, na.rm = TRUE),
            q75 = quantile(step_m, .75, na.rm = TRUE),
            n   = n(), .groups = "drop")
print(med_step)

## filter the point to keep only location when they are separated of at least d. distance (e.i, 25m)
greedy_thin_idx <- function(xy, d) {
  keep <- logical(nrow(xy))
  kept <- matrix(nrow = 0, ncol = 2)
  for (i in seq_len(nrow(xy))) {
    if (nrow(kept) == 0 ||
        all(sqrt((kept[, 1] - xy[i, 1])^2 + (kept[, 2] - xy[i, 2])^2) >= d)) {
      keep[i] <- TRUE
      kept <- rbind(kept, xy[i, ])
    }
  }
  keep
}

d_thin <- median(med_step$median_step_m, na.rm = TRUE)   # e.g. override: d_thin <- 50

terr_sf$keep <- FALSE
for (id in unique(terr_sf$individual.local.identifier)) {
  ix <- which(terr_sf$individual.local.identifier == id)
  ix <- ix[order(terr_sf$timestamp[ix])]
  xy <- st_coordinates(terr_sf)[ix, , drop = FALSE]
  terr_sf$keep[ix] <- greedy_thin_idx(xy, d_thin)
}

terr_thin <- terr_sf %>% filter(keep)
message(sprintf("Thinning at d = %.0f m : %d -> %d points",
                d_thin, nrow(terr_sf), nrow(terr_thin)))


#' ### STEP 3b : buffer extraction ---------------------------------------------

## Buffer radius (metres)
buf_r   <- 100
buffers <- st_buffer(terr_thin, dist = buf_r)   # buffers are in crs_proj (metres)

## extract categorical and numerical values
extract_into <- function(buffers_proj, r, categorical = FALSE, fun = "mean") {
  b <- st_transform(buffers_proj, terra::crs(r))   # match each raster's CRS
  if (categorical) {
    exact_extract(r, b, "frac", progress = FALSE)
  } else {
    exact_extract(r, b, fun, progress = FALSE)
  }
}

terr_thin$elevation  <- extract_into(buffers, covar_1_elevation)
terr_thin$dist_ridgeline  <- extract_into(buffers, covar_2_distance_ridgeline)
terr_thin$tri             <- extract_into(buffers, covar_3_TRI)
terr_thin$slope           <- extract_into(buffers, covar_4_slope)
terr_thin$tpi             <- extract_into(buffers, covar_5_TPI)
terr_thin$dist_settlement <- extract_into(buffers, covar_6_distance_settlement)
terr_thin$dens_settlement <- extract_into(buffers, covar_7_density_settlement)
terr_thin$pop_dens <- extract_into(buffers, covar_8_density_pop_km2)

lc_frac <- extract_into(buffers, covar_9_landcover, categorical = TRUE)  # frac_* cols

covar_tbl <- terr_thin %>% st_drop_geometry() %>% bind_cols(lc_frac)
saveRDS(covar_tbl, file.path(output_dir, "terrestrial_covariates.rds"))


#' ### STEP 3c : exploratory visualisation -------------------------------------

## (1) RIDGELINES: full distribution shape, one row per individual, per variable.
##     Reveals skew / bimodality that a boxplot hides (e.g. a long right tail in
##     distance-to-settlement, or two modes = two roost "types").
human_vars <- c("dist_settlement", "dens_settlement", "pop_dens")
long_human <- long %>%
  dplyr::filter(variable %in% human_vars)
p_ridge_human <- ggplot(
  long_human,
  aes(x = value, y = individual.id, fill = variable)
) +
  geom_density_ridges(
    scale = 1.2,
    alpha = .6,
    rel_min_height = .01,
    quantile_lines = TRUE,
    quantiles = 2
  ) +
  facet_wrap(~ variable, scales = "free_x") +
  guides(fill = "none") +
  theme_minimal(base_size = 11) +
  labs(
    y = "individual",
    x = "human-related covariate value at terrestrial locations"
  )

p_ridge_human



## (2) COLLINEARITY: slope/TRI/TPI and the three human-footprint layers are very
##     likely correlated. Check before any modelling.
M <- cor(covar_tbl[, cont_vars], use = "pairwise.complete.obs", method = "spearman")
# corrplot::corrplot(M, method = "number", type = "upper")
print(round(M, 2))

matrix of correlation




## (3) LANDCOVER composition per individual (mean fractional cover).
##     Consider reclassifying CLC codes to level-1 / human-vs-natural first.
lc_long <- covar_tbl %>%
  dplyr::select(individual.local.identifier, dplyr::starts_with("frac_")) %>%
  dplyr::group_by(individual.local.identifier) %>%
  dplyr::summarise(
    dplyr::across(dplyr::starts_with("frac_"), ~ mean(.x, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  tidyr::pivot_longer(
    cols = -individual.local.identifier,
    names_to = "class",
    values_to = "frac"
  )

p_lc <- ggplot(lc_long, aes(factor(individual.local.identifier), frac, fill = class)) +
  geom_col() + coord_flip() +
  theme_minimal(base_size = 11) +
  labs(x = "individual", y = "mean fractional cover in buffer")




















# Step selection functions preparation #########################################

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
  dplyr::select(id, steps) %>%
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

#' **Philosophy:** use is defined when the end of a step correspond to a terrestrial 
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
covar_1_elevation <- terra::rast("C:/Users/lfaure7/OneDrive/MEMOIRE M2/eagle_projet/data/pretraitements/Region-Alpes-Dem/region-alpes-dem.tif")
covar_2_distance_ridgeline <- terra::rast("C:/Users/lfaure7/OneDrive/MEMOIRE M2/donnees/raster/topography/distance_to_ridge_line_complete_version.tif")
covar_3_TRI <- terra::rast("C:/Users/lfaure7/OneDrive/MEMOIRE M2/eagle_projet/data/pretraitements/DEM_25m/TRI/TRI.tif")
covar_4_slope <- terra::rast("C:/Users/lfaure7/OneDrive/MEMOIRE M2/eagle_projet/data/pretraitements/DEM_25m/slope/slope_25.tif")
covar_5_TPI <- terra::rast("C:/Users/lfaure7/OneDrive/MEMOIRE M2/eagle_projet/data/pretraitements/DEM_25m/TPI/Topographic Position Index.tif")
covar_6_distance_settlement <- terra::rast("C:/Users/lfaure7/Desktop/COUCHES QGIS/settlements/distance to settlements 30 m cropped/dist_to_settlements_30m_croped.tif")
covar_7_density_settlement <- terra::rast("C:/Users/lfaure7/Desktop/COUCHES QGIS/settlements/density_rad_495m.tif")
covar_8_density_pop_km2 <- terra::rast("C:/Users/lfaure7/Desktop/COUCHES QGIS/density-resident-pop-km2/density_pop_km2.tif")
covar_9_landcover <- terra::rast("C:/Users/lfaure7/OneDrive/THESE/Conférences/Sempach workshop/Donnees/Landscape layers/CLC_longlat_10m.tif")

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