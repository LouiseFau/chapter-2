#' ---
#' title: "Nest site location per individuals"
#' author: "Louise Faure"
#' date of creation: 22.06.2026
#' details: (1) identify the location of the nest site for each individuals 
#' (2) for each individuals, defined an anthropophobic score based on the 
#' proportion of terrestrial behaviors expressed in the vicinity of human-dominated
#' landscape
#' (3) plot the result with the variability associated to the variation in the 
#' nbr of fix per individuals  
#' ---   

# library 
library(sp)
library(purrr)
library(dplyr)
library(ggplot2)

# data
output_dir <- "C:/Users/lfaure7/Documents/git/chapter-2/preparation donnee aigle ssf/resultats-intermediaire"

# emigration dates and golden eagle data
emig_dates <- readRDS("C:/Users/lfaure7/Documents/git/chapter-2/preparation donnee aigle ssf/donnees/emigration dates/emigration_dates_20250417.rds")

rds_dir <- "C:/Users/lfaure7/Documents/git/chapter-2/preparation donnee aigle ssf/donnees/no_burst_GE/no_burst_GE"
rds_files <- tools::file_path_sans_ext(list.files(rds_dir, pattern = "\\.rds$", ignore.case = TRUE))
rds_ids   <- as.numeric(gsub("_gpsNoDup_moveObj", "", rds_files))

# Filters individuals to retain those that have an emigration date
emig_dates_filtered <- emig_dates[emig_dates$did_disperse == TRUE &
                                    emig_dates$individual.id %in% rds_ids, ]
rds_files_filtered  <- rds_files[rds_ids %in% emig_dates_filtered$individual.id]


#' #############################################################################
#' ### Step 1 : identify the location of the nesting site 
#' #############################################################################


# Number of first GPS fixes used to estimate the nest location
nest_n_first_fixes <- 30

# Output folder for nest-site locations
nest_out_dir <- "C:/Users/lfaure7/Documents/git/chapter-2/preparation donnee aigle ssf/donnees/nest site location/nest_site_location"
dir.create(nest_out_dir, recursive = TRUE, showWarnings = FALSE)

# Helper function to extract coordinates and timestamps from one Move object
extract_move_data <- function(file_stub, rds_dir) {
  
  rds_path <- file.path(rds_dir, paste0(file_stub, ".rds"))
  mv <- readRDS(rds_path)
  
  individual_id <- as.numeric(gsub("_gpsNoDup_moveObj", "", file_stub))
  
  coords <- sp::coordinates(mv)
  ts <- move::timestamps(mv)
  
  # Try to recover individual local identifier from move object metadata
  id_dat <- tryCatch(
    as.data.frame(move::idData(mv)),
    error = function(e) data.frame()
  )
  
  individual_local <- if ("individual.local.identifier" %in% names(id_dat)) {
    as.character(id_dat$individual.local.identifier[1])
  } else if ("individual.local.identifier" %in% names(mv@data)) {
    as.character(mv@data$individual.local.identifier[1])
  } else {
    as.character(file_stub)
  }
  
  # Try to recover CRS from move object
  mv_crs <- tryCatch(
    sf::st_crs(sp::proj4string(mv)),
    error = function(e) sf::st_crs(4326)
  )
  
  if (is.na(mv_crs)) {
    mv_crs <- sf::st_crs(4326)
  }
  
  df <- data.frame(
    individual.id = individual_id,
    individual.local.identifier = individual_local,
    timestamp = as.POSIXct(ts, tz = "UTC"),
    x = coords[, 1],
    y = coords[, 2]
  )
  
  pts <- sf::st_as_sf(
    df,
    coords = c("x", "y"),
    crs = mv_crs,
    remove = FALSE
  )
  
  return(pts)
}

# Estimate nest site from first GPS fixes
estimate_nest_site <- function(file_stub, rds_dir, n_first = 30) {
  
  pts <- extract_move_data(file_stub, rds_dir)
  
  pts_first <- pts %>%
    dplyr::arrange(timestamp) %>%
    dplyr::slice_head(n = n_first)
  
  if (nrow(pts_first) == 0) {
    return(NULL)
  }
  
  # Work in metric CRS for averaging coordinates
  pts_3035 <- sf::st_transform(pts_first, 3035)
  xy <- sf::st_coordinates(pts_3035)
  
  nest_x <- mean(xy[, "X"], na.rm = TRUE)
  nest_y <- mean(xy[, "Y"], na.rm = TRUE)
  
  nest_sf <- sf::st_sf(
    individual.id = unique(pts_first$individual.id)[1],
    individual.local.identifier = unique(pts_first$individual.local.identifier)[1],
    n_fixes_used = nrow(pts_first),
    first_timestamp = min(pts_first$timestamp, na.rm = TRUE),
    last_timestamp_used = max(pts_first$timestamp, na.rm = TRUE),
    geometry = sf::st_sfc(sf::st_point(c(nest_x, nest_y)), crs = 3035)
  )
  
  return(nest_sf)
}

# Apply to all retained individuals
nest_sites_sf <- purrr::map(
  rds_files_filtered,
  estimate_nest_site,
  rds_dir = rds_dir,
  n_first = nest_n_first_fixes
) %>%
  purrr::compact() %>%
  do.call(rbind, .)

# Add lon/lat coordinates as attributes
nest_sites_lonlat <- sf::st_transform(nest_sites_sf, 4326)
nest_lonlat_xy <- sf::st_coordinates(nest_sites_lonlat)

nest_sites_sf$nest_longitude <- nest_lonlat_xy[, "X"]
nest_sites_sf$nest_latitude  <- nest_lonlat_xy[, "Y"]

# Save as RDS, CSV and shapefile
saveRDS(
  nest_sites_sf,
  file.path(nest_out_dir, "nest_site_location.rds")
)

sf::st_write(
  nest_sites_sf,
  file.path(nest_out_dir, "nest_site_location.shp"),
  delete_dsn = TRUE
)

print(nest_sites_sf)



################################################################################
#' ### STEP 2 : anthropophobia score from terrestrial behavior near humans 
#' 
#' ** Philosophy ** : proportion of terrestrial behavior within 100 of anthropised
#' areas divided by all the behavioral fixes within 100m of anthropised areas. A 
#' high score indicate that individuals rarely expresses terrestrial behavior near
#' humans while a low score indicates that individuals often expresses terrestrial
#' behavior near humans. 
#' 
#' The uncertainty associated to the variation in the nbr of gps fixes within a 
#' 100m buffer from highly human dominated areas is ploted. 
################################################################################


# Behaviors considered as terrestrial
terrestrial_behaviors <- c(
  "overnight roosting",
  "short resting"
)

# Minimum number of retained GPS fixes within 100 m of anthropised areas
min_n_near100 <- 30

# Wilson confidence interval for a binomial proportion
wilson_ci <- function(x, n, conf = 0.95) {
  z <- qnorm(1 - (1 - conf) / 2)
  p <- x / n
  
  denom <- 1 + z^2 / n
  centre <- p + z^2 / (2 * n)
  half_width <- z * sqrt((p * (1 - p) + z^2 / (4 * n)) / n)
  
  lower <- (centre - half_width) / denom
  upper <- (centre + half_width) / denom
  
  data.frame(lower = lower, upper = upper)
}

# Compute score
anthropophobia_score <- dispersal_thinned %>%
  dplyr::filter(
    !is.na(behavior_refined)
  ) %>%
  dplyr::group_by(individual.id, individual.local.identifier) %>%
  dplyr::summarise(
    n_total_behavioral_fixes = dplyr::n(),
    
    n_near100 = sum(near_human_100m %in% TRUE, na.rm = TRUE),
    
    n_terrestrial_near100 = sum(
      near_human_100m %in% TRUE &
        behavior_refined %in% terrestrial_behaviors,
      na.rm = TRUE
    ),
    
    prop_near100 = n_near100 / n_total_behavioral_fixes,
    
    prop_terrestrial_near100 = n_terrestrial_near100 / n_near100,
    
    anthropophobia_score = 1 - prop_terrestrial_near100,
    
    .groups = "drop"
  ) %>%
  dplyr::filter(
    n_near100 >= min_n_near100,
    is.finite(prop_terrestrial_near100)
  )

# Add Wilson confidence intervals
ci_terr <- wilson_ci(
  x = anthropophobia_score$n_terrestrial_near100,
  n = anthropophobia_score$n_near100
)

anthropophobia_score <- anthropophobia_score %>%
  dplyr::mutate(
    prop_terrestrial_lwr = ci_terr$lower,
    prop_terrestrial_upr = ci_terr$upper,
    
    # Since anthropophobia = 1 - terrestrial proportion,
    # the CI is inverted.
    anthropophobia_lwr = 1 - prop_terrestrial_upr,
    anthropophobia_upr = 1 - prop_terrestrial_lwr,
    
    n_non_terrestrial_near100 = n_near100 - n_terrestrial_near100
  ) %>%
  dplyr::arrange(anthropophobia_score)

print(anthropophobia_score)

# export as a shapefile 
# Join anthropophobia score to nest-site locations
anthropophobia_score_sf <- nest_sites_sf %>%
  dplyr::left_join(
    anthropophobia_score,
    by = c("individual.id", "individual.local.identifier")
  )

# Save as GeoPackage: recommended, keeps full column names
sf::st_write(
  anthropophobia_score_sf,
  file.path(anthro_out_dir, "anthropophobia_score_nest_sites.gpkg"),
  delete_dsn = TRUE
)

# Shapefile-safe version: short column names, because shapefile truncates names
anthropophobia_score_shp <- anthropophobia_score_sf %>%
  dplyr::rename(
    ind_id   = individual.id,
    indiv    = individual.local.identifier,
    n_tot    = n_total_behavioral_fixes,
    n100     = n_near100,
    nterr100 = n_terrestrial_near100,
    pnear100 = prop_near100,
    pterr100 = prop_terrestrial_near100,
    anthro   = anthropophobia_score,
    terr_lwr = prop_terrestrial_lwr,
    terr_upr = prop_terrestrial_upr,
    anth_lwr = anthropophobia_lwr,
    anth_upr = anthropophobia_upr,
    nnonterr = n_non_terrestrial_near100
  )

sf::st_write(
  anthropophobia_score_shp,
  file.path(anthro_out_dir, "anthropophobia_score_nest_sites.shp"),
  delete_dsn = TRUE
)


# Plot uncertainty
# Clean individual names if needed
clean_individual_name <- function(x) {
  gsub("\\s*\\([^)]*\\)\\s*$", "", x)
}

anthropophobia_score <- anthropophobia_score %>%
  dplyr::mutate(
    individual_plot = factor(
      individual.local.identifier,
      levels = individual.local.identifier
    )
  )

p_anthropophobia <- ggplot(
  anthropophobia_score,
  aes(
    x = individual_plot,
    y = anthropophobia_score
  )
) +
  geom_point(aes(size = n_near100)) +
  geom_errorbar(
    aes(
      ymin = anthropophobia_lwr,
      ymax = anthropophobia_upr
    ),
    width = 0.2
  ) +
  coord_flip() +
  scale_x_discrete(labels = clean_individual_name) +
  labs(
    x = "Individual",
    y = "Anthropophobia score",
    size = "Number of fixes ≤ 100 m",
    title = "Individual variation in behavioral anthropophobia",
    subtitle = "Score = 1 - proportion of terrestrial behaviors within ≤ 100 m of anthropised areas"
  ) +
  theme_minimal(base_size = 11)

print(p_anthropophobia)
ggsave(
  filename = "C:/Users/lfaure7/Documents/git/chapter-2/preparation donnee aigle ssf/dossier de plots/resultats preliminaire 23.06/p_anthropophobia.png",
  plot = p_anthropophobia,
  width = 8,
  height = 5,
  dpi = 300
)
