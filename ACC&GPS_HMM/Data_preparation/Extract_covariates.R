#' ----------------------------------------------------------------------------- 
# Title: Extraction of covariates values one thinned dataset ----
#' Authors : Louise Faure
#' Date : 16.07.26
#' Info : this script follow the data_preparation.R script where data are thinned
#' **Purpose** : extract covariates below each gps points. Covariates all have the same 
#' crs and resolution. **Steps**:
#' (1) read the 60-min dataset;
#' (2) split the dataset by individual;
#' (3) calculate diel temporal covariates;
#' (4) create and project GPS points;
#' (5) extract raster values below GPS locations;
#' (6) recombine individuals;
#' (7) export the resulting dataframe.
#' -----------------------------------------------------------------------------


# library
library(move2)
library(sf)
library(terra)
library(dplyr)
library(tidyverse)
library(lubridate)
library(data.table)

# parameters
tz_loc <- "Europe/Zurich"
terra::terraOptions(threads = 5)

# golden eagle thinned data and rasters
# for acc data
regular_60_sf <- readRDS("/Users/louisefaure/Library/CloudStorage/OneDrive-Personnel/THESE/CHAPITRE 2/git/chapter-2/HMM/HMM on ACC-classified behaviors/donnees intermediaire (2)/GE_60_min_thinned_behavior_assigned2.rds") 
# for gps data 
regular_60_sf <- readRDS("/Users/louisefaure/Library/CloudStorage/OneDrive-Personnel/THESE/CHAPITRE 2/git/chapter-2/HMM/preparation HMM/donnees intermediaire/GE_60_min_thinned.rds")

# topographic and humans covariates
raster_directory <- "/Users/louisefaure/Desktop/dossier sans titre/Rasters"
raster_files <- list.files(raster_directory, pattern = "\\.(tif|tiff|vrt|img)$", full.names = TRUE, ignore.case = TRUE)
rasters <- lapply(raster_files,function(file){x <- rast(file)
  crs(x) <- "EPSG:3035"
  x})
names(rasters) <- basename(raster_files)

# names for each raster 
population_density <- rasters[["population_density_1km2_100m.tif"]]
settlement_density <- rasters[["settlement_density_1km2_100m.tif"]]
elevation <- rasters[["elevation_100m.tif"]]
ruggedness <- rasters[["ruggedness_100m.tif"]]
slope <- rasters[["slope100m.tif"]]
dist_ridgeline <- rasters[["distance_to_ridgeline_100m.tif"]]
landcover <- rasters[["landcover_100m.tif"]]

# continuous versus categorical values 
continuous_covariates <- list(
  population_density = population_density,
  settlement_density = settlement_density,
  elevation_100m = elevation,
  ruggedness_100m = ruggedness,
  slope_100m = slope,
  distance_to_ridgeline_100m = dist_ridgeline
)
names(continuous_covariates) <- c("population_density","settlement_density","elevation_100m","ruggedness_100m","slope_100m","distance_to_ridgeline_100m")


# 1.1 Convert to a standard dataframe ----
regular_60_df <- regular_60_sf %>%sf::st_drop_geometry() %>%as.data.frame()

# 1.2 Split the dataset by individual and prepare function to extract landcover percentage ----
regular_60_list <- split(regular_60_df,as.character(regular_60_df$individual.local.identifier))

# Extract land-cover proportions from the central raster cell
# and its four orthogonal neighbours.
extract_landcover_5cells <- function(landcover_raster,points_raster){
  central_cells <- terra::cellFromXY(landcover_raster,terra::crds(points_raster))
  central_row_col <- terra::rowColFromCell(landcover_raster,central_cells)
  offsets <- rbind(c(0,0),c(-1,0),c(1,0),c(0,-1),c(0,1))
  rows <- outer(central_row_col[,1],offsets[,1],"+")
  cols <- outer(central_row_col[,2],offsets[,2],"+")
  valid <- rows >= 1 & rows <= terra::nrow(landcover_raster) & cols >= 1 & cols <= terra::ncol(landcover_raster)
  cells_5 <- matrix(NA_real_,nrow = nrow(rows),ncol = ncol(rows))
  cells_5[valid] <- terra::cellFromRowCol(landcover_raster,rows[valid],cols[valid])
  values_5 <- matrix(NA_real_,nrow = nrow(cells_5),ncol = ncol(cells_5))
  valid_cells <- !is.na(cells_5)
  values_5[valid_cells] <- landcover_raster[as.vector(cells_5[valid_cells])][,1]
  complete_neighbourhood <- rowSums(!is.na(values_5)) == 5
  n_forest <- rowSums(values_5 == 2 | values_5 == 3 | values_5 == 4,na.rm = TRUE)
  n_low_vegetation <- rowSums(values_5 == 5 | values_5 == 6 | values_5 == 7,na.rm = TRUE)
  n_rocky_terrain <- rowSums(values_5 == 8 | values_5 == 9 | values_5 == 11,na.rm = TRUE)
  n_other <- rowSums(!is.na(values_5),na.rm = TRUE) - n_forest - n_low_vegetation - n_rocky_terrain
  data.frame(
    prop_forest_5cells = ifelse(complete_neighbourhood,n_forest / 5,NA_real_),
    prop_low_vegetation_5cells = ifelse(complete_neighbourhood,n_low_vegetation / 5,NA_real_),
    prop_rocky_terrain_5cells = ifelse(complete_neighbourhood,n_rocky_terrain / 5,NA_real_),
    prop_other_5cells = ifelse(complete_neighbourhood,n_other / 5,NA_real_)
  )
}

# 1.3 Process each individual ----
regular_60_annotated_list <- lapply(
  seq_along(regular_60_list),
  function(i) {
    
    individual_name <- names(regular_60_list)[i]
    
    cat(
      "Processing individual",
      i,
      "of",
      length(regular_60_list),
      ":",
      individual_name,
      "\n"
    )
    
    flush.console()
    
    
    # Select one individual
    df <- regular_60_list[[i]]
    
    
    # 1.3.1 Calculate diel temporal covariates using apparent solar time ----
    if (!all(c("cos_diel", "sin_time") %in% names(df))) {
      
      df <- df %>%
        dplyr::mutate(
          timestamp = as.POSIXct(
            timestamp,
            tz = "UTC"
          ),
          
          # Decimal UTC time
          utc_decimal_hour =
            lubridate::hour(timestamp) +
            lubridate::minute(timestamp) / 60 +
            lubridate::second(timestamp) / 3600,
          
          # Fractional year for equation of time
          days_in_year =
            dplyr::if_else(
              lubridate::leap_year(timestamp),
              366,
              365
            ),
          
          fractional_year =
            2 * pi / days_in_year *
            (
              lubridate::yday(timestamp) - 1 +
                (utc_decimal_hour - 12) / 24
            ),
          
          # Equation of time (minutes)
          equation_of_time_min =
            229.18 *
            (
              0.000075 +
                0.001868 * cos(fractional_year) -
                0.032077 * sin(fractional_year) -
                0.014615 * cos(2 * fractional_year) -
                0.040849 * sin(2 * fractional_year)
            ),
          
          # Apparent solar time (minutes)
          # Longitude positive eastward
          solar_time_min =
            (
              utc_decimal_hour * 60 +
                equation_of_time_min +
                4 * lon
            ) %% 1440,
          
          solar_decimal_hour =
            solar_time_min / 60,
          
          solar_time_angle =
            2 * pi * solar_decimal_hour / 24,
          
          cos_diel =
            cos(solar_time_angle),
          
          sin_time =
            sin(solar_time_angle)
        )
      
    } else {
      
      message(
        "cos_diel and sin_time already present: solar diel covariates not recalculated."
      )
      
    }
    
    
    # 1.3.2 Create GPS points in longitude-latitude coordinates ----
    pts <- terra::vect(df,geom = c("lon","lat"),crs = "EPSG:4326")
    pts_raster <- terra::project(pts,"EPSG:3035")
    
    # 1.3.3 extract numerical covariates 
    covariate_values <- lapply(continuous_covariates,function(x){
      terra::extract(x,pts_raster,ID = FALSE)[,1]
    })
    
    covariate_values <- as.data.frame(covariate_values)
    df <- dplyr::bind_cols(df,covariate_values)
    
    # 1.3.5 Extract land-cover proportions from five cells ----
    landcover_5cells <- extract_landcover_5cells(
      landcover_raster = landcover,
      points_raster = pts_raster)
    
    df <- dplyr::bind_cols(
      df,
      landcover_5cells)
    
    # Remove temporary spatial objects
    rm(pts,pts_raster)
    
    # Return the annotated dataframe
    df
  }
)


# Preserve list names
names(regular_60_annotated_list) <- names(regular_60_list)


# 1.5 Recombine all individuals ----
GE_60_min_covariates <- regular_60_annotated_list %>%
  dplyr::bind_rows() %>%
  dplyr::arrange(
    individual.local.identifier,
    timestamp
  ) %>%
  as.data.frame()

# 1.6 Retain relevant columns ----
columns_to_keep_60 <- c(
  "sensor_type_id",
  "individual.local.identifier",
  "eobs_horizontal_accuracy_estimate",
  "eobs_speed_accuracy_estimate",
  "eobs_type_of_fix",
  "gps_dop",
  "gps_satellite_count",
  "ground_speed",
  "height_above_ellipsoid",
  "timestamp",
  "lon",
  "lat",
  "dist.traveled",
  "rf8fitted", "behavior_binary",
  "behavior_assignment_method",
  "behavior_base",
  "behavior_reclassified",
  "behavior_reclassification_rule",
  "burst_n",
  "burst_id",
  "row_in_burst",
  
  # age 
  "age_since_emig_weeks", 
  "age_since_emig_days", 
  "distance_to_nest_km",
  
  # Diel temporal covariates
  "cos_diel",
  "sin_time",
  
  # Point-level environmental covariates
  "hfi_point",
  "dem_elevation",
  "ruggedness_100m",
  "slope_100m",
  "distance_to_ridgeline_100m",
  "elevation_100m",
  
  # Land-cover proportions from five raster cells
  "prop_forest_5cells",
  "prop_low_vegetation_5cells",
  "prop_rocky_terrain_5cells",
  "prop_other_5cells",
  
  # human covariates
  "population_density",
  "settlement_density")


GE_60_min_covariates <- GE_60_min_covariates %>%
  dplyr::select(
    dplyr::any_of(
      columns_to_keep_60))



# 1.7 Save the final dataset ----
# change the name for saving acc data
saveRDS(
  GE_60_min_covariates,
  file = "/Users/louisefaure/Desktop/dossier sans titre/donnees filtree/GE_acc_60_min_covariates_hfi(2).rds",
  compress = "gzip")
