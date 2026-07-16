#' ----------------------------------------------------------------------------- 
# Title: Extraction of covariates values one thinned dataset ----
#' Authors : Louise Faure
#' Date : 16.07.26
#' Info : this script follow the data_preparation.R script where data are thinned
#' Purpose : 
#' (1) extract covariates below each gps points. Covariates all have the same 
#' crs and resolution.
#' (2) extract the mean, q75 and maximum of HFI in two different buffers (500m 
#' and 1000m)
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
regular_60_sf <- readRDS("/Users/louisefaure/Desktop/dossier sans titre/donnees filtree/GE_60_min_thinned.rds")
human_footprint <- terra::rast("/Users/louisefaure/Desktop/dossier sans titre/human_footprint_index_building_pop_builtprop_100m.tif")
ruggedness <- terra::rast("/Users/louisefaure/Desktop/dossier sans titre/ruggedness_TRI_100m.tif")
slope <- terra::rast("/Users/louisefaure/Desktop/dossier sans titre/slope_100m.tif")
dist_ridgeline <- terra::rast("/Users/louisefaure/Desktop/dossier sans titre/distance_to_ridgeline_100m.tif")
landcover <- terra::rast("/Users/louisefaure/Desktop/dossier sans titre/landcover_100m.tif")




#'------------------------------------------------------------------------------
# STEP 1: extract point-level environmental covariates ----
#'
#' (1) read the 60-min dataset;
#' (2) split the dataset by individual;
#' (3) calculate diel temporal covariates;
#' (4) create and project GPS points;
#' (5) extract raster values below GPS locations;
#' (6) recombine individuals;
#' (7) export the resulting dataframe.



# 1.1 Convert to a standard dataframe ----
regular_60_df <- regular_60_sf %>%
  sf::st_drop_geometry() %>%
  as.data.frame()


# 1.2 Split the dataset by individual ----
regular_60_list <- split(
  regular_60_df,
  as.character(
    regular_60_df$individual.local.identifier))


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
    
    
    # 1.3.1 Calculate diel temporal covariates ----
    df <- df %>%
      dplyr::mutate(
        timestamp = as.POSIXct(
          timestamp,
          tz = "UTC"
        ),
        
        time_local = lubridate::with_tz(
          timestamp,
          tzone = tz_loc
        ),
        
        decimal_hour =
          lubridate::hour(time_local) +
          lubridate::minute(time_local) / 60 +
          lubridate::second(time_local) / 3600,
        
        time_angle =
          2 * pi * decimal_hour / 24,
        
        cos_diel =
          cos(time_angle),
        
        sin_time =
          sin(time_angle)
      )
    
    
    # Stop if no point is available
    if (nrow(df) == 0) {
      return(df)
    }
    
    
    # 1.3.2 Create GPS points in longitude-latitude coordinates ----
    pts <- terra::vect(
      df,
      geom = c(
        "location.long",
        "location.lat"
      ),
      crs = "EPSG:4326"
    )
    
    
    # 1.3.3 Project GPS points once to the common raster CRS ----
    pts_raster <- terra::project(pts, terra::crs(human_footprint))
    
    
    # 1.3.4 Extract raster values below GPS locations ----
    df$hfi_point <- terra::extract(
      human_footprint,
      pts_raster,
      ID = FALSE
    )[, 1]
    
    df$ruggedness_100m <- terra::extract(
      ruggedness,
      pts_raster,
      ID = FALSE
    )[, 1]
    
    df$slope_100m <- terra::extract(
      slope,
      pts_raster,
      ID = FALSE
    )[, 1]
    
    df$distance_to_ridgeline_100m <- terra::extract(
      dist_ridgeline,
      pts_raster,
      ID = FALSE
    )[, 1]
    
    df$landcover_100m <- terra::extract(
      landcover,
      pts_raster,
      ID = FALSE
    )[, 1]
    
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




#'------------------------------------------------------------------------------
# STEP 2: extract buffered HFI covariates by individual ----
#'
#' (1) split the annotated dataframe by individual;
#' (2) recreate and project GPS points;
#' (3) create 500-m and 1000-m buffers;
#' (4) extract HFI mean, maximum and Q75 in each buffer;
#' (5) recombine individuals;
#' (6) export the final dataframe.





# 2.1 Split the annotated dataset by individual ----
GE_60_min_covariates_list <- split(GE_60_min_covariates,
  as.character(
    GE_60_min_covariates$individual.local.identifier))


# 2.2 Extract buffered HFI covariates for each individual ----
GE_60_min_hfi_list <- lapply(
  seq_along(GE_60_min_covariates_list),
  function(i) {
    
    individual_name <- names(
      GE_60_min_covariates_list
    )[i]
    
    cat(
      "Processing HFI buffers for individual",
      i,
      "of",
      length(GE_60_min_covariates_list),
      ":",
      individual_name,
      "\n")
    
    flush.console()
    
    # Select one individual
    df <- GE_60_min_covariates_list[[i]]
    
    # Stop if no point is available
    if (nrow(df) == 0) {return(df)}
    
    # 2.2.1 Create GPS points in longitude-latitude coordinates ----
    pts <- terra::vect(
      df,
      geom = c(
        "location.long",
        "location.lat"
      ),
      crs = "EPSG:4326")
    
    
    # 2.2.2 Project GPS points once to the raster CRS ----
    pts_raster <- terra::project(pts,terra::crs(human_footprint))
    
    # 2.2.3 Create 500-m buffers ----
    buffers_500m <- terra::buffer(pts_raster, width = 500)
    
    # Mean HFI within 500 m
    df$hfi_mean_500m <- terra::extract(
      human_footprint,
      buffers_500m,
      fun = mean,
      na.rm = TRUE,
      ID = FALSE
    )[, 1]
    
    # Maximum HFI within 500 m
    df$hfi_max_500m <- terra::extract(
      human_footprint,
      buffers_500m,
      fun = max,
      na.rm = TRUE,
      ID = FALSE
    )[, 1]
    
    
    # Third quartile of HFI within 500 m
    df$hfi_q75_500m <- terra::extract(
      human_footprint,
      buffers_500m,
      fun = function(x, ...) {
        
        stats::quantile(
          x,
          probs = 0.75,
          na.rm = TRUE,
          names = FALSE,
          type = 7
        )
      },
      ID = FALSE
    )[, 1]
    
    # Remove 500-m buffers before creating the next ones
    rm(buffers_500m)
    gc()
    
    # 2.2.4 Create 1000-m buffers ----
    buffers_1000m <- terra::buffer(pts_raster, width = 1000)
    
    # Mean HFI within 1000 m
    df$hfi_mean_1000m <- terra::extract(
      human_footprint,
      buffers_1000m,
      fun = mean,
      na.rm = TRUE,
      ID = FALSE
    )[, 1]
    
    # Maximum HFI within 1000 m
    df$hfi_max_1000m <- terra::extract(
      human_footprint,
      buffers_1000m,
      fun = max,
      na.rm = TRUE,
      ID = FALSE
    )[, 1]
    
    # Third quartile of HFI within 1000 m
    df$hfi_q75_1000m <- terra::extract(
      human_footprint,
      buffers_1000m,
      fun = function(x, ...) {
        
        stats::quantile(
          x,
          probs = 0.75,
          na.rm = TRUE,
          names = FALSE,
          type = 7
        )
      },
      ID = FALSE
    )[, 1]
    
    # Remove temporary spatial objects
    rm(
      pts,
      pts_raster,
      buffers_1000m
    )
    gc()
    
    
    # Return the annotated dataframe
    df
  }
)


# Preserve list names
names(GE_60_min_hfi_list) <- names(GE_60_min_covariates_list)


# 2.3 Recombine all individuals ----
GE_60_min_covariates_hfi <- GE_60_min_hfi_list %>%
  dplyr::bind_rows() %>%
  dplyr::arrange(
    individual.local.identifier,
    timestamp
  ) %>%
  as.data.frame()


# 2.4 Save the final dataset ----
saveRDS(
  GE_60_min_covariates_hfi,
  file = "/Users/louisefaure/Desktop/dossier sans titre/donnees filtree/GE_60_min_covariates_hfi.rds",
  compress = "gzip")