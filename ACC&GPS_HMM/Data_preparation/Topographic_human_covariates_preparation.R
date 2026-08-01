#' -----------------------------------------------------------------------------
#' Title: Preparation of topographic and human raster layers ----
#' Authors: Louise Faure
#' Date: 31.07.26
#' Purpose:
#' (1) topographic covariates preparation
#' (2) human footprint index preparation



# STEP 1: prepare topographic covariates ----
#' **Note:** EU-DEM v1.1 tiles E30N20 and E40N20 have a native resolution of
#' 25 m and use the ETRS89-LAEA projection (EPSG:3035). The two tiles are
#' combined as a virtual raster and cropped to the study extent.
#'
#' Slope and Terrain Ruggedness Index are first calculated from the native DEM.
#' Both variables are then averaged within a circular window with a radius of
#' 20 native cells, corresponding to 500 m. TPI is calculated as elevation
#' minus mean surrounding elevation within a circular window with a radius of
#' 40 native cells, corresponding to 1000 m. The central cell is excluded from
#' the TPI neighborhood.
#'
#' Ridgelines are defined as cells with positive TPI that are also local TPI
#' maxima within a 3 x 3 neighborhood. Distance to the nearest ridgeline is
#' calculated after transferring ridgeline presence to the final 100 m grid.
#'
#' **Steps:**
#' (i) load DEM and calculate slope and ruggedness within a 20-cell radius
#' (ii) calculate TPI within a 40-cell radius and identify ridgelines
#' (iii) calculate distance to the nearest ridgeline
#' (iv) reduce resolution to 100 m and export all covariates as GeoTIFF files


library(terra)

# Options and paths ----
temporary_directory <- "/Users/louisefaure/Downloads/topography_temporary"
dir.create(temporary_directory, recursive = TRUE, showWarnings = FALSE)
terraOptions(tempdir = temporary_directory, memfrac = 0.15, memmax = 4, todisk = TRUE, progress = 1)

dem_file <- "/Users/louisefaure/Desktop/dossier sans titre/Rasters/dem.tif"
terrain_native_file <- file.path(temporary_directory, "terrain_native_25m.tif")
slope_25m_file <- file.path(temporary_directory, "slope_500m_25m.tif")
ruggedness_25m_file <- file.path(temporary_directory, "ruggedness_500m_25m.tif")
dem_mean_40cells_file <- file.path(temporary_directory, "dem_mean_1000m_25m.tif")
tpi_25m_file <- file.path(temporary_directory, "tpi_1000m_25m.tif")
tpi_max_file <- file.path(temporary_directory, "tpi_local_maximum_25m.tif")
ridgeline_25m_file <- file.path(temporary_directory, "ridgeline_25m.tif")
ridgeline_100m_file <- file.path(temporary_directory, "ridgeline_100m.tif")

slope_output_file <- "/Users/louisefaure/Downloads/slope_500m_100m.tif"
ruggedness_output_file <- "/Users/louisefaure/Downloads/ruggedness_500m_100m.tif"
tpi_output_file <- "/Users/louisefaure/Downloads/tpi_1000m_100m.tif"
ridgeline_distance_output_file <- "/Users/louisefaure/Downloads/distance_to_ridgeline_100m.tif"

write_options <- list(datatype = "FLT4S", gdal = c("TILED=YES", "COMPRESS=ZSTD", "PREDICTOR=3", "BIGTIFF=IF_SAFER"))

# Load the already merged DEM ----
dem <- rast(dem_file)
analysis_extent <- ext(3841300, 4846200, 2236600, 2860700)
tpi_radius_m <- 1000
analysis_extent_buffer <- ext(xmin(analysis_extent) - tpi_radius_m, xmax(analysis_extent) + tpi_radius_m, ymin(analysis_extent) - tpi_radius_m, ymax(analysis_extent) + tpi_radius_m)
dem_buffer <- crop(dem, analysis_extent_buffer, snap = "out")

# Circular focal windows ----
circular_window <- function(radius_cells) {
  offsets <- -radius_cells:radius_cells
  outer(offsets, offsets, function(y, x) ifelse(x^2 + y^2 <= radius_cells^2, 1, NA_real_))
}

window_40cells <- circular_window(40)
window_40cells[41, 41] <- NA_real_


# Tiled focal function ----
focal_tiled <- function(x, w, fun, filename, tile_size = 1000) {
  tile_directory <- tempfile("focal_tiles_", tmpdir = temporary_directory)
  dir.create(tile_directory)
  halo <- floor(nrow(w) / 2)
  core_extents <- getTileExtents(x, y = c(tile_size, tile_size), buffer = 0)
  buffered_extents <- getTileExtents(x, y = c(tile_size, tile_size), buffer = halo)
  tile_files <- file.path(tile_directory, sprintf("tile_%05d.tif", seq_len(nrow(core_extents))))
  for(i in seq_len(nrow(core_extents))) {
    message("Processing tile ", i, " of ", nrow(core_extents))
    input_tile <- crop(x, ext(as.numeric(buffered_extents[i, ])), snap = "out")
    result_tile <- focal(input_tile, w = w, fun = fun, na.rm = TRUE)
    result_tile <- crop(result_tile, ext(as.numeric(core_extents[i, ])), snap = "near")
    writeRaster(result_tile, tile_files[i], overwrite = TRUE, wopt = write_options)
    rm(input_tile, result_tile)
    gc()
    message("Finished tile ", i, " of ", nrow(core_extents))
  }
  result_vrt <- vrt(tile_files, filename = file.path(tile_directory, "result.vrt"), overwrite = TRUE)
  result <- writeRaster(result_vrt, filename, overwrite = TRUE, wopt = write_options)
  unlink(tile_directory, recursive = TRUE, force = TRUE)
  result
}


# Native slope and TRI using the eight adjacent cells ----
terrain_native <- terrain(dem_buffer, v = c("slope", "TRI"), neighbors = 8, unit = "degrees", filename = terrain_native_file, overwrite = TRUE, wopt = write_options)
slope_native <- terrain_native[["slope"]]
ruggedness_native <- terrain_native[["TRI"]]

# Reduce slope and TRI from 25 m to 100 m ----
aggregation_factor <- round(100 / res(dem_buffer)[1])

slope_100m_buffer <- aggregate(slope_native, fact = aggregation_factor, fun = "mean", na.rm = TRUE)
ruggedness_100m_buffer <- aggregate(ruggedness_native, fact = aggregation_factor, fun = "mean", na.rm = TRUE)

slope <- crop(slope_100m_buffer, analysis_extent, snap = "near", filename = slope_output_file, overwrite = TRUE, wopt = write_options)
ruggedness <- crop(ruggedness_100m_buffer, analysis_extent, snap = "near", filename = ruggedness_output_file, overwrite = TRUE, wopt = write_options)

# Reduce DEM to 100 m before calculating TPI ----
dem_100m_buffer_file <- file.path(temporary_directory, "dem_100m_buffer.tif")
dem_mean_1000m_file <- file.path(temporary_directory, "dem_mean_1000m_100m.tif")
tpi_z_1000m_file <- file.path(temporary_directory, "tpi_z_1000m_100m.tif")

dem_100m_buffer <- aggregate(dem_buffer, fact = aggregation_factor, fun = "mean", na.rm = TRUE, filename = dem_100m_buffer_file, overwrite = TRUE, wopt = write_options)

# Circular window with a 1,000-m radius at 100-m resolution ----
window_10cells <- circular_window(10)
window_10cells[11, 11] <- NA_real_

# Calculate TPI at 100 m ----
dem_mean_1000m <- focal_tiled(dem_100m_buffer, window_10cells, "mean", dem_mean_1000m_file, tile_size = 1000)
tpi_1000m_buffer <- writeRaster(dem_100m_buffer - dem_mean_1000m, tpi_output_file, overwrite = TRUE, wopt = write_options)

# Standardize TPI using the final analysis extent ----
tpi_analysis <- crop(tpi_1000m_buffer, analysis_extent, snap = "near")
tpi_statistics <- global(tpi_analysis, c("mean", "sd"), na.rm = TRUE)
tpi_mean <- tpi_statistics[1, "mean"]
tpi_sd <- tpi_statistics[1, "sd"]

tpi_z_1000m_buffer <- writeRaster((tpi_1000m_buffer - tpi_mean) / tpi_sd, tpi_z_1000m_file, overwrite = TRUE, wopt = write_options)

# Identify ridgelines following Weiss: TPI > one standard deviation ----
ridgeline_100m <- ifel(tpi_z_1000m_buffer > 1, 1, NA, filename = ridgeline_100m_file, overwrite = TRUE, datatype = "INT1U", gdal = c("TILED=YES", "COMPRESS=ZSTD", "BIGTIFF=IF_SAFER"))

# Calculate distance to ridgelines on the buffered raster ----
distance_to_ridgeline_buffer <- distance(ridgeline_100m, unit = "m")
distance_to_ridgeline <- crop(distance_to_ridgeline_buffer, analysis_extent, snap = "near", filename = ridgeline_distance_output_file, overwrite = TRUE, datatype = "FLT4S", gdal = c("TILED=YES", "COMPRESS=ZSTD", "PREDICTOR=3", "BIGTIFF=IF_SAFER"))

# Final raster stack ----
topographic_covariates <- rast(c(slope_output_file, ruggedness_output_file, ridgeline_distance_output_file))
names(topographic_covariates) <- c("slope_8_neighbors", "tri_8_neighbors", "distance_to_ridgeline")



# STEP 2: prepare settlement and population density rasters ----
#' **Note:** Building footprints were downloaded from the Overture Maps Buildings
#' dataset (https://docs.overturemaps.org/getting-data/) and previously rasterized
#' with GDAL as a binary 100 m res. raster.
#'
#' Resident population was downloaded from GHS-POP R2023A, epoch 2020, resolution
#' 100 m, World Mollweide projection
#' (https://human-settlement.emergency.copernicus.eu/ghs_pop2023.php), mosaicked
#' over the Alps and previously reprojected and aligned to the settlement raster
#' in EPSG:3035. 
#' 
#' **Aims:** calculate density of people and settlement within a circular window
#' containing 97 cells, corresponding to approximately 0.97 km2 and a radius of
#' approximately 564 m.
#'
#' **Steps:**
#' (i) calculate settlement density as the proportion of built cells within 0.97 km2
#' (ii) calculate population density in residents per km2 within the same window
#' (iii) save both rasters as GeoTIFF files in the Downloads directory



library(terra)

# Parameters ----
terraOptions(memfrac = 0.10, memmax = 2, todisk = TRUE, progress = 1)

# Data ----
pop <- rast("C:/Users/lfaure7/Downloads/HFI_100m/temporary/population_count_100m_3035.tif")
settlement <- rast("C:/Users/lfaure7/Downloads/HFI_100m/built_cells_100m.tif")

# Circular window of approximately 1 km2 ----
cell_offsets <- -6:6
radius_m <- sqrt(1e6 / pi)
weight_matrix <- outer(cell_offsets, cell_offsets, function(row_offset, column_offset) ifelse((row_offset * 100)^2 + (column_offset * 100)^2 <= radius_m^2, 1, NA_real_))

# Tiled focal calculation ----
focal_tiled <- function(x, w, filename, scale_factor = 1, tile_rows = 500) {
  tile_directory <- tempfile("terra_focal_")
  dir.create(tile_directory)
  on.exit(unlink(tile_directory, recursive = TRUE, force = TRUE), add = TRUE)
  halo <- floor(nrow(w) / 2)
  core_extents <- getTileExtents(x, y = c(tile_rows, ncol(x)), buffer = 0)
  buffered_extents <- getTileExtents(x, y = c(tile_rows, ncol(x)), buffer = halo)
  tile_files <- file.path(tile_directory, sprintf("tile_%03d.tif", seq_len(nrow(core_extents))))
  for (i in seq_len(nrow(core_extents))) {
    input_tile <- crop(x, ext(as.numeric(buffered_extents[i, ])), snap = "out")
    raw_file <- file.path(tile_directory, sprintf("raw_%03d.tif", i))
    output_tile <- focal(input_tile, w = w, fun = "mean", na.rm = TRUE, filename = raw_file, overwrite = TRUE, wopt = list(datatype = "FLT4S", gdal = c("TILED=YES", "COMPRESS=ZSTD", "PREDICTOR=3", "BIGTIFF=IF_SAFER")))
    output_tile <- crop(output_tile, ext(as.numeric(core_extents[i, ])), snap = "near") * scale_factor
    writeRaster(output_tile, tile_files[i], overwrite = TRUE, datatype = "FLT4S", gdal = c("TILED=YES", "COMPRESS=ZSTD", "PREDICTOR=3", "BIGTIFF=IF_SAFER"))
    unlink(raw_file, force = TRUE)
  }
  output_vrt <- vrt(tile_files, filename = file.path(tile_directory, "output.vrt"), overwrite = TRUE)
  writeRaster(output_vrt, filename, overwrite = TRUE, datatype = "FLT4S", gdal = c("TILED=YES", "COMPRESS=ZSTD", "PREDICTOR=3", "BIGTIFF=IF_SAFER"))
}

# Settlement density: proportion of built 100 m cells ----
settlement_density <- focal_tiled(settlement, weight_matrix, "C:/Users/lfaure7/Downloads/settlement_density_1km2_100m.tif")

# Population density: mean residents per 100 m cell divided by 0.01 km2 ----
population_density <- focal_tiled(pop, weight_matrix, "C:/Users/lfaure7/Downloads/population_density_1km2_100m.tif", scale_factor = 100)