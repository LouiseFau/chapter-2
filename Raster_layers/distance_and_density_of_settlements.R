#' ---
#' title: "Preparing Raster Data"
#' author: "Louise Faure"
#' date: 28.05.2026
#' details: (i) reduce the resolution of the input raster, (ii) calculate distance to each cells, (iii) compute density of building map
#' (iv) export each raster independently 
#' ---   

#' ## Preamble
# libraries
library(terra)
library(sf)

# raster data
settlement_raw <- rast("C:/Users/lfaure7/Desktop/COUCHES QGIS/settlements/buildings_raw_3035.tif")

# polygone for the alpine area
alpine_area <- vect("C:/Users/lfaure7/Desktop/COUCHES QGIS/alpine_area/alpine_area_corrected_geometry.shp")

################################################################################
#' Step 1 : proportion of built areas per km2 
################################################################################

# rasterize the alpine polygon to attribute NA to pixels outside the Alps
alps_mask <- rasterize(
  alpine_area,
  settlement_raw,
  field = 1,
  background = NA
)

# binary grid: NA outside Alps / 1 built cell /0 non-built cell within Alps
settlement_01 <- ifel(
  is.na(alps_mask),
  NA,
  ifel(is.na(settlement_raw), 0, 1)
)

# prepare a circular window of 1km2

res_m <- res(settlement_01)[1]

if (abs(res_m - 50) > 1e-6) {
  warning(paste("La résolution n'est pas exactement 50 m :", res_m))
}

radius_1km <- sqrt(1e6 / pi)

k <- ceiling(radius_1km / res_m)
d <- seq(-k, k) * res_m

w_1km_circle <- outer(
  d, d,
  function(x, y) ifelse(x^2 + y^2 <= radius_1km^2, 1, NA)
)

window_area_km2 <- sum(!is.na(w_1km_circle)) * res_m^2 / 1e6
print(window_area_km2)

# calculate the proportion of built cells per km2

built_prop_1km <- focal(
  settlement_01,
  w = w_1km_circle,
  fun = "mean",
  na.policy = "omit",
  na.rm = TRUE,
  filename = "C:/Users/lfaure7/Desktop/COUCHES QGIS/settlements/built_proportion_1km_window_50m.tif",
  overwrite = TRUE,
  wopt = list(
    gdal = c(
      "COMPRESS=ZSTD",
      "PREDICTOR=3",
      "ZSTD_LEVEL=9",
      "TILED=YES",
      "NUM_THREADS=ALL_CPUS",
      "BIGTIFF=YES"
    )
  )
)

################################################################################
#' ### Step 2 : Densité de bâtiments par fenêtre glissante de 1 km²
################################################################################

building_count_50m <- rast(
  "C:/Users/lfaure7/Desktop/COUCHES QGIS/settlements/building_count_50m.tif"
)

building_count_50m <- mask(building_count_50m, alps_mask)

building_sum_1km <- focal(
  building_count_50m,
  w = w_1km_circle,
  fun = "sum",
  na.policy = "omit",
  na.rm = TRUE,
  filename = "C:/Users/lfaure7/Desktop/COUCHES QGIS/settlements/building_count_sum_1km_window_50m.tif",
  overwrite = TRUE,
  wopt = list(
    gdal = c(
      "COMPRESS=ZSTD",
      "PREDICTOR=3",
      "ZSTD_LEVEL=9",
      "TILED=YES",
      "NUM_THREADS=ALL_CPUS",
      "BIGTIFF=YES"
    )
  )
)

building_density_1km <- building_sum_1km / window_area_km2

writeRaster(
  building_density_1km,
  "C:/Users/lfaure7/Desktop/COUCHES QGIS/settlements/building_density_count_per_km2_1km_window_50m.tif",
  overwrite = TRUE,
  wopt = list(
    gdal = c(
      "COMPRESS=ZSTD",
      "PREDICTOR=3",
      "ZSTD_LEVEL=9",
      "TILED=YES",
      "NUM_THREADS=ALL_CPUS",
      "BIGTIFF=YES"
    )
  )
)
