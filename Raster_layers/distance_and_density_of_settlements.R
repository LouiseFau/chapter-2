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
binary_settlement <- rast("C:/Users/lfaure7/Desktop/COUCHES QGIS/binary_mask.tif")

# polygone for the alpine area
alpine_area <- vect("C:/Users/lfaure7/Desktop/COUCHES QGIS/alpine_area/alpine_area_corrected_geometry.shp")

# Resolution reduction : for each 30m cells (env. 16 pixel), count the number of cell with a value 1  
built_count <- aggregate(binary_settlement, fact = 4, fun = "sum", na.rm = TRUE)

# Filter and obtain a binary raster : give a 1 value when the 30 m cell is >= 4 
binary_agg2 <- ifel(built_count >= 4, 1, NA)

# Calculate density over 9 neighboring cells 
nb <- focal(ifel(is.na(binary_agg2), 0, 1), w = 3, fun = "sum")

# Filter : remove the cell that have not at least 2 built neighbours, this allow to remove the very isolated settlements
binary_agg3 <- ifel(!is.na(binary_agg2) & nb >= 2, 1, NA)

# Calculate distance and clip to the alps
rdist <- distance(binary_agg3,
                  filename  = "dist_to_buildings.tif",
                  overwrite = TRUE)

rdist_clip <- crop(rdist, alpine_area, mask = TRUE)
plot(rdist_clip)

# checks the values 
quantile(values(rdist_clip, mat = FALSE),
         probs = c(.5, .75, .9, .95, .99, 1), na.rm = TRUE)

# export the raster croped to the alpine area
writeRaster(rdist_clip, "C:/Users/lfaure7/Desktop/COUCHES QGIS/settlements/distance to settlements 30 m cropped/dist_to_settlements_30m_croped.tif",
            overwrite = TRUE,
            datatype  = "INT4U",
            gdal = c(
              "COMPRESS=ZSTD",        
              "PREDICTOR=2",          
              "ZSTD_LEVEL=9",
              "TILED=YES",            
              "NUM_THREADS=ALL_CPUS",
              "BIGTIFF=YES"           
            )
)


#' ## Density of settlements map
alps_mask     <- rasterize(alpine_area, built_count, field = 1, background = NA)
settlement_01 <- ifel(is.na(alps_mask), NA,
                      ifel(is.na(binary_agg2), 0, 1))   # cohérent avec la distance

density_local <- focal(settlement_01, w = 17, fun = "mean",
                       na.policy = "omit", na.rm = TRUE,
                       filename = "C:/Users/lfaure7/Desktop/COUCHES QGIS/settlements/density_500m.tif",
                       overwrite = TRUE,
                       wopt = list(gdal = c("COMPRESS=ZSTD","PREDICTOR=3","ZSTD_LEVEL=9",
                                            "TILED=YES","NUM_THREADS=ALL_CPUS","BIGTIFF=YES")))

plot(density_local)
plot(alpine_area, add = TRUE, border = "white", lwd = 0.5)
quantile(values(density_local, mat = FALSE),
         probs = c(.5, .75, .9, .95, .99, 1), na.rm = TRUE)