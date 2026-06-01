#' ## Buildings from Overture Maps - distance pipeline
# First tests by Louise on chapter 2
# Packages requis (à installer si manquants) : duckdb, DBI

library(terra)
library(sf)
library(duckdb)
library(DBI)

# Réutilise alpine_area déjà chargé plus haut dans le script
alpine_area <- vect("C:/Users/lfaure7/Desktop/COUCHES QGIS/alpine_area/alpine_area_corrected_geometry.shp")
out_path <- vect("C:/Users/lfaure7/Desktop/COUCHES QGIS/settlements/overture_buildings_alps.gpkg")

# === 1. Régler la topologie de alpine_area ===
# Ton shapefile a un avertissement d'auto-intersection ; ça peut faire
# échouer les filtres spatiaux. On nettoie d'abord.
alpine_valid <- makeValid(alpine_area)
alpine_wgs84 <- project(alpine_valid, "EPSG:4326")

# === 2. Lecture FILTRÉE du GPKG ===
# terra::vect a un argument `filter` : ne charge en mémoire QUE les polygones
# intersectant le filtre. Utilise l'index spatial GDAL du GPKG -> rapide.
# buildings_alps_wgs84 <- vect(out_path, filter = alpine_wgs84)
length(buildings_alps_wgs84)   # tu devrais voir un nombre sensiblement plus petit

# === 3. Reprojection en EPSG:3035 ===
buildings_3035 <- project(buildings_alps_wgs84, "EPSG:3035")

# Optionnel : sauvegarde du subset propre pour ne pas refaire le filtre
writeVector(buildings_3035,
            "C:/Users/lfaure7/Desktop/COUCHES QGIS/settlements/overture_buildings_alps_3035.gpkg",
            overwrite = TRUE)

# === 4. Rasterisation 30 m, alignée sur alpine_area ===
template <- rast(ext(alpine_valid), resolution = 30, crs = "EPSG:3035")

buildings_rast <- rasterize(
  buildings_3035, template,
  field      = 1,
  background = NA,
  touches    = TRUE,                      # indispensable : la plupart des bâtiments < 30 m
  filename   = "C:/Users/lfaure7/Desktop/COUCHES QGIS/settlements/buildings_overture_30m.tif",
  overwrite  = TRUE,
  wopt = list(gdal = c("COMPRESS=ZSTD","PREDICTOR=2","ZSTD_LEVEL=9",
                       "TILED=YES","NUM_THREADS=ALL_CPUS","BIGTIFF=YES"))
)

# === 5. Distance (transformée euclidienne, rapide) ===
rdist <- distance(buildings_rast,
                  filename = "C:/Users/lfaure7/Desktop/COUCHES QGIS/settlements/dist_to_overture.tif",
                  overwrite = TRUE)
rdist_clip <- crop(rdist, alpine_valid, mask = TRUE)
plot(rdist_clip)
quantile(values(rdist_clip, mat = FALSE),
         probs = c(.5, .75, .9, .95, .99, 1), na.rm = TRUE)