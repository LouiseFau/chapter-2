#' title: "Distance from building from Overture Maps"
#' author: "Louise Faure"
#' date: 01.06.2026
#' details: (i) filter buildings

# packages
library(terra)
library(sf)
library(duckdb)
library(DBI)



# === Chemins ===
gpkg_in     <- "C:/Users/lfaure7/Desktop/COUCHES QGIS/settlements/overture_buildings_alps.gpkg"
gpkg_out    <- "C:/Users/lfaure7/Desktop/COUCHES QGIS/settlements/overture_buildings_alps_3035.gpkg"
alpine_path <- "C:/Users/lfaure7/Desktop/COUCHES QGIS/alpine_area/alpine_area_corrected_geometry.shp"

# === 1. Polygone alpin : validé et exporté en WGS84 pour DuckDB ===
alpine_area  <- vect(alpine_path)
alpine_valid <- makeValid(alpine_area)
alpine_wgs84 <- project(alpine_valid, "EPSG:4326")

tmp_alps <- tempfile(fileext = ".gpkg")
writeVector(alpine_wgs84, tmp_alps, overwrite = TRUE)

# === 2. DuckDB : filtre polygonal + reprojection, en une seule passe ===
con <- dbConnect(duckdb::duckdb())
dbExecute(con, "INSTALL spatial; LOAD spatial;")

query <- sprintf("
COPY (
  SELECT
    b.id,
    ST_Transform(b.geom, 'OGC:CRS84', 'EPSG:3035') AS geometry
  FROM ST_Read('%s') b, ST_Read('%s') a
  WHERE ST_Intersects(b.geom, ST_SetCRS(a.geom, 'OGC:CRS84'))
) TO '%s' WITH (FORMAT GDAL, DRIVER 'GPKG');
", gpkg_in, tmp_alps, gpkg_out)

system.time(dbExecute(con, query))
dbDisconnect(con, shutdown = TRUE)

# === 3. Chargement du GPKG filtré (petit, et déjà en EPSG:3035) ===
buildings_3035 <- vect(gpkg_out)
nrow(buildings_3035)                          # bien plus petit que la version bbox
crs(buildings_3035, describe = TRUE)$code     # "3035"

# === 4. Rasterisation 30 m, alignée sur alpine_area ===
template <- rast(ext(alpine_valid), resolution = 30, crs = "EPSG:3035")

buildings_rast <- rasterize(
  buildings_3035, template,
  field      = 1,
  background = NA,
  touches    = TRUE,
  filename   = "C:/Users/lfaure7/Desktop/COUCHES QGIS/settlements/buildings_overture_30m.tif",
  overwrite  = TRUE,
  wopt = list(gdal = c("COMPRESS=ZSTD","PREDICTOR=2","ZSTD_LEVEL=9",
                       "TILED=YES","NUM_THREADS=ALL_CPUS","BIGTIFF=YES"))
)

# === 5. Distance (EDT, rapide) ===
rdist <- distance(
  buildings_rast,
  filename = "C:/Users/lfaure7/Desktop/COUCHES QGIS/settlements/dist_to_overture.tif",
  overwrite = TRUE
)
rdist_clip <- crop(rdist, alpine_valid, mask = TRUE)
plot(rdist_clip)
quantile(values(rdist_clip, mat = FALSE),
         probs = c(.5, .75, .9, .95, .99, 1), na.rm = TRUE)