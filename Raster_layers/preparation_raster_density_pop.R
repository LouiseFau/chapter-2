library(terra)
library(sf)

# --- chargement ---
r1  <- terra::rast("C:/Users/lfaure7/Downloads/GHS_POP_E2020_GLOBE_R2023A_54009_100_V1_0_R4_C20/GHS_POP_E2020_GLOBE_R2023A_54009_100_V1_0_R4_C20.tif")
r2  <- terra::rast("C:/Users/lfaure7/Downloads/GHS_POP_E2020_GLOBE_R2023A_54009_100_V1_0_R4_C19/GHS_POP_E2020_GLOBE_R2023A_54009_100_V1_0_R4_C19.tif")
alpes <- sf::st_read("C:/Users/lfaure7/Desktop/COUCHES QGIS/alpine_area/alpine_area_corrected_geometry.shp")

# --- reprojeter le polygone dans le CRS des rasters si nécessaire ---
alpes_v <- terra::vect(alpes)
alpes_v <- terra::project(alpes_v, terra::crs(r1))

# --- fusionner les deux plaques ---
# merge() gère les zones qui se chevauchent (prend r1 en priorité)
# mosaic() avec fun="mean" si tu veux moyenner les overlaps
r_merged <- terra::merge(r1, r2)

# --- crop (bounding box) PUIS mask (forme exacte du polygone) ---
# mask() met les cellules hors polygone à NA, pas à 0
r_cropped <- terra::crop(r_merged, alpes_v)
r_masked  <- terra::mask(r_cropped, alpes_v)

# vérification : les cellules hors Alpes doivent être NA
# les cellules dans les Alpes avec population = 0 restent à 0
freq_check <- terra::freq(r_masked, value = 0)
print(freq_check)  # nombre de cellules à 0 (vraie valeur)

# --- export ---
terra::writeRaster(r_masked,
                   filename  = "C:/Users/lfaure7/Desktop/COUCHES QGIS/density-resident-pop-km2/density_pop_km2.tif",
                   datatype  = "INT4U",       # entier non signé 32 bits, adapté aux données pop
                   NAflag    = -9999,          # valeur sentinelle pour NA dans le fichier
                   overwrite = TRUE)


# calculate human density per km2
r_merged <- terra::merge(r1, r2)

# calculer la densité focale sur le raster COMPLET
w <- matrix(1, nrow = 11, ncol = 11) 
r_density_focal <- terra::focal(r_merged, w = w, fun = "sum", na.rm = TRUE)

# crop puis mask APRÈS le focal
r_density_cropped <- terra::crop(r_density_focal, alpes_v)
r_density_masked  <- terra::mask(r_density_cropped, alpes_v)

# export
terra::writeRaster(r_density_masked,
                   filename  = "C:/Users/lfaure7/Desktop/COUCHES QGIS/density-resident-pop-km2/density_pop_km2.tif",
                   datatype  = "FLT4S",   # float, pas INT4U : les valeurs focales ne sont plus entières
                   NAflag    = -9999,
                   overwrite = TRUE)




# control
vals <- terra::values(r_density_masked)
vals_nonzero <- vals[!is.na(vals) & vals > 0]

hist(log(vals_nonzero), breaks = 100,
     main = "Distribution log(pop) dans les Alpes",
     xlab = "log(population)")


# Plotting the CDF
graph_cdf <- plot(ecdf(log1p(vals_nonzero)),
     main = "CDF of log(population + 1)",
     xlab = "log(population + 1)",
     ylab = "Cumulative proportion",
     pch  = NA, col = "steelblue", lwd = 2)

# add reference lines for readable quantiles
qs <- quantile(vals_nonzero, c(0.50, 0.75, 0.90, 0.95))
abline(v = log1p(qs), col = "grey60", lty = 2)
text(log1p(qs), 0.05,
     labels = paste0(c(50,75,90,95), "%\n(", round(qs,1), ")"),
     cex = 0.75, adj = 0)
