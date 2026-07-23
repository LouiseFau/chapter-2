#' ---
#' title: "Preparing Raster Data"
#' author: "Louise Faure"
#' date: 28.05.2026
#' details: ce script s'applique à la couche de bâtiment du jeu de données 
#' Overture. Pour le calcul des densité issues du jeu de données Human Footprint 
#' Settlement 2019, se référer à une version précédente du code. Nous n'utilisons
#' plus le jeu de données Human Footprint Settlement 2019 car des pixels identifié
#' bâtiments sont présents dans des zones sans bâtiments par comparaison à OSM.
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
#' 
#' À partir du raster binaire corrigé du bâti à 50 m, nous avons d’abord utilisé 
#' GDAL pour agréger les cellules à 100 m par moyenne, de façon à obtenir pour 
#' chaque cellule de 100 m une proportion locale de bâti comprise entre 0 et 1. 
#' Nous avons ensuite importé ce raster dans GRASS GIS, défini la région de 
#' calcul avec g.region, puis appliqué r.neighbors avec method=average et un 
#' fichier de pondération circulaire de taille 13 × 13 cellules. Cette fenêtre 
#' représente approximativement un disque de 1 km² autour de chaque cellule 
#' focale. Le résultat exporté avec r.out.gdal est un raster à 100 m indiquant,
#'  pour chaque cellule, la proportion moyenne de bâti dans son voisinage 
#' d’environ 1 km².
################################################################################


################################################################################
#' ### Step 2 : Densité de bâtiments par fenêtre glissante de 1 km²
#' 
#' La densité de bâtiments a été calculée indépendamment de la taille des 
#' polygones bâtis afin que chaque bâtiment contribue de manière équivalente à 
#' la mesure d’anthropisation. Pour cela, chaque polygone de bâtiment issu 
#' d’Overture Maps a d’abord été converti en un point représentatif situé à 
#' l’intérieur de son emprise. Ces points ont ensuite été rasterisés sur une 
#' grille de 100 m en utilisant une opération additive, de sorte que chaque 
#' cellule contienne le nombre de bâtiments dont le point représentatif y tombe.
#'  Le raster de comptage a ensuite été importé dans GRASS GIS, où une somme 
#'  glissante a été calculée avec `r.neighbors` à l’aide d’une fenêtre circulaire 
#'  représentant approximativement 1 km². Afin de corriger les effets de bord du
#'   masque alpin, une seconde somme glissante a été appliquée à un raster 
#'   binaire indiquant les cellules valides à l’intérieur des Alpes. La densité 
#'   finale a été obtenue en divisant le nombre de bâtiments présents dans la 
#'   fenêtre par la surface alpine effectivement disponible dans cette même 
#'   fenêtre, exprimée en km².

################################################################################
