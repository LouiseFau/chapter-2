#'-------------------------------------------------------------------------------
#' Title: iSSF fitting and counterfactual map prep. ----
#' Authors : Louise Faure
#' Date : 04.08.26
#' 
#' **Info:** this script follow the Data_processing_annotation.R script and 
#' Covariates selection.R This script is based on E. Nourani 
#' 05_INLA_prediction_map.R 
#' (https://github.com/kamransafi/GoldenEagles/blob/main/WP3_Soaring_Ontogeny/MS2_flight_landscape/05_INLA_prediction_map.R)
#' 
#' **Purpose:** build a conterfactual scenario where settlement density is set to 
#' observed q05 values to highlight changes in relative selection strenght 
#' associated to human buildings. 
#' 
#' **Steps:**
#' (1) align raster and prepare prediction grid
#' (2) create a new dataset for prediction with  the values of the settlement 
#' density settled to 0 (=q5). 
#' (3) split the dataset into chunck that will be sent to raven. 

#' (3) control coefficient match between predicted and observed dataset
#' (4) use the models on the two dataset, extract probability of landing per cell
#' using predict() function
#' (3) export two rasters 
#' (4) prepare banane buffers
#' ------------------------------------------------------------------------------

library(tidyverse)
library(terra)
library(sf)
library(mapview)
library(glmmTMB)

# Topographic and human layers ----
settlement_density_60 <- terra::rast("/Users/louisefaure/Desktop/dossier sans titre/Rasters/settlement_density_1km2_100m.tif")
elevation_100m_60 <- terra::rast("/Users/louisefaure/Desktop/dossier sans titre/Rasters/elevation_100m.tif")
ruggedness_100m_60 <- terra::rast("/Users/louisefaure/Desktop/dossier sans titre/Rasters/ruggedness_100m.tif")
landcover_100m_60 <- terra::rast("/Users/louisefaure/Desktop/dossier sans titre/Rasters/landcover_100m.tif")

names(settlement_density_60) <- "settlement_density"
names(elevation_100m_60) <- "elevation_100m"
names(ruggedness_100m_60) <- "ruggedness_100m"
names(landcover_100m_60) <- "landcover_100m"

#-------------------------------------------------------------------------------STEP 1: align raster layers and prepare the prediction grid ----

# 1.1 Define the common extent and reference grid ----
common_extent_60 <- Reduce(terra::intersect,list(
  terra::ext(settlement_density_60),
  terra::ext(elevation_100m_60),
  terra::ext(ruggedness_100m_60),
  terra::ext(landcover_100m_60)
))

prediction_template_60 <- terra::crop(landcover_100m_60,common_extent_60,snap = "in")
names(prediction_template_60) <- "landcover_100m"

# 1.2 Align continuous rasters with the reference grid ----
settlement_density_60 <- terra::resample(terra::crop(settlement_density_60,common_extent_60,snap = "out"),prediction_template_60,method = "bilinear")
elevation_100m_60 <- terra::resample(terra::crop(elevation_100m_60,common_extent_60,snap = "out"),prediction_template_60,method = "bilinear")
ruggedness_100m_60 <- terra::resample(terra::crop(ruggedness_100m_60,common_extent_60,snap = "out"),prediction_template_60,method = "bilinear")
landcover_100m_60 <- prediction_template_60

names(settlement_density_60) <- "settlement_density"
names(elevation_100m_60) <- "elevation_100m"
names(ruggedness_100m_60) <- "ruggedness_100m"
names(landcover_100m_60) <- "landcover_100m"

# 1.3 Reconstitute the proportion of low-vegetation cells ----
low_vegetation_binary_60 <- terra::ifel(
  landcover_100m_60 == 5 |
    landcover_100m_60 == 6 |
    landcover_100m_60 == 7,
  1,
  0)

names(low_vegetation_binary_60) <- "low_vegetation_binary"

rook_window_5cells_60 <- matrix(c(
  NA,1,NA,
  1,1,1,
  NA,1,NA
),nrow = 3,byrow = TRUE)

prop_low_vegetation_5cells_60 <- terra::focal(low_vegetation_binary_60,w = rook_window_5cells_60,fun = "mean",na.rm = FALSE)
names(prop_low_vegetation_5cells_60) <- "prop_low_vegetation_5cells"

# 1.4 Combine prediction covariates ----
stck_covariates_60 <- c(settlement_density_60,elevation_100m_60,ruggedness_100m_60,prop_low_vegetation_5cells_60)
names(stck_covariates_60) <- c("settlement_density","elevation_100m","ruggedness_100m","prop_low_vegetation_5cells")

# 1.5 Prepare prediction data in spatial-table format ----
alps_df_60 <- terra::as.data.frame(stck_covariates_60,xy = TRUE,cells = TRUE,na.rm = FALSE) %>%
  dplyr::rename(x_3035 = x,y_3035 = y)

prediction_covariates_60 <- c("settlement_density","elevation_100m","ruggedness_100m","prop_low_vegetation_5cells")

alps_df_no_na_60 <- alps_df_60 %>%
  tidyr::drop_na(dplyr::all_of(prediction_covariates_60)) %>%
  as.data.frame()

saveRDS(alps_df_no_na_60, file = "/Users/louisefaure/Desktop/dossier sans titre/iSSF_intermediate_results/alps_grid.rds")



#------------------------------------------------------------------------------- STEP 2: create new datasets for prediction ----
all_data_60 <- readRDS("/Users/louisefaure/Library/CloudStorage/OneDrive-Personnel/THESE/CHAPITRE 2/git/chapter-2/ACC&GPS_HMM/Results/Intermediate_dataset/issf_model_data_60.rds")
alps_df_no_na_60 <- readRDS("/Users/louisefaure/Desktop/dossier sans titre/iSSF_intermediate_results/alps_grid.rds")

prediction_covariates_60 <- c("settlement_density","elevation_100m","ruggedness_100m","prop_low_vegetation_5cells")
prediction_covariates_z_60 <- paste0(prediction_covariates_60,"_z")

set.seed(500)
n_prediction_cells_60 <- nrow(alps_df_no_na_60)

# 2.1 Extract standardization parameters from the training dataset ----
prediction_standardization_60 <- tibble::tibble(
  variable = prediction_covariates_60,
  center = vapply(prediction_covariates_60,\(x) mean(all_data_60[[x]],na.rm = TRUE),numeric(1)),
  scale = vapply(prediction_covariates_60,\(x) stats::sd(all_data_60[[x]],na.rm = TRUE),numeric(1)))

# 2.2 Sample rows from the training dataset ----
prediction_seed_rows_60 <- all_data_60 %>%
  dplyr::group_by(stratum_ID) %>%
  dplyr::slice_sample(n = 1) %>%
  dplyr::ungroup() %>%
  dplyr::slice_sample(n = n_prediction_cells_60,replace = TRUE) %>%
  dplyr::select(-dplyr::any_of(c(prediction_covariates_60,prediction_covariates_z_60)))

# 2.3 Add the raster covariates ----
alps_data_observed_60 <- prediction_seed_rows_60 %>%
  dplyr::bind_cols(alps_df_no_na_60) %>%
  dplyr::mutate(used = NA_integer_,scenario = "observed") %>%
  as.data.frame()

# 2.4 Standardize raster covariates using the training dataset ----
for(i in seq_len(nrow(prediction_standardization_60))) {
  variable_60 <- prediction_standardization_60$variable[[i]]
  center_60 <- prediction_standardization_60$center[[i]]
  scale_60 <- prediction_standardization_60$scale[[i]]
  alps_data_observed_60[[paste0(variable_60,"_z")]] <- (alps_data_observed_60[[variable_60]] - center_60) / scale_60}

# 2.5 Calculate settlement-density Q05 ----
settlement_q05_raw_60 <- as.numeric(stats::quantile(all_data_60$settlement_density[all_data_60$used == 0L],probs = 0.05,na.rm = TRUE,names = FALSE))

settlement_standardization_60 <- prediction_standardization_60 %>% dplyr::filter(variable == "settlement_density")

settlement_q05_z_60 <- (settlement_q05_raw_60 - settlement_standardization_60$center) / settlement_standardization_60$scale

#------------------------------------------------------------------------------- STEP 3: prepare cluster input files ----

cluster_input_directory_60 <- "/Users/louisefaure/Desktop/dossier sans titre/iSSF_cluster/input"
cluster_chunk_directory_60 <- file.path(cluster_input_directory_60,"chunks")
dir.create(cluster_chunk_directory_60,recursive = TRUE,showWarnings = FALSE)

prediction_columns_60 <- c(
  "cell","stratum_ID","animal_ID",
  "step_length_km","log_step_length_km","cos_turning_angle",
  "settlement_density_z","elevation_100m_z",
  "prop_low_vegetation_5cells_z","ruggedness_100m_z"
)

alps_prediction_compact_60 <- alps_data_observed_60 %>%
  dplyr::select(dplyr::all_of(prediction_columns_60)) %>%
  as.data.frame()

chunk_size_60 <- 250000L
n_chunks_60 <- ceiling(nrow(alps_prediction_compact_60) / chunk_size_60)

for(i in seq_len(n_chunks_60)) {
  first_row_60 <- (i - 1L) * chunk_size_60 + 1L
  last_row_60 <- min(i * chunk_size_60,nrow(alps_prediction_compact_60))
  
  saveRDS(
    alps_prediction_compact_60[first_row_60:last_row_60,,drop = FALSE],
    file.path(cluster_chunk_directory_60,sprintf("prediction_chunk_%04d.rds",i)),
    compress = "gzip"
  )
}

counterfactual_parameters_60 <- list(
  settlement_q05_raw = settlement_q05_raw_60,
  settlement_q05_z = settlement_q05_z_60,
  n_prediction_cells = nrow(alps_prediction_compact_60),
  n_chunks = n_chunks_60,
  chunk_size = chunk_size_60
)

saveRDS(
  counterfactual_parameters_60,
  file.path(cluster_input_directory_60,"counterfactual_parameters_60.rds"),
  compress = "gzip"
)

file.copy(
  "/Users/louisefaure/Library/CloudStorage/OneDrive-Personnel/THESE/CHAPITRE 2/git/chapter-2/ACC&GPS_HMM/Results/Intermediate_dataset/issf_model_data_60.rds",
  file.path(cluster_input_directory_60,"issf_model_data_60.rds"),
  overwrite = TRUE
)

terra::writeRaster(
  prediction_template_60,
  file.path(cluster_input_directory_60,"prediction_template_60.tif"),
  overwrite = TRUE,
  wopt = list(gdal = c("COMPRESS=DEFLATE","BIGTIFF=YES"))
)

writeLines(
  as.character(n_chunks_60),
  file.path(cluster_input_directory_60,"n_chunks.txt")
)

print(c(
  n_cells = nrow(alps_prediction_compact_60),
  n_chunks = n_chunks_60,
  rows_per_chunk = chunk_size_60
))





#-------------------------------------------------------------------------------STEP 3: make sure model coefficients match the original model's ----
local_root_60 <- "/Users/louisefaure/Desktop/dossier sans titre/iSSF_cluster"

best_model_raven_60 <- readRDS(
  file.path(
    local_root_60,
    "output/model/best_issf_model_60.rds"
  )
)

summary(best_model_raven_60) # to check whether we have the same coefficient for settlement_density, and it is confirmed !

counterfactual_parameters_60 <- readRDS(
  file.path(
    local_root_60,
    "input/counterfactual_parameters_60.rds"
  )
)

counterfactual_parameters_60$settlement_q05_raw
counterfactual_parameters_60$settlement_q05_z


controle_q05_zero_60 <- isTRUE(
  all.equal(
    as.numeric(counterfactual_parameters_60$settlement_q05_raw),
    0
  )
)

print(controle_q05_zero_60)


# plot result cluster 
library(terra)

local_root_60 <- "/Users/louisefaure/Desktop/dossier sans titre/iSSF_cluster"

final_prediction_raster_60 <- terra::rast(
  file.path(
    local_root_60,
    "output/final/issf_counterfactual_predictions_60.tif"
  )
)

# Retain the two population-level scenarios
comparison_eta_60 <- final_prediction_raster_60[[
  c("eta_observed","eta_settlement_q05")
]]

names(comparison_eta_60) <- c(
  "Observed settlement density",
  "Settlement density fixed at Q05"
)

# Sample values only to define robust common display limits
controle_eta_sample_60 <- terra::spatSample(
  comparison_eta_60,
  size = 500000,
  method = "regular",
  values = TRUE,
  na.rm = TRUE,
  as.df = TRUE
)

common_eta_limits_60 <- as.numeric(
  stats::quantile(
    unlist(controle_eta_sample_60,use.names = FALSE),
    probs = c(0.01,0.99),
    na.rm = TRUE,
    names = FALSE
  )
)

print(common_eta_limits_60)

# Clip only for visualization; the original raster is not modified
comparison_eta_display_60 <- terra::clamp(
  comparison_eta_60,
  lower = common_eta_limits_60[1],
  upper = common_eta_limits_60[2],
  values = TRUE
)

map_palette_60 <- hcl.colors(
  100,
  palette = "Inferno",
  rev = FALSE
)

map_output_60 <- file.path(
  local_root_60,
  "output/final/issf_observed_vs_q05.png"
)

png(
  filename = map_output_60,
  width = 3600,
  height = 1800,
  res = 300
)

graphics::par(
  mfrow = c(1,2),
  mar = c(2.5,2.5,4,5)
)

terra::plot(
  comparison_eta_display_60[["Observed settlement density"]],
  col = map_palette_60,
  range = common_eta_limits_60,
  main = "Observed settlement density",
  axes = FALSE,
  maxcell = 1000000,
  plg = list(
    title = "Log relative\nselection strength"
  )
)

terra::plot(
  comparison_eta_display_60[["Settlement density fixed at Q05"]],
  col = map_palette_60,
  range = common_eta_limits_60,
  main = if(controle_q05_zero_60) {
    "Settlement density set to zero"
  } else {
    "Settlement density fixed at Q05"
  },
  axes = FALSE,
  maxcell = 1000000,
  plg = list(
    title = "Log relative\nselection strength"
  )
)

graphics::mtext(
  "Population-level relative selection strength for landing sites",
  outer = TRUE,
  line = -1.5,
  cex = 1.1
)

dev.off()

print(map_output_60)



# test 
final_prediction_raster_60 <- terra::rast('/Users/louisefaure/Desktop/dossier sans titre/iSSF_cluster/output/final/issf_counterfactual_predictions_60.tif')
reduction_rss_percent_60 <-
  100 * (
    1 -
      final_prediction_raster_60[["rss_observed_vs_q05"]]
  )

names(reduction_rss_percent_60) <-
  "reduction_rss_percent"
reduction_breaks_60 <- c(
  0, 1, 5, 10, 25, 50, 75, 90, 100
)

reduction_colors_60 <- hcl.colors(
  length(reduction_breaks_60) - 1,
  palette = "YlOrRd"
)

terra::plot(
  reduction_rss_percent_60,
  breaks = reduction_breaks_60,
  col = reduction_colors_60,
  axes = FALSE,
  main = paste0(
    "Relative reduction in landing-site selection\n",
    "associated with observed settlement density"
  ),
  plg = list(
    title = "Reduction in RSS (%)",
    shrink = 0.85
  ),
  maxcell = 2000000
)

output_reduction_raster_60 <- "/Users/louisefaure/Desktop/dossier sans titre/iSSF_cluster/output/final/reduction_rss_percent_60.tif"

terra::writeRaster(
  reduction_rss_percent_60,
  filename = output_reduction_raster_60,
  overwrite = TRUE,
  datatype = "FLT4S",
  NAflag = -9999,
  gdal = c(
    "COMPRESS=DEFLATE",
    "PREDICTOR=3",
    "BIGTIFF=YES"
  )
)

print(output_reduction_raster_60)
