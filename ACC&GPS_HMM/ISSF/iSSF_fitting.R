#'-------------------------------------------------------------------------------
#' Title: iSSF fitting in Switzerland and historical counterfactual maps ----
#' Author: Louise Faure
#' Date: 06.10.2026
#' **Purpose:** fit the landing iSSF within Switzerland using the RegBL 2026
#' built-cell density (25 m sub-cells aggregated at 100 m, ~1 km2 window), then
#' predict population-level relative selection strength (RSS) for ten historical
#' settlement scenarios (1918-2026), holding topography and land cover constant.
#' Based on E. Nourani 05_INLA_prediction_map.R
#' **Steps:**
#' (1) annotate used/available locations with the RegBL 2026 density, keep strata
#'     entirely inside Switzerland, refit the best iSSF;
#' (2) build Swiss covariate grids (elevation, ruggedness, open vegetation);
#' (3) compute population-level log-RSS maps for each date and RSS relative to a
#'     no-building scenario (habitat terms only, movement terms excluded);
#' (4) summarise the effective area lost per date, with uncertainty.
#' ------------------------------------------------------------------------------

library(tidyverse)
library(terra)
library(glmmTMB)

# Paths and parameters ----
raster_dir <- "/Users/louisefaure/Desktop/dossier sans titre/Rasters"
settlement_dir <- file.path(raster_dir,"historical_settlements_rasters")
output_dir_ch <- "/Users/louisefaure/Desktop/dossier sans titre/iSSF_switzerland"
dir.create(output_dir_ch,recursive=TRUE,showWarnings=FALSE)

annotated_file <- "/Users/louisefaure/Library/CloudStorage/OneDrive-Personnel/THESE/CHAPITRE 2/git/chapter-2/ACC&GPS_HMM/Results/Intermediate_dataset/issf_generated_observed_location_annotated(2).rds"
coordinate_columns <- c("destination_x_m","destination_y_m") # step end points (used and available), EPSG:3035
coordinate_crs <- "EPSG:3035"

dates <- c(1918,1945,1960,1970,1980,1990,2000,2010,2015,2026)
settlement_files <- file.path(settlement_dir,sprintf("built_25m_share_1km2_100m_%d.tif",dates))

minimum_landings <- 30L
fixed_stratum_sd <- 1e3
n_draws <- 1000L
model_control <- glmmTMBControl(optCtrl=list(iter.max=10000,eval.max=10000))
covariates <- c("settlement_density","elevation_100m","ruggedness_100m","prop_low_vegetation_5cells")
write_options <- list(datatype="FLT4S",gdal=c("COMPRESS=ZSTD","PREDICTOR=3","TILED=YES"))
set.seed(500)


#------------------------------------------------------------------------------- STEP 1: refit the iSSF in Switzerland ----
settlement_2026 <- rast(settlement_files[dates==2026])
annotated_data <- readRDS(annotated_file)

# 1.1 Replace the Overture settlement density by the RegBL 2026 density (NA outside Switzerland) ----
locations <- project(vect(as.matrix(as.data.frame(annotated_data)[,coordinate_columns]),crs=coordinate_crs),crs(settlement_2026))
annotated_data$settlement_density <- terra::extract(settlement_2026,locations,ID=FALSE)[[1]]

# 1.2 Keep complete choice sets entirely inside Switzerland ----
n_strata_all <- n_distinct(paste(annotated_data$individual.local.identifier,annotated_data$stratum))

data_ch <- annotated_data %>%
  mutate(used=as.integer(used),individual.local.identifier=as.character(individual.local.identifier),stratum=as.character(stratum)) %>%
  filter(!is.na(individual.local.identifier),!is.na(stratum),used %in% c(0L,1L)) %>%
  group_by(individual.local.identifier,stratum) %>%
  filter(all(!is.na(settlement_density))) %>% # every location of the stratum lies in Switzerland
  ungroup() %>%
  filter(step_length_km>0,if_all(all_of(c("step_length_km","turning_angle_rad",covariates)),~!is.na(.x) & is.finite(.x))) %>%
  group_by(individual.local.identifier,stratum) %>%
  filter(sum(used==1L)==1L,sum(used==0L)>=1L) %>%
  ungroup()

retained_individuals <- data_ch %>% filter(used==1L) %>% count(individual.local.identifier) %>%
  filter(n>=minimum_landings) %>% pull(individual.local.identifier)

data_model_ch <- data_ch %>%
  filter(individual.local.identifier %in% retained_individuals) %>%
  mutate(animal_ID=factor(individual.local.identifier),
         stratum_ID=interaction(individual.local.identifier,stratum,drop=TRUE,lex.order=TRUE),
         log_step_length_km=log(step_length_km),cos_turning_angle=cos(turning_angle_rad))

cat("Strata:",n_strata_all,"in total |",n_distinct(data_model_ch$stratum_ID),"retained in Switzerland | individuals:",
    n_distinct(data_model_ch$animal_ID),"\n")

# 1.3 Standardize covariates (parameters reused for all historical scenarios) ----
standardization <- tibble(variable=covariates,
                          center=map_dbl(covariates,~mean(data_model_ch[[.x]])),
                          scale=map_dbl(covariates,~sd(data_model_ch[[.x]])))
for(i in seq_len(nrow(standardization))){
  data_model_ch[[paste0(standardization$variable[i],"_z")]] <- (data_model_ch[[standardization$variable[i]]]-standardization$center[i])/standardization$scale[i]}

# 1.4 Fit the best iSSF (Poisson formulation, fixed stratum variance) ----
formula_best <- used ~ -1 + settlement_density_z + elevation_100m_z + prop_low_vegetation_5cells_z + ruggedness_100m_z +
  step_length_km + log_step_length_km + cos_turning_angle + (1|stratum_ID) + (0+settlement_density_z|animal_ID)

model_ch <- glmmTMB(formula_best,family=poisson(link="log"),data=data_model_ch,
                    map=list(theta=factor(c(NA,1L))),start=list(theta=c(log(fixed_stratum_sd),0)),control=model_control)
print(summary(model_ch))
saveRDS(list(model=model_ch,standardization=standardization),file.path(output_dir_ch,"issf_switzerland_regbl2026.rds"))

beta <- fixef(model_ch)$cond
beta_se <- sqrt(diag(vcov(model_ch)$cond))


#------------------------------------------------------------------------------- STEP 2: Swiss covariate grids ----
grid_ch <- settlement_2026 # Swiss 100 m grid, NA outside Switzerland

align_to_grid <- function(x,method){
  x <- crop(x,grid_ch,snap="out")
  if(!compareGeom(x,grid_ch,stopOnError=FALSE)) x <- resample(x,grid_ch,method=method)
  mask(x,grid_ch)}

elevation_ch <- align_to_grid(rast(file.path(raster_dir,"elevation_100m.tif")),"bilinear")
ruggedness_ch <- align_to_grid(rast(file.path(raster_dir,"ruggedness_100m.tif")),"bilinear")

# open vegetation computed on a buffered extent so that border cells keep their neighbours
landcover <- crop(rast(file.path(raster_dir,"landcover_100m.tif")),ext(grid_ch)+500,snap="out")
low_vegetation <- ifel(landcover==5 | landcover==6 | landcover==7,1,0)
rook_window <- matrix(c(NA,1,NA,1,1,1,NA,1,NA),nrow=3,byrow=TRUE)
prop_low_vegetation_ch <- align_to_grid(focal(low_vegetation,w=rook_window,fun="mean",na.rm=FALSE),"near")


#------------------------------------------------------------------------------- STEP 3: log-RSS maps for each date ----
z <- function(x,v){s <- standardization[standardization$variable==v,]; (x-s$center)/s$scale}

# habitat part of the linear predictor (constant across dates; movement terms excluded)
eta_habitat <- beta[["elevation_100m_z"]]*z(elevation_ch,"elevation_100m") +
  beta[["prop_low_vegetation_5cells_z"]]*z(prop_low_vegetation_ch,"prop_low_vegetation_5cells") +
  beta[["ruggedness_100m_z"]]*z(ruggedness_ch,"ruggedness_100m")

settlement_stack <- rast(settlement_files)
names(settlement_stack) <- paste0("y",dates)
z_settlement <- mask(z(settlement_stack,"settlement_density"),eta_habitat)
z_zero <- z(0,"settlement_density")
b_settlement <- beta[["settlement_density_z"]]

eta_dates <- eta_habitat + b_settlement*z_settlement             # population-level log-RSS
rss_vs_zero <- exp(b_settlement*(z_settlement-z_zero))           # 1 = no effect of buildings
names(eta_dates) <- names(rss_vs_zero) <- paste0("y",dates)

writeRaster(eta_dates,file.path(output_dir_ch,"issf_log_rss_by_date.tif"),overwrite=TRUE,wopt=write_options)
writeRaster(rss_vs_zero,file.path(output_dir_ch,"issf_rss_vs_no_buildings_by_date.tif"),overwrite=TRUE,wopt=write_options)


#------------------------------------------------------------------------------- STEP 4: effective area lost per date ----
cell_area_km2 <- prod(res(grid_ch))/1e6
w_zero <- global(exp(eta_habitat+b_settlement*z_zero),"sum",na.rm=TRUE)[1,1]
b_draws <- rnorm(n_draws,b_settlement,beta_se[["settlement_density_z"]])

summary_dates <- map_dfr(seq_along(dates),function(i){
  d <- values(z_settlement[[i]]-z_zero,mat=FALSE)
  tab <- tibble(d=round(d[!is.na(d)],3)) %>% count(d) # binning makes the uncertainty draws fast
  area_lost <- function(b) sum(tab$n*(1-exp(b*tab$d)))*cell_area_km2
  lost_draws <- map_dbl(b_draws,area_lost)
  tibble(date=dates[i],
         area_total_km2=sum(tab$n)*cell_area_km2,
         effective_area_lost_km2=area_lost(b_settlement),
         lost_ci_low=quantile(lost_draws,0.025),
         lost_ci_high=quantile(lost_draws,0.975),
         selection_weighted_loss_pct=100*(1-global(exp(eta_dates[[i]]),"sum",na.rm=TRUE)[1,1]/w_zero))}) %>%
  mutate(lost_since_1918_km2=effective_area_lost_km2-first(effective_area_lost_km2))

print(summary_dates)
write_csv(summary_dates,file.path(output_dir_ch,"effective_area_lost_by_date.csv"))
