#'-------------------------------------------------------------------------------
#' Title: iSSF on acc and gps classified behaviors ----
#' Author: Louise Faure
#' Date: 09.10.2026
#' **Purpose:** this script follow Data_processing_annotation.R. The script
#' compare the acc and gps classified dataset, and obtain behavior specific
#' coefficients of relative selection strength.
#' **Steps:**
#' (1) comparison of acc and gps classified dataset and obtention of a common plot
#' with the selection coefficient for each variables with their interval of
#' confidence (visualisation n°1)
#' (2) prediction at the alpine scale with the generic acc-based model: raster of
#' the change in relative selection strength between the observed settlement
#' density (Overture) and a settlement density set to 0.
#' (3) table of the settlement density selection coefficients for the GPS model,
#' the generic ACC model, and resting and feeding from the interaction model
#' (visualisation n°2).
#' ------------------------------------------------------------------------------


# libraries
library(tidyverse)
library(terra)
library(glmmTMB)
library(gt)

# Paths
intermediate_dataset_directory <- "/Users/louisefaure/Library/CloudStorage/OneDrive-Personnel/THESE/CHAPITRE 2/git/chapter-2/ACC&GPS_HMM/Results/Intermediate_dataset"
issf_directory <- "/Users/louisefaure/Library/CloudStorage/OneDrive-Personnel/THESE/CHAPITRE 2/git/chapter-2/ACC&GPS_HMM/ISSF"
gps_model_file <- file.path(issf_directory,"issf_best_model_60.rds")                                         # exported by iSSF_covariate_selection.R
acc_annotated_file <- file.path(intermediate_dataset_directory,"issf_generated_observed_location_annotated_acc.rds")
overture_file <- "/Users/louisefaure/Library/CloudStorage/OneDrive-Personnel/THESE/COUCHES QGIS/COUCHES QGIS/settlements/overture_built_25m_share_1km2_100m.tif"
alps_perimeter_file <- "/Users/louisefaure/Library/CloudStorage/OneDrive-Personnel/THESE/COUCHES QGIS/COUCHES QGIS/alpine_area/alpine_area_corrected_geometry.shp"
raster_dir <- "/Users/louisefaure/Desktop/dossier sans titre/Rasters"
static_files <- c(elevation_100m=file.path(raster_dir,"elevation_100m.tif"),
                  ruggedness_100m=file.path(raster_dir,"ruggedness_100m.tif"))
landcover_file <- file.path(raster_dir,"landcover_100m.tif")
output_dir_alps <- "/Users/louisefaure/Desktop/dossier sans titre/iSSF_alps"
dir.create(output_dir_alps,recursive=TRUE,showWarnings=FALSE)

# Parameters ----
minimum_landings <- 30L
fixed_stratum_sd <- 1e3
model_control <- glmmTMBControl(optCtrl=list(iter.max=10000,eval.max=10000))
covariates <- c("settlement_density","elevation_100m","ruggedness_100m","prop_low_vegetation_5cells")
habitat_terms <- paste0(covariates,"_z")
movement_terms <- c("step_length_km","log_step_length_km","cos_turning_angle")
low_vegetation_classes <- c(5,6,7)
prediction_chunk_rows <- 500000L  # rows per predict() call (glmmTMB rebuilds its objects for each call: chunks keep memory low)
n_check_cells <- 5000L
write_options <- list(datatype="FLT4S",gdal=c("COMPRESS=ZSTD","PREDICTOR=3","TILED=YES"))
font_family <- "Baskerville"
clr_gps <- "grey25"               # dark grey for the GPS model
clr_acc <- "slateblue2"           # purple for the ACC models
set.seed(500)
stopifnot(file.exists(c(gps_model_file,acc_annotated_file,overture_file,alps_perimeter_file,static_files,landcover_file)))

# Formulas: generic (as selected on the GPS dataset) and with the feeding interaction ----
formula_generic <- used ~ -1 + settlement_density_z + elevation_100m_z + prop_low_vegetation_5cells_z + ruggedness_100m_z +
  step_length_km + log_step_length_km + cos_turning_angle + (1|stratum_ID) + (0+settlement_density_z|animal_ID)
# feeding is 0/1 (constant within a stratum): no main effect, one deviation per habitat term, delta kept fixed
formula_interaction <- update(formula_generic,. ~ . + feeding:settlement_density_z + feeding:elevation_100m_z +
                                feeding:prop_low_vegetation_5cells_z + feeding:ruggedness_100m_z)

fit_issf <- function(model_formula,data){
  glmmTMB(model_formula,family=poisson(link="log"),data=data,
          map=list(theta=factor(c(NA,1L))),start=list(theta=c(log(fixed_stratum_sd),0)),control=model_control)}


#------------------------------------------------------------------------------- STEP 1: comparison of the GPS and ACC models ----
#' **Steps:**
#' (i) load the selected GPS model and its standardization parameters;
#' (ii) prepare the ACC dataset as in iSSF_covariate_selection.R (complete strata,
#'      individuals with at least 30 landings, standardization on the ACC data);
#' (iii) fit the generic ACC model with the GPS formula;
#' (iv) plot the coefficients of both models on one figure (visualisation n°1).

# 1.1 GPS model ----
gps_export <- readRDS(gps_model_file)
model_gps <- gps_export$model
standardization_gps <- gps_export$standardization
data_gps <- model_gps$frame %>%
  mutate(settlement_density=settlement_density_z*standardization_gps$scale[standardization_gps$variable=="settlement_density"]+
           standardization_gps$center[standardization_gps$variable=="settlement_density"])   # raw density, for the common contrast
cat("GPS model:",gps_export$model_name,"| strata:",n_distinct(data_gps$stratum_ID),"| individuals:",n_distinct(data_gps$animal_ID),"\n")

# 1.2 ACC dataset ----
acc_annotated <- readRDS(acc_annotated_file)

data_complete_acc <- acc_annotated %>%
  mutate(used=as.integer(used),individual.local.identifier=as.character(individual.local.identifier),stratum=as.character(stratum)) %>%
  filter(!is.na(individual.local.identifier),!is.na(stratum),used %in% c(0L,1L),step_length_km>0,
         if_all(all_of(c("step_length_km","turning_angle_rad",covariates)),~!is.na(.x) & is.finite(.x))) %>%
  group_by(individual.local.identifier,stratum) %>%
  filter(sum(used==1L)==1L,sum(used==0L)>=1L) %>%
  ungroup()

retained_individuals_acc <- data_complete_acc %>% filter(used==1L) %>% count(individual.local.identifier) %>%
  filter(n>=minimum_landings) %>% pull(individual.local.identifier)

data_acc <- data_complete_acc %>%
  filter(individual.local.identifier %in% retained_individuals_acc) %>%
  mutate(animal_ID=factor(individual.local.identifier),
         stratum_ID=interaction(individual.local.identifier,stratum,drop=TRUE,lex.order=TRUE),
         log_step_length_km=log(step_length_km),cos_turning_angle=cos(turning_angle_rad),
         feeding=as.integer(behavior_detail_end=="feeding"))   # terrestrial behavior following the landing: feeding = 1, resting = 0

# Standardization on the ACC tracking data (centre and scale reused for the alpine predictions)
standardization_acc <- tibble(variable=covariates,
                              center=map_dbl(covariates,~mean(data_acc[[.x]])),
                              scale=map_dbl(covariates,~sd(data_acc[[.x]])))
for(i in seq_len(nrow(standardization_acc))){
  data_acc[[paste0(standardization_acc$variable[i],"_z")]] <- (data_acc[[standardization_acc$variable[i]]]-standardization_acc$center[i])/standardization_acc$scale[i]}

# control: strata and individuals, by terrestrial behavior
control_acc_strata <- data_acc %>% filter(used==1L) %>%
  group_by(behavior_detail_end) %>%
  summarise(n_strata=n(),n_individuals=n_distinct(animal_ID),
            median_strata_per_individual=median(table(animal_ID)[table(animal_ID)>0]),.groups="drop")
print(control_acc_strata)
cat("ACC dataset: strata",n_distinct(data_acc$stratum_ID),"| individuals",n_distinct(data_acc$animal_ID),"\n")

# 1.3 Generic ACC model ----
model_acc <- fit_issf(formula_generic,data_acc)
print(summary(model_acc))
saveRDS(list(model=model_acc,standardization=standardization_acc),file.path(output_dir_alps,"issf_acc_generic.rds"))

# VISUALISATION n°1: coefficients of the GPS and ACC models ----
coefficient_labels <- c(
  settlement_density_z="Settlement density",
  step_length_km="Step length",
  log_step_length_km="Log step length",
  cos_turning_angle="Cosine turning angle",
  elevation_100m_z="Elevation",
  prop_low_vegetation_5cells_z="Low vegetation",
  ruggedness_100m_z="Ruggedness")

coefficient_table <- function(model,dataset){
  confint(model) %>% as.data.frame() %>% tibble::rownames_to_column("Factor") %>%
    filter(Factor %in% names(coefficient_labels)) %>%
    rename(Lower=2,Upper=3) %>%
    mutate(dataset=dataset)}

graph_comparison <- bind_rows(coefficient_table(model_gps,"GPS"),coefficient_table(model_acc,"ACC")) %>%
  mutate(Factor=factor(Factor,levels=rev(names(coefficient_labels))),
         dataset=factor(dataset,levels=c("GPS","ACC")))

coefs_comparison <- ggplot(graph_comparison,aes(x=Estimate,y=Factor,colour=dataset))+
  geom_vline(xintercept=0,linetype="dashed",color="gray",linewidth=0.5)+
  geom_linerange(aes(xmin=Lower,xmax=Upper),linewidth=0.8,position=position_dodge(width=0.5))+
  geom_point(size=2.4,position=position_dodge(width=0.5))+
  scale_colour_manual(values=c(GPS=clr_gps,ACC=clr_acc),name=NULL)+
  scale_y_discrete(labels=coefficient_labels)+
  labs(x=NULL,y=NULL)+
  theme_minimal()+
  theme(text=element_text(family=font_family,size=14),
        axis.text.x=element_text(family=font_family,size=14),
        axis.text.y=element_text(family=font_family,size=15),
        legend.position="top")

print(coefs_comparison)
ggsave(coefs_comparison,filename=file.path(issf_directory,"coeffs_gps_vs_acc.svg"),width=7,height=3,dpi=400)


#------------------------------------------------------------------------------- STEP 2: alpine raster of the RSS change due to buildings (Nourani) ----
#' **Steps:**
#' (i) alpine grid: the Overture raster (100 m) cropped and masked by the Alpine
#'     area polygon; elevation and ruggedness resampled on it, low vegetation
#'     recomputed from the land cover as in the annotation (Nourani STEP 0A);
#' (ii) grid as a data frame (Nourani STEP 1), standardized with the mean and sd
#'      of the ACC tracking data, movement terms at their sample mean, grouping
#'      variables NA (Nourani STEP 2 / TMB_alps.R);
#' (iii) two predictions with predict(): the observed settlement density, and the
#'       same grid with the settlement density set to 0 (counterfactual);
#' (iv) change in RSS = 1 - exp(eta_observed - eta_no_building): 0 = no building,
#'      1 = selection entirely lost; written as a raster for QGIS;
#' (v) check: the change equals its closed form 1 - exp(beta x (z_observed - z_0)).

# 2.1 Alpine area and grid ----
overture <- rast(overture_file)
alps_area <- vect(alps_perimeter_file)
if(crs(alps_area)=="") stop("The Alpine area polygon has no CRS: set it (e.g. crs(alps_area) <- 'EPSG:4326') before projecting")
alps_area <- aggregate(project(alps_area,crs(overture)))   # one (multi)polygon in the raster CRS
if(is.null(intersect(ext(alps_area),ext(overture)))) stop("The Alpine area polygon and the Overture raster do not overlap: check the CRS of the polygon")
cat("Alpine area:",round(expanse(alps_area,unit="km")),"km2 |",
    "share inside the Overture raster:",round(100*expanse(crop(alps_area,ext(overture)),unit="km")/expanse(alps_area,unit="km"),1),"%\n")

grid_alps <- mask(crop(overture,alps_area,snap="out"),alps_area)
names(grid_alps) <- "settlement_density"

to_alps_grid <- function(file,method){
  x <- rast(file)
  x <- if(same.crs(x,grid_alps)) crop(x,grid_alps,snap="out") else project(x,grid_alps,method=method)
  if(!compareGeom(x,grid_alps,stopOnError=FALSE)){ message(basename(file),": resampled (",method,")"); x <- resample(x,grid_alps,method=method)}
  mask(x,grid_alps)}

static_stack <- rast(lapply(static_files,to_alps_grid,method="bilinear"))
names(static_stack) <- names(static_files)

# Proportion of low vegetation: raster equivalent of extract_landcover_5cells() (classes 5, 6, 7;
# central cell + 4 rook neighbours; NA if any of the five is NA), on the land cover aligned on the grid first
landcover_template <- extend(grid_alps,5)
landcover <- rast(landcover_file)
landcover <- if(same.crs(landcover,grid_alps)) resample(crop(landcover,landcover_template,snap="out"),landcover_template,method="near") else project(landcover,landcover_template,method="near")
low_vegetation <- ifel(landcover %in% low_vegetation_classes,1,0)
rook_window <- matrix(c(NA,1,NA,1,1,1,NA,1,NA),nrow=3,byrow=TRUE)
prop_low_vegetation <- focal(low_vegetation,w=rook_window,fun="mean",na.rm=FALSE) %>% crop(grid_alps,snap="out") %>% mask(grid_alps)
static_stack <- c(static_stack,prop_low_vegetation,grid_alps)
names(static_stack) <- covariates[c(2,3,4,1)]

# 2.2 Grid as a data frame (Nourani STEP 1), standardized with the ACC parameters (Nourani STEP 2) ----
alps_df <- as.data.frame(static_stack,xy=TRUE,cells=TRUE,na.rm=TRUE) %>% drop_na()
cat("Alpine cells with complete covariates:",nrow(alps_df),"(",round(nrow(alps_df)*prod(res(grid_alps))/1e6),"km2 )\n")

z_score_acc <- function(x,variable){s <- standardization_acc[standardization_acc$variable==variable,]; (x-s$center)/s$scale}
for(v in covariates) alps_df[[paste0(v,"_z")]] <- z_score_acc(alps_df[[v]],v)
# Movement terms at their sample mean (TMB_alps.R: step_length = mean(step_length)); grouping variables NA
# for a population-level prediction (TMB_alps.R: stratum_ID = NA)
alps_df <- alps_df %>%
  mutate(step_length_km=mean(data_acc$step_length_km),log_step_length_km=mean(data_acc$log_step_length_km),
         cos_turning_angle=mean(data_acc$cos_turning_angle),stratum_ID=NA,animal_ID=NA)

# Counterfactual grid: the same cells with the settlement density set to 0 (Nourani: only the week changes)
alps_df_no_building <- alps_df %>% mutate(settlement_density=0,settlement_density_z=z_score_acc(0,"settlement_density"))

# 2.3 Predictions with predict(), population level, by chunks of rows ----
predict_grid <- function(model,newdata){
  chunks <- split(seq_len(nrow(newdata)),ceiling(seq_len(nrow(newdata))/prediction_chunk_rows))
  eta <- numeric(nrow(newdata))
  for(i in seq_along(chunks)){
    cat("  predict(): chunk",i,"of",length(chunks),"\n"); flush.console()
    eta[chunks[[i]]] <- predict(model,newdata=newdata[chunks[[i]],],type="link",re.form=NA)}
  eta}

cat("Observed settlement density\n");  eta_observed <- predict_grid(model_acc,alps_df)
cat("Settlement density set to 0\n");  eta_no_building <- predict_grid(model_acc,alps_df_no_building)

# 2.4 Change in RSS between the two scenarios, written as a raster ----
rss_ratio <- exp(eta_observed-eta_no_building)   # RSS observed / RSS no building (1 = no building)
rss_change <- 1-rss_ratio                        # share of the RSS lost to buildings (0 = no building, 1 = entirely lost)

to_raster <- function(values,layer_name){
  r <- rast(grid_alps,nlyrs=1,vals=NA_real_,names=layer_name)
  r[alps_df$cell] <- values
  r}
rss_change_alps <- c(to_raster(rss_change,"rss_change_due_to_buildings"),to_raster(rss_ratio,"rss_observed_vs_no_building"))
writeRaster(rss_change_alps,file.path(output_dir_alps,"issf_alps_acc_rss_change_buildings.tif"),overwrite=TRUE,wopt=write_options)



#------------------------------------------------------------------------------- STEP 3: settlement coefficients by behavior (visualisation n°2) ----
#' **Steps:**
#' (i) fit the interaction model: beta = selection for landings followed by resting,
#'     delta = deviation for landings followed by feeding, beta + delta = feeding;
#' (ii) one common raw Q05-Q95 contrast of settlement density, supported by both
#'      datasets (largest Q05, smallest Q95 of the available destinations), converted
#'      to each model's standardized scale;
#' (iii) table: settlement coefficient [95% CI] and Q95-Q05 log-RSS [95% CI] for the
#'       GPS model, the generic ACC model, and resting and feeding from the
#'       interaction model.

# 3.1 Interaction model ----
model_interaction <- fit_issf(formula_interaction,data_acc)
print(summary(model_interaction))
saveRDS(list(model=model_interaction,standardization=standardization_acc),file.path(output_dir_alps,"issf_acc_interaction.rds"))

fixed_interaction <- fixef(model_interaction)$cond
vcov_interaction <- vcov(model_interaction)$cond
settlement_interaction_term <- grep("feeding.*settlement_density_z|settlement_density_z.*feeding",names(fixed_interaction),value=TRUE)
stopifnot(length(settlement_interaction_term)==1L)

beta_resting <- fixed_interaction[["settlement_density_z"]]
se_resting <- sqrt(vcov_interaction["settlement_density_z","settlement_density_z"])
delta_feeding <- fixed_interaction[[settlement_interaction_term]]
se_delta <- sqrt(vcov_interaction[settlement_interaction_term,settlement_interaction_term])
beta_feeding <- beta_resting+delta_feeding
se_feeding <- sqrt(se_resting^2+se_delta^2+2*vcov_interaction["settlement_density_z",settlement_interaction_term])
cat(sprintf("Settlement density: resting %.3f (SE %.3f) | feeding %.3f (SE %.3f) | deviation %.3f (SE %.3f, p = %.3f)\n",
            beta_resting,se_resting,beta_feeding,se_feeding,delta_feeding,se_delta,2*pnorm(-abs(delta_feeding/se_delta))))

# 3.2 Common raw Q05-Q95 contrast of settlement density ----
# Largest Q05 and smallest Q95 over the available destinations of both datasets, so that the
# contrast is supported by both; converted to each model's standardized scale
available_quantiles <- bind_rows(
  GPS=tibble(q05=quantile(data_gps$settlement_density[data_gps$used==0L],0.05),q95=quantile(data_gps$settlement_density[data_gps$used==0L],0.95)),
  ACC=tibble(q05=quantile(data_acc$settlement_density[data_acc$used==0L],0.05),q95=quantile(data_acc$settlement_density[data_acc$used==0L],0.95)),.id="dataset")
common_contrast <- tibble(q05_raw=max(available_quantiles$q05),q95_raw=min(available_quantiles$q95))
contrast_z <- function(standardization){
  s <- standardization[standardization$variable=="settlement_density",]
  (common_contrast$q95_raw-common_contrast$q05_raw)/s$scale}
contrast_z_gps <- contrast_z(standardization_gps)
contrast_z_acc <- contrast_z(standardization_acc)
print(available_quantiles); print(common_contrast)

# 3.3 Rows of the table: coefficient, CI, Q95-Q05 log-RSS, CI ----
settlement_row <- function(estimate,se,contrast,model_group,response){
  tibble(model_group=model_group,response=response,
         coefficient=estimate,coefficient_low=estimate-1.96*se,coefficient_high=estimate+1.96*se,
         log_rss=estimate*contrast,log_rss_low=(estimate-1.96*se)*contrast,log_rss_high=(estimate+1.96*se)*contrast)}
se_of <- function(model,term) sqrt(diag(vcov(model)$cond))[[term]]

settlement_table <- bind_rows(
  settlement_row(fixef(model_gps)$cond[["settlement_density_z"]],se_of(model_gps,"settlement_density_z"),contrast_z_gps,
                 "GPS model","Landing site selection"),
  settlement_row(fixef(model_acc)$cond[["settlement_density_z"]],se_of(model_acc,"settlement_density_z"),contrast_z_acc,
                 "ACC models","Landing site selection"),
  settlement_row(beta_resting,se_resting,contrast_z_acc,"ACC models","Landing site selection | resting"),
  settlement_row(beta_feeding,se_feeding,contrast_z_acc,"ACC models","Landing site selection | feeding"))
print(settlement_table,width=Inf)
write_csv(settlement_table,file.path(issf_directory,"settlement_coefficients_gps_acc.csv"))

# VISUALISATION n°2: table of the settlement density coefficients ----
# VISUALISATION n°2: table of the settlement density coefficients ----
# same layout and style as the GLMM table (Covariates_HFI_Selection.R, section 7), with the
# coefficient and its confidence interval in one cell
settlement_summary_table <- settlement_table %>%
  dplyr::mutate(
    hfi_coefficient_CI = sprintf("%.3f [%.3f; %.3f]",coefficient,coefficient_low,coefficient_high),
    Q95_Q05_log_rss_CI = sprintf("%.3f [%.3f; %.3f]",log_rss,log_rss_low,log_rss_high)) %>%
  dplyr::select(dataset = model_group,response,hfi_coefficient_CI,Q95_Q05_log_rss_CI)

settlement_gt <- settlement_summary_table %>%
  gt::gt(groupname_col = "dataset") %>%
  gt::tab_header(
    title = gt::md("**Settlement density influence on landing site selection in GPS and ACC models**"),
    subtitle = paste0(
      "Settlement density; common raw contrast: Q05 = ",
      sprintf("%.3f",common_contrast$q05_raw),
      "; Q95 = ",
      sprintf("%.3f",common_contrast$q95_raw))) %>%
  gt::cols_label(
    response = "Model response",
    hfi_coefficient_CI = "Settlement coefficient [95% CI]",
    Q95_Q05_log_rss_CI = "Q95−Q05 log-RSS [95% CI]") %>%
  gt::cols_align(
    align = "left",
    columns = response) %>%
  gt::cols_align(
    align = "center",
    columns = c(hfi_coefficient_CI,Q95_Q05_log_rss_CI)) %>%
  gt::tab_style(
    style = list(
      gt::cell_fill(color = "#ECECF8"),
      gt::cell_text(weight = "bold")),
    locations = gt::cells_column_labels()) %>%
  gt::tab_style(
    style = list(
      gt::cell_fill(color = "#F4F4FA"),
      gt::cell_text(weight = "bold")),
    locations = gt::cells_row_groups()) %>%
  gt::tab_options(
    table.width = gt::pct(100),
    table.font.size = gt::px(14),
    heading.title.font.size = gt::px(19),
    heading.subtitle.font.size = gt::px(13),
    column_labels.font.weight = "bold",
    data_row.padding = gt::px(7)) %>%
  gt::tab_source_note(
    source_note = gt::md(
      paste0(
        "A negative coefficient indicates avoidance of settlements when selecting a landing site. ",
        "The feeding row is the selection of landing sites followed by feeding (interaction model). ",
        "For the feeding versus resting row, the coefficient is the deviation of feeding landings ",
        "from resting landings: a negative value indicates stronger avoidance when landing to feed.")))

print(settlement_gt)
gt::gtsave(settlement_gt,file.path(issf_directory,"settlement_coefficients_gps_acc.html"),inline_css = TRUE)
