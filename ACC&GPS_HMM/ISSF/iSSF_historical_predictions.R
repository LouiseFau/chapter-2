#'-------------------------------------------------------------------------------
#' Title: iSSF historical predictions in Switzerland (RegBL, 1918-2026) ----
#' Author: Louise Faure
#' Date: 08.10.2026
#' **Purpose:** apply the best landing iSSF selected in iSSF_covariate_selection.R
#' to a Swiss 100 m prediction grid for ten historical settlement scenarios
#' (RegBL 1918-2026), holding topography, land cover and movement terms constant,
#' and map the population-level log-RSS for each date. Based on E. Nourani
#' 05_INLA_prediction_map.R (grid preparation) and
#' cluster_prep/glmmTMB_preds/TMB_alps.R (one prediction per period).
#' **Steps:**
#' (1) prepare the Swiss prediction grid: crop the Alpine covariate layers to
#'     Switzerland, convert the grid to a data frame, standardize with the mean
#'     and sd of the tracking data, check the covariate ranges;
#' (2) predict the population-level log-RSS for each date and write the rasters;
#' (3) check that the coefficients are those of the selected model and that the
#'     manual linear predictor equals predict().
#' ------------------------------------------------------------------------------

library(tidyverse)
library(terra)
library(glmmTMB)

# Paths and parameters ----
model_file <- "/Users/louisefaure/Library/CloudStorage/OneDrive-Personnel/THESE/CHAPITRE 2/git/chapter-2/ACC&GPS_HMM/ISSF/issf_best_model_60.rds" # exported by iSSF_covariate_selection.R
raster_dir <- "/Users/louisefaure/Desktop/dossier sans titre/Rasters"
settlement_dir <- file.path(raster_dir,"historical_settlements_rasters")
output_dir_ch <- "/Users/louisefaure/Desktop/dossier sans titre/iSSF_switzerland"
dir.create(output_dir_ch,recursive=TRUE,showWarnings=FALSE)

# Alpine covariate layers used to annotate the tracking data (Data_processing_annotation.R, STEP 6).
# Elevation and ruggedness are cropped as they are. The proportion of low vegetation is recomputed from
# landcover_100m.tif with the same rule as the annotation (classes 5, 6, 7; central cell + 4 rook neighbours).
static_files <- c(elevation_100m=file.path(raster_dir,"elevation_100m.tif"),
                  ruggedness_100m=file.path(raster_dir,"ruggedness_100m.tif"))
landcover_file <- file.path(raster_dir,"landcover_100m.tif")
low_vegetation_classes <- c(5,6,7)
static_covariates <- c(names(static_files),"prop_low_vegetation_5cells")

dates <- c(1918,1945,1960,1970,1980,1990,2000,2010,2015,2026)
settlement_files <- file.path(settlement_dir,sprintf("built_25m_share_1km2_100m_%d.tif",dates))

n_check_cells <- 5000L  # cells used to compare predict() with the manual linear predictor
write_options <- list(datatype="FLT4S",gdal=c("COMPRESS=ZSTD","PREDICTOR=3","TILED=YES"))
set.seed(500)
stopifnot(file.exists(c(model_file,static_files,landcover_file,settlement_files)))

# Selected model, its data and its standardization parameters ----
model_export <- readRDS(model_file)
model_best <- model_export$model
standardization <- model_export$standardization
model_frame <- model_best$frame                      # data the model was fitted on (z-scores, movement terms, IDs)
fixed_effects <- fixef(model_best)$cond
cat("Model:",model_export$model_name,"| strata:",n_distinct(model_frame$stratum_ID),"| individuals:",n_distinct(model_frame$animal_ID),"\n")
print(round(fixed_effects,4))

z_score <- function(x,variable){s <- standardization[standardization$variable==variable,]; (x-s$center)/s$scale}


#------------------------------------------------------------------------------- STEP 1: Swiss prediction grid (Nourani STEP 0A to 2) ----
# 1.1 The RegBL 2026 layer defines the grid: extent, 100 m resolution and mask (NA outside Switzerland) ----
settlement_stack <- rast(settlement_files)
names(settlement_stack) <- paste0("settlement_",dates)
grid_ch <- settlement_stack[["settlement_2026"]]

# 1.2 Crop the Alpine layers to the Swiss grid ----
# The Alpine layers contain Switzerland, so cropping is enough. A layer is resampled only
# if it is not on the RegBL grid (other origin or CRS), and this is reported.
crop_to_grid <- function(file){
  x <- crop(rast(file),grid_ch,snap="out")
  if(!compareGeom(x,grid_ch,stopOnError=FALSE)){
    message(basename(file),": not on the RegBL grid, resampled (bilinear)")
    x <- resample(x,grid_ch,method="bilinear")}
  mask(x,grid_ch)}

static_stack <- rast(lapply(static_files,crop_to_grid))
names(static_stack) <- names(static_files)

# 1.2b Proportion of low vegetation over the central cell and its four rook neighbours ----
# Raster equivalent of extract_landcover_5cells() in Data_processing_annotation.R:
# - low vegetation = landcover classes 5, 6, 7 (1), anything else (0);
# - mean over 5 cells: the window has NA in the corners so that terra ignores the diagonal cells (rook);
# - na.rm = FALSE: NA as soon as one of the five cells is NA (complete_neighbourhood);
# - computed on an extent buffered by 500 m so that Swiss border cells keep their neighbours.
landcover <- crop(rast(landcover_file),ext(grid_ch)+500,snap="out")
low_vegetation <- ifel(landcover %in% low_vegetation_classes,1,0)
rook_window <- matrix(c(NA,1,NA,1,1,1,NA,1,NA),nrow=3,byrow=TRUE)
prop_low_vegetation <- focal(low_vegetation,w=rook_window,fun="mean",na.rm=FALSE)
prop_low_vegetation <- crop(prop_low_vegetation,grid_ch,snap="out")
if(!compareGeom(prop_low_vegetation,grid_ch,stopOnError=FALSE)){
  message("landcover: not on the RegBL grid, resampled (near)")
  prop_low_vegetation <- resample(prop_low_vegetation,grid_ch,method="near")}
static_stack <- c(static_stack,mask(prop_low_vegetation,grid_ch))
names(static_stack) <- static_covariates

# 1.3 Grid as a data frame (Nourani STEP 1): one row per Swiss cell, complete for all layers ----
# Cells missing in any layer (static or any date) are dropped so that all dates share the same domain.
swiss_df <- as.data.frame(c(static_stack,settlement_stack),xy=TRUE,cells=TRUE,na.rm=TRUE) %>% drop_na()
cat("Swiss cells:",global(not.na(grid_ch),"sum")[1,1],"| cells with complete covariates:",nrow(swiss_df),"\n")

# 1.4 Standardize with the mean and sd of the tracking data (Nourani STEP 2) ----
# Movement terms are fixed at their sample mean (Nourani TMB_alps.R: step_length = mean(step_length)):
# a constant shifts all cells equally and cancels in the RSS ratio between scenarios.
# NA grouping variables = population-level prediction (glmmTMB recommendation).
for(v in static_covariates) swiss_df[[paste0(v,"_z")]] <- z_score(swiss_df[[v]],v)
swiss_df <- swiss_df %>%
  mutate(step_length_km=mean(model_frame$step_length_km),
         log_step_length_km=mean(model_frame$log_step_length_km),
         cos_turning_angle=mean(model_frame$cos_turning_angle),
         stratum_ID=NA,animal_ID=NA)

# 1.5 Covariate ranges: tracking data versus prediction grid, on the z scale (Nourani STEP 4 opening check) ----
# Nourani compares max(new$dem_200) with max(og$dem_200) and greys out the cells outside the training
# range. Here: quantiles of each covariate at the eagle locations (training) and over the Swiss grid,
# and the share of grid cells outside the training range, i.e. where the model extrapolates.
quantile_text <- function(x) sprintf("%.2f / %.2f / %.2f",quantile(x,0.05),median(x),quantile(x,0.95))
range_check <- function(variable,training,grid) tibble(variable=variable,
                                                       training_q05_q50_q95=quantile_text(training),grid_q05_q50_q95=quantile_text(grid),
                                                       pct_cells_outside_training=round(100*mean(grid<min(training) | grid>max(training)),2))

covariate_ranges <- bind_rows(
  map_dfr(static_covariates,~range_check(paste0(.x,"_z"),model_frame[[paste0(.x,"_z")]],swiss_df[[paste0(.x,"_z")]])),
  map_dfr(dates,~range_check(paste0("settlement_",.x,"_z"),model_frame$settlement_density_z,
                             z_score(swiss_df[[paste0("settlement_",.x)]],"settlement_density"))))
print(covariate_ranges,n=Inf)
write_csv(covariate_ranges,file.path(output_dir_ch,"covariate_ranges_training_vs_grid.csv"))
saveRDS(swiss_df,file.path(output_dir_ch,"swiss_prediction_grid_100m.rds"))


#------------------------------------------------------------------------------- STEP 2: predictions for each date (Nourani TMB_alps.R STEP 2) ----
# 2.1 Population-level linear predictor = sum of coefficient x covariate ----
# Identical to predict(model, newdata, type = "link", re.form = NA) (verified in STEP 3), computed as a
# matrix product so that 4 million cells x 11 scenarios run in seconds on a laptop, without the cluster.
predictor_columns <- names(fixed_effects)
linear_predictor <- function(newdata) as.vector(as.matrix(newdata[,predictor_columns]) %*% fixed_effects)

# 2.2 One scenario per date: a new dataset where only the settlement density changes (Nourani: only the week changes) ----
eta_dates <- sapply(dates,function(d) linear_predictor(mutate(swiss_df,settlement_density_z=z_score(.data[[paste0("settlement_",d)]],"settlement_density"))))
colnames(eta_dates) <- paste0("y",dates)

# 2.3 Back to rasters through the cell numbers ----
to_raster <- function(values_matrix){
  r <- rast(grid_ch,nlyrs=ncol(values_matrix),vals=NA_real_,names=colnames(values_matrix))
  r[swiss_df$cell] <- values_matrix
  r}

log_rss_by_date <- to_raster(eta_dates)              # population-level log-RSS, one layer per date
writeRaster(log_rss_by_date,file.path(output_dir_ch,"issf_log_rss_by_date.tif"),overwrite=TRUE,wopt=write_options)
writeRaster(exp(log_rss_by_date),file.path(output_dir_ch,"issf_rss_by_date.tif"),overwrite=TRUE,wopt=write_options)


#------------------------------------------------------------------------------- STEP 3: animation through time ----
# Layout of each frame:
#   - centre: Switzerland. Background = hillshade (black to white) under a semi-transparent elevation
#     layer (light to dark brown); over it, the RSS reduction due to buildings,
#       rss_change = 1 - RSS_date / RSS_reference = 1 - exp(eta_date - eta_reference),
#     in grey with an opacity equal to the reduction (0 transparent, 1 opaque); then lakes, main rivers,
#     the Swiss border, the nests and the main cities;
#   - right: a single bar, the area (km2) of cells with rss_change > `bar_threshold`, growing over time;
#   - bottom: a timeline with one tick per date and a cursor on the current date.
# Reference scenario (`reference_scenario`):
#   "no_building": the same landscape with the settlement density set to 0 everywhere. The 1918 frame then
#                  already shows the cells lost to the buildings that existed in 1918, and the later
#                  frames add the losses of 1918-1945, 1945-1960, ... on top of them (cumulative);
#   "1918":        the 1918 landscape itself. The 1918 frame is empty and the frames show the loss since 1918.
# Frames are saved as PNG, one per date, then assembled into a GIF and an MP4.
library(patchwork)
library(gifski)
library(av)
library(sf)
library(ggnewscale)
library(shadowtext)   # text with a white halo
select <- dplyr::select; filter <- dplyr::filter   # terra masks these two dplyr verbs

# Parameters ----
reference_scenario <- "no_building"   # "no_building" or "1918" (see above)
bar_threshold <- 0.7              # RSS change above which a cell counts in the bar
legend_position <- c(0.99,0.02)   # colour bar anchor (its bottom-right corner) inside the map panel, 0-1 from bottom-left
reference_text <- if(reference_scenario=="no_building") "relative to a landscape without buildings" else sprintf("since the reference year %d",dates[1])
legend_title <- sprintf("Reduction in relative selection strength (%%)\n%s",reference_text)
bar_title <- sprintf("Cumulative area (km²)\nwhere relative selection\nstrength has fallen by\nmore than %d %%,\n%s",round(100*bar_threshold),reference_text)
include_reference_frame <- TRUE   # TRUE: first frame = 1918; FALSE: start at 1945
display_factor <- 2               # aggregation factor for drawing only (2 = 200 m); 1 draws at 100 m, slower
frame_delay <- 1.2                # seconds per frame
font_family <- "Times New Roman"  # all text; "Times" also works on macOS
clr_loss <- "grey45"              # medium grey for the RSS reduction on the map (transparent -> this grey)
clr_bar <- "grey35"               # medium-dark grey for the bar and the timeline cursor
clr_text <- "grey25"              # dark grey for all text
clr_marks <- "grey55"             # lighter grey for the scale bar and the north arrow
clr_timeline <- "grey50"          # grey for the timeline and the outline of the bar
clr_water <- "#9DC3E6"; clr_water_dark <- "#2F5597"; water_alpha <- 0.55   # rivers, lakes and river names
clr_nest <- "#7B1E3A"; nest_alpha <- 0.75   # burgundy nests, thin black outline
loss_alpha_max <- 0.9             # opacity of a 100 % RSS reduction (1 = fully opaque)
elevation_colours <- c("#F3E4CF","#DDB98A","#C08A4E","#8F5A2A","#5A3416")   # light to dark brown
elevation_alpha <- 0.7            # opacity of the elevation layer over the hillshade

hillshade_file <- "/Users/louisefaure/Desktop/PROJET QGIS/Carte Poster ECBB/Couches/ombrage.tif"
tlmregio_dir <- "/Users/louisefaure/Downloads/swisstlmregio_2025_2056.gpkg"
tlmregio_file <- file.path(tlmregio_dir,"swissTLMRegio_Product_LV95.gpkg")
tlmregio_boundaries_file <- file.path(tlmregio_dir,"swissTLMRegio_BOUNDARIES_LV95.gpkg")
nest_file <- "/Users/louisefaure/Library/CloudStorage/OneDrive-Personnel/THESE/CHAPITRE 2/git/chapter-2/DONNEES AIGLES/nest site location/nest_site_location/nest_site_location.shp"
# Rivers drawn, by their swissTLMRegio names (column namn; the Rhône and the Rhine carry French names on some reaches)
main_rivers <- c("Aare","Inn","Limmat","Reuss","Rhein","Le Rhin","Thur","Ticino","Le Rhône","Rotten")
# River names on the map: position (LV95) and text angle, placed by hand along a straight reach of each river
river_labels <- tribble(~river,~x,~y,~angle,
                        "Rhône",   2610000,1126000, 35,
                        "Rhein",   2760000,1215000, 80,
                        "Aare",    2636000,1240000, 25,
                        "Reuss",   2668000,1233000, 85,
                        "Limmat",  2673000,1254000,-30,
                        "Thur",    2714000,1272000,-10,
                        "Ticino",  2719000,1127000,-75,
                        "Inn",     2792000,1166000, 35)
scale_km <- 25                    # length of the scale bar
min_lake_km2 <- 8                 # lakes smaller than this are not drawn
main_cities <- tribble(~name,~x,~y,          # LV95 (EPSG:2056)
                       "Zürich",2683000,1248000,  "Genève",2500000,1118000,  "Basel",2611500,1267500,
                       "Bern",2600000,1200000,    "Lausanne",2538000,1152500, "Luzern",2666000,1211500,
                       "Lugano",2717500,1096000,  "St. Gallen",2746000,1254500,"Chur",2759500,1190500,
                       "Sion",2594000,1120000)
stopifnot(file.exists(c(hillshade_file,tlmregio_file,tlmregio_boundaries_file,nest_file)))
frames_dir <- file.path(output_dir_ch,"animation_frames")
dir.create(frames_dir,showWarnings=FALSE)
map_crs <- crs(grid_ch)

# 3.1 RSS reduction relative to the reference scenario and area above the bar threshold, full 100 m resolution ----
# Reference log-RSS: the 1918 layer, or the 1918 layer with its building term replaced by that of a
# density of zero (only the settlement term differs between scenarios, so it is the only term to change)
log_rss_reference <- if(reference_scenario=="no_building"){
  b_settlement <- fixed_effects[["settlement_density_z"]]
  z_1918 <- mask(z_score(settlement_stack[["settlement_1918"]],"settlement_density"),log_rss_by_date[["y1918"]])
  log_rss_by_date[["y1918"]]-b_settlement*(z_1918-z_score(0,"settlement_density"))
}else log_rss_by_date[["y1918"]]
rss_change_vs_1918 <- 1-exp(log_rss_by_date-log_rss_reference)   # name kept for the rest of the script
names(rss_change_vs_1918) <- paste0("y",dates)
writeRaster(rss_change_vs_1918,file.path(output_dir_ch,sprintf("issf_rss_change_vs_%s_by_date.tif",reference_scenario)),
            overwrite=TRUE,wopt=write_options)

cell_area_km2 <- prod(res(grid_ch))/1e6
area_by_date <- tibble(date=dates,
                       area_lost_km2=as.numeric(global(rss_change_vs_1918>bar_threshold,"sum",na.rm=TRUE)[,1])*cell_area_km2,
                       area_grid_km2=nrow(swiss_df)*cell_area_km2) %>%
  mutate(pct_grid_lost=100*area_lost_km2/area_grid_km2)
print(area_by_date)
write_csv(area_by_date,file.path(output_dir_ch,sprintf("area_lost_vs_%s_by_date.csv",reference_scenario)))

# 3.2 Vector context from swissTLMRegio (LV95), reprojected to the raster CRS ----
# Layer and column names are found by pattern because they differ between releases: the matched
# layer is printed, and st_layers(tlmregio_file) lists all layers if a pattern has to be adjusted.
find_layer <- function(file,pattern){
  layers <- st_layers(file)$name
  hits <- grep(pattern,layers,ignore.case=TRUE,value=TRUE)
  if(length(hits)==0) stop("No layer matches '",pattern,"' in ",basename(file),". Layers: ",paste(layers,collapse=", "))
  message("'",pattern,"' -> ",basename(file),": ",hits[1]); hits[1]}
read_tlm <- function(pattern,file=tlmregio_file) st_read(file,layer=find_layer(file,pattern),quiet=TRUE) %>% st_transform(map_crs)
text_matches <- function(x,pattern){   # TRUE for rows where any character column matches the pattern
  chr <- Filter(is.character,st_drop_geometry(x))
  if(length(chr)==0) return(rep(FALSE,nrow(x)))
  Reduce(`|`,lapply(chr,function(v) grepl(pattern,v,ignore.case=TRUE)))}

# Swiss border: swissTLMRegio also covers the neighbouring countries, so Switzerland is selected by code or name
countries <- read_tlm("landesgebiet",tlmregio_boundaries_file)
is_switzerland <- text_matches(countries,"^CH$|^CHE$|Schweiz|Switzerland|Suisse|Svizzera")
if(!any(is_switzerland)) is_switzerland <- rep(TRUE,nrow(countries))
swiss_border <- countries[is_switzerland,] %>% st_union() %>% st_sf()
sf_use_s2(FALSE)   # planar operations on projected data
# Clip to Switzerland and keep only the geometry type asked for: st_intersection can return mixed
# collections (a river touching the border yields points, drawn as dots by geom_sf), which are dropped here.
clip_to_switzerland <- function(x,type){
  suppressWarnings(st_intersection(x,swiss_border)) %>% st_collection_extract(type) %>% filter(!st_is_empty(.))}

# Lakes clipped to Switzerland. Rivers are clipped to Switzerland buffered by 2 km: a river that forms the
# border (the Rhine from Basel to Lake Constance, for instance) would otherwise be cut along the border
# line and the grey border would show through it.
lakes <- read_tlm("lake") %>% filter(as.numeric(st_area(.))/1e6>=min_lake_km2) %>% clip_to_switzerland("POLYGON")
swiss_border_buffered <- st_buffer(swiss_border,2000)
rivers <- read_tlm("flowing") %>% filter(namn %in% main_rivers) %>%
  st_intersection(swiss_border_buffered) %>% suppressWarnings() %>% st_collection_extract("LINESTRING") %>%
  filter(!st_is_empty(.))
river_labels <- st_as_sf(river_labels,coords=c("x","y"),crs=2056) %>% st_transform(map_crs) %>%
  mutate(x=st_coordinates(.)[,1],y=st_coordinates(.)[,2]) %>% st_drop_geometry()
cat("Lakes:",nrow(lakes),"| river names found:",paste(sort(unique(rivers$namn)),collapse=", "),"\n")

cities <- st_as_sf(main_cities,coords=c("x","y"),crs=2056) %>% st_transform(map_crs) %>%
  mutate(x=st_coordinates(.)[,1],y=st_coordinates(.)[,2])

# Nest sites of the Swiss-born individuals: those located inside Switzerland
nests <- st_read(nest_file,quiet=TRUE) %>% st_transform(map_crs) %>% st_filter(swiss_border) %>%
  mutate(x=st_coordinates(.)[,1],y=st_coordinates(.)[,2])
cat("Nests in Switzerland:",nrow(nests),"\n")

# North arrow: an arrowhead only (filled triangle) with an N above it, in the same grey as the scale bar
north_grob <- grid::grobTree(
  grid::polygonGrob(x=c(0.5,0.22,0.78),y=c(0.68,0.08,0.08),gp=grid::gpar(fill=clr_marks,col=NA)),
  grid::textGrob("N",x=0.5,y=0.86,gp=grid::gpar(fontfamily=font_family,fontsize=8,col=clr_marks)))

# Map limits: Switzerland plus a margin on every side for the border and the neighbouring country names
swiss_bbox <- st_bbox(swiss_border)
margin <- 0.05
w <- unname(swiss_bbox$xmax-swiss_bbox$xmin); h <- unname(swiss_bbox$ymax-swiss_bbox$ymin)   # unname: st_bbox values carry names
map_xlim <- unname(c(swiss_bbox$xmin-margin*w,swiss_bbox$xmax+margin*w))
map_ylim <- unname(c(swiss_bbox$ymin-margin*h,swiss_bbox$ymax+margin*h))
# Scale bar: bottom-left of the map, label centred above the bar, north arrow centred above the label (map units)
scale_bar <- tibble(x0=map_xlim[1]+0.03*w,x1=map_xlim[1]+0.03*w+scale_km*1000,y=map_ylim[1]+0.035*h)
scale_centre <- unname((scale_bar$x0+scale_bar$x1)/2)
north_width <- 0.035*w   # width of the north arrow (map units), centred above the scale label
north_box <- c(xmin=scale_centre-north_width/2,xmax=scale_centre+north_width/2,
               ymin=unname(scale_bar$y)+0.035*h,ymax=unname(scale_bar$y)+0.035*h+north_width*1.3)
map_aspect <- diff(map_ylim)/diff(map_xlim)   # used to size the figure so that the map fills its panel

# 3.3 Raster background: hillshade and elevation, cropped to Switzerland, aggregated for drawing ----
to_display <- function(x) if(display_factor>1) aggregate(x,fact=display_factor,fun="mean",na.rm=TRUE) else x
swiss_vect <- vect(swiss_border)
to_swiss_grid <- function(file){   # crop to the Swiss extent, align on the RegBL grid, mask with the border
  x <- rast(file)
  x <- crop(x,project(vect(ext(grid_ch),crs=map_crs),crs(x)),snap="out")
  x <- if(same.crs(x,grid_ch)) resample(x,grid_ch,method="bilinear") else project(x,grid_ch,method="bilinear")
  mask(x,swiss_vect)}
hillshade_df <- as.data.frame(to_display(to_swiss_grid(hillshade_file)),xy=TRUE,na.rm=TRUE) %>% rename(shade=3)
elevation_df <- as.data.frame(to_display(to_swiss_grid(static_files[["elevation_100m"]])),xy=TRUE,na.rm=TRUE) %>% rename(elevation=3)

change_display <- to_display(mask(rss_change_vs_1918,swiss_vect))
change_df <- map(seq_along(dates),~as.data.frame(change_display[[.x]],xy=TRUE,na.rm=TRUE) %>%
                   rename(change=3) %>%
                   filter(change>0))   # gains (negative change) and unchanged cells are left transparent

# 3.4 Panels: map, bar, timeline ----
loss_gradient <- alpha(clr_loss,seq(0,loss_alpha_max,length.out=100))
legend_theme <- if(packageVersion("ggplot2")>="3.5.0"){
  theme(legend.position="inside",legend.position.inside=legend_position,legend.justification=c(1,0))
}else theme(legend.position=legend_position,legend.justification=c(1,0))

map_plot <- function(i){
  ggplot()+
    # background: hillshade, then semi-transparent elevation
    geom_raster(data=hillshade_df,aes(x,y,fill=shade))+
    scale_fill_gradient(low="black",high="white",guide="none")+
    new_scale_fill()+
    geom_raster(data=elevation_df,aes(x,y,fill=elevation),alpha=elevation_alpha)+
    scale_fill_gradientn(colours=elevation_colours,guide="none")+
    new_scale_fill()+
    # RSS change
    geom_raster(data=change_df[[i]],aes(x,y,fill=change))+
    scale_fill_gradientn(colours=loss_gradient,limits=c(0,1),breaks=c(0,0.5,1),labels=c("0 %","50 %","100 %"),
                         oob=scales::squish,name=legend_title,
                         guide=guide_colourbar(direction="horizontal",title.position="top",
                                               barwidth=unit(2.6,"cm"),barheight=unit(0.25,"cm"),
                                               frame.colour="grey40",ticks.colour="grey40"))+
    # border, then rivers under lakes (both clipped to Switzerland), nests, cities
    geom_sf(data=swiss_border,fill=NA,colour="grey20",linewidth=0.45)+
    geom_sf(data=rivers,colour=alpha(clr_water_dark,water_alpha),linewidth=0.4)+
    geom_sf(data=lakes,fill=alpha(clr_water,water_alpha),colour=alpha(clr_water_dark,water_alpha),linewidth=0.25)+
    geom_text(data=river_labels,aes(x,y,label=river,angle=angle),family=font_family,fontface="italic",size=2.8,
              colour=alpha(clr_water_dark,water_alpha))+
    geom_point(data=nests,aes(x,y),shape=21,size=2,fill=alpha(clr_nest,nest_alpha),colour=alpha("black",nest_alpha),stroke=0.5)+
    geom_point(data=cities,aes(x,y),shape=21,size=1.9,fill=clr_text,colour="white",stroke=0.5)+
    geom_shadowtext(data=cities,aes(x,y,label=name),family=font_family,size=3,colour=clr_text,
                    bg.colour=alpha("white",0.7),bg.r=0.15,nudge_y=0.012*h,vjust=0)+
    # scale bar, label centred above
    annotate("segment",x=scale_bar$x0,xend=scale_bar$x1,y=scale_bar$y,yend=scale_bar$y,colour=clr_marks,linewidth=0.7)+
    annotate("segment",x=c(scale_bar$x0,scale_bar$x1),xend=c(scale_bar$x0,scale_bar$x1),
             y=scale_bar$y-0.006*h,yend=scale_bar$y+0.006*h,colour=clr_marks,linewidth=0.7)+
    annotate("text",x=scale_centre,y=scale_bar$y+0.012*h,label=paste(scale_km,"km"),
             family=font_family,size=2.8,colour=clr_marks,vjust=0)+
    annotation_custom(north_grob,xmin=north_box[["xmin"]],xmax=north_box[["xmax"]],ymin=north_box[["ymin"]],ymax=north_box[["ymax"]])+
    coord_sf(xlim=map_xlim,ylim=map_ylim,expand=FALSE,crs=map_crs,datum=NA)+
    theme_void(base_family=font_family)+
    theme(legend.title=element_text(size=8,colour=clr_text),legend.text=element_text(size=8,colour=clr_text),
          legend.title.align=0)+
    legend_theme}

bar_max <- max(area_by_date$area_lost_km2)*1.12
bar_plot <- function(i){
  area <- area_by_date$area_lost_km2[i]
  ggplot(tibble(x="",y=area),aes(x,y))+
    geom_col(fill=clr_bar,colour=clr_timeline,linewidth=0.4,width=0.45)+
    annotate("text",x=1,y=area,label=if(area>0) paste(format(round(area),big.mark=" "),"km²") else "",
             vjust=-0.6,size=3.8,family=font_family,colour=clr_text)+
    scale_y_continuous(limits=c(0,bar_max),expand=expansion(mult=0))+
    theme_void(base_family=font_family)}

# Bar title: drawn without clipping so that a line wider than the narrow panel is not cut
bar_title_plot <- ggplot()+
  annotate("text",x=0,y=1,label=bar_title,family=font_family,size=2.5,hjust=0.5,vjust=1,lineheight=0.95,colour=clr_text)+
  scale_x_continuous(limits=c(-1,1))+scale_y_continuous(limits=c(0,1))+
  coord_cartesian(clip="off")+
  theme_void(base_family=font_family)

timeline_plot <- function(year,cursor_colour=clr_bar){
  ggplot()+
    annotate("segment",x=min(dates),xend=max(dates),y=0,yend=0,colour=clr_timeline,linewidth=0.8)+
    annotate("segment",x=dates,xend=dates,y=-0.15,yend=0.15,colour=clr_timeline,linewidth=0.5)+
    annotate("text",x=dates,y=-0.55,label=dates,size=3,colour=clr_timeline,family=font_family)+
    annotate("point",x=year,y=0,shape=21,size=4.5,fill=cursor_colour,colour="white",stroke=1)+
    scale_x_continuous(limits=c(min(dates)-4,max(dates)+4),expand=expansion(mult=0))+
    scale_y_continuous(limits=c(-1,0.6),expand=expansion(mult=0))+
    theme_void(base_family=font_family)}

# 3.5 One frame per date ----
# Map (A) and bar (B) share rows 1-3, so the bar's base is level with the map's lower edge; the timeline (C)
# and the bar title (D) share row 4. The figure height is derived from the map extent so that the map fills
# its panel exactly (no letterboxing), which keeps the alignment of A and B.
frame_layout <- "
AAAAB
AAAAB
AAAAB
CCCCD
"
layout_widths <- c(1,1,1,1,0.45); layout_heights <- c(1,1,1,0.24)
fig_width <- 12
map_panel_width <- fig_width*sum(layout_widths[1:4])/sum(layout_widths)
fig_height <- map_panel_width*map_aspect*sum(layout_heights)/sum(layout_heights[1:3])

make_frame <- function(i){
  year <- dates[i]
  panels <- list(map_plot(i),bar_plot(i),timeline_plot(year),bar_title_plot)   # order = design letters A, B, C, D
  frame <- wrap_plots(panels,design=frame_layout,widths=layout_widths,heights=layout_heights)
  file <- file.path(frames_dir,sprintf("frame_%02d_%d.png",i,year))
  ggsave(file,frame,width=fig_width,height=fig_height,dpi=200,bg="white")
  file}

# Check one frame before building them all: open the PNG itself, not the RStudio viewer, whose small
# window squeezes the side panels when the map keeps its fixed aspect ratio
check_frame <- make_frame(2); system(paste("open",shQuote(check_frame)))

frame_index <- if(include_reference_frame) seq_along(dates) else seq_along(dates)[-1]
frame_files <- map_chr(frame_index,make_frame)

# 3.6 Assemble the animation (last frame held twice as long) ----
# GIF for the web; MP4 for Keynote, PowerPoint and QuickTime. macOS Preview lists GIF frames instead of
# playing them: open the GIF with Quick Look (space bar in Finder) or a browser, or use the MP4.
animation_frames <- c(frame_files,frame_files[length(frame_files)])
frame_size <- dim(png::readPNG(frame_files[1]))[2:1]
gifski(animation_frames,gif_file=file.path(output_dir_ch,"issf_rss_change_1918_2026.gif"),
       width=frame_size[1],height=frame_size[2],delay=frame_delay)
av_encode_video(animation_frames,output=file.path(output_dir_ch,"issf_rss_change_1918_2026.mp4"),
                framerate=1/frame_delay)

#------------------------------------------------------------------------------- STEP 4: area still available through time ----
# Mirror of STEP 3. Instead of the area lost, each frame shows the area still available to the eagles for
# terrestrial behaviours (landing), which shrinks over time:
#   available_t = cells counted as habitat and whose RSS has not fallen by more than `bar_threshold`
#                 relative to the reference scenario of STEP 3 (rss_change <= bar_threshold).
# With `habitat_quantile = 0` every Swiss cell is habitat and the available area is exactly the complement
# of the lost area of STEP 3 (total grid area minus area lost); with the "no_building" reference the 1918
# bar is already below the total. With `habitat_quantile = 0.5`, for instance, only the better half of the
# reference landscape (log-RSS above the median) counts as habitat.
# Map: same hillshade with a more opaque elevation layer; urbanisation drawn as white patches whose opacity
# grows with the RSS reduction; water and cities lightened; nests purple with a black outline. The colour
# bar is drawn by hand (bottom-right): the map colour fading to white, with a thin frame and no background.
# Reuses the objects of STEP 3 (rss_change_vs_1918, hillshade_df, elevation_df, change_df, lakes, rivers,
# river_labels, cities, nests, scale_bar, north_grob, layout, timeline_plot).

# Parameters ----
habitat_quantile <- 0             # 0 = all cells are habitat in 1918; 0.5 = cells above the 1918 median log-RSS
elevation_alpha_4 <- 0.85         # elevation layer opacity (more opaque than in STEP 3)
clr_urban <- "white"
urban_alpha_max <- 0.95           # opacity of a 100 % RSS reduction
water_alpha_4 <- 0.55             # opacity of rivers, lakes and river names
clr_nest_4 <- "#7B1E3A"; nest_alpha_4 <- 0.75   # burgundy nests, black outline, slightly transparent
clr_bar_4 <- elevation_colours[4] # brown for the bar and the timeline cursor
legend_title_4 <- sprintf("Habitat loss: reduction in relative selection\nstrength (%%) %s",reference_text)
bar_title_4 <- sprintf("Area (km²) still available\nfor landing: relative selection\nstrength reduced by less than\n%d %%, %s",round(100*bar_threshold),reference_text)
frames_dir_4 <- file.path(output_dir_ch,"animation_frames_available")
dir.create(frames_dir_4,showWarnings=FALSE)

# 4.1 Available area per date, full 100 m resolution ----
habitat_1918 <- if(habitat_quantile>0){
  log_rss_reference>=global(log_rss_reference,fun=function(v) quantile(v,habitat_quantile,na.rm=TRUE))[1,1]
}else not.na(log_rss_reference)
available_cells <- habitat_1918 & (rss_change_vs_1918<=bar_threshold)
available_by_date <- tibble(date=dates,
                            area_habitat_1918_km2=global(habitat_1918,"sum",na.rm=TRUE)[1,1]*cell_area_km2,
                            area_available_km2=as.numeric(global(available_cells,"sum",na.rm=TRUE)[,1])*cell_area_km2) %>%
  mutate(pct_of_1918_remaining=100*area_available_km2/area_habitat_1918_km2)
print(available_by_date)
write_csv(available_by_date,file.path(output_dir_ch,"area_available_by_date.csv"))

# 4.2 Map ----
urban_gradient <- alpha(clr_urban,seq(0,urban_alpha_max,length.out=100))
# Hand-drawn colour bar, bottom-right corner of the map (map units): the plateau colour fading to white,
# i.e. what the map does as a cell is lost; a thin grey frame keeps the white end visible on the page
legend_ramp <- colorRampPalette(c(elevation_colours[2],clr_urban))(100)
legend_box <- c(x0=map_xlim[2]-0.19*w,x1=map_xlim[2]-0.03*w,y0=map_ylim[1]+0.035*h,y1=map_ylim[1]+0.05*h)
legend_ticks <- tibble(x=legend_box[["x0"]]+c(0,0.5,1)*(legend_box[["x1"]]-legend_box[["x0"]]),label=c("0 %","50 %","100 %"))

map_plot_4 <- function(i){
  ggplot()+
    geom_raster(data=hillshade_df,aes(x,y,fill=shade))+
    scale_fill_gradient(low="black",high="white",guide="none")+
    new_scale_fill()+
    geom_raster(data=elevation_df,aes(x,y,fill=elevation),alpha=elevation_alpha_4)+
    scale_fill_gradientn(colours=elevation_colours,guide="none")+
    new_scale_fill()+
    geom_raster(data=change_df[[i]],aes(x,y,fill=change))+
    scale_fill_gradientn(colours=urban_gradient,limits=c(0,1),oob=scales::squish,guide="none")+
    geom_sf(data=swiss_border,fill=NA,colour="grey20",linewidth=0.45)+
    geom_sf(data=rivers,colour=alpha(clr_water_dark,water_alpha_4),linewidth=0.4)+
    geom_sf(data=lakes,fill=alpha(clr_water,water_alpha_4),colour=alpha(clr_water_dark,water_alpha_4),linewidth=0.25)+
    geom_text(data=river_labels,aes(x,y,label=river,angle=angle),family=font_family,fontface="italic",size=2.8,
              colour=alpha(clr_water_dark,water_alpha_4))+
    geom_point(data=nests,aes(x,y),shape=21,size=2,fill=alpha(clr_nest_4,nest_alpha_4),colour=alpha("black",nest_alpha_4),stroke=0.6)+
    geom_point(data=cities,aes(x,y),shape=21,size=1.9,fill=clr_text,colour="white",stroke=0.5)+
    geom_shadowtext(data=cities,aes(x,y,label=name),family=font_family,size=3,colour=clr_text,
                    bg.colour=alpha("white",0.7),bg.r=0.15,nudge_y=0.012*h,vjust=0)+
    # colour bar
    annotation_raster(matrix(legend_ramp,nrow=1),xmin=legend_box[["x0"]],xmax=legend_box[["x1"]],
                      ymin=legend_box[["y0"]],ymax=legend_box[["y1"]],interpolate=TRUE)+
    annotate("rect",xmin=legend_box[["x0"]],xmax=legend_box[["x1"]],ymin=legend_box[["y0"]],ymax=legend_box[["y1"]],
             fill=NA,colour="grey40",linewidth=0.3)+
    annotate("segment",x=legend_ticks$x,xend=legend_ticks$x,y=legend_box[["y0"]],yend=legend_box[["y0"]]-0.006*h,
             colour=clr_text,linewidth=0.4)+
    annotate("text",x=legend_ticks$x,y=legend_box[["y0"]]-0.01*h,label=legend_ticks$label,family=font_family,
             size=2.6,colour=clr_text,vjust=1)+
    annotate("text",x=legend_box[["x0"]],y=legend_box[["y1"]]+0.012*h,label=legend_title_4,family=font_family,
             size=2.6,colour=clr_text,hjust=0,vjust=0,lineheight=0.9)+
    # scale bar and north
    annotate("segment",x=scale_bar$x0,xend=scale_bar$x1,y=scale_bar$y,yend=scale_bar$y,colour=clr_marks,linewidth=0.7)+
    annotate("segment",x=c(scale_bar$x0,scale_bar$x1),xend=c(scale_bar$x0,scale_bar$x1),
             y=scale_bar$y-0.006*h,yend=scale_bar$y+0.006*h,colour=clr_marks,linewidth=0.7)+
    annotate("text",x=scale_centre,y=scale_bar$y+0.012*h,label=paste(scale_km,"km"),
             family=font_family,size=2.8,colour=clr_marks,vjust=0)+
    annotation_custom(north_grob,xmin=north_box[["xmin"]],xmax=north_box[["xmax"]],ymin=north_box[["ymin"]],ymax=north_box[["ymax"]])+
    coord_sf(xlim=map_xlim,ylim=map_ylim,expand=FALSE,crs=map_crs,datum=NA)+
    theme_void(base_family=font_family)}

# 4.3 Bar: available area, starting full in 1918 and shrinking ----
bar_max_4 <- max(available_by_date$area_available_km2)*1.12
bar_plot_4 <- function(i){
  area <- available_by_date$area_available_km2[i]
  ggplot(tibble(x="",y=area),aes(x,y))+
    geom_col(fill=clr_bar_4,colour=clr_timeline,linewidth=0.4,width=0.45)+
    annotate("text",x=1,y=area,label=paste(format(round(area),big.mark=" "),"km²"),
             vjust=-0.6,size=3.8,family=font_family,colour=clr_text)+
    scale_y_continuous(limits=c(0,bar_max_4),expand=expansion(mult=0))+
    theme_void(base_family=font_family)}

bar_title_plot_4 <- ggplot()+
  annotate("text",x=0,y=1,label=bar_title_4,family=font_family,size=2.5,hjust=0.5,vjust=1,lineheight=0.95,colour=clr_text)+
  scale_x_continuous(limits=c(-1,1))+scale_y_continuous(limits=c(0,1))+
  coord_cartesian(clip="off")+
  theme_void(base_family=font_family)

# 4.4 Frames and animation ----
make_frame_4 <- function(i){
  year <- dates[i]
  panels <- list(map_plot_4(i),bar_plot_4(i),timeline_plot(year,cursor_colour=clr_bar_4),bar_title_plot_4)
  frame <- wrap_plots(panels,design=frame_layout,widths=layout_widths,heights=layout_heights)
  file <- file.path(frames_dir_4,sprintf("frame_%02d_%d.png",i,year))
  ggsave(file,frame,width=fig_width,height=fig_height,dpi=200,bg="white")
  file}

check_frame_4 <- make_frame_4(2); system(paste("open",shQuote(check_frame_4)))

frame_files_4 <- map_chr(frame_index,make_frame_4)
animation_frames_4 <- c(frame_files_4,frame_files_4[length(frame_files_4)])
gifski(animation_frames_4,gif_file=file.path(output_dir_ch,"issf_available_area_1918_2026.gif"),
       width=frame_size[1],height=frame_size[2],delay=frame_delay)
av_encode_video(animation_frames_4,output=file.path(output_dir_ch,"issf_available_area_1918_2026.mp4"),
                framerate=1/frame_delay)


