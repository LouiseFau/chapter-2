#'-------------------------------------------------------------------------------
#' Title: Historical built-cell density in Switzerland ----
#' Author: Louise Faure
#' Date: 06.10.2026
#' **Description:** this script produces cumulative built-cell density rasters
#' for Switzerland at ten dates (1918, 1945, 1960, 1970, 1980, 1990, 2000, 2010,
#' 2015, 2026). It uses the Federal Register of Buildings and Dwellings (RegBL),
#' public data as of 30.09.2026, downloaded on 01.10.2026 from
#' https://www.housing-stat.ch/fr/data/supply/public.html
#' **Steps:**
#' (1) Data preparation:
#'     (a) filter the RegBL: remove buildings without construction period (GBAUP)
#'         and planned or unrealised buildings (GSTAT);
#'     (b) for each building, derive its year of appearance (end of its GBAUP
#'         class) and keep its year of demolition (GABBJ). Nb. demolition dates 
#'         are mostly missing before 2000.
#'     (c) project building coordinates from EPSG:2056 to EPSG:3035;
#'     (d) build a 100 m grid aligned with the other covariates and cropped to
#'         Switzerland, with cells outside Switzerland set to NA.
#' (2) Loop over the ten dates:
#'     (a) keep buildings appearing on or before the date and not demolished
#'         before it;
#'     (b) rasterize them as built (1) / not built (0) cells;
#'     (c) compute the proportion of built cells within a circular window of
#'         about 1 km2, ignoring cells outside Switzerland.
#' (3) Save one GeoTIFF per date.
#' ------------------------------------------------------------------------------

library(data.table)
library(terra)

# Paths ----
regbl_dir <- "/Users/louisefaure/Downloads/ch(1)"
boundary_file <- "/Users/louisefaure/Desktop/dossier sans titre/Rasters/historical_settlements_rasters/swissBOUNDARIES3D_1_5_LV95_LN02.gpkg"
output_dir <- "/Users/louisefaure/Desktop/dossier sans titre/Rasters/historical_settlements_rasters"
dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)

dates <- c(1918,1945,1960,1970,1980,1990,2000,2010,2015,2026)
write_options <- list(datatype="FLT4S",gdal=c("COMPRESS=ZSTD","PREDICTOR=3","TILED=YES"))


#------------------------------------------------------------------------------- STEP 1: data preparation ----
# (a) Filter the RegBL ----
bat <- fread(file.path(regbl_dir,"gebaeude_batiment_edificio.csv"),
             select=c("EGID","GKODE","GKODN","GSTAT","GBAUP","GABBJ"))
print(bat[,.N,by=GSTAT][order(GSTAT)])

# keep existing (1004), unusable (1005) and demolished (1007) buildings with a period and coordinates
bat <- bat[GSTAT %in% c(1004,1005,1007) & !is.na(GBAUP) & !is.na(GKODE) & !is.na(GKODN)]
cat("Demolished buildings without GABBJ removed:",bat[GSTAT==1007 & is.na(GABBJ),.N],"\n")
bat <- bat[!(GSTAT==1007 & is.na(GABBJ))]

# (b) Year of appearance = end of the GBAUP class ----
period_end <- data.table(GBAUP=8011:8023,year_built=c(1918,1945,1960,1970,1980,seq(1985,2015,5),2026))
bat <- period_end[bat,on="GBAUP"]

# (c) Project buildings to EPSG:3035 ----
xy <- crds(project(vect(bat,geom=c("GKODE","GKODN"),crs="EPSG:2056"),"EPSG:3035"))
bat[,`:=`(x=xy[,1],y=xy[,2])]

# (d) Swiss boundary and 100 m grid ----
layer_name <- grep("landesgebiet",vector_layers(boundary_file),ignore.case=TRUE,value=TRUE)[1]
ch <- vect(boundary_file,layer=layer_name)
icc_col <- grep("^icc$",names(ch),ignore.case=TRUE,value=TRUE)
ch <- project(aggregate(ch[values(ch)[[icc_col]]=="CH"]),"EPSG:3035") # remove Liechtenstein

grid <- rast(ext(3841300,4846200,2236600,2860700),res=100,crs="EPSG:3035")
grid <- crop(grid,ch,snap="out")
swiss_mask <- rasterize(ch,grid,field=0,touches=TRUE) # 0 inside Switzerland, NA outside

# (e) 100 m cell and 25 m sub-cell of each building, computed once ----
bat[,cell:=cellFromXY(grid,cbind(x,y))]
bat <- bat[!is.na(cell)]
bat[,sub25:=((ymax(grid)-y)%/%25)*1e6+((x-xmin(grid))%/%25)] # unique id of the 25 m sub-cell

# circular window of about 1 km2 (97 cells of 100 m, radius about 564 m)
offsets <- -6:6
radius_m <- sqrt(1e6/pi)
window_1km2 <- outer(offsets,offsets,function(r,c) ifelse((r*100)^2+(c*100)^2<=radius_m^2,1,NA))


# ------------------------------------------------------------------------------ STEP 2: loop over dates ----
mask_values <- values(swiss_mask,mat=FALSE)
summary_table <- data.table()

for(d in dates){
  # (a) buildings existing at date d
  alive <- bat[year_built<=d & (is.na(GABBJ) | GABBJ>d)]
  
  # (b) share of built 25 m sub-cells in each 100 m cell (0-1), NA outside Switzerland
  share_d <- alive[,.(share=uniqueN(sub25)/16),by=cell]
  share_d <- share_d[!is.na(mask_values[cell])]
  v <- mask_values
  v[share_d$cell] <- share_d$share
  built <- setValues(swiss_mask,v)
  
  # (c) % of built 25 m sub-cells within about 1 km2
  density <- mask(focal(built,w=window_1km2,fun="mean",na.rm=TRUE)*100,swiss_mask)
  
  # (3) save
  writeRaster(density,file.path(output_dir,sprintf("built_25m_share_1km2_100m_%d.tif",d)),
              overwrite=TRUE,wopt=write_options)
  
  summary_table <- rbind(summary_table,data.table(date=d,buildings=nrow(alive),
                                                  built_cells_100m=nrow(share_d),built_subcells_25m=sum(share_d$share*16)))
  message("Date ",d," done")
}

print(summary_table)
fwrite(summary_table,file.path(output_dir,"summary_buildings_built_cells_25m.csv"))

