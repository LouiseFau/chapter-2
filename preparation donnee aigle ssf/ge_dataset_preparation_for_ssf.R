#' ---
#' title: "Preparing Golden Eagle dataset for ssf"
#' author: "Louise Faure"
#' date: 02.06.2026
#' details: (i) select all the golden eagles gps values for the first 15 weeks of the dispersal phase, 
#' (ii) define a threshold value to keep the location close to the group
#' (iii) identify the time at which each location become independent    
#' ---   
#' 
#' #library
library(dplyr)





#' # STEP 0 : load the data ----------------------------------------------------
#' # set output directory 
output_dir <- "C:/Users/lfaure7/Documents/git/chapter-2/preparation donnee aigle ssf/resultats-intermediaire"

# Emigration date
emig_dates <- readRDS("C:/Users/lfaure7/Documents/git/chapter-2/preparation donnee aigle ssf/donnees/emigration dates/emigration_dates_20250417.rds")

# golden eagles with a date of dispersal
rds_dir <- "C:/Users/lfaure7/Documents/git/chapter-2/preparation donnee aigle ssf/donnees/no_burst_GE/no_burst_GE"
rds_files <- tools::file_path_sans_ext(list.files(rds_dir, pattern = "\\.rds$", ignore.case = TRUE))
rds_ids   <- as.numeric(gsub("_gpsNoDup_moveObj", "", rds_files))

# Keep only individuals with a confirmed dispersal date AND an RDS file
emig_dates_filtered <- emig_dates[emig_dates$did_disperse == TRUE &
                                    emig_dates$individual.id %in% rds_ids, ]
rds_files_filtered  <- rds_files[rds_ids %in% emig_dates_filtered$individual.id]

cat("Individuals with confirmed dispersal date and RDS file:", nrow(emig_dates_filtered), "\n")

#' # STEP 1 : filter the golden eagles data ------------------------------------
#' (i) keep only GPS locations within the first 15 weeks after emigration date
#' (ii) inspect speed and height to define filtering thresholds

#' # STEP 1 : filter the golden eagles data ------------------------------------

dispersal_data <- lapply(rds_files_filtered, function(f) {
  
  id  <- as.numeric(gsub("_gpsNoDup_moveObj", "", f))
  obj <- readRDS(file.path(rds_dir, paste0(f, ".rds")))
  df  <- as.data.frame(obj)
  
  # --- Filtrage des colonnes DÈS la lecture (comme ton ancien code) ---
  df <- df[, c("individual.local.identifier", "timestamp",
               "location.long", "location.lat",
               "height.above.ellipsoid", "ground.speed", 
               "eobs.speed.accuracy.estimate", "eobs.horizontal.accuracy.estimate", 
               "gps.dop", "height.above.ellipsoid", "vert.speed",
               "turn.angle", "step.length", "gr.speed")]
  
  df$timestamp <- as.POSIXct(df$timestamp, tz = "UTC")
  emig_dt <- emig_dates_filtered$dispersal_date[emig_dates_filtered$individual.id == id]
  
  # Garde les 15 premières semaines (105 jours) après émigration
  df <- df[df$timestamp >= emig_dt &
             df$timestamp <= emig_dt + 105 * 24 * 3600, ]
  
  df$individual.id   <- id
  df$dispersal_date  <- emig_dt
  df$days_since_emig <- as.numeric(difftime(df$timestamp, emig_dt, units = "days"))
  
  df
})

dispersal_data <- bind_rows(dispersal_data)
rownames(dispersal_data) <- NULL

