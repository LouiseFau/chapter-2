#'------------------------------------------------------------------------------
# Title : download and cleaned data from Movebank ----
#'Author : Louise Faure
#'Date: 21.07.2026
#'Adapted from Hester Bronnvik 00_access_data.R
#'Purpose: download the gps burst data from movebank for the individuals that 
#'emigrate and which behaviors have been classified by Julia Hatzl and Louise 
#'Faure (who reuse and adpated J.H random forest script). These data will be 
#'associated to the acc classified beahviors in "prepare_data.R" script.
#'------------------------------------------------------------------------------

#'libraries
library(move)
library(tidyverse)
library(lubridate)
library(move2)

#' Movebank parameters
study_id <- 282734839
speed_max_kmh <- 70

#' Movebank connection
movebank_username <- readline(
  prompt = "Nom d'utilisateur Movebank : "
)

movebank_connection <- movebank_handle(
  username = movebank_username
)

#' output dir
output_dir <- paste0("/Users/louisefaure/Desktop/dossier sans titre/", "donnees aigles gps burst")
dir.create(output_dir,recursive = TRUE,showWarnings = FALSE)

#' golden eagle dataset
emig_dat <- readRDS('/Users/louisefaure/Library/CloudStorage/OneDrive-Personnel/THESE/CHAPITRE 2/git/chapter-2/DONNEES AIGLES/emigration dates/emigration_dates_20250417.rds')
non_classified_birds_dir <- file.path(
  "/Users/louisefaure/Library/CloudStorage/OneDrive-Personnel",
  "THESE/CHAPITRE 2/git/chapter-2/DONNEES AIGLES",
  "Individus non classifies/rf_assigned"
)
julia_classified_birds_file <- read.csv("/Users/louisefaure/Library/CloudStorage/OneDrive-Personnel/THESE/CHAPITRE 2/git/chapter-2/DONNEES AIGLES/classified_acc_data/2024_01_24_alldata_allbirds_merged_rf_raw.csv")

#' general parameters
gps_dop_max <- 10 # Maximum accepted GPS dilution of precision
speed_spike_max_ms <- 50 # Max speed
max_spike_iterations <- 10L # Maximum number of successive spike-removal iterations
tz_loc <- "Europe/Zurich"



# 1. Obtain the list of individuals ----
classified_names <- julia_classified_birds_file |> transmute( individual.local.identifier = as.character(individualID))

non_classified_names <- list.files(
  non_classified_birds_dir,
  pattern = "\\.csv$",
  full.names = TRUE
) |>
  map_dfr(
    \(file) read_csv(file, show_col_types = FALSE) |>
      transmute(
        individual.local.identifier =
          as.character(individualID)))

individuals <- bind_rows(
  classified_names,
  non_classified_names
) |>
  filter(
    !is.na(individual.local.identifier),
    individual.local.identifier != ""
  ) |>
  distinct(individual.local.identifier) |>
  arrange(individual.local.identifier)

stopifnot(nrow(individuals) == 66L)


# 2. Define the temporal extraction period ----
extraction_periods <- individuals |>
  left_join(
    emig_dat |>
      transmute(
        individual.local.identifier =
          as.character(individual.local.identifier),
        extraction_start = as_date(dispersal_date)
      ) |>
      distinct(individual.local.identifier, .keep_all = TRUE),
    by = "individual.local.identifier"
  ) |>
  mutate(
    extraction_end = extraction_start + days(105)
  )

stopifnot(!anyNA(extraction_periods$extraction_start))


#' 3. Download GPS data for one individual and one period
download_individual_gps <- function(
    individual_id,
    extraction_start,
    extraction_end,
    handle,
    max_attempts = 5L
) {
  
  message(
    "Downloading : ",
    individual_id
  )
  
  extraction_start <- as_datetime(
    extraction_start,
    tz = tz_loc
  )
  
  extraction_end <- as_datetime(
    extraction_end,
    tz = tz_loc
  )
  
  attempt <- 1L
  
  
  repeat {
    
    gps_data <- tryCatch(
      
      movebank_download_study(
        study_id = study_id,
        sensor_type_id = "gps",
        individual_local_identifier = individual_id,
        timestamp_start = extraction_start,
        timestamp_end = extraction_end,
        remove_movebank_outliers = TRUE,
        handle = handle
      ),
      
      error = identity
      
    )
    
    
    if (!inherits(gps_data, "error")) {
      break
    }
    
    
    message(
      "Error for ",
      individual_id,
      " : ",
      conditionMessage(gps_data)
    )
    
    
    if (attempt >= max_attempts) {
      
      warning(
        "Skipping individual : ",
        individual_id
      )
      
      return(NULL)
      
    }
    
    
    Sys.sleep(
      min(
        30 * 2^(attempt - 1L),
        300
      )
    )
    
    
    attempt <- attempt + 1L
    
  }
  
  
  if (nrow(gps_data) == 0L) {
    
    warning(
      "No GPS data for ",
      individual_id
    )
    
    return(NULL)
    
  }
  
  
  gps_data <- gps_data[
    !sf::st_is_empty(gps_data) &
      mt_time(gps_data) >= extraction_start &
      mt_time(gps_data) < extraction_end,
  ]
  
  
  if (nrow(gps_data) == 0L) {
    
    warning(
      "No GPS data after date filtering for ",
      individual_id
    )
    
    return(NULL)
    
  }
  
  
  gps_data <- gps_data[
    order(
      mt_time(gps_data)
    ),
  ]
  
  
  message(
    "Completed : ",
    individual_id,
    " (",
    nrow(gps_data),
    " locations)"
  )
  
  
  gps_data
  
}



#' 4. Download all selected individuals
gps_bursts_raw_move2 <- extraction_periods |>
  transmute(
    individual_id =
      individual.local.identifier,
    extraction_start,
    extraction_end
  ) |>
  pmap(
    \(individual_id,
      extraction_start,
      extraction_end) {
      
      download_individual_gps(
        individual_id = individual_id,
        extraction_start = extraction_start,
        extraction_end = extraction_end,
        handle = movebank_connection
      )
      
    },
    .progress = "GPS Movebank"
  ) |>
  discard(is.null) |>
  mt_stack()



#' 5. Save raw downloaded GPS data
saveRDS(
  gps_bursts_raw_move2,
  file.path(
    output_dir,
    "gps_bursts_raw_move2.rds"
  )
)


# 6. Clean data
library(tidyverse)
library(sf)
library(geosphere)


# 1. Split GPS data by individual
gps_ls <- split(
  gps_bursts_raw_move2,
  gps_bursts_raw_move2$individual_local_identifier
)


# 2. Clean GPS data by individual
progress_bar <- txtProgressBar(
  min = 0,
  max = length(gps_ls),
  style = 3
)

clean_gps <- lapply(
  seq_along(gps_ls),
  function(i) {
    
    setTxtProgressBar(
      progress_bar,
      i
    )
    
    ind <- gps_ls[[i]]
    coords <- sf::st_coordinates(ind)
    
    ind <- ind |>
      mutate(
        id = names(gps_ls)[i],
        lon = coords[, "X"],
        lat = coords[, "Y"],
        index = row_number()
      ) |>
      as.data.frame() |>
      drop_na(lon, lat) |>
      arrange(timestamp, index)
    
    
    # Retain one location per timestamp:
    # 1. non-missing eobs_status;
    # 2. non-missing GPS DOP;
    # 3. smallest GPS DOP;
    # 4. first original row.
    ind <- ind |>
      arrange(
        timestamp,
        is.na(eobs_status),
        is.na(gps_dop),
        as.numeric(gps_dop),
        index
      ) |>
      distinct(
        timestamp,
        .keep_all = TRUE
      )
    
    
    # GPS quality filtering
    ind <- ind |>
      mutate(
        gps_dop = as.numeric(gps_dop),
        gps_dop_available = !is.na(gps_dop)
      ) |>
      select(
        -index
      ) |>
      filter(
        eobs_status == "A",
        is.na(gps_dop) | gps_dop <= gps_dop_max
      ) |>
      ungroup()
    
    
    if (nrow(ind) == 0L) {
      return(NULL)
    }
    
    
    # Distance and speed between consecutive locations
    distance_traveled <- c(
      NA_real_,
      geosphere::distHaversine(
        cbind(
          ind$lon,
          ind$lat
        )
      )
    )
    
    ind |>
      mutate(
        dist.traveled = distance_traveled,
        time.lag = as.numeric(
          difftime(
            timestamp,
            lag(timestamp),
            units = "secs"
          )
        ),
        gr.speed = if_else(
          time.lag > 0,
          dist.traveled / time.lag,
          NA_real_
        )
      ) |>
      as.data.frame()
  }
)

close(progress_bar)


# 3. Assign individual names before removing empty elements
names(clean_gps) <- names(gps_ls)


# 4. Summarize the cleaning process
cleaning_summary <- tibble(
  individual.local.identifier = names(gps_ls),
  number.raw = map_int(
    gps_ls,
    NROW
  ),
  number.clean = map_int(
    clean_gps,
    NROW
  )
) |>
  mutate(
    number.removed = number.raw - number.clean,
    proportion.retained = number.clean / number.raw,
    individual.lost = number.clean == 0L
  )


# 5. Identify individuals lost during cleaning
lost_individuals <- cleaning_summary |>
  filter(individual.lost) |>
  select(
    individual.local.identifier,
    number.raw,
    number.clean,
    number.removed
  )


# 6. Remove empty individuals while retaining their names
clean_gps <- compact(clean_gps)


# 7. Global summary
cleaning_overview <- cleaning_summary |>
  summarise(
    number.initial.individuals = n(),
    number.retained.individuals =
      sum(!individual.lost),
    number.lost.individuals =
      sum(individual.lost),
    number.raw.locations =
      sum(number.raw),
    number.clean.locations =
      sum(number.clean),
    number.removed.locations =
      sum(number.removed)
  )

print(cleaning_overview)
print(lost_individuals)
# 8. Print the total number of removed locations
print(
  cleaning_overview$number.removed.locations
)
# 9. Print removed locations by individual
cleaning_summary |>
  arrange(
    desc(number.removed)
  ) |>
  select(
    individual.local.identifier,
    number.raw,
    number.clean,
    number.removed,
    proportion.retained
  ) |>
  print(
    n = Inf
  )


# 7. Global summary
cleaning_overview <- cleaning_summary |>
  summarise(
    number.initial.individuals = n(),
    number.retained.individuals =
      sum(!individual.lost),
    number.lost.individuals =
      sum(individual.lost),
    number.raw.locations =
      sum(number.raw),
    number.clean.locations =
      sum(number.clean),
    number.removed.locations =
      sum(number.removed)
  )

print(cleaning_overview)


# 8. Print the total number of removed locations
print(
  cleaning_overview$number.removed.locations
)


# 9. Print removed locations by individual
cleaning_summary |>
  arrange(
    desc(number.removed)
  ) |>
  select(
    individual.local.identifier,
    number.raw,
    number.clean,
    number.removed,
    proportion.retained
  ) |>
  print(
    n = Inf
  )


# 10. Print lost individuals
print(
  lost_individuals,
  n = Inf
)

# 8. Save cleaned data and cleaning summary
saveRDS(
  clean_gps,
  file.path(
    output_dir,
    "gps_clean_by_individual.rds"
  )
)
