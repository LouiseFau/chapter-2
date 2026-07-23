


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                        script for extracting summary statistics from ACC data
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~






#rm(list=ls()) # clear memory


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# load packages
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(move)
library(reshape)
library(lubridate)
library(tidyverse)
library(moments)
library(stats)
library(tseries)
library(zoo)
# library(animalTrack)
library(plyr)
library(moveACC)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# import functions ====
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

source("C:/Users/lfaure7/Downloads/swisstransfer_39b015db-7422-43f4-a71e-1e9281382302/Louise_PhD/code_acc_classification/1_1__function_extracting_ss.R")
#################################################################################s



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### extract summary statistics from wild golden eagle data ====
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# make list of names in folder of working directory. following loop will be applied to all dataset listed there. make sure all of the datafiles are stored in that folder
eaglelist <- list.files(path = "C:/Users/lfaure7/Desktop/aigle non classifié", pattern = "*.csv", full.names = F)

# Function obtained from https://gitlab.com/api/v4/projects/anneks%2FmoveACC/repository/files/R%2FACCwave.R/raw?ref=master
# roll ##############
roll <- # from library "animalTrack" currently not available on cran
  function(y,z){
    roll = atan2(y,z)
    return(roll)  
  }

# pitch ###############
pitch <- # from library "animalTrack" currently not available on cran
  function(x, y, z){  
    pitch90 = -atan(x/sqrt(y^2 + z^2))  
    return(pitch90)  
  }

# loop over all data files that are mentioned in eaglelist
for (i in eaglelist){
  df <- read.csv(file = file.path("C:/Users/lfaure7/Desktop/aigle non classifié", i)) 
  df$nrsamples <- unlist(lapply(strsplit(as.character(df$eobs.accelerations.raw),' '), length))
  df <- subset(df, nrsamples == "792")
  sensitiv <- data.frame(TagID=c(5861, 7000, 7001, 5941, 7101, 7039, 7006, 7009, 5704, 7040, 7034, 5942, 7588, 7035, 5860, 6483, 6462, 5943, 5796, 6225, 7003, 7005, 7049, 5940, 7002, 7041, 7010, 5859, 4570, 7014, 7015, 6226), sensitivity="low") # tag number generation 3 logger, sensitivity automatically low
  testDFtot <- TransformRawACC(df, sensitivity.settings=sensitiv, # transform rawACC data into "g"
                                units="g", showProgress=T)
  waveDFtot <- ACCwavexyz(testDFtot,transformedData=T, showProgress=T)
  str(waveDFtot)

  write.csv(waveDFtot, row.names = FALSE, file = file.path("C:/Users/lfaure7/Desktop/aigle non classifié/ss_extract2", paste0("ss_", i)))
}
#################################################################################s
