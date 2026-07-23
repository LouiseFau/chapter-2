# rm(list=ls()) # clear memory




# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                        script for random forest behavioural classification
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


library("randomForest")
library("dplyr")
library("caret")





# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### # # Data import ====
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

df1 <- read.csv("C:/Users/lfaure7/Downloads/swisstransfer_39b015db-7422-43f4-a71e-1e9281382302/Louise_PhD/data_acc_classification/2__random_forest/validation_month_wildgoldeneagles.csv", sep=",",row.names = NULL)
df1$behaviour <- as.factor(df1$behaviour)
str(df1)


# remove the additional columns of the dataset which contain information about the individual, tag, burst etc. This ensures that the dataset afterwards only contains the behavioural information (behaviour) and the summary statistics. For detailed description of the abbrevations please refer to the vignette of the moveACC pagage and the supplementary material of this paper (table S1)


df1 <- df1 %>% dplyr::select(-c(individualID, tagID, burstID, timestamp, event.id))
str(df1)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### Splitting training testing ====
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# splitting dataset into training and testing with maintaining ratio for behaviours
set.seed(24)
summary(df1$behaviour.num)
trainindexland <- caret::createDataPartition(df1$behaviour, p=.8, list= F, times = 1)
trainDatafin <- df1[trainindexland,]
testDatafin <- df1[-trainindexland,]
str(trainDatafin)








# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Random Forest fitting ====
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

###############################################################################################################d
tunegrid <- expand.grid(.mtry=c(3:15))
metric <- "Accuracy"
summary(trainDatafin)

fitControl <- trainControl(
  method="repeatedcv", number=10,
  repeats = 3,
  search="grid",
  savePredictions = 'final',       # saves predictions for optimal tuning parameter
  classProbs = T,                  # should class probabilities be returned
  summaryFunction=multiClassSummary,  # results summary function
  sampling="up"
)

set.seed(55)
model_fitted8behuni2 = train(behaviour~., data=trainDatafin, method="rf", metric=metric, tuneGrid=tunegrid, trControl=fitControl, ntree=2500, importance=TRUE)
model_fitted8behuni2
plot(model_fitted8behuni2)

impvar <- varImp(model_fitted8behuni2, scale=FALSE, n.var=58)
varImpPlot(model_fitted8behuni2$finalModel, n.var=58)
varimpdf
plot(impvar)
importance(model_fitted8behuni2$finalModel, type=1)


set.seed(58)
predictedmodel_fitted8behuni <- predict(model_fitted8behuni2, testDatafin)
confusionMatrix(reference = testDatafin$behaviour, data = predictedmodel_fitted8behuni)

levels(testDatafin$behaviour)
table(testDatafin$behaviour)
table(predictedmodel_fitted8behuni)


###########################################################d






# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Prediction of wild golden eagle data ====
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~




# make list of names in folder of working directory. following loop will be applied to all dataset listed there. make sure all of the datafiles are stored in that folder, not any others

eaglelist <- list.files(path = "C:/Users/lfaure7/Desktop/aigle non classifié/ss_extract2", pattern = "*.csv", full.names = F) # import extracted summary data


# define colums to be selected for prediction. need to be the same columns of summary statistics as when training the random forest
for.columns <- c("varWaveWingBeat", "varRestWaves", "varOrigWave", "eigenValue1", "beatsSec", "amplitude",
                 "propExplPC1", "propExplPC2", "propExplPC3", "odbaAvg", "odbasd",
                 "burstmeanx", "burstmeany", "burstmeanz", "sdx", "sdy", "sdz",
                 "skewx", "skewy", "skwz", "kurtx", "kurty", "kurtz",
                 "xmax", "ymax", "zmax", "qmax", "xmin", "ymin", "zmin",
                 "xycorr", "xzcorr", "yzcorr", "trendx", "trendy", "trendz",
                 "dynxmax", "dynymax", "dynzmax", "xycross", "xzcross", "yzcross",
                 "rollanimaltrack", "pitchanimaltrack",
                 "autcorxmean", "autcorymean", "autcorzmean",
                 "autcorxsd", "autcorysd", "autcorzsd",
                 "dynxmean", "dynymean", "dynzmean",
                 "dynxmin", "dynymin", "dynzmin",
                 "dynxsd", "dynysd", "dynzsd")

for.supp <- c("individualID", "tagID", "burstID", "timestamp", "event.id", "odbaAvg", "rollanimaltrack", "pitchanimaltrack","burstmeanx", "burstmeany", "burstmeanz")

# loop over all data files that are mentioned in eaglelist
for (i in eaglelist){
  for.dat <- read.csv(file = file.path("C:/Users/lfaure7/Desktop/aigle non classifié/ss_extract2", i)) # read in csv 
  dfpp <- for.dat [,for.columns] # make subset of dataframe according to the vector defined outside the loop
  for.dat <- for.dat [,for.supp] # make subset of dataframe according to the vector defined outside the loop

 
  uni_rf8fitted_pre <- predict(model_fitted8behuni2, dfpp, type="prob")
  for.dat$pro_rf8fitted <- apply(uni_rf8fitted_pre,1,function(x)max(x))
  for.dat$rf8fitted <- predict(model_fitted8behuni2, dfpp) # get class

  # save new file in different folder
  write.csv(for.dat, row.names = FALSE, file = file.path("C:/Users/lfaure7/Desktop/aigle non classifié/rf_assigned2", paste0("rf_", i)))
}



# autre boucle chat

out.dir <- "C:/Users/lfaure7/Desktop/aigle non classifié/rf_assigned2"

if (!dir.exists(out.dir)) {
  dir.create(out.dir, recursive = TRUE)
}

for (i in eaglelist) {
  
  message("Processing: ", i)
  
  for.dat.full <- read.csv(
    file = file.path("C:/Users/lfaure7/Desktop/aigle non classifié/ss_extract2", i)
  )
  
  # Vérifier les colonnes manquantes
  missing.cols <- setdiff(for.columns, names(for.dat.full))
  
  if (length(missing.cols) > 0) {
    stop(
      paste0(
        "Dans le fichier ", i, ", colonnes manquantes : ",
        paste(missing.cols, collapse = ", ")
      )
    )
  }
  
  dfpp <- for.dat.full[, for.columns, drop = FALSE]
  for.dat <- for.dat.full[, for.supp, drop = FALSE]
  
  # Identifier les lignes utilisables
  ok.rows <- complete.cases(dfpp) & rowSums(!is.finite(as.matrix(dfpp))) == 0
  
  # Diagnostic si certaines lignes sont exclues
  if (any(!ok.rows)) {
    bad.rows <- which(!ok.rows)
    
    message(
      "Attention : ", length(bad.rows),
      " ligne(s) avec NA/NaN/Inf dans ", i,
      ". Elles seront gardées dans le fichier de sortie, mais sans prédiction."
    )
    
    bad.cells <- which(
      is.na(dfpp) | !is.finite(as.matrix(dfpp)),
      arr.ind = TRUE
    )
    
    print(
      data.frame(
        file = i,
        row = bad.cells[, 1],
        column = names(dfpp)[bad.cells[, 2]]
      )
    )
  }
  
  # Initialiser les colonnes de sortie avec NA
  for.dat$pro_rf8fitted <- NA_real_
  for.dat$rf8fitted <- NA_character_
  
  # Prédire seulement les lignes complètes
  uni_rf8fitted_pre <- predict(
    model_fitted8behuni2,
    dfpp[ok.rows, , drop = FALSE],
    type = "prob"
  )
  
  for.dat$pro_rf8fitted[ok.rows] <- apply(
    uni_rf8fitted_pre,
    1,
    max
  )
  
  for.dat$rf8fitted[ok.rows] <- as.character(
    predict(
      model_fitted8behuni2,
      dfpp[ok.rows, , drop = FALSE]
    )
  )
  
  write.csv(
    for.dat,
    row.names = FALSE,
    file = file.path(out.dir, paste0("rf_", i))
  )
}


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### Validation month prediction ====
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

val <- read.csv("C:/Users/lfaure7/Downloads/swisstransfer_39b015db-7422-43f4-a71e-1e9281382302/Louise_PhD/data_acc_classification/2__random_forest/validation_month_wildgoldeneagles.csv")

str(val)
val$individualID <- as.factor(val$individualID)
val$behaviour <- as.factor(val$behaviour)


# make subset of this dataset of only summary statistics
dfrf_val <- val %>% 
  dplyr::select(varWaveWingBeat, varRestWaves, varOrigWave, eigenValue1, beatsSec, amplitude, propExplPC1, propExplPC2, propExplPC3, odbaAvg, odbasd, burstmeanx, burstmeany, burstmeanz, sdx, sdy, sdz, skewx, skewy, skwz, kurtx, kurty, kurtz, xmax, ymax, zmax, qmax, xmin, ymin, zmin,
                xycorr, xzcorr, yzcorr, trendx, trendy, trendz, dynxmax, dynymax, dynzmax,  xycross, xzcross, yzcross, rollanimaltrack, pitchanimaltrack, autcorxmean, autcorymean, autcorzmean, autcorxsd, autcorysd, autcorzsd, dynxmean, dynymean, dynzmean, dynxmin, dynymin, dynzmin, dynxsd, dynysd, dynzsd)


df_supp_val <- val %>%
  dplyr::select(individualID, tagID, burstID, timestamp, event.id, odbaAvg, rollanimaltrack, pitchanimaltrack, burstmeanx, burstmeany, burstmeanz)
dput(colnames(df_supp_val))



dfpp<- dfrf_val

set.seed(123)



df_supp_val$rf8fitted <- predict(model_fitted8behuni2, dfpp)

uni_rf8fitted_pre <- predict(model_fitted8behuni2, dfpp, type="prob")
df_supp_val$rf8fitted_pre <- apply(uni_rf8fitted_pre,1,function(x)max(x))

uni_rf8fitted_class <- predict(model_fitted8behuni2, dfpp)


str(df_supp_val)
df_supp_val$ref_behaviour <- val$behaviour
levels(df_supp_val$ref_behaviour)
levels(df_supp_val$rf8fitted)


table(df_supp_val$rf8fitted)
table(df_supp_val$ref_behaviour)

confusionMatrix(reference = df_supp_val$ref_behaviour, data = df_supp_val$rf8fitted)

