


# clear memory
# rm(list=ls()) 
#necessary packages
library(moments)  #to calculate skewness and kurtosis of axis
library(stats)
library(tseries)
library("move")
# library("moveACC")
library("reshape")
library("tidyr")
library("lubridate")
library("dplyr")
library("ggplot2")
library("tidyverse")
library("data.table")
library(zoo)
# library(animalTrack)
library(plyr)


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

## this section is unchanged from the function "waveACC" implement in the package moveACC, version 0.1

    ACCwavexyz <- function(df,transformedData=c(FALSE,TRUE),showProgress=TRUE){ 
      whichaxes <- grep(pattern="acceleration.axes|acceleration_axes|eobs.acceleration.axes",names(df),value=T)
      if(length(whichaxes)>1){stop("Data contains more than one column referring to 'acceleration.axes'")}
      if(!transformedData){axescol <- grep(pattern="accelerations.raw|accelerations_raw|eobs.accelerations.raw",names(df),value=T)
      if(length(axescol)>1){stop("Data contains more than one column referring to 'accelerations.raw'")}}
      if(transformedData){axescol <- "accelerationTransformed"}
      indv <- grep(pattern="individual.local.identifier|individual_local_identifier",names(df),value=T)
      if(length(indv)==0){indv <- "local_identifier"}
      sampfreq <- grep(pattern="acceleration.sampling.frequency.per.axis|acceleration_sampling_frequency_per_axis|eobs.acceleration.sampling.frequency.per.axis",names(df),value=T)
      if(length(sampfreq)>1){stop("Data contains more than one column referring to 'acceleration.sampling.frequency.per.axis'")}
      df$NumberOfAxis<-(nchar(as.character(df[,whichaxes])))
      df$TotalNumberSamples <- unlist(lapply(strsplit(as.character(df[,axescol]),' '), length))
      df$BurstDurationSecs <- round((df$TotalNumberSamples/df$NumberOfAxis)/df[,sampfreq],3)
      tagid <- grep(pattern="tag.local.identifier|tag_local_identifier",names(df),value=T)
      eventId <- grep(pattern="event_id|event.id",names(df),value=T)
      bird_id <- grep(pattern="bird_id|individual.local.identifier",names(df),value=T) 
      timestamp_UTC <- grep(pattern="timestamp_UTC", names(df), value=T)
      

      df[,indv] <- as.character(df[,indv])
      df <- df[complete.cases(df[,indv]),]
      DFlist <- split(df, f = df[,indv])
      
      waveFx <- function(accDF){
        accDF$burstID <- c(0:(nrow(accDF)-1))
        accDFL <- split(accDF, 1:nrow(accDF))
        
        axisDFL <- lapply(accDFL, function(accDFLrow){
          accAxesDf <- as.numeric(unlist(strsplit(as.character(accDFLrow[,axescol]), split=" ", fixed=T)))
          ax1t<-accAxesDf[seq(1,accDFLrow$TotalNumberSamples,3)]
          ax2t<-accAxesDf[seq(2,accDFLrow$TotalNumberSamples,3)]
          ax3t<-accAxesDf[seq(3,accDFLrow$TotalNumberSamples,3)]
    
    
          # magnitude q ################ -> has been used in nathan 2012: using tri-axial acceleration
          # Definition: sqrt of the sum of squares of the three axes ( i.e. the length of the diagonal of the  x-y-z volume)
          qax1t<-numeric(length = length(ax1t))   # getting same values if i use vector or dataframe
          for (i in 1:length(qax1t)) {
            qax1t[i] <- (ax1t[i] - mean(ax1t))^2
          }
          
          qax2t<-numeric(length = length(ax2t))
          for (i in 1:length(qax2t)) {
            qax2t[i] <- (ax2t[i] - mean(ax2t))^2
          }
          
          qax3t<-numeric(length = length(ax3t))
          for (i in 1:length(qax3t)) {
            qax3t[i] <- (ax3t[i] - mean(ax3t))^2
          }
          # square root of sum of squares of sum of squares of xyz
          q <- numeric(length = length(qax1t))
          for (i in 1:length(q)) {
            q[i] <- sqrt(qax1t[i] + qax2t[i] + qax3t[i])
          }
          
          
          ##
          
          # mean  ###################

          ax1mean <- sum(ax1t)/length(ax1t)  
          ax2mean <- sum(ax2t)/length(ax2t)
          ax3mean <- sum(ax3t)/length(ax3t)
          qmean <- sum(q)/length(q)
          
          ##
          
          
          # create dataframe ###############  
          burst <- data.frame(ax1t, ax2t, ax3t, q)
          
          
          
          # standard diviation ###############
          
          ax1sd <- sd(ax1t)     
          ax2sd <- sd(ax2t)
          ax3sd <- sd(ax3t)
          qsd <- sd(q)
          
          
          ##
          # skewness ##################
          #distribution of mean in relation to median, requires package "moments"
          
          ax1skw <- skewness(ax1t)           
          ax2skw <- skewness(ax2t)
          ax3skw <- skewness(ax3t)
          qskw <- skewness(q)
          
          ##
          
          
          # kurtosis ########################
          # sharpness of peak in distribution
          
          ax1krt <- kurtosis(ax1t)           
          ax2krt <- kurtosis(ax2t)
          ax3krt <- kurtosis(ax3t)
          qkrt <- kurtosis(q)
          
          
          # maximum value #############
          
          ax1max <- max(ax1t)           
          ax2max <- max(ax2t)
          ax3max <- max(ax3t)
          qmax <- max(q)
          
          
          # minimum value ############
          
          ax1min <- min(ax1t)           
          ax2min <- min(ax2t)
          ax3min <- min(ax3t)
          qmin <- min(q)
          
          
          # median value ####
          
          ax1medi <- median(ax1t)
          ax2medi <- median(ax2t)
          ax3medi <- median(ax3t)
          qmedi <- median(q)
          
          ##
          
          # dynamic acceleration ################
          # dynamic acceleration = ODBA for each axis
          
          # mean dyn
          dynx <-abs(burst$ax1t-ax1mean)
          dyny <- abs(burst$ax2t-ax2mean)
          dynz <- abs(burst$ax3t-ax3mean)
          
          dynxmax <- max(burst$ax1t-ax1mean)
          dynymax <- max(burst$ax2t-ax2mean)
          dynzmax <- max(burst$ax3t-ax3mean)
          
          dynxmean <- mean(burst$ax1t-ax1mean)          
          dynymean <- mean(burst$ax2t-ax2mean)
          dynzmean <- mean(burst$ax3t-ax3mean)
          
          dynxmin <- min(burst$ax1t-ax1mean)          
          dynymin <- min(burst$ax2t-ax2mean)
          dynzmin <- min(burst$ax3t-ax3mean)
          
          # max dyn 
          dynxmax <- max(burst$ax1t-ax1mean)           
          dynymax <- max(burst$ax2t-ax2mean)
          dynzmax <- max(burst$ax3t-ax3mean)
          
          # median dyn
          dynxmedi <- median(burst$ax1t-ax1mean)          
          dynymedi <- median(burst$ax2t-ax2mean)
          dynzmedi <- median(burst$ax3t-ax3mean)
          
          # sd dyn 
          dynxsd <- sd(burst$ax1t-ax1mean)           
          dynysd <- sd(burst$ax2t-ax2mean)
          dynzsd <- sd(burst$ax3t-ax3mean)
          
          # ODBA #############
          #: dyn - static acceleration
          
          odba <- abs(burst$ax1t-ax1mean)+abs(burst$ax2t-ax2mean)+abs(burst$ax3t-ax3mean)
          odbamean <- mean(odba)
          odbamedi <- median(odba)
          odbasd <- sd(odba)
          
          
          
          # pairwise correlation ###########################
          xycorr <- cor(ax1t, ax2t)            
          xzcorr <- cor(ax1t, ax3t)
          yzcorr <- cor(ax2t, ax3t)
          
          
          ##
          
          # autocorrelation  ################ 
          # lag <- df$TotalNumberSamples/3
          autcorx <- acf(ax1t,  plot=F)
          autcory <- acf(ax2t, plot=F)
          autcorz <- acf(ax3t, plot=F)
          
          
          
          
          # roll ##############
          
          rollxaxis <- -ax2t # in eobs data loggers axis are different, normally x front, y side, z up. eobs: y back, x side, z down! important anjusting y and x!!!! consult eobs manual
          rollyaxis <- ax1t
          rollzaxis <- -ax3t
          
          rollanimaltrack <- roll(rollyaxis, rollzaxis) # roll how it is caluculated in the animaltrack package
          
          
          # pitch ###############
          # same needs to be done to calculate pitch -> adjust axis
          pitchanimaltrack <- pitch(rollxaxis, rollyaxis, rollzaxis)
          
          
          #trend ################
          # 
          burst$id <- 1:length(ax1t)
          trendx <- lm(ax1t~id, data=burst)$coefficients[2]
          trendy <- lm(ax2t~id, data=burst)$coefficients[2]
          trendz <- lm(ax3t~id, data=burst)$coefficients[2]
          trendq <- lm(q~id, data=burst)$coefficients[2]
          
          
          # FFT ###################
          # fast furier transorfmation of axis 
          pc <- prcomp(cbind(ax1t,ax2t,ax3t), scale. = FALSE) #principal component analysisa
          
          fftpc1 <- fft(pc$x[,1])   
          fftpc2 <- fft(pc$x[,2])
          fftpc3 <- fft(pc$x[,3])
          m1 <- Mod(fftpc1)
          m2 <- Mod(fftpc2)
          m3 <- Mod(fftpc3)
          
          peaks1 <- which.max(head(m1, round(length(m1)/2)))
          ffilt <- fftpc1
          ffilt[-peaks1] <- 0
          ffilt2 <- Re(fft(ffilt, inverse=TRUE))/264
          
          eigen1 <- eigen(cov(cbind(m1, m2, m3)))$values[1]
          
          # maximum frequency values
          psdmaxx <- max(m1)
          psdmaxy <- max(m2)
          psdmaxz <- max(m3)
          
          # minimum frequency values
          psdminx <- min(m1)
          psdminy <- min(m2)
          psdminz <- min(m3)
          
          # mean frequency values
          psdmeanx <- mean(m1)
          psdmeany <- mean(m2)
          psdmeanz <- mean(m3)
          
          
          # calculate line crossing #### 
          # how often do axis cross each other 
          
          nSignChanges <- function(x, y, z) {
            signsxy <- sign(x-y)
            signsxz <- sign(x-z)
            signsyz <- sign(y-z)
            list(xycross=sum(signsxy[-1] != signsxy[-length(x)]),
                 xzcross=sum(signsxz[-1] != signsxz[-length(z)]),
                 yzcross=sum(signsyz[-1] != signsyz[-length(z)]))
            
          }
          
          crossraw <- nSignChanges(ax1t, ax2t, ax3t)
          
          
          # include output of wanted summary statistics #########

          
          
          axisDFr <- data.frame(individualID=accDFLrow[,indv], tagID=accDFLrow[,tagid], burstID=accDFLrow$burstID, 
                                timestamp=accDFLrow$timestamp, event.id=accDFLrow[,eventId], accelerationTransformed=accDFLrow[,axescol],
                                
                                bird_id=accDFLrow[,bird_id], timestamp_UTC=accDFLrow[,timestamp_UTC], 

                                
                                beatsSec=peaks1/accDFLrow$BurstDurationSecs, 
                                amplitude= max(m1)/(length(m1)/2),  propExplPC1=summary(pc)$importance[2,1],
                                propExplPC2=summary(pc)$importance[2,2], propExplPC3=summary(pc)$importance[2,3],
                                odbaAvg=mean(odba),odbaMedian=median(odba), odbasd=odbasd,
                                varWaveWingBeat=var(ffilt2), varRestWaves=var(pc$x[,1]-ffilt2),
                                varOrigWave=var(pc$x[,1]), eigenValue1=eigen1, accAxes=accDFLrow[,whichaxes],
                                numberSamplesPerAxis=accDFLrow$TotalNumberSamples/accDFLrow$NumberOfAxis,
                                samplingFreqPerAxis=accDFLrow[,sampfreq], burstDurationSecs=accDFLrow$BurstDurationSecs,
                                
                                burstmeanx=ax1mean, burstmeany=ax2mean, burstmeanz=ax3mean, burstmeanq=qmean,
                                dynxmean=dynxmean, dynymean=dynymean, dynzmean=dynzmean,
                                dynxmax=dynxmax, dynymax=dynymax, dynzmax=dynzmax,
                                dynxmin=dynxmin, dynymin=dynymin, dynzmin=dynzmin,
                                dynxmedi=dynxmedi, dynymedi=dynymedi, dynzmedi=dynzmedi,
                                dynxsd=dynxsd, dynysd=dynysd, dynzsd=dynzsd,
                                
                                psdmaxx = psdmaxx, psdmaxy = psdmaxy, psdmaxz = psdmaxz,
                                psdminx=psdminx, psdminy=psdminy, psdminz=psdminz,
                                psdmeanx=psdmeanx, psdmeany=psdmeany, psdmeanz=psdmeanz,
                                
                                xycross=crossraw$xycross, xzcross=crossraw$xzcross, yzcross=crossraw$yzcross,
                                
                                rollanimaltrack=mean(rollanimaltrack), 
                                
                                pitchanimaltrack=mean(pitchanimaltrack),
                                
                                autcorxmean=mean(autcorx$acf), autcorymean=mean(autcory$acf), autcorzmean=mean(autcorz$acf),
                                autcorxsd=sd(autcorx$acf), autcorysd=sd(autcory$acf), autcorzsd=sd(autcorz$acf),

                                
                                skewx=ax1skw, skewy=ax2skw, skwz=ax3skw, skewq=qskw,
                                kurtx=ax1krt, kurty=ax2krt, kurtz=ax3krt, kurtq=qkrt, 
                                
                                xycorr=xycorr, xzcorr=xzcorr, yzcorr=yzcorr, 
                                trendx=trendx, trendy=trendy, trendz=trendz, trendq=trendq, 
                                
                                sdx=ax1sd, sdy=ax2sd, sdz=ax3sd, sdq=qsd,
                                xmin=ax1min, ymin=ax2min, zmin=ax3min, qmin=qmin,
                                xmax=ax1max, ymax=ax2max, zmax=ax3max, qmax=qmax,
                                medix=ax1medi, mediy=ax2medi, mediz=ax3medi, qmedi=qmedi) 
          return(axisDFr)
        })
        axisDF <- do.call("rbind",axisDFL)
        return(axisDF)
      }
      if(showProgress){waveDFL <- llply(DFlist,waveFx,.progress=progress_text(char = "."))}
      if(!showProgress){waveDFL <- llply(DFlist,waveFx)}
      
      waveDF <- do.call("rbind",waveDFL)
      return(waveDF)
    }
    