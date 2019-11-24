################  GRF of Fish and Salamander Forelimbs ########################
## Code for evaluating force data collected from the Blob lab force plate #####
# Abbreviations: Vert = Vertical, ML = mediolateral, Hz = Horizontal (refers to anteroposterior direction)

# 10-14-19: Made edits, so that the AP component is not converted by a negative one.  This is because we rotate the force plate, which results
# in the sign of the AP component changing (e.g., a posterior reading would be negative for the oppossum work, but positive for the salamanders.

# Clear everything in the R workspace, so nothing gets mixed up between different trials
rm(list=ls(all=TRUE))

# Determining the date you're running these analyses
today <- Sys.Date()
SaveDate <- format(today, format="%y%m%d")

#### LOAD THE LIBRARIES ####
# if (!require(c("devtools", "signal"))) {
#   install.packages(c("devtools", "signal"), dependencies = TRUE)
#   library(c(devtools, signal))
# }
library(devtools)

# install_github("MorphoFun/kraken") # uncomment if this is the first time using kraken or if you need to update
library(kraken)

#### LOAD THE CALIBRATION FILE ####
CalibFile <- data.frame(read.csv("./dataraw/FinLimbGRFs_Calibs.csv", header=TRUE))

#### LOAD THE VIDEO INFO FILE ####
VideoFile <- data.frame(read.csv("./dataraw/FinLimbGRFs_VideoInfo.csv", header=TRUE))


#### LOADING THE DATA ####
setwd("./dataraw/force_Raw")

#setwd("./dataraw/force_Raw/problemTrials/pec")

fileList <- list.files(pattern = ".txt", full.names = FALSE)
myFiles <- lapply(fileList, FUN = read.table)
names(myFiles) <- lapply(fileList, FUN = function(x) substring(x, 1,7))

myData <- lapply(myFiles, function(x) {
  # remove unnecessary rows
  x <- x[,c(1:12)] 
  
  # add row for the sweep number
  x <- cbind(x, 1:nrow(x))  
  
  # rename column names
  colnames <- c("light_Volts", "Vert1.Volts", "Vert2.Volts", "Vert3.Volts", "Vert4.Volts", "VertSum.Volts", "ML1.Volts", "ML2.Volts", "MLSum.Volts", "Hz1.Volts", "Hz2.Volts", "HzSum.Volts", "Sweep")
  setNames(x, colnames)
  }
)

#### QUANTIFY GROUND REACTION FORCES ####

GRFanalysis <- function(myData) {
  Trial <- names(myData)
  
  ## LOOKING UP THE VIDEO INFO ##
  # Appendages listed as "Both" have both pectoral and pelvic appendage data
  VideoInfo <- VideoFile[VideoFile$File.name %in% Trial,]
  # VideoInfo$Pectoral.Start.Frame <- as.numeric(VideoInfo$Pectoral.Start.Frame)
  # VideoInfo$Pectoral.End.Frame <- as.numeric(VideoInfo$Pectoral.End.Frame)
  # VideoInfo$Pelvic.Start.Frame <- as.numeric(VideoInfo$Pelvic.Start.Frame)
  # VideoInfo$Pelvic.End.Frame <- as.numeric(VideoInfo$Pelvic.End.Frame)
  Date <- format(as.Date(VideoInfo$Date.Filmed, format = "%m/%d/%y"), format="%y%m%d")
  
  ## LOOKING UP THE CALIB INFO ##
  CalibInfo <- CalibFile[CalibFile$Date %in% Date,]
  

  #### PECTORAL APPENDAGE ####
  
  ## Selecting pectoral trials
  Pec_VideoInfo <- VideoInfo[VideoInfo$TrialsToUse == "both"|VideoInfo$TrialsToUse == "pec",]
  Pec_Data <- myData[names(myData) %in% Pec_VideoInfo$File.name]

  
  Pec_VideoInfo_Trial <- NULL
  Pec_CalibInfo_Trial <- NULL
  Pec_GRFs <- NULL
  saveAs <- NULL
  Pec_GRFs_Filtered <- NULL
  Pec_GRFs_Filtered_dataset <- NULL
  BW_N <- NULL
  Pec_GRFs_Filtered_dataset_noOverlap <- NULL
  Pec_GRFs_Filtered_dataset_noOverlap_Peak <- NULL
  pdfSave <- NULL
  GraphTitle <- NULL
  
  for (i in 1:length(Pec_Data)) {
    
    Pec_VideoInfo_Trial[[i]] <- data.frame(Pec_VideoInfo[Pec_VideoInfo$File.name == names(Pec_Data)[i],])
    Pec_CalibInfo_Trial[[i]] <- CalibInfo[as.character(CalibInfo$Filename) == as.character(Pec_VideoInfo_Trial[[i]]$Calib_Force),]

    
    ## Converting LabView output to GRFs ##
    Pec_GRFs[[i]] <- voltToForce(Pec_Data[[i]], calib = Pec_CalibInfo_Trial[[i]][3:5], zeroStart = Pec_VideoInfo_Trial[[i]]$ForceZeroStart, lightStartFrame = Pec_VideoInfo_Trial[[i]]$Light.Start, 
                                 startFrame = Pec_VideoInfo_Trial[[i]]$Pectoral.Start.Frame, endFrame = Pec_VideoInfo_Trial[[i]]$Pectoral.End.Frame, filename = Pec_VideoInfo_Trial[[i]]$File.name, BW = Pec_VideoInfo_Trial[[i]]$Body.Weight.kg)
    names(Pec_GRFs)[[i]] <- names(Pec_Data)[[i]]
    
    ## Filter the data
    saveAs[i] <- paste(Pec_VideoInfo_Trial[[i]]$File.name, "_Pec_Filtered.pdf", sep = "")
    Pec_GRFs_Filtered[[i]] <- butterFilteR(Pec_GRFs[[i]], saveAs = saveAs[[i]], saveGraph = "yes")  
    names(Pec_GRFs_Filtered)[[i]] <- names(Pec_Data)[[i]]
    
    ## Calculating angles of GRF orientation
    Pec_GRFs_Filtered_dataset[[i]] <- GRFAngles(Pec_GRFs_Filtered[[i]]$GRF0Sum_filter_interp)
    names(Pec_GRFs_Filtered_dataset)[[i]] <- names(Pec_Data)[[i]]
    
    ## Converting to units of body weight
    BW_N[i] <- Pec_VideoInfo_Trial[[i]]$Body.Weight.kg*9.8 # converting animal's body weight from kilograms to Newtons
    Pec_GRFs_Filtered_dataset[[i]]$InterpV_BW <- Pec_GRFs_Filtered_dataset[[i]]$InterpV_N/BW_N[i]
    Pec_GRFs_Filtered_dataset[[i]]$InterpML_BW <- Pec_GRFs_Filtered_dataset[[i]]$InterpML_N/BW_N[i]
    Pec_GRFs_Filtered_dataset[[i]]$InterpAP_BW <- Pec_GRFs_Filtered_dataset[[i]]$InterpAP_N/BW_N[i]
    Pec_GRFs_Filtered_dataset[[i]]$NetGRF_BW <- Pec_GRFs_Filtered_dataset[[i]]$NetGRF_N/BW_N[i]
    
    ## Remove data when pectoral and pelvic appendages overlap
    Pec_GRFs_Filtered_dataset_noOverlap[[i]] <- removeOverlaps(Pec_GRFs_Filtered_dataset[[i]], Pec_VideoInfo_Trial[[i]][2:3], Pec_VideoInfo_Trial[[i]][4:5], Pec_VideoInfo_Trial[[i]]$Filming.Rate.Hz)
    names(Pec_GRFs_Filtered_dataset_noOverlap)[[i]] <- names(Pec_Data)[[i]]
    
    ## Evaluating at peak/max net GRF
    ## Also making sure that the peak does not occur during portions where the limbs overlap on the force plate
    Pec_GRFs_Filtered_dataset_noOverlap_Peak[[i]] <- Pec_GRFs_Filtered_dataset_noOverlap[[i]][which.max(Pec_GRFs_Filtered_dataset_noOverlap[[i]]$NetGRF_BW),]
    Pec_GRFs_Filtered_dataset_noOverlap_Peak[[i]]$filename <- names(Pec_Data)[[i]]
    
    ## Saving the graphs
    pdfSave[i] <- paste(Pec_VideoInfo_Trial[[i]]$File.name, "_Pec_PostProcess.pdf", sep = "")
    pdf(pdfSave[[i]])

    par(mfrow=c(2,2), oma = c(3, 0, 2, 0))  # oma = outer margin with 2 lines above the top of the graphs

    # Vertical component of GRF graph
    plot(Pec_GRFs_Filtered_dataset[[i]]$PercentStance, Pec_GRFs_Filtered_dataset[[i]]$InterpV_BW, xlab='Percent Stance', ylab='GRF - Vertical (BW)', main='Zeroed GRF (Vertical) Force', type="l", col="black")
    lines(Pec_GRFs_Filtered_dataset_noOverlap[[i]]$PercentStance, Pec_GRFs_Filtered_dataset_noOverlap[[i]]$InterpV_BW, type="l", col="blue", lwd = 4)
    abline(v=Pec_GRFs_Filtered_dataset_noOverlap_Peak[[i]]$PercentStance, col='grey30', lty=2, lwd=2) # Plotting vertical line at Peak Net GRF

    # Mediolateral component of GRF graph
    plot(Pec_GRFs_Filtered_dataset[[i]]$PercentStance, Pec_GRFs_Filtered_dataset[[i]]$InterpML_BW, xlab='Percent Stance', ylab='GRF - Mediolateral (BW)', main='Zeroed GRF (Mediolateral) Force', type="l", col="black")
    lines(Pec_GRFs_Filtered_dataset_noOverlap[[i]]$PercentStance, Pec_GRFs_Filtered_dataset_noOverlap[[i]]$InterpML_BW, type="l", col="red", lwd = 4)
    abline(v=Pec_GRFs_Filtered_dataset_noOverlap_Peak[[i]]$PercentStance, col='grey30', lty=2, lwd=2)

    # Horizontal (Anteroposterior) component of GRF graph
    plot(Pec_GRFs_Filtered_dataset[[i]]$PercentStance, Pec_GRFs_Filtered_dataset[[i]]$InterpAP_BW, xlab='Percent Stance', ylab='GRF - Anteroposterior (BW)', main='Zeroed GRF (Anteroposterior) Force', type="l", col="black")
    lines(Pec_GRFs_Filtered_dataset_noOverlap[[i]]$PercentStance, Pec_GRFs_Filtered_dataset_noOverlap[[i]]$InterpAP_BW, type="l", col="forestgreen", lwd = 4)
    abline(v=Pec_GRFs_Filtered_dataset_noOverlap_Peak[[i]]$PercentStance, col='grey30', lty=2, lwd=2)

    # Net GRF graph
    plot(Pec_GRFs_Filtered_dataset[[i]]$PercentStance, Pec_GRFs_Filtered_dataset[[i]]$NetGRF_BW, xlab='Percent Stance', ylab='Net GRF (BW)', main='Zeroed Net GRF Force', type="l", col="black")
    lines(Pec_GRFs_Filtered_dataset_noOverlap[[i]]$PercentStance, Pec_GRFs_Filtered_dataset_noOverlap[[i]]$NetGRF_BW, type="l", col="purple", lwd = 4)
    abline(v=Pec_GRFs_Filtered_dataset_noOverlap_Peak[[i]]$PercentStance, col='grey30', lty=2, lwd=2)

    GraphTitle[i] <- pdfSave[i]
    mtext(GraphTitle[i], line=0.5, outer=TRUE)  # writes an overall title over the graphs
    mtext('Dashed grey line = % Stance for Peak Net GRF', side=1, outer=TRUE, col = 'grey30')
    mtext('Black lines indicate times with more than 1 structure on plate', side=1, line=1.5, outer=TRUE, col = 'black')
    dev.off()

  }
  
  ## Collapsing data at the peak net GRF into a single data.frame
  Pec_GRFs_Filtered_dataset_noOverlap_Peak <- do.call(rbind, Pec_GRFs_Filtered_dataset_noOverlap_Peak)
  
  
  #### PELVIC APPENDAGE ####
  
  ## Selecting pelvic trials
  Pel_VideoInfo <- VideoInfo[VideoInfo$TrialsToUse == "both"|VideoInfo$TrialsToUse == "pel",]
  Pel_Data <- myData[names(myData) %in% Pel_VideoInfo$File.name]
  
  Pel_VideoInfo_Trial <- NULL
  Pel_CalibInfo_Trial <- NULL
  Pel_GRFs <- NULL
  Pel_saveAs <- NULL
  Pel_GRFs_Filtered <- NULL
  Pel_GRFs_Filtered_dataset <- NULL
  Pel_BW_N <- NULL
  Pel_GRFs_Filtered_dataset_noOverlap <- NULL
  Pel_GRFs_Filtered_dataset_noOverlap_Peak <- NULL
  
  for (i in 1:length(Pel_Data)) {
    Pel_VideoInfo_Trial[[i]] <- data.frame(Pel_VideoInfo[Pel_VideoInfo$File.name == names(Pel_Data)[i],])
    Pel_CalibInfo_Trial[[i]] <- CalibInfo[as.character(CalibInfo$Filename) == as.character(Pel_VideoInfo_Trial[[i]]$Calib_Force),]
    
    ## Converting LabView output to GRFs ##
    Pel_GRFs[[i]] <- voltToForce(Pel_Data[[i]], Pel_CalibInfo_Trial[[i]][3:5], zeroStart = Pel_VideoInfo_Trial[[i]]$ForceZeroStart, lightStartFrame = Pel_VideoInfo_Trial[[i]]$Light.Start, 
                                 startFrame = Pel_VideoInfo_Trial[[i]]$Pelvic.Start.Frame, endFrame = Pel_VideoInfo_Trial[[i]]$Pelvic.End.Frame, filename = Pel_VideoInfo_Trial[[i]]$File.name, BW = Pel_VideoInfo_Trial[[i]]$Body.Weight.kg)
    names(Pel_GRFs)[[i]] <- names(Pel_Data)[[i]]
    
    ## Filter the data
    Pel_saveAs[i] <- paste(Pel_VideoInfo_Trial[[i]]$File.name, "_Pel_Filtered.pdf", sep = "")
    Pel_GRFs_Filtered[[i]] <- butterFilteR(Pel_GRFs[[i]], saveAs = Pel_saveAs[[i]], saveGraph = "yes")  
    names(Pel_GRFs_Filtered)[[i]] <- names(Pel_Data)[[i]]
    
    ## Calculating angles of GRF orientation
    Pel_GRFs_Filtered_dataset[[i]] <- GRFAngles(Pel_GRFs_Filtered[[i]]$GRF0Sum_filter_interp)
    names(Pel_GRFs_Filtered_dataset)[[i]] <- names(Pel_Data)[[i]]
    
    ## Converting to units of body weight
    Pel_BW_N[i] <- Pel_VideoInfo_Trial[[i]]$Body.Weight.kg*9.8 # converting animal's body weight from kilograms to Newtons
    Pel_GRFs_Filtered_dataset[[i]]$InterpV_BW <- Pel_GRFs_Filtered_dataset[[i]]$InterpV_N/BW_N[i]
    Pel_GRFs_Filtered_dataset[[i]]$InterpML_BW <- Pel_GRFs_Filtered_dataset[[i]]$InterpML_N/BW_N[i]
    Pel_GRFs_Filtered_dataset[[i]]$InterpAP_BW <- Pel_GRFs_Filtered_dataset[[i]]$InterpAP_N/BW_N[i]
    Pel_GRFs_Filtered_dataset[[i]]$NetGRF_BW <- Pel_GRFs_Filtered_dataset[[i]]$NetGRF_N/BW_N[i]
    
    ## Remove data when pectoral and pelvic appendages overlap
    Pel_GRFs_Filtered_dataset_noOverlap[[i]] <- removeOverlaps(Pel_GRFs_Filtered_dataset[[i]], Pel_VideoInfo_Trial[[i]][4:5], Pel_VideoInfo_Trial[[i]][2:3], Pel_VideoInfo_Trial[[i]]$Filming.Rate.Hz)
    names(Pel_GRFs_Filtered_dataset_noOverlap)[[i]] <- names(Pel_Data)[[i]]

    
    ## Evaluating at peak/max net GRF
    ## Also making sure that the peak does not occur during portions where the limbs overlap on the force plate
    Pel_GRFs_Filtered_dataset_noOverlap_Peak[[i]] <- Pel_GRFs_Filtered_dataset_noOverlap[[i]][which.max(Pel_GRFs_Filtered_dataset_noOverlap[[i]]$NetGRF_BW),]
    Pel_GRFs_Filtered_dataset_noOverlap_Peak[[i]]$filename <- names(Pel_Data)[[i]]
    
    ## Saving the graphs
    pdfSave[i] <- paste(Pel_VideoInfo_Trial[[i]]$File.name, "_Pel_PostProcess.pdf", sep = "")
    pdf(pdfSave[[i]])

    par(mfrow=c(2,2), oma = c(3, 0, 2, 0))  # oma = outer margin with 2 lines above the top of the graphs

    # Vertical component of GRF graph
    plot(Pel_GRFs_Filtered_dataset[[i]]$PercentStance, Pel_GRFs_Filtered_dataset[[i]]$InterpV_BW, xlab='Percent Stance', ylab='GRF - Vertical (BW)', main='Zeroed GRF (Vertical) Force', type="l", col="black")
    lines(Pel_GRFs_Filtered_dataset_noOverlap[[i]]$PercentStance, Pel_GRFs_Filtered_dataset_noOverlap[[i]]$InterpV_BW, type="l", col="blue", lwd = 4)
    abline(v=Pel_GRFs_Filtered_dataset_noOverlap_Peak[[i]]$PercentStance, col='grey30', lty=2, lwd=2) # Plotting vertical line at Peak Net GRF

    # Mediolateral component of GRF graph
    plot(Pel_GRFs_Filtered_dataset[[i]]$PercentStance, Pel_GRFs_Filtered_dataset[[i]]$InterpML_BW, xlab='Percent Stance', ylab='GRF - Mediolateral (BW)', main='Zeroed GRF (Mediolateral) Force', type="l", col="black")
    lines(Pel_GRFs_Filtered_dataset_noOverlap[[i]]$PercentStance, Pel_GRFs_Filtered_dataset_noOverlap[[i]]$InterpML_BW, type="l", col="red", lwd = 4)
    abline(v=Pel_GRFs_Filtered_dataset_noOverlap_Peak[[i]]$PercentStance, col='grey30', lty=2, lwd=2)

    # Horizontal (Anteroposterior) component of GRF graph
    plot(Pel_GRFs_Filtered_dataset[[i]]$PercentStance, Pel_GRFs_Filtered_dataset[[i]]$InterpAP_BW, xlab='Percent Stance', ylab='GRF - Anteroposterior (BW)', main='Zeroed GRF (Anteroposterior) Force', type="l", col="black")
    lines(Pel_GRFs_Filtered_dataset_noOverlap[[i]]$PercentStance, Pel_GRFs_Filtered_dataset_noOverlap[[i]]$InterpAP_BW, type="l", col="forestgreen", lwd = 4)
    abline(v=Pel_GRFs_Filtered_dataset_noOverlap_Peak[[i]]$PercentStance, col='grey30', lty=2, lwd=2)

    # Net GRF graph
    plot(Pel_GRFs_Filtered_dataset[[i]]$PercentStance, Pel_GRFs_Filtered_dataset[[i]]$NetGRF_BW, xlab='Percent Stance', ylab='Net GRF (BW)', main='Zeroed Net GRF Force', type="l", col="black")
    lines(Pel_GRFs_Filtered_dataset_noOverlap[[i]]$PercentStance, Pel_GRFs_Filtered_dataset_noOverlap[[i]]$NetGRF_BW, type="l", col="purple", lwd = 4)
    abline(v=Pel_GRFs_Filtered_dataset_noOverlap_Peak[[i]]$PercentStance, col='grey30', lty=2, lwd=2)

    GraphTitle[i] <- pdfSave[i]
    mtext(GraphTitle[i], line=0.5, outer=TRUE)  # writes an overall title over the graphs
    mtext('Dashed grey line = % Stance for Peak Net GRF', side=1, outer=TRUE, col = 'grey30')
    mtext('Black lines indicate times with more than 1 structure on plate', side=1, line=1.5, outer=TRUE, col = 'black')
    dev.off()

  }

  ## Collapsing data at the peak net GRF into a single data.frame
  Pel_GRFs_Filtered_dataset_noOverlap_Peak <- do.call(rbind, Pel_GRFs_Filtered_dataset_noOverlap_Peak)
  
  
  #### GENERATE OUTPUT #### 
  Pec_output <- list(
    Pec_GRFs = Pec_GRFs,
    Pec_GRFs_Filtered = Pec_GRFs_Filtered,
    Pec_GRFs_Filtered_dataset = Pec_GRFs_Filtered_dataset,
    Pec_GRFs_Filtered_dataset_noOverlap = Pec_GRFs_Filtered_dataset_noOverlap,
    Pec_GRFs_Filtered_PeakNet = Pec_GRFs_Filtered_dataset_noOverlap_Peak
  )
  
  Pel_output <- list(
    Pel_GRFs = Pel_GRFs,
    Pel_GRFs_Filtered = Pel_GRFs_Filtered,
    Pel_GRFs_Filtered_dataset = Pel_GRFs_Filtered_dataset,
    Pel_GRFs_Filtered_dataset_noOverlap = Pel_GRFs_Filtered_dataset_noOverlap,
    Pel_GRFs_Filtered_PeakNet = Pel_GRFs_Filtered_dataset_noOverlap_Peak
  )
  
  output <- list(
    VideoInfo = VideoInfo,
    CalibInfo = CalibInfo,
    Pectoral = Pec_output,
    Pelvic = Pel_output
  )
  
  return(output)
}

GRFs <- GRFanalysis(myData)




#### SAVING THE DATA ####
## go to the parent directory then save output in 'output' folder
setwd('..')
setwd('..')
setwd('./output')

## Pectoral data

  ## Save the dataset that was filtered and had areas of overlap excluded
  Save_FilterAll_Pec <- NULL
  for (i in 1:length(GRFs$Pectoral$Pec_GRFs_Filtered_dataset)) {
    Save_FilterAll_Pec[[i]] <- paste(names(GRFs$Pectoral$Pec_GRFs_Filtered_dataset)[[i]], "_Pec_Filtered_",SaveDate, ".csv", sep="")    
    write.table(GRFs$Pectoral$Pec_GRFs_Filtered_dataset[[i]], file = Save_FilterAll_Pec[[i]], sep =",", row.names=FALSE)
  }
  
  ## Save the dataset that was filtered and had areas of overlap excluded
  Save_FilterNoOverlap_Pec <- NULL
  for (i in 1:length(GRFs$Pectoral$Pec_GRFs_Filtered_dataset_noOverlap)) {
    Save_FilterNoOverlap_Pec[[i]] <- paste(names(GRFs$Pectoral$Pec_GRFs_Filtered_dataset_noOverlap)[[i]], "_Pec_Filtered_noOverlap_",SaveDate, ".csv", sep="")    
    write.table(GRFs$Pectoral$Pec_GRFs_Filtered_dataset_noOverlap[[i]], file = Save_FilterNoOverlap_Pec[[i]], sep =",", row.names=FALSE)
  }
  
  ## Save the data taken at the peak Net GRF
  Save_FilterNoOverlap_PeakNet_Pec <- paste("FinLimbs_Pec_Filtered_noOverlap_Peak_",SaveDate, ".csv", sep="")
  write.table(GRFs$Pectoral$Pec_GRFs_Filtered_PeakNet, file = Save_FilterNoOverlap_PeakNet_Pec, sep =",", row.names=FALSE)


## Pelvic data
  ## Save the dataset that was filtered and had areas of overlap excluded
  Save_FilterAll_Pel <- NULL
  for (i in 1:length(GRFs$Pelvic$Pel_GRFs_Filtered_dataset)) {
    Save_FilterAll_Pel[[i]] <- paste(names(GRFs$Pelvic$Pel_GRFs_Filtered_dataset)[[i]], "_Pel_Filtered_",SaveDate, ".csv", sep="")    
    write.table(GRFs$Pelvic$Pel_GRFs_Filtered_dataset[[i]], file = Save_FilterAll_Pel[[i]], sep =",", row.names=FALSE)
  }
  
  ## Save the dataset that was filtered and had areas of overlap excluded
  Save_FilterNoOverlap_Pel <- NULL
  for (i in 1:length(GRFs$Pelvic$Pel_GRFs_Filtered_dataset_noOverlap)) {
    Save_FilterNoOverlap_Pel[[i]] <- paste(names(GRFs$Pelvic$Pel_GRFs_Filtered_dataset_noOverlap)[[i]], "_Pel_Filtered_noOverlap_",SaveDate, ".csv", sep="")    
    write.table(GRFs$Pelvic$Pel_GRFs_Filtered_dataset_noOverlap[[i]], file = Save_FilterNoOverlap_Pel[[i]], sep =",", row.names=FALSE)
  }
  
  ## Save the data taken at the peak Net GRF
  Save_FilterNoOverlap_PeakNet_Pel <- paste("FinLimbs_Pel_Filtered_noOverlap_Peak_",SaveDate, ".csv", sep="")
  write.table(GRFs$Pelvic$Pel_GRFs_Filtered_PeakNet, file = Save_FilterNoOverlap_PeakNet_Pel, sep =",", row.names=FALSE)
  

  #### DISCRIMINANT FUNCTION ANALYSIS ####
