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
if (!require(c("devtools", "signal"))) {
  install.packages(c("devtools", "signal"), dependencies = TRUE)
  library(c(devtools, signal))
}
library(devtools)
install_github("MorphoFun/kraken")
library(kraken)

#### LOAD THE CALIBRATION FILE ####
CalibFile <- data.frame(read.csv("./dataraw/FinLimbGRFs_Calibs.csv", header=TRUE))

#### LOAD THE VIDEO INFO FILE ####
VideoFile <- data.frame(read.csv("./dataraw/FinLimbGRFs_VideoInfo.csv", header=TRUE))


#### LOADING THE DATA ####
setwd("./dataraw")

myFile <- file.choose()
myData <- read.table(myFile, header=FALSE)
myData <- myData[,c(1:12)] # Last 4 columns/channels were unused, so only subsetting what I need
names(myData) <- c("light_Volts", "Vert1.Volts", "Vert2.Volts", "Vert3.Volts", "Vert4.Volts", "VertSum.Volts", "ML1.Volts", "ML2.Volts", "MLSum.Volts", "Hz1.Volts", "Hz2.Volts", "HzSum.Volts")
myData$Sweep <- 1:nrow(myData)

# Determining trial name
Trial <-  substring(myFile, nchar(myFile)-10, nchar(myFile)-4)


#### LOOKING UP THE VIDEO INFO ####
# Appendages listed as "Both" have both pectoral and pelvic appendage data
VideoInfo <- VideoFile[VideoFile$File.name %in% Trial,]
VideoInfo$Pectoral.Start.Frame <- as.numeric(VideoInfo$Pectoral.Start.Frame)
VideoInfo$Pectoral.End.Frame <- as.numeric(VideoInfo$Pectoral.End.Frame)
VideoInfo$Pelvic.Start.Frame <- as.numeric(VideoInfo$Pelvic.Start.Frame)
VideoInfo$Pelvic.End.Frame <- as.numeric(VideoInfo$Pelvic.End.Frame)
Date <- format(as.Date(VideoInfo$Date.Filmed, format = "%m/%d/%y"), format="%y%m%d")


#### LOOKING UP THE CALIB INFO ####
CalibInfo <- CalibFile[CalibFile$Date %in% Date,]


#### CONVERTING LABVIEW OUTPUT TO GRFS ####

if (!VideoInfo$Appendages == 'Pelvic') Pec_GRFs <- voltToForce(myData, CalibInfo[,3:5], zeroStart = VideoInfo$Pectoral.Start.Frame, lightStartFrame = VideoInfo$Light.Start, startFrame = VideoInfo$Pectoral.Start.Frame, endFrame = VideoInfo$Pectoral.End.Frame, filename = Trial, BW = VideoInfo$Body.Weight.kg)
if (!VideoInfo$Appendages == 'Pectoral') Pel_GRFs <- voltToForce(myData, CalibInfo[,3:5], VideoInfo$Pectoral.Start.Frame, VideoInfo$Light.Start, VideoInfo$Pelvic.Start.Frame, VideoInfo$Pelvic.End.Frame, filename = Trial, BW = VideoInfo$Body.Weight.kg)


#### FILTERING THE DATA ####
if (!VideoInfo$Appendages == 'Pelvic') {
  saveAs <- paste(Trial, "_Pec_Filtered.pdf", sep = "")
  Pec_GRFs_Filtered <- butterFilteR(Pec_GRFs, saveAs = saveAs, saveGraph = "yes")
}

if (!VideoInfo$Appendages == 'Pectoral') {
  saveAs <- paste(Trial, "_Pel_Filtered.pdf", sep = "")
  Pel_GRFs_Filtered <- butterFilteR(Pel_GRFs, saveAs = saveAs, saveGraph = "yes")
} 


#### CALCULATING ANGLES OF GRF ORIENTATION ####
if (!VideoInfo$Appendages == 'Pelvic') Pec_GRFs_Filtered_dataset <- GRFAngles(Pec_GRFs_Filtered$GRF0Sum_filter_interp)
if (!VideoInfo$Appendages == 'Pectoral') Pel_GRFs_Filtered_dataset <- GRFAngles(Pel_GRFs_Filtered$GRF0Sum_filter_interp)

#### CONVERTING TO UNITS OF BODY WEIGHT ####
BW_N <- VideoInfo$Body.Weight.kg*9.8 # converting animal's body weight from kilograms to Newtons
if (!VideoInfo$Appendages == 'Pelvic') {
  Pec_GRFs_Filtered_dataset$InterpV_BW <- Pec_GRFs_Filtered_dataset$InterpV_N/BW_N
  Pec_GRFs_Filtered_dataset$InterpML_BW <- Pec_GRFs_Filtered_dataset$InterpML_N/BW_N
  Pec_GRFs_Filtered_dataset$InterpAP_BW <- Pec_GRFs_Filtered_dataset$InterpAP_N/BW_N
  Pec_GRFs_Filtered_dataset$NetGRF_BW <- Pec_GRFs_Filtered_dataset$NetGRF_N/BW_N
}

if (!VideoInfo$Appendages == 'Pectoral') {
  Pel_GRFs_Filtered_dataset$InterpV_BW <- Pel_GRFs_Filtered_dataset$InterpV_N/BW_N
  Pel_GRFs_Filtered_dataset$InterpML_BW <- Pel_GRFs_Filtered_dataset$InterpML_N/BW_N
  Pel_GRFs_Filtered_dataset$InterpAP_BW <- Pel_GRFs_Filtered_dataset$InterpAP_N/BW_N
  Pel_GRFs_Filtered_dataset$NetGRF_BW <- Pel_GRFs_Filtered_dataset$NetGRF_N/BW_N
}


#### REMOVING AREAS OF OVERLAP ####
## This is done to exclude the parts of stance where two appendages were contacting the plate at the same time

Pec_GRFs_Filtered_dataset_noOverlap <- removeOverlaps(Pec_GRFs_Filtered_dataset, VideoInfo[,2:3], VideoInfo[,4:5], VideoInfo$Filming.Rate.Hz)
Pel_GRFs_Filtered_dataset_noOverlap <- removeOverlaps(Pel_GRFs_Filtered_dataset, VideoInfo[,4:5], VideoInfo[,2:3], VideoInfo$Filming.Rate.Hz)


#### DETERMINING VALUES AT THE PEAK NET GRF ####
## Evaluating at peak/max net GRF
## Also making sure that the peak does not occur during portions where the limbs overlap on the force plate
Pec_GRFs_Filtered_dataset_noOverlap_Peak <- Pec_GRFs_Filtered_dataset_noOverlap[which.max(Pec_GRFs_Filtered_dataset_noOverlap$NetGRF_BW),]
Pel_GRFs_Filtered_dataset_noOverlap_Peak <- Pel_GRFs_Filtered_dataset_noOverlap[which.max(Pel_GRFs_Filtered_dataset_noOverlap$NetGRF_BW),]

#### SAVING THE GRAPHS ####

### Graphing post-processed data
## Pectoral appendage
if (!VideoInfo$Appendages == 'Pelvic') {
  pdfSave <- paste(Trial, "_Pec_PostProcess.pdf", sep = "")
  pdf(pdfSave)
  
  par(mfrow=c(2,2), oma = c(3, 0, 2, 0))  # oma = outer margin with 2 lines above the top of the graphs
  
  # Vertical component of GRF graph
  plot(Pec_GRFs_Filtered_dataset$PercentStance, Pec_GRFs_Filtered_dataset$InterpV_BW, xlab='Percent Stance', ylab='GRF - Vertical (BW)', main='Zeroed GRF (Vertical) Force', type="l", col="black")
  if (VideoInfo$Appendages == "Both")
  {
    lines(Pec_GRFs_Filtered_dataset_noOverlap$PercentStance, Pec_GRFs_Filtered_dataset_noOverlap$InterpV_BW, type="l", col="blue", lwd = 4)
  }
  abline(v=Pec_GRFs_Filtered_dataset_noOverlap_Peak$PercentStance, col='grey30', lty=2, lwd=2) # Plotting vertical line at Peak Net GRF
  
  # Mediolateral component of GRF graph
  plot(Pec_GRFs_Filtered_dataset$PercentStance, Pec_GRFs_Filtered_dataset$InterpML_BW, xlab='Percent Stance', ylab='GRF - Mediolateral (BW)', main='Zeroed GRF (Mediolateral) Force', type="l", col="black")
  if (VideoInfo$Appendages == "Both")
  {
    lines(Pec_GRFs_Filtered_dataset_noOverlap$PercentStance, Pec_GRFs_Filtered_dataset_noOverlap$InterpML_BW, type="l", col="red", lwd = 4)
  }
  abline(v=Pec_GRFs_Filtered_dataset_noOverlap_Peak$PercentStance, col='grey30', lty=2, lwd=2)
  
  # Horizontal (Anteroposterior) component of GRF graph
  plot(Pec_GRFs_Filtered_dataset$PercentStance, Pec_GRFs_Filtered_dataset$InterpAP_BW, xlab='Percent Stance', ylab='GRF - Anteroposterior (BW)', main='Zeroed GRF (Anteroposterior) Force', type="l", col="black")
  if (VideoInfo$Appendages == "Both")
  {
    lines(Pec_GRFs_Filtered_dataset_noOverlap$PercentStance, Pec_GRFs_Filtered_dataset_noOverlap$InterpAP_BW, type="l", col="forestgreen", lwd = 4)
  }
  abline(v=Pec_GRFs_Filtered_dataset_noOverlap_Peak$PercentStance, col='grey30', lty=2, lwd=2)
  
  # Net GRF graph
  plot(Pec_GRFs_Filtered_dataset$PercentStance, Pec_GRFs_Filtered_dataset$NetGRF_BW, xlab='Percent Stance', ylab='Net GRF (BW)', main='Zeroed Net GRF Force', type="l", col="black")
  if (VideoInfo$Appendages == "Both")
  {
    lines(Pec_GRFs_Filtered_dataset_noOverlap$PercentStance, Pec_GRFs_Filtered_dataset_noOverlap$NetGRF_BW, type="l", col="purple", lwd = 4)
  }
  abline(v=Pec_GRFs_Filtered_dataset_noOverlap_Peak$PercentStance, col='grey30', lty=2, lwd=2)
  GraphTitle <- pdfSave
  mtext(GraphTitle, line=0.5, outer=TRUE)  # writes an overall title over the graphs
  mtext('Dashed pink line = % Stance for Peak Net GRF', side=1, outer=TRUE, col = 'magenta')
  mtext('Black lines indicate times with more than 1 structure on plate', side=1, line=1.5, outer=TRUE, col = 'black')
  dev.off()
}

## Pelvic Appendage
if (!VideoInfo$Appendages == 'Pectoral') {
  pdfSave <- paste(Trial, "_Pel_PostProcess.pdf", sep = "")
  pdf(pdfSave)
  
  par(mfrow=c(2,2), oma = c(3, 0, 2, 0))  # oma = outer margin with 2 lines above the top of the graphs
  
  # Vertical component of GRF graph
  plot(Pel_GRFs_Filtered_dataset$PercentStance, Pel_GRFs_Filtered_dataset$InterpV_BW, xlab='Percent Stance', ylab='GRF - Vertical (BW)', main='Zeroed GRF (Vertical) Force', type="l", col="black")
  if (VideoInfo$Appendages == "Both")
  {
    lines(Pel_GRFs_Filtered_dataset_noOverlap$PercentStance, Pel_GRFs_Filtered_dataset_noOverlap$InterpV_BW, type="l", col="blue", lwd = 4)
  }
  abline(v=Pel_GRFs_Filtered_dataset_noOverlap_Peak$PercentStance, col='grey30', lty=2, lwd=2) # Plotting vertical line at Peak Net GRF
  
  # Mediolateral component of GRF graph
  plot(Pel_GRFs_Filtered_dataset$PercentStance, Pel_GRFs_Filtered_dataset$InterpML_BW, xlab='Percent Stance', ylab='GRF - Mediolateral (BW)', main='Zeroed GRF (Mediolateral) Force', type="l", col="black")
  if (VideoInfo$Appendages == "Both")
  {
    lines(Pel_GRFs_Filtered_dataset_noOverlap$PercentStance, Pel_GRFs_Filtered_dataset_noOverlap$InterpML_BW, type="l", col="red", lwd = 4)
  }
  abline(v=Pel_GRFs_Filtered_dataset_noOverlap_Peak$PercentStance, col='grey30', lty=2, lwd=2)
  
  # Horizontal (Anteroposterior) component of GRF graph
  plot(Pel_GRFs_Filtered_dataset$PercentStance, Pel_GRFs_Filtered_dataset$InterpAP_BW, xlab='Percent Stance', ylab='GRF - Anteroposterior (BW)', main='Zeroed GRF (Anteroposterior) Force', type="l", col="black")
  if (VideoInfo$Appendages == "Both")
  {
    lines(Pel_GRFs_Filtered_dataset_noOverlap$PercentStance, Pel_GRFs_Filtered_dataset_noOverlap$InterpAP_BW, type="l", col="forestgreen", lwd = 4)
  }
  abline(v=Pel_GRFs_Filtered_dataset_noOverlap_Peak$PercentStance, col='grey30', lty=2, lwd=2)
  
  # Net GRF graph
  plot(Pel_GRFs_Filtered_dataset$PercentStance, Pel_GRFs_Filtered_dataset$NetGRF_BW, xlab='Percent Stance', ylab='Net GRF (BW)', main='Zeroed Net GRF Force', type="l", col="black")
  if (VideoInfo$Appendages == "Both")
  {
    lines(Pel_GRFs_Filtered_dataset_noOverlap$PercentStance, Pel_GRFs_Filtered_dataset_noOverlap$NetGRF_BW, type="l", col="purple", lwd = 4)
  }
  abline(v=Pel_GRFs_Filtered_dataset_noOverlap_Peak$PercentStance, col='magenta', lty=2, lwd=2)
  GraphTitle <- pdfSave
  mtext(GraphTitle, line=0.5, outer=TRUE)  # writes an overall title over the graphs
  mtext('Dashed pink line = % Stance for Peak Net GRF', side=1, outer=TRUE, col = 'magenta')
  mtext('Black lines indicate times with more than 1 structure on plate', side=1, line=1.5, outer=TRUE, col = 'black')
  dev.off()
}


#### SAVING THE DATA ####

## Pectoral data
if (!VideoInfo$Appendages == 'Pelvic') {
  ## Save the dataset that was filtered and had areas of overlap excluded
  Save_FilterAll_Pec <- paste(Trial,"_Pec_Filtered_All_",SaveDate, ".csv", sep="")
  write.table(Pec_GRFs_Filtered_dataset, file = Save_FilterAll_Pec, sep =",", row.names=FALSE)
  
  ## Save the dataset that was filtered and had areas of overlap excluded
  Save_FilterNoOverlap_Pec <- paste(Trial,"_Pec_Filtered_noOverlap_",SaveDate, ".csv", sep="")
  write.table(Pec_GRFs_Filtered_dataset_noOverlap, file = Save_FilterNoOverlap_Pec, sep =",", row.names=FALSE)
  
  ## Save the data taken at the peak Net GRF
  Save_FilterNoOverlap_PeakNet_Pec <- paste(Trial,"_Pec_Filtered_noOverlap_Peak",SaveDate, ".csv", sep="")
  write.table(Pec_GRFs_Filtered_dataset_noOverlap_Peak, file = Save_FilterNoOverlap_PeakNet_Pec, sep =",", row.names=FALSE)
}

## Pelvic data
if (!VideoInfo$Appendages == 'Pectoral') {
  ## Save the dataset that was filtered and had areas of overlap excluded
  Save_FilterAll_Pel <- paste(Trial,"_Pel_Filtered_All_",SaveDate, ".csv", sep="")
  write.table(Pel_GRFs_Filtered_dataset, file = Save_FilterAll_Pel, sep =",", row.names=FALSE)
  
  ## Save the dataset that was filtered and had areas of overlap excluded
  Save_FilterNoOverlap_Pel <- paste(Trial,"_Pel_Filtered_noOverlap_",SaveDate, ".csv", sep="")
  write.table(Pel_GRFs_Filtered_dataset_noOverlap, file = Save_FilterNoOverlap_Pel, sep =",", row.names=FALSE)
  
  ## Save the data taken at the peak Net GRF
  Save_FilterNoOverlap_PeakNet_Pel <- paste(Trial,"_Pel_Filtered_noOverlap_Peak",SaveDate, ".csv", sep="")
  write.table(Pel_GRFs_Filtered_dataset_noOverlap_Peak, file = Save_FilterNoOverlap_PeakNet_Pel, sep =",", row.names=FALSE)
}



