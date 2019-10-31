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
if (!VideoInfo$Appendages == 'Pelvic') Pec_GRFs_Filtered <- butterFilteR(Pec_GRFs)
if (!VideoInfo$Appendages == 'Pectoral') Pel_GRFs_Filtered <- butterFilteR(Pel_GRFs)


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

###################  PECTORAL APPENDAGE BEFORE PELVIC APPENDAGE FOR PECTORAL FILES #################

GRFoverlaps <- function(df, primary, secondary, filmRate, ...) {
  ### isolated footfalls
  if (is.na(secondary[,1])==TRUE) {
    output <- df
  }
  
  ### For files that have the primary appendage hitting the plate 1st
  if (is.na(secondary[,1])==FALSE & as.numeric(primary[,1])<as.numeric(secondary[,1])) {
  
    LimbCycleLength <- primary[,2]-primary[,1]
    OverlapStart <- floor(((secondary[,1]-primary[,1])/LimbCycleLength)*filmRate)  # want everything before overlap starts

    # Taking data points at every 5% of stance
    UsableRows <- df[[1]][1:(OverlapStart+1)]+1 # hindlimb overlap tends to occur towards the end
    df_noOverlap <- df[UsableRows,]

    # Making cells within the overlap as "NA" so the data.frame maintains the same dimensions with all 0 -> 101 pts at 5% increments included
    UnusableRows <- df[[1]][-c(1:(OverlapStart+1))]+1

    NA.Fillers <- data.frame(UnusableRows-1)
    NA.Fillers[ , paste("x", 1:(length(df_noOverlap)-1), sep = "")] <- NA
    names(NA.Fillers) <- names(df_noOverlap)
    
    # Adding filler NA's to actual data to keep/analyze
    df_WFill <- rbind(df_noOverlap, NA.Fillers)
    
    output <- df_WFill
  }
  
  ### For files that have the secondary appendage hitting the plate 1st
  if (is.na(secondary[,1])==FALSE & as.numeric(primary[,1])>as.numeric(secondary[,1]))
  {
    LimbCycleLength <- primary[,2]-primary[,1]
    OverlapEnd <- ceiling(((secondary[,2]-primary[,1])/LimbCycleLength)*filmRate) # Rounding up because want data after overlap is done

    # Taking data points at every 5% of stance
    UsableRows <- df[[1]][-c(0:OverlapEnd)]+1 # overlap tends to occur towards the beginning # adding 1 b/c data rows don't start until row 2
    
    df_noOverlap <- df[UsableRows,]
    
    # Making cells within the overlap as "NA" so the data.frame maintains the same dimensions with all 0 -> 101 pts at 5% increments included
    UnusableRows <- df[[1]][c(0:(OverlapEnd))]+1
    NA.Fillers <- data.frame(UnusableRows-1)
    NA.Fillers[ , paste("x", 1:(length(df_noOverlap)-1), sep = "")] <- NA
    names(NA.Fillers) <- names(df_noOverlap) # need to have the same variable names to rbind
    
    # Adding filler NA's to actual data to keep/analyze
    df_WFill <- rbind(NA.Fillers, df_noOverlap)
    
    output <- df_WFill
    
  }
  
  return(output)
}
    
    


#### SAVING THE DATA ####
SaveAllDataName <- paste(Trial,"_AllPrep_",SaveDate, ".csv", sep="")
write.table(myData_Forces$allData, file=SaveAllDataName, sep =",", row.names=FALSE)

# Saving the filter prep data
# Butterworth filtering will be conducted using another set of R code
# First, code the date that the file is being saved
if (!VideoInfo$Appendages == 'Pelvic') SaveFileName.Pec <- paste(Trial,"_FilterPrepPec_", SaveDate, ".csv", sep="")
if (!VideoInfo$Appendages == 'Pectoral') SaveFileName.Pel <- paste(Trial,"_FilterPrepPel_", SaveDate, ".csv", sep="")

if (!VideoInfo$Appendages == 'Pelvic') write.table(Pec_Forces$filterPrep, file=SaveFileName.Pec, sep =",", row.names=FALSE)
if (!VideoInfo$Appendages == 'Pectoral') write.table(Pel_Forces$filterPrep, file=SaveFileName.Pel, sep =",", row.names=FALSE)


#### 