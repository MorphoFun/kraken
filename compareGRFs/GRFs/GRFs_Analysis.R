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


CalibInfo <- CalibFile[CalibFile$Date %in% Date,]

#### CONVERTING LABVIEW OUTPUT TO GRFS ####

if (!VideoInfo$Appendages == 'Pelvic') Pec_GRFs <- voltToForce(myData, CalibInfo[,3:5], VideoInfo$Pectoral.Start.Frame, VideoInfo$Light.Start, VideoInfo$Pectoral.Start.Frame, VideoInfo$Pectoral.End.Frame, filename = Trial, BW = VideoInfo$Body.Weight.kg)
if (!VideoInfo$Appendages == 'Pectoral') Pel_GRFs <- voltToForce(myData, CalibInfo[,3:5], VideoInfo$Pectoral.Start.Frame, VideoInfo$Light.Start, VideoInfo$Pelvic.Start.Frame, VideoInfo$Pelvic.End.Frame, filename = Trial, BW = VideoInfo$Body.Weight.kg)


#### FILTERING THE DATA ####
if (!VideoInfo$Appendages == 'Pelvic') Pec_GRFs_Filtered <- butterFilteR(Pec_GRFs)
if (!VideoInfo$Appendages == 'Pectoral') Pel_GRFs_Filtered <- butterFilteR(Pel_GRFs)


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
