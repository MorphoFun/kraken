############### Comparing two padding methods for filtering data ##############
### Our previous R code would pad the beginning and end of stance with the mean value 
### calculated at the beginning and end of stance, respectively, and replicate that value
### so it extends the data 25% before and after stance.
### However, Chris Richards tends to mirror the data at the beginning and end of stance instead.
### Chris' method is a lot easier to code and seems to produce results that are less sensitive
### to noise at the edges of the data. 
### Otherwise, the traces are pretty much exactly the same. 



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


# Determining what sweep number the light is turned on, so we can sync with video frames
lightSwitch <- myData[which(myData$light_Volts<0),]
lightOnset <- lightSwitch[which(lightSwitch$sweep == min(lightSwitch$sweep)),]

# Synching video frames with force sweep numbers
lightStartFrame <- VideoInfo$Light.Start
startFrame <- VideoInfo$Pectoral.Start.Frame
endFrame <- VideoInfo$Pectoral.End.Frame
zeroStart <- VideoInfo$Pectoral.Start.Frame
forceHz <- 5000
videoHz <- 100

startSweep <- lightOnset$sweep-((lightStartFrame-startFrame)*(forceHz/videoHz))
endSweep <- lightOnset$sweep-((lightStartFrame-endFrame)*(forceHz/videoHz))
zeroSweep <- lightOnset$sweep-((lightStartFrame-zeroStart)*(forceHz/videoHz))

# Correcting for the offset from the baseline (zero)
offsetCalcStart <- zeroSweep-2000
offsetCalcEnd <- zeroSweep-1000
GRF_vertSumCalib_N_Offset <- mean(myData$GRF_vertSumCalib_N[offsetCalcStart:offsetCalcEnd])
GRF_MLSumCalib_N_Offset <- mean(myData$GRF_MLSumCalib_N[offsetCalcStart:offsetCalcEnd])
GRF_APSumCalib_N_Offset <- mean(myData$GRF_APSumCalib_N[offsetCalcStart:offsetCalcEnd])



#### CONVERTING LABVIEW OUTPUT TO GRFS ####

if (!VideoInfo$Appendages == 'Pelvic') Pec_GRFs <- voltToForce(myData, CalibInfo[,3:5], zeroStart = VideoInfo$Pectoral.Start.Frame, lightStartFrame = VideoInfo$Light.Start, startFrame = VideoInfo$Pectoral.Start.Frame, endFrame = VideoInfo$Pectoral.End.Frame, filename = Trial, BW = VideoInfo$Body.Weight.kg)
if (!VideoInfo$Appendages == 'Pectoral') Pel_GRFs <- voltToForce(myData, CalibInfo[,3:5], VideoInfo$Pectoral.Start.Frame, VideoInfo$Light.Start, VideoInfo$Pelvic.Start.Frame, VideoInfo$Pelvic.End.Frame, filename = Trial, BW = VideoInfo$Body.Weight.kg)


#### FILTERING THE DATA ####
if (!VideoInfo$Appendages == 'Pelvic') Pec_GRFs_Filtered <- butterFilteR(Pec_GRFs)
if (!VideoInfo$Appendages == 'Pectoral') Pel_GRFs_Filtered <- butterFilteR(Pel_GRFs)

#### LOADING THE DATA PROCESSED WITH THE OLD PADDING METHOD ####
old_Pec <- read.csv("af01f18_FilteredPecAll_120803.csv", header = TRUE)
old_Pel <- read.csv("af01f18_FilteredPelAll_120803.csv", header = TRUE)


#### PLOTTING THE DIFFERENCES ####
raw <- list("data" = myData[,c(13,17:19)], "filename" = "af01f18")
raw$data$Time_s <- c(0, cumsum(rep(1/forceHz, (nrow(raw$data)-1))))

## Pectoral appendage
quartz()
plot(raw$data$sweep[25000:33000], raw$data$GRF_vertSumCalib_N[25000:33000] - GRF_vertSumCalib_N_Offset, type = "l")
lines(Pec_GRFs_Filtered$GRF0Sum_filter_sweeps$sweep, Pec_GRFs_Filtered$GRF0Sum_filter_sweeps$filterVN, type = "l", col = "blue", lwd = 5)
lines(old_Pec$Sweep, old_Pec$FilterVN, type = "l", col = "red", lwd = 4)

quartz()
plot(raw$data$sweep[25000:33000], raw$data$GRF_MLSumCalib_N[25000:33000] - GRF_MLSumCalib_N_Offset, type = "l")
lines(Pec_GRFs_Filtered$GRF0Sum_filter_sweeps$sweep, Pec_GRFs_Filtered$GRF0Sum_filter_sweeps$filterMLN, type = "l", col = "blue", lwd = 5)
lines(old_Pec$Sweep, old_Pec$FilterMLN, type = "l", col = "red", lwd = 4)

quartz()
plot(raw$data$sweep[25000:33000], -(raw$data$GRF_APSumCalib_N[25000:33000] - GRF_APSumCalib_N_Offset), type = "l")
lines(Pec_GRFs_Filtered$GRF0Sum_filter_sweeps$sweep, Pec_GRFs_Filtered$GRF0Sum_filter_sweeps$filterAPN, type = "l", col = "blue", lwd = 5)
lines(old_Pec$Sweep, old_Pec$FilterHzN, type = "l", col = "red", lwd = 4)


#### COMPARING THE GRF ANGLES ####
if (!VideoInfo$Appendages == 'Pelvic') Pec_GRFs_Filtered_dataset <- GRFAngles(Pec_GRFs_Filtered$GRF0Sum_filter_interp)
if (!VideoInfo$Appendages == 'Pectoral') Pel_GRFs_Filtered_dataset <- GRFAngles(Pel_GRFs_Filtered$GRF0Sum_filter_interp)

old_Pec_angles <- read.csv("af01f18_Pec_NoIncrement_120804.csv", header = TRUE)
old_Pel_angles <- read.csv("af01f18_Pel_NoIncrement_120804.csv", header = TRUE)


#### PLOTTING DIFFERENCES IN THE GRF ANGLES ####
quartz()
plot(Pec_GRFs_Filtered_dataset[,1], Pec_GRFs_Filtered_dataset$MLAngle_Convert_deg, type = "l", lwd = 4)
lines(old_Pec_angles[,1], old_Pec_angles$MLAngle.Convert, col = "blue")

quartz()
plot(Pec_GRFs_Filtered_dataset[,1], Pec_GRFs_Filtered_dataset$APAngle_Convert_deg, type = "l", lwd = 4)
lines(old_Pec_angles[,1], old_Pec_angles$APAngle.Convert, col = "blue")
