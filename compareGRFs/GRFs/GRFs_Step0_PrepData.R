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

#### Loading libraries ####
if (!require("devtools")) {
  install.packages("devtools", dependencies = TRUE)
  library(devtools)
}
install_github("MorphoFun/kraken")
library(kraken)

# Read in values to calibrate each GRF component by; uses most recent Calibration Overview file
# when using list.files, need to use setwd as an assigned variable and also by itlself (see below) or it won't work
CalibFile <- data.frame(read.csv("./dataraw/FinLimbGRFs_Calibs.csv", header=TRUE))

# Read the most recent Video Info file
VideoFile <- data.frame(read.csv("./dataraw/FinLimbGRFs_VideoInfo.csv", header=TRUE))

setwd('./dataraw/Force_LabView_Output')

myFile <- file.choose()
myData <- read.table(myFile, header=FALSE)
myData <- myData[,c(1:12)] # Last 4 columns/channels were unused, so only subsetting what I need
names(myData) <- c("Light.Volts", "Vert1.Volts", "Vert2.Volts", "Vert3.Volts", "Vert4.Volts", "VertSum.Volts", "ML1.Volts", "ML2.Volts", "MLSum.Volts", "Hz1.Volts", "Hz2.Volts", "HzSum.Volts")
myData$Sweep <- 1:nrow(myData)

# Determining trial name
Trial <-  substring(myFile, nchar(myFile)-10, nchar(myFile)-4)

# Looking up video info
# Appendages listed as "Both" have both pectoral and pelvic appendage data
VideoInfo <- VideoFile[VideoFile$File.name %in% Trial,]
VideoInfo$Pectoral.Start.Frame <- as.numeric(VideoInfo$Pectoral.Start.Frame)
VideoInfo$Pectoral.End.Frame <- as.numeric(VideoInfo$Pectoral.End.Frame)
VideoInfo$Pelvic.Start.Frame <- as.numeric(VideoInfo$Pelvic.Start.Frame)
VideoInfo$Pelvic.End.Frame <- as.numeric(VideoInfo$Pelvic.End.Frame)
Date <- format(as.Date(VideoInfo$Date.Filmed, format = "%m/%d/%y"), format="%y%m%d")

## Converting the voltage readings from the force plate to measures of force

forces <- voltToForce

 
# Calibrating the raw force data and converting to newtons (if needed)
CalibInfo <- CalibFile[CalibFile$Date %in% Date,]
myData$VertSumCalib.N <- myData$VertSum.Volts*CalibInfo$Vert.Calib
myData$MLSumCalib.N <- myData$MLSum.Volts*CalibInfo$ML.Calib
myData$HzSumCalib.N <- myData$HzSum.Volts*CalibInfo$Hz.Calib

# Putting forces in terms of GRF (which is opposite in direction to the force produced by the limb onto the force plate)
myData$GRF.VertSumCalib.N <- myData$VertSumCalib.N # Already made negative based on the calibration calculations conducted earlier in the excel calibration files
myData$GRF.MLSumCalib.N <- -myData$MLSumCalib.N
myData$GRF.HzSumCalib.N <- myData$HzSumCalib.N # Rotating the force plate results in the signage changes, so don't need to multiply by -1

# Determining what sweep number the light is turned on, so I can sync with video frames
LightSwitch <- myData[which(myData$Light.Volts<0),]
LightOnset <- LightSwitch[which(LightSwitch$Sweep == min(LightSwitch$Sweep)),]

# Synching video frames with force sweep numbers
# *50 is due to the conversion from sweeps to frames (5000 Hz for force/100 Hz for video)
if (!is.na(VideoInfo$Pectoral.Start.Frame)) PectoralStartSweep <- LightOnset$Sweep-((VideoInfo$Light.Start-VideoInfo$Pectoral.Start.Frame)*50)
if (!is.na(VideoInfo$Pectoral.End.Frame)) PectoralEndSweep <- LightOnset$Sweep-((VideoInfo$Light.Start-VideoInfo$Pectoral.End.Frame)*50)
if (!is.na(VideoInfo$Pelvic.Start.Frame)) PelvicStartSweep <- LightOnset$Sweep-((VideoInfo$Light.Start-VideoInfo$Pelvic.Start.Frame)*50)
if (!is.na(VideoInfo$Pelvic.End.Frame)) PelvicEndSweep <- LightOnset$Sweep-((VideoInfo$Light.Start-VideoInfo$Pelvic.End.Frame)*50)

# Calculating the difference between the force trace zero and what should really be zero on the force trace (i.e., calculating horizontal offset)
# First, determine whether the pectoral or pelvic appendages appear first on the force trace
if (VideoInfo$Appendages == 'Pectoral') First <- data.frame(PectoralStartSweep)
if (VideoInfo$Appendages == 'Pelvic') First <- data.frame(PelvicStartSweep)
if (VideoInfo$Appendages == 'Both') First <- data.frame(PectoralStartSweep,PelvicStartSweep)
FirstTrace <- First[which(First==min(First))]
OffsetCalcStart <- (FirstTrace[1,1])-2000
OffsetCalcEnd <- (FirstTrace[1,1])-1000
GRF.VertSumCalib.N.Offset <- mean(myData$GRF.VertSumCalib.N[OffsetCalcStart:OffsetCalcEnd])
GRF.MLSumCalib.N.Offset <- mean(myData$GRF.MLSumCalib.N[OffsetCalcStart:OffsetCalcEnd])
GRF.HzSumCalib.N.Offset <- mean(myData$GRF.HzSumCalib.N[OffsetCalcStart:OffsetCalcEnd])

# Zeroing the force trace data using the horizontal offset value
myData$GRF.VertSumCalib.N.Zero <- myData$GRF.VertSumCalib.N-GRF.VertSumCalib.N.Offset
myData$GRF.MLSumCalib.N.Zero <- myData$GRF.MLSumCalib.N-GRF.MLSumCalib.N.Offset
myData$GRF.HzSumCalib.N.Zero <- myData$GRF.HzSumCalib.N-GRF.HzSumCalib.N.Offset

# Plotting force traces
if (VideoInfo$Appendages == 'Pectoral') Last <- data.frame(PectoralEndSweep)
if (VideoInfo$Appendages == 'Pelvic') Last <- data.frame(PelvicEndSweep)
if (VideoInfo$Appendages == 'Both') Last <- data.frame(PectoralEndSweep, PelvicEndSweep)
LastTrace <- Last[which(Last==max(Last))]
PlotStart <- FirstTrace[1,1]-1000
PlotEnd <- LastTrace[1,1]+1000


if (!is.na(VideoInfo$Pectoral.Start.Frame) & !is.na(VideoInfo$Pelvic.Start.Frame)) {
ImpPointsX <- data.frame(PectoralStartSweep, PectoralEndSweep,PelvicStartSweep, PelvicEndSweep)
names(ImpPointsX) <- c('Pectoral Start', 'Pectoral End','Pelvic Start', 'Pelvic End') }

if (!is.na(VideoInfo$Pectoral.Start.Frame) & is.na(VideoInfo$Pelvic.Start.Frame)) {
ImpPointsX <- data.frame(PectoralStartSweep, PectoralEndSweep)
names(ImpPointsX) <- c('Pectoral Start', 'Pectoral End') }

if (is.na(VideoInfo$Pectoral.Start.Frame) & !is.na(VideoInfo$Pelvic.Start.Frame)) {
ImpPointsX <- data.frame(PelvicStartSweep, PelvicEndSweep)
names(ImpPointsX) <- c('Pelvic Start', 'Pelvic End') }
  
attach(myData)
if (!is.na(VideoInfo$Pectoral.Start.Frame) & is.na(VideoInfo$Pelvic.Start.Frame)) {
ImpPoints.GRFVert <- data.frame(GRF.VertSumCalib.N.Zero[PectoralStartSweep], GRF.VertSumCalib.N.Zero[PectoralEndSweep])
ImpPoints.GRFML <- data.frame(GRF.MLSumCalib.N.Zero[PectoralStartSweep], GRF.MLSumCalib.N.Zero[PectoralEndSweep])
ImpPoints.GRFHz <- data.frame(GRF.HzSumCalib.N.Zero[PectoralStartSweep], GRF.HzSumCalib.N.Zero[PectoralEndSweep])
}

if (is.na(VideoInfo$Pectoral.Start.Frame) & !is.na(VideoInfo$Pelvic.Start.Frame)) {
ImpPoints.GRFVert <- data.frame(GRF.VertSumCalib.N.Zero[PelvicStartSweep], GRF.VertSumCalib.N.Zero[PelvicEndSweep])
ImpPoints.GRFML <- data.frame(GRF.MLSumCalib.N.Zero[PelvicStartSweep], GRF.MLSumCalib.N.Zero[PelvicEndSweep])
ImpPoints.GRFHz <- data.frame(GRF.HzSumCalib.N.Zero[PelvicStartSweep], GRF.HzSumCalib.N.Zero[PelvicEndSweep])
}

if (!is.na(VideoInfo$Pectoral.Start.Frame) & !is.na(VideoInfo$Pelvic.Start.Frame)) {
ImpPoints.GRFVert <- data.frame(GRF.VertSumCalib.N.Zero[PectoralStartSweep], GRF.VertSumCalib.N.Zero[PectoralEndSweep], GRF.VertSumCalib.N.Zero[PelvicStartSweep], GRF.VertSumCalib.N.Zero[PelvicEndSweep])
ImpPoints.GRFML <- data.frame(GRF.MLSumCalib.N.Zero[PectoralStartSweep], GRF.MLSumCalib.N.Zero[PectoralEndSweep], GRF.MLSumCalib.N.Zero[PelvicStartSweep], GRF.MLSumCalib.N.Zero[PelvicEndSweep])
ImpPoints.GRFHz <- data.frame(GRF.HzSumCalib.N.Zero[PectoralStartSweep], GRF.HzSumCalib.N.Zero[PectoralEndSweep], GRF.HzSumCalib.N.Zero[PelvicStartSweep], GRF.HzSumCalib.N.Zero[PelvicEndSweep])
}
detach(myData)

quartz(width=10)
par(mfrow=c(1,3), oma = c(0, 0, 2, 0))  # oma = outer margin with 2 lines above the top of the graphs
# Vertical component of GRF graph
plot(myData$Sweep[PlotStart:PlotEnd], myData$GRF.VertSumCalib.N.Zero[PlotStart:PlotEnd], xlab='Sweep', ylab='GRF - Vertical (N)', main='Zeroed GRF (Vertical) Force', type="l", col="blue")
points(ImpPointsX[1,], ImpPoints.GRFVert[1,], type='p', pch='O', col='cyan')
text(ImpPointsX[1,], ImpPoints.GRFVert[1,], labels=names(ImpPointsX), pos=3, font=2) # pos: 1 = below, 2 = left, 3 = above, 4 = right

# Mediolateral component of GRF graph
plot(myData$Sweep[PlotStart:PlotEnd], myData$GRF.MLSumCalib.N.Zero[PlotStart:PlotEnd], xlab='Sweep', ylab='GRF - Mediolateral (N)', main='Zeroed GRF (Mediolateral) Force', type="l", col="red")
points(ImpPointsX[1,], ImpPoints.GRFML[1,], type='p', pch='O', col='cyan')
text(ImpPointsX[1,], ImpPoints.GRFML[1,], labels=names(ImpPointsX), pos=3, font=2) # pos: 1 = below, 2 = left, 3 = above, 4 = right

# Horizontal (Anteroposterior) component of GRF graph
plot(myData$Sweep[PlotStart:PlotEnd], myData$GRF.HzSumCalib.N.Zero[PlotStart:PlotEnd], xlab='Sweep', ylab='GRF - Horizontal (N)', main='Zeroed GRF (Horizontal) Force', type="l", col="forestgreen")
points(ImpPointsX[1,], ImpPoints.GRFHz[1,], type='p', pch='O', col='cyan')
text(ImpPointsX[1,], ImpPoints.GRFHz[1,], labels=names(ImpPointsX), pos=3, font=2) # pos: 1 = below, 2 = left, 3 = above, 4 = right

mtext(Trial, line=0.5, outer=TRUE)  # writes an overall title over the graphs


# For saving as pdf
PdfSave <- paste(substring(myFile, 1, 74), "R Analysis/Step 1 Calibrate and Organize Data/Save All Data/", Trial, "_Filter_", SaveDate, ".pdf", sep="")
pdf(PdfSave, width=11)

par(mfrow=c(1,3), oma = c(0, 0, 2, 0))  # oma = outer margin with 2 lines above the top of the graphs
# Vertical component of GRF graph
plot(myData$Sweep[PlotStart:PlotEnd], myData$GRF.VertSumCalib.N.Zero[PlotStart:PlotEnd], xlab='Sweep', ylab='GRF - Vertical (N)', main='Zeroed GRF (Vertical) Force', type="l", col="blue")
points(ImpPointsX[1,], ImpPoints.GRFVert[1,], type='p', pch='O', col='cyan')
text(ImpPointsX[1,], ImpPoints.GRFVert[1,], labels=names(ImpPointsX), pos=3, font=2) # pos: 1 = below, 2 = left, 3 = above, 4 = right

# Mediolateral component of GRF graph
plot(myData$Sweep[PlotStart:PlotEnd], myData$GRF.MLSumCalib.N.Zero[PlotStart:PlotEnd], xlab='Sweep', ylab='GRF - Mediolateral (N)', main='Zeroed GRF (Mediolateral) Force', type="l", col="red")
points(ImpPointsX[1,], ImpPoints.GRFML[1,], type='p', pch='O', col='cyan')
text(ImpPointsX[1,], ImpPoints.GRFML[1,], labels=names(ImpPointsX), pos=3, font=2) # pos: 1 = below, 2 = left, 3 = above, 4 = right

# Horizontal (Anteroposterior) component of GRF graph
plot(myData$Sweep[PlotStart:PlotEnd], myData$GRF.HzSumCalib.N.Zero[PlotStart:PlotEnd], xlab='Sweep', ylab='GRF - Horizontal (N)', main='Zeroed GRF (Horizontal) Force', type="l", col="forestgreen")
points(ImpPointsX[1,], ImpPoints.GRFHz[1,], type='p', pch='O', col='cyan')
text(ImpPointsX[1,], ImpPoints.GRFHz[1,], labels=names(ImpPointsX), pos=3, font=2) # pos: 1 = below, 2 = left, 3 = above, 4 = right

mtext(Trial, line=0.5, outer=TRUE)  # writes an overall title over the graphs

dev.off()  # for use with pdf()



# Preparing the data to be filtered
# Pectoral appendage stance
if (!VideoInfo$Appendages == 'Pelvic') {

# Calculating the range of sweeps for your stance
StanceSweeps.Pec <- myData$Sweep[PectoralStartSweep:PectoralEndSweep]
# Leaving extra sweeps in front and back of trace to reduce likelihood of edge effects (Extra sweeps = 25% of stance duration)
Extra.Pec <- round(.25*length(StanceSweeps.Pec))
FilterPrep.Sweeps.Pec <- myData$Sweep[c(PectoralStartSweep-Extra.Pec):c(PectoralEndSweep+Extra.Pec)]
FilterGRFVert.Pec <- myData$GRF.VertSumCalib.N.Zero[c(PectoralStartSweep-Extra.Pec):c(PectoralEndSweep+Extra.Pec)]
FilterGRFML.Pec <- myData$GRF.MLSumCalib.N.Zero[c(PectoralStartSweep-Extra.Pec):c(PectoralEndSweep+Extra.Pec)]
FilterGRFHz.Pec <- myData$GRF.HzSumCalib.N.Zero[c(PectoralStartSweep-Extra.Pec):c(PectoralEndSweep+Extra.Pec)]
FilterPrep.Pec <- data.frame(FilterPrep.Sweeps.Pec,FilterGRFVert.Pec, FilterGRFML.Pec, FilterGRFHz.Pec)
names(FilterPrep.Pec) <- c('Sweep', 'GRF0SumVN', 'GRF0SumMLN', 'GRF0SumHzN')
}

if (!VideoInfo$Appendages == 'Pectoral') {
StanceSweeps.Pel <- myData$Sweep[PelvicStartSweep:PelvicEndSweep]
Extra.Pel <- round(.25*length(StanceSweeps.Pel))
FilterPrep.Sweeps.Pel <- myData$Sweep[c(PelvicStartSweep-Extra.Pel):c(PelvicEndSweep+Extra.Pel)]
FilterGRFVert.Pel <- myData$GRF.VertSumCalib.N.Zero[c(PelvicStartSweep-Extra.Pel):c(PelvicEndSweep+Extra.Pel)]
FilterGRFML.Pel <- myData$GRF.MLSumCalib.N.Zero[c(PelvicStartSweep-Extra.Pel):c(PelvicEndSweep+Extra.Pel)]
FilterGRFHz.Pel <- myData$GRF.HzSumCalib.N.Zero[c(PelvicStartSweep-Extra.Pel):c(PelvicEndSweep+Extra.Pel)]
FilterPrep.Pel <- data.frame(FilterPrep.Sweeps.Pel,FilterGRFVert.Pel, FilterGRFML.Pel, FilterGRFHz.Pel)
names(FilterPrep.Pel) <- c('Sweep', 'GRF0SumVN', 'GRF0SumMLN', 'GRF0SumHzN')
}

# Saving all of the data
setwd('/Users/SandyMKawano/Desktop/Research/Forelimbs on Force Plates/Force Data/R Analysis/Step 1 Calibrate and Organize Data/Save All Data')
SaveAllDataName <- paste(Trial,"_AllPrep_",SaveDate, ".csv", sep="")
write.table(myData, file=SaveAllDataName, sep =",", row.names=FALSE)

# Saving the filter prep data
# Butterworth filtering will be conducted using another set of R code
# First, code the date that the file is being saved
if (!VideoInfo$Appendages == 'Pelvic') SaveFileName.Pec <- paste(Trial,"_FilterPrepPec_", SaveDate, ".csv", sep="")
if (!VideoInfo$Appendages == 'Pectoral') SaveFileName.Pel <- paste(Trial,"_FilterPrepPel_", SaveDate, ".csv", sep="")

# Changing directory
setwd('/Users/SandyMKawano/Desktop/Research/Forelimbs on Force Plates/Force Data/R Analysis/Step 1 Calibrate and Organize Data/Prep Data for Filtering')
if (!VideoInfo$Appendages == 'Pelvic') write.table(FilterPrep.Pec, file=SaveFileName.Pec, sep =",", row.names=FALSE)
if (!VideoInfo$Appendages == 'Pectoral') write.table(FilterPrep.Pel, file=SaveFileName.Pel, sep =",", row.names=FALSE)


# Changing directory for saving the figure of the graphs
setwd('/Users/SandyMKawano/Desktop/Research/Forelimbs on Force Plates/Force Data/R Analysis/Step 1 Calibrate and Organize Data/Save All Data')

