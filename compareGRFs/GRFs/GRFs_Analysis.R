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
    
    
    if (!is.na(VideoInfo$Body.1.Start))
    {
      if (VideoInfo$Body.1.Start>=VideoInfo$Pectoral.Start.Frame & VideoInfo$Body.1.Start<=VideoInfo$Pectoral.End.Frame)
      {
        Body1Start <- round(((VideoInfo$Body.1.Start-VideoInfo$Pectoral.Start.Frame)/LimbCycleLength)*VideoInfo$Filming.Rate.Hz)
        ImpPoints.Body1Start <- data.frame(Body1Start, df$InterpV_BW[Body1Start+1], df$InterpML_BW[Body1Start+1], df$InterpHz_BW[Body1Start+1], df$NetGRF_BW[Body1Start+1], row.names='Body 1 Start')
        names(ImpPoints.Body1Start) <- c('Percent Stance', 'Vertical (BW)', 'Mediolateral (BW)', 'Horizontal (BW)', 'Net GRF (BW)')
      }
    }
    if (!is.na(VideoInfo$Body.1.End))
    {
      if (VideoInfo$Body.1.End>=VideoInfo$Pectoral.Start.Frame & VideoInfo$Body.1.End<=VideoInfo$Pectoral.End.Frame)
      {
        Body1End <- round(((VideoInfo$Body.1.End-VideoInfo$Pectoral.Start.Frame)/LimbCycleLength)*VideoInfo$Filming.Rate.Hz)
        ImpPoints.Body1End <- data.frame(Body1End, df$InterpV_BW[Body1End+1], df$InterpML_BW[Body1End+1], df$InterpHz_BW[Body1End+1], df$NetGRF_BW[Body1End+1], row.names='Body 1 End')
        names(ImpPoints.Body1End) <- c('Percent Stance', 'Vertical (BW)', 'Mediolateral (BW)', 'Horizontal (BW)', 'Net GRF (BW)')
      }
    }
    
    if (!is.na(VideoInfo$Body.2.Start))
    {
      if (VideoInfo$Body.2.Start>=VideoInfo$Pectoral.Start.Frame & VideoInfo$Body.2.Start<=VideoInfo$Pectoral.End.Frame)
      {
        Body2Start <- round(((VideoInfo$Body.2.Start-VideoInfo$Pectoral.Start.Frame)/LimbCycleLength)*VideoInfo$Filming.Rate.Hz)
        ImpPoints.Body2Start <- data.frame(Body2Start, df$InterpV_BW[Body2Start+1], df$InterpML_BW[Body2Start+1], df$InterpHz_BW[Body2Start+1], df$NetGRF_BW[Body2Start+1], row.names='Body 2 Start')
        names(ImpPoints.Body2Start) <- c('Percent Stance', 'Vertical (BW)', 'Mediolateral (BW)', 'Horizontal (BW)', 'Net GRF (BW)')
      }
    }
    if (!is.na(VideoInfo$Body.2.End))
    {
      if (VideoInfo$Body.2.End>=VideoInfo$Pectoral.Start.Frame & VideoInfo$Body.2.End<=VideoInfo$Pectoral.End.Frame)
      {
        Body2End <- round(((VideoInfo$Body.2.End-VideoInfo$Pectoral.Start.Frame)/LimbCycleLength)*VideoInfo$Filming.Rate.Hz)
        ImpPoints.Body2End <- data.frame(Body2End, df$InterpV_BW[Body2End+1], df$InterpML_BW[Body2End+1], df$InterpHz_BW[Body2End+1], df$NetGRF_BW[Body2End+1], row.names='Body 2 End')
        names(ImpPoints.Body2End) <- c('Percent Stance', 'Vertical (BW)', 'Mediolateral (BW)', 'Horizontal (BW)', 'Net GRF (BW)')
      }
    }
    
    if (exists("ImpPoints.Body1Start") & exists("ImpPoints.Body1End") & exists("ImpPoints.Body2Start") & exists("ImpPoints.Body2End"))
    {
      ImpPoints <- rbind(ImpPoints.Body1Start, ImpPoints.Body1End, ImpPoints.Body2Start, ImpPoints.Body2End)
    }
    if (exists("ImpPoints.Body1Start") & exists("ImpPoints.Body1End") & exists("ImpPoints.Body2Start") & exists("ImpPoints.Body2End")=='FALSE')
    {
      ImpPoints <- rbind(ImpPoints.Body1Start, ImpPoints.Body1End, ImpPoints.Body2Start)
    }
    if (exists("ImpPoints.Body1Start")=='FALSE' & exists("ImpPoints.Body1End") & exists("ImpPoints.Body2Start") & exists("ImpPoints.Body2End"))
    {
      ImpPoints <- rbind(ImpPoints.Body1End, ImpPoints.Body2Start, ImpPoints.Body2End)
    }
    if (exists("ImpPoints.Body1Start")=='FALSE' & exists("ImpPoints.Body1End") & exists("ImpPoints.Body2Start") & exists("ImpPoints.Body2End"))
    {
      ImpPoints <- rbind(ImpPoints.Body1End, ImpPoints.Body2Start)
    }
    if (exists("ImpPoints.Body1Start") & exists("ImpPoints.Body1End") & exists("ImpPoints.Body2Start")=='FALSE' & exists("ImpPoints.Body2End")=='FALSE')
    {
      ImpPoints <- rbind(ImpPoints.Body1Start, ImpPoints.Body1End)
    }
    if (exists("ImpPoints.Body1Start")=='FALSE' & exists("ImpPoints.Body1End")=='FALSE' & exists("ImpPoints.Body2Start") & exists("ImpPoints.Body2End"))
    {
      ImpPoints <- rbind(ImpPoints.Body2Start, ImpPoints.Body2End)
    }
    if (exists("ImpPoints.Body1Start")=='FALSE' & exists("ImpPoints.Body1End")=='FALSE' & exists("ImpPoints.Body2Start")=='FALSE' & exists("ImpPoints.Body2End"))
    {
      ImpPoints <- rbind(ImpPoints.Body2End)
    }
    if (exists("ImpPoints.Body1Start")=='FALSE' & exists("ImpPoints.Body1End")=='FALSE' & exists("ImpPoints.Body2Start") & exists("ImpPoints.Body2End")=='FALSE')
    {
      ImpPoints <- rbind(ImpPoints.Body2Start)
    }
    if (exists("ImpPoints.Body1Start") & exists("ImpPoints.Body1End")=='FALSE' & exists("ImpPoints.Body2Start")=='FALSE' & exists("ImpPoints.Body2End")=='FALSE')
    {
      ImpPoints <- rbind(ImpPoints.Body1Start)
    }
    if (exists("ImpPoints.Body1Start")=='FALSE' & exists("ImpPoints.Body1End") & exists("ImpPoints.Body2Start")=='FALSE' & exists("ImpPoints.Body2End")=='FALSE')
    {
      ImpPoints <- rbind(ImpPoints.Body1End)
    }
    ifelse ((exists("ImpPoints")=='TRUE'), ImpPoints <-rbind(ImpPoints.Overlap,ImpPoints), ImpPoints <-rbind(ImpPoints.Overlap))
    
       
  }
    

  
  ## pectoral before pelvic
  if (!VideoInfo$Appendages == 'Pelvic' & is.na(VideoInfo$Pelvic.Start.Frame)==FALSE & as.numeric(VideoInfo$Pectoral.Start.Frame)<as.numeric(VideoInfo$Pelvic.Start.Frame)) status <- "pec_pecToPel"
  
  ## pelbic before pectoral
  if (!VideoInfo$Appendages == 'Pelvic' & is.na(VideoInfo$Pelvic.Start.Frame)==FALSE & as.numeric(VideoInfo$Pectoral.Start.Frame)>as.numeric(VideoInfo$Pelvic.Start.Frame)) status <- "pec_pelToPec"
  
  ### Pelvic trials  
  ## only pelvic
  if (!VideoInfo$Appendages == 'Pectoral' & is.na(VideoInfo$Pectoral.Start.Frame)==TRUE) status <- "pel_pelOnly"

  ## pectoral before pelvic
  if (!VideoInfo$Appendages == 'Pectoral' & !is.na(VideoInfo$Pectoral.Start.Frame)==TRUE & as.numeric(VideoInfo$Pectoral.Start.Frame)<as.numeric(VideoInfo$Pelvic.Start.Frame)) status <- "pel_pecToPel"
    
  
}

if (!VideoInfo$Appendages == 'Pelvic' & is.na(VideoInfo$Pelvic.Start.Frame)==FALSE & as.numeric(VideoInfo$Pectoral.Start.Frame)<as.numeric(VideoInfo$Pelvic.Start.Frame)) # For files that have the pectoral appendage hitting the plate 1st
{
  
  LimbCycleLength <- VideoInfo$Pectoral.End.Frame-VideoInfo$Pectoral.Start.Frame
  OverlapStart <- floor(((VideoInfo$Pelvic.Start.Frame-VideoInfo$Pectoral.Start.Frame)/LimbCycleLength)*VideoInfo$Filming.Rate.Hz)  # want everything before overlap starts
  ImpPoints.Overlap <- data.frame(OverlapStart, FilterInterp$InterpV.BW[OverlapStart+1], FilterInterp$InterpML.BW[OverlapStart+1], FilterInterp$InterpHz.BW[OverlapStart+1], FilterInterp$NetGRF.BW[OverlapStart+1], row.names='Overlap Start')
  # need to include +1 within square brackets b/c obs rows are 1 value greater than the % stance (i.e., 0% stance = row 1)
  names(ImpPoints.Overlap) <- c('Percent Stance', 'Vertical (BW)', 'Mediolateral (BW)', 'Horizontal (BW)', 'Net GRF (BW)')
  
  if (!is.na(VideoInfo$Body.1.Start))
  {
    if (VideoInfo$Body.1.Start>=VideoInfo$Pectoral.Start.Frame & VideoInfo$Body.1.Start<=VideoInfo$Pectoral.End.Frame)
    {
      Body1Start <- round(((VideoInfo$Body.1.Start-VideoInfo$Pectoral.Start.Frame)/LimbCycleLength)*VideoInfo$Filming.Rate.Hz)
      ImpPoints.Body1Start <- data.frame(Body1Start, FilterInterp$InterpV.BW[Body1Start+1], FilterInterp$InterpML.BW[Body1Start+1], FilterInterp$InterpHz.BW[Body1Start+1], FilterInterp$NetGRF.BW[Body1Start+1], row.names='Body 1 Start')
      names(ImpPoints.Body1Start) <- c('Percent Stance', 'Vertical (BW)', 'Mediolateral (BW)', 'Horizontal (BW)', 'Net GRF (BW)')
    }
  }
  if (!is.na(VideoInfo$Body.1.End))
  {
    if (VideoInfo$Body.1.End>=VideoInfo$Pectoral.Start.Frame & VideoInfo$Body.1.End<=VideoInfo$Pectoral.End.Frame)
    {
      Body1End <- round(((VideoInfo$Body.1.End-VideoInfo$Pectoral.Start.Frame)/LimbCycleLength)*VideoInfo$Filming.Rate.Hz)
      ImpPoints.Body1End <- data.frame(Body1End, FilterInterp$InterpV.BW[Body1End+1], FilterInterp$InterpML.BW[Body1End+1], FilterInterp$InterpHz.BW[Body1End+1], FilterInterp$NetGRF.BW[Body1End+1], row.names='Body 1 End')
      names(ImpPoints.Body1End) <- c('Percent Stance', 'Vertical (BW)', 'Mediolateral (BW)', 'Horizontal (BW)', 'Net GRF (BW)')
    }
  }
  
  if (!is.na(VideoInfo$Body.2.Start))
  {
    if (VideoInfo$Body.2.Start>=VideoInfo$Pectoral.Start.Frame & VideoInfo$Body.2.Start<=VideoInfo$Pectoral.End.Frame)
    {
      Body2Start <- round(((VideoInfo$Body.2.Start-VideoInfo$Pectoral.Start.Frame)/LimbCycleLength)*VideoInfo$Filming.Rate.Hz)
      ImpPoints.Body2Start <- data.frame(Body2Start, FilterInterp$InterpV.BW[Body2Start+1], FilterInterp$InterpML.BW[Body2Start+1], FilterInterp$InterpHz.BW[Body2Start+1], FilterInterp$NetGRF.BW[Body2Start+1], row.names='Body 2 Start')
      names(ImpPoints.Body2Start) <- c('Percent Stance', 'Vertical (BW)', 'Mediolateral (BW)', 'Horizontal (BW)', 'Net GRF (BW)')
    }
  }
  if (!is.na(VideoInfo$Body.2.End))
  {
    if (VideoInfo$Body.2.End>=VideoInfo$Pectoral.Start.Frame & VideoInfo$Body.2.End<=VideoInfo$Pectoral.End.Frame)
    {
      Body2End <- round(((VideoInfo$Body.2.End-VideoInfo$Pectoral.Start.Frame)/LimbCycleLength)*VideoInfo$Filming.Rate.Hz)
      ImpPoints.Body2End <- data.frame(Body2End, FilterInterp$InterpV.BW[Body2End+1], FilterInterp$InterpML.BW[Body2End+1], FilterInterp$InterpHz.BW[Body2End+1], FilterInterp$NetGRF.BW[Body2End+1], row.names='Body 2 End')
      names(ImpPoints.Body2End) <- c('Percent Stance', 'Vertical (BW)', 'Mediolateral (BW)', 'Horizontal (BW)', 'Net GRF (BW)')
    }
  }
  
  if (exists("ImpPoints.Body1Start") & exists("ImpPoints.Body1End") & exists("ImpPoints.Body2Start") & exists("ImpPoints.Body2End"))
  {
    ImpPoints <- rbind(ImpPoints.Body1Start, ImpPoints.Body1End, ImpPoints.Body2Start, ImpPoints.Body2End)
  }
  if (exists("ImpPoints.Body1Start") & exists("ImpPoints.Body1End") & exists("ImpPoints.Body2Start") & exists("ImpPoints.Body2End")=='FALSE')
  {
    ImpPoints <- rbind(ImpPoints.Body1Start, ImpPoints.Body1End, ImpPoints.Body2Start)
  }
  if (exists("ImpPoints.Body1Start")=='FALSE' & exists("ImpPoints.Body1End") & exists("ImpPoints.Body2Start") & exists("ImpPoints.Body2End"))
  {
    ImpPoints <- rbind(ImpPoints.Body1End, ImpPoints.Body2Start, ImpPoints.Body2End)
  }
  if (exists("ImpPoints.Body1Start")=='FALSE' & exists("ImpPoints.Body1End") & exists("ImpPoints.Body2Start") & exists("ImpPoints.Body2End"))
  {
    ImpPoints <- rbind(ImpPoints.Body1End, ImpPoints.Body2Start)
  }
  if (exists("ImpPoints.Body1Start") & exists("ImpPoints.Body1End") & exists("ImpPoints.Body2Start")=='FALSE' & exists("ImpPoints.Body2End")=='FALSE')
  {
    ImpPoints <- rbind(ImpPoints.Body1Start, ImpPoints.Body1End)
  }
  if (exists("ImpPoints.Body1Start")=='FALSE' & exists("ImpPoints.Body1End")=='FALSE' & exists("ImpPoints.Body2Start") & exists("ImpPoints.Body2End"))
  {
    ImpPoints <- rbind(ImpPoints.Body2Start, ImpPoints.Body2End)
  }
  if (exists("ImpPoints.Body1Start")=='FALSE' & exists("ImpPoints.Body1End")=='FALSE' & exists("ImpPoints.Body2Start")=='FALSE' & exists("ImpPoints.Body2End"))
  {
    ImpPoints <- rbind(ImpPoints.Body2End)
  }
  if (exists("ImpPoints.Body1Start")=='FALSE' & exists("ImpPoints.Body1End")=='FALSE' & exists("ImpPoints.Body2Start") & exists("ImpPoints.Body2End")=='FALSE')
  {
    ImpPoints <- rbind(ImpPoints.Body2Start)
  }
  if (exists("ImpPoints.Body1Start") & exists("ImpPoints.Body1End")=='FALSE' & exists("ImpPoints.Body2Start")=='FALSE' & exists("ImpPoints.Body2End")=='FALSE')
  {
    ImpPoints <- rbind(ImpPoints.Body1Start)
  }
  if (exists("ImpPoints.Body1Start")=='FALSE' & exists("ImpPoints.Body1End") & exists("ImpPoints.Body2Start")=='FALSE' & exists("ImpPoints.Body2End")=='FALSE')
  {
    ImpPoints <- rbind(ImpPoints.Body1End)
  }
  ifelse ((exists("ImpPoints")=='TRUE'), ImpPoints <-rbind(ImpPoints.Overlap,ImpPoints), ImpPoints <-rbind(ImpPoints.Overlap))
  
  # Taking data points at every 5% of stance
  IncrementRows <- seq(1,101,5) # row 1 = 0% stance, so essentially need to do row#-1 for each percentage of stance
  UsableRows <- FilterInterp$PercentStance[1:(OverlapStart+1)]+1 # hindlimb overlap tends to occur towards the end
  StanceWithoutOverlap.Increment <- IncrementRows[IncrementRows %in% UsableRows] 
  
  FilterInterp.NoIncrement <- FilterInterp[UsableRows,c(1,10,12:16)]
  FilterInterp.Increment <- FilterInterp[StanceWithoutOverlap.Increment,c(1,10,12:16)]
  
  # Making cells within the overlap as "NA" so the data.frame maintains the same dimensions with all 0 -> 101 pts at 5% increments included
  UnusableRows <- FilterInterp$PercentStance[-c(1:(OverlapStart+1))]+1
  StanceWithOverlap.Increment <- IncrementRows[IncrementRows %in% UnusableRows]
  NA.Increment <- rep("NA", length(StanceWithOverlap.Increment))
  StanceWithOverlap.IncrementPercent <- StanceWithOverlap.Increment-1
  NA.Fillers <- data.frame(cbind(StanceWithOverlap.IncrementPercent, NA.Increment, NA.Increment, NA.Increment, NA.Increment, NA.Increment, NA.Increment))
  names(NA.Fillers) <- names(FilterInterp.Increment) # need to have the same variable names to rbind
  
  NA.NoIncrement <- rep("NA", length(seq(1,101,1)))  
  StanceWithOverlap.Percent <- UnusableRows-1
  NA.Fillers.NoIncrement <- data.frame(cbind(StanceWithOverlap.Percent, NA.Increment, NA.Increment, NA.Increment, NA.Increment, NA.Increment, NA.Increment))
  names(NA.Fillers.NoIncrement) <- names(FilterInterp)
  
  # Adding filler NA's to actual data to keep/analyze
  FilterInterp.IncrementWFill <- rbind(FilterInterp.Increment, NA.Fillers)
  FilterInterp.NoIncrementWFill <- rbind(FilterINterp.NoIncrement, NA.Fillers.NoIncrement)
  
  # Evaluating at peak/max net GRF
  # Also making sure that I'm only including data points that occurred before overlap with other structures started
  PeakNetGRF <- FilterInterp[which(FilterInterp$NetGRF.BW==max(FilterInterp$NetGRF.BW[1:OverlapStart+1])),]
  
  # Creating file names for saving as .csv
  SaveFilteredNoIncrement <- paste(Trial,"_Pec_NoIncrement_", SaveDate, ".csv", sep="")
  SaveFilteredIncrement <- paste(Trial,"_Pec_Increment_", SaveDate, ".csv", sep="")
  SaveFilteredPeakNet<- paste(Trial,"_Pec_PeakNetGRF_", SaveDate, ".csv", sep="")
  SaveFiltered101 <- paste(Trial, "_Pec_101_", SaveDate, ".txt", sep="")
  
  # Set the directory to save the data files
  setwd('/Users/SandyMKawano/Desktop/Research/Forelimbs on Force Plates/Force Data/R Analysis/Step 3 Post Filter Processing/Filtered at Increments')
  write.table(FilterInterp.IncrementWFill, file=SaveFilteredIncrement, sep =",", row.names=FALSE)
  
  # setwd('/Users/SandyMKawano/Desktop/Research/Forelimbs on Force Plates/Force Data/R Analysis/Step 3 Post Filter Processing/Filtered without Increments')
  # write.table(FilterInterp.NoIncrement, file=SaveFilteredNoIncrement, sep =",", row.names=FALSE)
  # 
  setwd('/Users/SandyMKawano/Desktop/Research/Forelimbs on Force Plates/Force Data/R Analysis/Step 3 Post Filter Processing/Peak net GRF')
  write.table(PeakNetGRF, file=SaveFilteredPeakNet, sep =",", row.names=FALSE)
  
  setwd('/Users/SandyMKawano/Desktop/Research/Forelimbs on Force Plates/Force Data/R Analysis/Step 3 Post Filter Processing/101 points')
  write.table(FilterInterp.NoIncrementWFill, file=SaveFiltered101, sep ="\t", row.names=FALSE)
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