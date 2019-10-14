################  GRF of Fish and Salamander Appendages ########################
## Code for evaluating force data collected from the Blob lab force plate #####
# Abbreviations: Vert = Vertical, ML = mediolateral, Hz = Horizontal (refers to anteroposterior direction)

# 7-23-11: Changed it so that the unconverted AP angle is saved instead for further analysis (converted ML angle still being used)

# 10-12-11: Changed it so that the net GRF data in units of BW is also saved in output

# 12-29-11: Changed code again, so it uses both the converted AP and converted ML angles

######################  POST-FILTERING PROCESSING #############################

# Clear everything in the R workspace, so nothing gets mixed up between different trials
rm(list=ls(all=TRUE))

# Read the most recent Video Info file

setwd('/Users/SandyMKawano/Desktop/Research/Forelimbs on Force Plates/Video/Video Info')
VideoPath <- setwd('/Users/SandyMKawano/Desktop/Research/Forelimbs on Force Plates/Video/Video Info')
Video.Ind <- list.files(VideoPath, pattern=".csv", full=TRUE)
Video.Dates <-substring(Video.Ind, 110, 115)
Video.Newest <- max(as.numeric(Video.Dates))
Video.Use <- Video.Ind[substring(Video.Ind, 110, 115) %in% Video.Newest]
VideoFile <- data.frame(read.csv(Video.Use, header=TRUE))
# Had issues when I was using read.table(); wasn't loading all the data correctly once I got to a certain size

# Set the directory, to facilitate choosing file
setwd('/Users/SandyMKawano/Desktop/Research/Forelimbs on Force Plates/Force Data/R Analysis/Step 2 Filter Data/Filtered Interpolated')

myFile <- file.choose()
myData <- read.csv(myFile, header=TRUE)

# Determining trial name and data collectiond date
Trial <- substring(myFile, 127, 133)

# Stamping with today's date
today <- Sys.Date()
SaveDate <- format(today, format="%y%m%d")


# Looking up video info
# Appendages listed as "Both" have both pectoral and pelvic appendage data
VideoInfo <- as.data.frame(VideoFile[VideoFile$File.name %in% Trial,])
VideoInfo$Pectoral.Start.Frame <- as.numeric(VideoInfo$Pectoral.Start.Frame)
VideoInfo$Pectoral.End.Frame <- as.numeric(VideoInfo$Pectoral.End.Frame)
VideoInfo$Pelvic.Start.Frame <- as.numeric(VideoInfo$Pelvic.Start.Frame)
VideoInfo$Pelvic.End.Frame <- as.numeric(VideoInfo$Pelvic.End.Frame)

# Identifying whether it's all of the filtered data or the filtered data interpolated to 101 points (0%-100% of stance)
# The interpolated data will always have 102 rows (101 points + 1 row for variable names), so if the data set has more than 102
# rows, then it must be the data set for all of the filtered data for that trial
ifelse ((nrow(myData) > 102), FilterAll <- myData, FilterInterp <- myData)


#######  For files that had the data interpolated to 101 points ################
  names(FilterInterp) <- c('PercentStance', 'InterpV.N', 'InterpML.N', 'InterpHz.N')

  FilterInterp$SquaredV.N <- FilterInterp$InterpV.N^2
  FilterInterp$SquaredML.N <- FilterInterp$InterpML.N^2
  FilterInterp$SquaredHz.N <- FilterInterp$InterpHz.N^2
  FilterInterp$NetGRF.N <- sqrt(FilterInterp$SquaredV.N+FilterInterp$SquaredML.N+FilterInterp$SquaredHz.N)
  FilterInterp$APAngle <- (acos(FilterInterp$InterpHz.N/(sqrt(FilterInterp$SquaredHz.N+FilterInterp$SquaredV.N))))*(180/pi)
  FilterInterp$APAngle.Convert <- 90-FilterInterp$APAngle
  FilterInterp$MLAngle <- (acos(FilterInterp$InterpML.N/(sqrt(FilterInterp$SquaredML.N+FilterInterp$SquaredV.N))))*(180/pi)
  FilterInterp$MLAngle.Convert <- 90-FilterInterp$MLAngle

  # Converting some data in Newtons to Body Weight (to help size standardize force data)
  BW.N <- VideoInfo$Body.Weight.kg*9.8 # converting animal's body weight from kilograms to Newtons
  FilterInterp$InterpV.BW <- FilterInterp$InterpV.N/BW.N
  FilterInterp$InterpML.BW <- FilterInterp$InterpML.N/BW.N
  FilterInterp$InterpHz.BW <- FilterInterp$InterpHz.N/BW.N
  FilterInterp$NetGRF.BW <- FilterInterp$NetGRF.N/BW.N

   ###################  ONLY PECTORAL APPENDAGE HITTING #################
  if (substring(myFile,143,145) == 'Pec' & is.na(VideoInfo$Pelvic.Start.Frame)==TRUE) # For files that have the pectoral appendage hitting the plate 1st
  {
    LimbCycleLength <- VideoInfo$Pectoral.End.Frame-VideoInfo$Pectoral.Start.Frame
      if (!is.na(VideoInfo$Body.1.Start))
        {
        if (VideoInfo$Body.1.Start>=VideoInfo$Pectoral.Start.Frame & VideoInfo$Body.1.Start<=VideoInfo$Pectoral.End.Frame)
          {
          Body1Start <- round(((VideoInfo$Body.1.Start-VideoInfo$Pectoral.Start.Frame)/LimbCycleLength)*VideoInfo$Filming.Rate.Hz)
          ImpPoints.Body1Start <- data.frame(Body1Start, FilterInterp$InterpV.BW[Body1Start+1], FilterInterp$InterpML.BW[Body1Start+1], FilterInterp$InterpHz.BW[Body1Start+1], FilterInterp$NetGRF.BW[Body1Start+1], row.names='Body 1 Start')
          # need to include +1 within square brackets b/c obs rows are 1 value greater than the % stance (i.e., 0% stance = row 1)
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
    # Taking data points at every 5% of stance
    IncrementRows <- seq(1,101,5) # row 1 = 0% stance, so essentially need to do row#-1 for each percentage of stance
    FilterInterp.NoIncrement <- FilterInterp[,c(1,10,12:16)]
    FilterInterp.Increment <- FilterInterp[IncrementRows,c(1,10,12:16)]

    # Evaluating at peak/max net GRF
    # Also making sure that I'm only including data points that occurred before overlap with other structures started
    PeakNetGRF <- FilterInterp[which(FilterInterp$NetGRF.BW==max(FilterInterp$NetGRF.BW)),]

    # Creating file names for saving as .csv
    SaveFilteredNoIncrement <- paste(Trial,"_Pec_NoIncrement_", SaveDate, ".csv", sep="")
    SaveFilteredIncrement <- paste(Trial,"_Pec_Increment_", SaveDate, ".csv", sep="")
    SaveFilteredPeakNet<- paste(Trial,"_Pec_PeakNetGRF_", SaveDate, ".csv", sep="")
    SaveFiltered101 <- paste(Trial, "_Pec_101_", SaveDate, ".txt", sep="")
    
    # Set the directory to save the data files
setwd('/Users/SandyMKawano/Desktop/Research/Forelimbs on Force Plates/Force Data/R Analysis/Step 3 Post Filter Processing/Filtered at Increments')
write.table(FilterInterp.Increment, file=SaveFilteredIncrement, sep =",", row.names=FALSE)

setwd('/Users/SandyMKawano/Desktop/Research/Forelimbs on Force Plates/Force Data/R Analysis/Step 3 Post Filter Processing/Filtered without Increments')
write.table(FilterInterp.NoIncrement, file=SaveFilteredNoIncrement, sep =",", row.names=FALSE)

setwd('/Users/SandyMKawano/Desktop/Research/Forelimbs on Force Plates/Force Data/R Analysis/Step 3 Post Filter Processing/Peak net GRF')
write.table(PeakNetGRF, file=SaveFilteredPeakNet, sep =",", row.names=FALSE)
    
    setwd('/Users/SandyMKawano/Desktop/Research/Forelimbs on Force Plates/Force Data/R Analysis/Step 3 Post Filter Processing/101 points')
    write.table(FilterInterp, file=SaveFiltered101, sep ="\t", row.names=FALSE)
  }
  
  
  

  ###################  PECTORAL APPENDAGE BEFORE PELVIC APPENDAGE FOR PECTORAL FILES #################

  if (substring(myFile,143,145) == 'Pec' & is.na(VideoInfo$Pelvic.Start.Frame)==FALSE & as.numeric(VideoInfo$Pectoral.Start.Frame)<as.numeric(VideoInfo$Pelvic.Start.Frame)) # For files that have the pectoral appendage hitting the plate 1st
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
  
  

###################  PECTORAL APPENDAGE AFTER PELVIC APPENDAGE FOR PECTORAL FILES #################

  if (substring(myFile,143,145) == 'Pec' & is.na(VideoInfo$Pelvic.Start.Frame)==FALSE & as.numeric(VideoInfo$Pectoral.Start.Frame)>as.numeric(VideoInfo$Pelvic.Start.Frame))
  {
      as.numeric(VideoInfo$Pectoral.Start.Frame)
      as.numeric(VideoInfo$Pectoral.End.Frame)
      as.numeric(VideoInfo$Pelvic.Start.Frame)
      as.numeric(VideoInfo$Pelvic.End.Frame)
      LimbCycleLength <- VideoInfo$Pectoral.End.Frame-VideoInfo$Pectoral.Start.Frame
      OverlapEnd <- ceiling(((VideoInfo$Pelvic.End.Frame-VideoInfo$Pectoral.Start.Frame)/LimbCycleLength)*VideoInfo$Filming.Rate.Hz) # Rounding up because want data after overlap is done
      ImpPoints.Overlap <- data.frame(OverlapEnd, FilterInterp$InterpV.BW[OverlapEnd+1], FilterInterp$InterpML.BW[OverlapEnd+1], FilterInterp$InterpHz.BW[OverlapEnd+1], FilterInterp$NetGRF.BW[OverlapEnd+1], row.names='Overlap End')
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
    UsableRows <- FilterInterp$PercentStance[-c(0:OverlapEnd)]+1 # overlap tends to occur towards the beginning # adding 1 b/c data rows don't start until row 2
    StanceWithoutOverlap.Increment <- IncrementRows[IncrementRows %in% UsableRows]

    FilterInterp.NoIncrement <- FilterInterp[UsableRows,c(1,10,12:16)]
    FilterInterp.Increment <- FilterInterp[StanceWithoutOverlap.Increment,c(1,10,12:16)]
    
    # Making cells within the overlap as "NA" so the data.frame maintains the same dimensions with all 0 -> 101 pts at 5% increments included
    UnusableRows <- FilterInterp$PercentStance[c(0:(OverlapEnd))]+1
    StanceWithOverlap.Increment <- IncrementRows[IncrementRows %in% UnusableRows]
    NA.Increment <- rep("NA", length(StanceWithOverlap.Increment))
    StanceWithOverlap.IncrementPercent <- StanceWithOverlap.Increment-1
    NA.Fillers <- data.frame(cbind(StanceWithOverlap.IncrementPercent, NA.Increment, NA.Increment, NA.Increment, NA.Increment, NA.Increment, NA.Increment))
    names(NA.Fillers) <- names(FilterInterp.Increment) # need to have the same variable names to rbind

    # Adding filler NA's to actual data to keep/analyze
  FilterInterp.IncrementWFill <- t(cbind(t(NA.Fillers), t(FilterInterp.Increment))) # simply using rbind didn't work

    # Evaluating at peak/max net GRF
    # Also making sure that I'm only including data points that occurred before overlap with other structures started
    PeakNetGRF <- FilterInterp[which(FilterInterp$NetGRF.BW==max(FilterInterp$NetGRF.BW[-c(1:OverlapEnd)])),]


    # Creating file names for saving as .csv
    SaveFilteredNoIncrement <- paste(Trial,"_Pec_NoIncrement_", SaveDate, ".csv", sep="")
    SaveFilteredIncrement <- paste(Trial,"_Pec_Increment_", SaveDate, ".csv", sep="")
    SaveFilteredPeakNet<- paste(Trial,"_Pec_PeakNetGRF_", SaveDate, ".csv", sep="")
    SaveFiltered101 <- paste(Trial, "_Pec_101_", SaveDate, ".txt", sep="")
    
    # Set the directory to save the data files
setwd('/Users/SandyMKawano/Desktop/Research/Forelimbs on Force Plates/Force Data/R Analysis/Step 3 Post Filter Processing/Filtered at Increments')
write.table(FilterInterp.IncrementWFill, file=SaveFilteredIncrement, sep =",", row.names=FALSE)

setwd('/Users/SandyMKawano/Desktop/Research/Forelimbs on Force Plates/Force Data/R Analysis/Step 3 Post Filter Processing/Filtered without Increments')
write.table(FilterInterp.NoIncrement, file=SaveFilteredNoIncrement, sep =",", row.names=FALSE)

setwd('/Users/SandyMKawano/Desktop/Research/Forelimbs on Force Plates/Force Data/R Analysis/Step 3 Post Filter Processing/Peak net GRF')
write.table(PeakNetGRF, file=SaveFilteredPeakNet, sep =",", row.names=FALSE)
      
      setwd('/Users/SandyMKawano/Desktop/Research/Forelimbs on Force Plates/Force Data/R Analysis/Step 3 Post Filter Processing/101 points')
      write.table(FilterInterp, file=SaveFiltered101, sep ="\t", row.names=FALSE)
  }
  
  
  

     ###################  ONLY PELVIC APPENDAGE HITTING PLATE #################
  if (substring(myFile,143,145) == 'Pel' & is.na(VideoInfo$Pectoral.Start.Frame==TRUE))
  {
    LimbCycleLength <- VideoInfo$Pelvic.End.Frame-VideoInfo$Pelvic.Start.Frame
      if (!is.na(VideoInfo$Body.1.Start))
        {
        if (VideoInfo$Body.1.Start>=VideoInfo$Pelvic.Start.Frame & VideoInfo$Body.1.Start<=VideoInfo$Pelvic.End.Frame)
          {
          Body1Start <- round(((VideoInfo$Body.1.Start-VideoInfo$Pelvic.Start.Frame)/LimbCycleLength)*VideoInfo$Filming.Rate.Hz)
          ImpPoints.Body1Start <- data.frame(Body1Start, FilterInterp$InterpV.BW[Body1Start+1], FilterInterp$InterpML.BW[Body1Start+1], FilterInterp$InterpHz.BW[Body1Start+1], FilterInterp$NetGRF.BW[Body1Start+1], row.names='Body 1 Start')
          names(ImpPoints.Body1Start) <- c('Percent Stance', 'Vertical (BW)', 'Mediolateral (BW)', 'Horizontal (BW)', 'Net GRF (BW)')
          }
        }
      if (!is.na(VideoInfo$Body.1.End))
        {
        if (VideoInfo$Body.1.End>=VideoInfo$Pelvic.Start.Frame & VideoInfo$Body.1.End<=VideoInfo$Pelvic.End.Frame)
          {
          Body1End <- round(((VideoInfo$Body.1.End-VideoInfo$Pelvic.Start.Frame)/LimbCycleLength)*VideoInfo$Filming.Rate.Hz)
          ImpPoints.Body1End <- data.frame(Body1End, FilterInterp$InterpV.BW[Body1End+1], FilterInterp$InterpML.BW[Body1End+1], FilterInterp$InterpHz.BW[Body1End+1], FilterInterp$NetGRF.BW[Body1End+1], row.names='Body 1 End')
          names(ImpPoints.Body1End) <- c('Percent Stance', 'Vertical (BW)', 'Mediolateral (BW)', 'Horizontal (BW)', 'Net GRF (BW)')
          }
        }

      if (!is.na(VideoInfo$Body.2.Start))
        {
        if (VideoInfo$Body.2.Start>=VideoInfo$Pelvic.Start.Frame & VideoInfo$Body.2.Start<=VideoInfo$Pelvic.End.Frame)
          {
          Body2Start <- round(((VideoInfo$Body.2.Start-VideoInfo$Pelvic.Start.Frame)/LimbCycleLength)*VideoInfo$Filming.Rate.Hz)
          ImpPoints.Body2Start <- data.frame(Body2Start, FilterInterp$InterpV.BW[Body2Start+1], FilterInterp$InterpML.BW[Body2Start+1], FilterInterp$InterpHz.BW[Body2Start+1], FilterInterp$NetGRF.BW[Body2Start+1], row.names='Body 2 Start')
          names(ImpPoints.Body2Start) <- c('Percent Stance', 'Vertical (BW)', 'Mediolateral (BW)', 'Horizontal (BW)', 'Net GRF (BW)')
          }
        }
      if (!is.na(VideoInfo$Body.2.End))
        {
        if (VideoInfo$Body.2.End>=VideoInfo$Pelvic.Start.Frame & VideoInfo$Body.2.End<=VideoInfo$Pelvic.End.Frame)
          {
          Body2End <- round(((VideoInfo$Body.2.End-VideoInfo$Pelvic.Start.Frame)/LimbCycleLength)*VideoInfo$Filming.Rate.Hz)
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


    # Taking data points at every 5% of stance
    IncrementRows <- seq(1,101,5) # row 1 = 0% stance, so essentially need to do row#-1 for each percentage of stance
    FilterInterp.NoIncrement <- FilterInterp[,c(1,10,12:16)]
    FilterInterp.Increment <- FilterInterp[IncrementRows,c(1,10,12:16)]

    # Evaluating at peak/max net GRF
    # Also making sure that I'm only including data points that occurred before overlap with other structures started
    PeakNetGRF <- FilterInterp[which(FilterInterp$NetGRF.BW==max(FilterInterp$NetGRF.BW)),]

    # Creating file names for saving as .csv
    SaveFilteredNoIncrement <- paste(Trial,"_Pel_NoIncrement_", SaveDate, ".csv", sep="")
    SaveFilteredIncrement <- paste(Trial,"_Pel_Increment_", SaveDate, ".csv", sep="")
    SaveFilteredPeakNet<- paste(Trial,"_Pel_PeakNetGRF_", SaveDate, ".csv", sep="")
    SaveFiltered101 <- paste(Trial, "_Pel_101_", SaveDate, ".txt", sep="")
    
    # Set the directory to save the data files
setwd('/Users/SandyMKawano/Desktop/Research/Forelimbs on Force Plates/Force Data/R Analysis/Step 3 Post Filter Processing/Filtered at Increments')
write.table(FilterInterp.Increment, file=SaveFilteredIncrement, sep =",", row.names=FALSE)

setwd('/Users/SandyMKawano/Desktop/Research/Forelimbs on Force Plates/Force Data/R Analysis/Step 3 Post Filter Processing/Filtered without Increments')
write.table(FilterInterp.NoIncrement, file=SaveFilteredNoIncrement, sep =",", row.names=FALSE)

setwd('/Users/SandyMKawano/Desktop/Research/Forelimbs on Force Plates/Force Data/R Analysis/Step 3 Post Filter Processing/Peak net GRF')
write.table(PeakNetGRF, file=SaveFilteredPeakNet, sep =",", row.names=FALSE)

    setwd('/Users/SandyMKawano/Desktop/Research/Forelimbs on Force Plates/Force Data/R Analysis/Step 3 Post Filter Processing/101 points')
    write.table(FilterInterp, file=SaveFiltered101, sep ="\t", row.names=FALSE)
  }

    ###################  PECTORAL APPENDAGE BEFORE PELVIC APPENDAGE FOR PELVIC TRIALS #################

  if (substring(myFile,143,145) == 'Pel' & !is.na(VideoInfo$Pectoral.Start.Frame)==TRUE & as.numeric(VideoInfo$Pectoral.Start.Frame)<as.numeric(VideoInfo$Pelvic.Start.Frame))
  {
      LimbCycleLength <- VideoInfo$Pelvic.End.Frame-VideoInfo$Pelvic.Start.Frame
      OverlapEnd <- ceiling(((VideoInfo$Pectoral.End.Frame-VideoInfo$Pelvic.Start.Frame)/LimbCycleLength)*VideoInfo$Filming.Rate.Hz) # Rounding up because want data after overlap is done
      ImpPoints.Overlap <- data.frame(OverlapEnd, FilterInterp$InterpV.BW[OverlapEnd+1], FilterInterp$InterpML.BW[OverlapEnd+1], FilterInterp$InterpHz.BW[OverlapEnd+1], FilterInterp$NetGRF.BW[OverlapEnd+1], row.names='Overlap End')
      names(ImpPoints.Overlap) <- c('Percent Stance', 'Vertical (BW)', 'Mediolateral (BW)', 'Horizontal (BW)', 'Net GRF (BW)')

      if (!is.na(VideoInfo$Body.1.Start))
        {
        if (VideoInfo$Body.1.Start>=VideoInfo$Pelvic.Start.Frame & VideoInfo$Body.1.Start<=VideoInfo$Pelvic.End.Frame)
          {
          Body1Start <- round(((VideoInfo$Body.1.Start-VideoInfo$Pelvic.Start.Frame)/LimbCycleLength)*VideoInfo$Filming.Rate.Hz)
          ImpPoints.Body1Start <- data.frame(Body1Start, FilterInterp$InterpV.BW[Body1Start+1], FilterInterp$InterpML.BW[Body1Start+1], FilterInterp$InterpHz.BW[Body1Start+1], FilterInterp$NetGRF.BW[Body1Start+1], row.names='Body 1 Start')
          names(ImpPoints.Body1Start) <- c('Percent Stance', 'Vertical (BW)', 'Mediolateral (BW)', 'Horizontal (BW)', 'Net GRF (BW)')
          }
        }
       if (!is.na(VideoInfo$Body.1.End))
        {
        if (VideoInfo$Body.1.End>=VideoInfo$Pelvic.Start.Frame & VideoInfo$Body.1.End<=VideoInfo$Pelvic.End.Frame)
          {
          Body1End <- round(((VideoInfo$Body.1.End-VideoInfo$Pelvic.Start.Frame)/LimbCycleLength)*VideoInfo$Filming.Rate.Hz)
          ImpPoints.Body1End <- data.frame(Body1End, FilterInterp$InterpV.BW[Body1End+1], FilterInterp$InterpML.BW[Body1End+1], FilterInterp$InterpHz.BW[Body1End+1], FilterInterp$NetGRF.BW[Body1End+1], row.names='Body 1 End')
          names(ImpPoints.Body1End) <- c('Percent Stance', 'Vertical (BW)', 'Mediolateral (BW)', 'Horizontal (BW)', 'Net GRF (BW)')
          }
        }
       if (!is.na(VideoInfo$Body.2.Start))
        {
        if (VideoInfo$Body.2.Start>=VideoInfo$Pelvic.Start.Frame & VideoInfo$Body.2.Start<=VideoInfo$Pelvic.End.Frame)
          {
          Body2Start <- round(((VideoInfo$Body.2.Start-VideoInfo$Pelvic.Start.Frame)/LimbCycleLength)*VideoInfo$Filming.Rate.Hz)
          ImpPoints.Body2Start <- data.frame(Body2Start, FilterInterp$InterpV.BW[Body2Start+1], FilterInterp$InterpML.BW[Body2Start+1], FilterInterp$InterpHz.BW[Body2Start+1], FilterInterp$NetGRF.BW[Body2Start+1], row.names='Body 2 Start')
          names(ImpPoints.Body2Start) <- c('Percent Stance', 'Vertical (BW)', 'Mediolateral (BW)', 'Horizontal (BW)', 'Net GRF (BW)')
          }
        }
       if (!is.na(VideoInfo$Body.2.End))
        {
        if (VideoInfo$Body.2.End>=VideoInfo$Pelvic.Start.Frame & VideoInfo$Body.2.End<=VideoInfo$Pelvic.End.Frame)
          {
          Body2End <- round(((VideoInfo$Body.2.End-VideoInfo$Pelvic.Start.Frame)/LimbCycleLength)*VideoInfo$Filming.Rate.Hz)
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
    UsableRows <- FilterInterp$PercentStance[-c(0:OverlapEnd)]+1 # forelimb overlap tends to occur towards the beginning # adding 1 b/c data rows don't start until row 2
    StanceWithoutOverlap.Increment <- IncrementRows[IncrementRows %in% UsableRows]
    FilterInterp.NoIncrement <- FilterInterp[UsableRows,c(1,10,12:16)]
    FilterInterp.Increment <- FilterInterp[StanceWithoutOverlap.Increment,c(1,10,12:16)]
    
    # Making cells within the overlap as "NA" so the data.frame maintains the same dimensions with all 0 -> 101 pts at 5% increments included
    UnusableRows <- FilterInterp$PercentStance[c(0:(OverlapEnd))]+1
    StanceWithOverlap.Increment <- IncrementRows[IncrementRows %in% UnusableRows]
    NA.Increment <- rep("NA", length(StanceWithOverlap.Increment))
    StanceWithOverlap.IncrementPercent <- StanceWithOverlap.Increment-1
    NA.Fillers <- data.frame(cbind(StanceWithOverlap.IncrementPercent, NA.Increment, NA.Increment, NA.Increment, NA.Increment, NA.Increment, NA.Increment))
    names(NA.Fillers) <- names(FilterInterp.Increment) # need to have the same variable names to rbind

    # Adding filler NA's to actual data to keep/analyze
    FilterInterp.IncrementWFill <- t(cbind(t(NA.Fillers), t(FilterInterp.Increment))) # simply using rbind didn't work

    # Evaluating at peak/max net GRF
    # Also making sure that I'm only including data points that occurred before overlap with other structures started
    PeakNetGRF <- FilterInterp[which(FilterInterp$NetGRF.BW==max(FilterInterp$NetGRF.BW[-c(1:OverlapEnd)])),]


    # Creating file names for saving as .csv
    SaveFilteredNoIncrement <- paste(Trial,"_Pel_NoIncrement_", SaveDate, ".csv", sep="")
    SaveFilteredIncrement <- paste(Trial,"_Pel_Increment_", SaveDate, ".csv", sep="")
    SaveFilteredPeakNet<- paste(Trial,"_Pel_PeakNetGRF_", SaveDate, ".csv", sep="")
    SaveFiltered101 <- paste(Trial, "_Pel_101_", SaveDate, ".txt", sep="")
    
    # Set the directory to save the data files
setwd('/Users/SandyMKawano/Desktop/Research/Forelimbs on Force Plates/Force Data/R Analysis/Step 3 Post Filter Processing/Filtered at Increments')
write.table(FilterInterp.IncrementWFill, file=SaveFilteredIncrement, sep =",", row.names=FALSE)

setwd('/Users/SandyMKawano/Desktop/Research/Forelimbs on Force Plates/Force Data/R Analysis/Step 3 Post Filter Processing/Filtered without Increments')
write.table(FilterInterp.NoIncrement, file=SaveFilteredNoIncrement, sep =",", row.names=FALSE)

setwd('/Users/SandyMKawano/Desktop/Research/Forelimbs on Force Plates/Force Data/R Analysis/Step 3 Post Filter Processing/Peak net GRF')
write.table(PeakNetGRF, file=SaveFilteredPeakNet, sep =",", row.names=FALSE)
      
setwd('/Users/SandyMKawano/Desktop/Research/Forelimbs on Force Plates/Force Data/R Analysis/Step 3 Post Filter Processing/101 points')
write.table(FilterInterp, file=SaveFiltered101, sep ="\t", row.names=FALSE)
  }





    # Plotting the data
    quartz()
     par(mfrow=c(2,2), oma = c(3, 0, 2, 0))  # oma = outer margin with 2 lines above the top of the graphs
      # Vertical component of GRF graph
      plot(FilterInterp$PercentStance, FilterInterp$InterpV.BW, xlab='Percent Stance', ylab='GRF - Vertical (BW)', main='Zeroed GRF (Vertical) Force', type="l", col="blue")
      if (exists("ImpPoints"))
      {
        points(ImpPoints[,1], ImpPoints[,2], type='p', pch='O', col='cyan')
        text(ImpPoints[,1], ImpPoints[,2], labels=row.names(ImpPoints), pos=3, font=2) # pos: 1 = below, 2 = left, 3 = above, 4 = right
      }
      abline(v=PeakNetGRF$PercentStance, col='magenta', lty=2, lwd=2) # Plotting vertical line at Peak Net GRF
      # Mediolateral component of GRF graph
      plot(FilterInterp$PercentStance, FilterInterp$InterpML.BW, xlab='Percent Stance', ylab='GRF - Mediolateral (BW)', main='Zeroed GRF (Mediolateral) Force', type="l", col="red")
      if (exists("ImpPoints"))
      {
        points(ImpPoints[,1], ImpPoints[,3], type='p', pch='O', col='cyan')
        text(ImpPoints[,1], ImpPoints[,3], labels=row.names(ImpPoints), pos=3, font=2) # pos: 1 = below, 2 = left, 3 = above, 4 = right
      }
      abline(v=PeakNetGRF$PercentStance, col='magenta', lty=2, lwd=2)
      # Horizontal (Anteroposterior) component of GRF graph
      plot(FilterInterp$PercentStance, FilterInterp$InterpHz.BW, xlab='Percent Stance', ylab='GRF - Horizontal (BW)', main='Zeroed GRF (Horizontal) Force', type="l", col="forestgreen")
      if (exists("ImpPoints"))
      {
        points(ImpPoints[,1], ImpPoints[,4], type='p', pch='O', col='cyan')
        text(ImpPoints[,1], ImpPoints[,4], labels=row.names(ImpPoints), pos=3, font=2) # pos: 1 = below, 2 = left, 3 = above, 4 = right
      }
      abline(v=PeakNetGRF$PercentStance, col='magenta', lty=2, lwd=2)
      # Net GRF graph
      plot(FilterInterp$PercentStance, FilterInterp$NetGRF.BW, xlab='Percent Stance', ylab='Net GRF (BW)', main='Zeroed Net GRF Force', type="l", col="purple")
      if (exists("ImpPoints"))
      {
        points(ImpPoints[,1], ImpPoints[,5], type='p', pch='O', col='cyan')
        text(ImpPoints[,1], ImpPoints[,5], labels=row.names(ImpPoints), pos=3, font=2) # pos: 1 = below, 2 = left, 3 = above, 4 = right
      }
      abline(v=PeakNetGRF$PercentStance, col='magenta', lty=2, lwd=2)
      GraphTitle <- paste(Trial, "-", substring(myFile, 143, 145), sep=" ")
      mtext(GraphTitle, line=0.5, outer=TRUE)  # writes an overall title over the graphs
      mtext('Dashed pink line = % Stance for Peak Net GRF', side=1, outer=TRUE, col = 'magenta')
      mtext('Blue circles = points of overlap with other structures', side=1, line=1.5, outer=TRUE, col = 'cyan')

    # Saving the graphs
    PdfSave <- paste(substring(myFile, 1, 85), "Step 3 Post Filter Processing/Graphs/", Trial, "_", substring(myFile, 143, 145), "_PostFilter_", SaveDate, ".pdf", sep="")
    pdf(PdfSave)
    par(mfrow=c(2,2), oma = c(3, 0, 2, 0))  # oma = outer margin with 2 lines above the top of the graphs
      # Vertical component of GRF graph
      plot(FilterInterp$PercentStance, FilterInterp$InterpV.BW, xlab='Percent Stance', ylab='GRF - Vertical (BW)', main='Zeroed GRF (Vertical) Force', type="l", col="blue")
      if (exists("ImpPoints"))
      {
        points(ImpPoints[,1], ImpPoints[,2], type='p', pch='O', col='cyan')
        text(ImpPoints[,1], ImpPoints[,2], labels=row.names(ImpPoints), pos=3, font=2) # pos: 1 = below, 2 = left, 3 = above, 4 = right
      }
      abline(v=PeakNetGRF$PercentStance, col='magenta', lty=2, lwd=2) # Plotting vertical line at Peak Net GRF
      # Mediolateral component of GRF graph
      plot(FilterInterp$PercentStance, FilterInterp$InterpML.BW, xlab='Percent Stance', ylab='GRF - Mediolateral (BW)', main='Zeroed GRF (Mediolateral) Force', type="l", col="red")
      if (exists("ImpPoints"))
      {
        points(ImpPoints[,1], ImpPoints[,3], type='p', pch='O', col='cyan')
        text(ImpPoints[,1], ImpPoints[,3], labels=row.names(ImpPoints), pos=3, font=2) # pos: 1 = below, 2 = left, 3 = above, 4 = right
      }
      abline(v=PeakNetGRF$PercentStance, col='magenta', lty=2, lwd=2)
      # Horizontal (Anteroposterior) component of GRF graph
      plot(FilterInterp$PercentStance, FilterInterp$InterpHz.BW, xlab='Percent Stance', ylab='GRF - Horizontal (BW)', main='Zeroed GRF (Horizontal) Force', type="l", col="forestgreen")
      if (exists("ImpPoints"))
      {
        points(ImpPoints[,1], ImpPoints[,4], type='p', pch='O', col='cyan')
        text(ImpPoints[,1], ImpPoints[,4], labels=row.names(ImpPoints), pos=3, font=2) # pos: 1 = below, 2 = left, 3 = above, 4 = right
      }
      abline(v=PeakNetGRF$PercentStance, col='magenta', lty=2, lwd=2)
      # Net GRF graph
      plot(FilterInterp$PercentStance, FilterInterp$NetGRF.BW, xlab='Percent Stance', ylab='Net GRF (BW)', main='Zeroed Net GRF Force', type="l", col="purple")
      if (exists("ImpPoints"))
      {
        points(ImpPoints[,1], ImpPoints[,5], type='p', pch='O', col='cyan')
        text(ImpPoints[,1], ImpPoints[,5], labels=row.names(ImpPoints), pos=3, font=2) # pos: 1 = below, 2 = left, 3 = above, 4 = right
      }
      abline(v=PeakNetGRF$PercentStance, col='magenta', lty=2, lwd=2)
      GraphTitle <- paste(Trial, "-", substring(myFile, 143, 145), sep=" ")
      mtext(GraphTitle, line=0.5, outer=TRUE)  # writes an overall title over the graphs
      mtext('Dashed pink line = % Stance for Peak Net GRF', side=1, outer=TRUE, col = 'magenta')
      mtext('Blue circles = points of overlap with other structures', side=1, line=1.5, outer=TRUE, col = 'cyan')
      dev.off()








