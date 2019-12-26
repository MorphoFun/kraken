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
library(lme4) # for lmer()
library(MuMIn) # for r.squaredGLMM()
library(emmeans) # for emmean() and contrast()
library(piecewiseSEM) # for rsquared()
library(gridExtra) # for grid.arrange()
library(car) # for qqPlot()
library(ggplot2) # for ggplot()
library(cowplot) # for ggdraw()
library(qqplotr) # for stat_qq_band()
library(robustlmm) # for rlmer()
library(nlme) # for lme() to allow unequal variances across groups in the random effects


library(EMAtools) # for lme.dscore(); although, may not use this
library(effsize) # for cohen.d(); although, may not use this


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

  
  ### Combine the peak net GRF data
  GRFs$Pelvic$Pel_GRFs_Filtered_PeakNet$appendage <- "pelvic"

  GRFs$Pectoral$Pec_GRFs_Filtered_PeakNet$appendage <- "pectoral"  

  
  peakNetGRF <- rbind(GRFs$Pelvic$Pel_GRFs_Filtered_PeakNet, GRFs$Pectoral$Pec_GRFs_Filtered_PeakNet)
  peakNetGRF$group <- paste(substring(peakNetGRF$filename, 1, 2), substring(peakNetGRF$appendage, 1, 3), sep = "_")
  peakNetGRF$individual <- substring(peakNetGRF$filename, 1, 4)
  
  ## Subsetting data into groups to analyze
  pec_peakNetGRFs <- subset(peakNetGRF, appendage == "pectoral")
  pel_peakNetGRFs <- subset(peakNetGRF, appendage == "pelvic")
  sal_peakNetGRFs <- subset(peakNetGRF, !substring(filename, 1, 2) == "pb")
  
  
  #### PeakNetGRF: summary stats ####
  variablesToAnalyze <- (c("PercentStance", "InterpV_BW", "InterpML_BW", "InterpAP_BW", "NetGRF_BW", "MLAngle_Convert_deg", "APAngle_Convert_deg", "group", "individual", "appendage"))
  aggregate(. ~ group, data = peakNetGRF[,variablesToAnalyze[1:(length(variablesToAnalyze)-2)]], FUN = function(x) c(mean = mean(x), sd = sd(x), n = length(x)))
  nVars <- 7
  
  
  #### YANK ####
  yank_pec <- list(
    vertical = lapply(GRFs$Pectoral$Pec_GRFs_Filtered_dataset_noOverlap, function(x) yank(x$PercentStance, x$InterpV_BW)),
    medioateral = lapply(GRFs$Pectoral$Pec_GRFs_Filtered_dataset_noOverlap, function(x) yank(x$PercentStance, x$InterpML_BW)),
    anteroposterior = lapply(GRFs$Pectoral$Pec_GRFs_Filtered_dataset_noOverlap, function(x) yank(x$PercentStance, x$InterpAP_BW)),
    net = lapply(GRFs$Pectoral$Pec_GRFs_Filtered_dataset_noOverlap, function(x) yank(x$PercentStance, x$NetGRF_BW))
  )
  
  yank_pel <- list(
    vertical = lapply(GRFs$Pelvic$Pec_GRFs_Filtered_dataset_noOverlap, function(x) yank(x$PercentStance, x$InterpV_BW)),
    medioateral = lapply(GRFs$Pelvic$Pec_GRFs_Filtered_dataset_noOverlap, function(x) yank(x$PercentStance, x$InterpML_BW)),
    anteroposterior = lapply(GRFs$Pelvic$Pec_GRFs_Filtered_dataset_noOverlap, function(x) yank(x$PercentStance, x$InterpAP_BW)),
    net = lapply(GRFs$Pelvic$Pec_GRFs_Filtered_dataset_noOverlap, function(x) yank(x$PercentStance, x$NetGRF_BW))
  )
  

  #### REMOVING OUTLIERS ####
  
  ## ggplot2 calculates more outliers bc baseplot::boxplot doesn't actually calculate the 1st and 3rd quantiles with even n
  # https://stackoverflow.com/questions/21793715/why-geom-boxplot-identify-more-outliers-than-base-boxplot

  # This function returns errors about 'subscript out of bounds' when the outlierRange is different from 1.5
removeOutliers <- function(df, yName, group, labelName, outlierRange, ... ){
  bp <- vector("list", nVars)
  outliers <- vector("list", nVars)
  usableData <- list()
  for(i in 1:length(yName)){
    # changing the range affects how far beyond the IQR is considered an outlier (1.5 set as default)
    bp[[i]] <- car::Boxplot(df[,yName[i]] ~ group, id.method = labelName, data = df, ylab = yName[i], range = outlierRange)
    outliers[[i]] <- df[bp[[i]],]
  }
  outliersCombined <- data.frame(do.call("rbind", outliers))
  outliersUnique <- unique(outliersCombined)
  usableData <- df[ ! df[,labelName] %in% outliersUnique$filename, ]

  output <- list(
    outliers = outliersUnique,
    usableData = usableData
  )
  return(output)
}

### identifying the outliers across all of the variables for the pectoral peak net GRF dataset
## changing outlier range to 2 so 'outliers' are points falling 2x away from the IQR
pec_peakNetGRF_noOutliers <- removeOutliers(pec_peakNetGRFs, variablesToAnalyze[1:7], "group", "filename", outlierRange = 2)



  #### LINEAR MIXED EFFECTS MODELS ####
  ## This will be used to compare the means between groups while accounting for the repeated trials within individuals
  ## Since there are no pelvic data for Periophthalmus, we'll have three models:
  ## 1) pectoral comparison between Periophthalmus, Ambystoma, and Pleurodeles
  ## 2) pelvic comparison between Amvbystoma and Pleurodeles
  ## 3) pectoral versus pelvic for Ambystoma and Pleurodeles
  ## The effect sizes of individual independent variables can be assessed through the fixed effects: 
  ## https://stat.ethz.ch/pipermail/r-sig-mixed-models/2013q4/021102.html
  ## Or, could consider the f2 value (Aiken and West 1991): https://largescaleassessmentsineducation.springeropen.com/articles/10.1186/s40536-018-0061-2
  
  #### LMER with random intercepts: Peak net GRF - Pectoral ####
  
  ### a) Run lmers
  
  ## with outliers and without outliers
  modelFormulae <- list()
  pec_LMM <- list()
  pec_LMM_noOutliers <- list()
  pec_LMM_residuals <- list()
  pec_LMM_noOutliers_residuals <- list()
  for (i in 1:nVars) {
    modelFormulae[[i]] <- as.formula(paste(variablesToAnalyze[i], "~group+(1|individual)", sep = ""))
    pec_LMM[[i]] <- lmer(modelFormulae[[i]], data = pec_peakNetGRFs)
    pec_LMM_residuals[[i]] <- resid(pec_LMM[[i]])
    pec_LMM_noOutliers[[i]] <- lmer(modelFormulae[[i]], data = data.frame(pec_peakNetGRF_noOutliers$usableData))
    pec_LMM_noOutliers_residuals[[i]] <- resid(pec_LMM_noOutliers[[i]])
  }
  names(pec_LMM) <- modelFormulae
  names(pec_LMM_noOutliers) <- modelFormulae
  names(pec_LMM_residuals) <- variablesToAnalyze[1:7]
  names(pec_LMM_noOutliers_residuals) <- variablesToAnalyze[1:7]

  pec_LMM_resids <- data.frame(do.call("cbind", pec_LMM_residuals))
  # problems with this one because there are different number of observations in each variable
  pec_LMM_noOutliers_resids <- data.frame(do.call("cbind", pec_LMM_noOutliers_residuals))
  

  #### Testing the assumptions ####
  ## Don't need to test for linearity of data because the predictors are categorical

  ## b) evaluating the normality of the residuals
  # the null of the Shapiro-Wilk test is that the input (e.g., residuals of data) are normal
  
  # For data with outliers
  pec_LMM_shapiro <- list()
  for (i in 1:nVars) {
    pec_LMM_shapiro[[i]] <- shapiro.test(resid(pec_LMM[[i]]))
    qqPlot(resid(pec_LMM[[i]]), ylab = paste(names(pec_LMM)[[i]], " residuals"))
  }
  #do.call(grid.arrange, grobs = list(qqplots)) # can't use grid.arrange on car plots
  names(pec_LMM_shapiro) <- modelFormulae
  # Only InterpAP_BW met the assumption of normality based on the Shapiro-Wilk tests
  # However, the graphs show that the values tended to deviate from normality mainly because
  # the smallest points tended to underestimate the fitted values, which is more conservative
  # than having the largest values overestimate the qq-line
  
  # For data without outliers
  pec_LMM_noOutliers_shapiro <- list()
  for (i in 1:nVars) {
    pec_LMM_noOutliers_shapiro[[i]] <- shapiro.test(resid(pec_LMM_noOutliers[[i]]))
    qqPlot(resid(pec_LMM_noOutliers[[i]]), ylab = paste(names(pec_LMM_noOutliers)[[i]], " residuals"))
  }
  names(pec_LMM_noOutliers_shapiro) <- modelFormulae
  # removing the outliers made it so more variables met the assumption of normal residuals
  # only Percent Stance and ML_BW were not normal

  
  ## c) Testing homogeneity of variances
  # the Bartlett's test is more sensitive to non-normal data so people often use Levene's
  # more info here; http://www.sthda.com/english/wiki/compare-multiple-sample-variances-in-r
  # The Fligner-Killeen test can be used for non-normal data  
  
  # with full data set
  apply(pec_peakNetGRFs[,variablesToAnalyze[1:7]],2,function(x) {leveneTest(x ~ as.factor(pec_peakNetGRFs$group))})
  apply(pec_peakNetGRFs[,variablesToAnalyze[1:7]],2,function(x) {fligner.test(x ~ as.factor(pec_peakNetGRFs$group))})
  
  # can also evaluate the homogeneity of variances graphically
  pp <- list()
  for (i in 1:nVars) {
    pp[[i]] <- plot(pec_LMM[[i]])
  }
  names(pp) <- names(pec_LMM)
  # there does not appear to be any observable pattern in the residuals vs. fitted plots, which suggests
  # that the variances are homogeneous
  
  # Testing assumption when the outliers were removed
  apply(pec_peakNetGRF_noOutliers$usableData[,variablesToAnalyze[1:7]],2,function(x) {leveneTest(x ~ as.factor(pec_peakNetGRF_noOutliers$usableData$group))})
  apply(pec_peakNetGRF_noOutliers$usableData[,variablesToAnalyze[1:7]],2,function(x) {fligner.test(x ~ as.factor(pec_peakNetGRF_noOutliers$usableData$group))})
  
  ## d) Testing the normality of the random effects
  # following suggestions from: https://stats.stackexchange.com/questions/117170/testing-whether-random-effects-are-normally-distributed-in-r
  
  ## random intercepts model
  r_int<- ranef(pec_LMM[[1]])$individual$`(Intercept)`
  qqnorm(r_int)
  qqline(r_int)
  shapiro.test(r_int)

  
  ### NOTE: LMMs with random intercepts and slopes were attempted but the sample size of the pectoral data set were not large enough to handle the more 
  ### complex random effects structure for certain variable, so only a random intercepts LMM was used to make the comparisons consistent

  
  
  ## Zu omega squared to assess the 'goodness of fit' for the entire model
  # Xu's omega method: http://onlinelibrary.wiley.com/doi/10.1002/sim.1572/abstract
  ## or through the performance package: (got the same exact results as my code)
  # performance::r2_xu(pec_LMM)
  
  Xu_omega2 <- function(lmm, ...) {
    1-var(residuals(lmm))/(var(model.response(model.frame(lmm))))
  }
  
  pec_LMM_omega2 <- lapply(pec_LMM, FUN = function(x) Xu_omega2(x))

   
  ## Can also do post-hoc pair-wise comparisons for the fixed effects: https://stats.stackexchange.com/questions/237512/how-to-perform-post-hoc-test-on-lmer-model
  ## This describes that the differences between the post-hoc options: https://stats.stackexchange.com/questions/204741/which-multiple-comparison-method-to-use-for-a-lmer-model-lsmeans-or-glht
  ## The main difference is how they calculate the p-value, which doesn't matter to LMMs. 
  ## multcomp::glht() can produces shorter CIs and small p-values (which may be due to assumming infinite dfs),
  ## so might be better to use lsmeans::lsmeans() to be more conservative. 
  ## lsmeans also handles interactions better. 
  ## lsmeans has migrated to emmeans. This describes how to get the contrasts: 
  ## https://stats.stackexchange.com/questions/331238/post-hoc-pairwise-comparison-of-interaction-in-mixed-effects-lmer-model
  ## and https://stats.stackexchange.com/questions/355611/pairwise-comparisons-with-emmeans-for-a-mixed-three-way-interaction-in-a-linear
  # pec_EMM <- emmeans(pec_LMM[[1]], ~ group)
  # contrast(pec_EMM[[1]], interaction = "pairwise")
  
  pec_EMM <- lapply(pec_LMM, FUN = function(x) emmeans(x, ~ group))
  pec_EMM_contrasts <- lapply(pec_EMM, FUN = function(x) contrast(x, interaction = "pairwise"))
  
  ## Cohen's f2 to estimate the effect size of a fixed effect: https://github.com/finch-f/effect-size-in-Linear-mixed-model/blob/master/How%20to%20calculate%20effect%20size%20of%20a%20fixed%20effect.pdf
  Cohen_f2 <- function(y, FE = NULL, RE = NULL, df = NULL, ...) {
    LMM_formula <- as.formula(paste(y, " ~ ", FE, " + ", RE))
    LMM <- lmer(LMM_formula, data = df)
    LMM_null_formula <- as.formula(paste(y, " ~ 1 + ", RE))
    LMM_null <- lmer(LMM_null_formula, data = df)
    LMM_pseudoR2 <- piecewiseSEM::rsquared(LMM)
    LMM_null_pseudoR2 <- piecewiseSEM::rsquared(LMM_null)
    LMM_f2 <- (LMM_pseudoR2[,5:6] - LMM_null_pseudoR2[,5:6])/(1 - LMM_pseudoR2[,5:6])
    # double-check last line bc I got a conditional f2 that was lower than the marginal, which shouldn't be possible
    # maybe that's because I should only be calculating the f2 for the marginal values
  }
  
  pec_LMM_null <- lmer(PercentStance ~ 1 + (1 |individual), data = subset(peakNetGRF, appendage == "pectoral"), REML = TRUE)
  pec_LMM_pseudoR2 <- piecewiseSEM::rsquared(pec_LMM)
  pec_LMM_null_pseudoR2 <- piecewiseSEM::rsquared(pec_LMM_null)
  pec_LMM_F2 <- (pec_LMM_null_pseudoR2 - pec_LMM_pseudoR2)/(1 - pec_LMM_null_pseudoR2)
  # got a lot of warnings about "In Ops.factor(left, right) : ‘-’ not meaningful for factors"
  pec_LMM_F2 <- (pec_LMM_pseudoR2[,5:6] - pec_LMM_null_pseudoR2[,5:6])/(1 - pec_LMM_pseudoR2[,5:6])
  # this worked, though. Values were different from MuMIn::r.squaredGLMM(pec_LMM)
  # using the marginal values since Cohen's f2 is based on the fixed effects, and the conditional value includes both the fixed and random effects.
  
  ## when plotting the vertical and net GRF data, maybe use filled regions to help highlight the areas of overlap 
  
  
  ### Calculating the Inter Class Correlation ###
  ## Based on: https://www.ssc.wisc.edu/sscc/pubs/MM/MM_DiagInfer.html
  ## This represents the variation due to the random effect
  ## see info for "sc": https://www.rdocumentation.org/packages/lme4/versions/1.1-21/topics/VarCorr
  ## or use: https://www.rdocumentation.org/packages/sjstats/versions/0.17.4/topics/icc
  ## looks like sjstats::icc has been changed to performance::icc
  ## ICC might already be reported in the MuMIn pseudo-R2 code, though
  
  r1Var <- as.numeric(VarCorr(pec_LMM[[1]])[["individual"]])
  residVar <- attr(VarCorr(pec_LMM[[1]]), "sc")^2
  r1Var
  residVar
  r1Var / (r1Var + residVar) # same as the adjusted ICC from performance::icc()
  
  ## performance::icc(pec_LMM[[1]])
  
  
  #### Robust LMM ####
  
  pec_RLMM <- list()
  pec_RLMM_noOutliers <- list()
  pec_RLMM_residuals <- list()
  pec_RLMM_noOutliers_residuals <- list()
  for (i in 1:nVars) {
    pec_RLMM[[i]] <- rlmer(modelFormulae[[i]], data = pec_peakNetGRFs)
    pec_RLMM_residuals[[i]] <- resid(pec_RLMM[[i]])
    pec_RLMM_noOutliers[[i]] <- rlmer(modelFormulae[[i]], data = data.frame(pec_peakNetGRF_noOutliers$usableData))
    pec_RLMM_noOutliers_residuals[[i]] <- resid(pec_RLMM_noOutliers[[i]])
  }
  names(pec_RLMM) <- modelFormulae
  names(pec_RLMM_noOutliers) <- modelFormulae
  names(pec_RLMM_residuals) <- variablesToAnalyze[1:7]
  names(pec_RLMM_noOutliers_residuals) <- variablesToAnalyze[1:7]
  
  
  # to improve the efficiency of the random effect
  
  pec_RLMM2 <- list()
  pec_RLMM2_noOutliers <- list()
  pec_RLMM2_residuals <- list()
  pec_RLMM2_noOutliers_residuals <- list()
  pec_RLMM_compare <- list()
  for (i in 1:nVars) {
    pec_RLMM2[[i]] <- update(pec_RLMM[[i]], rho.sigma.e = psi2propII(smoothPsi, k = 2.28),rho.sigma.b = psi2propII(smoothPsi, k = 2.28))
    pec_RLMM2_residuals[[i]] <- resid(pec_RLMM[[i]])
    pec_RLMM2_noOutliers[[i]] <- update(pec_RLMM_noOutliers[[i]], rho.sigma.e = psi2propII(smoothPsi, k = 2.28),rho.sigma.b = psi2propII(smoothPsi, k = 2.28))
    pec_RLMM2_noOutliers_residuals[[i]] <- resid(pec_RLMM2_noOutliers[[i]])
    pec_RLMM_compare[[i]] <- compare(pec_LMM[[i]], pec_RLMM[[i]], pec_RLMM_noOutliers[[i]],  pec_RLMM2[[i]],pec_RLMM2_noOutliers[[i]], show.rho.functions = FALSE)
  }
  names(pec_RLMM2) <- modelFormulae
  names(pec_RLMM2_noOutliers) <- modelFormulae
  names(pec_RLMM2_residuals) <- variablesToAnalyze[1:7]
  names(pec_RLMM2_noOutliers_residuals) <- variablesToAnalyze[1:7]
  names(pec_RLMM_compare) <- modelFormulae
  
  
  savePlots_RLMM <- function(rlmm, fileName, width = 10, height = 6, units = "in", res = 300, nGraphs = 3, ...) {
    saveName <- list("_FittedVsResiduals.jpg",
                     "_NormalQQVsResiduals.jpg",
                     "_NormalQQVsREs.jpg", 
                     "_REsForGroup.jpg"
    )
    for (i in 1:nGraphs) {
     robustlmm::plot(rlmm, which = i)
     ggsave(paste(fileName, saveName[[i]], sep = ""))
   }
  }

  
  for (i in 1:7) { 
    savePlots_RLMM(pec_RLMM[[i]], paste("pec_RLMM_", variablesToAnalyze[i], sep = ""))
  }


  
  
  # For data with outliers
  pec_RLMM_shapiro <- list()
  pec_RLMM2_shapiro <- list()
  for (i in 1:nVars) {
    pec_RLMM_shapiro[[i]] <- shapiro.test(resid(pec_RLMM[[i]]))
    pec_RLMM2_shapiro[[i]] <- shapiro.test(resid(pec_RLMM2[[i]]))
  }

  names(pec_RLMM_shapiro) <- modelFormulae
  # Only InterpAP_BW met the assumption of normality based on the Shapiro-Wilk tests
  
  names(pec_RLMM2_shapiro) <- modelFormulae
  
  
  
  # For data without outliers
  pec_RLMM_noOutliers_shapiro <- list()
  pec_RLMM2_noOutliers_shapiro <- list()
  for (i in 1:nVars) {
    pec_RLMM_noOutliers_shapiro[[i]] <- shapiro.test(resid(pec_RLMM_noOutliers[[i]]))
    pec_RLMM2_noOutliers_shapiro[[i]] <- shapiro.test(resid(pec_RLMM2_noOutliers[[i]]))
  }
  names(pec_RLMM_noOutliers_shapiro) <- modelFormulae
  names(pec_RLMM2_noOutliers_shapiro) <- modelFormulae
  # removing the outliers made it so more variables met the assumption of normal residuals
  # only V_BW and NetGRF_BW were normal
  
  
  

  #### NLME with random intercepts: peak net GRF - Pectoral ####
  ## based on: https://stats.stackexchange.com/questions/77891/checking-assumptions-lmer-lme-mixed-models-in-r
  
  lm.base2 <- lme(PercentStance ~ group, random= ~1|individual, method="ML", data= pec_peakNetGRF_noOutliers$usableData)
  plot(lm.base2)
  qqnorm(resid(lm.base2))
  
  
  ### a) Run lmes
  
  ## with outliers and without outliers
  # make sure the iterations are high enough to allow the models to converge
  lmeControl(msMaxIter = 50)
  
  modelFormulae2 <- list()
  pec_LMM2 <- list()
  pec_LMM2_noOutliers <- list()
  pec_LMM2_residuals <- list()
  pec_LMM2_noOutliers_residuals <- list()
  for (i in 1:nVars) {
    modelFormulae2[[i]] <- as.formula(paste(variablesToAnalyze[i], "~group", sep = ""))
    pec_LMM2[[i]] <- lme(modelFormulae2[[i]], random = ~1|individual, weights = varIdent(form = ~ 1 | group), data = pec_peakNetGRFs)
    pec_LMM2_residuals[[i]] <- residuals(pec_LMM2[[i]])
    pec_LMM2_noOutliers[[i]] <- lme(modelFormulae2[[i]], random = ~1|individual, weights = varIdent(form = ~ 1 | group), data = data.frame(pec_peakNetGRF_noOutliers$usableData))
    pec_LMM2_noOutliers_residuals[[i]] <- residuals(pec_LMM2_noOutliers[[i]])
  }
  names(pec_LMM2) <- modelFormulae2
  names(pec_LMM2_noOutliers) <- modelFormulae2
  names(pec_LMM2_residuals) <- variablesToAnalyze[1:7]
  names(pec_LMM2_noOutliers_residuals) <- variablesToAnalyze[1:7]
  
  pec_LMM2_resids <- data.frame(do.call("cbind", pec_LMM2_residuals))
  pec_LMM2_noOutliers_resids <- data.frame(do.call("cbind", pec_LMM2_noOutliers_residuals))
  
  
  
  #### Testing the assumptions ####
  ## Don't need to test for linearity of data because the predictors are categorical
  
  ## b) evaluating the normality of the residuals
  # the null of the Shapiro-Wilk test is that the input (e.g., residuals of data) are normal
  
  # For data with outliers
  pec_LMM2_shapiro <- list()
  for (i in 1:nVars) {
    pec_LMM2_shapiro[[i]] <- shapiro.test(resid(pec_LMM2[[i]]))
    qqPlot(resid(pec_LMM2[[i]]), ylab = paste(names(pec_LMM2)[[i]], " residuals"))
  }
 
  names(pec_LMM2_shapiro) <- modelFormulae2
  # Only InterpAP_BW met the assumption of normality based on the Shapiro-Wilk tests
  
  # For data without outliers
  pec_LMM2_noOutliers_shapiro <- list()
  for (i in 1:nVars) {
    pec_LMM2_noOutliers_shapiro[[i]] <- shapiro.test(resid(pec_LMM2_noOutliers[[i]]))
    qqPlot(resid(pec_LMM2_noOutliers[[i]]), ylab = paste(names(pec_LMM2_noOutliers)[[i]], " residuals"))
  }
  names(pec_LMM2_noOutliers_shapiro) <- modelFormulae2
  # removing the outliers made it so more variables met the assumption of normal residuals
  # however, still not as many variables as the lmer models
  
  
  ## c) Testing homogeneity of variances
  # the Bartlett's test is more sensitive to non-normal data so people often use Levene's
  # more info here; http://www.sthda.com/english/wiki/compare-multiple-sample-variances-in-r
  # The Fligner-Killeen test can be used for non-normal data  
  
  # with full data set (same as lmer)
  apply(pec_peakNetGRFs[,variablesToAnalyze[1:7]],2,function(x) {leveneTest(x ~ as.factor(pec_peakNetGRFs$group))})
  apply(pec_peakNetGRFs[,variablesToAnalyze[1:7]],2,function(x) {fligner.test(x ~ as.factor(pec_peakNetGRFs$group))})
  
  # can also evaluate the homogeneity of variances graphically
  pp2 <- list()
  for (i in 1:nVars) {
    pp2[[i]] <- plot(pec_LMM2[[i]])
  }
  names(pp2) <- names(pec_LMM2)
  # there does not appear to be any observable pattern in the residuals vs. fitted plots, which suggests
  # that the variances are homogeneous
  
  # Testing assumption when the outliers were removed (same as lmer)
  apply(pec_peakNetGRF_noOutliers$usableData[,variablesToAnalyze[1:7]],2,function(x) {leveneTest(x ~ as.factor(pec_peakNetGRF_noOutliers$usableData$group))})
  apply(pec_peakNetGRF_noOutliers$usableData[,variablesToAnalyze[1:7]],2,function(x) {fligner.test(x ~ as.factor(pec_peakNetGRF_noOutliers$usableData$group))})
  
  # can also evaluate the homogeneity of variances graphically
  pp2_noOutliers <- list()
  for (i in 1:nVars) {
    pp2_noOutliers[[i]] <- plot(pec_LMM2_noOutliers[[i]])
  }
  names(pp2_noOutliers) <- names(pec_LMM2_noOutliers)
  
  
  ## d) Testing the normality of the random effects
  # following suggestions from: https://stats.stackexchange.com/questions/117170/testing-whether-random-effects-are-normally-distributed-in-r
  
  # ## random intercepts model
  # r_int2<- ranef(pec_LMM2[[1]])$individual$`(Intercept)`
  # qqnorm(r_int2)
  # qqline(r_int2)
  # shapiro.test(r_int2)
  # 
  # ## checking for equal variances across the random effects
  # jpeg("pec_peakNetGRF_PercentStance_REs.jpg", width = 10, height = 6, units = "in", res = 300)
  # plot( lm.base2, resid(., type = "p") ~ fitted(.) | individual,
  #       id = 0.05, adj = -0.3 )
  # dev.off()
  
  
  
  
  #### FIGURES - STAT ASSUMPTIONS ####
  
  #### Testing assumptions of LMER (full vs. no outlier) ####
  ## (export as 500 width x 600 height)
  
  ## PercentStance assumptions
  pec_LMM_PercentStance_QQ <- ggplot(data = pec_LMM_resids, mapping = aes(sample = PercentStance)) +
    stat_qq_band() +
    stat_qq_line() +
    stat_qq_point() +
    labs(x = "Theoretical Quantiles", y = "Residual Quantiles") + 
    theme_classic()
  
  pec_LMM_PercentStance_Fitted <- plot(pec_LMM[[1]])
  
  pec_LMM_noOutliers_PercentStance_QQ <- ggplot(data = pec_LMM_noOutliers_resids, mapping = aes(sample = PercentStance)) +
    stat_qq_band() +
    stat_qq_line() +
    stat_qq_point() +
    labs(x = "Theoretical Quantiles", y = "Residual Quantiles") + 
    theme_classic()
  
  pec_LMM_noOutliers_PercentStance_Fitted <- plot(pec_LMM_noOutliers[[1]])
  
  jpeg("pec_peakNetGRF_percentStance_stats.jpg", width = 10, height = 6, units = "in", res = 300)
  cowplot::plot_grid(pec_LMM_PercentStance_QQ, pec_LMM_PercentStance_Fitted, pec_LMM_noOutliers_PercentStance_QQ, pec_LMM_noOutliers_PercentStance_Fitted, labels = c("a", "b", "c", "d"))
  dev.off()
  
  ## Vertical GRF assumptions
  pec_LMM_VBW_QQ <- ggplot(data = pec_LMM_resids, mapping = aes(sample = InterpV_BW)) +
    stat_qq_band() +
    stat_qq_line() +
    stat_qq_point() +
    labs(x = "Theoretical Quantiles", y = "Residual Quantiles") + 
    theme_classic()
  
  pec_LMM_VBW_Fitted <- plot(pec_LMM[[2]])
  
  pec_LMM_noOutliers_VBW_QQ <- ggplot(data = pec_LMM_noOutliers_resids, mapping = aes(sample = InterpV_BW)) +
    stat_qq_band() +
    stat_qq_line() +
    stat_qq_point() +
    labs(x = "Theoretical Quantiles", y = "Residual Quantiles") + 
    theme_classic()
  
  pec_LMM_noOutliers_VBW_Fitted <- plot(pec_LMM_noOutliers[[2]])
  
  jpeg("pec_peakNetGRF_VBW_stats.jpg", width = 10, height = 6, units = "in", res = 300)
  cowplot::plot_grid(pec_LMM_VBW_QQ, pec_LMM_VBW_Fitted, pec_LMM_noOutliers_VBW_QQ, pec_LMM_noOutliers_VBW_Fitted, labels = c("a", "b", "c", "d"))
  dev.off()
  
  ## Mediolateral GRF assumptions
  pec_LMM_MLBW_QQ <- ggplot(data = pec_LMM_resids, mapping = aes(sample = InterpML_BW)) +
    stat_qq_band() +
    stat_qq_line() +
    stat_qq_point() +
    labs(x = "Theoretical Quantiles", y = "Residual Quantiles") + 
    theme_classic()
  
  pec_LMM_MLBW_Fitted <- plot(pec_LMM[[3]])
  
  pec_LMM_noOutliers_MLBW_QQ <- ggplot(data = pec_LMM_noOutliers_resids, mapping = aes(sample = InterpML_BW)) +
    stat_qq_band() +
    stat_qq_line() +
    stat_qq_point() +
    labs(x = "Theoretical Quantiles", y = "Residual Quantiles") + 
    theme_classic()
  
  pec_LMM_noOutliers_MLBW_Fitted <- plot(pec_LMM_noOutliers[[3]])
  
  jpeg("pec_peakNetGRF_MLBW_stats.jpg", width = 10, height = 6, units = "in", res = 300)
  cowplot::plot_grid(pec_LMM_MLBW_QQ, pec_LMM_MLBW_Fitted, pec_LMM_noOutliers_MLBW_QQ, pec_LMM_noOutliers_MLBW_Fitted, labels = c("a", "b", "c", "d"))
  dev.off()
  
  ## Anteroposterior GRF assumptions
  pec_LMM_APBW_QQ <- ggplot(data = pec_LMM_resids, mapping = aes(sample = InterpAP_BW)) +
    stat_qq_band() +
    stat_qq_line() +
    stat_qq_point() +
    labs(x = "Theoretical Quantiles", y = "Residual Quantiles") + 
    theme_classic()
  
  pec_LMM_APBW_Fitted <- plot(pec_LMM[[4]])
  
  pec_LMM_noOutliers_APBW_QQ <- ggplot(data = pec_LMM_noOutliers_resids, mapping = aes(sample = InterpAP_BW)) +
    stat_qq_band() +
    stat_qq_line() +
    stat_qq_point() +
    labs(x = "Theoretical Quantiles", y = "Residual Quantiles") + 
    theme_classic()
  
  pec_LMM_noOutliers_APBW_Fitted <- plot(pec_LMM_noOutliers[[4]])
  
  jpeg("pec_peakNetGRF_APBW_stats.jpg", width = 10, height = 6, units = "in", res = 300)
  cowplot::plot_grid(pec_LMM_APBW_QQ, pec_LMM_APBW_Fitted, pec_LMM_noOutliers_APBW_QQ, pec_LMM_noOutliers_APBW_Fitted, labels = c("a", "b", "c", "d"))
  dev.off()
  
  ## Net GRF assumptions
  pec_LMM_NetGRFBW_QQ <- ggplot(data = pec_LMM_resids, mapping = aes(sample = NetGRF_BW)) +
    stat_qq_band() +
    stat_qq_line() +
    stat_qq_point() +
    labs(x = "Theoretical Quantiles", y = "Residual Quantiles") + 
    theme_classic()
  
  pec_LMM_NetGRFBW_Fitted <- plot(pec_LMM[[5]])
  
  pec_LMM_noOutliers_NetGRFBW_QQ <- ggplot(data = pec_LMM_noOutliers_resids, mapping = aes(sample = NetGRF_BW)) +
    stat_qq_band() +
    stat_qq_line() +
    stat_qq_point() +
    labs(x = "Theoretical Quantiles", y = "Residual Quantiles") + 
    theme_classic()
  
  pec_LMM_noOutliers_NetGRFBW_Fitted <- plot(pec_LMM_noOutliers[[5]])
  
  jpeg("pec_peakNetGRF_NetGRFBW_stats.jpg", width = 10, height = 6, units = "in", res = 300)
  cowplot::plot_grid(pec_LMM_NetGRFBW_QQ, pec_LMM_NetGRFBW_Fitted, pec_LMM_noOutliers_NetGRFBW_QQ, pec_LMM_noOutliers_NetGRFBW_Fitted, labels = c("a", "b", "c", "d"))
  dev.off()
  
  ## ML Angle assumptions
  pec_LMM_MLAngleConvert_QQ <- ggplot(data = pec_LMM_resids, mapping = aes(sample = MLAngle_Convert_deg)) +
    stat_qq_band() +
    stat_qq_line() +
    stat_qq_point() +
    labs(x = "Theoretical Quantiles", y = "Residual Quantiles") + 
    theme_classic()
  
  pec_LMM_MLAngleConvert_Fitted <- plot(pec_LMM[[6]])
  
  pec_LMM_noOutliers_MLAngleConvert_QQ <- ggplot(data = pec_LMM_noOutliers_resids, mapping = aes(sample = MLAngle_Convert_deg)) +
    stat_qq_band() +
    stat_qq_line() +
    stat_qq_point() +
    labs(x = "Theoretical Quantiles", y = "Residual Quantiles") + 
    theme_classic()
  
  pec_LMM_noOutliers_MLAngleConvert_Fitted <- plot(pec_LMM_noOutliers[[6]])
  
  jpeg("pec_peakNetGRF_MLAngleConvert_stats.jpg", width = 10, height = 6, units = "in", res = 300)
  cowplot::plot_grid(pec_LMM_MLAngleConvert_QQ, pec_LMM_MLAngleConvert_Fitted, pec_LMM_noOutliers_MLAngleConvert_QQ, pec_LMM_noOutliers_MLAngleConvert_Fitted, labels = c("a", "b", "c", "d"))
  dev.off()
  
  ## AP Angle assumptions
  pec_LMM_APAngleConvert_QQ <- ggplot(data = pec_LMM_resids, mapping = aes(sample = APAngle_Convert_deg)) +
    stat_qq_band() +
    stat_qq_line() +
    stat_qq_point() +
    labs(x = "Theoretical Quantiles", y = "Residual Quantiles") + 
    theme_classic()
  
  pec_LMM_APAngleConvert_Fitted <- plot(pec_LMM[[7]])
  
  pec_LMM_noOutliers_APAngleConvert_QQ <- ggplot(data = pec_LMM_noOutliers_resids, mapping = aes(sample = APAngle_Convert_deg)) +
    stat_qq_band() +
    stat_qq_line() +
    stat_qq_point() +
    labs(x = "Theoretical Quantiles", y = "Residual Quantiles") + 
    theme_classic()
  
  pec_LMM_noOutliers_APAngleConvert_Fitted <- plot(pec_LMM_noOutliers[[7]])
  
  jpeg("pec_peakNetGRF_APAngleConvert_stats.jpg", width = 10, height = 6, units = "in", res = 300)
  cowplot::plot_grid(pec_LMM_APAngleConvert_QQ, pec_LMM_APAngleConvert_Fitted, pec_LMM_noOutliers_APAngleConvert_QQ, pec_LMM_noOutliers_APAngleConvert_Fitted, labels = c("a", "b", "c", "d"))
  dev.off()
  
  
  #### Testing assumptions of LMM (no outliers) ####
  jpeg("pec_peakNetGRF_noOutlier_stats1.jpg", width = 8.5, height = 11, units = "in", res = 300)
  cowplot::plot_grid(
    pec_LMM_noOutliers_PercentStance_QQ, pec_LMM_noOutliers_PercentStance_Fitted,
    pec_LMM_noOutliers_VBW_QQ, pec_LMM_noOutliers_VBW_Fitted,
    pec_LMM_noOutliers_MLBW_QQ, pec_LMM_noOutliers_MLBW_Fitted,
    pec_LMM_noOutliers_APBW_QQ, pec_LMM_noOutliers_APBW_Fitted,
    ncol = 2,
    labels = c("a", "b", "c", "d", "e", "f", "g", "h")
    )
  dev.off()
  
  jpeg("pec_peakNetGRF_noOutlier_stats2.jpg", width = 8.5, height = 11, units = "in", res = 300)
  cowplot::plot_grid(
  pec_LMM_noOutliers_NetGRFBW_QQ, pec_LMM_noOutliers_NetGRFBW_Fitted,
  pec_LMM_noOutliers_MLAngleConvert_QQ, pec_LMM_noOutliers_MLAngleConvert_Fitted,
  pec_LMM_noOutliers_APAngleConvert_QQ, pec_LMM_noOutliers_APAngleConvert_Fitted,
  ncol = 2,
  labels = c("i", "j", "k", "l", "m", "n")
  )
  dev.off()
  
  
  
  #### Testing assumptions of LME (full vs. no outlier) ####
  ## (export as 500 width x 600 height)
  
  ## PercentStance assumptions
  pec_LMM2_PercentStance_QQ <- ggplot(data = pec_LMM2_resids, mapping = aes(sample = PercentStance)) +
    stat_qq_band() +
    stat_qq_line() +
    stat_qq_point() +
    labs(x = "Theoretical Quantiles", y = "Residual Quantiles") + 
    theme_classic()
  
  pec_LMM2_PercentStance_Fitted <- plot(pec_LMM2[[1]])
  
  pec_LMM2_noOutliers_PercentStance_QQ <- ggplot(data = pec_LMM2_noOutliers_resids, mapping = aes(sample = PercentStance)) +
    stat_qq_band() +
    stat_qq_line() +
    stat_qq_point() +
    labs(x = "Theoretical Quantiles", y = "Residual Quantiles") + 
    theme_classic()
  
  pec_LMM2_noOutliers_PercentStance_Fitted <- plot(pec_LMM2_noOutliers[[1]])
  
  jpeg("pec_peakNetGRF_percentStance_stats2.jpg", width = 10, height = 6, units = "in", res = 300)
  cowplot::plot_grid(pec_LMM2_PercentStance_QQ, pec_LMM2_PercentStance_Fitted, pec_LMM2_noOutliers_PercentStance_QQ, pec_LMM2_noOutliers_PercentStance_Fitted, labels = c("a", "b", "c", "d"))
  dev.off()
  
  ## Vertical GRF assumptions
  pec_LMM2_VBW_QQ <- ggplot(data = pec_LMM2_resids, mapping = aes(sample = InterpV_BW)) +
    stat_qq_band() +
    stat_qq_line() +
    stat_qq_point() +
    labs(x = "Theoretical Quantiles", y = "Residual Quantiles") + 
    theme_classic()
  
  pec_LMM2_VBW_Fitted <- plot(pec_LMM2[[2]])
  
  pec_LMM2_noOutliers_VBW_QQ <- ggplot(data = pec_LMM2_noOutliers_resids, mapping = aes(sample = InterpV_BW)) +
    stat_qq_band() +
    stat_qq_line() +
    stat_qq_point() +
    labs(x = "Theoretical Quantiles", y = "Residual Quantiles") + 
    theme_classic()
  
  pec_LMM2_noOutliers_VBW_Fitted <- plot(pec_LMM2_noOutliers[[2]])
  
  jpeg("pec_peakNetGRF_VBW_stats2.jpg", width = 10, height = 6, units = "in", res = 300)
  cowplot::plot_grid(pec_LMM2_VBW_QQ, pec_LMM2_VBW_Fitted, pec_LMM2_noOutliers_VBW_QQ, pec_LMM2_noOutliers_VBW_Fitted, labels = c("a", "b", "c", "d"))
  dev.off()
  
  ## Mediolateral GRF assumptions
  pec_LMM2_MLBW_QQ <- ggplot(data = pec_LMM2_resids, mapping = aes(sample = InterpML_BW)) +
    stat_qq_band() +
    stat_qq_line() +
    stat_qq_point() +
    labs(x = "Theoretical Quantiles", y = "Residual Quantiles") + 
    theme_classic()
  
  pec_LMM2_MLBW_Fitted <- plot(pec_LMM2[[3]])
  
  pec_LMM2_noOutliers_MLBW_QQ <- ggplot(data = pec_LMM2_noOutliers_resids, mapping = aes(sample = InterpML_BW)) +
    stat_qq_band() +
    stat_qq_line() +
    stat_qq_point() +
    labs(x = "Theoretical Quantiles", y = "Residual Quantiles") + 
    theme_classic()
  
  pec_LMM2_noOutliers_MLBW_Fitted <- plot(pec_LMM2_noOutliers[[3]])
  
  jpeg("pec_peakNetGRF_MLBW_stats2.jpg", width = 10, height = 6, units = "in", res = 300)
  cowplot::plot_grid(pec_LMM2_MLBW_QQ, pec_LMM2_MLBW_Fitted, pec_LMM2_noOutliers_MLBW_QQ, pec_LMM2_noOutliers_MLBW_Fitted, labels = c("a", "b", "c", "d"))
  dev.off()
  
  ## Anteroposterior GRF assumptions
  pec_LMM2_APBW_QQ <- ggplot(data = pec_LMM2_resids, mapping = aes(sample = InterpAP_BW)) +
    stat_qq_band() +
    stat_qq_line() +
    stat_qq_point() +
    labs(x = "Theoretical Quantiles", y = "Residual Quantiles") + 
    theme_classic()
  
  pec_LMM2_APBW_Fitted <- plot(pec_LMM2[[4]])
  
  pec_LMM2_noOutliers_APBW_QQ <- ggplot(data = pec_LMM2_noOutliers_resids, mapping = aes(sample = InterpAP_BW)) +
    stat_qq_band() +
    stat_qq_line() +
    stat_qq_point() +
    labs(x = "Theoretical Quantiles", y = "Residual Quantiles") + 
    theme_classic()
  
  pec_LMM2_noOutliers_APBW_Fitted <- plot(pec_LMM2_noOutliers[[4]])
  
  jpeg("pec_peakNetGRF_APBW_stats2.jpg", width = 10, height = 6, units = "in", res = 300)
  cowplot::plot_grid(pec_LMM2_APBW_QQ, pec_LMM2_APBW_Fitted, pec_LMM2_noOutliers_APBW_QQ, pec_LMM2_noOutliers_APBW_Fitted, labels = c("a", "b", "c", "d"))
  dev.off()
  
  ## Net GRF assumptions
  pec_LMM2_NetGRFBW_QQ <- ggplot(data = pec_LMM2_resids, mapping = aes(sample = NetGRF_BW)) +
    stat_qq_band() +
    stat_qq_line() +
    stat_qq_point() +
    labs(x = "Theoretical Quantiles", y = "Residual Quantiles") + 
    theme_classic()
  
  pec_LMM2_NetGRFBW_Fitted <- plot(pec_LMM2[[5]])
  
  pec_LMM2_noOutliers_NetGRFBW_QQ <- ggplot(data = pec_LMM2_noOutliers_resids, mapping = aes(sample = NetGRF_BW)) +
    stat_qq_band() +
    stat_qq_line() +
    stat_qq_point() +
    labs(x = "Theoretical Quantiles", y = "Residual Quantiles") + 
    theme_classic()
  
  pec_LMM2_noOutliers_NetGRFBW_Fitted <- plot(pec_LMM2_noOutliers[[5]])
  
  jpeg("pec_peakNetGRF_NetGRFBW_stats2.jpg", width = 10, height = 6, units = "in", res = 300)
  cowplot::plot_grid(pec_LMM2_NetGRFBW_QQ, pec_LMM2_NetGRFBW_Fitted, pec_LMM2_noOutliers_NetGRFBW_QQ, pec_LMM2_noOutliers_NetGRFBW_Fitted, labels = c("a", "b", "c", "d"))
  dev.off()
  
  ## ML Angle assumptions
  pec_LMM2_MLAngleConvert_QQ <- ggplot(data = pec_LMM2_resids, mapping = aes(sample = MLAngle_Convert_deg)) +
    stat_qq_band() +
    stat_qq_line() +
    stat_qq_point() +
    labs(x = "Theoretical Quantiles", y = "Residual Quantiles") + 
    theme_classic()
  
  pec_LMM2_MLAngleConvert_Fitted <- plot(pec_LMM2[[6]])
  
  pec_LMM2_noOutliers_MLAngleConvert_QQ <- ggplot(data = pec_LMM2_noOutliers_resids, mapping = aes(sample = MLAngle_Convert_deg)) +
    stat_qq_band() +
    stat_qq_line() +
    stat_qq_point() +
    labs(x = "Theoretical Quantiles", y = "Residual Quantiles") + 
    theme_classic()
  
  pec_LMM2_noOutliers_MLAngleConvert_Fitted <- plot(pec_LMM2_noOutliers[[6]])
  
  jpeg("pec_peakNetGRF_MLAngleConvert_stats2.jpg", width = 10, height = 6, units = "in", res = 300)
  cowplot::plot_grid(pec_LMM2_MLAngleConvert_QQ, pec_LMM2_MLAngleConvert_Fitted, pec_LMM2_noOutliers_MLAngleConvert_QQ, pec_LMM2_noOutliers_MLAngleConvert_Fitted, labels = c("a", "b", "c", "d"))
  dev.off()
  
  ## AP Angle assumptions
  pec_LMM2_APAngleConvert_QQ <- ggplot(data = pec_LMM2_resids, mapping = aes(sample = APAngle_Convert_deg)) +
    stat_qq_band() +
    stat_qq_line() +
    stat_qq_point() +
    labs(x = "Theoretical Quantiles", y = "Residual Quantiles") + 
    theme_classic()
  
  pec_LMM2_APAngleConvert_Fitted <- plot(pec_LMM2[[7]])
  
  pec_LMM2_noOutliers_APAngleConvert_QQ <- ggplot(data = pec_LMM2_noOutliers_resids, mapping = aes(sample = APAngle_Convert_deg)) +
    stat_qq_band() +
    stat_qq_line() +
    stat_qq_point() +
    labs(x = "Theoretical Quantiles", y = "Residual Quantiles") + 
    theme_classic()
  
  pec_LMM2_noOutliers_APAngleConvert_Fitted <- plot(pec_LMM2_noOutliers[[7]])
  
  jpeg("pec_peakNetGRF_APAngleConvert_stats3.jpg", width = 10, height = 6, units = "in", res = 300)
  cowplot::plot_grid(pec_LMM2_APAngleConvert_QQ, pec_LMM2_APAngleConvert_Fitted, pec_LMM2_noOutliers_APAngleConvert_QQ, pec_LMM2_noOutliers_APAngleConvert_Fitted, labels = c("a", "b", "c", "d"))
  dev.off()
  
  
  #### Testing assumptions of LMM (no outliers) ####
  jpeg("pec_peakNetGRF_noOutlier_stats2.jpg", width = 8.5, height = 11, units = "in", res = 300)
  cowplot::plot_grid(
    pec_LMM2_noOutliers_PercentStance_QQ, pec_LMM2_noOutliers_PercentStance_Fitted,
    pec_LMM2_noOutliers_VBW_QQ, pec_LMM2_noOutliers_VBW_Fitted,
    pec_LMM2_noOutliers_MLBW_QQ, pec_LMM2_noOutliers_MLBW_Fitted,
    pec_LMM2_noOutliers_APBW_QQ, pec_LMM2_noOutliers_APBW_Fitted,
    ncol = 2,
    labels = c("a", "b", "c", "d", "e", "f", "g", "h")
  )
  dev.off()
  
  jpeg("pec_peakNetGRF_noOutlier_stats4.jpg", width = 8.5, height = 11, units = "in", res = 300)
  cowplot::plot_grid(
    pec_LMM2_noOutliers_NetGRFBW_QQ, pec_LMM2_noOutliers_NetGRFBW_Fitted,
    pec_LMM2_noOutliers_MLAngleConvert_QQ, pec_LMM2_noOutliers_MLAngleConvert_Fitted,
    pec_LMM2_noOutliers_APAngleConvert_QQ, pec_LMM2_noOutliers_APAngleConvert_Fitted,
    ncol = 2,
    labels = c("i", "j", "k", "l", "m", "n")
  )
  dev.off()
  
  
  
  #### FIGURES - DATA ####
  
  ## Create function to produce box plots with jitter points and marginal densities on the sides
  boxWithDensityPlot <- function(df, xName, yName, xLabel, yLabel) {
    # create the plot
    original_plot <- df %>% 
      ggplot(aes_string(x = xName, y = yName)) + 
      geom_boxplot(aes_string(color = xName), show.legend = FALSE) +
      geom_jitter(position=position_jitter(0.2), alpha = 0.5, aes_string(color = xName), show.legend = FALSE) +
      xlab(paste("\n", xLabel)) +
      ylab(paste(yLabel, "\n")) +
      theme_pubr() + 
      border()      
    
    y_density <- axis_canvas(original_plot, axis = "y", coord_flip = TRUE) +
      geom_density(data = df, aes_string(x = yName, fill = xName), color = NA, alpha = 0.5) +
      coord_flip()
    
    # create the combined plot
    #combined_plot %<>% insert_yaxis_grob(., y_density, position = "right")
    combined_plot <- insert_yaxis_grob(original_plot, y_density, position = "right")
    
    # plot the resulting combined plot
    ggdraw(combined_plot)
  }
  
  
  ### With all data 
  pec_peakNetGRFs$species <- substring(pec_peakNetGRFs$group, 1, 2)

  pec_peakNetGRFs_VBW_plot <- boxWithDensityPlot(pec_peakNetGRFs, "species", "InterpV_BW", "", "GRF - vertical")
  pec_peakNetGRFs_MLBW_plot <- boxWithDensityPlot(pec_peakNetGRFs, "species", "InterpML_BW", "", "GRF - mediolateral")
  pec_peakNetGRFs_APBW_plot <- boxWithDensityPlot(pec_peakNetGRFs, "species", "InterpAP_BW", "", "GRF - anteroposterior")
  pec_peakNetGRFs_netBW_plot <- boxWithDensityPlot(pec_peakNetGRFs, "species", "NetGRF_BW", "", "GRF - net")
  pec_peakNetGRFs_MLangle_plot <- boxWithDensityPlot(pec_peakNetGRFs, "species", "MLAngle_Convert_deg", "", "mediolateral angle")
  pec_peakNetGRFs_APangle_plot <- boxWithDensityPlot(pec_peakNetGRFs, "species", "APAngle_Convert_deg", "", "anteroposterior angle")
  
  jpeg("pec_peakNetGRF_plots.jpg", width = 8.5, height = 11, units = "in", res = 300)
  cowplot::plot_grid(pec_peakNetGRFs_VBW_plot, 
                     pec_peakNetGRFs_MLBW_plot, 
                     pec_peakNetGRFs_APBW_plot,
                     pec_peakNetGRFs_netBW_plot, 
                     pec_peakNetGRFs_MLangle_plot,
                     pec_peakNetGRFs_APangle_plot,
                     ncol = 3,
                     #labels = c("a", "b", "c", "d", "e", "f")
                     labels = "AUTO"
                     )
  dev.off()
  
  
  ### With Outliers removed
 
  #pec_peakNetGRFs$species <- substring(pec_peakNetGRFs$group, 1, 2)
  
  pec_peakNetGRFs_noOutliers_VBW_plot <- boxWithDensityPlot(pec_peakNetGRF_noOutliers$usableData, "species", "InterpV_BW", "", "GRF - vertical")
  pec_peakNetGRFs_noOutliers_MLBW_plot <- boxWithDensityPlot(pec_peakNetGRF_noOutliers$usableData, "species", "InterpML_BW", "", "GRF - mediolateral")
  pec_peakNetGRFs_noOutliers_APBW_plot <- boxWithDensityPlot(pec_peakNetGRF_noOutliers$usableData, "species", "InterpAP_BW", "", "GRF - anteroposterior")
  pec_peakNetGRFs_noOutliers_netBW_plot <- boxWithDensityPlot(pec_peakNetGRF_noOutliers$usableData, "species", "NetGRF_BW", "", "GRF - net")
  pec_peakNetGRFs_noOutliers_MLangle_plot <- boxWithDensityPlot(pec_peakNetGRF_noOutliers$usableData, "species", "MLAngle_Convert_deg", "", "mediolateral angle")
  pec_peakNetGRFs_noOutliers_APangle_plot <- boxWithDensityPlot(pec_peakNetGRF_noOutliers$usableData, "species", "APAngle_Convert_deg", "", "anteroposterior angle")
  
  jpeg("pec_peakNetGRF_noOutliers_plots.jpg", width = 8.5, height = 11, units = "in", res = 300)
  cowplot::plot_grid(pec_peakNetGRFs_noOutliers_VBW_plot, 
                     pec_peakNetGRFs_noOutliers_MLBW_plot, 
                     pec_peakNetGRFs_noOutliers_APBW_plot,
                     pec_peakNetGRFs_noOutliers_netBW_plot, 
                     pec_peakNetGRFs_noOutliers_MLangle_plot,
                     pec_peakNetGRFs_noOutliers_APangle_plot,
                     ncol = 3,
                     #labels = c("a", "b", "c", "d", "e", "f")
                     labels = "AUTO"
  )
  dev.off()
  
  

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
  
  
  ############## UNUSED: Discriminant Function Analysis   ##########
  ## There doesn't seem to be a DFA with mixed effects option: https://stats.stackexchange.com/questions/33372/discriminant-analysis-with-random-effects
  ## so the best option is to collapse data to individual means, so the sample size will be based
  ## on total number of individuals rather than trials. 
  ## However, I may not have a large enough sample size bc I'd have more traits than individuals per group.
  ## For example, max # of independent variables is n - 2 where n = sample size: http://userwww.sfsu.edu/efc/classes/biol710/discrim/discrim.pdf
  ## Given that, it doesn't make sense to run a discriminant function analysis on these data
  
  ## calculating Cohen's d for each fixed effect
  ## https://stackoverflow.com/questions/57566566/calculating-confidence-intervals-around-lme-dscore-cohens-d-for-mixed-effect-mo
  # EMAtools::lme.dscore(pec_LME, data = subset(peakNetGRF, appendage == "pectoral"), type="lme4")
  
  # ## calculating confidence intervals for each fixed effect
  # confint(pec_EMM, method="Wald")
  
  # ## pseudo-R2 for mixed effects models: https://ecologyforacrowdedplanet.wordpress.com/2013/08/27/r-squared-in-mixed-models-the-easy-way/
  # ## "Marginal R2GLMM represents the variance explained by the fixed effects
  # ## Conditional R2GLMM is interpreted as a variance explained by the entire model, including both
  # ## fixed and random effects."
  # r.squaredGLMM(pec_LMM)
  # 
  # ## From the MuMIN docs: "R2GLMM can be calculated also for fixed-effect models. In the simpliest case of OLS it reduces to
  # ## var(fitted) / (var(fitted) + deviance / 2). Unlike likelihood-ratio based R2 for OLS, value of this statistic differs from that of the classical R2."
  # 
  ## Testing assumption of normality with base plots
  # qqnorm(resid(pec_LMM[[i]]))
  # qqline(resid(pec_LMM[[i]]))
  

  # # UNUSED: PeakNetGRF - Pectoral: plot raw data ####
  # ## evaluate whether there are major outliers using boxplots and violin plots
  # 
  # ## need to change this so it uses the subsetted data for the pectoral trials
  # ## need to get rid of the redudant legends
  # # this shows how to get a common legend for combined plots: https://www.datanovia.com/en/lessons/combine-multiple-ggplots-into-a-figure/
  # p <- list()
  # for(i in 1:nVars){
  #   p[[i]] <- ggplot(pec_peakNetGRFs, aes_string(x = variablesToAnalyze[8], y = variablesToAnalyze[i], fill = variablesToAnalyze[8])) + 
  #     geom_violin(trim = FALSE) + 
  #     geom_boxplot(width= 0.1, fill = "white") +
  #     labs(x = "appendage", y = variablesToAnalyze[i]) + 
  #     theme(legend.position = "none") + 
  #     theme_classic()
  # }
  # do.call(grid.arrange, p)
  # 
  # histos <- list()
  # for(i in 1:nVars){
  #   histos[[i]] <- ggplot(pec_peakNetGRFs, aes_string(x = variablesToAnalyze[i], color = group, fill = group)) + 
  #     geom_histogram(aes(y=..density..), alpha = 0.2, colour="black") +
  #     geom_density(alpha = 0.5) +
  #     theme_classic()
  # }
  # 
  # ggplot(pec_peakNetGRFs, aes(x = PercentStance)) + 
  #   geom_histogram(aes(y=..density..), colour="black", fill="white") +
  #   geom_density(alpha=.2, fill="#FF6666") +
  #   theme_classic()


  ## When removing outliers for each variable separately and then figuring which are the common ('unique') outliers across the variables
  # pec_peakNetGRF_noOutliers <- removeOutliers(pec_peakNetGRFs, variablesToAnalyze[1:7], "group", "filename", outlierRange = 1.5)
  # pec_peakNetGRF_allOutliers <- do.call("rbind", pec_peakNetGRF_noOutliers$outliers)
  # unique(pec_peakNetGRF_allOutliers$filename)
  
  
  # # Checking if taking the log helps PercentStance meet normality
  # shapiro.test(resid(pec_PercentStance_log_LMM))
  # qqPlot(resid(pec_PercentStance_log_LMM), ylab = "log(PercentStance) residuals") # that may be worse than no log
  # 
  # shapiro.test(resid(pec_PercentStance_log_noOutliers_LMM))
  # qqPlot(resid(pec_PercentStance_log_noOutliers_LMM), ylab = "log(PercentStance) residuals") 
  # # no, it does not. Regardless of whether the outliers were removed
  # 
