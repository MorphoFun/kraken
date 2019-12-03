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


#### YANK ####



  ## Need to add duty factor and appendage frequency data
  
  ### Combine the peak net GRF data
  GRFs$Pelvic$Pel_GRFs_Filtered_PeakNet$appendage <- "pelvic"

  GRFs$Pectoral$Pec_GRFs_Filtered_PeakNet$appendage <- "pectoral"  

  
  peakNetGRF <- rbind(GRFs$Pelvic$Pel_GRFs_Filtered_PeakNet, GRFs$Pectoral$Pec_GRFs_Filtered_PeakNet)
  peakNetGRF$group <- paste(substring(peakNetGRF$filename, 1, 2), substring(peakNetGRF$appendage, 1, 3), sep = "_")
  peakNetGRF$individual <- substring(peakNetGRF$filename, 1, 4)
  
  #### PeakNetGRF: summary stats ####
  variablesToAnalyze <- (c("PercentStance", "InterpV_BW", "InterpML_BW", "InterpAP_BW", "NetGRF_BW", "MLAngle_Convert_deg", "APAngle_Convert_deg", "group", "individual", "appendage"))
  aggregate(. ~ group, data = peakNetGRF[,variablesToAnalyze[1:(length(variablesToAnalyze)-2)]], FUN = function(x) c(mean = mean(x), sd = sd(x), n = length(x)))
  nVars <- 7
  
  #### PeakNetGRF: plot raw data ####
  ## evaluate whether there are major outliers using boxplots and violin plots
  
  ## need to get rid of the redudant legends
  # this shows how to get a common legend for combined plots: https://www.datanovia.com/en/lessons/combine-multiple-ggplots-into-a-figure/
  p <- list()
  for(i in 1:nVars){
    p[[i]] <- ggplot(peakNetGRF, aes_string(x = variablesToAnalyze[8], y = variablesToAnalyze[i], fill = variablesToAnalyze[8])) + 
      geom_violin(trim = FALSE) + 
      geom_boxplot(width= 0.1, fill = "white") +
      labs(x = "appendage", y = "variable") + 
      theme(legend.position = "none") + 
      theme_classic()
  }
  do.call(grid.arrange, p)
  
  
  ## ggplot2 calculates more outliers bc baseplot::boxplot doesn't actually calculate the 1st and 3rd quantiles with even n
  # https://stackoverflow.com/questions/21793715/why-geom-boxplot-identify-more-outliers-than-base-boxplot
  boxplot_wOutliers <- function(df, x, y, group = NULL, label = NULL, ...) {
    # calculate boxplot object
    g <- ggplot(df, aes_string(x = x, y = y, fill = group)) + geom_boxplot() + guides(fill=FALSE)
    
    # get list of outliers 
    out <- ggplot_build(g)[["data"]][[1]][["outliers"]]
    
    # label list elements with factor levels
    names(out) <- levels(as.factor(group))
    
    # convert to tidy data
    tidyout <- purrr::map_df(out, tibble::as_tibble, .id = "category")
    tidyout$label <- df[,label][tidyout$value]
    tidyout$group <- levels(as.factor(df[,group]))
    
    # plot boxplots with labels
    g + geom_text(data = tidyout, aes(group, value, label = label),  hjust = -0.1)  +
      theme(legend.position = "none") +
      theme_classic()
  }
  
  boxplot_wOutliers(peakNetGRF, x = "group", y = "PercentStance", group = "group", label = "filename")
  
removeOutliers <- function(df, x, y, ...) {
  out <- boxplot(df[,y] ~ df[,x], col = rainbow(5), xlab = x, ylab = y)$out
  outliers <- df[out,]
  usableData <- df[-c(out),]
  output <- list(
    outliers = outliers,
    usableData = usableData
  )
  return(output)
}

peakNetGRF_noOutliers <- removeOutliers(peakNetGRF, "group", "PercentStance")

## can't seem to get rid of the trials that are labeled as outliers
# the outliers that are determined by boxplot don't seem to be matching up with the plot
boxplot_wOutliers(peakNetGRF_noOutliers$usableData, x = "group", y = "PercentStance", group = "group", label = "filename")


  boxplot_wOutliers(peakNetGRF_noOutliers, x = "group", y = "PercentStance", group = "group", label = "filename")
  
  
  
  
  #### LINEAR MIXED EFFECTS MODELS ####
  ## This will be used to compare the means between groups while accounting for the repeated trials within individuals
  ## Since there are no pelvic data for Periophthalmus, we'll have three models:
  ## 1) pectoral comparison between Periophthalmus, Ambystoma, and Pleurodeles
  ## 2) pelvic comparison between Amvbystoma and Pleurodeles
  ## 3) pectoral versus pelvic for Ambystoma and Pleurodeles
  ## The effect sizes of individual independent variables can be assessed through the fixed effects: 
  ## https://stat.ethz.ch/pipermail/r-sig-mixed-models/2013q4/021102.html
  ## Or, could consider the f2 value (Aiken and West 1991): https://largescaleassessmentsineducation.springeropen.com/articles/10.1186/s40536-018-0061-2

  
  #### LMM: Pec Comparison ####
  #pec_LMM <- lmer(. ~ group + (1 |individual), data = subset(peakNetGRF[,variablesToAnalyze], appendage == "pectoral"), REML = TRUE)
  #summary(pec_LMM)
  
  modelFormulae <- list()
  pec_LMM <- list()
  for (i in 1:nVars) {
    modelFormulae[[i]] <- as.formula(paste(variablesToAnalyze[i], "~group+(1|individual)", sep = ""))
    pec_LMM[[i]] <- lmer(modelFormulae[[i]], data = subset(peakNetGRF[,variablesToAnalyze], appendage == "pectoral"))

  }
  names(pec_LMM) <- modelFormulae
  
  
  ## trying to log-transform to help with normality but can't take the log of negative numbers
  modelFormulae_log <- list()
  pec_LMM_log <- list()
  for (i in 1:nVars) {
    modelFormulae_log[[i]] <- as.formula(paste("log(", variablesToAnalyze[i], ")~group+(1|individual)", sep = ""))
    pec_LMM_log[[i]] <- lmer(modelFormulae_log[[i]], data = subset(peakNetGRF[,variablesToAnalyze], appendage == "pectoral"))
    
  }
  names(pec_LMM_log) <- modelFormulae_log


  ### Testing the assumptions
  ## Looking at residuals vs. fitted values to evaluate whether there are any nonlinear relationships
  ## gunshot pattern suggests that there aren't any nonlinear patterns, so okay to use a LMM vs. GLMM
  plots <- list()
  for (i in 1:nVars) {
    plots[[i]] <- plot(pec_LMM[[i]])
  }
  do.call(grid.arrange, plots)
  
  ## evaluating the normality of the residuals
  qqplots <- list()
  pec_LMM_shapiro <- list()
  for (i in 1:nVars) {
    qqplots[[i]] <-   qqnorm(resid(pec_LMM[[i]]))
    #qqline(resid(pec_LMM[[i]]))
    pec_LMM_shapiro[[i]] <- shapiro.test(resid(pec_LMM[[i]]))
  }
  do.call(grid.arrange, qqplots)
  names(pec_LMM_shapiro) <- modelFormulae
  
  ## following up 
  
  ## testing out the car::qqPlot() function bc it provides confidence intervals around the qqline
  qqPlot(resid(pec_LMM[[7]]))
  

  
  
  qqnorm(resid(pec_LMM[[i]]))
  qqline(resid(pec_LMM[[i]]))
  
  
  ### Continuing on with statistical analyses since data seem to meet assumptions of LMM
  
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
  
