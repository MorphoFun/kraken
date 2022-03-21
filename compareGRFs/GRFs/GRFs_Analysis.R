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
library(reshape2) # for melt()
library(MuMIn) # for r.squaredGLMM()
library(emmeans) # for emmeans() and pairs(); contrast() could be an option, too
library(gridExtra) # for grid.arrange()
library(grid) # for textGrob()
library(car) # for qqPlot()
library(ggplot2) # for ggplot()
library(cowplot) # for plot_grid() and ggdraw()
library(performance) # r2_nakagawa()
library(qqplotr) # for stat_qq_band()
library(ggpubr) # for theme_pubr

# check all objects in kraken package
ls("package:kraken")

forcePath = "./force_rawFiles"
setwd(forcePath)

# #### LOAD THE CALIBRATION FILE ####
# CalibFile <- data.frame(read.csv("./FinLimbGRFs_Calibs.csv", header=TRUE))
# save(CalibFile, file = "./data/FinLimbGRFs_Calibs.rda")
# 
# #### LOAD THE VIDEO INFO FILE ####
# VideoFile <- data.frame(read.csv("./compareGRFs/FinLimbGRFs_VideoInfo.csv", header=TRUE))
# save(CalibFile, file = "./data/FinLimbGRFs_VideoInfo.rda")


#### LOADING THE DATA ####


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
  VideoInfo$ForceZeroStart <- 1
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
peakNetGRF$speciesCode <- substring(peakNetGRF$filename, 1, 2)

peakNetGRF$species <- gsub('pb', 'P. barbarus', peakNetGRF$species)
peakNetGRF$species <- gsub('pw', 'P. waltl', peakNetGRF$species)
peakNetGRF$species <- gsub('af', 'A. tigrinum', peakNetGRF$species)

## Subsetting peak net GRF data into groups to analyze
pec_peakNetGRFs <- subset(peakNetGRF, appendage == "pectoral")
pel_peakNetGRFs <- subset(peakNetGRF, appendage == "pelvic")
sal_peakNetGRFs <- subset(peakNetGRF, !substring(filename, 1, 2) == "pb")


#### GRF - PROFILE PLOTS ####

## Choosing a color palette that is friendly to color blindness
#brown, light blue, green
cbPalette <- c("#661100", "#56B4E9", "#009E73")

#brown, green
cbPalette_brgr <- c("#661100", "#009E73")


# Create common label for x-axis
x.grob <- textGrob("\n Percent Stance", 
                   gp=gpar(fontface="bold", col="black", fontsize=12))

# produce common legend
get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}


### Compiling pectoral data
pec_GRFs <- list(
  vertical = lapply(GRFs$Pectoral$Pec_GRFs_Filtered_dataset_noOverlap, function(x) data.frame(percentStance = seq(0,100), InterpV_BW = x$InterpV_BW)),
  mediolateral = lapply(GRFs$Pectoral$Pec_GRFs_Filtered_dataset_noOverlap, function(x) data.frame(percentStance = seq(0,100), InterpML_BW = x$InterpML_BW)),
  anteroposterior = lapply(GRFs$Pectoral$Pec_GRFs_Filtered_dataset_noOverlap, function(x) data.frame(percentStance = seq(0,100), InterpAP_BW = x$InterpAP_BW)),
  net = lapply(GRFs$Pectoral$Pec_GRFs_Filtered_dataset_noOverlap, function(x) data.frame(percentStance = seq(0,100), net = x$NetGRF_BW)),
  ml_angle =   lapply(GRFs$Pectoral$Pec_GRFs_Filtered_dataset_noOverlap, function(x) data.frame(percentStance = seq(0,100), MLAngle_Convert_deg = x$MLAngle_Convert_deg)),
  ap_angle = lapply(GRFs$Pectoral$Pec_GRFs_Filtered_dataset_noOverlap, function(x) data.frame(percentStance = seq(0,100), APAngle_Convert_deg = x$APAngle_Convert_deg))
)

pec_GRFs_combined <- list(
  vertical = melt(pec_GRFs$vertical, id.vars = "percentStance", value.name = "vertical_BW"),
  mediolateral = melt(pec_GRFs$mediolateral, id.vars = "percentStance", value.name = "mediolateral_BW"),
  anteroposterior = melt(pec_GRFs$anteroposterior, id.vars = "percentStance", value.name = "anteroposterior_BW"),
  net = melt(pec_GRFs$net, id.vars = "percentStance", value.name = "net_BW"), 
  ml_ang = melt(pec_GRFs$ml_angle, id.vars = "percentStance", value.name = "mlang_deg"),
  ap_ang = melt(pec_GRFs$ap_angle, id.vars = "percentStance", value.name = "apang_deg")
)

for (i in 1:length(pec_GRFs_combined)) {
  pec_GRFs_combined[[i]]$speciesCode = substring(pec_GRFs_combined[[i]][,4], 1, 2)
  pec_GRFs_combined[[i]]$species <- ifelse(pec_GRFs_combined[[i]]$speciesCode == "pb", "P. barbarus", ifelse(pec_GRFs_combined[[i]]$speciesCode == "af", "A. tigrinum", "P. waltl"))
}


## Pectoral - GRF plots (in units of BW per percent of stance)
pec_GRF_v <- profilePlotR(pec_GRFs_combined$vertical, "percentStance", "vertical_BW", "species", "Percent Stance", "GRF - vertical (BW)", colorPalette = cbPalette, yrange = c(0, 0.5))
pec_GRF_ml <- profilePlotR(pec_GRFs_combined$mediolateral, "percentStance", "mediolateral_BW", "species", "Percent Stance", "GRF - mediolateral (BW)", colorPalette = cbPalette, yrange = c(-0.2, 0.2))
pec_GRF_ap <- profilePlotR(pec_GRFs_combined$anteroposterior, "percentStance", "anteroposterior_BW", "species", "Percent Stance", "GRF - anteroposterior (BW)", colorPalette = cbPalette, yrange = c(-0.2, 0.2))
pec_GRF_net <- profilePlotR(pec_GRFs_combined$net, "percentStance", "net_BW", "species", "Percent Stance", "GRF - Net (BW)", colorPalette = cbPalette, yrange = c(0, 0.5))
pec_GRF_mlang <- profilePlotR(pec_GRFs_combined$ml_ang, "percentStance", "mlang_deg", "species", "Percent Stance", "GRF - mediolateral angle (deg)", colorPalette = cbPalette, yrange = c(-50, 50))
pec_GRF_apang <- profilePlotR(pec_GRFs_combined$ap_ang, "percentStance", "apang_deg", "species", "Percent Stance", "GRF - anteroposterior angle (deg)", colorPalette = cbPalette, yrange = c(-50, 50))



pec_GRF_prow <- cowplot::plot_grid(
  pec_GRF_net + theme(axis.title.x = element_blank(),
                      legend.position = "none" ),
  pec_GRF_v + theme(axis.title.x = element_blank(),
                      legend.position = "none"), 
  pec_GRF_ml + theme(axis.title.x = element_blank(),
                       legend.position = "none" ), 
  pec_GRF_ap + theme(axis.title.x = element_blank(),
                       legend.position = "none" ),
  pec_GRF_mlang + theme(axis.title.x = element_blank(),
                      legend.position = "none" ),
  pec_GRF_apang + theme(axis.title.x = element_blank(),
                      legend.position = "none" ),
  ncol = 2,
  labels = "auto")


pec_GRF_legend <- get_legend(pec_GRF_v)

# Produce plot with insets and common x-axis label
jpeg("pec_profile_plots.jpg", width = 8.5, height = 11, units = "in", res = 600)
  grid.arrange(arrangeGrob(pec_GRF_prow, bottom = x.grob), pec_GRF_legend, heights = c(1, .2))
dev.off()

pdf("pec_profile_plots.pdf", width = 8.5, height = 11)
  grid.arrange(arrangeGrob(pec_GRF_prow, bottom = x.grob), pec_GRF_legend, heights = c(1, .2), top = textGrob("Figure 1: Pectoral appendages \n",gp=gpar(fontsize=20)))
dev.off()


### Compiling pelvic data
pel_GRFs <- list(
  vertical = lapply(GRFs$Pelvic$Pel_GRFs_Filtered_dataset_noOverlap, function(x) data.frame(percentStance = seq(0,100), InterpV_BW = x$InterpV_BW)),
  mediolateral = lapply(GRFs$Pelvic$Pel_GRFs_Filtered_dataset_noOverlap, function(x) data.frame(percentStance = seq(0,100), InterpML_BW = x$InterpML_BW)),
  anteroposterior = lapply(GRFs$Pelvic$Pel_GRFs_Filtered_dataset_noOverlap, function(x) data.frame(percentStance = seq(0,100), InterpAP_BW = x$InterpAP_BW)),
  net = lapply(GRFs$Pelvic$Pel_GRFs_Filtered_dataset_noOverlap, function(x) data.frame(percentStance = seq(0,100), net = x$NetGRF_BW)),
  ml_angle =   lapply(GRFs$Pelvic$Pel_GRFs_Filtered_dataset_noOverlap, function(x) data.frame(percentStance = seq(0,100), MLAngle_Convert_deg = x$MLAngle_Convert_deg)),
  ap_angle = lapply(GRFs$Pelvic$Pel_GRFs_Filtered_dataset_noOverlap, function(x) data.frame(percentStance = seq(0,100), APAngle_Convert_deg = x$APAngle_Convert_deg))
)

pel_GRFs_combined <- list(
  vertical = melt(pel_GRFs$vertical, id.vars = "percentStance", value.name = "vertical_BW"),
  mediolateral = melt(pel_GRFs$mediolateral, id.vars = "percentStance", value.name = "mediolateral_BW"),
  anteroposterior = melt(pel_GRFs$anteroposterior, id.vars = "percentStance", value.name = "anteroposterior_BW"),
  net = melt(pel_GRFs$net, id.vars = "percentStance", value.name = "net_BW"), 
  ml_ang = melt(pel_GRFs$ml_angle, id.vars = "percentStance", value.name = "mlang_deg"),
  ap_ang = melt(pel_GRFs$ap_angle, id.vars = "percentStance", value.name = "apang_deg")
)

for (i in 1:length(pel_GRFs_combined)) {
  pel_GRFs_combined[[i]]$speciesCode = substring(pel_GRFs_combined[[i]][,4], 1, 2)
  pel_GRFs_combined[[i]]$species <- ifelse(pel_GRFs_combined[[i]]$speciesCode == "af", "A. tigrinum", "P. waltl")
}

## Pelvic - GRF plots (in units of BW per percent of stance)
pel_GRF_v <- profilePlotR(pel_GRFs_combined$vertical, "percentStance", "vertical_BW", "species", "Percent Stance", "GRF - vertical (BW)", cbPalette_brgr, yrange = c(0, 0.5))
pel_GRF_ml <- profilePlotR(pel_GRFs_combined$mediolateral, "percentStance", "mediolateral_BW", "species", "Percent Stance", "GRF - mediolateral (BW)", cbPalette_brgr, yrange = c(-0.2, 0.2))
pel_GRF_ap <- profilePlotR(pel_GRFs_combined$anteroposterior, "percentStance", "anteroposterior_BW", "species", "Percent Stance", "GRF - anteroposterior (BW)", cbPalette_brgr, yrange = c(-0.2, 0.2))
pel_GRF_net <- profilePlotR(pel_GRFs_combined$net, "percentStance", "net_BW", "species", "Percent Stance", "GRF - Net (BW)", cbPalette_brgr, yrange = c(0, 0.5))
pel_GRF_mlang <- profilePlotR(pel_GRFs_combined$ml_ang, "percentStance", "mlang_deg", "species", "Percent Stance", "GRF - mediolateral angle (deg)", cbPalette_brgr, yrange = c(-50, 50))
pel_GRF_apang <- profilePlotR(pel_GRFs_combined$ap_ang, "percentStance", "apang_deg", "species", "Percent Stance", "GRF - anteroposterior angle (deg)", cbPalette_brgr, yrange = c(-50, 50))



pel_GRF_prow <- cowplot::plot_grid(
  pel_GRF_net + theme(axis.title.x = element_blank(),
                      legend.position = "none" ),
  pel_GRF_v + theme(axis.title.x = element_blank(),
                    legend.position = "none"), 
  pel_GRF_ml + theme(axis.title.x = element_blank(),
                     legend.position = "none" ), 
  pel_GRF_ap + theme(axis.title.x = element_blank(),
                     legend.position = "none" ),
  pel_GRF_mlang + theme(axis.title.x = element_blank(),
                        legend.position = "none" ),
  pel_GRF_apang + theme(axis.title.x = element_blank(),
                        legend.position = "none" ),
  ncol = 2,
  labels = "auto")


pel_GRF_legend <- get_legend(pel_GRF_v)

# Produce plot with insets and common x-axis label
jpeg("pel_profile_plots.jpg", width = 8.5, height = 11, units = "in", res = 600)
grid.arrange(arrangeGrob(pel_GRF_prow, bottom = x.grob), pel_GRF_legend, heights = c(1, .2))
dev.off()

pdf("pel_profile_plots.pdf", width = 8.5, height = 11)
grid.arrange(arrangeGrob(pel_GRF_prow, bottom = x.grob), pel_GRF_legend, heights = c(1, .2), top = textGrob("Figure 3: Pelvic appendages \n",gp=gpar(fontsize=20)))
dev.off()


#### PeakNetGRF: summary stats ####
variablesToAnalyze <- (c("PercentStance", "InterpV_BW", "InterpML_BW", "InterpAP_BW", "NetGRF_BW", "MLAngle_Convert_deg", "APAngle_Convert_deg", "group", "individual", "appendage"))

## summarize the mean, sd, and n (sample size) for each variable
aggregate(. ~ group, data = peakNetGRF[,variablesToAnalyze[1:(length(variablesToAnalyze)-2)]], FUN = function(x) c(mean = mean(x), sd = sd(x), n = length(x)))

## identify number of dependent variables to analyze
# "group", "individual", "appendage" are the independent variables
nVars <- length(variablesToAnalyze) - 3


#### PEAK NET GRF - PLOTS ####

## Create function to produce box plots with jitter points and marginal densities on the sides
boxWithDensityPlot <- function(df, xName, yName, xLabel, yLabel, grouplevels = NULL, colorPalette = NULL, xTickAngle = 0, ...) {
  # create the plot
  if(is.null(grouplevels)) {grouplevels = unique(df[,xName])}
  original_plot <- df %>% 
    ggplot(aes_string(x = xName, y = yName)) + 
    geom_boxplot(aes_string(color = xName), show.legend = FALSE) +
    geom_jitter(position=position_jitter(0.2), alpha = 0.5, aes_string(color = xName), show.legend = FALSE) +
    scale_color_manual(name = xName, # changing legend title
                       labels = grouplevels, # Changing legend labels
                       values = colorPalette) +
    scale_fill_manual(name = xName,
                      labels = grouplevels,
                      values = colorPalette) +
    xlab(paste("\n", xLabel)) +
    ylab(paste(yLabel, "\n")) +
    theme_pubr(x.text.angle = xTickAngle) + 
    border()      
  
  y_density <- axis_canvas(original_plot, axis = "y", coord_flip = TRUE) +
    geom_density(data = df, aes_string(x = yName, fill = xName), color = NA, alpha = 0.5) +
    scale_fill_manual(name = xName,
                      labels = grouplevels,
                      values = colorPalette) +
    coord_flip()
  
  # create the combined plot
  #combined_plot %<>% insert_yaxis_grob(., y_density, position = "right")
  combined_plot <- insert_yaxis_grob(original_plot, y_density, position = "right")
  
  # plot the resulting combined plot
  ggdraw(combined_plot)
}

### REORDERING SPECIES ###
pec_peakNetGRFs$species <- factor(pec_peakNetGRFs$species , levels=c("A. tigrinum", "P. waltl", "P. barbarus"))

## Choosing a color palette that is friendly to color blindness
#brown, green,  light blue
cbPalette_brblgr <- c("#661100", "#009E73", "#56B4E9")

### PECTORAL - ALL DATA
#pec_peakNetGRFs$species <- substring(pec_peakNetGRFs$group, 1, 2)

pec_peakNetGRFs_VBW_plot <- boxWithDensityPlot(pec_peakNetGRFs, "species", "InterpV_BW", "", "GRF - vertical", colorPalette = cbPalette_brblgr, xTickAngle = 45)
pec_peakNetGRFs_MLBW_plot <- boxWithDensityPlot(pec_peakNetGRFs, "species", "InterpML_BW", "", "GRF - mediolateral", colorPalette = cbPalette_brblgr, xTickAngle = 45)
pec_peakNetGRFs_APBW_plot <- boxWithDensityPlot(pec_peakNetGRFs, "species", "InterpAP_BW", "", "GRF - anteroposterior", colorPalette = cbPalette_brblgr, xTickAngle = 45)
pec_peakNetGRFs_netBW_plot <- boxWithDensityPlot(pec_peakNetGRFs, "species", "NetGRF_BW", "", "GRF - net", colorPalette = cbPalette_brblgr, xTickAngle = 45)
pec_peakNetGRFs_MLangle_plot <- boxWithDensityPlot(pec_peakNetGRFs, "species", "MLAngle_Convert_deg", "", "mediolateral angle", colorPalette = cbPalette_brblgr, xTickAngle = 45)
pec_peakNetGRFs_APangle_plot <- boxWithDensityPlot(pec_peakNetGRFs, "species", "APAngle_Convert_deg", "", "anteroposterior angle", colorPalette = cbPalette_brblgr, xTickAngle = 45)

jpeg("pec_peakNetGRF_plots.jpg", width = 8.5, height = 11, units = "in", res = 300)
cowplot::plot_grid(pec_peakNetGRFs_netBW_plot,
                   pec_peakNetGRFs_VBW_plot, 
                   pec_peakNetGRFs_MLBW_plot, 
                   pec_peakNetGRFs_APBW_plot,
                   pec_peakNetGRFs_MLangle_plot,
                   pec_peakNetGRFs_APangle_plot,
                   ncol = 3,
                   #labels = c("a", "b", "c", "d", "e", "f")
                   labels = "AUTO"
)
dev.off()

Fig2_title <- ggdraw() + 
  draw_label(
    "Figure 2",
    fontface = 'bold',
    x = 0,
    hjust = 0
  ) 
pdf("pec_peakNetGRF_plots.pdf", width = 8.5, height = 11, title = "Figure 4")
cowplot::plot_grid(Fig2_title, plot_grid(pec_peakNetGRFs_netBW_plot,
                   pec_peakNetGRFs_VBW_plot, 
                   pec_peakNetGRFs_MLBW_plot, 
                   pec_peakNetGRFs_APBW_plot,
                   pec_peakNetGRFs_MLangle_plot,
                   pec_peakNetGRFs_APangle_plot,
                   ncol = 3,
                   #labels = c("a", "b", "c", "d", "e", "f")
                   labels = "AUTO"
), ncol = 1, rel_heights = c(0.1, 1))
dev.off()


### PELVIC - ALL DATA
#pel_peakNetGRFs$species <- substring(pel_peakNetGRFs$group, 1, 2)

pel_peakNetGRFs_VBW_plot <- boxWithDensityPlot(pel_peakNetGRFs, "species", "InterpV_BW", "", "GRF - vertical", colorPalette = cbPalette_brgr, xTickAngle = 45)
pel_peakNetGRFs_MLBW_plot <- boxWithDensityPlot(pel_peakNetGRFs, "species", "InterpML_BW", "", "GRF - mediolateral", colorPalette = cbPalette_brgr, xTickAngle = 45)
pel_peakNetGRFs_APBW_plot <- boxWithDensityPlot(pel_peakNetGRFs, "species", "InterpAP_BW", "", "GRF - anteroposterior", colorPalette = cbPalette_brgr, xTickAngle = 45)
pel_peakNetGRFs_netBW_plot <- boxWithDensityPlot(pel_peakNetGRFs, "species", "NetGRF_BW", "", "GRF - net", colorPalette = cbPalette_brgr, xTickAngle = 45)
pel_peakNetGRFs_MLangle_plot <- boxWithDensityPlot(pel_peakNetGRFs, "species", "MLAngle_Convert_deg", "", "mediolateral angle", colorPalette = cbPalette_brgr, xTickAngle = 45)
pel_peakNetGRFs_APangle_plot <- boxWithDensityPlot(pel_peakNetGRFs, "species", "APAngle_Convert_deg", "", "anteroposterior angle", colorPalette = cbPalette_brgr, xTickAngle = 45)

jpeg("pel_peakNetGRF_plots.jpg", width = 8.5, height = 11, units = "in", res = 300)
cowplot::plot_grid(pel_peakNetGRFs_netBW_plot, 
                   pel_peakNetGRFs_VBW_plot, 
                   pel_peakNetGRFs_MLBW_plot, 
                   pel_peakNetGRFs_APBW_plot,
                   pel_peakNetGRFs_MLangle_plot,
                   pel_peakNetGRFs_APangle_plot,
                   ncol = 3,
                   #labels = c("a", "b", "c", "d", "e", "f")
                   labels = "AUTO"
)
dev.off()


Fig4_title <- ggdraw() + 
  draw_label(
    "Figure 4",
    fontface = 'bold',
    x = 0,
    hjust = 0
  ) 
pdf("pel_peakNetGRF_plots.pdf", width = 8.5, height = 11)
cowplot::plot_grid(Fig4_title, plot_grid(pel_peakNetGRFs_netBW_plot, 
                   pel_peakNetGRFs_VBW_plot, 
                   pel_peakNetGRFs_MLBW_plot, 
                   pel_peakNetGRFs_APBW_plot,
                   pel_peakNetGRFs_MLangle_plot,
                   pel_peakNetGRFs_APangle_plot,
                   ncol = 3,
                   #labels = c("a", "b", "c", "d", "e", "f")
                   labels = "AUTO"
), ncol = 1, rel_heights = c(0.1, 1))
dev.off()


#### PEAK NET GRF - LMM ####

### PECTORAL
# vertical 
pec_peakNetGRF_v_lmm <- lmer(InterpV_BW ~ species + (1|individual), data = pec_peakNetGRFs)
pec_peakNetGRF_v_emm <- emmeans(pec_peakNetGRF_v_lmm, "species")
pairs(pec_peakNetGRF_v_emm)
pec_peakNetGRF_v_lmm_omega2 <- performance::r2_xu(pec_peakNetGRF_v_lmm) # 0.2239348
performance::r2_nakagawa(pec_peakNetGRF_v_lmm) # c = 0.214, m = 0.154

# mediolateral
pec_peakNetGRF_ml_lmm <- lmer(InterpML_BW ~ species + (1|individual), data = pec_peakNetGRFs)
pec_peakNetGRF_ml_emm <- emmeans(pec_peakNetGRF_ml_lmm, "species")
pairs(pec_peakNetGRF_ml_emm)
pec_peakNetGRF_ml_lmm_omega2 <- performance::r2_xu(pec_peakNetGRF_ml_lmm) # 0.5055679
performance::r2_nakagawa(pec_peakNetGRF_ml_lmm) # c = 0.491, m = 0.260

# anteroposterior
pec_peakNetGRF_ap_lmm <- lmer(InterpAP_BW ~ species + (1|individual), data = pec_peakNetGRFs)
pec_peakNetGRF_ap_emm <- emmeans(pec_peakNetGRF_ap_lmm, "species")
pairs(pec_peakNetGRF_ap_emm)
pec_peakNetGRF_ap_lmm_omega2 <- performance::r2_xu(pec_peakNetGRF_ap_lmm) # 0.670572
performance::r2_nakagawa(pec_peakNetGRF_ap_lmm) # c = 0.651, m = 0.515

# net
pec_peakNetGRF_net_lmm <- lmer(NetGRF_BW ~ species + (1|individual), data = pec_peakNetGRFs)
pec_peakNetGRF_net_emm <- emmeans(pec_peakNetGRF_net_lmm, "species")
pairs(pec_peakNetGRF_net_emm)
pec_peakNetGRF_net_lmm_omega2 <- performance::r2_xu(pec_peakNetGRF_net_lmm) # 0.2309856
performance::r2_nakagawa(pec_peakNetGRF_net_lmm) # c = 0.225, m = 0.110

# mediolateral angle
pec_peakNetGRF_mlang_lmm <- lmer(MLAngle_Convert_deg ~ species + (1|individual), data = pec_peakNetGRFs)
pec_peakNetGRF_mlang_emm <- emmeans(pec_peakNetGRF_mlang_lmm, "species")
pairs(pec_peakNetGRF_mlang_emm)
pec_peakNetGRF_mlang_lmm_omega2 <- performance::r2_xu(pec_peakNetGRF_mlang_lmm) # 0.427933
performance::r2_nakagawa(pec_peakNetGRF_mlang_lmm) # c = 0.421, m = 0.258

# anteroposterior angle
pec_peakNetGRF_apang_lmm <- lmer(APAngle_Convert_deg ~ species + (1|individual), data = pec_peakNetGRFs)
pec_peakNetGRF_apang_emm <- emmeans(pec_peakNetGRF_apang_lmm, "species")
pairs(pec_peakNetGRF_apang_emm)
pec_peakNetGRF_apang_lmm_omega2 <- performance::r2_xu(pec_peakNetGRF_apang_lmm) # 0.6657909
performance::r2_nakagawa(pec_peakNetGRF_apang_lmm) # c = 0.645, m = 0.543

# percent stance
pec_peakNetGRF_percentstance_lmm <- lmer(PercentStance ~ species + (1|individual), data = pec_peakNetGRFs)
pec_peakNetGRF_percentstance_emm <- emmeans(pec_peakNetGRF_percentstance_lmm, "species")
pairs(pec_peakNetGRF_percentstance_emm)
pec_peakNetGRF_percentstance_lmm_omega2 <- performance::r2_xu(pec_peakNetGRF_percentstance_lmm) # 0.4574002
performance::r2_nakagawa(pec_peakNetGRF_percentstance_lmm) # c = 0.462, m = 0.226


### PELVIC
# vertical 
pel_peakNetGRF_v_lmm <- lmer(InterpV_BW ~ species + (1|individual), data = pel_peakNetGRFs)
pel_peakNetGRF_v_emm <- emmeans(pel_peakNetGRF_v_lmm, "species")
pairs(pel_peakNetGRF_v_emm)
pel_peakNetGRF_v_lmm_omega2 <- performance::r2_xu(pel_peakNetGRF_v_lmm) # 0.750374
performance::r2_nakagawa(pel_peakNetGRF_v_lmm) # c = 0.734, m = 0.558

# mediolateral
pel_peakNetGRF_ml_lmm <- lmer(InterpML_BW ~ species + (1|individual), data = pel_peakNetGRFs)
pel_peakNetGRF_ml_emm <- emmeans(pel_peakNetGRF_ml_lmm, "species")
pairs(pel_peakNetGRF_ml_emm)
pel_peakNetGRF_ml_lmm_omega2 <- performance::r2_xu(pel_peakNetGRF_ml_lmm) # 0.3266032 
performance::r2_nakagawa(pel_peakNetGRF_ml_lmm) # c = 0.350, m = 0.056

# anteroposterior
pel_peakNetGRF_ap_lmm <- lmer(InterpAP_BW ~ species + (1|individual), data = pel_peakNetGRFs)
pel_peakNetGRF_ap_emm <- emmeans(pel_peakNetGRF_ap_lmm, "species")
pairs(pel_peakNetGRF_ap_emm)
pel_peakNetGRF_ap_lmm_omega2 <- performance::r2_xu(pel_peakNetGRF_ap_lmm) # 0.4262423
performance::r2_nakagawa(pel_peakNetGRF_ap_lmm) # c = 0.418, m = 0.102

# net
pel_peakNetGRF_net_lmm <- lmer(NetGRF_BW ~ species + (1|individual), data = pel_peakNetGRFs)
pel_peakNetGRF_net_emm <- emmeans(pel_peakNetGRF_net_lmm, "species")
pairs(pel_peakNetGRF_net_emm)
pel_peakNetGRF_net_lmm_omega2 <- performance::r2_xu(pel_peakNetGRF_net_lmm) # 0.7405002
performance::r2_nakagawa(pel_peakNetGRF_net_lmm) # c = 0.715, m = 0.585

# mediolateral angle
pel_peakNetGRF_mlang_lmm <- lmer(MLAngle_Convert_deg ~ species + (1|individual), data = pel_peakNetGRFs)
pel_peakNetGRF_mlang_emm <- emmeans(pel_peakNetGRF_mlang_lmm, "species")
pairs(pel_peakNetGRF_mlang_emm)
pel_peakNetGRF_mlang_lmm_omega2 <- performance::r2_xu(pel_peakNetGRF_mlang_lmm) # 0.4554476
performance::r2_nakagawa(pel_peakNetGRF_mlang_lmm) # c = 0.483, m = 0.215 

# anteroposterior angle
pel_peakNetGRF_apang_lmm <- lmer(APAngle_Convert_deg ~ species + (1|individual), data = pel_peakNetGRFs)
pel_peakNetGRF_apang_emm <- emmeans(pel_peakNetGRF_apang_lmm, "species")
pairs(pel_peakNetGRF_apang_emm)
pel_peakNetGRF_apang_lmm_omega2 <- performance::r2_xu(pel_peakNetGRF_apang_lmm) # 0.3052985
performance::r2_nakagawa(pel_peakNetGRF_apang_lmm) # c = 0.315, m = 0.001

# percent stance
pel_peakNetGRF_percentstance_lmm <- lmer(PercentStance ~ species + (1|individual), data = pel_peakNetGRFs)
pel_peakNetGRF_percentstance_emm <- emmeans(pel_peakNetGRF_percentstance_lmm, "species")
pairs(pel_peakNetGRF_percentstance_emm)
pel_peakNetGRF_percentstance_lmm_omega2 <- performance::r2_xu(pel_peakNetGRF_percentstance_lmm) # 0.1924656
performance::r2_nakagawa(pel_peakNetGRF_percentstance_lmm) # c = 0.183, m = 0.019



#### PEAK NET GRF - QQ PLOTS ####

plotQQ <- function(lmm, ...) {
  lmm_resids <- data.frame(resids = residuals(lmm))
  ggplot(data = lmm_resids, mapping = aes(sample = resids)) +
    qqplotr::stat_qq_band() +
    qqplotr::stat_qq_line() +
    qqplotr::stat_qq_point() +
    labs(x = "Theoretical Quantiles", y = "Residual Quantiles") + 
    theme_classic()
}

#### PEAK NET GRF - PLOT ASSUMPTIONS ####
## (export as 500 width x 600 height)

### PECTORAL 
# vertical
pec_peakNetGRF_v_lmm_QQ <- plotQQ(pec_peakNetGRF_v_lmm)
pec_peakNetGRF_v_lmm_Fitted <- plot(pec_peakNetGRF_v_lmm)

# mediolateral
pec_peakNetGRF_ml_lmm_QQ <- plotQQ(pec_peakNetGRF_ml_lmm)
pec_peakNetGRF_ml_lmm_Fitted <- plot(pec_peakNetGRF_ml_lmm)

# anteroposterior
pec_peakNetGRF_ap_lmm_QQ <- plotQQ(pec_peakNetGRF_ap_lmm)
pec_peakNetGRF_ap_lmm_Fitted <- plot(pec_peakNetGRF_ap_lmm)

# net
pec_peakNetGRF_net_lmm_QQ <- plotQQ(pec_peakNetGRF_net_lmm)
pec_peakNetGRF_net_lmm_Fitted <- plot(pec_peakNetGRF_net_lmm)

# mediolateral angle
pec_peakNetGRF_mlang_lmm_QQ <- plotQQ(pec_peakNetGRF_mlang_lmm)
pec_peakNetGRF_mlang_lmm_Fitted <- plot(pec_peakNetGRF_mlang_lmm)

# anteroposterior angle
pec_peakNetGRF_apang_lmm_QQ <- plotQQ(pec_peakNetGRF_apang_lmm)
pec_peakNetGRF_apang_lmm_Fitted <- plot(pec_peakNetGRF_apang_lmm)

jpeg("pec_peaknetGRF_QQ.jpg", width = 8, height = 10, units = "in", res = 300)
cowplot::plot_grid(pec_peakNetGRF_net_lmm_QQ, pec_peakNetGRF_v_lmm_QQ, pec_peakNetGRF_ml_lmm_QQ,
                   pec_peakNetGRF_ap_lmm_QQ, pec_peakNetGRF_mlang_lmm_QQ, pec_peakNetGRF_apang_lmm_QQ, ncol = 2, labels = letters[1:6])
dev.off()

### PELVIC  
# vertical
pel_peakNetGRF_v_lmm_QQ <- plotQQ(pel_peakNetGRF_v_lmm)
pel_peakNetGRF_v_lmm_Fitted <- plot(pel_peakNetGRF_v_lmm)

# mediolateral
pel_peakNetGRF_ml_lmm_QQ <- plotQQ(pel_peakNetGRF_ml_lmm)
pel_peakNetGRF_ml_lmm_Fitted <- plot(pel_peakNetGRF_ml_lmm)

# anteroposterior
pel_peakNetGRF_ap_lmm_QQ <- plotQQ(pel_peakNetGRF_ap_lmm)
pel_peakNetGRF_ap_lmm_Fitted <- plot(pel_peakNetGRF_ap_lmm)

# net
pel_peakNetGRF_net_lmm_QQ <- plotQQ(pel_peakNetGRF_net_lmm)
pel_peakNetGRF_net_lmm_Fitted <- plot(pel_peakNetGRF_net_lmm)

# mediolateral angle
pel_peakNetGRF_mlang_lmm_QQ <- plotQQ(pel_peakNetGRF_mlang_lmm)
pel_peakNetGRF_mlang_lmm_Fitted <- plot(pel_peakNetGRF_mlang_lmm)

# anteroposterior angle
pel_peakNetGRF_apang_lmm_QQ <- plotQQ(pel_peakNetGRF_apang_lmm)
pel_peakNetGRF_apang_lmm_Fitted <- plot(pel_peakNetGRF_apang_lmm)

jpeg("pel_peaknetGRF_QQ.jpg", width = 8, height = 10, units = "in", res = 300)
cowplot::plot_grid(pel_peakNetGRF_net_lmm_QQ, pel_peakNetGRF_v_lmm_QQ, pel_peakNetGRF_ml_lmm_QQ,
                   pel_peakNetGRF_ap_lmm_QQ, pel_peakNetGRF_mlang_lmm_QQ, pel_peakNetGRF_apang_lmm_QQ, ncol = 2, labels = letters[1:6])
dev.off()


#### PEAK NET GRF - REMOVING OUTLIERS ####
## outliers are set as those points falling outside 2x the interquantile range
## only doing this for the variables that seemed to deviate from normality

# pec - ml
pec_peakNetGRF_ml_bp <-  boxplot(InterpML_BW ~ species, data = pec_peakNetGRFs, range = 2)
pec_peakNetGRF_ml_noOutliers <- subset(pec_peakNetGRFs, !InterpML_BW %in% c(pec_peakNetGRF_ml_bp$out))
pec_peakNetGRF_ml_noOutliers_lmm <- lmer(InterpML_BW ~ species + (1|individual), data = pec_peakNetGRF_ml_noOutliers)
plotQQ(pec_peakNetGRF_ml_noOutliers_lmm)
shapiro.test(resid(pec_peakNetGRF_ml_noOutliers_lmm))
# the QQ plot is improved by removing the one outlier
pec_peakNetGRF_ml_noOutliers_emm <- emmeans(pec_peakNetGRF_ml_noOutliers_lmm, "species")

# pec - ap 
pec_peakNetGRF_ap_bp <-  boxplot(InterpAP_BW ~ species, data = pec_peakNetGRFs, range = 2)
pec_peakNetGRF_ap_noOutliers <- subset(pec_peakNetGRFs, !InterpAP_BW %in% c(pec_peakNetGRF_ap_bp$out))
pec_peakNetGRF_ap_noOutliers_lmm <- lmer(InterpAP_BW ~ species + (1|individual), data = pec_peakNetGRF_ap_noOutliers)
plotQQ(pec_peakNetGRF_ap_noOutliers_lmm)
shapiro.test(resid(pec_peakNetGRF_ap_noOutliers_lmm))
pec_peakNetGRF_ap_noOutliers_emm <- emmeans(pec_peakNetGRF_ap_noOutliers_lmm, "species")
# note: no outliers, so no change in the estimates

# pec - ml angle
pec_peakNetGRF_mlang_bp <-  boxplot(MLAngle_Convert_deg ~ species, data = pec_peakNetGRFs, range = 2)
pec_peakNetGRF_mlang_noOutliers <- subset(pec_peakNetGRFs, !MLAngle_Convert_deg %in% c(pec_peakNetGRF_mlang_bp$out))
pec_peakNetGRF_mlang_noOutliers_lmm <- lmer(MLAngle_Convert_deg ~ species + (1|individual), data = pec_peakNetGRF_mlang_noOutliers)
plotQQ(pec_peakNetGRF_mlang_noOutliers_lmm)
shapiro.test(resid(pec_peakNetGRF_mlang_noOutliers_lmm))
# the QQ plot is improved by removing the two outliers
pec_peakNetGRF_mlang_noOutliers_emm <- emmeans(pec_peakNetGRF_mlang_noOutliers_lmm, "species")

# pec - ap angle
pec_peakNetGRF_apang_bp <-  boxplot(APAngle_Convert_deg ~ species, data = pec_peakNetGRFs, range = 2)
pec_peakNetGRF_apang_noOutliers <- subset(pec_peakNetGRFs, !APAngle_Convert_deg %in% c(pec_peakNetGRF_apang_bp$out))
pec_peakNetGRF_apang_noOutliers_lmm <- lmer(APAngle_Convert_deg ~ species + (1|individual), data = pec_peakNetGRF_apang_noOutliers)
plotQQ(pec_peakNetGRF_apang_noOutliers_lmm)
shapiro.test(resid(pec_peakNetGRF_apang_noOutliers_lmm))
# the QQ plot is improved by removing the one outlier
pec_peakNetGRF_apang_noOutliers_emm <- emmeans(pec_peakNetGRF_apang_noOutliers_lmm, "species")


# pel - ap angle
pel_peakNetGRF_apang_bp <-  boxplot(APAngle_Convert_deg ~ species, data = pel_peakNetGRFs, range = 2)
pel_peakNetGRF_apang_noOutliers <- subset(pel_peakNetGRFs, !APAngle_Convert_deg %in% c(pel_peakNetGRF_apang_bp$out))
pel_peakNetGRF_apang_noOutliers_lmm <- lmer(APAngle_Convert_deg ~ species + (1|individual), data = pel_peakNetGRF_apang_noOutliers)
plotQQ(pel_peakNetGRF_apang_noOutliers_lmm)
shapiro.test(resid(pel_peakNetGRF_apang_noOutliers_lmm))
# the QQ plot is improved by removing the one outlier
pel_peakNetGRF_apang_noOutliers_emm <- emmeans(pel_peakNetGRF_apang_noOutliers_lmm, "species")




#### YANK - CALCULATIONS ####

pec_yank <- list(
  vertical = lapply(GRFs$Pectoral$Pec_GRFs_Filtered_dataset_noOverlap, function(x) data.frame(percentStance = seq(0,100), yank = yank(x$PercentStance, x$InterpV_BW)[,4])),
  mediolateral = lapply(GRFs$Pectoral$Pec_GRFs_Filtered_dataset_noOverlap, function(x) data.frame(percentStance = seq(0,100), yank = yank(x$PercentStance, x$InterpML_BW)[,4])),
  anteroposterior = lapply(GRFs$Pectoral$Pec_GRFs_Filtered_dataset_noOverlap, function(x) data.frame(percentStance = seq(0,100), yank = yank(x$PercentStance, x$InterpAP_BW)[,4])),
  net = lapply(GRFs$Pectoral$Pec_GRFs_Filtered_dataset_noOverlap, function(x) data.frame(percentStance = seq(0,100), yank = yank(x$PercentStance, x$NetGRF_BW)[,4]))
)

pel_yank <- list(
  vertical = lapply(GRFs$Pelvic$Pel_GRFs_Filtered_dataset_noOverlap, function(x) data.frame(percentStance = seq(0,100), yank = yank(x$PercentStance, x$InterpV_BW)[,4])),
  mediolateral = lapply(GRFs$Pelvic$Pel_GRFs_Filtered_dataset_noOverlap, function(x) data.frame(percentStance = seq(0,100), yank = yank(x$PercentStance, x$InterpML_BW)[,4])),
  anteroposterior = lapply(GRFs$Pelvic$Pel_GRFs_Filtered_dataset_noOverlap, function(x) data.frame(percentStance = seq(0,100), yank = yank(x$PercentStance, x$InterpAP_BW)[,4])),
  net = lapply(GRFs$Pelvic$Pel_GRFs_Filtered_dataset_noOverlap, function(x) data.frame(percentStance = seq(0,100), yank = yank(x$PercentStance, x$NetGRF_BW)[,4]))
)
  

pec_yank_combined <- list(
  vertical = melt(pec_yank$vertical, id.vars = "percentStance", value.name = "yank"),
  mediolateral = melt(pec_yank$mediolateral, id.vars = "percentStance", value.name = "yank"),
  anteroposterior = melt(pec_yank$anteroposterior, id.vars = "percentStance", value.name = "yank"),
  net = melt(pec_yank$net, id.vars = "percentStance", value.name = "yank")
)

for (i in 1:length(pec_yank_combined)) {
  pec_yank_combined[[i]]$speciesCode = substring(pec_yank_combined[[i]][,4], 1, 2)
  pec_yank_combined[[i]]$species <- ifelse(pec_yank_combined[[i]]$speciesCode == "pb", "P. barbarus", ifelse(pec_yank_combined[[i]]$speciesCode == "af", "A. tigrinum", "P. waltl"))
}

pel_yank_combined <- list(
  vertical = melt(pel_yank$vertical, id.vars = "percentStance", value.name = "yank"),
  mediolateral = melt(pel_yank$mediolateral, id.vars = "percentStance", value.name = "yank"),
  anteroposterior = melt(pel_yank$anteroposterior, id.vars = "percentStance", value.name = "yank"),
  net = melt(pel_yank$net, id.vars = "percentStance", value.name = "yank")
)

for (i in 1:length(pel_yank_combined)) {
  pel_yank_combined[[i]]$speciesCode = substring(pel_yank_combined[[i]][,4], 1, 2)
  pel_yank_combined[[i]]$species <- ifelse(pel_yank_combined[[i]]$speciesCode == "pb", "P. barbarus", ifelse(pel_yank_combined[[i]]$speciesCode == "af", "A. tigrinum", "P. waltl"))
}

#### YANK - PLOTS ####


## Pectoral - yank plots (in units of BW per percent of stance)
pec_yank_pv <- profilePlotR(pec_yank_combined$vertical, "percentStance", "yank", "species", "Percent Stance", "Yank - vertical", colorPalette = cbPalette, yrange = c(-0.02, 0.02))
pec_yank_pml <- profilePlotR(pec_yank_combined$mediolateral, "percentStance", "yank", "species", "Percent Stance", "Yank - mediolateral", colorPalette = cbPalette, yrange = c(-0.02, 0.02))
pec_yank_pap <- profilePlotR(pec_yank_combined$anteroposterior, "percentStance", "yank", "species", "Percent Stance", "Yank - anteroposterior", colorPalette = cbPalette, yrange = c(-0.02, 0.02))
pec_yank_pnet <- profilePlotR(pec_yank_combined$net, "percentStance", "yank", "species", "Percent Stance", "Yank - Net", colorPalette = cbPalette, yrange = c(-0.02, 0.02))

pec_yank_prow <- cowplot::plot_grid(
                  pec_yank_pnet + theme(axis.title.x = element_blank(),
                                      legend.position = "none" ),                 
                  pec_yank_pv + theme(axis.title.x = element_blank(),
                                       legend.position = "none"), 
                   pec_yank_pml + theme(axis.title.x = element_blank(),
                                       legend.position = "none" ), 
                   pec_yank_pap + theme(axis.title.x = element_blank(),
                                       legend.position = "none" ),
                   nrow = 2,
                   labels = "auto")


pec_yank_legend <- get_legend(pec_yank_pv)

# Produce plot with insets and common x-axis label
jpeg("compareGRFs_pec_yank.jpg", width = 8, height = 10, units = "in", res = 300)
grid.arrange(arrangeGrob(pec_yank_prow, bottom = x.grob), pec_yank_legend, heights = c(1, .2), top = textGrob("Pectoral appendages \n",gp=gpar(fontsize=20)))
dev.off()

pdf("compareGRFs_pec_yank.pdf", width = 8, height = 10)
grid.arrange(arrangeGrob(pec_yank_prow, bottom = x.grob), pec_yank_legend, heights = c(1, .2), top = textGrob("Figure 5: Pectoral appendages \n",gp=gpar(fontsize=20)))
dev.off()

## Pelvic - yank plots
pel_yank_pv <- profilePlotR(pel_yank_combined$vertical, "percentStance", "yank", "species", "Percent Stance", "Yank - vertical", colorPalette = cbPalette_brgr, yrange = c(-0.02, 0.02))
pel_yank_pml <- profilePlotR(pel_yank_combined$mediolateral, "percentStance", "yank", "species", "Percent Stance", "Yank - mediolateral", colorPalette = cbPalette_brgr, yrange = c(-0.02, 0.02))
pel_yank_pap <- profilePlotR(pel_yank_combined$anteroposterior, "percentStance", "yank", "species", "Percent Stance", "Yank - anteroposterior", colorPalette = cbPalette_brgr, yrange = c(-0.02, 0.02))
pel_yank_pnet <- profilePlotR(pel_yank_combined$net, "percentStance", "yank", "species", "Percent Stance", "Yank - Net", colorPalette = cbPalette_brgr, yrange = c(-0.02, 0.02))

pel_yank_prow <- cowplot::plot_grid(
  pel_yank_pnet + theme(axis.title.x = element_blank(),
                        legend.position = "none" ),
  pel_yank_pv + theme(axis.title.x = element_blank(),
                      legend.position = "none"), 
  pel_yank_pml + theme(axis.title.x = element_blank(),
                       legend.position = "none" ), 
  pel_yank_pap + theme(axis.title.x = element_blank(),
                       legend.position = "none" ),
  nrow = 2,
  labels = "auto")


pel_yank_legend <- get_legend(pel_yank_pv)

# Produce plot with insets and common x-axis label
jpeg("compareGRFs_pel_yank.jpg", width = 8, height = 10, units = "in", res = 300)
grid.arrange(arrangeGrob(pel_yank_prow, bottom = x.grob), pel_yank_legend, heights = c(1, .2), top = textGrob("Pelvic appendages \n",gp=gpar(fontsize=20)))
dev.off()

pdf("compareGRFs_pel_yank.pdf", width = 8, height = 10)
grid.arrange(arrangeGrob(pel_yank_prow, bottom = x.grob), pel_yank_legend, heights = c(1, .2), top = textGrob("Figure 6: Pelvic appendages \n",gp=gpar(fontsize=20)))
dev.off()



#### YANK - SUMMARY (pooled data) #####
## remove the scientific notation
options(scipen=999)


### PECTORAL 
## calculating values of the mean and standard deviation for each percent of stance between the species
# vertical
pec_yank_sumStats_v <- aggregate(yank~percentStance*species, data = pec_yank_combined$vertical, function(x) c(mean = mean(x), sd = sd(x)))

# mediolateral
pec_yank_sumStats_ml <- aggregate(yank~percentStance*species, data = pec_yank_combined$mediolateral, function(x) c(mean = mean(x), sd = sd(x)))

# anteroposterior
pec_yank_sumStats_ap <- aggregate(yank~percentStance*species, data = pec_yank_combined$anteroposterior, function(x) c(mean = mean(x), sd = sd(x)))

# net
pec_yank_sumStats_net <- aggregate(yank~percentStance*species, data = pec_yank_combined$net, function(x) c(mean = mean(x), sd = sd(x)))



## calculating maximum yank from the average values
# don't need to substract 1 from which.max output for some reason

# vertical
aggregate(yank[,1]~species, data = pec_yank_sumStats_v, function(x) max(x, na.rm = TRUE))
# which percent of stance this occurs 
aggregate(yank[,1]~species, data = pec_yank_sumStats_v, function(x) which.max(x))


# mediolateral
aggregate(yank[,1]~species, data = pec_yank_sumStats_ml, function(x) max(x, na.rm = TRUE))
# which percent of stance this occurs 
aggregate(yank[,1]~species, data = pec_yank_sumStats_ml, function(x) which.max(x))

# anteroposterior
aggregate(yank[,1]~species, data = pec_yank_sumStats_ap, function(x) max(x, na.rm = TRUE))
# which percent of stance this occurs 
aggregate(yank[,1]~species, data = pec_yank_sumStats_ap, function(x) which.max(x))


# net
aggregate(yank[,1]~species, data = pec_yank_sumStats_net, function(x) max(x, na.rm = TRUE))
# which percent of stance this occurs 
aggregate(yank[,1]~species, data = pec_yank_sumStats_net, function(x) which.max(x))


## calculating minimum yank from the average values
# don't need to substract 1 from which.max output for some reason

# vertical
aggregate(yank[,1]~species, data = pec_yank_sumStats_v, function(x) min(x, na.rm = TRUE))
# which percent of stance this occurs 
aggregate(yank[,1]~species, data = pec_yank_sumStats_v, function(x) which.min(x))

# mediolateral
aggregate(yank[,1]~species, data = pec_yank_sumStats_ml, function(x) min(x, na.rm = TRUE))
# which percent of stance this occurs 
aggregate(yank[,1]~species, data = pec_yank_sumStats_ml, function(x) which.min(x))

# anteroposterior
aggregate(yank[,1]~species, data = pec_yank_sumStats_ap, function(x) min(x, na.rm = TRUE))
# which percent of stance this occurs 
aggregate(yank[,1]~species, data = pec_yank_sumStats_ap, function(x) which.min(x))

# net
aggregate(yank[,1]~species, data = pec_yank_sumStats_net, function(x) min(x, na.rm = TRUE))
# which percent of stance this occurs 
aggregate(yank[,1]~species, data = pec_yank_sumStats_ap, function(x) which.min(x))



### PELVIC
## calculating values of the mean and standard deviation for each percent of stance between the species
# vertical
pel_yank_sumStats_v <- aggregate(yank~percentStance*species, data = pel_yank_combined$vertical, function(x) c(mean = mean(x), sd = sd(x)))

# mediolateral
pel_yank_sumStats_ml <- aggregate(yank~percentStance*species, data = pel_yank_combined$mediolateral, function(x) c(mean = mean(x), sd = sd(x)))

# anteroposterior
pel_yank_sumStats_ap <- aggregate(yank~percentStance*species, data = pel_yank_combined$anteroposterior, function(x) c(mean = mean(x), sd = sd(x)))

# net
pel_yank_sumStats_net <- aggregate(yank~percentStance*species, data = pel_yank_combined$net, function(x) c(mean = mean(x), sd = sd(x)))



## calculating maximum yank from the average values
# don't need to substract 1 from which.max output for some reason

# vertical
aggregate(yank[,1]~species, data = pel_yank_sumStats_v, function(x) max(x, na.rm = TRUE))
# which percent of stance this occurs 
aggregate(yank[,1]~species, data = pel_yank_sumStats_v, function(x) which.max(x))

# mediolateral
aggregate(yank[,1]~species, data = pel_yank_sumStats_ml, function(x) max(x, na.rm = TRUE))
# which percent of stance this occurs 
aggregate(yank[,1]~species, data = pel_yank_sumStats_ml, function(x) which.max(x))

# anteroposterior
aggregate(yank[,1]~species, data = pel_yank_sumStats_ap, function(x) max(x, na.rm = TRUE))
# which percent of stance this occurs 
aggregate(yank[,1]~species, data = pel_yank_sumStats_ap, function(x) which.max(x))

# net
aggregate(yank[,1]~species, data = pel_yank_sumStats_net, function(x) max(x, na.rm = TRUE))
# which percent of stance this occurs 
aggregate(yank[,1]~species, data = pel_yank_sumStats_net, function(x) which.max(x))


## calculating minimum yank from the average values
# don't need to substract 1 from which.max output for some reason

# vertical
aggregate(yank[,1]~species, data = pel_yank_sumStats_v, function(x) min(x, na.rm = TRUE))
# which percent of stance this occurs 
aggregate(yank[,1]~species, data = pel_yank_sumStats_v, function(x) which.min(x))

# mediolateral
aggregate(yank[,1]~species, data = pel_yank_sumStats_ml, function(x) min(x, na.rm = TRUE))
# which percent of stance this occurs 
aggregate(yank[,1]~species, data = pel_yank_sumStats_ml, function(x) which.min(x))

# anteroposterior
aggregate(yank[,1]~species, data = pel_yank_sumStats_ap, function(x) min(x, na.rm = TRUE))
# which percent of stance this occurs 
aggregate(yank[,1]~species, data = pel_yank_sumStats_ap, function(x) which.min(x))

# net
aggregate(yank[,1]~species, data = pel_yank_sumStats_net, function(x) min(x, na.rm = TRUE))
# which percent of stance this occurs 
aggregate(yank[,1]~species, data = pel_yank_sumStats_net, function(x) which.min(x))


#### YANK - SUMMARY (unpooled data) ####
## this calculates the maximum and minimum yank values within each trial, rather than from the averaged data

### PECTORAL
## substrating one from which.max because the first observation is 0% of stance

## Maximum
# vertical 
pec_yank_v_max <- data.frame(yank_max = unlist(cbind(lapply(pec_yank$vertical, FUN = function(x) max(x[,2], na.rm = TRUE)))))
pec_yank_v_max$species <- substring(names(pec_yank$vertical), 1, 2)
pec_yank_v_max$individual <- substring(names(pec_yank$vertical), 1, 4)
aggregate(yank_max~species, data = pec_yank_v_max, function(x) c(mean = mean(x), sd = sd(x)))

pec_yank_v_maxwhich <- data.frame(yank_maxwhich = unlist(cbind(lapply(pec_yank$vertical, FUN = function(x) which.max(x[,2])-1))))
pec_yank_v_maxwhich$species <- substring(names(pec_yank$vertical), 1, 2)
aggregate(yank_maxwhich~species, data = pec_yank_v_maxwhich, function(x) c(mean = mean(x), sd = sd(x)))

# mediolateral 
pec_yank_ml_max <- data.frame(yank_max = unlist(cbind(lapply(pec_yank$mediolateral, FUN = function(x) max(x[,2], na.rm = TRUE)))))
pec_yank_ml_max$species <- substring(names(pec_yank$mediolateral), 1, 2)
pec_yank_ml_max$individual <- substring(names(pec_yank$mediolateral), 1, 4)
aggregate(yank_max~species, data = pec_yank_ml_max, function(x) c(mean = mean(x), sd = sd(x)))

pec_yank_ml_maxwhich <- data.frame(yank_maxwhich = unlist(cbind(lapply(pec_yank$mediolateral, FUN = function(x) which.max(x[,2])-1))))
pec_yank_ml_maxwhich$species <- substring(names(pec_yank$mediolateral), 1, 2)
aggregate(yank_maxwhich~species, data = pec_yank_ml_maxwhich, function(x) c(mean = mean(x), sd = sd(x)))

# anteroposterior
pec_yank_ap_max <- data.frame(yank_max = unlist(cbind(lapply(pec_yank$anteroposterior, FUN = function(x) max(x[,2], na.rm = TRUE)))))
pec_yank_ap_max$species <- substring(names(pec_yank$anteroposterior), 1, 2)
pec_yank_ap_max$individual <- substring(names(pec_yank$anteroposterior), 1, 4)
aggregate(yank_max~species, data = pec_yank_ap_max, function(x) c(mean = mean(x), sd = sd(x)))

pec_yank_ap_maxwhich <- data.frame(yank_maxwhich = unlist(cbind(lapply(pec_yank$anteroposterior, FUN = function(x) which.max(x[,2])-1))))
pec_yank_ap_maxwhich$species <- substring(names(pec_yank$anteroposterior), 1, 2)
aggregate(yank_maxwhich~species, data = pec_yank_ap_maxwhich, function(x) c(mean = mean(x), sd = sd(x)))

# net
pec_yank_net_max <- data.frame(yank_max = unlist(cbind(lapply(pec_yank$net, FUN = function(x) max(x[,2], na.rm = TRUE)))))
pec_yank_net_max$species <- substring(names(pec_yank$net), 1, 2)
pec_yank_net_max$individual <- substring(names(pec_yank$net), 1, 4)
aggregate(yank_max~species, data = pec_yank_net_max, function(x) c(mean = mean(x), sd = sd(x)))

pec_yank_net_maxwhich <- data.frame(yank_maxwhich = unlist(cbind(lapply(pec_yank$net, FUN = function(x) which.max(x[,2])-1))))
pec_yank_net_maxwhich$species <- substring(names(pec_yank$net), 1, 2)
aggregate(yank_maxwhich~species, data = pec_yank_net_maxwhich, function(x) c(mean = mean(x), sd = sd(x)))


## Minimum
# vertical 
pec_yank_v_min <- data.frame(yank_min = unlist(cbind(lapply(pec_yank$vertical, FUN = function(x) min(x[,2], na.rm = TRUE)))))
pec_yank_v_min$species <- substring(names(pec_yank$vertical), 1, 2)
pec_yank_v_min$individual <- substring(names(pec_yank$vertical), 1, 4)
aggregate(yank_min~species, data = pec_yank_v_min, function(x) c(mean = mean(x), sd = sd(x)))

pec_yank_v_minwhich <- data.frame(yank_minwhich = unlist(cbind(lapply(pec_yank$vertical, FUN = function(x) which.min(x[,2])-1))))
pec_yank_v_minwhich$species <- substring(names(pec_yank$vertical), 1, 2)
aggregate(yank_minwhich~species, data = pec_yank_v_minwhich, function(x) c(mean = mean(x), sd = sd(x)))

# mediolateral 
pec_yank_ml_min <- data.frame(yank_min = unlist(cbind(lapply(pec_yank$mediolateral, FUN = function(x) min(x[,2], na.rm = TRUE)))))
pec_yank_ml_min$species <- substring(names(pec_yank$mediolateral), 1, 2)
pec_yank_ml_min$individual <- substring(names(pec_yank$mediolateral), 1, 4)
aggregate(yank_min~species, data = pec_yank_ml_min, function(x) c(mean = mean(x), sd = sd(x)))

pec_yank_ml_minwhich <- data.frame(yank_minwhich = unlist(cbind(lapply(pec_yank$mediolateral, FUN = function(x) which.min(x[,2])-1))))
pec_yank_ml_minwhich$species <- substring(names(pec_yank$mediolateral), 1, 2)
aggregate(yank_minwhich~species, data = pec_yank_ml_minwhich, function(x) c(mean = mean(x), sd = sd(x)))

# anteroposterior
pec_yank_ap_min <- data.frame(yank_min = unlist(cbind(lapply(pec_yank$anteroposterior, FUN = function(x) min(x[,2], na.rm = TRUE)))))
pec_yank_ap_min$species <- substring(names(pec_yank$anteroposterior), 1, 2)
pec_yank_ap_min$individual <- substring(names(pec_yank$anteroposterior), 1, 4)
aggregate(yank_min~species, data = pec_yank_ap_min, function(x) c(mean = mean(x), sd = sd(x)))

pec_yank_ap_minwhich <- data.frame(yank_minwhich = unlist(cbind(lapply(pec_yank$anteroposterior, FUN = function(x) which.min(x[,2])-1))))
pec_yank_ap_minwhich$species <- substring(names(pec_yank$anteroposterior), 1, 2)
aggregate(yank_minwhich~species, data = pec_yank_ap_minwhich, function(x) c(mean = mean(x), sd = sd(x)))

# net
pec_yank_net_min <- data.frame(yank_min = unlist(cbind(lapply(pec_yank$net, FUN = function(x) min(x[,2], na.rm = TRUE)))))
pec_yank_net_min$species <- substring(names(pec_yank$net), 1, 2)
pec_yank_net_min$individual <- substring(names(pec_yank$net), 1, 4)
aggregate(yank_min~species, data = pec_yank_net_min, function(x) c(mean = mean(x), sd = sd(x)))

pec_yank_net_minwhich <- data.frame(yank_minwhich = unlist(cbind(lapply(pec_yank$net, FUN = function(x) which.min(x[,2])-1))))
pec_yank_net_minwhich$species <- substring(names(pec_yank$net), 1, 2)
aggregate(yank_minwhich~species, data = pec_yank_net_minwhich, function(x) c(mean = mean(x), sd = sd(x)))


### PELVIC
## substrating one from which.max because the first observation is 0% of stance

## Maximum
# vertical 
pel_yank_v_max <- data.frame(yank_max = unlist(cbind(lapply(pel_yank$vertical, FUN = function(x) max(x[,2], na.rm = TRUE)))))
pel_yank_v_max$species <- substring(names(pel_yank$vertical), 1, 2)
pel_yank_v_max$individual <- substring(names(pel_yank$vertical), 1, 4)
aggregate(yank_max~species, data = pel_yank_v_max, function(x) c(mean = mean(x), sd = sd(x)))

pel_yank_v_maxwhich <- data.frame(yank_maxwhich = unlist(cbind(lapply(pel_yank$vertical, FUN = function(x) which.max(x[,2])-1))))
pel_yank_v_maxwhich$species <- substring(names(pel_yank$vertical), 1, 2)
aggregate(yank_maxwhich~species, data = pel_yank_v_maxwhich, function(x) c(mean = mean(x), sd = sd(x)))

# mediolateral 
pel_yank_ml_max <- data.frame(yank_max = unlist(cbind(lapply(pel_yank$mediolateral, FUN = function(x) max(x[,2], na.rm = TRUE)))))
pel_yank_ml_max$species <- substring(names(pel_yank$mediolateral), 1, 2)
pel_yank_ml_max$individual <- substring(names(pel_yank$mediolateral), 1, 4)
aggregate(yank_max~species, data = pel_yank_ml_max, function(x) c(mean = mean(x), sd = sd(x)))

pel_yank_ml_maxwhich <- data.frame(yank_maxwhich = unlist(cbind(lapply(pel_yank$mediolateral, FUN = function(x) which.max(x[,2])-1))))
pel_yank_ml_maxwhich$species <- substring(names(pel_yank$mediolateral), 1, 2)
aggregate(yank_maxwhich~species, data = pel_yank_ml_maxwhich, function(x) c(mean = mean(x), sd = sd(x)))

# anteroposterior
pel_yank_ap_max <- data.frame(yank_max = unlist(cbind(lapply(pel_yank$anteroposterior, FUN = function(x) max(x[,2], na.rm = TRUE)))))
pel_yank_ap_max$species <- substring(names(pel_yank$anteroposterior), 1, 2)
pel_yank_ap_max$individual <- substring(names(pel_yank$anteroposterior), 1, 4)
aggregate(yank_max~species, data = pel_yank_ap_max, function(x) c(mean = mean(x), sd = sd(x)))

pel_yank_ap_maxwhich <- data.frame(yank_maxwhich = unlist(cbind(lapply(pel_yank$anteroposterior, FUN = function(x) which.max(x[,2])-1))))
pel_yank_ap_maxwhich$species <- substring(names(pel_yank$anteroposterior), 1, 2)
aggregate(yank_maxwhich~species, data = pel_yank_ap_maxwhich, function(x) c(mean = mean(x), sd = sd(x)))

# net
pel_yank_net_max <- data.frame(yank_max = unlist(cbind(lapply(pel_yank$net, FUN = function(x) max(x[,2], na.rm = TRUE)))))
pel_yank_net_max$species <- substring(names(pel_yank$net), 1, 2)
pel_yank_net_max$individual <- substring(names(pel_yank$net), 1, 4)
aggregate(yank_max~species, data = pel_yank_net_max, function(x) c(mean = mean(x), sd = sd(x)))

pel_yank_net_maxwhich <- data.frame(yank_maxwhich = unlist(cbind(lapply(pel_yank$net, FUN = function(x) which.max(x[,2])-1))))
pel_yank_net_maxwhich$species <- substring(names(pel_yank$net), 1, 2)
aggregate(yank_maxwhich~species, data = pel_yank_net_maxwhich, function(x) c(mean = mean(x), sd = sd(x)))

## Minimum
# vertical 
pel_yank_v_min <- data.frame(yank_min = unlist(cbind(lapply(pel_yank$vertical, FUN = function(x) min(x[,2], na.rm = TRUE)))))
pel_yank_v_min$species <- substring(names(pel_yank$vertical), 1, 2)
pel_yank_v_min$individual <- substring(names(pel_yank$vertical), 1, 4)
aggregate(yank_min~species, data = pel_yank_v_min, function(x) c(mean = mean(x), sd = sd(x)))

pel_yank_v_minwhich <- data.frame(yank_minwhich = unlist(cbind(lapply(pel_yank$vertical, FUN = function(x) which.min(x[,2])-1))))
pel_yank_v_minwhich$species <- substring(names(pel_yank$vertical), 1, 2)
aggregate(yank_minwhich~species, data = pel_yank_v_minwhich, function(x) c(mean = mean(x), sd = sd(x)))

# mediolateral 
pel_yank_ml_min <- data.frame(yank_min = unlist(cbind(lapply(pel_yank$mediolateral, FUN = function(x) min(x[,2], na.rm = TRUE)))))
pel_yank_ml_min$species <- substring(names(pel_yank$mediolateral), 1, 2)
pel_yank_ml_min$individual <- substring(names(pel_yank$mediolateral), 1, 4)
aggregate(yank_min~species, data = pel_yank_ml_min, function(x) c(mean = mean(x), sd = sd(x)))

pel_yank_ml_minwhich <- data.frame(yank_minwhich = unlist(cbind(lapply(pel_yank$mediolateral, FUN = function(x) which.min(x[,2])-1))))
pel_yank_ml_minwhich$species <- substring(names(pel_yank$mediolateral), 1, 2)
aggregate(yank_minwhich~species, data = pel_yank_ml_minwhich, function(x) c(mean = mean(x), sd = sd(x)))

# anteroposterior
pel_yank_ap_min <- data.frame(yank_min = unlist(cbind(lapply(pel_yank$anteroposterior, FUN = function(x) min(x[,2], na.rm = TRUE)))))
pel_yank_ap_min$species <- substring(names(pel_yank$anteroposterior), 1, 2)
pel_yank_ap_min$individual <- substring(names(pel_yank$anteroposterior), 1, 4)
aggregate(yank_min~species, data = pel_yank_ap_min, function(x) c(mean = mean(x), sd = sd(x)))

pel_yank_ap_minwhich <- data.frame(yank_minwhich = unlist(cbind(lapply(pel_yank$anteroposterior, FUN = function(x) which.min(x[,2])-1))))
pel_yank_ap_minwhich$species <- substring(names(pel_yank$anteroposterior), 1, 2)
aggregate(yank_minwhich~species, data = pel_yank_ap_minwhich, function(x) c(mean = mean(x), sd = sd(x)))

# net
pel_yank_net_min <- data.frame(yank_min = unlist(cbind(lapply(pel_yank$net, FUN = function(x) min(x[,2], na.rm = TRUE)))))
pel_yank_net_min$species <- substring(names(pel_yank$net), 1, 2)
pel_yank_net_min$individual <- substring(names(pel_yank$net), 1, 4)
aggregate(yank_min~species, data = pel_yank_net_min, function(x) c(mean = mean(x), sd = sd(x)))

pel_yank_net_minwhich <- data.frame(yank_minwhich = unlist(cbind(lapply(pel_yank$net, FUN = function(x) which.min(x[,2])-1))))
pel_yank_net_minwhich$species <- substring(names(pel_yank$net), 1, 2)
aggregate(yank_minwhich~species, data = pel_yank_net_minwhich, function(x) c(mean = mean(x), sd = sd(x)))


## Zu omega squared to assess the 'goodness of fit' for the entire model
# Xu's omega method: http://onlinelibrary.wiley.com/doi/10.1002/sim.1572/abstract
## or through the performance package: (got the same exact results as my code)
# performance::r2_xu(pec_LMM)

# Xu_omega2 <- function(lmm, ...) {
#   1-var(residuals(lmm))/(var(model.response(model.frame(lmm))))
# }



#### YANK - PEC - LMM ####

## double-check number of trials in each species
table(factor(pec_yank_v_max$species, levels = c("af", "pw", "pb"))) 

### Maximum - magnitude 
# vertical 
pec_yank_v_max_lmm <- lmer(yank_max ~ species + (1|individual), data = pec_yank_v_max)
pec_yank_v_max_emm <- emmeans(pec_yank_v_max_lmm, "species")
pairs(pec_yank_v_max_emm)
pec_yank_v_max_lmm_omega2 <- performance::r2_xu(pec_yank_v_max_lmm) # 0.2424792
performance::r2_nakagawa(pec_yank_v_max_lmm) # c = 0.218, m = 0.083

# medioateral
pec_yank_ml_max_lmm <- lmer(yank_max ~ species + (1|individual), data = pec_yank_ml_max)
pec_yank_ml_max_emm <- emmeans(pec_yank_ml_max_lmm, "species")
pairs(pec_yank_ml_max_emm)
pec_yank_ml_max_lmm_omega2 <- performance::r2_xu(pec_yank_ml_max_lmm) # 0.418855
performance::r2_nakagawa(pec_yank_ml_max_lmm) # c = 0.391, m = 0.122

# anteroposterior
pec_yank_ap_max_lmm <- lmer(yank_max ~ species + (1|individual), data = pec_yank_ap_max)
pec_yank_ap_max_emm <- emmeans(pec_yank_ap_max_lmm, "species")
pairs(pec_yank_ap_max_emm)
pec_yank_ap_max_lmm_omega2 <- performance::r2_xu(pec_yank_ap_max_lmm) # 0.3946568
performance::r2_nakagawa(pec_yank_ap_max_lmm) # c = 0.366, m = 0.243

# net
pec_yank_net_max_lmm <- lmer(yank_max ~ species + (1|individual), data = pec_yank_net_max)
pec_yank_net_max_emm <- emmeans(pec_yank_net_max_lmm, "species")
pairs(pec_yank_net_max_emm)
pec_yank_net_max_lmm_omega2 <- performance::r2_xu(pec_yank_net_max_lmm) # 0.2503454
performance::r2_nakagawa(pec_yank_net_max_lmm) # c = 0.230, m = 0.060


### Minimum - magnitude

# vertical
pec_yank_v_min_lmm <- lmer(yank_min ~ species + (1|individual), data = pec_yank_v_min)
pec_yank_v_min_emm <- emmeans(pec_yank_v_min_lmm, "species")
pairs(pec_yank_v_min_emm)
pec_yank_v_min_lmm_omega2 <- performance::r2_xu(pec_yank_v_min_lmm) # 0.1692728
performance::r2_nakagawa(pec_yank_v_min_lmm) # c = 0.152, 0.109

# mediolateral
pec_yank_ml_min_lmm <- lmer(yank_min ~ species + (1|individual), data = pec_yank_ml_min)
pec_yank_ml_min_emm <- emmeans(pec_yank_ml_min_lmm, "species")
pairs(pec_yank_ml_min_emm)
pec_yank_ml_min_lmm_omega2 <- performance::r2_xu(pec_yank_ml_min_lmm) # 0.3740738 
performance::r2_nakagawa(pec_yank_ml_min_lmm) # c = 0.363, m = 0.118

# anteroposterior
pec_yank_ap_min_lmm <- lmer(yank_min ~ species + (1|individual), data = pec_yank_ap_min)
pec_yank_ap_min_emm <- emmeans(pec_yank_ap_min_lmm, "species")
pairs(pec_yank_ap_min_emm)
pec_yank_ap_min_lmm_omega2 <- performance::r2_xu(pec_yank_ap_min_lmm) # 0.3512727
performance::r2_nakagawa(pec_yank_ap_min_lmm) # c = 0.335, m = 0.165

# net
pec_yank_net_min_lmm <- lmer(yank_min ~ species + (1|individual), data = pec_yank_net_min)
pec_yank_net_min_emm <- emmeans(pec_yank_net_min_lmm, "species")
pairs(pec_yank_net_min_emm)
pec_yank_net_min_lmm_omega2 <- performance::r2_xu(pec_yank_net_min_lmm) # 0.1636566
performance::r2_nakagawa(pec_yank_net_min_lmm) # c = 0.145, m = 0.081



#### YANK - PEL - LMM ####

## double-check number of trials in each species
table(factor(pel_yank_v_max$species, levels = c("af", "pw", "pb"))) 

### Maximum - magnitude 
# vertical 
pel_yank_v_max_lmm <- lmer(yank_max ~ species + (1|individual), data = pel_yank_v_max)
pel_yank_v_max_emm <- emmeans(pel_yank_v_max_lmm, "species")
pairs(pel_yank_v_max_emm)
pel_yank_v_max_lmm_omega2 <- performance::r2_xu(pel_yank_v_max_lmm) # 0.2254354
performance::r2_nakagawa(pel_yank_v_max_lmm) # c = 0.204, m = 0.106

# medioateral
pel_yank_ml_max_lmm <- lmer(yank_max ~ species + (1|individual), data = pel_yank_ml_max)
pel_yank_ml_max_emm <- emmeans(pel_yank_ml_max_lmm, "species")
pairs(pel_yank_ml_max_emm)
pel_yank_ml_max_lmm_omega2 <- performance::r2_xu(pel_yank_ml_max_lmm) # .3096737
performance::r2_nakagawa(pel_yank_ml_max_lmm) # c = 0.314,  m = 0.001

# anteroposterior
pel_yank_ap_max_lmm <- lmer(yank_max ~ species + (1|individual), data = pel_yank_ap_max)
pel_yank_ap_max_emm <- emmeans(pel_yank_ap_max_lmm, "species")
pairs(pel_yank_ap_max_emm)
pel_yank_ap_max_lmm_omega2 <- performance::r2_xu(pel_yank_ap_max_lmm) # 0.2067791
performance::r2_nakagawa(pel_yank_ap_max_lmm) # c = 0.171, m = 0.083

# net
pel_yank_net_max_lmm <- lmer(yank_max ~ species + (1|individual), data = pel_yank_net_max)
pel_yank_net_max_emm <- emmeans(pel_yank_net_max_lmm, "species")
pairs(pel_yank_net_max_emm)
pel_yank_net_max_lmm_omega2 <- performance::r2_xu(pel_yank_net_max_lmm) # 0.1992901
performance::r2_nakagawa(pel_yank_net_max_lmm) # c = 0.180, m = 0.102


### Minimum - magnitude

# vertical
pel_yank_v_min_lmm <- lmer(yank_min ~ species + (1|individual), data = pel_yank_v_min)
pel_yank_v_min_emm <- emmeans(pel_yank_v_min_lmm, "species")
pairs(pel_yank_v_min_emm)
pel_yank_v_min_lmm_omega2 <- performance::r2_xu(pel_yank_v_min_lmm) # 0.6515398
performance::r2_nakagawa(pel_yank_v_min_lmm) # c = 0.629, m = 0.401

# mediolateral
pel_yank_ml_min_lmm <- lmer(yank_min ~ species + (1|individual), data = pel_yank_ml_min)
pel_yank_ml_min_emm <- emmeans(pel_yank_ml_min_lmm, "species")
pairs(pel_yank_ml_min_emm)
pel_yank_ml_min_lmm_omega2 <- performance::r2_xu(pel_yank_ml_min_lmm) # 0.07954376
performance::r2_nakagawa(pel_yank_ml_min_lmm) # c = 0.063, m = 0.001

# anteroposterior
pel_yank_ap_min_lmm <- lmer(yank_min ~ species + (1|individual), data = pel_yank_ap_min)
pel_yank_ap_min_emm <- emmeans(pel_yank_ap_min_lmm, "species")
pairs(pel_yank_ap_min_emm)
pel_yank_ap_min_lmm_omega2 <- performance::r2_xu(pel_yank_ap_min_lmm) # 0.2469142
performance::r2_nakagawa(pel_yank_ap_min_lmm) # c = 0.208, m = 0.035

# net
pel_yank_net_min_lmm <- lmer(yank_min ~ species + (1|individual), data = pel_yank_net_min)
pel_yank_net_min_emm <- emmeans(pel_yank_net_min_lmm, "species")
pairs(pel_yank_net_min_emm)
pel_yank_net_min_lmm_omega2 <- performance::r2_xu(pel_yank_net_min_lmm) # 0.6509524
performance::r2_nakagawa(pel_yank_net_min_lmm) # c = 0.619, m = 0.454


### YANK - SHAPIRO-WILK ####

## a) Don't need to test for linearity of data because the predictors are categorical

## b) evaluating the normality of the residuals
# the null of the Shapiro-Wilk test is that the input (e.g., residuals of data) are normal

## Pec - max
shapiro.test(resid(pec_yank_v_max_lmm)) # p-value = 0.871
shapiro.test(resid(pec_yank_ml_max_lmm)) # p-value = 0.0000052
shapiro.test(resid(pec_yank_ap_max_lmm)) # p-value =  0.0000001552
shapiro.test(resid(pec_yank_net_max_lmm)) # p-value = 0.866

## Pec - min
shapiro.test(resid(pec_yank_v_min_lmm)) # p-value = 0.4643
shapiro.test(resid(pec_yank_ml_min_lmm)) # p-value = 0.0006225
shapiro.test(resid(pec_yank_ap_min_lmm)) # p-value =  0.0006075
shapiro.test(resid(pec_yank_net_min_lmm)) # p-value = 0.4183

## Pel - max
shapiro.test(resid(pel_yank_v_max_lmm)) # p-value = 0.00000358
shapiro.test(resid(pel_yank_ml_max_lmm)) # p-value = 0.02667
shapiro.test(resid(pel_yank_ap_max_lmm)) # p-value = 0.000001058
shapiro.test(resid(pel_yank_net_max_lmm)) # p-value = 0.000011

## Pel - min
shapiro.test(resid(pel_yank_v_min_lmm)) # p-value = 0.07706
shapiro.test(resid(pel_yank_ml_min_lmm)) # p-value = 0.00008501
shapiro.test(resid(pel_yank_ap_min_lmm)) # p-value = 0.001612
shapiro.test(resid(pel_yank_net_min_lmm)) # p-value = 0.0585


#### YANK - QQ PLOT ASSUMPTIONS ####
## (export as 500 width x 600 height)

### PECTORAL - MAX
pec_yank_v_max_lmm_QQ <- plotQQ(pec_yank_v_max_lmm)
pec_yank_v_max_lmm_Fitted <- plot(pec_yank_v_max_lmm)

pec_yank_ml_max_lmm_QQ <- plotQQ(pec_yank_ml_max_lmm)
pec_yank_ml_max_lmm_Fitted <- plot(pec_yank_ml_max_lmm)

pec_yank_ap_max_lmm_QQ <- plotQQ(pec_yank_ap_max_lmm)
pec_yank_ap_max_lmm_Fitted <- plot(pec_yank_ap_max_lmm)

pec_yank_net_max_lmm_QQ <- plotQQ(pec_yank_net_max_lmm)
pec_yank_net_max_lmm_Fitted <- plot(pec_yank_net_max_lmm)


### PECTORAL - MIN
pec_yank_v_min_lmm_QQ <- plotQQ(pec_yank_v_min_lmm)
pec_yank_v_min_lmm_Fitted <- plot(pec_yank_v_min_lmm)

pec_yank_ml_min_lmm_QQ <- plotQQ(pec_yank_ml_min_lmm)
pec_yank_ml_min_lmm_Fitted <- plot(pec_yank_ml_min_lmm)

pec_yank_ap_min_lmm_QQ <- plotQQ(pec_yank_ap_min_lmm)
pec_yank_ap_min_lmm_Fitted <- plot(pec_yank_ap_min_lmm)

pec_yank_net_min_lmm_QQ <- plotQQ(pec_yank_net_min_lmm)
pec_yank_net_min_lmm_Fitted <- plot(pec_yank_net_min_lmm)


### PELVIC - MAX
pel_yank_v_max_lmm_QQ <- plotQQ(pel_yank_v_max_lmm)
pel_yank_v_max_lmm_Fitted <- plot(pel_yank_v_max_lmm)

pel_yank_ml_max_lmm_QQ <- plotQQ(pel_yank_ml_max_lmm)
pel_yank_ml_max_lmm_Fitted <- plot(pel_yank_ml_max_lmm)

pel_yank_ap_max_lmm_QQ <- plotQQ(pel_yank_ap_max_lmm)
pel_yank_ap_max_lmm_Fitted <- plot(pel_yank_ap_max_lmm)

pel_yank_net_max_lmm_QQ <- plotQQ(pel_yank_net_max_lmm)
pel_yank_net_max_lmm_Fitted <- plot(pel_yank_net_max_lmm)


### PELVIC - MIN
pel_yank_v_min_lmm_QQ <- plotQQ(pel_yank_v_min_lmm)
pel_yank_v_min_lmm_Fitted <- plot(pel_yank_v_min_lmm)

pel_yank_ml_min_lmm_QQ <- plotQQ(pel_yank_ml_min_lmm)
pel_yank_ml_min_lmm_Fitted <- plot(pel_yank_ml_min_lmm)

pel_yank_ap_min_lmm_QQ <- plotQQ(pel_yank_ap_min_lmm)
pel_yank_ap_min_lmm_Fitted <- plot(pel_yank_ap_min_lmm)

pel_yank_net_min_lmm_QQ <- plotQQ(pel_yank_net_min_lmm)
pel_yank_net_min_lmm_Fitted <- plot(pel_yank_net_min_lmm)


jpeg("pec_yank_max_QQ.jpg", width = 10, height = 6, units = "in", res = 300)
cowplot::plot_grid(pec_yank_net_max_lmm_QQ, pec_yank_v_max_lmm_QQ, pec_yank_ml_max_lmm_QQ, pec_yank_ap_max_lmm_QQ, labels = c("a", "b", "c", "d"))
dev.off()

jpeg("pec_yank_min_QQ.jpg", width = 10, height = 6, units = "in", res = 300)
cowplot::plot_grid(pec_yank_net_min_lmm_QQ, pec_yank_v_min_lmm_QQ, pec_yank_ml_min_lmm_QQ, pec_yank_ap_min_lmm_QQ, labels = c("a", "b", "c", "d"))
dev.off()

jpeg("pel_yank_max_QQ.jpg", width = 10, height = 6, units = "in", res = 300)
cowplot::plot_grid(pel_yank_net_max_lmm_QQ, pel_yank_v_max_lmm_QQ, pel_yank_ml_max_lmm_QQ, pel_yank_ap_max_lmm_QQ, labels = c("a", "b", "c", "d"))
dev.off()

jpeg("pel_yank_min_QQ.jpg", width = 10, height = 6, units = "in", res = 300)
cowplot::plot_grid(pel_yank_net_min_lmm_QQ, pel_yank_v_min_lmm_QQ, pel_yank_ml_min_lmm_QQ, pel_yank_ap_min_lmm_QQ, labels = c("a", "b", "c", "d"))
dev.off()


#### YANK - REMOVING OUTLIERS ####
## outliers are set as those points falling outside 2x the interquantile range
## only doing this for the variables that seemed to deviate from normality

# pec - ml - max
pec_yank_ml_max_bp <-  boxplot(yank_max ~ species, data = pec_yank_ml_max, range = 2)
pec_yank_ml_max_noOutliers <- subset(pec_yank_ml_max, !yank_max %in% c(pec_yank_ml_max_bp$out))
pec_yank_ml_max_noOutliers_lmm <- lmer(yank_max ~ species + (1|individual), data = pec_yank_ml_max_noOutliers)
plotQQ(pec_yank_ml_max_noOutliers_lmm)
shapiro.test(resid(pec_yank_ml_max_noOutliers_lmm))
# the QQ plot is improved by removing the two outliers
pec_yank_ml_max_noOutliers_emm <- emmeans(pec_yank_ml_max_noOutliers_lmm, "species")

# pec - ap - max
pec_yank_ap_max_bp <-  boxplot(yank_max ~ species, data = pec_yank_ap_max, range = 2)
pec_yank_ap_max_noOutliers <- subset(pec_yank_ap_max, !yank_max %in% c(pec_yank_ap_max_bp$out))
pec_yank_ap_max_noOutliers_lmm <- lmer(yank_max ~ species + (1|individual), data = pec_yank_ap_max_noOutliers)
plotQQ(pec_yank_ap_max_noOutliers_lmm)
shapiro.test(resid(pec_yank_ap_max_noOutliers_lmm))
# the QQ plot is improved by removing the one outlier
pec_yank_ap_max_noOutliers_emm <- emmeans(pec_yank_ap_max_noOutliers_lmm, "species")

# pec - ml - min
pec_yank_ml_min_bp <-  boxplot(yank_min ~ species, data = pec_yank_ml_min, range = 2)
pec_yank_ml_min_noOutliers <- subset(pec_yank_ml_min, !yank_min %in% c(pec_yank_ml_min_bp$out))
pec_yank_ml_min_noOutliers_lmm <- lmer(yank_min ~ species + (1|individual), data = pec_yank_ml_min_noOutliers)
plotQQ(pec_yank_ml_min_noOutliers_lmm)
shapiro.test(resid(pec_yank_ml_min_noOutliers_lmm))
# the QQ plot is improved by removing the three outliers
pec_yank_ml_min_noOutliers_emm <- emmeans(pec_yank_ml_min_noOutliers_lmm, "species")

# pec - ap - min
pec_yank_ap_min_bp <-  boxplot(yank_min ~ species, data = pec_yank_ap_min, range = 2)
pec_yank_ap_min_noOutliers <- subset(pec_yank_ap_min, !yank_min %in% c(pec_yank_ap_min_bp$out))
pec_yank_ap_min_noOutliers_lmm <- lmer(yank_min ~ species + (1|individual), data = pec_yank_ap_min_noOutliers)
plotQQ(pec_yank_ap_min_noOutliers_lmm)
shapiro.test(resid(pec_yank_ap_min_noOutliers_lmm))
# the QQ plot is improved by removing the two outliers
pec_yank_ap_min_noOutliers_emm <- emmeans(pec_yank_ap_min_noOutliers_lmm, "species")


# pel - v - max
pel_yank_v_max_bp <-  boxplot(yank_max ~ species, data = pel_yank_v_max, range = 2)
pel_yank_v_max_noOutliers <- subset(pel_yank_v_max, !yank_max %in% c(pel_yank_v_max_bp$out))
pel_yank_v_max_noOutliers_lmm <- lmer(yank_max ~ species + (1|individual), data = pel_yank_v_max_noOutliers)
plotQQ(pel_yank_v_max_noOutliers_lmm)
shapiro.test(resid(pel_yank_v_max_noOutliers_lmm))
# the QQ plot didn't really change
pel_yank_v_max_noOutliers_emm <- emmeans(pel_yank_v_max_noOutliers_lmm, "species")


# pel - ap - max
pel_yank_ap_max_bp <-  boxplot(yank_max ~ species, data = pel_yank_ap_max, range = 2)
pel_yank_ap_max_noOutliers <- subset(pel_yank_ap_max, !yank_max %in% c(pel_yank_ap_max_bp$out))
pel_yank_ap_max_noOutliers_lmm <- lmer(yank_max ~ species + (1|individual), data = pel_yank_ap_max_noOutliers)
plotQQ(pel_yank_ap_max_noOutliers_lmm)
shapiro.test(resid(pel_yank_ap_max_noOutliers_lmm))
# the QQ plot is improved by removing the 6 outliers
pel_yank_ap_max_noOutliers_emm <- emmeans(pel_yank_ap_max_noOutliers_lmm, "species")


# pel - net - max
pel_yank_net_max_bp <-  boxplot(yank_max ~ species, data = pel_yank_net_max, range = 2)
pel_yank_net_max_noOutliers <- subset(pel_yank_net_max, !yank_max %in% c(pel_yank_net_max_bp$out))
pel_yank_net_max_noOutliers_lmm <- lmer(yank_max ~ species + (1|individual), data = pel_yank_net_max_noOutliers)
plotQQ(pel_yank_net_max_noOutliers_lmm)
shapiro.test(resid(pel_yank_net_max_noOutliers_lmm))
# the QQ plot didn't really change
pel_yank_net_max_noOutliers_emm <- emmeans(pel_yank_net_max_noOutliers_lmm, "species")


# pel - ml - min
pel_yank_ml_min_bp <-  boxplot(yank_min ~ species, data = pel_yank_ml_min, range = 2)
pel_yank_ml_min_noOutliers <- subset(pel_yank_ml_min, !yank_min %in% c(pel_yank_ml_min_bp$out))
pel_yank_ml_min_noOutliers_lmm <- lmer(yank_min ~ species + (1|individual), data = pel_yank_ml_min_noOutliers)
plotQQ(pel_yank_ml_min_noOutliers_lmm)
shapiro.test(resid(pel_yank_ml_min_noOutliers_lmm))
# the QQ plot is improved by removing the two outliers
pel_yank_ml_min_noOutliers_emm <- emmeans(pel_yank_ml_min_noOutliers_lmm, "species")


# pel - ap - min
pel_yank_ap_min_bp <-  boxplot(yank_min ~ species, data = pel_yank_ap_min, range = 2)
pel_yank_ap_min_noOutliers <- subset(pel_yank_ap_min, !yank_min %in% c(pel_yank_ap_min_bp$out))
pel_yank_ap_min_noOutliers_lmm <- lmer(yank_min ~ species + (1|individual), data = pel_yank_ap_min_noOutliers)
plotQQ(pel_yank_ap_min_noOutliers_lmm)
shapiro.test(resid(pel_yank_ap_min_noOutliers_lmm))
# the QQ plot is improved by removing the two outliers
pel_yank_ap_min_noOutliers_emm <- emmeans(pel_yank_ap_min_noOutliers_lmm, "species")



  #### SAVING THE DATA ####
  ## go to the parent directory then save output in 'output' folder
  setwd('..')
  setwd('..')
  setwd('./output')
  
  ## Pectoral data
  
  ## Save the dataset that was filtered and had areas of overlap included
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
  
  
  
  ##### COMPARING STANCE DURATIONS ####
  
  VideoInfo <- GRFs$VideoInfo
  names(VideoInfo)[1] <- "filename"
  VideoInfo$individual <- substring(VideoInfo$filename, 1,4)
  
  # Separate info between pectoral vs. pelvic trials
  Pec_VideoInfo <- VideoInfo[VideoInfo$TrialsToUse == "both"|VideoInfo$TrialsToUse == "pec",]
  Pel_VideoInfo <- VideoInfo[VideoInfo$TrialsToUse == "both"|VideoInfo$TrialsToUse == "pel",]
  
  Pec_VideoInfo$stance_s <- (Pec_VideoInfo$Pectoral.End.Frame - Pec_VideoInfo$Pectoral.Start.Frame)/Pec_VideoInfo$Filming.Rate.Hz
  Pel_VideoInfo$stance_s <- (Pel_VideoInfo$Pelvic.End.Frame - Pel_VideoInfo$Pelvic.Start.Frame)/Pel_VideoInfo$Filming.Rate.Hz
  
  
  ## double-check number of trials in each species
  table(factor(Pec_VideoInfo$Species, levels = c("Ambystoma_tigrinum", "Pleurodeles_waltl", "Periophthalmus_barbarus"))) 
  table(factor(Pel_VideoInfo$Species, levels = c("Ambystoma_tigrinum", "Pleurodeles_waltl", "Periophthalmus_barbarus"))) 
  
  
  ### stance comparisons
  pec_stance_lmm <- lmer(stance_s ~ Species + (1|individual), data = Pec_VideoInfo)
  pec_stance_emm <- emmeans(pec_stance_lmm, "Species")
  pairs(pec_stance_emm)
  pec_stance_lmm_omega2 <- performance::r2_xu(pec_stance_lmm) # 0.4148974
  performance::r2_nakagawa(pec_stance_lmm) # c = 0.392, m = 0.271
  (summary(pec_stance_lmm)$coefficients[1,1]/(summary(pec_stance_lmm)$coefficients[2,1]+summary(pec_stance_lmm)$coefficients[1,1]))*100
  # stance duration is only ~35.6% different between af and pb
  (summary(pec_stance_lmm)$coefficients[1,1]/(summary(pec_stance_lmm)$coefficients[3,1]+summary(pec_stance_lmm)$coefficients[1,1]))*100
  # stance duration is only ~3.5% different between af and pw
  ((summary(pec_stance_lmm)$coefficients[3,1]+summary(pec_stance_lmm)$coefficients[1,1])/(summary(pec_stance_lmm)$coefficients[2,1]+summary(pec_stance_lmm)$coefficients[1,1]))*100
  # stance duration is only ~40.5% different between af and pb
  
  pel_stance_lmm <- lmer(stance_s ~ Species + (1|individual), data = Pel_VideoInfo)
  pel_stance_emm <- emmeans(pel_stance_lmm, "Species")
  pairs(pel_stance_emm)
  pel_stance_lmm_omega2 <- performance::r2_xu(pel_stance_lmm) # 0.1784667
  performance::r2_nakagawa(pel_stance_lmm) # c = 0.157, m = 0.002
 (summary(pel_stance_lmm)$coefficients[1,1]/(summary(pel_stance_lmm)$coefficients[2,1]+summary(pel_stance_lmm)$coefficients[1,1]))*100
  # stance duration is only ~2.6% different between the species
  
  
  #### RUNNING STANCE DURATION AS COVARIATE IN LMMS - PEAK NET GRF ####
  ## Using likelihood ratio tests to determine whether stance needs to be included in the model
  
  pec_peakNetGRFs_wStance <- merge(pec_peakNetGRFs, Pec_VideoInfo, by = c("filename", "individual"))

  
  ### PECTORAL
  # vertical 
  pec_peakNetGRF_v_stance_lmm <- lmer(InterpV_BW ~ species*stance_s + (1|individual), data = pec_peakNetGRFs_wStance)
  pec_peakNetGRF_v_stance_lmm_ML <- lmer(InterpV_BW ~ species*stance_s + (1|individual), data = pec_peakNetGRFs_wStance, REML = FALSE)
  pec_peakNetGRF_v_lmm_ML <- lmer(InterpV_BW ~ species + (1|individual), data = pec_peakNetGRFs_wStance, REML = FALSE)
  anova(pec_peakNetGRF_v_lmm_ML, pec_peakNetGRF_v_stance_lmm_ML) # p-value = 0.4905
  
  pec_peakNetGRF_v_stance_emm <- emmeans(pec_peakNetGRF_v_stance_lmm, ~species|stance_s)
  pairs(pec_peakNetGRF_v_stance_emm)
  performance::r2_nakagawa(pec_peakNetGRF_v_stance_lmm) # c = 0.197, m = 0.148
  summary(pec_peakNetGRF_v_stance_lmm)
  Anova(pec_peakNetGRF_v_stance_lmm)
  
  # mediolateral
  pec_peakNetGRF_ml_stance_lmm <- lmer(InterpML_BW ~ species*stance_s + (1|individual), data = pec_peakNetGRFs_wStance)
  pec_peakNetGRF_ml_stance_lmm_ML <- lmer(InterpML_BW ~ species*stance_s + (1|individual), data = pec_peakNetGRFs_wStance, REML = FALSE)
  pec_peakNetGRF_ml_lmm_ML <- lmer(InterpML_BW ~ stance_s + (1|individual), data = pec_peakNetGRFs_wStance, REML = FALSE)
  anova(pec_peakNetGRF_ml_lmm_ML, pec_peakNetGRF_ml_stance_lmm_ML) # p-value 0.0004
  
  pec_peakNetGRF_ml_stance_emm <- emmeans(pec_peakNetGRF_ml_stance_lmm, ~species|stance_s)
  pairs(pec_peakNetGRF_ml_stance_emm)
  performance::r2_nakagawa(pec_peakNetGRF_ml_stance_lmm) # c = 0.529, m = 0.305
  summary(pec_peakNetGRF_ml_stance_lmm)
  Anova(pec_peakNetGRF_ml_stance_lmm)
  
  
  # anteroposterior
  pec_peakNetGRF_ap_stance_lmm <- lmer(InterpAP_BW ~ species*stance_s + (1|individual), data = pec_peakNetGRFs_wStance)
  pec_peakNetGRF_ap_stance_lmm_ML <- lmer(InterpAP_BW ~ species*stance_s + (1|individual), data = pec_peakNetGRFs_wStance, REML = FALSE)
  pec_peakNetGRF_ap_lmm_ML <- lmer(InterpAP_BW ~ stance_s + (1|individual), data = pec_peakNetGRFs_wStance, REML = FALSE)
  anova(pec_peakNetGRF_ap_lmm_ML, pec_peakNetGRF_ap_stance_lmm_ML) # p-value <0.0001
  
  pec_peakNetGRF_ap_stance_emm <- emmeans(pec_peakNetGRF_ap_stance_lmm, ~species|stance_s)
  pairs(pec_peakNetGRF_ap_stance_emm)
  performance::r2_nakagawa(pec_peakNetGRF_ap_stance_lmm) # c = 0.709, m = 0.594
  summary(pec_peakNetGRF_ap_stance_lmm)
  Anova(pec_peakNetGRF_ap_stance_lmm)
  
  # net
  pec_peakNetGRF_net_stance_lmm <- lmer(NetGRF_BW ~ species*stance_s + (1|individual), data = pec_peakNetGRFs_wStance)
  pec_peakNetGRF_net_stance_lmm_ML <- lmer(NetGRF_BW ~ species*stance_s + (1|individual), data = pec_peakNetGRFs_wStance, REML = FALSE)
  pec_peakNetGRF_net_lmm_ML <- lmer(NetGRF_BW ~ stance_s + (1|individual), data = pec_peakNetGRFs_wStance, REML = FALSE)
  anova(pec_peakNetGRF_net_lmm_ML, pec_peakNetGRF_net_stance_lmm_ML) # p-value = 0.09234
  
  pec_peakNetGRF_net_stance_emm <- emmeans(pec_peakNetGRF_net_stance_lmm, ~species|stance_s)
  pairs(pec_peakNetGRF_net_stance_emm)
  performance::r2_nakagawa(pec_peakNetGRF_net_stance_lmm) # c = 0.211, m = 0.113
  summary(pec_peakNetGRF_net_stance_lmm)
  Anova(pec_peakNetGRF_net_stance_lmm)
  
  # mediolateral angle
  pec_peakNetGRF_mlang_stance_lmm <- lmer(MLAngle_Convert_deg ~ species*stance_s + (1|individual), data = pec_peakNetGRFs_wStance)
  pec_peakNetGRF_mlang_stance_lmm_ML <- lmer(MLAngle_Convert_deg ~ species*stance_s + (1|individual), data = pec_peakNetGRFs_wStance, REML = FALSE)
  pec_peakNetGRF_mlang_lmm_ML <- lmer(MLAngle_Convert_deg ~ stance_s + (1|individual), data = pec_peakNetGRFs_wStance, REML = FALSE)
  anova(pec_peakNetGRF_mlang_lmm_ML, pec_peakNetGRF_mlang_stance_lmm_ML) # p-value = 0.0005425 
  
  pec_peakNetGRF_mlang_stance_emm <- emmeans(pec_peakNetGRF_mlang_stance_lmm, ~species|stance_s)
  pairs(pec_peakNetGRF_mlang_stance_emm)
  performance::r2_nakagawa(pec_peakNetGRF_mlang_stance_lmm) # c = 0.426, m = 0.295
  summary(pec_peakNetGRF_mlang_stance_lmm)
  Anova(pec_peakNetGRF_mlang_stance_lmm)
  
  # anteroposterior angle
  pec_peakNetGRF_apang_stance_lmm <- lmer(APAngle_Convert_deg ~ species*stance_s + (1|individual), data = pec_peakNetGRFs_wStance)
  pec_peakNetGRF_apang_stance_lmm_ML <- lmer(APAngle_Convert_deg ~ species*stance_s + (1|individual), data = pec_peakNetGRFs_wStance, REML = FALSE)
  pec_peakNetGRF_apang_lmm_ML <- lmer(APAngle_Convert_deg ~ stance_s + (1|individual), data = pec_peakNetGRFs_wStance, REML = FALSE)
  anova(pec_peakNetGRF_apang_lmm_ML, pec_peakNetGRF_apang_stance_lmm_ML) # p-value = 0.0000000006941
  
  pec_peakNetGRF_apang_stance_emm <- emmeans(pec_peakNetGRF_apang_stance_lmm, ~species|stance_s)
  pairs(pec_peakNetGRF_apang_stance_emm)
  performance::r2_nakagawa(pec_peakNetGRF_apang_stance_lmm) # c = 0.732, m = 0.653
  summary(pec_peakNetGRF_apang_stance_lmm)
  Anova(pec_peakNetGRF_apang_stance_lmm)
  
  # percent stance
  pec_peakNetGRF_percentstance_stance_lmm <- lmer(PercentStance ~ species*stance_s + (1|individual), data = pec_peakNetGRFs_wStance)
  pec_peakNetGRF_percentstance_stance_lmm_ML <- lmer(PercentStance ~ species*stance_s + (1|individual), data = pec_peakNetGRFs_wStance, REML = FALSE)
  pec_peakNetGRF_percentstance_lmm_ML <- lmer(PercentStance ~ stance_s + (1|individual), data = pec_peakNetGRFs_wStance, REML = FALSE)
  anova(pec_peakNetGRF_percentstance_lmm_ML, pec_peakNetGRF_percentstance_stance_lmm_ML) # p-value = 0.01524
  
  pec_peakNetGRF_percentstance_stance_emm <- emmeans(pec_peakNetGRF_percentstance_stance_lmm, species*stance_s)
  pairs(pec_peakNetGRF_percentstance_stance_emm)
  performance::r2_nakagawa(pec_peakNetGRF_percentstance_stance_lmm) # c = 0.469, m = 0.250
  summary(pec_peakNetGRF_percentstance_stance_lmm)
  Anova(pec_peakNetGRF_percentstance_stance_lmm)
  
  
  
  
  
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
