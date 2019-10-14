################  GRF of Fish and Salamander Appendages ########################
## Code for evaluating force data collected from the Blob lab force plate #####
# Abbreviations: Vert = Vertical, ML = mediolateral, Hz = Horizontal (refers to anteroposterior direction)
# Sometimes this file acts a little funny, so may need to select all text and then click the "Run selection" button

# 12-29-12: Editing code to produce prettier plots with standard error using ggplot2 package

######################  ORGANIZING POST-FILTERED RESULTS #############################

# Clear everything in the R workspace, so nothing gets mixed up between different trials
rm(list=ls(all=TRUE))

today <- Sys.Date()
SaveDate <- format(today, format="%y%m%d")


#####################  PEAK NET GRF #############################################
setwd('/Users/SandyMKawano/Desktop/Research/Forelimbs on Force Plates/Force Data/R Analysis/Step 3 Post Filter Processing/Peak net GRF')


PeakNetGRFPath <- setwd('/Users/SandyMKawano/Desktop/Research/Forelimbs on Force Plates/Force Data/R Analysis/Step 3 Post Filter Processing/Peak net GRF')
PeakNetGRF.Ind <- list.files(PeakNetGRFPath, pattern=".csv", full=TRUE)

# Loading all files
PeakNetGRF.All <- lapply(PeakNetGRF.Ind, read.csv, header=TRUE)

# Combining all the data into a single dataframe
PeakNetGRF.Combined <- data.frame(do.call("rbind", PeakNetGRF.All))

# Adding file names and appendage type as additional variables
PeakNetGRF.Combined$File.name <- substring(PeakNetGRF.Ind, 129, 135)
PeakNetGRF.Combined$Appendage <- substring(PeakNetGRF.Ind, 137, 139)
PeakNetGRF.Combined$Species <- substring(PeakNetGRF.Ind, 129,130)

# Dividing into pelvic vs. pectoral files
PeakNetGRF.Pec <- PeakNetGRF.Combined[grep('Pec', PeakNetGRF.Combined$Appendage),]
PeakNetGRF.Pel <- PeakNetGRF.Combined[grep('Pel', PeakNetGRF.Combined$Appendage),]

# Choosing only those that I've decided to keep
# Loading spreadsheet of which trials were kept

# Read the most recent Video Info file
setwd('/Users/SandyMKawano/Desktop/Research/Forelimbs on Force Plates/Trials Kept')
TrialsPath <- setwd('/Users/SandyMKawano/Desktop/Research/Forelimbs on Force Plates/Trials Kept')
Trials.Ind <- list.files(TrialsPath, pattern=".csv", full=TRUE)
Trials.Dates <-substring(Trials.Ind, 88, 93)
Trials.Newest <- max(as.numeric(Trials.Dates))
Trials.Use <- Trials.Ind[substring(Trials.Ind, 88, 93) %in% Trials.Newest]
TrialsNotes <- data.frame(read.csv(Trials.Use, header=TRUE))

# Subsetting data to organize data into trials kept by appendage
PecFiles <- TrialsNotes[which(TrialsNotes$Appendage=='Pectoral'),]
PecFilesKept <- PecFiles[which(PecFiles$Use=='Yes'),]
PecKept <- PeakNetGRF.Pec[which(PeakNetGRF.Pec$File.name %in% PecFilesKept$File.name),]

PelFiles <- TrialsNotes[which(TrialsNotes$Appendage=='Pelvic'),]
PelFilesKept <- PelFiles[which(PelFiles$Use=='Yes'),]
PelKept <- PeakNetGRF.Pel[which(PeakNetGRF.Pel$File.name %in% PelFilesKept$File.name),]

# Read the most recent Video Info file
setwd('/Users/SandyMKawano/Desktop/Research/Forelimbs on Force Plates/Video/Video Info')
VideoPath <- setwd('/Users/SandyMKawano/Desktop/Research/Forelimbs on Force Plates/Video/Video Info')
Video.Ind <- list.files(VideoPath, pattern=".csv", full=TRUE)
Video.Dates <-substring(Video.Ind, 110, 115)
Video.Newest <- max(as.numeric(Video.Dates))
Video.Use <- Video.Ind[substring(Video.Ind, 110, 115) %in% Video.Newest]
VideoFile <- data.frame(read.csv(Video.Use, header=TRUE))

Pec.Af.Kept <- subset(PecKept, substring(PecKept$File.name, 1,2)=="af")
Pec.Af.Kept.VideoInfo <- VideoFile[which(VideoFile$File.name %in% Pec.Af.Kept$File.name),]

Pec.Pb.Kept <- subset(PecKept, substring(PecKept$File.name, 1,2)=="pb")
Pec.Pb.Kept.VideoInfo <- VideoFile[which(VideoFile$File.name %in% Pec.Pb.Kept$File.name),]

Pec.Pw.Kept <- subset(PecKept, substring(PecKept$File.name, 1,2)=="pw")
Pec.Pw.Kept.VideoInfo <- VideoFile[which(VideoFile$File.name %in% Pec.Pw.Kept$File.name),]

Pel.Af.Kept <- subset(PelKept, substring(PelKept$File.name, 1,2)=="af")
Pel.Af.Kept.VideoInfo <- VideoFile[which(VideoFile$File.name %in% Pel.Af.Kept$File.name),]

Pel.Pw.Kept <- subset(PelKept, substring(PelKept$File.name, 1,2)=="pw")
Pel.Pw.Kept.VideoInfo <- VideoFile[which(VideoFile$File.name %in% Pel.Pw.Kept$File.name),]

Pec.Sl.Kept <- subset(PecKept, substring(PecKept$File.name, 1,2)=="sl")
Pec.Sl.Kept.VideoInfo <- VideoFile[which(VideoFile$File.name %in% Pec.Sl.Kept$File.name),]


####### ASSESSING DIFFERENCES IN STANCE DURATIONS ######
#####  THIS IS THE ACTUAL METHOD THAT WE USE FOR OUR ANALYSES
##  Evaluating significance
library(lme4)
library(sciplot) # for calculating standard errors
library(boot) # for bootstrapping
library(lsmeans) #for Tukey post hoc comparison
library(pbkrtest)
# use the p-values in pMCMC (posterior distribution from a Markov Chain Monte Carlo)
# The intercept is the baseline (usually the group name starting with a letter earlier in the alphabet)
# The way you interpret the output is the the pMCMC value will tell you whether that group (e.g., pw)
# is significantly different from the baseline or intercept


###  The MCMC sampling from before isn't acceptable anymore, so going to use alternative approaches
##  Bolker et al. (2008) described how likelihood ratio tests are not recommended for GLMMs
## Although I'm not using the results from the LRT, I'm running it just to see the results
## Tukey post hoc comparisons on the least squares means are then being use to look at pairwise comparisons
## Although REML is less biased, you need to have Maximum Likelihood to calculate your AIC and BIC values and to conduct the LRT

Pec.Pb.Kept.VideoInfo$Group <- "Pec.Pb"
Pec.Af.Kept.VideoInfo$Group <- "Pec.Af"
Pel.Af.Kept.VideoInfo$Group <- "Pel.Af"
Pec.Pw.Kept.VideoInfo$Group <- "Pec.Pw"
Pel.Pw.Kept.VideoInfo$Group <- "Pel.Pw"

Pec.Pb.Kept.VideoInfo$Ind <- substring(Pec.Pb.Kept.VideoInfo$File.name, 1,4)
Pec.Af.Kept.VideoInfo$Ind <- substring(Pec.Af.Kept.VideoInfo$File.name, 1,4)
Pel.Af.Kept.VideoInfo$Ind <- substring(Pel.Af.Kept.VideoInfo$File.name, 1,4)
Pec.Pw.Kept.VideoInfo$Ind <- substring(Pec.Pw.Kept.VideoInfo$File.name, 1,4)
Pel.Pw.Kept.VideoInfo$Ind <- substring(Pel.Pw.Kept.VideoInfo$File.name, 1,4)


StanceDurationData <- rbind(Pec.Pb.Kept.VideoInfo, Pec.Af.Kept.VideoInfo, Pel.Af.Kept.VideoInfo, Pec.Pw.Kept.VideoInfo, Pel.Pw.Kept.VideoInfo) 

for (i in 1:nrow(StanceDurationData)) {
  ifelse(StanceDurationData$Appendages[[i]]==unique("Pelvic"), StanceDurationData$Stance[[i]] <- StanceDurationData$Pelvic.End.Frame[[i]]-StanceDurationData$Pelvic.Start.Frame[[i]],
         StanceDurationData$Stance[[i]] <- StanceDurationData$Pectoral.End.Frame[[i]]-StanceDurationData$Pectoral.Start.Frame[[i]])
         
}

Stance.LMER <- lmer(Stance~Group+(1|Ind), StanceDurationData, verbose=T, REML=F)
Stance.LMER.Null <- lmer(Stance~1+(1|Ind), StanceDurationData, REML=F)
# Likelihood ratio tests cannot take REML products; it will autotmatically convert it to Maximum Likehood before analysis
Stance.LMER.Compare <- anova(Stance.LMER, Stance.LMER.Null)
Stance.Tukey <- lsmeans(Stance.LMER, pairwise~as.factor(Group), adjust="tukey")
Stance.Means <- aggregate(StanceDurationData$Stance, by=list(StanceDurationData$Group), FUN=mean)[2]
Stance.SE <- aggregate(StanceDurationData$Stance, by=list(StanceDurationData$Group), FUN=se)[2]
StanceSummaryStats <- data.frame(unique(StanceDurationData$Group), Stance.Means, Stance.SE)
names(StanceSummaryStats) <- c("group", "mean", "se")



# Organizing data
# Categorizing by individual and species
PecKept$Ind <- substring(PecKept$File.name, 1, 4)
PecKept$Species <- substring(PecKept$File.name, 1, 2)      # note: af refers to Ambystoma tigrinum forelimb, so instead of at, using af, to differentiate from Megan's dataset which used "at"

PelKept$Ind <- substring(PelKept$File.name, 1, 4)
PelKept$Species <- substring(PelKept$File.name, 1, 2)

# Subsetting species
PecKept.Af <- PecKept[grep('af', PecKept$Species),];  PecKept.Pb <- PecKept[grep('pb', PecKept$Species),]; PecKept.Pw <- PecKept[grep('pw', PecKept$Species),]
PelKept.Af <- PelKept[grep('af', PelKept$Species),];  PelKept.Pw <- PelKept[grep('pw', PelKept$Species),]; PecKept.Sl <- PecKept[grep('sl', PecKept$Species),]


# Calculating means for all variables for those trials that were kept
PeakNetGRF.Pec.Af.Means <- sapply(PecKept.Af[,c(1:16)], FUN=function(x) mean(x))
PeakNetGRF.Pel.Af.Means <- sapply(PelKept.Af[,c(1:16)], FUN=function(x) mean(x))
PeakNetGRF.Pec.Pb.Means <- sapply(PecKept.Pb[,c(1:16)], FUN=function(x) mean(x))
PeakNetGRF.Pec.Pw.Means <- sapply(PecKept.Pw[,c(1:16)], FUN=function(x) mean(x))
PeakNetGRF.Pel.Pw.Means <- sapply(PelKept.Pw[,c(1:16)], FUN=function(x) mean(x))
PeakNetGRF.Pec.Sl.Means <- sapply(PecKept.Sl[,c(1:16)], FUN=function(x) mean(x))

# Standard errors
library(sciplot)
PeakNetGRF.Pec.Af.SE  <- sapply(PecKept.Af[,c(1:16)], FUN=function(x) se(x, na.rm=TRUE))
PeakNetGRF.Pel.Af.SE  <- sapply(PelKept.Af[,c(1:16)], FUN=function(x) se(x, na.rm=TRUE))
PeakNetGRF.Pec.Pb.SE  <- sapply(PecKept.Pb[,c(1:16)], FUN=function(x) se(x, na.rm=TRUE))
PeakNetGRF.Pec.Pw.SE  <- sapply(PecKept.Pw[,c(1:16)], FUN=function(x) se(x, na.rm=TRUE))
PeakNetGRF.Pel.Pw.SE  <- sapply(PelKept.Pw[,c(1:16)], FUN=function(x) se(x, na.rm=TRUE))
PeakNetGRF.Pec.Sl.SE  <- sapply(PecKept.Sl[,c(1:16)], FUN=function(x) se(x, na.rm=TRUE))

# Putting together the summary info
PeakNetGRF.Pec.Af.Sum <- rbind(PeakNetGRF.Pec.Af.Means, PeakNetGRF.Pec.Af.SE); PeakNetGRF.Pec.Af.Sum <- cbind(rbind("Mean", "Standard Error"), PeakNetGRF.Pec.Af.Sum)
PeakNetGRF.Pel.Af.Sum <- rbind(PeakNetGRF.Pel.Af.Means, PeakNetGRF.Pel.Af.SE); PeakNetGRF.Pel.Af.Sum <- cbind(rbind("Mean", "Standard Error"), PeakNetGRF.Pel.Af.Sum)
PeakNetGRF.Pec.Pb.Sum <- rbind(PeakNetGRF.Pec.Pb.Means, PeakNetGRF.Pec.Pb.SE); PeakNetGRF.Pec.Pb.Sum <- cbind(rbind("Mean", "Standard Error"), PeakNetGRF.Pec.Pb.Sum)
PeakNetGRF.Pec.Pw.Sum <- rbind(PeakNetGRF.Pec.Pw.Means, PeakNetGRF.Pec.Pw.SE); PeakNetGRF.Pec.Pw.Sum <- cbind(rbind("Mean", "Standard Error"), PeakNetGRF.Pec.Pw.Sum)
PeakNetGRF.Pel.Pw.Sum <- rbind(PeakNetGRF.Pel.Pw.Means, PeakNetGRF.Pel.Pw.SE); PeakNetGRF.Pel.Pw.Sum <- cbind(rbind("Mean", "Standard Error"), PeakNetGRF.Pel.Pw.Sum)
PeakNetGRF.Pec.Sl.Sum <- rbind(PeakNetGRF.Pec.Sl.Means, PeakNetGRF.Pec.Sl.SE); PeakNetGRF.Pec.Sl.Sum <- cbind(rbind("Mean", "Standard Error"), PeakNetGRF.Pec.Sl.Sum)


# Saving data
setwd('/Users/SandyMKawano/Desktop/Research/Forelimbs on Force Plates/Force Data/R Analysis/Step 4 Results/Peak Net GRF')

SavePecAfKeep <- paste("Pec_Af_PeakNetGRF_", SaveDate, ".csv", sep="")
SavePelAfKeep <- paste("Pel_Af_PeakNetGRF_", SaveDate, ".csv", sep="")
SavePecPbKeep <- paste("Pec_Pb_PeakNetGRF_", SaveDate, ".csv", sep="")
SavePecPwKeep <- paste("Pec_Pw_PeakNetGRF_", SaveDate, ".csv", sep="")
SavePelPwKeep <- paste("Pel_Pw_PeakNetGRF_", SaveDate, ".csv", sep="")
SavePecSlKeep <- paste("Pec_Sl_PeakNetGRF_", SaveDate, ".csv", sep="")

write.table(PecKept.Af, file=SavePecAfKeep, sep =",", row.names=FALSE)
write.table(PelKept.Af, file=SavePelAfKeep, sep =",", row.names=FALSE)
write.table(PecKept.Pb, file=SavePecPbKeep, sep =",", row.names=FALSE)
write.table(PecKept.Pw, file=SavePecPwKeep, sep =",", row.names=FALSE)
write.table(PelKept.Pw, file=SavePelPwKeep, sep =",", row.names=FALSE)
write.table(PecKept.Sl, file=SavePecSlKeep, sep =",", row.names=FALSE)

SavePecAfKeepSum <- paste("Pec_Af_PeakNetGRFSum_", SaveDate, ".csv", sep="")
SavePelAfKeepSum <- paste("Pel_Af_PeakNetGRFSum_", SaveDate, ".csv", sep="")
SavePecPbKeepSum <- paste("Pec_Pb_PeakNetGRFSum_", SaveDate, ".csv", sep="")
SavePecPwKeepSum <- paste("Pec_Pw_PeakNetGRFSum_", SaveDate, ".csv", sep="")
SavePelPwKeepSum <- paste("Pel_Pw_PeakNetGRFSum_", SaveDate, ".csv", sep="")
SavePecSlKeepSum <- paste("Pec_Sl_PeakNetGRFSum_", SaveDate, ".csv", sep="")

write.table(PeakNetGRF.Pec.Af.Sum, file=SavePecAfKeepSum, sep =",", row.names=FALSE)
write.table(PeakNetGRF.Pel.Af.Sum, file=SavePelAfKeepSum, sep =",", row.names=FALSE)
write.table(PeakNetGRF.Pec.Pb.Sum, file=SavePecPbKeepSum, sep =",", row.names=FALSE)
write.table(PeakNetGRF.Pec.Pw.Sum, file=SavePecPwKeepSum, sep =",", row.names=FALSE)
write.table(PeakNetGRF.Pel.Pw.Sum, file=SavePelPwKeepSum, sep =",", row.names=FALSE)
write.table(PeakNetGRF.Pec.Sl.Sum, file=SavePecSlKeepSum, sep =",", row.names=FALSE)


######################  GRF DATA AT INCREMENTS #################################
setwd('/Users/SandyMKawano/Desktop/Research/Forelimbs on Force Plates/Force Data/R Analysis/Step 3 Post Filter Processing/Filtered at Increments')
GRFPath <- setwd('/Users/SandyMKawano/Desktop/Research/Forelimbs on Force Plates/Force Data/R Analysis/Step 3 Post Filter Processing/Filtered at Increments')
GRF.Ind <- list.files(GRFPath, pattern=".csv", full=TRUE)   # Only want .csv files

File.name <- substring(GRF.Ind, 139, 145)
Appendage <- substring(GRF.Ind, 147, 149)
Species <- substring(GRF.Ind, 139,140)

PercentStance <- seq(0,100,5)

# Loading all files
GRF.All <- lapply(GRF.Ind, read.csv, header=TRUE)
#names(GRF.All) <- File.name   # want to be able to subset by appendage type


# Combining all the data into a single dataframe
# Need to create empty matrix to input data from the loop into
#GRF.APAngle <- rep(NA, length(GRF.All))
#for (i in 1:length(GRF.All)) {GRF.APAngle[i] <- GRF.All[[i]][2]}
#GRF.APAngle.Combined <- as.data.frame(do.call("rbind", GRF.APAngle))
#names(GRF.APAngle.Combined) <- PercentStance
#GRF.APAngle.Combined$File.name <- File.name
#GRF.APAngle.Combined$Appendage <- Appendage

# If do.call produces a smaller number of rows than in GRF.All, look in GRF.APAngle.Convert, and see if any are illogical

## For converted AP Angle
GRF.APAngle.Convert <- rep(NA, length(GRF.All))
for (i in 1:length(GRF.All)) {GRF.APAngle.Convert[i] <- GRF.All[[i]][2]}
GRF.APAngle.Convert.Combined <- as.data.frame(do.call("rbind", GRF.APAngle.Convert))
names(GRF.APAngle.Convert.Combined) <- PercentStance
GRF.APAngle.Convert.Combined$File.name <- File.name
GRF.APAngle.Convert.Combined$Appendage <- Appendage
GRF.APAngle.Convert.Combined$Species <- Species


GRF.MLAngle.Convert <- rep(NA, length(GRF.All))
for (i in 1:length(GRF.All)) {GRF.MLAngle.Convert[i] <- GRF.All[[i]][3]}
GRF.MLAngle.Convert.Combined <- as.data.frame(do.call("rbind", GRF.MLAngle.Convert))
names(GRF.MLAngle.Convert.Combined) <- PercentStance
GRF.MLAngle.Convert.Combined$File.name <- File.name
GRF.MLAngle.Convert.Combined$Appendage <- Appendage
GRF.MLAngle.Convert.Combined$Species <- Species

GRF.InterpV.BW <- rep(NA, length(GRF.All))
for (i in 1:length(GRF.All)) {GRF.InterpV.BW[i] <- GRF.All[[i]][4]}
GRF.InterpV.BW.Combined <- as.data.frame(do.call("rbind", GRF.InterpV.BW))
names(GRF.InterpV.BW.Combined) <- PercentStance
GRF.InterpV.BW.Combined$File.name <- File.name
GRF.InterpV.BW.Combined$Appendage <- Appendage
GRF.InterpV.BW.Combined$Species <- Species
  
GRF.InterpML.BW <- rep(NA, length(GRF.All))
for (i in 1:length(GRF.All)) {GRF.InterpML.BW[i] <- GRF.All[[i]][5]}
GRF.InterpML.BW.Combined <- as.data.frame(do.call("rbind", GRF.InterpML.BW))
names(GRF.InterpML.BW.Combined) <- PercentStance
GRF.InterpML.BW.Combined$File.name <- File.name
GRF.InterpML.BW.Combined$Appendage <- Appendage
GRF.InterpML.BW.Combined$Species <- Species

GRF.InterpHz.BW <- rep(NA, length(GRF.All))
for (i in 1:length(GRF.All)) {GRF.InterpHz.BW[i] <- GRF.All[[i]][6]}
GRF.InterpHz.BW.Combined <- as.data.frame(do.call("rbind", GRF.InterpHz.BW))
names(GRF.InterpHz.BW.Combined) <- PercentStance
GRF.InterpHz.BW.Combined$File.name <- File.name
GRF.InterpHz.BW.Combined$Appendage <- Appendage
GRF.InterpHz.BW.Combined$Species <- Species

GRF.NetGRF.BW <- rep(NA, length(GRF.All))
for (i in 1:length(GRF.All)) {GRF.NetGRF.BW[i] <- GRF.All[[i]][7]}
GRF.NetGRF.BW.Combined <- as.data.frame(do.call("rbind", GRF.NetGRF.BW))
names(GRF.NetGRF.BW.Combined) <- PercentStance
GRF.NetGRF.BW.Combined$File.name <- File.name
GRF.NetGRF.BW.Combined$Appendage <- Appendage
GRF.NetGRF.BW.Combined$Species <- Species


# Subsetting into Pectoral and Pelvic files
#GRF.APAngle.Combined.Pec <- GRF.APAngle.Combined[which(GRF.APAngle.Combined$Appendage=="Pec"),]
GRF.APAngle.Convert.Combined.Pec <- GRF.APAngle.Convert.Combined[which(GRF.APAngle.Convert.Combined$Appendage=="Pec"),]
GRF.MLAngle.Convert.Combined.Pec <- GRF.MLAngle.Convert.Combined[which(GRF.MLAngle.Convert.Combined$Appendage=="Pec"),]
GRF.InterpV.BW.Combined.Pec <- GRF.InterpV.BW.Combined[which(GRF.InterpV.BW.Combined$Appendage=="Pec"),]
GRF.InterpML.BW.Combined.Pec <- GRF.InterpML.BW.Combined[which(GRF.InterpML.BW.Combined$Appendage=="Pec"),]
GRF.InterpHz.BW.Combined.Pec <- GRF.InterpHz.BW.Combined[which(GRF.InterpHz.BW.Combined$Appendage=="Pec"),]
GRF.NetGRF.BW.Combined.Pec <- GRF.NetGRF.BW.Combined[which(GRF.NetGRF.BW.Combined$Appendage=="Pec"),]

#GRF.APAngle.Combined.Pel <- GRF.APAngle.Combined[which(GRF.APAngle.Combined$Appendage=="Pel"),]
GRF.APAngle.Convert.Combined.Pel <- GRF.APAngle.Convert.Combined[which(GRF.APAngle.Convert.Combined$Appendage=="Pel"),]
GRF.MLAngle.Convert.Combined.Pel <- GRF.MLAngle.Convert.Combined[which(GRF.MLAngle.Convert.Combined$Appendage=="Pel"),]
GRF.InterpV.BW.Combined.Pel <- GRF.InterpV.BW.Combined[which(GRF.InterpV.BW.Combined$Appendage=="Pel"),]
GRF.InterpML.BW.Combined.Pel <- GRF.InterpML.BW.Combined[which(GRF.InterpML.BW.Combined$Appendage=="Pel"),]
GRF.InterpHz.BW.Combined.Pel <- GRF.InterpHz.BW.Combined[which(GRF.InterpHz.BW.Combined$Appendage=="Pel"),]
GRF.NetGRF.BW.Combined.Pel <- GRF.NetGRF.BW.Combined[which(GRF.NetGRF.BW.Combined$Appendage=="Pel"),]


# Subsetting only files I decided to keep
#GRF.APAngle.Combined.Pec.Keep <- GRF.APAngle.Combined.Pec[which(GRF.APAngle.Combined.Pec$File.name %in% PecFilesKept$File.name),]
GRF.APAngle.Convert.Combined.Pec.Keep <- GRF.APAngle.Convert.Combined.Pec[which(GRF.APAngle.Convert.Combined.Pec$File.name %in% PecFilesKept$File.name),]
GRF.MLAngle.Convert.Combined.Pec.Keep <- GRF.MLAngle.Convert.Combined.Pec[which(GRF.MLAngle.Convert.Combined.Pec$File.name %in% PecFilesKept$File.name),]
GRF.InterpV.BW.Combined.Pec.Keep <- GRF.InterpV.BW.Combined.Pec[which(GRF.InterpV.BW.Combined.Pec$File.name %in% PecFilesKept$File.name),]
GRF.InterpML.BW.Combined.Pec.Keep <- GRF.InterpML.BW.Combined.Pec[which(GRF.InterpML.BW.Combined.Pec$File.name %in% PecFilesKept$File.name),]
GRF.InterpHz.BW.Combined.Pec.Keep <- GRF.InterpHz.BW.Combined.Pec[which(GRF.InterpHz.BW.Combined.Pec$File.name %in% PecFilesKept$File.name),]
GRF.NetGRF.BW.Combined.Pec.Keep <- GRF.NetGRF.BW.Combined.Pec[which(GRF.NetGRF.BW.Combined.Pec$File.name %in% PecFilesKept$File.name),]


#GRF.APAngle.Combined.Pel.Keep <- GRF.APAngle.Combined.Pel[which(GRF.APAngle.Combined.Pel$File.name %in% PelFilesKept$File.name),]
GRF.APAngle.Convert.Combined.Pel.Keep <- GRF.APAngle.Convert.Combined.Pel[which(GRF.APAngle.Convert.Combined.Pel$File.name %in% PelFilesKept$File.name),]
GRF.MLAngle.Convert.Combined.Pel.Keep <- GRF.MLAngle.Convert.Combined.Pel[which(GRF.MLAngle.Convert.Combined.Pel$File.name %in% PelFilesKept$File.name),]
GRF.InterpV.BW.Combined.Pel.Keep <- GRF.InterpV.BW.Combined.Pel[which(GRF.InterpV.BW.Combined.Pel$File.name %in% PelFilesKept$File.name),]
GRF.InterpML.BW.Combined.Pel.Keep <- GRF.InterpML.BW.Combined.Pel[which(GRF.InterpML.BW.Combined.Pel$File.name %in% PelFilesKept$File.name),]
GRF.InterpHz.BW.Combined.Pel.Keep <- GRF.InterpHz.BW.Combined.Pel[which(GRF.InterpHz.BW.Combined.Pel$File.name %in% PelFilesKept$File.name),]
GRF.NetGRF.BW.Combined.Pel.Keep <- GRF.NetGRF.BW.Combined.Pel[which(GRF.NetGRF.BW.Combined.Pel$File.name %in% PelFilesKept$File.name),]


# Dividing into Ambystoma, Periophthalmus, Pleurodeles, and Siren files
#Af.GRF.APAngle.Combined.Pec <- GRF.APAngle.Combined.Pec.Keep[which(substring(GRF.APAngle.Combined.Pec.Keep$File.name, 1,2)=="af"),]
Af.GRF.APAngle.Convert.Combined.Pec <- GRF.APAngle.Convert.Combined.Pec.Keep[which(substring(GRF.APAngle.Convert.Combined.Pec.Keep$File.name, 1,2)=="af"),]
Af.GRF.MLAngle.Convert.Combined.Pec <- GRF.MLAngle.Convert.Combined.Pec.Keep[which(substring(GRF.MLAngle.Convert.Combined.Pec.Keep$File.name, 1,2)=="af"),]
Af.GRF.InterpV.BW.Combined.Pec <- GRF.InterpV.BW.Combined.Pec.Keep[which(substring(GRF.InterpV.BW.Combined.Pec.Keep$File.name, 1,2)=="af"),]
Af.GRF.InterpML.BW.Combined.Pec <- GRF.InterpML.BW.Combined.Pec.Keep[which(substring(GRF.InterpML.BW.Combined.Pec.Keep$File.name, 1,2)=="af"),]
Af.GRF.InterpHz.BW.Combined.Pec <- GRF.InterpHz.BW.Combined.Pec.Keep[which(substring(GRF.InterpHz.BW.Combined.Pec.Keep$File.name, 1,2)=="af"),]
Af.GRF.NetGRF.BW.Combined.Pec <- GRF.NetGRF.BW.Combined.Pec.Keep[which(substring(GRF.NetGRF.BW.Combined.Pec.Keep$File.name, 1,2)=="af"),]

#Af.GRF.APAngle.Combined.Pel <- GRF.APAngle.Combined.Pel.Keep[which(substring(GRF.APAngle.Combined.Pel.Keep$File.name, 1,2)=="af"),]
Af.GRF.APAngle.Convert.Combined.Pel <- GRF.APAngle.Convert.Combined.Pel.Keep[which(substring(GRF.APAngle.Convert.Combined.Pel.Keep$File.name, 1,2)=="af"),]
Af.GRF.MLAngle.Convert.Combined.Pel <- GRF.MLAngle.Convert.Combined.Pel.Keep[which(substring(GRF.MLAngle.Convert.Combined.Pel.Keep$File.name, 1,2)=="af"),]
Af.GRF.InterpV.BW.Combined.Pel <- GRF.InterpV.BW.Combined.Pel.Keep[which(substring(GRF.InterpV.BW.Combined.Pel.Keep$File.name, 1,2)=="af"),]
Af.GRF.InterpML.BW.Combined.Pel <- GRF.InterpML.BW.Combined.Pel.Keep[which(substring(GRF.InterpML.BW.Combined.Pel.Keep$File.name, 1,2)=="af"),]
Af.GRF.InterpHz.BW.Combined.Pel <- GRF.InterpHz.BW.Combined.Pel.Keep[which(substring(GRF.InterpHz.BW.Combined.Pel.Keep$File.name, 1,2)=="af"),]
Af.GRF.NetGRF.BW.Combined.Pel <- GRF.NetGRF.BW.Combined.Pel.Keep[which(substring(GRF.NetGRF.BW.Combined.Pel.Keep$File.name, 1,2)=="af"),]

#Pb.GRF.APAngle.Combined.Pec <- GRF.APAngle.Combined.Pec.Keep[which(substring(GRF.APAngle.Combined.Pec.Keep$File.name, 1,2)=="pb"),]
Pb.GRF.APAngle.Convert.Combined.Pec <- GRF.APAngle.Convert.Combined.Pec.Keep[which(substring(GRF.APAngle.Convert.Combined.Pec.Keep$File.name, 1,2)=="pb"),]
Pb.GRF.MLAngle.Convert.Combined.Pec <- GRF.MLAngle.Convert.Combined.Pec.Keep[which(substring(GRF.MLAngle.Convert.Combined.Pec.Keep$File.name, 1,2)=="pb"),]
Pb.GRF.InterpV.BW.Combined.Pec <- GRF.InterpV.BW.Combined.Pec.Keep[which(substring(GRF.InterpV.BW.Combined.Pec.Keep$File.name, 1,2)=="pb"),]
Pb.GRF.InterpML.BW.Combined.Pec <- GRF.InterpML.BW.Combined.Pec.Keep[which(substring(GRF.InterpML.BW.Combined.Pec.Keep$File.name, 1,2)=="pb"),]
Pb.GRF.InterpHz.BW.Combined.Pec <- GRF.InterpHz.BW.Combined.Pec.Keep[which(substring(GRF.InterpHz.BW.Combined.Pec.Keep$File.name, 1,2)=="pb"),]
Pb.GRF.NetGRF.BW.Combined.Pec <- GRF.NetGRF.BW.Combined.Pec.Keep[which(substring(GRF.NetGRF.BW.Combined.Pec.Keep$File.name, 1,2)=="pb"),]

#Pw.GRF.APAngle.Combined.Pec <- GRF.APAngle.Combined.Pec.Keep[which(substring(GRF.APAngle.Combined.Pec.Keep$File.name, 1,2)=="pw"),]
Pw.GRF.APAngle.Convert.Combined.Pec <- GRF.APAngle.Convert.Combined.Pec.Keep[which(substring(GRF.APAngle.Convert.Combined.Pec.Keep$File.name, 1,2)=="pw"),]
Pw.GRF.MLAngle.Convert.Combined.Pec <- GRF.MLAngle.Convert.Combined.Pec.Keep[which(substring(GRF.MLAngle.Convert.Combined.Pec.Keep$File.name, 1,2)=="pw"),]
Pw.GRF.InterpV.BW.Combined.Pec <- GRF.InterpV.BW.Combined.Pec.Keep[which(substring(GRF.InterpV.BW.Combined.Pec.Keep$File.name, 1,2)=="pw"),]
Pw.GRF.InterpML.BW.Combined.Pec <- GRF.InterpML.BW.Combined.Pec.Keep[which(substring(GRF.InterpML.BW.Combined.Pec.Keep$File.name, 1,2)=="pw"),]
Pw.GRF.InterpHz.BW.Combined.Pec <- GRF.InterpHz.BW.Combined.Pec.Keep[which(substring(GRF.InterpHz.BW.Combined.Pec.Keep$File.name, 1,2)=="pw"),]
Pw.GRF.NetGRF.BW.Combined.Pec <- GRF.NetGRF.BW.Combined.Pec.Keep[which(substring(GRF.NetGRF.BW.Combined.Pec.Keep$File.name, 1,2)=="pw"),]

#Pw.GRF.APAngle.Combined.Pel <- GRF.APAngle.Combined.Pel.Keep[which(substring(GRF.APAngle.Combined.Pel.Keep$File.name, 1,2)=="pw"),]
Pw.GRF.APAngle.Convert.Combined.Pel <- GRF.APAngle.Convert.Combined.Pel.Keep[which(substring(GRF.APAngle.Convert.Combined.Pel.Keep$File.name, 1,2)=="pw"),]
Pw.GRF.MLAngle.Convert.Combined.Pel <- GRF.MLAngle.Convert.Combined.Pel.Keep[which(substring(GRF.MLAngle.Convert.Combined.Pel.Keep$File.name, 1,2)=="pw"),]
Pw.GRF.InterpV.BW.Combined.Pel <- GRF.InterpV.BW.Combined.Pel.Keep[which(substring(GRF.InterpV.BW.Combined.Pel.Keep$File.name, 1,2)=="pw"),]
Pw.GRF.InterpML.BW.Combined.Pel <- GRF.InterpML.BW.Combined.Pel.Keep[which(substring(GRF.InterpML.BW.Combined.Pel.Keep$File.name, 1,2)=="pw"),]
Pw.GRF.InterpHz.BW.Combined.Pel <- GRF.InterpHz.BW.Combined.Pel.Keep[which(substring(GRF.InterpHz.BW.Combined.Pel.Keep$File.name, 1,2)=="pw"),]
Pw.GRF.NetGRF.BW.Combined.Pel <- GRF.NetGRF.BW.Combined.Pel.Keep[which(substring(GRF.NetGRF.BW.Combined.Pel.Keep$File.name, 1,2)=="pw"),]

#Sl.GRF.APAngle.Combined.Pec <- GRF.APAngle.Combined.Pec.Keep[which(substring(GRF.APAngle.Combined.Pec.Keep$File.name, 1,2)=="sl"),]
Sl.GRF.APAngle.Convert.Combined.Pec <- GRF.APAngle.Convert.Combined.Pec.Keep[which(substring(GRF.APAngle.Convert.Combined.Pec.Keep$File.name, 1,2)=="sl"),]
Sl.GRF.MLAngle.Convert.Combined.Pec <- GRF.MLAngle.Convert.Combined.Pec.Keep[which(substring(GRF.MLAngle.Convert.Combined.Pec.Keep$File.name, 1,2)=="sl"),]
Sl.GRF.InterpV.BW.Combined.Pec <- GRF.InterpV.BW.Combined.Pec.Keep[which(substring(GRF.InterpV.BW.Combined.Pec.Keep$File.name, 1,2)=="sl"),]
Sl.GRF.InterpML.BW.Combined.Pec <- GRF.InterpML.BW.Combined.Pec.Keep[which(substring(GRF.InterpML.BW.Combined.Pec.Keep$File.name, 1,2)=="sl"),]
Sl.GRF.InterpHz.BW.Combined.Pec <- GRF.InterpHz.BW.Combined.Pec.Keep[which(substring(GRF.InterpHz.BW.Combined.Pec.Keep$File.name, 1,2)=="sl"),]
Sl.GRF.NetGRF.BW.Combined.Pec <- GRF.NetGRF.BW.Combined.Pec.Keep[which(substring(GRF.NetGRF.BW.Combined.Pec.Keep$File.name, 1,2)=="sl"),]


setwd('/Users/SandyMKawano/Desktop/Research/Forelimbs on Force Plates/Force Data/R Analysis/Step 4 Results/Angle AP')
#Save.Af.APAngle.Pec <- paste ("Pec_Af_APAngle_", SaveDate, ".csv", sep="")
#write.table(Af.GRF.APAngle.Combined.Pec, file=Save.Af.APAngle.Pec, sep =",", row.names=FALSE)
#Save.Pb.APAngle.Pec <- paste ("Pec_Pb_APAngle_", SaveDate, ".csv", sep="")
#write.table(Pb.GRF.APAngle.Combined.Pec, file=Save.Pb.APAngle.Pec, sep =",", row.names=FALSE)
#Save.Af.APAngle.Pel <- paste ("Pel_Af_APAngle_", SaveDate, ".csv", sep="")
#write.table(Af.GRF.APAngle.Combined.Pel, file=Save.Af.APAngle.Pel, sep =",", row.names=FALSE)

Save.Af.APAngle.Pec <- paste ("Pec_Af_APAngle_", SaveDate, ".csv", sep="")
write.table(Af.GRF.APAngle.Convert.Combined.Pec, file=Save.Af.APAngle.Pec, sep =",", row.names=FALSE)
Save.Pb.APAngle.Pec <- paste ("Pec_Pb_APAngle_", SaveDate, ".csv", sep="")
write.table(Pb.GRF.APAngle.Convert.Combined.Pec, file=Save.Pb.APAngle.Pec, sep =",", row.names=FALSE)
Save.Af.APAngle.Pel <- paste ("Pel_Af_APAngle_", SaveDate, ".csv", sep="")
write.table(Af.GRF.APAngle.Convert.Combined.Pel, file=Save.Af.APAngle.Pel, sep =",", row.names=FALSE)
Save.Pw.APAngle.Pec <- paste ("Pec_Pw_APAngle_", SaveDate, ".csv", sep="")
write.table(Pw.GRF.APAngle.Convert.Combined.Pec, file=Save.Pw.APAngle.Pec, sep =",", row.names=FALSE)
Save.Pw.APAngle.Pel <- paste ("Pel_Pw_APAngle_", SaveDate, ".csv", sep="")
write.table(Pw.GRF.APAngle.Convert.Combined.Pel, file=Save.Pw.APAngle.Pel, sep =",", row.names=FALSE)
#

setwd('/Users/SandyMKawano/Desktop/Research/Forelimbs on Force Plates/Force Data/R Analysis/Step 4 Results/Angle ML')
Save.Af.MLAngle.Pec <- paste ("Pec_Af_MLAngle_", SaveDate, ".csv", sep="")
write.table(Af.GRF.MLAngle.Convert.Combined.Pec, file=Save.Af.MLAngle.Pec, sep =",", row.names=FALSE)
Save.Pb.MLAngle.Pec <- paste ("Pec_Pb_MLAngle_", SaveDate, ".csv", sep="")
write.table(Pb.GRF.MLAngle.Convert.Combined.Pec, file=Save.Pb.MLAngle.Pec, sep =",", row.names=FALSE)
Save.Af.MLAngle.Pel <- paste ("Pel_Af_MLAngle_", SaveDate, ".csv", sep="")
write.table(Af.GRF.MLAngle.Convert.Combined.Pel, file=Save.Af.MLAngle.Pel, sep =",", row.names=FALSE)
Save.Pw.MLAngle.Pec <- paste ("Pec_Pw_MLAngle_", SaveDate, ".csv", sep="")
write.table(Pw.GRF.MLAngle.Convert.Combined.Pec, file=Save.Pw.MLAngle.Pec, sep =",", row.names=FALSE)
Save.Pw.MLAngle.Pel <- paste ("Pel_Pw_MLAngle_", SaveDate, ".csv", sep="")
write.table(Pw.GRF.MLAngle.Convert.Combined.Pel, file=Save.Pw.MLAngle.Pel, sep =",", row.names=FALSE)


setwd('/Users/SandyMKawano/Desktop/Research/Forelimbs on Force Plates/Force Data/R Analysis/Step 4 Results/GRF V')
Save.Af.GRFV.BW.Pec <- paste ("Pec_Af_GRFV_BW_", SaveDate, ".csv", sep="")
write.table(Af.GRF.InterpV.BW.Combined.Pec, file=Save.Af.GRFV.BW.Pec, sep =",", row.names=FALSE)
Save.Pb.GRFV.BW.Pec <- paste ("Pec_Pb_GRFV_BW_", SaveDate, ".csv", sep="")
write.table(Pb.GRF.InterpV.BW.Combined.Pec, file=Save.Pb.GRFV.BW.Pec, sep =",", row.names=FALSE)
Save.Af.GRFV.BW.Pel <- paste ("Pel_Af_GRFV_BW_", SaveDate, ".csv", sep="")
write.table(Af.GRF.InterpV.BW.Combined.Pel, file=Save.Af.GRFV.BW.Pel, sep =",", row.names=FALSE)
Save.Pw.GRFV.BW.Pec <- paste ("Pec_Pw_GRFV_BW_", SaveDate, ".csv", sep="")
write.table(Pw.GRF.InterpV.BW.Combined.Pec, file=Save.Pw.GRFV.BW.Pec, sep =",", row.names=FALSE)
Save.Pw.GRFV.BW.Pel <- paste ("Pel_Pw_GRFV_BW_", SaveDate, ".csv", sep="")
write.table(Pw.GRF.InterpV.BW.Combined.Pel, file=Save.Pw.GRFV.BW.Pel, sep =",", row.names=FALSE)


setwd('/Users/SandyMKawano/Desktop/Research/Forelimbs on Force Plates/Force Data/R Analysis/Step 4 Results/GRF ML')
Save.Af.GRFML.BW.Pec <- paste ("Pec_Af_GRFML_BW_", SaveDate, ".csv", sep="")
write.table(Af.GRF.InterpML.BW.Combined.Pec, file=Save.Af.GRFML.BW.Pec, sep =",", row.names=FALSE)
Save.Pb.GRFML.BW.Pec <- paste ("Pec_Pb_GRFML_BW_", SaveDate, ".csv", sep="")
write.table(Pb.GRF.InterpML.BW.Combined.Pec, file=Save.Pb.GRFML.BW.Pec, sep =",", row.names=FALSE)
Save.Af.GRFML.BW.Pel <- paste ("Pel_Af_GRFML_BW_", SaveDate, ".csv", sep="")
write.table(Af.GRF.InterpML.BW.Combined.Pel, file=Save.Af.GRFML.BW.Pel, sep =",", row.names=FALSE)
Save.Pw.GRFML.BW.Pec <- paste ("Pec_Pw_GRFML_BW_", SaveDate, ".csv", sep="")
write.table(Pw.GRF.InterpML.BW.Combined.Pec, file=Save.Pw.GRFML.BW.Pec, sep =",", row.names=FALSE)
Save.Pw.GRFML.BW.Pel <- paste ("Pel_Pw_GRFML_BW_", SaveDate, ".csv", sep="")
write.table(Pw.GRF.InterpML.BW.Combined.Pel, file=Save.Pw.GRFML.BW.Pel, sep =",", row.names=FALSE)


setwd('/Users/SandyMKawano/Desktop/Research/Forelimbs on Force Plates/Force Data/R Analysis/Step 4 Results/GRF Hz')
Save.Af.GRFHz.BW.Pec <- paste ("Pec_Af_GRFHz_BW_", SaveDate, ".csv", sep="")
write.table(Af.GRF.InterpHz.BW.Combined.Pec, file=Save.Af.GRFHz.BW.Pec, sep =",", row.names=FALSE)
Save.Pb.GRFHz.BW.Pec <- paste ("Pec_Pb_GRFHz_BW_", SaveDate, ".csv", sep="")
write.table(Pb.GRF.InterpHz.BW.Combined.Pec, file=Save.Pb.GRFHz.BW.Pec, sep =",", row.names=FALSE)
Save.Af.GRFHz.BW.Pel <- paste ("Pel_Af_GRFHz_BW_", SaveDate, ".csv", sep="")
write.table(Af.GRF.InterpHz.BW.Combined.Pel, file=Save.Af.GRFHz.BW.Pel, sep =",", row.names=FALSE)
Save.Pw.GRFHz.BW.Pec <- paste ("Pec_Pw_GRFHz_BW_", SaveDate, ".csv", sep="")
write.table(Pw.GRF.InterpHz.BW.Combined.Pec, file=Save.Pw.GRFHz.BW.Pec, sep =",", row.names=FALSE)
Save.Pw.GRFHz.BW.Pel <- paste ("Pel_Pw_GRFHz_BW_", SaveDate, ".csv", sep="")
write.table(Pw.GRF.InterpHz.BW.Combined.Pel, file=Save.Pw.GRFHz.BW.Pel, sep =",", row.names=FALSE)


setwd('/Users/SandyMKawano/Desktop/Research/Forelimbs on Force Plates/Force Data/R Analysis/Step 4 Results/Net GRF')
Save.Af.NetGRF.BW.Pec <- paste ("Pec_Af_NetGRF_BW_", SaveDate, ".csv", sep="")
write.table(Af.GRF.NetGRF.BW.Combined.Pec, file=Save.Af.NetGRF.BW.Pec, sep =",", row.names=FALSE)
Save.Pb.NetGRF.BW.Pec <- paste ("Pec_Pb_NetGRF_BW_", SaveDate, ".csv", sep="")
write.table(Pb.GRF.NetGRF.BW.Combined.Pec, file=Save.Pb.NetGRF.BW.Pec, sep =",", row.names=FALSE)
Save.Af.NetGRF.BW.Pel <- paste ("Pel_Af_NetGRF_BW_", SaveDate, ".csv", sep="")
write.table(Af.GRF.NetGRF.BW.Combined.Pel, file=Save.Af.NetGRF.BW.Pel, sep =",", row.names=FALSE)
Save.Pw.NetGRF.BW.Pec <- paste ("Pec_Pw_NetGRF_BW_", SaveDate, ".csv", sep="")
write.table(Pw.GRF.NetGRF.BW.Combined.Pec, file=Save.Pw.NetGRF.BW.Pec, sep =",", row.names=FALSE)
Save.Pw.NetGRF.BW.Pel <- paste ("Pel_Pw_NetGRF_BW_", SaveDate, ".csv", sep="")
write.table(Pw.GRF.NetGRF.BW.Combined.Pel, file=Save.Pw.NetGRF.BW.Pel, sep =",", row.names=FALSE)



####################  PLOTTING GRF DATA ####################################

# Calculating standard errors at each 5% increment of stance
library(sciplot) # used for calculating SE
#Af.GRF.APAngle.Combined.Pec.SE  <- sapply(Af.GRF.APAngle.Combined.Pec[,c(1:21)], FUN=function(x) se(x, na.rm=TRUE))
Af.GRF.APAngle.Convert.Combined.Pec.SE  <- sapply(Af.GRF.APAngle.Convert.Combined.Pec[,c(1:21)], FUN=function(x) se(x, na.rm=TRUE))
Af.GRF.MLAngle.Convert.Combined.Pec.SE  <- sapply(Af.GRF.MLAngle.Convert.Combined.Pec[,c(1:21)], FUN=function(x) se(x, na.rm=TRUE))
Af.GRF.InterpV.BW.Combined.Pec.SE  <- sapply(Af.GRF.InterpV.BW.Combined.Pec[,c(1:21)], FUN=function(x) se(x, na.rm=TRUE))
Af.GRF.InterpML.BW.Combined.Pec.SE  <- sapply(Af.GRF.InterpML.BW.Combined.Pec[,c(1:21)], FUN=function(x) se(x, na.rm=TRUE))
Af.GRF.InterpHz.BW.Combined.Pec.SE  <- sapply(Af.GRF.InterpHz.BW.Combined.Pec[,c(1:21)], FUN=function(x) se(x, na.rm=TRUE))
Af.GRF.NetGRF.BW.Combined.Pec.SE  <- sapply(Af.GRF.NetGRF.BW.Combined.Pec[,c(1:21)], FUN=function(x) se(x, na.rm=TRUE))

#Af.GRF.APAngle.Combined.Pel.SE  <- sapply(Af.GRF.APAngle.Combined.Pel[,c(1:21)], FUN=function(x) se(x, na.rm=TRUE))
Af.GRF.APAngle.Convert.Combined.Pel.SE  <- sapply(Af.GRF.APAngle.Convert.Combined.Pel[,c(1:21)], FUN=function(x) se(x, na.rm=TRUE))
Af.GRF.MLAngle.Convert.Combined.Pel.SE  <- sapply(Af.GRF.MLAngle.Convert.Combined.Pel[,c(1:21)], FUN=function(x) se(x, na.rm=TRUE))
Af.GRF.InterpV.BW.Combined.Pel.SE  <- sapply(Af.GRF.InterpV.BW.Combined.Pel[,c(1:21)], FUN=function(x) se(x, na.rm=TRUE))
Af.GRF.InterpML.BW.Combined.Pel.SE  <- sapply(Af.GRF.InterpML.BW.Combined.Pel[,c(1:21)], FUN=function(x) se(x, na.rm=TRUE))
Af.GRF.InterpHz.BW.Combined.Pel.SE  <- sapply(Af.GRF.InterpHz.BW.Combined.Pel[,c(1:21)], FUN=function(x) se(x, na.rm=TRUE))
Af.GRF.NetGRF.BW.Combined.Pel.SE  <- sapply(Af.GRF.NetGRF.BW.Combined.Pel[,c(1:21)], FUN=function(x) se(x, na.rm=TRUE))

#Pb.GRF.APAngle.Combined.Pec.SE  <- sapply(Pb.GRF.APAngle.Combined.Pec[,c(1:21)], FUN=function(x) se(x, na.rm=TRUE))
Pb.GRF.APAngle.Convert.Combined.Pec.SE  <- sapply(Pb.GRF.APAngle.Convert.Combined.Pec[,c(1:21)], FUN=function(x) se(x, na.rm=TRUE))
Pb.GRF.MLAngle.Convert.Combined.Pec.SE  <- sapply(Pb.GRF.MLAngle.Convert.Combined.Pec[,c(1:21)], FUN=function(x) se(x, na.rm=TRUE))
Pb.GRF.InterpV.BW.Combined.Pec.SE  <- sapply(Pb.GRF.InterpV.BW.Combined.Pec[,c(1:21)], FUN=function(x) se(x, na.rm=TRUE))
Pb.GRF.InterpML.BW.Combined.Pec.SE  <- sapply(Pb.GRF.InterpML.BW.Combined.Pec[,c(1:21)], FUN=function(x) se(x, na.rm=TRUE))
Pb.GRF.InterpHz.BW.Combined.Pec.SE  <- sapply(Pb.GRF.InterpHz.BW.Combined.Pec[,c(1:21)], FUN=function(x) se(x, na.rm=TRUE))
Pb.GRF.NetGRF.BW.Combined.Pec.SE  <- sapply(Pb.GRF.NetGRF.BW.Combined.Pec[,c(1:21)], FUN=function(x) se(x, na.rm=TRUE))

#Pw.GRF.APAngle.Combined.Pec.SE  <- sapply(Pw.GRF.APAngle.Combined.Pec[,c(1:21)], FUN=function(x) se(x, na.rm=TRUE))
Pw.GRF.APAngle.Convert.Combined.Pec.SE  <- sapply(Pw.GRF.APAngle.Convert.Combined.Pec[,c(1:21)], FUN=function(x) se(x, na.rm=TRUE))
Pw.GRF.MLAngle.Convert.Combined.Pec.SE  <- sapply(Pw.GRF.MLAngle.Convert.Combined.Pec[,c(1:21)], FUN=function(x) se(x, na.rm=TRUE))
Pw.GRF.InterpV.BW.Combined.Pec.SE  <- sapply(Pw.GRF.InterpV.BW.Combined.Pec[,c(1:21)], FUN=function(x) se(x, na.rm=TRUE))
Pw.GRF.InterpML.BW.Combined.Pec.SE  <- sapply(Pw.GRF.InterpML.BW.Combined.Pec[,c(1:21)], FUN=function(x) se(x, na.rm=TRUE))
Pw.GRF.InterpHz.BW.Combined.Pec.SE  <- sapply(Pw.GRF.InterpHz.BW.Combined.Pec[,c(1:21)], FUN=function(x) se(x, na.rm=TRUE))
Pw.GRF.NetGRF.BW.Combined.Pec.SE  <- sapply(Pw.GRF.NetGRF.BW.Combined.Pec[,c(1:21)], FUN=function(x) se(x, na.rm=TRUE))

#Pw.GRF.APAngle.Combined.Pel.SE  <- sapply(Pw.GRF.APAngle.Combined.Pel[,c(1:21)], FUN=function(x) se(x, na.rm=TRUE))
Pw.GRF.APAngle.Convert.Combined.Pel.SE  <- sapply(Pw.GRF.APAngle.Convert.Combined.Pel[,c(1:21)], FUN=function(x) se(x, na.rm=TRUE))
Pw.GRF.MLAngle.Convert.Combined.Pel.SE  <- sapply(Pw.GRF.MLAngle.Convert.Combined.Pel[,c(1:21)], FUN=function(x) se(x, na.rm=TRUE))
Pw.GRF.InterpV.BW.Combined.Pel.SE  <- sapply(Pw.GRF.InterpV.BW.Combined.Pel[,c(1:21)], FUN=function(x) se(x, na.rm=TRUE))
Pw.GRF.InterpML.BW.Combined.Pel.SE  <- sapply(Pw.GRF.InterpML.BW.Combined.Pel[,c(1:21)], FUN=function(x) se(x, na.rm=TRUE))
Pw.GRF.InterpHz.BW.Combined.Pel.SE  <- sapply(Pw.GRF.InterpHz.BW.Combined.Pel[,c(1:21)], FUN=function(x) se(x, na.rm=TRUE))
Pw.GRF.NetGRF.BW.Combined.Pel.SE  <- sapply(Pw.GRF.NetGRF.BW.Combined.Pel[,c(1:21)], FUN=function(x) se(x, na.rm=TRUE))

#Sl.GRF.APAngle.Combined.Pec.SE  <- sapply(Sl.GRF.APAngle.Combined.Pec[,c(1:21)], FUN=function(x) se(x, na.rm=TRUE))
Sl.GRF.APAngle.Convert.Combined.Pec.SE  <- sapply(Sl.GRF.APAngle.Convert.Combined.Pec[,c(1:21)], FUN=function(x) se(x, na.rm=TRUE))
Sl.GRF.MLAngle.Convert.Combined.Pec.SE  <- sapply(Sl.GRF.MLAngle.Convert.Combined.Pec[,c(1:21)], FUN=function(x) se(x, na.rm=TRUE))
Sl.GRF.InterpV.BW.Combined.Pec.SE  <- sapply(Sl.GRF.InterpV.BW.Combined.Pec[,c(1:21)], FUN=function(x) se(x, na.rm=TRUE))
Sl.GRF.InterpML.BW.Combined.Pec.SE  <- sapply(Sl.GRF.InterpML.BW.Combined.Pec[,c(1:21)], FUN=function(x) se(x, na.rm=TRUE))
Sl.GRF.InterpHz.BW.Combined.Pec.SE  <- sapply(Sl.GRF.InterpHz.BW.Combined.Pec[,c(1:21)], FUN=function(x) se(x, na.rm=TRUE))
Sl.GRF.NetGRF.BW.Combined.Pec.SE  <- sapply(Sl.GRF.NetGRF.BW.Combined.Pec[,c(1:21)], FUN=function(x) se(x, na.rm=TRUE))

#Af.APAngle.Pec.Mean <- sapply(Af.GRF.APAngle.Combined.Pec[,c(1:21)], FUN=function(x) mean(x, na.rm=TRUE))
Af.APAngle.Pec.Mean <- sapply(Af.GRF.APAngle.Convert.Combined.Pec[,c(1:21)], FUN=function(x) mean(x, na.rm=TRUE))
Af.MLAngle.Pec.Mean <- sapply(Af.GRF.MLAngle.Convert.Combined.Pec[,c(1:21)], FUN=function(x) mean(x, na.rm=TRUE))
Af.InterpV.Pec.Mean <- sapply(Af.GRF.InterpV.BW.Combined.Pec[,c(1:21)], FUN=function(x) mean(x, na.rm=TRUE))
Af.InterpML.Pec.Mean <- sapply(Af.GRF.InterpML.BW.Combined.Pec[,c(1:21)], FUN=function(x) mean(x, na.rm=TRUE))
Af.InterpHz.Pec.Mean <- sapply(Af.GRF.InterpHz.BW.Combined.Pec[,c(1:21)], FUN=function(x) mean(x, na.rm=TRUE))
Af.NetGRF.Pec.Mean <- sapply(Af.GRF.NetGRF.BW.Combined.Pec[,c(1:21)], FUN=function(x) mean(x, na.rm=TRUE))

#Af.APAngle.Pel.Mean <- sapply(Af.GRF.APAngle.Combined.Pel[,c(1:21)], FUN=function(x) mean(x, na.rm=TRUE))
Af.APAngle.Pel.Mean <- sapply(Af.GRF.APAngle.Convert.Combined.Pel[,c(1:21)], FUN=function(x) mean(x, na.rm=TRUE))
Af.MLAngle.Pel.Mean <- sapply(Af.GRF.MLAngle.Convert.Combined.Pel[,c(1:21)], FUN=function(x) mean(x, na.rm=TRUE))
Af.InterpV.Pel.Mean <- sapply(Af.GRF.InterpV.BW.Combined.Pel[,c(1:21)], FUN=function(x) mean(x, na.rm=TRUE))
Af.InterpML.Pel.Mean <- sapply(Af.GRF.InterpML.BW.Combined.Pel[,c(1:21)], FUN=function(x) mean(x, na.rm=TRUE))
Af.InterpHz.Pel.Mean <- sapply(Af.GRF.InterpHz.BW.Combined.Pel[,c(1:21)], FUN=function(x) mean(x, na.rm=TRUE))
Af.NetGRF.Pel.Mean <- sapply(Af.GRF.NetGRF.BW.Combined.Pel[,c(1:21)], FUN=function(x) mean(x, na.rm=TRUE))

#Pb.APAngle.Pec.Mean <- sapply(Pb.GRF.APAngle.Combined.Pec[,c(1:21)], FUN=function(x) mean(x, na.rm=TRUE))
Pb.APAngle.Pec.Mean <- sapply(Pb.GRF.APAngle.Convert.Combined.Pec[,c(1:21)], FUN=function(x) mean(x, na.rm=TRUE))
Pb.MLAngle.Pec.Mean <- sapply(Pb.GRF.MLAngle.Convert.Combined.Pec[,c(1:21)], FUN=function(x) mean(x, na.rm=TRUE))
Pb.InterpV.Pec.Mean <- sapply(Pb.GRF.InterpV.BW.Combined.Pec[,c(1:21)], FUN=function(x) mean(x, na.rm=TRUE))
Pb.InterpML.Pec.Mean <- sapply(Pb.GRF.InterpML.BW.Combined.Pec[,c(1:21)], FUN=function(x) mean(x, na.rm=TRUE))
Pb.InterpHz.Pec.Mean <- sapply(Pb.GRF.InterpHz.BW.Combined.Pec[,c(1:21)], FUN=function(x) mean(x, na.rm=TRUE))
Pb.NetGRF.Pec.Mean <- sapply(Pb.GRF.NetGRF.BW.Combined.Pec[,c(1:21)], FUN=function(x) mean(x, na.rm=TRUE))

#Pw.APAngle.Pec.Mean <- sapply(Pw.GRF.APAngle.Combined.Pec[,c(1:21)], FUN=function(x) mean(x, na.rm=TRUE))
Pw.APAngle.Pec.Mean <- sapply(Pw.GRF.APAngle.Convert.Combined.Pec[,c(1:21)], FUN=function(x) mean(x, na.rm=TRUE))
Pw.MLAngle.Pec.Mean <- sapply(Pw.GRF.MLAngle.Convert.Combined.Pec[,c(1:21)], FUN=function(x) mean(x, na.rm=TRUE))
Pw.InterpV.Pec.Mean <- sapply(Pw.GRF.InterpV.BW.Combined.Pec[,c(1:21)], FUN=function(x) mean(x, na.rm=TRUE))
Pw.InterpML.Pec.Mean <- sapply(Pw.GRF.InterpML.BW.Combined.Pec[,c(1:21)], FUN=function(x) mean(x, na.rm=TRUE))
Pw.InterpHz.Pec.Mean <- sapply(Pw.GRF.InterpHz.BW.Combined.Pec[,c(1:21)], FUN=function(x) mean(x, na.rm=TRUE))
Pw.NetGRF.Pec.Mean <- sapply(Pw.GRF.NetGRF.BW.Combined.Pec[,c(1:21)], FUN=function(x) mean(x, na.rm=TRUE))

#Pw.APAngle.Pel.Mean <- sapply(Pw.GRF.APAngle.Combined.Pel[,c(1:21)], FUN=function(x) mean(x, na.rm=TRUE))
Pw.APAngle.Pel.Mean <- sapply(Pw.GRF.APAngle.Convert.Combined.Pel[,c(1:21)], FUN=function(x) mean(x, na.rm=TRUE))
Pw.MLAngle.Pel.Mean <- sapply(Pw.GRF.MLAngle.Convert.Combined.Pel[,c(1:21)], FUN=function(x) mean(x, na.rm=TRUE))
Pw.InterpV.Pel.Mean <- sapply(Pw.GRF.InterpV.BW.Combined.Pel[,c(1:21)], FUN=function(x) mean(x, na.rm=TRUE))
Pw.InterpML.Pel.Mean <- sapply(Pw.GRF.InterpML.BW.Combined.Pel[,c(1:21)], FUN=function(x) mean(x, na.rm=TRUE))
Pw.InterpHz.Pel.Mean <- sapply(Pw.GRF.InterpHz.BW.Combined.Pel[,c(1:21)], FUN=function(x) mean(x, na.rm=TRUE))
Pw.NetGRF.Pel.Mean <- sapply(Pw.GRF.NetGRF.BW.Combined.Pel[,c(1:21)], FUN=function(x) mean(x, na.rm=TRUE))

#Sl.APAngle.Pec.Mean <- sapply(Sl.GRF.APAngle.Combined.Pec[,c(1:21)], FUN=function(x) mean(x, na.rm=TRUE))
Sl.APAngle.Pec.Mean <- sapply(Sl.GRF.APAngle.Convert.Combined.Pec[,c(1:21)], FUN=function(x) mean(x, na.rm=TRUE))
Sl.MLAngle.Pec.Mean <- sapply(Sl.GRF.MLAngle.Convert.Combined.Pec[,c(1:21)], FUN=function(x) mean(x, na.rm=TRUE))
Sl.InterpV.Pec.Mean <- sapply(Sl.GRF.InterpV.BW.Combined.Pec[,c(1:21)], FUN=function(x) mean(x, na.rm=TRUE))
Sl.InterpML.Pec.Mean <- sapply(Sl.GRF.InterpML.BW.Combined.Pec[,c(1:21)], FUN=function(x) mean(x, na.rm=TRUE))
Sl.InterpHz.Pec.Mean <- sapply(Sl.GRF.InterpHz.BW.Combined.Pec[,c(1:21)], FUN=function(x) mean(x, na.rm=TRUE))
Sl.NetGRF.Pec.Mean <- sapply(Sl.GRF.NetGRF.BW.Combined.Pec[,c(1:21)], FUN=function(x) mean(x, na.rm=TRUE))


# Saving summary data
#Af.APAngle.Pec.Summary <- rbind(Af.APAngle.Pec.Mean, Af.GRF.APAngle.Combined.Pec.SE); Af.APAngle.Pec.Summary <- cbind(rbind("Mean", "Standard Error"), Af.APAngle.Pec.Summary)
Af.APAngle.Pec.Summary <- rbind(Af.APAngle.Pec.Mean, Af.GRF.APAngle.Convert.Combined.Pec.SE); Af.APAngle.Pec.Summary <- cbind(rbind("Mean", "Standard Error"), Af.APAngle.Pec.Summary)
Af.MLAngle.Pec.Summary <- rbind(Af.MLAngle.Pec.Mean, Af.GRF.MLAngle.Convert.Combined.Pec.SE); Af.MLAngle.Pec.Summary <- cbind(rbind("Mean", "Standard Error"), Af.MLAngle.Pec.Summary)
Af.InterpV.Pec.Summary <- rbind(Af.InterpV.Pec.Mean, Af.GRF.InterpV.BW.Combined.Pec.SE); Af.InterpV.Pec.Summary <- cbind(rbind("Mean", "Standard Error"), Af.InterpV.Pec.Summary)
Af.InterpML.Pec.Summary <- rbind(Af.InterpML.Pec.Mean, Af.GRF.InterpML.BW.Combined.Pec.SE); Af.InterpML.Pec.Summary <- cbind(rbind("Mean", "Standard Error"), Af.InterpML.Pec.Summary)
Af.InterpHz.Pec.Summary <- rbind(Af.InterpHz.Pec.Mean, Af.GRF.InterpHz.BW.Combined.Pec.SE); Af.InterpHz.Pec.Summary <- cbind(rbind("Mean", "Standard Error"), Af.InterpHz.Pec.Summary)
Af.NetGRF.Pec.Summary <- rbind(Af.NetGRF.Pec.Mean, Af.GRF.NetGRF.BW.Combined.Pec.SE); Af.NetGRF.Pec.Summary <- cbind(rbind("Mean", "Standard Error"), Af.NetGRF.Pec.Summary)

#Af.APAngle.Pel.Summary <- rbind(Af.APAngle.Pel.Mean, Af.GRF.APAngle.Combined.Pel.SE); Af.APAngle.Pel.Summary <- cbind(rbind("Mean", "Standard Error"), Af.APAngle.Pel.Summary)
Af.APAngle.Pel.Summary <- rbind(Af.APAngle.Pel.Mean, Af.GRF.APAngle.Convert.Combined.Pel.SE); Af.APAngle.Pel.Summary <- cbind(rbind("Mean", "Standard Error"), Af.APAngle.Pel.Summary)
Af.MLAngle.Pel.Summary <- rbind(Af.MLAngle.Pel.Mean, Af.GRF.MLAngle.Convert.Combined.Pel.SE); Af.MLAngle.Pel.Summary <- cbind(rbind("Mean", "Standard Error"), Af.MLAngle.Pel.Summary)
Af.InterpV.Pel.Summary <- rbind(Af.InterpV.Pel.Mean, Af.GRF.InterpV.BW.Combined.Pel.SE); Af.InterpV.Pel.Summary <- cbind(rbind("Mean", "Standard Error"), Af.InterpV.Pel.Summary)
Af.InterpML.Pel.Summary <- rbind(Af.InterpML.Pel.Mean, Af.GRF.InterpML.BW.Combined.Pel.SE); Af.InterpML.Pel.Summary <- cbind(rbind("Mean", "Standard Error"), Af.InterpML.Pel.Summary)
Af.InterpHz.Pel.Summary <- rbind(Af.InterpHz.Pel.Mean, Af.GRF.InterpHz.BW.Combined.Pel.SE); Af.InterpHz.Pel.Summary <- cbind(rbind("Mean", "Standard Error"), Af.InterpHz.Pel.Summary)
Af.NetGRF.Pel.Summary <- rbind(Af.NetGRF.Pel.Mean, Af.GRF.NetGRF.BW.Combined.Pel.SE); Af.NetGRF.Pel.Summary <- cbind(rbind("Mean", "Standard Error"), Af.NetGRF.Pel.Summary)

#Pb.APAngle.Pec.Summary <- rbind(Pb.APAngle.Pec.Mean, Pb.GRF.APAngle.Combined.Pec.SE); Pb.APAngle.Pec.Summary <- cbind(rbind("Mean", "Standard Error"), Pb.APAngle.Pec.Summary)
Pb.APAngle.Pec.Summary <- rbind(Pb.APAngle.Pec.Mean, Pb.GRF.APAngle.Convert.Combined.Pec.SE); Pb.APAngle.Pec.Summary <- cbind(rbind("Mean", "Standard Error"), Pb.APAngle.Pec.Summary)
Pb.MLAngle.Pec.Summary <- rbind(Pb.MLAngle.Pec.Mean, Pb.GRF.MLAngle.Convert.Combined.Pec.SE); Pb.MLAngle.Pec.Summary <- cbind(rbind("Mean", "Standard Error"), Pb.MLAngle.Pec.Summary)
Pb.InterpV.Pec.Summary <- rbind(Pb.InterpV.Pec.Mean, Pb.GRF.InterpV.BW.Combined.Pec.SE); Pb.InterpV.Pec.Summary <- cbind(rbind("Mean", "Standard Error"), Pb.InterpV.Pec.Summary)
Pb.InterpML.Pec.Summary <- rbind(Pb.InterpML.Pec.Mean, Pb.GRF.InterpML.BW.Combined.Pec.SE); Pb.InterpML.Pec.Summary <- cbind(rbind("Mean", "Standard Error"), Pb.InterpML.Pec.Summary)
Pb.InterpHz.Pec.Summary <- rbind(Pb.InterpHz.Pec.Mean, Pb.GRF.InterpHz.BW.Combined.Pec.SE); Pb.InterpHz.Pec.Summary <- cbind(rbind("Mean", "Standard Error"), Pb.InterpHz.Pec.Summary)
Pb.NetGRF.Pec.Summary <- rbind(Pb.NetGRF.Pec.Mean, Pb.GRF.NetGRF.BW.Combined.Pec.SE); Pb.NetGRF.Pec.Summary <- cbind(rbind("Mean", "Standard Error"), Pb.NetGRF.Pec.Summary)

#Pw.APAngle.Pec.Summary <- rbind(Pw.APAngle.Pec.Mean, Pw.GRF.APAngle.Combined.Pec.SE); Pw.APAngle.Pec.Summary <- cbind(rbind("Mean", "Standard Error"), Pw.APAngle.Pec.Summary)
Pw.APAngle.Pec.Summary <- rbind(Pw.APAngle.Pec.Mean, Pw.GRF.APAngle.Convert.Combined.Pec.SE); Pw.APAngle.Pec.Summary <- cbind(rbind("Mean", "Standard Error"), Pw.APAngle.Pec.Summary)
Pw.MLAngle.Pec.Summary <- rbind(Pw.MLAngle.Pec.Mean, Pw.GRF.MLAngle.Convert.Combined.Pec.SE); Pw.MLAngle.Pec.Summary <- cbind(rbind("Mean", "Standard Error"), Pw.MLAngle.Pec.Summary)
Pw.InterpV.Pec.Summary <- rbind(Pw.InterpV.Pec.Mean, Pw.GRF.InterpV.BW.Combined.Pec.SE); Pw.InterpV.Pec.Summary <- cbind(rbind("Mean", "Standard Error"), Pw.InterpV.Pec.Summary)
Pw.InterpML.Pec.Summary <- rbind(Pw.InterpML.Pec.Mean, Pw.GRF.InterpML.BW.Combined.Pec.SE); Pw.InterpML.Pec.Summary <- cbind(rbind("Mean", "Standard Error"), Pw.InterpML.Pec.Summary)
Pw.InterpHz.Pec.Summary <- rbind(Pw.InterpHz.Pec.Mean, Pw.GRF.InterpHz.BW.Combined.Pec.SE); Pw.InterpHz.Pec.Summary <- cbind(rbind("Mean", "Standard Error"), Pw.InterpHz.Pec.Summary)
Pw.NetGRF.Pec.Summary <- rbind(Pw.NetGRF.Pec.Mean, Pw.GRF.NetGRF.BW.Combined.Pec.SE); Pw.NetGRF.Pec.Summary <- cbind(rbind("Mean", "Standard Error"), Pw.NetGRF.Pec.Summary)

#Pw.APAngle.Pel.Summary <- rbind(Pw.APAngle.Pel.Mean, Pw.GRF.APAngle.Combined.Pel.SE); Pw.APAngle.Pel.Summary <- cbind(rbind("Mean", "Standard Error"), Pw.APAngle.Pel.Summary)
Pw.APAngle.Pel.Summary <- rbind(Pw.APAngle.Pel.Mean, Pw.GRF.APAngle.Convert.Combined.Pel.SE); Pw.APAngle.Pel.Summary <- cbind(rbind("Mean", "Standard Error"), Pw.APAngle.Pel.Summary)
Pw.MLAngle.Pel.Summary <- rbind(Pw.MLAngle.Pel.Mean, Pw.GRF.MLAngle.Convert.Combined.Pel.SE); Pw.MLAngle.Pel.Summary <- cbind(rbind("Mean", "Standard Error"), Pw.MLAngle.Pel.Summary)
Pw.InterpV.Pel.Summary <- rbind(Pw.InterpV.Pel.Mean, Pw.GRF.InterpV.BW.Combined.Pel.SE); Pw.InterpV.Pel.Summary <- cbind(rbind("Mean", "Standard Error"), Pw.InterpV.Pel.Summary)
Pw.InterpML.Pel.Summary <- rbind(Pw.InterpML.Pel.Mean, Pw.GRF.InterpML.BW.Combined.Pel.SE); Pw.InterpML.Pel.Summary <- cbind(rbind("Mean", "Standard Error"), Pw.InterpML.Pel.Summary)
Pw.InterpHz.Pel.Summary <- rbind(Pw.InterpHz.Pel.Mean, Pw.GRF.InterpHz.BW.Combined.Pel.SE); Pw.InterpHz.Pel.Summary <- cbind(rbind("Mean", "Standard Error"), Pw.InterpHz.Pel.Summary)
Pw.NetGRF.Pel.Summary <- rbind(Pw.NetGRF.Pel.Mean, Pw.GRF.NetGRF.BW.Combined.Pel.SE); Pw.NetGRF.Pel.Summary <- cbind(rbind("Mean", "Standard Error"), Pw.NetGRF.Pel.Summary)

#Sl.APAngle.Pec.Summary <- rbind(Sl.APAngle.Pec.Mean, Sl.GRF.APAngle.Combined.Pec.SE); Sl.APAngle.Pec.Summary <- cbind(rbind("Mean", "Standard Error"), Sl.APAngle.Pec.Summary)
Sl.APAngle.Pec.Summary <- rbind(Sl.APAngle.Pec.Mean, Sl.GRF.APAngle.Convert.Combined.Pec.SE); Sl.APAngle.Pec.Summary <- cbind(rbind("Mean", "Standard Error"), Sl.APAngle.Pec.Summary)
Sl.MLAngle.Pec.Summary <- rbind(Sl.MLAngle.Pec.Mean, Sl.GRF.MLAngle.Convert.Combined.Pec.SE); Sl.MLAngle.Pec.Summary <- cbind(rbind("Mean", "Standard Error"), Sl.MLAngle.Pec.Summary)
Sl.InterpV.Pec.Summary <- rbind(Sl.InterpV.Pec.Mean, Sl.GRF.InterpV.BW.Combined.Pec.SE); Sl.InterpV.Pec.Summary <- cbind(rbind("Mean", "Standard Error"), Sl.InterpV.Pec.Summary)
Sl.InterpML.Pec.Summary <- rbind(Sl.InterpML.Pec.Mean, Sl.GRF.InterpML.BW.Combined.Pec.SE); Sl.InterpML.Pec.Summary <- cbind(rbind("Mean", "Standard Error"), Sl.InterpML.Pec.Summary)
Sl.InterpHz.Pec.Summary <- rbind(Sl.InterpHz.Pec.Mean, Sl.GRF.InterpHz.BW.Combined.Pec.SE); Sl.InterpHz.Pec.Summary <- cbind(rbind("Mean", "Standard Error"), Sl.InterpHz.Pec.Summary)
Sl.NetGRF.Pec.Summary <- rbind(Sl.NetGRF.Pec.Mean, Sl.GRF.NetGRF.BW.Combined.Pec.SE); Sl.NetGRF.Pec.Summary <- cbind(rbind("Mean", "Standard Error"), Sl.NetGRF.Pec.Summary)


###### CREATING GRAPHS ######
# library(plotrix) # allows you to easily plot means + error bars
# 
# ForelimbColor <- "blue"
# HindlimbColor <- "green"
# PecFinColor <- "purple"

library(ggplot2)

# Reproducing colors used in default color palette in ggplot2
ggplotColours <- function(n=6, h=c(0, 360) +15){
  if ((diff(h)%%360) < 1) h[2] <- h[2] - 360/n
  hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}

Stance <- seq(0,100,5)

##  Saving graphs
setwd('/Users/SandyMKawano/Desktop/Research/Forelimbs on Force Plates/Force Data/R Analysis/Step 4 Results/Graphs')


## AP Angle 
Af.APAngle.Pec.Mean.SE <- data.frame(Af.APAngle.Pec.Mean, Af.GRF.APAngle.Convert.Combined.Pec.SE, Stance)
Af.APAngle.Pec.Mean.SE$Type <- "Terrestrial FL"
names(Af.APAngle.Pec.Mean.SE) <- c("Mean", "SE", "Stance", "Type")
Af.APAngle.Pel.Mean.SE <- data.frame(Af.APAngle.Pel.Mean, Af.GRF.APAngle.Convert.Combined.Pel.SE, Stance)
Af.APAngle.Pel.Mean.SE$Type <- "Terrestrial HL"
names(Af.APAngle.Pel.Mean.SE) <- c("Mean", "SE", "Stance", "Type")
Pb.APAngle.Pec.Mean.SE <- data.frame(Pb.APAngle.Pec.Mean, Pb.GRF.APAngle.Convert.Combined.Pec.SE, Stance)
Pb.APAngle.Pec.Mean.SE$Type <- "Terrestrial PF"
names(Pb.APAngle.Pec.Mean.SE) <- c("Mean", "SE", "Stance", "Type")
Pw.APAngle.Pec.Mean.SE <- data.frame(Pw.APAngle.Pec.Mean, Pw.GRF.APAngle.Convert.Combined.Pec.SE, Stance)
Pw.APAngle.Pec.Mean.SE$Type <- "Semi-Aquatic FL"
names(Pw.APAngle.Pec.Mean.SE) <- c("Mean", "SE", "Stance", "Type")
Pw.APAngle.Pel.Mean.SE <- data.frame(Pw.APAngle.Pel.Mean, Pw.GRF.APAngle.Convert.Combined.Pel.SE, Stance)
Pw.APAngle.Pel.Mean.SE$Type <- "Semi-Aquatic HL"
names(Pw.APAngle.Pel.Mean.SE) <- c("Mean", "SE", "Stance", "Type")
Sl.APAngle.Pec.Mean.SE <- data.frame(Sl.APAngle.Pec.Mean, Sl.GRF.APAngle.Convert.Combined.Pec.SE, Stance)
Sl.APAngle.Pec.Mean.SE$Type <- "Aquatic FL"
names(Sl.APAngle.Pec.Mean.SE) <- c("Mean", "SE", "Stance", "Type")


APAngle.Mean.SE <- rbind(Af.APAngle.Pec.Mean.SE, Af.APAngle.Pel.Mean.SE, Pb.APAngle.Pec.Mean.SE, Pw.APAngle.Pec.Mean.SE, Pw.APAngle.Pel.Mean.SE)
#APAngle.Mean.SE <- rbind(Af.APAngle.Pec.Mean.SE, Af.APAngle.Pel.Mean.SE, Pb.APAngle.Pec.Mean.SE, Pw.APAngle.Pec.Mean.SE, Pw.APAngle.Pel.Mean.SE, Sl.APAngle.Pec.Mean.SE)

APAngle.MaxMin <- aes(ymax=APAngle.Mean.SE$Mean+APAngle.Mean.SE$SE, ymin=APAngle.Mean.SE$Mean-APAngle.Mean.SE$SE)



tiff(filename=paste("APAngle.Compare_", SaveDate, ".tif", sep=""), width=1000, height=1000)
ggplot(data=APAngle.Mean.SE, aes(x=APAngle.Mean.SE$Stance, y=APAngle.Mean.SE$Mean, fill=Type, linetype=Type))+
  scale_y_continuous("Anteroposterior Angle (degrees)\n")+
  scale_x_continuous("\nStance (%)\n")+
  geom_line(size=1, alpha=0.75)+
  geom_ribbon(APAngle.MaxMin, alpha=0.5)+
  scale_colour_manual(values=c("ivory4", "ivory4", "ivory4", "ivory4", "ivory4"))+
  #scale_fill_manual(values=c("ivory4", "ivory4", "ivory4", "ivory4", "ivory4"))+
  scale_fill_manual(name="Appendage:", 
                    labels=c("Semi-quatic FL    ", "Semi-quatic HL    ", "Terrestrial FL    ", "Terrestrial HL  ", "Terrestrial PF  "),
                    #values=c("ivory4", "ivory4", "ivory4"))+
                    values=c("lightblue3", "steelblue3", "burlywood", "rosybrown4", "seagreen"))+
  scale_linetype_manual(values=c("dotted", "dashed", "dotted", "solid", "dashed"))+
  theme(axis.title.x=element_text(colour="black", size=30))+ # vjust=0 puts a little more spacing btwn the axis text and label
  theme(axis.title.y=element_text(colour='black', size=30))+
  theme(axis.text.x=element_text(colour='black', size=15))+
  theme(axis.text.y=element_text(colour='black', size=15))+
  theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank())+ # get rid of gridlines
  theme(panel.background=element_blank())+ # make background white
  theme(axis.line=element_line(colour="black", linetype="solid"))+ # put black lines for axes
  theme(legend.position="none")
dev.off()


# ML Angle
Af.MLAngle.Pec.Mean.SE <- data.frame(Af.MLAngle.Pec.Mean, Af.GRF.MLAngle.Convert.Combined.Pec.SE, Stance)
Af.MLAngle.Pec.Mean.SE$Type <- "Terrestrial FL"
names(Af.MLAngle.Pec.Mean.SE) <- c("Mean", "SE", "Stance", "Type")
Af.MLAngle.Pel.Mean.SE <- data.frame(Af.MLAngle.Pel.Mean, Af.GRF.MLAngle.Convert.Combined.Pel.SE, Stance)
Af.MLAngle.Pel.Mean.SE$Type <- "Terrestrial HL"
names(Af.MLAngle.Pel.Mean.SE) <- c("Mean", "SE", "Stance", "Type")
Pb.MLAngle.Pec.Mean.SE <- data.frame(Pb.MLAngle.Pec.Mean, Pb.GRF.MLAngle.Convert.Combined.Pec.SE, Stance)
Pb.MLAngle.Pec.Mean.SE$Type <- "Terrestrial PF"
names(Pb.MLAngle.Pec.Mean.SE) <- c("Mean", "SE", "Stance", "Type")
Pw.MLAngle.Pec.Mean.SE <- data.frame(Pw.MLAngle.Pec.Mean, Pw.GRF.MLAngle.Convert.Combined.Pec.SE, Stance)
Pw.MLAngle.Pec.Mean.SE$Type <- "Semi-aquatic FL"
names(Pw.MLAngle.Pec.Mean.SE) <- c("Mean", "SE", "Stance", "Type")
Pw.MLAngle.Pel.Mean.SE <- data.frame(Pw.MLAngle.Pel.Mean, Pw.GRF.MLAngle.Convert.Combined.Pel.SE, Stance)
Pw.MLAngle.Pel.Mean.SE$Type <- "Semi-aquatic HL"
names(Pw.MLAngle.Pel.Mean.SE) <- c("Mean", "SE", "Stance", "Type")
Sl.MLAngle.Pec.Mean.SE <- data.frame(Sl.MLAngle.Pec.Mean, Sl.GRF.MLAngle.Convert.Combined.Pec.SE, Stance)
Sl.MLAngle.Pec.Mean.SE$Type <- "Aquatic FL"
names(Sl.MLAngle.Pec.Mean.SE) <- c("Mean", "SE", "Stance", "Type")

MLAngle.Mean.SE <- rbind(Af.MLAngle.Pec.Mean.SE, Af.MLAngle.Pel.Mean.SE, Pb.MLAngle.Pec.Mean.SE, Pw.MLAngle.Pec.Mean.SE, Pw.MLAngle.Pel.Mean.SE)
#MLAngle.Mean.SE <- rbind(Af.MLAngle.Pec.Mean.SE, Af.MLAngle.Pel.Mean.SE, Pb.MLAngle.Pec.Mean.SE, Pw.MLAngle.Pec.Mean.SE, Pw.MLAngle.Pel.Mean.SE, Sl.MLAngle.Pec.Mean.SE)

MLAngle.MaxMin <- aes(ymax=MLAngle.Mean.SE$Mean+MLAngle.Mean.SE$SE, ymin=MLAngle.Mean.SE$Mean-MLAngle.Mean.SE$SE)


tiff(filename=paste("MLAngle.Compare_", SaveDate, ".tif", sep=""), width=1000, height=1000)
ggplot(data=MLAngle.Mean.SE, aes(x=MLAngle.Mean.SE$Stance, y=MLAngle.Mean.SE$Mean, fill=Type, linetype=Type))+
  scale_y_continuous("Mediolateral Angle (degrees)\n")+
  scale_x_continuous("\nStance (%)\n")+
  geom_line(size=1, alpha=0.75)+
  geom_ribbon(MLAngle.MaxMin, alpha=0.5)+
  scale_colour_manual(values=c("ivory4", "ivory4", "ivory4", "ivory4", "ivory4"))+
  #scale_fill_manual(values=c("ivory4", "ivory4", "ivory4"))+
  scale_fill_manual(name="Appendage:", 
                    labels=c("Semi-aquatic FL   ", "Semi-aquatic HL    ", "Terrestrial FL    ", "Terrestrial HL  ", "Terrestrial PF  "),
                    #values=c("ivory4", "ivory4", "ivory4"))+
                    values=c("lightblue3", "steelblue3", "burlywood", "rosybrown4", "seagreen"))+
  scale_linetype_manual(values=c("dashed", "dotted", "dashed", "dotted", "solid"))+
  theme(axis.title.x=element_text(colour="black", size=30))+ # vjust=0 puts a little more spacing btwn the axis text and label
  theme(axis.title.y=element_text(colour='black', size=30))+
  theme(axis.text.x=element_text(colour='black', size=15))+
  theme(axis.text.y=element_text(colour='black', size=15))+
  theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank())+ # get rid of gridlines
  theme(panel.background=element_blank())+ # make background white
  theme(axis.line=element_line(colour="black", linetype="solid"))+ # put black lines for axes
  theme(legend.position="none")
dev.off()


#Vertical GRF
Af.InterpV.Pec.Mean.SE <- data.frame(Af.InterpV.Pec.Mean, Af.GRF.InterpV.BW.Combined.Pec.SE, Stance)
Af.InterpV.Pec.Mean.SE$Type <- "Terrestrial FL"
names(Af.InterpV.Pec.Mean.SE) <- c("Mean", "SE", "Stance", "Type")
Af.InterpV.Pel.Mean.SE <- data.frame(Af.InterpV.Pel.Mean, Af.GRF.InterpV.BW.Combined.Pel.SE, Stance)
Af.InterpV.Pel.Mean.SE$Type <- "Terrestrial HL"
names(Af.InterpV.Pel.Mean.SE) <- c("Mean", "SE", "Stance", "Type")
Pb.InterpV.Pec.Mean.SE <- data.frame(Pb.InterpV.Pec.Mean, Pb.GRF.InterpV.BW.Combined.Pec.SE, Stance)
Pb.InterpV.Pec.Mean.SE$Type <- "Terrestrial PF"
names(Pb.InterpV.Pec.Mean.SE) <- c("Mean", "SE", "Stance", "Type")
Pw.InterpV.Pec.Mean.SE <- data.frame(Pw.InterpV.Pec.Mean, Pw.GRF.InterpV.BW.Combined.Pec.SE, Stance)
Pw.InterpV.Pec.Mean.SE$Type <- "Semi-aquatic FL"
names(Pw.InterpV.Pec.Mean.SE) <- c("Mean", "SE", "Stance", "Type")
Pw.InterpV.Pel.Mean.SE <- data.frame(Pw.InterpV.Pel.Mean, Pw.GRF.InterpV.BW.Combined.Pel.SE, Stance)
Pw.InterpV.Pel.Mean.SE$Type <- "Semi-aquatic HL"
names(Pw.InterpV.Pel.Mean.SE) <- c("Mean", "SE", "Stance", "Type")
Sl.InterpV.Pec.Mean.SE <- data.frame(Sl.InterpV.Pec.Mean, Sl.GRF.InterpV.BW.Combined.Pec.SE, Stance)
Sl.InterpV.Pec.Mean.SE$Type <- "Aquatic FL"
names(Sl.InterpV.Pec.Mean.SE) <- c("Mean", "SE", "Stance", "Type")


InterpV.Mean.SE <- rbind(Af.InterpV.Pec.Mean.SE, Af.InterpV.Pel.Mean.SE, Pb.InterpV.Pec.Mean.SE, Pw.InterpV.Pec.Mean.SE, Pw.InterpV.Pel.Mean.SE)
#InterpV.Mean.SE <- rbind(Af.InterpV.Pec.Mean.SE, Af.InterpV.Pel.Mean.SE, Pb.InterpV.Pec.Mean.SE, Pw.InterpV.Pec.Mean.SE, Pw.InterpV.Pel.Mean.SE, Sl.InterpV.Pec.Mean.SE)

InterpV.MaxMin <- aes(ymax=InterpV.Mean.SE$Mean+InterpV.Mean.SE$SE, ymin=InterpV.Mean.SE$Mean-InterpV.Mean.SE$SE)


tiff(filename=paste("VerticalBW.Compare_", SaveDate, ".tif", sep=""), width=1000, height=1000)
ggplot(data=InterpV.Mean.SE, aes(x=InterpV.Mean.SE$Stance, y=InterpV.Mean.SE$Mean, fill=Type, linetype=Type))+
  scale_y_continuous("Vertical GRF (BW)\n", limits=c(0,0.5))+
  scale_x_continuous("\nStance (%)\n")+
  geom_line(size=1, alpha=0.75)+
  geom_ribbon(InterpV.MaxMin, alpha=0.5)+
  scale_colour_manual(values=c("ivory4", "ivory4", "ivory4", "ivory4", "ivory4"))+
  #scale_fill_manual(values=c("ivory4", "ivory4", "ivory4"))+
  scale_fill_manual(name="Appendage:", 
                    labels=c("Semi-aquatic FL    ", "Semi-aquatic HL    ", "Terrestrial FL    ", "Terrestrial HL  ", "Terrestrial PF  "),
                    #values=c("ivory4", "ivory4", "ivory4"))+
                    values=c("lightblue3", "steelblue3", "burlywood", "rosybrown4", "seagreen"))+
  scale_linetype_manual(values=c("dashed", "dotted", "dashed", "dotted", "solid"))+
  theme(axis.title.x=element_text(colour="black", size=30))+ # vjust=0 puts a little more spacing btwn the axis text and label
  theme(axis.title.y=element_text(colour='black', size=30))+
  theme(axis.text.x=element_text(colour='black', size=15))+
  theme(axis.text.y=element_text(colour='black', size=15))+
  theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank())+ # get rid of gridlines
  theme(panel.background=element_blank())+ # make background white
  theme(axis.line=element_line(colour="black", linetype="solid"))+ # put black lines for axes
  theme(legend.position="none")
dev.off()



#Mediolateral GRF
Af.InterpML.Pec.Mean.SE <- data.frame(Af.InterpML.Pec.Mean, Af.GRF.InterpML.BW.Combined.Pec.SE, Stance)
Af.InterpML.Pec.Mean.SE$Type <- "Terrestrial FL"
names(Af.InterpML.Pec.Mean.SE) <- c("Mean", "SE", "Stance", "Type")
Af.InterpML.Pel.Mean.SE <- data.frame(Af.InterpML.Pel.Mean, Af.GRF.InterpML.BW.Combined.Pel.SE, Stance)
Af.InterpML.Pel.Mean.SE$Type <- "Terrestrial HL"
names(Af.InterpML.Pel.Mean.SE) <- c("Mean", "SE", "Stance", "Type")
Pb.InterpML.Pec.Mean.SE <- data.frame(Pb.InterpML.Pec.Mean, Pb.GRF.InterpML.BW.Combined.Pec.SE, Stance)
Pb.InterpML.Pec.Mean.SE$Type <- "Terrestrial PF"
names(Pb.InterpML.Pec.Mean.SE) <- c("Mean", "SE", "Stance", "Type")
Pw.InterpML.Pec.Mean.SE <- data.frame(Pw.InterpML.Pec.Mean, Pw.GRF.InterpML.BW.Combined.Pec.SE, Stance)
Pw.InterpML.Pec.Mean.SE$Type <- "Semi-aquatic FL"
names(Pw.InterpML.Pec.Mean.SE) <- c("Mean", "SE", "Stance", "Type")
Pw.InterpML.Pel.Mean.SE <- data.frame(Pw.InterpML.Pel.Mean, Pw.GRF.InterpML.BW.Combined.Pel.SE, Stance)
Pw.InterpML.Pel.Mean.SE$Type <- "Semi-aquatic HL"
names(Pw.InterpML.Pel.Mean.SE) <- c("Mean", "SE", "Stance", "Type")
Sl.InterpML.Pec.Mean.SE <- data.frame(Sl.InterpML.Pec.Mean, Sl.GRF.InterpML.BW.Combined.Pec.SE, Stance)
Sl.InterpML.Pec.Mean.SE$Type <- "Aquatic FL"
names(Sl.InterpML.Pec.Mean.SE) <- c("Mean", "SE", "Stance", "Type")

InterpML.Mean.SE <- rbind(Af.InterpML.Pec.Mean.SE, Af.InterpML.Pel.Mean.SE, Pb.InterpML.Pec.Mean.SE, Pw.InterpML.Pec.Mean.SE, Pw.InterpML.Pel.Mean.SE)
#InterpML.Mean.SE <- rbind(Af.InterpML.Pec.Mean.SE, Af.InterpML.Pel.Mean.SE, Pb.InterpML.Pec.Mean.SE, Pw.InterpML.Pec.Mean.SE, Pw.InterpML.Pel.Mean.SE, Sl.InterpML.Pec.Mean.SE)

InterpML.MaxMin <- aes(ymax=InterpML.Mean.SE$Mean+InterpML.Mean.SE$SE, ymin=InterpML.Mean.SE$Mean-InterpML.Mean.SE$SE)


tiff(filename=paste("MediolateralBW.Compare_", SaveDate, ".tif", sep=""), width=1000, height=1000)
ggplot(data=InterpML.Mean.SE, aes(x=InterpML.Mean.SE$Stance, y=InterpML.Mean.SE$Mean, fill=Type, linetype=Type))+
  scale_y_continuous("Mediolateral GRF (BW)\n")+
  scale_x_continuous("\nStance (%)\n")+
  geom_line(size=1, alpha=0.75)+
  geom_ribbon(InterpML.MaxMin, alpha=0.5)+
  scale_colour_manual(values=c("ivory4", "ivory4", "ivory4", "ivory4", "ivory4"))+
  #scale_fill_manual(values=c("ivory4", "ivory4", "ivory4"))+
  scale_fill_manual(name="Appendage:", 
                    labels=c("Semi-aquatic FL    ", "Semi-aquatic HL    ", "Terrestrial FL    ", "Terrestrial HL  ", "Terrestrial PF  "),
                    #values=c("ivory4", "ivory4", "ivory4"))+
                    values=c("lightblue3", "steelblue3", "burlywood", "rosybrown4", "seagreen"))+
  scale_linetype_manual(values=c("dashed", "dotted", "dashed", "dotted", "solid"))+
  theme(axis.title.x=element_text(colour="black", size=30))+ # vjust=0 puts a little more spacing btwn the axis text and label
  theme(axis.title.y=element_text(colour='black', size=30))+
  theme(axis.text.x=element_text(colour='black', size=15))+
  theme(axis.text.y=element_text(colour='black', size=15))+
  theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank())+ # get rid of gridlines
  theme(panel.background=element_blank())+ # make background white
  theme(axis.line=element_line(colour="black", linetype="solid"))+ # put black lines for axes
  theme(legend.position="none")
dev.off()



#Anteroposterior GRF
Af.InterpHz.Pec.Mean.SE <- data.frame(Af.InterpHz.Pec.Mean, Af.GRF.InterpHz.BW.Combined.Pec.SE, Stance)
Af.InterpHz.Pec.Mean.SE$Type <- "Terrestrial FL"
names(Af.InterpHz.Pec.Mean.SE) <- c("Mean", "SE", "Stance", "Type")
Af.InterpHz.Pel.Mean.SE <- data.frame(Af.InterpHz.Pel.Mean, Af.GRF.InterpHz.BW.Combined.Pel.SE, Stance)
Af.InterpHz.Pel.Mean.SE$Type <- "Terrestrial HL"
names(Af.InterpHz.Pel.Mean.SE) <- c("Mean", "SE", "Stance", "Type")
Pb.InterpHz.Pec.Mean.SE <- data.frame(Pb.InterpHz.Pec.Mean, Pb.GRF.InterpHz.BW.Combined.Pec.SE, Stance)
Pb.InterpHz.Pec.Mean.SE$Type <- "Terrestrial PF"
names(Pb.InterpHz.Pec.Mean.SE) <- c("Mean", "SE", "Stance", "Type")
Pw.InterpHz.Pec.Mean.SE <- data.frame(Pw.InterpHz.Pec.Mean, Pw.GRF.InterpHz.BW.Combined.Pec.SE, Stance)
Pw.InterpHz.Pec.Mean.SE$Type <- "Semi-aquatic FL"
names(Pw.InterpHz.Pec.Mean.SE) <- c("Mean", "SE", "Stance", "Type")
Pw.InterpHz.Pel.Mean.SE <- data.frame(Pw.InterpHz.Pel.Mean, Pw.GRF.InterpHz.BW.Combined.Pel.SE, Stance)
Pw.InterpHz.Pel.Mean.SE$Type <- "Semi-aquatic HL"
names(Pw.InterpHz.Pel.Mean.SE) <- c("Mean", "SE", "Stance", "Type")
Sl.InterpHz.Pec.Mean.SE <- data.frame(Sl.InterpHz.Pec.Mean, Sl.GRF.InterpHz.BW.Combined.Pec.SE, Stance)
Sl.InterpHz.Pec.Mean.SE$Type <- "Aquatic FL"
names(Sl.InterpHz.Pec.Mean.SE) <- c("Mean", "SE", "Stance", "Type")

InterpHz.Mean.SE <- rbind(Af.InterpHz.Pec.Mean.SE, Af.InterpHz.Pel.Mean.SE, Pb.InterpHz.Pec.Mean.SE, Pw.InterpHz.Pec.Mean.SE, Pw.InterpHz.Pel.Mean.SE)
#InterpHz.Mean.SE <- rbind(Af.InterpHz.Pec.Mean.SE, Af.InterpHz.Pel.Mean.SE, Pb.InterpHz.Pec.Mean.SE, Pw.InterpHz.Pec.Mean.SE, Pw.InterpHz.Pel.Mean.SE, Sl.InterpHz.Pec.Mean.SE)

InterpHz.MaxMin <- aes(ymax=InterpHz.Mean.SE$Mean+InterpHz.Mean.SE$SE, ymin=InterpHz.Mean.SE$Mean-InterpHz.Mean.SE$SE)


tiff(filename=paste("AnteroposteriorBW.Compare_", SaveDate, ".tif", sep=""), width=1000, height=1000)
ggplot(data=InterpHz.Mean.SE, aes(x=InterpHz.Mean.SE$Stance, y=InterpHz.Mean.SE$Mean, fill=Type, linetype=Type))+
  scale_y_continuous("Anteroposterior GRF (BW)\n")+
  scale_x_continuous("\nStance (%)\n")+
  geom_line(size=1, alpha=0.75)+
  geom_ribbon(InterpHz.MaxMin, alpha=0.5)+
  scale_colour_manual(values=c("ivory4", "ivory4", "ivory4", "ivory4"))+
  #scale_fill_manual(values=c("ivory4", "ivory4", "ivory4"))+
  scale_fill_manual(name="Appendage:", 
                    labels=c("Semi-aquatic FL    ", "Semi-aquatic HL    ", "Terrestrial FL    ", "Terrestrial HL  ", "Terrestrial PF  "),
                    #values=c("ivory4", "ivory4", "ivory4"))+
                    values=c("lightblue3", "steelblue3", "burlywood", "rosybrown4", "seagreen"))+
  scale_linetype_manual(values=c("dashed", "dotted", "dashed", "dotted", "solid"))+
  theme(axis.title.x=element_text(colour="black", size=30))+ # vjust=0 puts a little more spacing btwn the axis text and label
  theme(axis.title.y=element_text(colour='black', size=30))+
  theme(axis.text.x=element_text(colour='black', size=15))+
  theme(axis.text.y=element_text(colour='black', size=15))+
  theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank())+ # get rid of gridlines
  theme(panel.background=element_blank())+ # make background white
  theme(axis.line=element_line(colour="black", linetype="solid"))+ # put black lines for axes
  theme(legend.position="none")
dev.off()

#Net GRF
Af.NetGRF.Pec.Mean.SE <- data.frame(Af.NetGRF.Pec.Mean, Af.GRF.NetGRF.BW.Combined.Pec.SE, Stance)
Af.NetGRF.Pec.Mean.SE$Type <- "Terrestrial FL"
names(Af.NetGRF.Pec.Mean.SE) <- c("Mean", "SE", "Stance", "Type")
Af.NetGRF.Pel.Mean.SE <- data.frame(Af.NetGRF.Pel.Mean, Af.GRF.NetGRF.BW.Combined.Pel.SE, Stance)
Af.NetGRF.Pel.Mean.SE$Type <- "Terrestrial HL"
names(Af.NetGRF.Pel.Mean.SE) <- c("Mean", "SE", "Stance", "Type")
Pb.NetGRF.Pec.Mean.SE <- data.frame(Pb.NetGRF.Pec.Mean, Pb.GRF.NetGRF.BW.Combined.Pec.SE, Stance)
Pb.NetGRF.Pec.Mean.SE$Type <- "Terrestrial PF"
names(Pb.NetGRF.Pec.Mean.SE) <- c("Mean", "SE", "Stance", "Type")
Pw.NetGRF.Pec.Mean.SE <- data.frame(Pw.NetGRF.Pec.Mean, Pw.GRF.NetGRF.BW.Combined.Pec.SE, Stance)
Pw.NetGRF.Pec.Mean.SE$Type <- "Semi-aquatic FL"
names(Pw.NetGRF.Pec.Mean.SE) <- c("Mean", "SE", "Stance", "Type")
Pw.NetGRF.Pel.Mean.SE <- data.frame(Pw.NetGRF.Pel.Mean, Pw.GRF.NetGRF.BW.Combined.Pel.SE, Stance)
Pw.NetGRF.Pel.Mean.SE$Type <- "Semi-aquatic HL"
names(Pw.NetGRF.Pel.Mean.SE) <- c("Mean", "SE", "Stance", "Type")
Sl.NetGRF.Pec.Mean.SE <- data.frame(Sl.NetGRF.Pec.Mean, Sl.GRF.NetGRF.BW.Combined.Pec.SE, Stance)
Sl.NetGRF.Pec.Mean.SE$Type <- "Aquatic FL"
names(Sl.NetGRF.Pec.Mean.SE) <- c("Mean", "SE", "Stance", "Type")


NetGRF.Mean.SE <- rbind(Af.NetGRF.Pec.Mean.SE, Af.NetGRF.Pel.Mean.SE, Pb.NetGRF.Pec.Mean.SE, Pw.NetGRF.Pec.Mean.SE, Pw.NetGRF.Pel.Mean.SE)
#NetGRF.Mean.SE <- rbind(Af.NetGRF.Pec.Mean.SE, Af.NetGRF.Pel.Mean.SE, Pb.NetGRF.Pec.Mean.SE, Pw.NetGRF.Pec.Mean.SE, Pw.NetGRF.Pel.Mean.SE, Sl.NetGRF.Pec.Mean.SE)

NetGRF.MaxMin <- aes(ymax=NetGRF.Mean.SE$Mean+NetGRF.Mean.SE$SE, ymin=NetGRF.Mean.SE$Mean-NetGRF.Mean.SE$SE)

tiff(filename=paste("NetBW.Compare_", SaveDate, ".tif", sep=""), width=1000, height=1000)
ggplot(data=NetGRF.Mean.SE, aes(x=NetGRF.Mean.SE$Stance, y=NetGRF.Mean.SE$Mean, fill=Type, linetype=Type))+
  scale_y_continuous("Net GRF (BW)\n", limits=c(0.00,0.5))+
  scale_x_continuous("\nStance (%)\n")+
  geom_line(size=1, alpha=0.75)+
  geom_ribbon(NetGRF.MaxMin, alpha=0.5)+
#   scale_colour_manual(name="Appendage", 
#                       labels=c("Semi-aquatic FL  ", "Semi-aquatic HL  ", "Terrestrial FL  ", "Terrestrial HL  ", "Terrestrial PF  "),
#                       values=c("ivory4", "ivory4", "ivory4", "ivory4", "ivory4"))+
#   scale_fill_manual(name="Appendage", 
#                     labels=c("Forelimb (FL)", "Hind limb (HL)", "Pectoral Fin (PF)"),
#                     values=c("ivory4", "ivory4", "ivory4"))+
  scale_fill_manual(name="Appendage:", 
                    labels=c("Semi-aquatic FL  ", "Semi-aquatic HL  ", "Terrestrial FL  ", "Terrestrial HL  ", "Terrestrial PF  "),
                    #values=c("ivory4", "ivory4", "ivory4"))+
                    values=c("lightblue3", "steelblue3", "burlywood", "rosybrown4", "seagreen"))+
  scale_linetype_manual(name="Appendage", 
                        labels=c("Semi-aquatic FL  ", "Semi-aquatic HL  ", "Terrestrial FL  ", "Terrestrial HL  ", "Terrestrial PF  "),
                        values=c("dashed", "dotted", "dashed", "dotted", "solid"))+
  theme(axis.title.x=element_text(colour="black", size=30))+ # vjust=0 puts a little more spacing btwn the axis text and label
  theme(axis.title.y=element_text(colour='black', size=30))+
  theme(axis.text.x=element_text(colour='black', size=15))+
  theme(axis.text.y=element_text(colour='black', size=15))+
  theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank())+ # get rid of gridlines
  theme(panel.background=element_blank())+ # make background white
  theme(axis.line=element_line(colour="black", linetype="solid"))+ # put black lines for axes
  theme(legend.position="bottom")
dev.off()

##  Saving graphs
setwd('/Users/SandyMKawano/Desktop/Research/Forelimbs on Force Plates/Force Data/R Analysis/Step 4 Results/Graphs')

library(gridExtra) # for using grid.arrange()
# PdfSave <- paste("Appendage Force Comparison_", SaveDate, ".pdf", sep="")
# pdf(PdfSave)

APAngle.Plot <- ggplot(data=APAngle.Mean.SE, aes(x=APAngle.Mean.SE$Stance, y=APAngle.Mean.SE$Mean, fill=Type, linetype=Type))+
  scale_y_continuous("Anteroposterior Angle (degrees)\n")+
  #scale_x_continuous("\nStance (%)\n")+
  scale_x_continuous(element_blank())+
  geom_rect(aes(xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=0), fill="gray95")+
  geom_vline(aes(xintercept=mean(Pec.Af.Kept$PercentStance)), colour="burlywood", linetype="dashed", cex=1, alpha=0.25)+
  geom_vline(aes(xintercept=mean(Pel.Af.Kept$PercentStance)), colour="rosybrown4", linetype="dotted", cex=1, alpha=0.25)+
  geom_vline(aes(xintercept=mean(Pec.Pb.Kept$PercentStance)), colour="seagreen", linetype="solid", cex=1, alpha=0.25)+
  geom_vline(aes(xintercept=mean(Pec.Pw.Kept$PercentStance)), colour="lightblue3", linetype="dashed", cex=1, alpha=0.25)+
  geom_vline(aes(xintercept=mean(Pel.Pw.Kept$PercentStance)), colour="steelblue3", linetype="dotted", cex=1, alpha=0.25)+
  geom_vline(aes(xintercept=mean(Pec.Sl.Kept$PercentStance)), colour="purple", linetype="dashed", cex=1, alpha=0.25)+
  geom_line(size=1, alpha=0.75)+
  geom_ribbon(APAngle.MaxMin, alpha=0.5)+
  scale_colour_manual(values=c("ivory4", "ivory4", "ivory4", "ivory4", "ivory4"))+
  #scale_fill_manual(values=c("ivory4", "ivory4", "ivory4"))+
  scale_fill_manual(values=c("lightblue3", "steelblue3", "burlywood", "rosybrown4", "seagreen", "purple"))+ # make sure types are alphabetically arranged
  scale_linetype_manual(values=c("dashed", "dotted", "dashed", "dotted", "solid", "dashed"))+
  theme(axis.title.x=element_text(colour="black"))+ # vjust=0 puts a little more spacing btwn the axis text and label
  theme(axis.title.y=element_text(colour='black'))+
  theme(axis.text.x=element_text(colour='black'))+
  theme(axis.text.y=element_text(colour='black'))+
  theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank())+ # get rid of gridlines
  theme(panel.background=element_blank())+ # make background white
  theme(axis.line=element_line(colour="black", linetype="solid"))+ # put black lines for axes
  theme(legend.position="none")+
  #ggtitle("FL vs. HL: 109.75° |  FL vs. PF: 108.97° |  HL vs. PF: 24.09°\n")+
  #geom_hline(aes(yintercept=0, colour="aliceblue"))+
  theme(plot.title=element_text(size=8))

#theme(legend.key.width=unit(5,"lines")) # changes width of legend

MLAngle.Plot <- ggplot(data=MLAngle.Mean.SE, aes(x=MLAngle.Mean.SE$Stance, y=MLAngle.Mean.SE$Mean, fill=Type, linetype=Type))+
  scale_y_continuous("Mediolateral Angle (degrees)\n")+
  #scale_x_continuous("\nStance (%)\n")+
  scale_x_continuous(element_blank())+
  geom_rect(aes(xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=0), fill="gray95")+
  geom_vline(aes(xintercept=mean(Pec.Af.Kept$PercentStance)), colour="burlywood", linetype="dashed", cex=1, alpha=0.25)+
  geom_vline(aes(xintercept=mean(Pel.Af.Kept$PercentStance)), colour="rosybrown4", linetype="dotted", cex=1, alpha=0.25)+
  geom_vline(aes(xintercept=mean(Pec.Pb.Kept$PercentStance)), colour="seagreen", linetype="solid", cex=1, alpha=0.25)+
  geom_vline(aes(xintercept=mean(Pec.Pw.Kept$PercentStance)), colour="lightblue3", linetype="dashed", cex=1, alpha=0.25)+
  geom_vline(aes(xintercept=mean(Pel.Pw.Kept$PercentStance)), colour="steelblue3", linetype="dotted", cex=1, alpha=0.25)+
  geom_line(size=1, alpha=0.75)+
  geom_ribbon(MLAngle.MaxMin, alpha=0.5)+
  scale_colour_manual(values=c("ivory4", "ivory4", "ivory4", "ivory4", "ivory4"))+
  #scale_fill_manual(values=c("ivory4", "ivory4", "ivory4"))+
  scale_fill_manual(values=c("lightblue3", "steelblue3", "burlywood", "rosybrown4", "seagreen"))+
  scale_linetype_manual(values=c("dashed", "dotted", "dashed", "dotted", "solid"))+
  theme(axis.title.x=element_text(colour="black"))+ # vjust=0 puts a little more spacing btwn the axis text and label
  theme(axis.title.y=element_text(colour='black'))+
  theme(axis.text.x=element_text(colour='black'))+
  theme(axis.text.y=element_text(colour='black'))+
  theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank())+ # get rid of gridlines
  theme(panel.background=element_blank())+ # make background white
  theme(axis.line=element_line(colour="black", linetype="solid"))+ # put black lines for axes
  theme(legend.position="none")+
  #ggtitle("FL vs. HL: 21.41° |  FL vs. PF: 27.33° |  HL vs. PF: 38.88°\n")+
  #geom_hline(aes(yintercept=0, colour="ivory3"))+
  theme(plot.title=element_text(size=8))


VBW.Plot <- ggplot(data=InterpV.Mean.SE, aes(x=InterpV.Mean.SE$Stance, y=InterpV.Mean.SE$Mean, fill=Type, linetype=Type))+
  scale_y_continuous("Vertical GRF (BW)\n", limits=c(0.05,0.5))+
  #scale_x_continuous("\nStance (%)\n")+
  scale_x_continuous(element_blank())+
  geom_rect(aes(xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=0), fill="gray95")+
  geom_vline(aes(xintercept=mean(Pec.Af.Kept$PercentStance)), colour="burlywood", linetype="dashed", cex=1, alpha=0.25)+
  geom_vline(aes(xintercept=mean(Pel.Af.Kept$PercentStance)), colour="rosybrown4", linetype="dotted", cex=1, alpha=0.25)+
  geom_vline(aes(xintercept=mean(Pec.Pb.Kept$PercentStance)), colour="seagreen", linetype="solid", cex=1, alpha=0.25)+
  geom_vline(aes(xintercept=mean(Pec.Pw.Kept$PercentStance)), colour="lightblue3", linetype="dashed", cex=1, alpha=0.25)+
  geom_vline(aes(xintercept=mean(Pel.Pw.Kept$PercentStance)), colour="steelblue3", linetype="dotted", cex=1, alpha=0.25)+
  geom_line(size=1, alpha=0.75)+
  geom_ribbon(InterpV.MaxMin, alpha=0.5)+
  scale_colour_manual(values=c("ivory4", "ivory4", "ivory4", "ivory4", "ivory4"))+
  #scale_fill_manual(values=c("ivory4", "ivory4", "ivory4"))+
  scale_fill_manual(values=c("lightblue3", "steelblue3", "burlywood", "rosybrown4", "seagreen"))+
  scale_linetype_manual(values=c("dashed", "dotted", "dashed", "dotted", "solid"))+
  theme(axis.title.x=element_text(colour="black"))+ # vjust=0 puts a little more spacing btwn the axis text and label
  theme(axis.title.y=element_text(colour='black'))+
  theme(axis.text.x=element_text(colour='black'))+
  theme(axis.text.y=element_text(colour='black'))+
  theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank())+ # get rid of gridlines
  theme(panel.background=element_blank())+ # make background white
  theme(axis.line=element_line(colour="black", linetype="solid"))+ # put black lines for axes
  theme(legend.position="none")+
  #ggtitle("FL vs. HL: 25.79° |  FL vs. PF: 9.72° | HL vs. PF: 16.28°\n")+
  #geom_hline(aes(yintercept=0, colour="ivory3"))+
  theme(plot.title=element_text(size=8))

MLBW.Plot <- ggplot(data=InterpML.Mean.SE, aes(x=InterpML.Mean.SE$Stance, y=InterpML.Mean.SE$Mean, fill=Type, linetype=Type))+
  scale_y_continuous("Mediolateral GRF (BW)\n")+
  #scale_x_continuous("\nStance (%)\n")+
  scale_x_continuous(element_blank())+
  geom_rect(aes(xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=0), fill="gray95")+
  geom_vline(aes(xintercept=mean(Pec.Af.Kept$PercentStance)), colour="burlywood", linetype="dashed", cex=1, alpha=0.25)+
  geom_vline(aes(xintercept=mean(Pel.Af.Kept$PercentStance)), colour="rosybrown4", linetype="dotted", cex=1, alpha=0.25)+
  geom_vline(aes(xintercept=mean(Pec.Pb.Kept$PercentStance)), colour="seagreen", linetype="solid", cex=1, alpha=0.25)+
  geom_vline(aes(xintercept=mean(Pec.Pw.Kept$PercentStance)), colour="lightblue3", linetype="dashed", cex=1, alpha=0.25)+
  geom_vline(aes(xintercept=mean(Pel.Pw.Kept$PercentStance)), colour="steelblue3", linetype="dotted", cex=1, alpha=0.25)+
  geom_line(size=1, alpha=0.75)+
  geom_ribbon(InterpML.MaxMin, alpha=0.5)+
  scale_colour_manual(values=c("ivory4", "ivory4", "ivory4", "ivory4", "ivory4"))+
  #scale_fill_manual(values=c("ivory4", "ivory4", "ivory4"))+
  scale_fill_manual(values=c("lightblue3", "steelblue3", "burlywood", "rosybrown4", "seagreen"))+
  scale_linetype_manual(values=c("dashed", "dotted", "dashed", "dotted", "solid"))+
  theme(axis.title.x=element_text(colour="black"))+ # vjust=0 puts a little more spacing btwn the axis text and label
  theme(axis.title.y=element_text(colour='black'))+
  theme(axis.text.x=element_text(colour='black'))+
  theme(axis.text.y=element_text(colour='black'))+
  theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank())+ # get rid of gridlines
  theme(panel.background=element_blank())+ # make background white
  theme(axis.line=element_line(colour="black", linetype="solid"))+ # put black lines for axes
  theme(legend.position="none")+
  #ggtitle("FL vs. HL: 28.18° |  FL vs. PF: 14.89° | HL vs. PF: 40.05°\n")+
  #geom_hline(aes(yintercept=0, colour="ivory3"))+
  theme(plot.title=element_text(size=8))

APBW.Plot <- ggplot(data=InterpHz.Mean.SE, aes(x=InterpHz.Mean.SE$Stance, y=InterpHz.Mean.SE$Mean, fill=Type, linetype=Type))+
  scale_y_continuous("Anteroposterior GRF (BW)\n", breaks=c(-0.1,0, 0.1,0.2), limits=c(-0.15,0.20))+
  #scale_x_continuous("\nStance (%)\n")+
  scale_x_continuous(element_blank())+
  geom_rect(aes(xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=0), fill="gray95")+
  geom_vline(aes(xintercept=mean(Pec.Af.Kept$PercentStance)), colour="burlywood", linetype="dashed", cex=1, alpha=0.25)+
  geom_vline(aes(xintercept=mean(Pel.Af.Kept$PercentStance)), colour="rosybrown4", linetype="dotted", cex=1, alpha=0.25)+
  geom_vline(aes(xintercept=mean(Pec.Pb.Kept$PercentStance)), colour="seagreen", linetype="solid", cex=1, alpha=0.25)+
  geom_vline(aes(xintercept=mean(Pec.Pw.Kept$PercentStance)), colour="lightblue3", linetype="dashed", cex=1, alpha=0.25)+
  geom_vline(aes(xintercept=mean(Pel.Pw.Kept$PercentStance)), colour="steelblue3", linetype="dotted", cex=1, alpha=0.25)+
  geom_line(size=1, alpha=0.75)+
  geom_ribbon(InterpHz.MaxMin, alpha=0.5)+
  scale_colour_manual(values=c("ivory4", "ivory4", "ivory4", "ivory4", "ivory4"))+
  #scale_fill_manual(values=c("ivory4", "ivory4", "ivory4"))+
  scale_fill_manual(values=c("lightblue3", "steelblue3", "burlywood", "rosybrown4", "seagreen"))+
  scale_linetype_manual(values=c("dashed", "dotted", "dashed", "dotted", "solid"))+
  theme(axis.title.x=element_text(colour="black"))+ # vjust=0 puts a little more spacing btwn the axis text and label
  theme(axis.title.y=element_text(colour='black'))+
  theme(axis.text.x=element_text(colour='black'))+
  theme(axis.text.y=element_text(colour='black'))+
  theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank())+ # get rid of gridlines
  theme(panel.background=element_blank())+ # make background white
  theme(axis.line=element_line(colour="black", linetype="solid"))+ # put black lines for axes
  theme(legend.position="none")+
  #ggtitle("FL vs. HL: 117.18° | FL vs. PF: 114.47° | HL vs. PF: 19.28°\n")+
  #geom_hline(aes(yintercept=0, colour="ivory3"))+
  theme(plot.title=element_text(size=8))

NetBW.Plot <- ggplot(data=NetGRF.Mean.SE, aes(x=NetGRF.Mean.SE$Stance, y=NetGRF.Mean.SE$Mean, fill=Type, linetype=Type))+
  scale_y_continuous("Net GRF (BW)\n", limits=c(0.05,0.5))+
  #scale_x_continuous("\nStance (%)\n")+
  scale_x_continuous(element_blank())+
  geom_rect(aes(xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=0), fill="gray95")+
  geom_vline(aes(xintercept=mean(Pec.Af.Kept$PercentStance)), colour="burlywood", linetype="dashed", cex=1, alpha=0.25)+
  geom_vline(aes(xintercept=mean(Pel.Af.Kept$PercentStance)), colour="rosybrown4", linetype="dotted", cex=1, alpha=0.25)+
  geom_vline(aes(xintercept=mean(Pec.Pb.Kept$PercentStance)), colour="seagreen", linetype="solid", cex=1, alpha=0.25)+
  geom_vline(aes(xintercept=mean(Pec.Pw.Kept$PercentStance)), colour="lightblue3", linetype="dashed", cex=1, alpha=0.25)+
  geom_vline(aes(xintercept=mean(Pel.Pw.Kept$PercentStance)), colour="steelblue3", linetype="dotted", cex=1, alpha=0.25)+
  geom_line(size=1, alpha=0.75)+
  geom_ribbon(NetGRF.MaxMin, alpha=0.5)+
  scale_colour_manual(name="Appendage:", # changing legend title
                      labels=c("Semi-Aquatic FL  ", "Semi-Aquatic HL  ", "Terrestrial FL  ", "Terrestrial HL  ", "Terrestrial PF  "), # Changing legend labels
                      values=c("ivory4", "ivory4", "ivory4", "ivory4", "ivory4"))+
  scale_fill_manual(name="Appendage:", 
                    labels=c("Semi-Aquatic FL  ", "Semi-Aquatic HL  ", "Terrestrial FL  ", "Terrestrial HL  ", "Terrestrial PF  "),
                    #values=c("ivory4", "ivory4", "ivory4"))+
                    values=c("lightblue3", "steelblue3", "burlywood", "rosybrown4", "seagreen"))+
  scale_linetype_manual(name="Appendage:", 
                        labels=c("Semi-Aquatic FL  ", "Semi-Aquatic HL  ", "Terrestrial FL  ", "Terrestrial HL  ", "Terrestrial PF  "),
                        values=c("dashed", "dotted", "dashed", "dotted", "solid"))+
  theme(axis.title.x=element_text(colour="black"))+ # vjust=0 puts a little more spacing btwn the axis text and label
  theme(axis.title.y=element_text(colour='black'))+
  theme(axis.text.x=element_text(colour='black'))+
  theme(axis.text.y=element_text(colour='black'))+
  theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank())+ # get rid of gridlines
  theme(panel.background=element_blank())+ # make background white
  theme(axis.line=element_line(colour="black", linetype="solid"))+ # put black lines for axes
  #theme(legend.position="none")+
  #ggtitle("FL vs. HL: 24.27° |  FL vs. PF: 8.97° | HL vs. PF: 15.58°\n")+
  #geom_hline(aes(yintercept=0))+
  theme(legend.position="bottom", legend.direction="horizontal")+
  theme(plot.title=element_text(size=8))


## Extracting the legend
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

mylegend<-g_legend(NetBW.Plot)

library(gridExtra) # for grid.arrange
# Putting all the profile plots onto a single frame, and saving as a tiff file
tiff(filename=paste("AllPlots.Compare_", SaveDate, ".tif", sep=""), width=8.08661*300, height=9.5*300, res=300)
AllPlots <- grid.arrange(arrangeGrob((NetBW.Plot + theme(legend.position="none")), VBW.Plot, APBW.Plot, MLBW.Plot, APAngle.Plot, MLAngle.Plot, sub="Stance (%)"), mylegend, nrow=2,heights=c(10, 1))
dev.off()

# PdfSave <- paste("Appendage Force Comparison_", SaveDate, ".pdf", sep="")
# pdf(PdfSave)


## AP Angle 
Af.APAngle.Pec.Mean.SE <- data.frame(Af.APAngle.Pec.Mean, Af.GRF.APAngle.Convert.Combined.Pec.SE, Stance)
Af.APAngle.Pec.Mean.SE$Type <- "Terrestrial FL"
names(Af.APAngle.Pec.Mean.SE) <- c("Mean", "SE", "Stance", "Type")
Af.APAngle.Pel.Mean.SE <- data.frame(Af.APAngle.Pel.Mean, Af.GRF.APAngle.Convert.Combined.Pel.SE, Stance)
Af.APAngle.Pel.Mean.SE$Type <- "Terrestrial HL"
names(Af.APAngle.Pel.Mean.SE) <- c("Mean", "SE", "Stance", "Type")
Pb.APAngle.Pec.Mean.SE <- data.frame(Pb.APAngle.Pec.Mean, Pb.GRF.APAngle.Convert.Combined.Pec.SE, Stance)
Pb.APAngle.Pec.Mean.SE$Type <- "Terrestrial PF"
names(Pb.APAngle.Pec.Mean.SE) <- c("Mean", "SE", "Stance", "Type")
Pw.APAngle.Pec.Mean.SE <- data.frame(Pw.APAngle.Pec.Mean, Pw.GRF.APAngle.Convert.Combined.Pec.SE, Stance)
Pw.APAngle.Pec.Mean.SE$Type <- "Aquatic FL"
names(Pw.APAngle.Pec.Mean.SE) <- c("Mean", "SE", "Stance", "Type")
Pw.APAngle.Pel.Mean.SE <- data.frame(Pw.APAngle.Pel.Mean, Pw.GRF.APAngle.Convert.Combined.Pel.SE, Stance)
Pw.APAngle.Pel.Mean.SE$Type <- "Aquatic HL"
names(Pw.APAngle.Pel.Mean.SE) <- c("Mean", "SE", "Stance", "Type")


APAngle.Mean.SE <- rbind(Af.APAngle.Pec.Mean.SE, Af.APAngle.Pel.Mean.SE, Pb.APAngle.Pec.Mean.SE, Pw.APAngle.Pec.Mean.SE, Pw.APAngle.Pel.Mean.SE)

APAngle.MaxMin <- aes(ymax=APAngle.Mean.SE$Mean+APAngle.Mean.SE$SE, ymin=APAngle.Mean.SE$Mean-APAngle.Mean.SE$SE)

ggplot(data=APAngle.Mean.SE, aes(x=APAngle.Mean.SE$Stance, y=APAngle.Mean.SE$Mean, fill=Type))+
  scale_y_continuous("Anteroposterior Angle (degrees)")+
  scale_x_continuous("Stance")+
  geom_line(size=2, alpha=0.75)+
  geom_ribbon(APAngle.MaxMin, alpha=0.5)



# ML Angle
Af.MLAngle.Pec.Mean.SE <- data.frame(Af.MLAngle.Pec.Mean, Af.GRF.MLAngle.Convert.Combined.Pec.SE, Stance)
Af.MLAngle.Pec.Mean.SE$Type <- "Terrestrial FL"
names(Af.MLAngle.Pec.Mean.SE) <- c("Mean", "SE", "Stance", "Type")
Af.MLAngle.Pel.Mean.SE <- data.frame(Af.MLAngle.Pel.Mean, Af.GRF.MLAngle.Convert.Combined.Pel.SE, Stance)
Af.MLAngle.Pel.Mean.SE$Type <- "Terrestrial HL"
names(Af.MLAngle.Pel.Mean.SE) <- c("Mean", "SE", "Stance", "Type")
Pb.MLAngle.Pec.Mean.SE <- data.frame(Pb.MLAngle.Pec.Mean, Pb.GRF.MLAngle.Convert.Combined.Pec.SE, Stance)
Pb.MLAngle.Pec.Mean.SE$Type <- "Terrestrial PF"
names(Pb.MLAngle.Pec.Mean.SE) <- c("Mean", "SE", "Stance", "Type")
Pw.MLAngle.Pec.Mean.SE <- data.frame(Pw.MLAngle.Pec.Mean, Pw.GRF.MLAngle.Convert.Combined.Pec.SE, Stance)
Pw.MLAngle.Pec.Mean.SE$Type <- "Aquatic FL"
names(Pw.MLAngle.Pec.Mean.SE) <- c("Mean", "SE", "Stance", "Type")
Pw.MLAngle.Pel.Mean.SE <- data.frame(Pw.MLAngle.Pel.Mean, Pw.GRF.MLAngle.Convert.Combined.Pel.SE, Stance)
Pw.MLAngle.Pel.Mean.SE$Type <- "Aquatic HL"
names(Pw.MLAngle.Pel.Mean.SE) <- c("Mean", "SE", "Stance", "Type")


MLAngle.Mean.SE <- rbind(Af.MLAngle.Pec.Mean.SE, Af.MLAngle.Pel.Mean.SE, Pb.MLAngle.Pec.Mean.SE, Pw.MLAngle.Pec.Mean.SE, Pw.MLAngle.Pel.Mean.SE)

MLAngle.MaxMin <- aes(ymax=MLAngle.Mean.SE$Mean+MLAngle.Mean.SE$SE, ymin=MLAngle.Mean.SE$Mean-MLAngle.Mean.SE$SE)

ggplot(data=MLAngle.Mean.SE, aes(x=MLAngle.Mean.SE$Stance, y=MLAngle.Mean.SE$Mean, fill=Type))+
  scale_y_continuous("Mediolateral Angle (degrees)")+
  scale_x_continuous("Stance")+
  geom_line(size=2, alpha=0.75)+
  geom_ribbon(MLAngle.MaxMin, alpha=0.5)


#Vertical GRF
Af.InterpV.Pec.Mean.SE <- data.frame(Af.InterpV.Pec.Mean, Af.GRF.InterpV.BW.Combined.Pec.SE, Stance)
Af.InterpV.Pec.Mean.SE$Type <- "Terrestrial FL"
names(Af.InterpV.Pec.Mean.SE) <- c("Mean", "SE", "Stance", "Type")
Af.InterpV.Pel.Mean.SE <- data.frame(Af.InterpV.Pel.Mean, Af.GRF.InterpV.BW.Combined.Pel.SE, Stance)
Af.InterpV.Pel.Mean.SE$Type <- "Terrestrial HL"
names(Af.InterpV.Pel.Mean.SE) <- c("Mean", "SE", "Stance", "Type")
Pb.InterpV.Pec.Mean.SE <- data.frame(Pb.InterpV.Pec.Mean, Pb.GRF.InterpV.BW.Combined.Pec.SE, Stance)
Pb.InterpV.Pec.Mean.SE$Type <- "Terrestrial PF"
names(Pb.InterpV.Pec.Mean.SE) <- c("Mean", "SE", "Stance", "Type")
Pw.InterpV.Pec.Mean.SE <- data.frame(Pw.InterpV.Pec.Mean, Pw.GRF.InterpV.BW.Combined.Pec.SE, Stance)
Pw.InterpV.Pec.Mean.SE$Type <- "Aquatic FL"
names(Pw.InterpV.Pec.Mean.SE) <- c("Mean", "SE", "Stance", "Type")
Pw.InterpV.Pel.Mean.SE <- data.frame(Pw.InterpV.Pel.Mean, Pw.GRF.InterpV.BW.Combined.Pel.SE, Stance)
Pw.InterpV.Pel.Mean.SE$Type <- "Aquatic HL"
names(Pw.InterpV.Pel.Mean.SE) <- c("Mean", "SE", "Stance", "Type")

InterpV.Mean.SE <- rbind(Af.InterpV.Pec.Mean.SE, Af.InterpV.Pel.Mean.SE, Pb.InterpV.Pec.Mean.SE, Pw.InterpV.Pec.Mean.SE, Pw.InterpV.Pel.Mean.SE)

InterpV.MaxMin <- aes(ymax=InterpV.Mean.SE$Mean+InterpV.Mean.SE$SE, ymin=InterpV.Mean.SE$Mean-InterpV.Mean.SE$SE)

ggplot(data=InterpV.Mean.SE, aes(x=InterpV.Mean.SE$Stance, y=InterpV.Mean.SE$Mean, fill=Type))+
  scale_y_continuous("Vertical (BW)")+
  scale_x_continuous("Stance")+
  geom_line(size=2, alpha=0.75)+
  geom_ribbon(InterpV.MaxMin, alpha=0.5)



#Mediolateral GRF
Af.InterpML.Pec.Mean.SE <- data.frame(Af.InterpML.Pec.Mean, Af.GRF.InterpML.BW.Combined.Pec.SE, Stance)
Af.InterpML.Pec.Mean.SE$Type <- "Terrestrial FL"
names(Af.InterpML.Pec.Mean.SE) <- c("Mean", "SE", "Stance", "Type")
Af.InterpML.Pel.Mean.SE <- data.frame(Af.InterpML.Pel.Mean, Af.GRF.InterpML.BW.Combined.Pel.SE, Stance)
Af.InterpML.Pel.Mean.SE$Type <- "Terrestrial HL"
names(Af.InterpML.Pel.Mean.SE) <- c("Mean", "SE", "Stance", "Type")
Pb.InterpML.Pec.Mean.SE <- data.frame(Pb.InterpML.Pec.Mean, Pb.GRF.InterpML.BW.Combined.Pec.SE, Stance)
Pb.InterpML.Pec.Mean.SE$Type <- "Terrestrial PF"
names(Pb.InterpML.Pec.Mean.SE) <- c("Mean", "SE", "Stance", "Type")
Pw.InterpML.Pec.Mean.SE <- data.frame(Pw.InterpML.Pec.Mean, Pw.GRF.InterpML.BW.Combined.Pec.SE, Stance)
Pw.InterpML.Pec.Mean.SE$Type <- "Aquatic FL"
names(Pw.InterpML.Pec.Mean.SE) <- c("Mean", "SE", "Stance", "Type")
Pw.InterpML.Pel.Mean.SE <- data.frame(Pw.InterpML.Pel.Mean, Pw.GRF.InterpML.BW.Combined.Pel.SE, Stance)
Pw.InterpML.Pel.Mean.SE$Type <- "Aquatic HL"
names(Pw.InterpML.Pel.Mean.SE) <- c("Mean", "SE", "Stance", "Type")


InterpML.Mean.SE <- rbind(Af.InterpML.Pec.Mean.SE, Af.InterpML.Pel.Mean.SE, Pb.InterpML.Pec.Mean.SE, Pw.InterpML.Pec.Mean.SE, Pw.InterpML.Pel.Mean.SE)

InterpML.MaxMin <- aes(ymax=InterpML.Mean.SE$Mean+InterpML.Mean.SE$SE, ymin=InterpML.Mean.SE$Mean-InterpML.Mean.SE$SE)

ggplot(data=InterpML.Mean.SE, aes(x=InterpML.Mean.SE$Stance, y=InterpML.Mean.SE$Mean, fill=Type))+
  scale_y_continuous("Mediolateral (BW)")+
  scale_x_continuous("Stance")+
  geom_line(size=2, alpha=0.75)+
  geom_ribbon(InterpML.MaxMin, alpha=0.5)


#Anteroposterior GRF
Af.InterpHz.Pec.Mean.SE <- data.frame(Af.InterpHz.Pec.Mean, Af.GRF.InterpHz.BW.Combined.Pec.SE, Stance)
Af.InterpHz.Pec.Mean.SE$Type <- "Terrestrial FL"
names(Af.InterpHz.Pec.Mean.SE) <- c("Mean", "SE", "Stance", "Type")
Af.InterpHz.Pel.Mean.SE <- data.frame(Af.InterpHz.Pel.Mean, Af.GRF.InterpHz.BW.Combined.Pel.SE, Stance)
Af.InterpHz.Pel.Mean.SE$Type <- "Terrestrial HL"
names(Af.InterpHz.Pel.Mean.SE) <- c("Mean", "SE", "Stance", "Type")
Pb.InterpHz.Pec.Mean.SE <- data.frame(Pb.InterpHz.Pec.Mean, Pb.GRF.InterpHz.BW.Combined.Pec.SE, Stance)
Pb.InterpHz.Pec.Mean.SE$Type <- "Terrestrial PF"
names(Pb.InterpHz.Pec.Mean.SE) <- c("Mean", "SE", "Stance", "Type")
Pw.InterpHz.Pec.Mean.SE <- data.frame(Pw.InterpHz.Pec.Mean, Pw.GRF.InterpHz.BW.Combined.Pec.SE, Stance)
Pw.InterpHz.Pec.Mean.SE$Type <- "Aquatic FL"
names(Pw.InterpHz.Pec.Mean.SE) <- c("Mean", "SE", "Stance", "Type")
Pw.InterpHz.Pel.Mean.SE <- data.frame(Pw.InterpHz.Pel.Mean, Pw.GRF.InterpHz.BW.Combined.Pel.SE, Stance)
Pw.InterpHz.Pel.Mean.SE$Type <- "Aquatic HL"
names(Pw.InterpHz.Pel.Mean.SE) <- c("Mean", "SE", "Stance", "Type")


InterpHz.Mean.SE <- rbind(Af.InterpHz.Pec.Mean.SE, Af.InterpHz.Pel.Mean.SE, Pb.InterpHz.Pec.Mean.SE, Pw.InterpHz.Pec.Mean.SE, Pw.InterpHz.Pel.Mean.SE)

InterpHz.MaxMin <- aes(ymax=InterpHz.Mean.SE$Mean+InterpHz.Mean.SE$SE, ymin=InterpHz.Mean.SE$Mean-InterpHz.Mean.SE$SE)

ggplot(data=InterpHz.Mean.SE, aes(x=InterpHz.Mean.SE$Stance, y=InterpHz.Mean.SE$Mean, fill=Type))+
  scale_y_continuous("Anteroposterior (BW)")+
  scale_x_continuous("Stance")+
  geom_line(size=2, alpha=0.75)+
  geom_ribbon(InterpHz.MaxMin, alpha=0.5)

#Net GRF
Af.NetGRF.Pec.Mean.SE <- data.frame(Af.NetGRF.Pec.Mean, Af.GRF.NetGRF.BW.Combined.Pec.SE, Stance)
Af.NetGRF.Pec.Mean.SE$Type <- "Terrestrial FL"
names(Af.NetGRF.Pec.Mean.SE) <- c("Mean", "SE", "Stance", "Type")
Af.NetGRF.Pel.Mean.SE <- data.frame(Af.NetGRF.Pel.Mean, Af.GRF.NetGRF.BW.Combined.Pel.SE, Stance)
Af.NetGRF.Pel.Mean.SE$Type <- "Terrestrial HL"
names(Af.NetGRF.Pel.Mean.SE) <- c("Mean", "SE", "Stance", "Type")
Pb.NetGRF.Pec.Mean.SE <- data.frame(Pb.NetGRF.Pec.Mean, Pb.GRF.NetGRF.BW.Combined.Pec.SE, Stance)
Pb.NetGRF.Pec.Mean.SE$Type <- "Terrestrial PF"
names(Pb.NetGRF.Pec.Mean.SE) <- c("Mean", "SE", "Stance", "Type")
Pw.NetGRF.Pec.Mean.SE <- data.frame(Pw.NetGRF.Pec.Mean, Pw.GRF.NetGRF.BW.Combined.Pec.SE, Stance)
Pw.NetGRF.Pec.Mean.SE$Type <- "Aquatic FL"
names(Pw.NetGRF.Pec.Mean.SE) <- c("Mean", "SE", "Stance", "Type")
Pw.NetGRF.Pel.Mean.SE <- data.frame(Pw.NetGRF.Pel.Mean, Pw.GRF.NetGRF.BW.Combined.Pel.SE, Stance)
Pw.NetGRF.Pel.Mean.SE$Type <- "Aquatic HL"
names(Pw.NetGRF.Pel.Mean.SE) <- c("Mean", "SE", "Stance", "Type")


NetGRF.Mean.SE <- rbind(Af.NetGRF.Pec.Mean.SE, Af.NetGRF.Pel.Mean.SE, Pb.NetGRF.Pec.Mean.SE, Pw.NetGRF.Pec.Mean.SE, Pw.NetGRF.Pel.Mean.SE)

NetGRF.MaxMin <- aes(ymax=NetGRF.Mean.SE$Mean+NetGRF.Mean.SE$SE, ymin=NetGRF.Mean.SE$Mean-NetGRF.Mean.SE$SE)

ggplot(data=NetGRF.Mean.SE, aes(x=NetGRF.Mean.SE$Stance, y=NetGRF.Mean.SE$Mean, fill=Type))+
  scale_y_continuous("Net GRF (BW)")+
  scale_x_continuous("Stance")+
  geom_line(size=2, alpha=0.75)+
  geom_ribbon(NetGRF.MaxMin, alpha=0.5)

dev.off()


### HINDLIMB VS. FORELIMB COMPARISON ONLY
## AP Angle 
Af.APAngle.Pec.Mean.SE <- data.frame(Af.APAngle.Pec.Mean, Af.GRF.APAngle.Convert.Combined.Pec.SE, Stance)
Af.APAngle.Pec.Mean.SE$Type <- "Terrestrial FL"
names(Af.APAngle.Pec.Mean.SE) <- c("Mean", "SE", "Stance", "Type")
Af.APAngle.Pel.Mean.SE <- data.frame(Af.APAngle.Pel.Mean, Af.GRF.APAngle.Convert.Combined.Pel.SE, Stance)
Af.APAngle.Pel.Mean.SE$Type <- "Terrestrial HL"
names(Af.APAngle.Pel.Mean.SE) <- c("Mean", "SE", "Stance", "Type")
Pw.APAngle.Pec.Mean.SE <- data.frame(Pw.APAngle.Pec.Mean, Pw.GRF.APAngle.Convert.Combined.Pec.SE, Stance)
Pw.APAngle.Pec.Mean.SE$Type <- "Aquatic FL"
names(Pw.APAngle.Pec.Mean.SE) <- c("Mean", "SE", "Stance", "Type")
Pw.APAngle.Pel.Mean.SE <- data.frame(Pw.APAngle.Pel.Mean, Pw.GRF.APAngle.Convert.Combined.Pel.SE, Stance)
Pw.APAngle.Pel.Mean.SE$Type <- "Aquatic HL"
names(Pw.APAngle.Pel.Mean.SE) <- c("Mean", "SE", "Stance", "Type")


Limb.APAngle.Mean.SE <- rbind(Af.APAngle.Pec.Mean.SE, Af.APAngle.Pel.Mean.SE, Pw.APAngle.Pec.Mean.SE, Pw.APAngle.Pel.Mean.SE)

Limb.APAngle.MaxMin <- aes(ymax=Limb.APAngle.Mean.SE$Mean+Limb.APAngle.Mean.SE$SE, ymin=Limb.APAngle.Mean.SE$Mean-Limb.APAngle.Mean.SE$SE)


ggplot(data=Limb.APAngle.Mean.SE, aes(x=Limb.APAngle.Mean.SE$Stance, y=Limb.APAngle.Mean.SE$Mean, fill=Type, linetype=Type))+
  scale_y_continuous("Anteroposterior Angle (degrees)")+
  scale_x_continuous("Stance")+
  geom_line(size=2, alpha=0.75)+
  geom_ribbon(Limb.APAngle.MaxMin, alpha=0.5)+
  #scale_fill_manual(values=ggplotColours(n=3)) # Although there are only 2 groups to compare, n=2 makes the hindlimb green
  scale_fill_manual(name="Appendage:", 
                    labels=c("Terrestrial FL    ", "Terrestrial HL    ", "Aquatic FL  ", "Aquatic HL  "),
                    #values=c("ivory4", "ivory4", "ivory4")+
                    values=c("lightblue3", "firebrick3", "green", "brown"))+
  scale_linetype_manual(name="Appendage:", 
                        labels=c("Terrestrial FL    ", "Terrestrial HL    ", "Aquatic FL  ", "Aquatic HL  "),
                        values=c("dashed", "dotted", "dashed", "dotted"))
  

# ML Angle
Af.MLAngle.Pec.Mean.SE <- data.frame(Af.MLAngle.Pec.Mean, Af.GRF.MLAngle.Convert.Combined.Pec.SE, Stance)
Af.MLAngle.Pec.Mean.SE$Type <- "Terrestrial FL"
names(Af.MLAngle.Pec.Mean.SE) <- c("Mean", "SE", "Stance", "Type")
Af.MLAngle.Pel.Mean.SE <- data.frame(Af.MLAngle.Pel.Mean, Af.GRF.MLAngle.Convert.Combined.Pel.SE, Stance)
Af.MLAngle.Pel.Mean.SE$Type <- "Terrestrial HL"
names(Af.MLAngle.Pel.Mean.SE) <- c("Mean", "SE", "Stance", "Type")
Pw.MLAngle.Pec.Mean.SE <- data.frame(Pw.MLAngle.Pec.Mean, Pw.GRF.MLAngle.Convert.Combined.Pec.SE, Stance)
Pw.MLAngle.Pec.Mean.SE$Type <- "Aquatic FL"
names(Pw.MLAngle.Pec.Mean.SE) <- c("Mean", "SE", "Stance", "Type")
Pw.MLAngle.Pel.Mean.SE <- data.frame(Pw.MLAngle.Pel.Mean, Pw.GRF.MLAngle.Convert.Combined.Pel.SE, Stance)
Pw.MLAngle.Pel.Mean.SE$Type <- "Aquatic HL"
names(Pw.MLAngle.Pel.Mean.SE) <- c("Mean", "SE", "Stance", "Type")


Limb.MLAngle.Mean.SE <- rbind(Af.MLAngle.Pec.Mean.SE, Af.MLAngle.Pel.Mean.SE, Pw.MLAngle.Pec.Mean.SE, Pw.MLAngle.Pel.Mean.SE)

Limb.MLAngle.MaxMin <- aes(ymax=Limb.MLAngle.Mean.SE$Mean+Limb.MLAngle.Mean.SE$SE, ymin=Limb.MLAngle.Mean.SE$Mean-Limb.MLAngle.Mean.SE$SE)

ggplot(data=Limb.MLAngle.Mean.SE, aes(x=Limb.MLAngle.Mean.SE$Stance, y=Limb.MLAngle.Mean.SE$Mean, fill=Type, linetype=Type))+
  scale_y_continuous("Mediolateral Angle (degrees)")+
  scale_x_continuous("Stance")+
  geom_line(size=2, alpha=0.75)+
  geom_ribbon(Limb.MLAngle.MaxMin, alpha=0.5)+
  #scale_fill_manual(values=ggplotColours(n=3)) # Although there are only 2 groups to compare, n=2 makes the hindlimb green
  scale_fill_manual(name="Appendage:", 
                    labels=c("Terrestrial FL    ", "Terrestrial HL    ", "Aquatic FL  ", "Aquatic HL  "),
                    #values=c("ivory4", "ivory4", "ivory4")+
                    values=c("lightblue3", "firebrick3", "green", "brown"))+
  scale_linetype_manual(name="Appendage:", 
                        labels=c("Terrestrial FL    ", "Terrestrial HL    ", "Aquatic FL  ", "Aquatic HL  "),
                        values=c("dashed", "dotted", "dashed", "dotted"))

#Vertical GRF
Af.InterpV.Pec.Mean.SE <- data.frame(Af.InterpV.Pec.Mean, Af.GRF.InterpV.BW.Combined.Pec.SE, Stance)
Af.InterpV.Pec.Mean.SE$Type <- "Terrestrial FL"
names(Af.InterpV.Pec.Mean.SE) <- c("Mean", "SE", "Stance", "Type")
Af.InterpV.Pel.Mean.SE <- data.frame(Af.InterpV.Pel.Mean, Af.GRF.InterpV.BW.Combined.Pel.SE, Stance)
Af.InterpV.Pel.Mean.SE$Type <- "Terrestrial HL"
names(Af.InterpV.Pel.Mean.SE) <- c("Mean", "SE", "Stance", "Type")
Pw.InterpV.Pec.Mean.SE <- data.frame(Pw.InterpV.Pec.Mean, Pw.GRF.InterpV.BW.Combined.Pec.SE, Stance)
Pw.InterpV.Pec.Mean.SE$Type <- "Aquatic FL"
names(Pw.InterpV.Pec.Mean.SE) <- c("Mean", "SE", "Stance", "Type")
Pw.InterpV.Pel.Mean.SE <- data.frame(Pw.InterpV.Pel.Mean, Pw.GRF.InterpV.BW.Combined.Pel.SE, Stance)
Pw.InterpV.Pel.Mean.SE$Type <- "Aquatic HL"
names(Pw.InterpV.Pel.Mean.SE) <- c("Mean", "SE", "Stance", "Type")


Limb.InterpV.Mean.SE <- rbind(Af.InterpV.Pec.Mean.SE, Af.InterpV.Pel.Mean.SE, Pw.InterpV.Pec.Mean.SE, Pw.InterpV.Pel.Mean.SE)

Limb.InterpV.MaxMin <- aes(ymax=Limb.InterpV.Mean.SE$Mean+Limb.InterpV.Mean.SE$SE, ymin=Limb.InterpV.Mean.SE$Mean-Limb.InterpV.Mean.SE$SE)

ggplot(data=Limb.InterpV.Mean.SE, aes(x=Limb.InterpV.Mean.SE$Stance, y=Limb.InterpV.Mean.SE$Mean, fill=Type, linetype=Type))+
  scale_y_continuous("Vertical (BW)")+
  scale_x_continuous("Stance")+
  geom_line(size=2, alpha=0.75)+
  geom_ribbon(Limb.InterpV.MaxMin, alpha=0.5)+
  #scale_fill_manual(values=ggplotColours(n=3)) # Although there are only 2 groups to compare, n=2 makes the hindlimb green
  scale_fill_manual(name="Appendage:", 
                    labels=c("Terrestrial FL    ", "Terrestrial HL    ", "Aquatic FL  ", "Aquatic HL  " ),
                    #values=c("ivory4", "ivory4", "ivory4")+
                    values=c("lightblue3", "firebrick3", "green", "brown"))+
  scale_linetype_manual(name="Appendage:", 
                        labels=c("Terrestrial FL    ", "Terrestrial HL    ", "Aquatic FL  ", "Aquatic HL  "),
                        values=c("dashed", "dotted", "dashed", "dotted"))



#Mediolateral GRF
Af.InterpML.Pec.Mean.SE <- data.frame(Af.InterpML.Pec.Mean, Af.GRF.InterpML.BW.Combined.Pec.SE, Stance)
Af.InterpML.Pec.Mean.SE$Type <- "Terrestrial FL"
names(Af.InterpML.Pec.Mean.SE) <- c("Mean", "SE", "Stance", "Type")
Af.InterpML.Pel.Mean.SE <- data.frame(Af.InterpML.Pel.Mean, Af.GRF.InterpML.BW.Combined.Pel.SE, Stance)
Af.InterpML.Pel.Mean.SE$Type <- "Terrestrial HL"
names(Af.InterpML.Pel.Mean.SE) <- c("Mean", "SE", "Stance", "Type")
Pw.InterpML.Pec.Mean.SE <- data.frame(Pw.InterpML.Pec.Mean, Pw.GRF.InterpML.BW.Combined.Pec.SE, Stance)
Pw.InterpML.Pec.Mean.SE$Type <- "Aquatic FL"
names(Pw.InterpML.Pec.Mean.SE) <- c("Mean", "SE", "Stance", "Type")
Pw.InterpML.Pel.Mean.SE <- data.frame(Pw.InterpML.Pel.Mean, Pw.GRF.InterpML.BW.Combined.Pel.SE, Stance)
Pw.InterpML.Pel.Mean.SE$Type <- "Aquatic HL"
names(Pw.InterpML.Pel.Mean.SE) <- c("Mean", "SE", "Stance", "Type")


Limb.InterpML.Mean.SE <- rbind(Af.InterpML.Pec.Mean.SE, Af.InterpML.Pel.Mean.SE, Pw.InterpML.Pec.Mean.SE, Pw.InterpML.Pel.Mean.SE)

Limb.InterpML.MaxMin <- aes(ymax=Limb.InterpML.Mean.SE$Mean+Limb.InterpML.Mean.SE$SE, ymin=Limb.InterpML.Mean.SE$Mean-Limb.InterpML.Mean.SE$SE)

ggplot(data=Limb.InterpML.Mean.SE, aes(x=Limb.InterpML.Mean.SE$Stance, y=Limb.InterpML.Mean.SE$Mean, fill=Type, linetype=Type))+
  scale_y_continuous("Mediolateral (BW)")+
  scale_x_continuous("Stance")+
  geom_line(size=2, alpha=0.75)+
  geom_ribbon(Limb.InterpML.MaxMin, alpha=0.5)+
  #scale_fill_manual(values=ggplotColours(n=3)) # Although there are only 2 groups to compare, n=2 makes the hindlimb green
  scale_fill_manual(name="Appendage:", 
                    labels=c("Terrestrial FL    ", "Terrestrial HL    ", "Aquatic FL  ", "Aquatic HL  "),
                    #values=c("ivory4", "ivory4", "ivory4")+
                    values=c("lightblue3", "firebrick3", "green", "brown"))+
  scale_linetype_manual(name="Appendage:", 
                        labels=c("Terrestrial FL    ", "Terrestrial HL   ", "Aquatic FL  ", "Aquatic HL  "),
                        values=c("dashed", "dotted", "dashed", "dotted"))


#Anteroposterior GRF
Af.InterpHz.Pec.Mean.SE <- data.frame(Af.InterpHz.Pec.Mean, Af.GRF.InterpHz.BW.Combined.Pec.SE, Stance)
Af.InterpHz.Pec.Mean.SE$Type <- "Terrestrial FL"
names(Af.InterpHz.Pec.Mean.SE) <- c("Mean", "SE", "Stance", "Type")
Af.InterpHz.Pel.Mean.SE <- data.frame(Af.InterpHz.Pel.Mean, Af.GRF.InterpHz.BW.Combined.Pel.SE, Stance)
Af.InterpHz.Pel.Mean.SE$Type <- "Terrestrial HL"
names(Af.InterpHz.Pel.Mean.SE) <- c("Mean", "SE", "Stance", "Type")
Pw.InterpHz.Pec.Mean.SE <- data.frame(Pw.InterpHz.Pec.Mean, Pw.GRF.InterpHz.BW.Combined.Pec.SE, Stance)
Pw.InterpHz.Pec.Mean.SE$Type <- "Aquatic FL"
names(Pw.InterpHz.Pec.Mean.SE) <- c("Mean", "SE", "Stance", "Type")
Pw.InterpHz.Pel.Mean.SE <- data.frame(Pw.InterpHz.Pel.Mean, Pw.GRF.InterpHz.BW.Combined.Pel.SE, Stance)
Pw.InterpHz.Pel.Mean.SE$Type <- "Aquatic HL"
names(Pw.InterpHz.Pel.Mean.SE) <- c("Mean", "SE", "Stance", "Type")


Limb.InterpHz.Mean.SE <- rbind(Af.InterpHz.Pec.Mean.SE, Af.InterpHz.Pel.Mean.SE, Pw.InterpHz.Pec.Mean.SE, Pw.InterpHz.Pel.Mean.SE)

Limb.InterpHz.MaxMin <- aes(ymax=Limb.InterpHz.Mean.SE$Mean+Limb.InterpHz.Mean.SE$SE, ymin=Limb.InterpHz.Mean.SE$Mean-Limb.InterpHz.Mean.SE$SE)

ggplot(data=Limb.InterpHz.Mean.SE, aes(x=Limb.InterpHz.Mean.SE$Stance, y=Limb.InterpHz.Mean.SE$Mean, fill=Type, linetype=Type))+
  scale_y_continuous("Anteroposterior (BW)")+
  scale_x_continuous("Stance")+
  geom_line(size=2, alpha=0.75)+
  geom_ribbon(Limb.InterpHz.MaxMin, alpha=0.5)+
  #scale_fill_manual(values=ggplotColours(n=3)) # Although there are only 2 groups to compare, n=2 makes the hindlimb green
  scale_fill_manual(name="Appendage:", 
                    labels=c("Terrestrial FL    ", "Terrestrial HL    ", "Aquatic FL  ", "Aquatic HL  "),
                    #values=c("ivory4", "ivory4", "ivory4")+
                    values=c("lightblue3", "firebrick3", "green", "brown"))+
  scale_linetype_manual(name="Appendage:", 
                        labels=c("Terrestrial FL    ", "Terrestrial HL    ", "Aquatic FL  ", "Aquatic HL  "),
                        values=c("dashed", "dotted", "dashed", "dotted"))

#Net GRF
Af.NetGRF.Pec.Mean.SE <- data.frame(Af.NetGRF.Pec.Mean, Af.GRF.NetGRF.BW.Combined.Pec.SE, Stance)
Af.NetGRF.Pec.Mean.SE$Type <- "Terrestrial FL"
names(Af.NetGRF.Pec.Mean.SE) <- c("Mean", "SE", "Stance", "Type")
Af.NetGRF.Pel.Mean.SE <- data.frame(Af.NetGRF.Pel.Mean, Af.GRF.NetGRF.BW.Combined.Pel.SE, Stance)
Af.NetGRF.Pel.Mean.SE$Type <- "Terrestrial HL"
names(Af.NetGRF.Pel.Mean.SE) <- c("Mean", "SE", "Stance", "Type")
Pw.NetGRF.Pec.Mean.SE <- data.frame(Pw.NetGRF.Pec.Mean, Pw.GRF.NetGRF.BW.Combined.Pec.SE, Stance)
Pw.NetGRF.Pec.Mean.SE$Type <- "Aquatic FL"
names(Pw.NetGRF.Pec.Mean.SE) <- c("Mean", "SE", "Stance", "Type")
Pw.NetGRF.Pel.Mean.SE <- data.frame(Pw.NetGRF.Pel.Mean, Pw.GRF.NetGRF.BW.Combined.Pel.SE, Stance)
Pw.NetGRF.Pel.Mean.SE$Type <- "Aquatic HL"
names(Pw.NetGRF.Pel.Mean.SE) <- c("Mean", "SE", "Stance", "Type")


Limb.NetGRF.Mean.SE <- rbind(Af.NetGRF.Pec.Mean.SE, Af.NetGRF.Pel.Mean.SE, Pw.NetGRF.Pec.Mean.SE, Pw.NetGRF.Pel.Mean.SE)

Limb.NetGRF.MaxMin <- aes(ymax=Limb.NetGRF.Mean.SE$Mean+Limb.NetGRF.Mean.SE$SE, ymin=Limb.NetGRF.Mean.SE$Mean-Limb.NetGRF.Mean.SE$SE)


tiff(filename=paste("Limb Net GRF_", SaveDate, ".tif", sep=""), width=9*600, height=5*600, res=600)
ggplot(data=Limb.NetGRF.Mean.SE, aes(x=Limb.NetGRF.Mean.SE$Stance, y=Limb.NetGRF.Mean.SE$Mean, fill=Type, linetype=Type))+
  scale_y_continuous("Net GRF (BW)")+
  scale_x_continuous("Stance")+
  geom_line(size=2, alpha=0.75)+
  geom_ribbon(Limb.NetGRF.MaxMin, alpha=0.5)+
  #scale_fill_manual(values=ggplotColours(n=3)) # Although there are only 2 groups to compare, n=2 makes the hindlimb green
  scale_fill_manual(name="Appendage:", 
                    labels=c("Terrestrial FL    ", "Terrestrial HL    ", "Aquatic FL  ", "Aquatic HL  "),
                    #values=c("ivory4", "ivory4", "ivory4")+
                    values=c("lightblue3", "firebrick3", "green", "brown"))+
  scale_linetype_manual(name="Appendage:", 
                        labels=c("Terrestrial FL    ", "Terrestrial HL    ", "Aquatic FL  ", "Aquatic HL  "),
                        values=c("dashed", "dotted", "dashed", "dotted"))
dev.off()



PdfSave.Limb <- paste("Limb Force Comparison_", SaveDate, ".pdf", sep="")
pdf(PdfSave.Limb)

## AP Angle 
Af.APAngle.Pec.Mean.SE <- data.frame(Af.APAngle.Pec.Mean, Af.GRF.APAngle.Convert.Combined.Pec.SE, Stance)
Af.APAngle.Pec.Mean.SE$Type <- "Terrestrial FL"
names(Af.APAngle.Pec.Mean.SE) <- c("Mean", "SE", "Stance", "Type")
Af.APAngle.Pel.Mean.SE <- data.frame(Af.APAngle.Pel.Mean, Af.GRF.APAngle.Convert.Combined.Pel.SE, Stance)
Af.APAngle.Pel.Mean.SE$Type <- "Terrestrial HL"
names(Af.APAngle.Pel.Mean.SE) <- c("Mean", "SE", "Stance", "Type")
Pw.APAngle.Pec.Mean.SE <- data.frame(Pw.APAngle.Pec.Mean, Pw.GRF.APAngle.Convert.Combined.Pec.SE, Stance)
Pw.APAngle.Pec.Mean.SE$Type <- "Aquatic FL"
names(Pw.APAngle.Pec.Mean.SE) <- c("Mean", "SE", "Stance", "Type")
Pw.APAngle.Pel.Mean.SE <- data.frame(Pw.APAngle.Pel.Mean, Pw.GRF.APAngle.Convert.Combined.Pel.SE, Stance)
Pw.APAngle.Pel.Mean.SE$Type <- "Aquatic HL"
names(Pw.APAngle.Pel.Mean.SE) <- c("Mean", "SE", "Stance", "Type")


Limb.APAngle.Mean.SE <- rbind(Af.APAngle.Pec.Mean.SE, Af.APAngle.Pel.Mean.SE, Pw.APAngle.Pec.Mean.SE, Pw.APAngle.Pel.Mean.SE)

Limb.APAngle.MaxMin <- aes(ymax=Limb.APAngle.Mean.SE$Mean+Limb.APAngle.Mean.SE$SE, ymin=Limb.APAngle.Mean.SE$Mean-Limb.APAngle.Mean.SE$SE)


ggplot(data=Limb.APAngle.Mean.SE, aes(x=Limb.APAngle.Mean.SE$Stance, y=Limb.APAngle.Mean.SE$Mean, fill=Type, linetype=Type))+
  scale_y_continuous("Anteroposterior Angle (degrees)")+
  scale_x_continuous("Stance")+
  geom_line(size=2, alpha=0.75)+
  geom_ribbon(Limb.APAngle.MaxMin, alpha=0.5)+
  #scale_fill_manual(values=ggplotColours(n=3)) # Although there are only 2 groups to compare, n=2 makes the hindlimb green
  scale_fill_manual(name="Appendage:", 
                    labels=c("Terrestrial FL    ", "Terrestrial HL    ", "Aquatic FL  ", "Aquatic HL  "),
                    #values=c("ivory4", "ivory4", "ivory4")+
                    values=c("lightblue3", "firebrick3", "green", "brown"))+
  scale_linetype_manual(name="Appendage:", 
                        labels=c("Terrestrial FL    ", "Terrestrial HL    ", "Aquatic FL  ", "Aquatic HL  "),
                        values=c("dashed", "dotted", "dashed", "dotted"))



# ML Angle
Af.MLAngle.Pec.Mean.SE <- data.frame(Af.MLAngle.Pec.Mean, Af.GRF.MLAngle.Convert.Combined.Pec.SE, Stance)
Af.MLAngle.Pec.Mean.SE$Type <- "Terrestrial FL"
names(Af.MLAngle.Pec.Mean.SE) <- c("Mean", "SE", "Stance", "Type")
Af.MLAngle.Pel.Mean.SE <- data.frame(Af.MLAngle.Pel.Mean, Af.GRF.MLAngle.Convert.Combined.Pel.SE, Stance)
Af.MLAngle.Pel.Mean.SE$Type <- "Terrestrial HL"
names(Af.MLAngle.Pel.Mean.SE) <- c("Mean", "SE", "Stance", "Type")
Pw.MLAngle.Pec.Mean.SE <- data.frame(Pw.MLAngle.Pec.Mean, Pw.GRF.MLAngle.Convert.Combined.Pec.SE, Stance)
Pw.MLAngle.Pec.Mean.SE$Type <- "Aquatic FL"
names(Pw.MLAngle.Pec.Mean.SE) <- c("Mean", "SE", "Stance", "Type")
Pw.MLAngle.Pel.Mean.SE <- data.frame(Pw.MLAngle.Pel.Mean, Pw.GRF.MLAngle.Convert.Combined.Pel.SE, Stance)
Pw.MLAngle.Pel.Mean.SE$Type <- "Aquatic HL"
names(Pw.MLAngle.Pel.Mean.SE) <- c("Mean", "SE", "Stance", "Type")


Limb.MLAngle.Mean.SE <- rbind(Af.MLAngle.Pec.Mean.SE, Af.MLAngle.Pel.Mean.SE, Pw.MLAngle.Pec.Mean.SE, Pw.MLAngle.Pel.Mean.SE)

Limb.MLAngle.MaxMin <- aes(ymax=Limb.MLAngle.Mean.SE$Mean+Limb.MLAngle.Mean.SE$SE, ymin=Limb.MLAngle.Mean.SE$Mean-Limb.MLAngle.Mean.SE$SE)

ggplot(data=Limb.MLAngle.Mean.SE, aes(x=Limb.MLAngle.Mean.SE$Stance, y=Limb.MLAngle.Mean.SE$Mean, fill=Type, linetype=Type))+
  scale_y_continuous("Mediolateral Angle (degrees)")+
  scale_x_continuous("Stance")+
  geom_line(size=2, alpha=0.75)+
  geom_ribbon(Limb.MLAngle.MaxMin, alpha=0.5)+
  #scale_fill_manual(values=ggplotColours(n=3)) # Although there are only 2 groups to compare, n=2 makes the hindlimb green
  scale_fill_manual(name="Appendage:", 
                    labels=c("Terrestrial FL    ", "Terrestrial HL    ", "Aquatic FL  ", "Aquatic HL  "),
                    #values=c("ivory4", "ivory4", "ivory4")+
                    values=c("lightblue3", "firebrick3", "green", "brown"))+
  scale_linetype_manual(name="Appendage:", 
                        labels=c("Terrestrial FL    ", "Terrestrial HL    ", "Aquatic FL  ", "Aquatic HL  "),
                        values=c("dashed", "dotted", "dashed", "dotted"))

#Vertical GRF
Af.InterpV.Pec.Mean.SE <- data.frame(Af.InterpV.Pec.Mean, Af.GRF.InterpV.BW.Combined.Pec.SE, Stance)
Af.InterpV.Pec.Mean.SE$Type <- "Terrestrial FL"
names(Af.InterpV.Pec.Mean.SE) <- c("Mean", "SE", "Stance", "Type")
Af.InterpV.Pel.Mean.SE <- data.frame(Af.InterpV.Pel.Mean, Af.GRF.InterpV.BW.Combined.Pel.SE, Stance)
Af.InterpV.Pel.Mean.SE$Type <- "Terrestrial HL"
names(Af.InterpV.Pel.Mean.SE) <- c("Mean", "SE", "Stance", "Type")
Pw.InterpV.Pec.Mean.SE <- data.frame(Pw.InterpV.Pec.Mean, Pw.GRF.InterpV.BW.Combined.Pec.SE, Stance)
Pw.InterpV.Pec.Mean.SE$Type <- "Aquatic FL"
names(Pw.InterpV.Pec.Mean.SE) <- c("Mean", "SE", "Stance", "Type")
Pw.InterpV.Pel.Mean.SE <- data.frame(Pw.InterpV.Pel.Mean, Pw.GRF.InterpV.BW.Combined.Pel.SE, Stance)
Pw.InterpV.Pel.Mean.SE$Type <- "Aquatic HL"
names(Pw.InterpV.Pel.Mean.SE) <- c("Mean", "SE", "Stance", "Type")


Limb.InterpV.Mean.SE <- rbind(Af.InterpV.Pec.Mean.SE, Af.InterpV.Pel.Mean.SE, Pw.InterpV.Pec.Mean.SE, Pw.InterpV.Pel.Mean.SE)

Limb.InterpV.MaxMin <- aes(ymax=Limb.InterpV.Mean.SE$Mean+Limb.InterpV.Mean.SE$SE, ymin=Limb.InterpV.Mean.SE$Mean-Limb.InterpV.Mean.SE$SE)

ggplot(data=Limb.InterpV.Mean.SE, aes(x=Limb.InterpV.Mean.SE$Stance, y=Limb.InterpV.Mean.SE$Mean, fill=Type, linetype=Type))+
  scale_y_continuous("Vertical (BW)")+
  scale_x_continuous("Stance")+
  geom_line(size=2, alpha=0.75)+
  geom_ribbon(Limb.InterpV.MaxMin, alpha=0.5)+
  #scale_fill_manual(values=ggplotColours(n=3)) # Although there are only 2 groups to compare, n=2 makes the hindlimb green
  scale_fill_manual(name="Appendage:", 
                    labels=c("Terrestrial FL    ", "Terrestrial HL    ", "Aquatic FL  ", "Aquatic HL  " ),
                    #values=c("ivory4", "ivory4", "ivory4")+
                    values=c("lightblue3", "firebrick3", "green", "brown"))+
  scale_linetype_manual(name="Appendage:", 
                        labels=c("Terrestrial FL    ", "Terrestrial HL    ", "Aquatic FL  ", "Aquatic HL  "),
                        values=c("dashed", "dotted", "dashed", "dotted"))



#Mediolateral GRF
Af.InterpML.Pec.Mean.SE <- data.frame(Af.InterpML.Pec.Mean, Af.GRF.InterpML.BW.Combined.Pec.SE, Stance)
Af.InterpML.Pec.Mean.SE$Type <- "Terrestrial FL"
names(Af.InterpML.Pec.Mean.SE) <- c("Mean", "SE", "Stance", "Type")
Af.InterpML.Pel.Mean.SE <- data.frame(Af.InterpML.Pel.Mean, Af.GRF.InterpML.BW.Combined.Pel.SE, Stance)
Af.InterpML.Pel.Mean.SE$Type <- "Terrestrial HL"
names(Af.InterpML.Pel.Mean.SE) <- c("Mean", "SE", "Stance", "Type")
Pw.InterpML.Pec.Mean.SE <- data.frame(Pw.InterpML.Pec.Mean, Pw.GRF.InterpML.BW.Combined.Pec.SE, Stance)
Pw.InterpML.Pec.Mean.SE$Type <- "Aquatic FL"
names(Pw.InterpML.Pec.Mean.SE) <- c("Mean", "SE", "Stance", "Type")
Pw.InterpML.Pel.Mean.SE <- data.frame(Pw.InterpML.Pel.Mean, Pw.GRF.InterpML.BW.Combined.Pel.SE, Stance)
Pw.InterpML.Pel.Mean.SE$Type <- "Aquatic HL"
names(Pw.InterpML.Pel.Mean.SE) <- c("Mean", "SE", "Stance", "Type")


Limb.InterpML.Mean.SE <- rbind(Af.InterpML.Pec.Mean.SE, Af.InterpML.Pel.Mean.SE, Pw.InterpML.Pec.Mean.SE, Pw.InterpML.Pel.Mean.SE)

Limb.InterpML.MaxMin <- aes(ymax=Limb.InterpML.Mean.SE$Mean+Limb.InterpML.Mean.SE$SE, ymin=Limb.InterpML.Mean.SE$Mean-Limb.InterpML.Mean.SE$SE)

ggplot(data=Limb.InterpML.Mean.SE, aes(x=Limb.InterpML.Mean.SE$Stance, y=Limb.InterpML.Mean.SE$Mean, fill=Type, linetype=Type))+
  scale_y_continuous("Mediolateral (BW)")+
  scale_x_continuous("Stance")+
  geom_line(size=2, alpha=0.75)+
  geom_ribbon(Limb.InterpML.MaxMin, alpha=0.5)+
  #scale_fill_manual(values=ggplotColours(n=3)) # Although there are only 2 groups to compare, n=2 makes the hindlimb green
  scale_fill_manual(name="Appendage:", 
                    labels=c("Terrestrial FL    ", "Terrestrial HL    ", "Aquatic FL  ", "Aquatic HL  "),
                    #values=c("ivory4", "ivory4", "ivory4")+
                    values=c("lightblue3", "firebrick3", "green", "brown"))+
  scale_linetype_manual(name="Appendage:", 
                        labels=c("Terrestrial FL    ", "Terrestrial HL   ", "Aquatic FL  ", "Aquatic HL  "),
                        values=c("dashed", "dotted", "dashed", "dotted"))


#Anteroposterior GRF
Af.InterpHz.Pec.Mean.SE <- data.frame(Af.InterpHz.Pec.Mean, Af.GRF.InterpHz.BW.Combined.Pec.SE, Stance)
Af.InterpHz.Pec.Mean.SE$Type <- "Terrestrial FL"
names(Af.InterpHz.Pec.Mean.SE) <- c("Mean", "SE", "Stance", "Type")
Af.InterpHz.Pel.Mean.SE <- data.frame(Af.InterpHz.Pel.Mean, Af.GRF.InterpHz.BW.Combined.Pel.SE, Stance)
Af.InterpHz.Pel.Mean.SE$Type <- "Terrestrial HL"
names(Af.InterpHz.Pel.Mean.SE) <- c("Mean", "SE", "Stance", "Type")
Pw.InterpHz.Pec.Mean.SE <- data.frame(Pw.InterpHz.Pec.Mean, Pw.GRF.InterpHz.BW.Combined.Pec.SE, Stance)
Pw.InterpHz.Pec.Mean.SE$Type <- "Aquatic FL"
names(Pw.InterpHz.Pec.Mean.SE) <- c("Mean", "SE", "Stance", "Type")
Pw.InterpHz.Pel.Mean.SE <- data.frame(Pw.InterpHz.Pel.Mean, Pw.GRF.InterpHz.BW.Combined.Pel.SE, Stance)
Pw.InterpHz.Pel.Mean.SE$Type <- "Aquatic HL"
names(Pw.InterpHz.Pel.Mean.SE) <- c("Mean", "SE", "Stance", "Type")


Limb.InterpHz.Mean.SE <- rbind(Af.InterpHz.Pec.Mean.SE, Af.InterpHz.Pel.Mean.SE, Pw.InterpHz.Pec.Mean.SE, Pw.InterpHz.Pel.Mean.SE)

Limb.InterpHz.MaxMin <- aes(ymax=Limb.InterpHz.Mean.SE$Mean+Limb.InterpHz.Mean.SE$SE, ymin=Limb.InterpHz.Mean.SE$Mean-Limb.InterpHz.Mean.SE$SE)

ggplot(data=Limb.InterpHz.Mean.SE, aes(x=Limb.InterpHz.Mean.SE$Stance, y=Limb.InterpHz.Mean.SE$Mean, fill=Type, linetype=Type))+
  scale_y_continuous("Anteroposterior (BW)")+
  scale_x_continuous("Stance")+
  geom_line(size=2, alpha=0.75)+
  geom_ribbon(Limb.InterpHz.MaxMin, alpha=0.5)+
  #scale_fill_manual(values=ggplotColours(n=3)) # Although there are only 2 groups to compare, n=2 makes the hindlimb green
  scale_fill_manual(name="Appendage:", 
                    labels=c("Terrestrial FL    ", "Terrestrial HL    ", "Aquatic FL  ", "Aquatic HL  "),
                    #values=c("ivory4", "ivory4", "ivory4")+
                    values=c("lightblue3", "firebrick3", "green", "brown"))+
  scale_linetype_manual(name="Appendage:", 
                        labels=c("Terrestrial FL    ", "Terrestrial HL    ", "Aquatic FL  ", "Aquatic HL  "),
                        values=c("dashed", "dotted", "dashed", "dotted"))

#Net GRF
Af.NetGRF.Pec.Mean.SE <- data.frame(Af.NetGRF.Pec.Mean, Af.GRF.NetGRF.BW.Combined.Pec.SE, Stance)
Af.NetGRF.Pec.Mean.SE$Type <- "Terrestrial FL"
names(Af.NetGRF.Pec.Mean.SE) <- c("Mean", "SE", "Stance", "Type")
Af.NetGRF.Pel.Mean.SE <- data.frame(Af.NetGRF.Pel.Mean, Af.GRF.NetGRF.BW.Combined.Pel.SE, Stance)
Af.NetGRF.Pel.Mean.SE$Type <- "Terrestrial HL"
names(Af.NetGRF.Pel.Mean.SE) <- c("Mean", "SE", "Stance", "Type")
Pw.NetGRF.Pec.Mean.SE <- data.frame(Pw.NetGRF.Pec.Mean, Pw.GRF.NetGRF.BW.Combined.Pec.SE, Stance)
Pw.NetGRF.Pec.Mean.SE$Type <- "Aquatic FL"
names(Pw.NetGRF.Pec.Mean.SE) <- c("Mean", "SE", "Stance", "Type")
Pw.NetGRF.Pel.Mean.SE <- data.frame(Pw.NetGRF.Pel.Mean, Pw.GRF.NetGRF.BW.Combined.Pel.SE, Stance)
Pw.NetGRF.Pel.Mean.SE$Type <- "Aquatic HL"
names(Pw.NetGRF.Pel.Mean.SE) <- c("Mean", "SE", "Stance", "Type")


Limb.NetGRF.Mean.SE <- rbind(Af.NetGRF.Pec.Mean.SE, Af.NetGRF.Pel.Mean.SE, Pw.NetGRF.Pec.Mean.SE, Pw.NetGRF.Pel.Mean.SE)

Limb.NetGRF.MaxMin <- aes(ymax=Limb.NetGRF.Mean.SE$Mean+Limb.NetGRF.Mean.SE$SE, ymin=Limb.NetGRF.Mean.SE$Mean-Limb.NetGRF.Mean.SE$SE)


ggplot(data=Limb.NetGRF.Mean.SE, aes(x=Limb.NetGRF.Mean.SE$Stance, y=Limb.NetGRF.Mean.SE$Mean, fill=Type, linetype=Type))+
  scale_y_continuous("Net GRF (BW)")+
  scale_x_continuous("Stance")+
  geom_line(size=2, alpha=0.75)+
  geom_ribbon(Limb.NetGRF.MaxMin, alpha=0.5)+
  #scale_fill_manual(values=ggplotColours(n=3)) # Although there are only 2 groups to compare, n=2 makes the hindlimb green
  scale_fill_manual(name="Appendage:", 
                    labels=c("Terrestrial FL    ", "Terrestrial HL    ", "Aquatic FL  ", "Aquatic HL  "),
                    #values=c("ivory4", "ivory4", "ivory4")+
                    values=c("lightblue3", "firebrick3", "green", "brown"))+
  scale_linetype_manual(name="Appendage:", 
                        labels=c("Terrestrial FL    ", "Terrestrial HL    ", "Aquatic FL  ", "Aquatic HL  "),
                        values=c("dashed", "dotted", "dashed", "dotted"))
dev.off()

######## PECTORAL APPENDAGE COMPARISON

## AP Angle 
Af.APAngle.Pec.Mean.SE <- data.frame(Af.APAngle.Pec.Mean, Af.GRF.APAngle.Convert.Combined.Pec.SE, Stance)
Af.APAngle.Pec.Mean.SE$Type <- "Terrestrial FL"
names(Af.APAngle.Pec.Mean.SE) <- c("Mean", "SE", "Stance", "Type")
Pb.APAngle.Pec.Mean.SE <- data.frame(Pb.APAngle.Pec.Mean, Pb.GRF.APAngle.Convert.Combined.Pec.SE, Stance)
Pb.APAngle.Pec.Mean.SE$Type <- "Terrestrial PF"
names(Pb.APAngle.Pec.Mean.SE) <- c("Mean", "SE", "Stance", "Type")
Pw.APAngle.Pec.Mean.SE <- data.frame(Pw.APAngle.Pec.Mean, Pw.GRF.APAngle.Convert.Combined.Pec.SE, Stance)
Pw.APAngle.Pec.Mean.SE$Type <- "Aquatic FL"
names(Pw.APAngle.Pec.Mean.SE) <- c("Mean", "SE", "Stance", "Type")


Pec.APAngle.Mean.SE <- rbind(Af.APAngle.Pec.Mean.SE, Pb.APAngle.Pec.Mean.SE, Pw.APAngle.Pec.Mean.SE)

Pec.APAngle.MaxMin <- aes(ymax=Pec.APAngle.Mean.SE$Mean+Pec.APAngle.Mean.SE$SE, ymin=Pec.APAngle.Mean.SE$Mean-Pec.APAngle.Mean.SE$SE)

ggplot(data=Pec.APAngle.Mean.SE, aes(x=Pec.APAngle.Mean.SE$Stance, y=Pec.APAngle.Mean.SE$Mean, fill=Type, linetype=Type))+
  scale_y_continuous("Anteroposterior Angle (degrees)")+
  scale_x_continuous("Stance")+
  geom_line(size=2, alpha=0.75)+
  geom_ribbon(Pec.APAngle.MaxMin, alpha=0.5)+
  #scale_fill_manual(values=ggplotColours(n=2)) 
  scale_fill_manual(name="Appendage:", 
                    labels=c("Terrestrial FL    ", "Terrestrial PF    ", "Aquatic FL  "),
                    #values=c("ivory4", "ivory4", "ivory4"))+
                    values=c("lightblue3", "darkorange1", "green"))+
  scale_linetype_manual(name="Appendage:", 
                        labels=c("Terrestrial FL    ", "Terrestrial PF    ", "Aquatic FL  "),
                        values=c("dashed", "solid", "dashed"))
  

# ML Angle
Af.MLAngle.Pec.Mean.SE <- data.frame(Af.MLAngle.Pec.Mean, Af.GRF.MLAngle.Convert.Combined.Pec.SE, Stance)
Af.MLAngle.Pec.Mean.SE$Type <- "Terrestrial FL"
names(Af.MLAngle.Pec.Mean.SE) <- c("Mean", "SE", "Stance", "Type")
Pb.MLAngle.Pec.Mean.SE <- data.frame(Pb.MLAngle.Pec.Mean, Pb.GRF.MLAngle.Convert.Combined.Pec.SE, Stance)
Pb.MLAngle.Pec.Mean.SE$Type <- "Terrestrial PF"
names(Pb.MLAngle.Pec.Mean.SE) <- c("Mean", "SE", "Stance", "Type")
Pw.MLAngle.Pec.Mean.SE <- data.frame(Pw.MLAngle.Pec.Mean, Pw.GRF.MLAngle.Convert.Combined.Pec.SE, Stance)
Pw.MLAngle.Pec.Mean.SE$Type <- "Aquatic FL"
names(Pw.MLAngle.Pec.Mean.SE) <- c("Mean", "SE", "Stance", "Type")


Pec.MLAngle.Mean.SE <- rbind(Af.MLAngle.Pec.Mean.SE, Pb.MLAngle.Pec.Mean.SE, Pw.MLAngle.Pec.Mean.SE)

Pec.MLAngle.MaxMin <- aes(ymax=Pec.MLAngle.Mean.SE$Mean+Pec.MLAngle.Mean.SE$SE, ymin=Pec.MLAngle.Mean.SE$Mean-Pec.MLAngle.Mean.SE$SE)

ggplot(data=Pec.MLAngle.Mean.SE, aes(x=Pec.MLAngle.Mean.SE$Stance, y=Pec.MLAngle.Mean.SE$Mean, fill=Type, linetype=Type))+
  scale_y_continuous("Mediolateral Angle (degrees)")+
  scale_x_continuous("Stance")+
  geom_line(size=2, alpha=0.75)+
  geom_ribbon(Pec.MLAngle.MaxMin, alpha=0.5)+
  #scale_fill_manual(values=ggplotColours(n=2)) 
  scale_fill_manual(name="Appendage:", 
                    labels=c("Terrestrial FL    ", "Terrestrial PF    ", "Aquatic FL  "),
                    #values=c("ivory4", "ivory4", "ivory4"))+
                    values=c("lightblue3", "darkorange1", "green"))+
  scale_linetype_manual(name="Appendage:", 
                        labels=c("Terrestrial FL    ", "Terrestrial PF    ", "Aquatic FL  "),
                        values=c("dashed", "solid", "dashed"))
  

#Vertical GRF
Af.InterpV.Pec.Mean.SE <- data.frame(Af.InterpV.Pec.Mean, Af.GRF.InterpV.BW.Combined.Pec.SE, Stance)
Af.InterpV.Pec.Mean.SE$Type <- "Terrestrial FL"
names(Af.InterpV.Pec.Mean.SE) <- c("Mean", "SE", "Stance", "Type")
Pb.InterpV.Pec.Mean.SE <- data.frame(Pb.InterpV.Pec.Mean, Pb.GRF.InterpV.BW.Combined.Pec.SE, Stance)
Pb.InterpV.Pec.Mean.SE$Type <- "Terrestrial PF"
names(Pb.InterpV.Pec.Mean.SE) <- c("Mean", "SE", "Stance", "Type")
Pw.InterpV.Pec.Mean.SE <- data.frame(Pw.InterpV.Pec.Mean, Pw.GRF.InterpV.BW.Combined.Pec.SE, Stance)
Pw.InterpV.Pec.Mean.SE$Type <- "Aquatic FL"
names(Pw.InterpV.Pec.Mean.SE) <- c("Mean", "SE", "Stance", "Type")


Pec.InterpV.Mean.SE <- rbind(Af.InterpV.Pec.Mean.SE, Pb.InterpV.Pec.Mean.SE, Pw.InterpV.Pec.Mean.SE)

Pec.InterpV.MaxMin <- aes(ymax=Pec.InterpV.Mean.SE$Mean+Pec.InterpV.Mean.SE$SE, ymin=Pec.InterpV.Mean.SE$Mean-Pec.InterpV.Mean.SE$SE)

ggplot(data=Pec.InterpV.Mean.SE, aes(x=Pec.InterpV.Mean.SE$Stance, y=Pec.InterpV.Mean.SE$Mean, fill=Type, linetype=Type))+
  scale_y_continuous("Vertical (BW)")+
  scale_x_continuous("Stance")+
  geom_line(size=2, alpha=0.75)+
  geom_ribbon(Pec.InterpV.MaxMin, alpha=0.5)+
  #scale_fill_manual(values=ggplotColours(n=2))
  scale_fill_manual(name="Appendage:", 
                    labels=c("Terrestrial FL    ", "Terrestrial PF    ", "Aquatic FL  "),
                    #values=c("ivory4", "ivory4", "ivory4"))+
                    values=c("lightblue3", "darkorange1", "green"))+
  scale_linetype_manual(name="Appendage:", 
                        labels=c("Terrestrial FL    ", "Terrestrial PF    ", "Aquatic FL"),
                        values=c("dashed", "solid", "dashed"))
  

#Mediolateral GRF
Af.InterpML.Pec.Mean.SE <- data.frame(Af.InterpML.Pec.Mean, Af.GRF.InterpML.BW.Combined.Pec.SE, Stance)
Af.InterpML.Pec.Mean.SE$Type <- "Terrestrial FL"
names(Af.InterpML.Pec.Mean.SE) <- c("Mean", "SE", "Stance", "Type")
Pb.InterpML.Pec.Mean.SE <- data.frame(Pb.InterpML.Pec.Mean, Pb.GRF.InterpML.BW.Combined.Pec.SE, Stance)
Pb.InterpML.Pec.Mean.SE$Type <- "Terrestrial PF"
names(Pb.InterpML.Pec.Mean.SE) <- c("Mean", "SE", "Stance", "Type")
Pw.InterpML.Pec.Mean.SE <- data.frame(Pw.InterpML.Pec.Mean, Pw.GRF.InterpML.BW.Combined.Pec.SE, Stance)
Pw.InterpML.Pec.Mean.SE$Type <- "Aquatic FL"
names(Pw.InterpML.Pec.Mean.SE) <- c("Mean", "SE", "Stance", "Type")


Pec.InterpML.Mean.SE <- rbind(Af.InterpML.Pec.Mean.SE, Pb.InterpML.Pec.Mean.SE, Pw.InterpML.Pec.Mean.SE)

Pec.InterpML.MaxMin <- aes(ymax=Pec.InterpML.Mean.SE$Mean+Pec.InterpML.Mean.SE$SE, ymin=Pec.InterpML.Mean.SE$Mean-Pec.InterpML.Mean.SE$SE)

ggplot(data=Pec.InterpML.Mean.SE, aes(x=Pec.InterpML.Mean.SE$Stance, y=Pec.InterpML.Mean.SE$Mean, fill=Type, linetype=Type))+
  scale_y_continuous("Mediolateral (BW)")+
  scale_x_continuous("Stance")+
  geom_line(size=2, alpha=0.75)+
  geom_ribbon(Pec.InterpML.MaxMin, alpha=0.5)+
  #scale_fill_manual(values=ggplotColours(n=2))
  scale_fill_manual(name="Appendage:", 
                    labels=c("Terrestrial FL    ", "Terrestrial PF    ", "Aquatic FL  "),
                    #values=c("ivory4", "ivory4", "ivory4"))+
                    values=c("lightblue3", "darkorange1", "green"))+
  scale_linetype_manual(name="Appendage:", 
                        labels=c("Terrestrial FL    ", "Terrestrial PF    ", "Aquatic FL  "),
                        values=c("dashed", "solid", "dashed"))
  

#Anteroposterior GRF
Af.InterpHz.Pec.Mean.SE <- data.frame(Af.InterpHz.Pec.Mean, Af.GRF.InterpHz.BW.Combined.Pec.SE, Stance)
Af.InterpHz.Pec.Mean.SE$Type <- "Terrestrial FL"
names(Af.InterpHz.Pec.Mean.SE) <- c("Mean", "SE", "Stance", "Type")
Pb.InterpHz.Pec.Mean.SE <- data.frame(Pb.InterpHz.Pec.Mean, Pb.GRF.InterpHz.BW.Combined.Pec.SE, Stance)
Pb.InterpHz.Pec.Mean.SE$Type <- "Terrestrial PF"
names(Pb.InterpHz.Pec.Mean.SE) <- c("Mean", "SE", "Stance", "Type")
Pw.InterpHz.Pec.Mean.SE <- data.frame(Pw.InterpHz.Pec.Mean, Pw.GRF.InterpHz.BW.Combined.Pec.SE, Stance)
Pw.InterpHz.Pec.Mean.SE$Type <- "Aquatic FL"
names(Pw.InterpHz.Pec.Mean.SE) <- c("Mean", "SE", "Stance", "Type")


Pec.InterpHz.Mean.SE <- rbind(Af.InterpHz.Pec.Mean.SE, Pb.InterpHz.Pec.Mean.SE, Pw.InterpHz.Pec.Mean.SE)

Pec.InterpHz.MaxMin <- aes(ymax=Pec.InterpHz.Mean.SE$Mean+Pec.InterpHz.Mean.SE$SE, ymin=Pec.InterpHz.Mean.SE$Mean-Pec.InterpHz.Mean.SE$SE)

ggplot(data=Pec.InterpHz.Mean.SE, aes(x=Pec.InterpHz.Mean.SE$Stance, y=Pec.InterpHz.Mean.SE$Mean, fill=Type, linetype=Type))+
  scale_y_continuous("Anteroposterior (BW)")+
  scale_x_continuous("Stance")+
  geom_line(size=2, alpha=0.75)+
  geom_ribbon(Pec.InterpHz.MaxMin, alpha=0.5)+
  #scale_fill_manual(values=ggplotColours(n=2))
  scale_fill_manual(name="Appendage:", 
                    labels=c("Terrestrial FL    ", "Terrestrial PF    ", "Aquatic FL  "),
                    #values=c("ivory4", "ivory4", "ivory4"))+
                    values=c("lightblue3", "darkorange1", "green"))+
  scale_linetype_manual(name="Appendage:", 
                        labels=c("Terrestrial FL    ", "Terrestrial PF    ", "Aquatic FL"),
                        values=c("dashed", "solid", "dashed"))
  
  
#Net GRF
Af.NetGRF.Pec.Mean.SE <- data.frame(Af.NetGRF.Pec.Mean, Af.GRF.NetGRF.BW.Combined.Pec.SE, Stance)
Af.NetGRF.Pec.Mean.SE$Type <- "Terrestrial FL"
names(Af.NetGRF.Pec.Mean.SE) <- c("Mean", "SE", "Stance", "Type")
Pb.NetGRF.Pec.Mean.SE <- data.frame(Pb.NetGRF.Pec.Mean, Pb.GRF.NetGRF.BW.Combined.Pec.SE, Stance)
Pb.NetGRF.Pec.Mean.SE$Type <- "Terrestrial PF"
names(Pb.NetGRF.Pec.Mean.SE) <- c("Mean", "SE", "Stance", "Type")
Pw.NetGRF.Pec.Mean.SE <- data.frame(Pw.NetGRF.Pec.Mean, Pw.GRF.NetGRF.BW.Combined.Pec.SE, Stance)
Pw.NetGRF.Pec.Mean.SE$Type <- "Aquatic FL"
names(Pw.NetGRF.Pec.Mean.SE) <- c("Mean", "SE", "Stance", "Type")


Pec.NetGRF.Mean.SE <- rbind(Af.NetGRF.Pec.Mean.SE, Pb.NetGRF.Pec.Mean.SE, Pw.NetGRF.Pec.Mean.SE)

Pec.NetGRF.MaxMin <- aes(ymax=Pec.NetGRF.Mean.SE$Mean+Pec.NetGRF.Mean.SE$SE, ymin=Pec.NetGRF.Mean.SE$Mean-Pec.NetGRF.Mean.SE$SE)

tiff(filename=paste("Pec Net GRF_", SaveDate, ".tif", sep=""), width=9*600, height=5*600, res=600)
ggplot(data=Pec.NetGRF.Mean.SE, aes(x=Pec.NetGRF.Mean.SE$Stance, y=Pec.NetGRF.Mean.SE$Mean, fill=Type, linetype=Type))+
  scale_y_continuous("Net GRF (BW)")+
  scale_x_continuous("Stance")+
  geom_line(size=2, alpha=0.75)+
  geom_ribbon(Pec.NetGRF.MaxMin, alpha=0.5)+
  #scale_fill_manual(values=ggplotColours(n=2))
  scale_fill_manual(name="Appendage:", 
                    labels=c("Terrestrial FL    ", "Terrestrial PF    ", "Aquatic FL  "),
                    #values=c("ivory4", "ivory4", "ivory4"))+
                    values=c("lightblue3", "darkorange1", "green"))+
  scale_linetype_manual(name="Appendage:", 
                        labels=c("Terrestrial FL    ", "Terrestrial PF    ", "Aquatic FL  "),
                        values=c("dashed", "solid", "dashed"))
dev.off()
  

Pec.PdfSave <- paste("Pectoral Force Comparison_", SaveDate, ".pdf", sep="")
pdf(Pec.PdfSave)

## AP Angle 
Af.APAngle.Pec.Mean.SE <- data.frame(Af.APAngle.Pec.Mean, Af.GRF.APAngle.Convert.Combined.Pec.SE, Stance)
Af.APAngle.Pec.Mean.SE$Type <- "Terrestrial FL"
names(Af.APAngle.Pec.Mean.SE) <- c("Mean", "SE", "Stance", "Type")
Pb.APAngle.Pec.Mean.SE <- data.frame(Pb.APAngle.Pec.Mean, Pb.GRF.APAngle.Convert.Combined.Pec.SE, Stance)
Pb.APAngle.Pec.Mean.SE$Type <- "Terrestrial PF"
names(Pb.APAngle.Pec.Mean.SE) <- c("Mean", "SE", "Stance", "Type")
Pw.APAngle.Pec.Mean.SE <- data.frame(Pw.APAngle.Pec.Mean, Pw.GRF.APAngle.Convert.Combined.Pec.SE, Stance)
Pw.APAngle.Pec.Mean.SE$Type <- "Aquatic FL"
names(Pw.APAngle.Pec.Mean.SE) <- c("Mean", "SE", "Stance", "Type")


Pec.APAngle.Mean.SE <- rbind(Af.APAngle.Pec.Mean.SE, Pb.APAngle.Pec.Mean.SE, Pw.APAngle.Pec.Mean.SE)

Pec.APAngle.MaxMin <- aes(ymax=Pec.APAngle.Mean.SE$Mean+Pec.APAngle.Mean.SE$SE, ymin=Pec.APAngle.Mean.SE$Mean-Pec.APAngle.Mean.SE$SE)

ggplot(data=Pec.APAngle.Mean.SE, aes(x=Pec.APAngle.Mean.SE$Stance, y=Pec.APAngle.Mean.SE$Mean, fill=Type, linetype=Type))+
  scale_y_continuous("Anteroposterior Angle (degrees)")+
  scale_x_continuous("Stance")+
  geom_line(size=2, alpha=0.75)+
  geom_ribbon(Pec.APAngle.MaxMin, alpha=0.5)+
  #scale_fill_manual(values=ggplotColours(n=2)) 
  scale_fill_manual(name="Appendage:", 
                    labels=c("Terrestrial FL    ", "Terrestrial PF    ", "Aquatic FL  "),
                    #values=c("ivory4", "ivory4", "ivory4"))+
                    values=c("lightblue3", "darkorange1", "green"))+
  scale_linetype_manual(name="Appendage:", 
                        labels=c("Terrestrial FL    ", "Terrestrial PF    ", "Aquatic FL  "),
                        values=c("dashed", "solid", "dashed"))


# ML Angle
Af.MLAngle.Pec.Mean.SE <- data.frame(Af.MLAngle.Pec.Mean, Af.GRF.MLAngle.Convert.Combined.Pec.SE, Stance)
Af.MLAngle.Pec.Mean.SE$Type <- "Terrestrial FL"
names(Af.MLAngle.Pec.Mean.SE) <- c("Mean", "SE", "Stance", "Type")
Pb.MLAngle.Pec.Mean.SE <- data.frame(Pb.MLAngle.Pec.Mean, Pb.GRF.MLAngle.Convert.Combined.Pec.SE, Stance)
Pb.MLAngle.Pec.Mean.SE$Type <- "Terrestrial PF"
names(Pb.MLAngle.Pec.Mean.SE) <- c("Mean", "SE", "Stance", "Type")
Pw.MLAngle.Pec.Mean.SE <- data.frame(Pw.MLAngle.Pec.Mean, Pw.GRF.MLAngle.Convert.Combined.Pec.SE, Stance)
Pw.MLAngle.Pec.Mean.SE$Type <- "Aquatic FL"
names(Pw.MLAngle.Pec.Mean.SE) <- c("Mean", "SE", "Stance", "Type")


Pec.MLAngle.Mean.SE <- rbind(Af.MLAngle.Pec.Mean.SE, Pb.MLAngle.Pec.Mean.SE, Pw.MLAngle.Pec.Mean.SE)

Pec.MLAngle.MaxMin <- aes(ymax=Pec.MLAngle.Mean.SE$Mean+Pec.MLAngle.Mean.SE$SE, ymin=Pec.MLAngle.Mean.SE$Mean-Pec.MLAngle.Mean.SE$SE)

ggplot(data=Pec.MLAngle.Mean.SE, aes(x=Pec.MLAngle.Mean.SE$Stance, y=Pec.MLAngle.Mean.SE$Mean, fill=Type, linetype=Type))+
  scale_y_continuous("Mediolateral Angle (degrees)")+
  scale_x_continuous("Stance")+
  geom_line(size=2, alpha=0.75)+
  geom_ribbon(Pec.MLAngle.MaxMin, alpha=0.5)+
  #scale_fill_manual(values=ggplotColours(n=2)) 
  scale_fill_manual(name="Appendage:", 
                    labels=c("Terrestrial FL    ", "Terrestrial PF    ", "Aquatic FL  "),
                    #values=c("ivory4", "ivory4", "ivory4"))+
                    values=c("lightblue3", "darkorange1", "green"))+
  scale_linetype_manual(name="Appendage:", 
                        labels=c("Terrestrial FL    ", "Terrestrial PF    ", "Aquatic FL  "),
                        values=c("dashed", "solid", "dashed"))


#Vertical GRF
Af.InterpV.Pec.Mean.SE <- data.frame(Af.InterpV.Pec.Mean, Af.GRF.InterpV.BW.Combined.Pec.SE, Stance)
Af.InterpV.Pec.Mean.SE$Type <- "Terrestrial FL"
names(Af.InterpV.Pec.Mean.SE) <- c("Mean", "SE", "Stance", "Type")
Pb.InterpV.Pec.Mean.SE <- data.frame(Pb.InterpV.Pec.Mean, Pb.GRF.InterpV.BW.Combined.Pec.SE, Stance)
Pb.InterpV.Pec.Mean.SE$Type <- "Terrestrial PF"
names(Pb.InterpV.Pec.Mean.SE) <- c("Mean", "SE", "Stance", "Type")
Pw.InterpV.Pec.Mean.SE <- data.frame(Pw.InterpV.Pec.Mean, Pw.GRF.InterpV.BW.Combined.Pec.SE, Stance)
Pw.InterpV.Pec.Mean.SE$Type <- "Aquatic FL"
names(Pw.InterpV.Pec.Mean.SE) <- c("Mean", "SE", "Stance", "Type")


Pec.InterpV.Mean.SE <- rbind(Af.InterpV.Pec.Mean.SE, Pb.InterpV.Pec.Mean.SE, Pw.InterpV.Pec.Mean.SE)

Pec.InterpV.MaxMin <- aes(ymax=Pec.InterpV.Mean.SE$Mean+Pec.InterpV.Mean.SE$SE, ymin=Pec.InterpV.Mean.SE$Mean-Pec.InterpV.Mean.SE$SE)

ggplot(data=Pec.InterpV.Mean.SE, aes(x=Pec.InterpV.Mean.SE$Stance, y=Pec.InterpV.Mean.SE$Mean, fill=Type, linetype=Type))+
  scale_y_continuous("Vertical (BW)")+
  scale_x_continuous("Stance")+
  geom_line(size=2, alpha=0.75)+
  geom_ribbon(Pec.InterpV.MaxMin, alpha=0.5)+
  #scale_fill_manual(values=ggplotColours(n=2))
  scale_fill_manual(name="Appendage:", 
                    labels=c("Terrestrial FL    ", "Terrestrial PF    ", "Aquatic FL  "),
                    #values=c("ivory4", "ivory4", "ivory4"))+
                    values=c("lightblue3", "darkorange1", "green"))+
  scale_linetype_manual(name="Appendage:", 
                        labels=c("Terrestrial FL    ", "Terrestrial PF    ", "Aquatic FL"),
                        values=c("dashed", "solid", "dashed"))


#Mediolateral GRF
Af.InterpML.Pec.Mean.SE <- data.frame(Af.InterpML.Pec.Mean, Af.GRF.InterpML.BW.Combined.Pec.SE, Stance)
Af.InterpML.Pec.Mean.SE$Type <- "Terrestrial FL"
names(Af.InterpML.Pec.Mean.SE) <- c("Mean", "SE", "Stance", "Type")
Pb.InterpML.Pec.Mean.SE <- data.frame(Pb.InterpML.Pec.Mean, Pb.GRF.InterpML.BW.Combined.Pec.SE, Stance)
Pb.InterpML.Pec.Mean.SE$Type <- "Terrestrial PF"
names(Pb.InterpML.Pec.Mean.SE) <- c("Mean", "SE", "Stance", "Type")
Pw.InterpML.Pec.Mean.SE <- data.frame(Pw.InterpML.Pec.Mean, Pw.GRF.InterpML.BW.Combined.Pec.SE, Stance)
Pw.InterpML.Pec.Mean.SE$Type <- "Aquatic FL"
names(Pw.InterpML.Pec.Mean.SE) <- c("Mean", "SE", "Stance", "Type")


Pec.InterpML.Mean.SE <- rbind(Af.InterpML.Pec.Mean.SE, Pb.InterpML.Pec.Mean.SE, Pw.InterpML.Pec.Mean.SE)

Pec.InterpML.MaxMin <- aes(ymax=Pec.InterpML.Mean.SE$Mean+Pec.InterpML.Mean.SE$SE, ymin=Pec.InterpML.Mean.SE$Mean-Pec.InterpML.Mean.SE$SE)

ggplot(data=Pec.InterpML.Mean.SE, aes(x=Pec.InterpML.Mean.SE$Stance, y=Pec.InterpML.Mean.SE$Mean, fill=Type, linetype=Type))+
  scale_y_continuous("Mediolateral (BW)")+
  scale_x_continuous("Stance")+
  geom_line(size=2, alpha=0.75)+
  geom_ribbon(Pec.InterpML.MaxMin, alpha=0.5)+
  #scale_fill_manual(values=ggplotColours(n=2))
  scale_fill_manual(name="Appendage:", 
                    labels=c("Terrestrial FL    ", "Terrestrial PF    ", "Aquatic FL  "),
                    #values=c("ivory4", "ivory4", "ivory4"))+
                    values=c("lightblue3", "darkorange1", "green"))+
  scale_linetype_manual(name="Appendage:", 
                        labels=c("Terrestrial FL    ", "Terrestrial PF    ", "Aquatic FL  "),
                        values=c("dashed", "solid", "dashed"))


#Anteroposterior GRF
Af.InterpHz.Pec.Mean.SE <- data.frame(Af.InterpHz.Pec.Mean, Af.GRF.InterpHz.BW.Combined.Pec.SE, Stance)
Af.InterpHz.Pec.Mean.SE$Type <- "Terrestrial FL"
names(Af.InterpHz.Pec.Mean.SE) <- c("Mean", "SE", "Stance", "Type")
Pb.InterpHz.Pec.Mean.SE <- data.frame(Pb.InterpHz.Pec.Mean, Pb.GRF.InterpHz.BW.Combined.Pec.SE, Stance)
Pb.InterpHz.Pec.Mean.SE$Type <- "Terrestrial PF"
names(Pb.InterpHz.Pec.Mean.SE) <- c("Mean", "SE", "Stance", "Type")
Pw.InterpHz.Pec.Mean.SE <- data.frame(Pw.InterpHz.Pec.Mean, Pw.GRF.InterpHz.BW.Combined.Pec.SE, Stance)
Pw.InterpHz.Pec.Mean.SE$Type <- "Aquatic FL"
names(Pw.InterpHz.Pec.Mean.SE) <- c("Mean", "SE", "Stance", "Type")


Pec.InterpHz.Mean.SE <- rbind(Af.InterpHz.Pec.Mean.SE, Pb.InterpHz.Pec.Mean.SE, Pw.InterpHz.Pec.Mean.SE)

Pec.InterpHz.MaxMin <- aes(ymax=Pec.InterpHz.Mean.SE$Mean+Pec.InterpHz.Mean.SE$SE, ymin=Pec.InterpHz.Mean.SE$Mean-Pec.InterpHz.Mean.SE$SE)

ggplot(data=Pec.InterpHz.Mean.SE, aes(x=Pec.InterpHz.Mean.SE$Stance, y=Pec.InterpHz.Mean.SE$Mean, fill=Type, linetype=Type))+
  scale_y_continuous("Anteroposterior (BW)")+
  scale_x_continuous("Stance")+
  geom_line(size=2, alpha=0.75)+
  geom_ribbon(Pec.InterpHz.MaxMin, alpha=0.5)+
  #scale_fill_manual(values=ggplotColours(n=2))
  scale_fill_manual(name="Appendage:", 
                    labels=c("Terrestrial FL    ", "Terrestrial PF    ", "Aquatic FL  "),
                    #values=c("ivory4", "ivory4", "ivory4"))+
                    values=c("lightblue3", "darkorange1", "green"))+
  scale_linetype_manual(name="Appendage:", 
                        labels=c("Terrestrial FL    ", "Terrestrial PF    ", "Aquatic FL"),
                        values=c("dashed", "solid", "dashed"))


#Net GRF
Af.NetGRF.Pec.Mean.SE <- data.frame(Af.NetGRF.Pec.Mean, Af.GRF.NetGRF.BW.Combined.Pec.SE, Stance)
Af.NetGRF.Pec.Mean.SE$Type <- "Terrestrial FL"
names(Af.NetGRF.Pec.Mean.SE) <- c("Mean", "SE", "Stance", "Type")
Pb.NetGRF.Pec.Mean.SE <- data.frame(Pb.NetGRF.Pec.Mean, Pb.GRF.NetGRF.BW.Combined.Pec.SE, Stance)
Pb.NetGRF.Pec.Mean.SE$Type <- "Terrestrial PF"
names(Pb.NetGRF.Pec.Mean.SE) <- c("Mean", "SE", "Stance", "Type")
Pw.NetGRF.Pec.Mean.SE <- data.frame(Pw.NetGRF.Pec.Mean, Pw.GRF.NetGRF.BW.Combined.Pec.SE, Stance)
Pw.NetGRF.Pec.Mean.SE$Type <- "Aquatic FL"
names(Pw.NetGRF.Pec.Mean.SE) <- c("Mean", "SE", "Stance", "Type")


Pec.NetGRF.Mean.SE <- rbind(Af.NetGRF.Pec.Mean.SE, Pb.NetGRF.Pec.Mean.SE, Pw.NetGRF.Pec.Mean.SE)

Pec.NetGRF.MaxMin <- aes(ymax=Pec.NetGRF.Mean.SE$Mean+Pec.NetGRF.Mean.SE$SE, ymin=Pec.NetGRF.Mean.SE$Mean-Pec.NetGRF.Mean.SE$SE)

ggplot(data=Pec.NetGRF.Mean.SE, aes(x=Pec.NetGRF.Mean.SE$Stance, y=Pec.NetGRF.Mean.SE$Mean, fill=Type, linetype=Type))+
  scale_y_continuous("Net GRF (BW)")+
  scale_x_continuous("Stance")+
  geom_line(size=2, alpha=0.75)+
  geom_ribbon(Pec.NetGRF.MaxMin, alpha=0.5)+
  #scale_fill_manual(values=ggplotColours(n=2))
  scale_fill_manual(name="Appendage:", 
                    labels=c("Terrestrial FL    ", "Terrestrial PF    ", "Aquatic FL  "),
                    #values=c("ivory4", "ivory4", "ivory4"))+
                    values=c("lightblue3", "darkorange1", "green"))+
  scale_linetype_manual(name="Appendage:", 
                        labels=c("Terrestrial FL    ", "Terrestrial PF    ", "Aquatic FL  "),
                        values=c("dashed", "solid", "dashed"))
dev.off()





setwd('/Users/SandyMKawano/Desktop/Research/Forelimbs on Force Plates/Force Data/R Analysis/Step 4 Results/Angle AP')
Save.Af.APAngle.Pec.Sum <- paste ("Pec_Af_APAngleSum_", SaveDate, ".csv", sep="")
write.table(Af.APAngle.Pec.Summary, file=Save.Af.APAngle.Pec.Sum, sep =",", row.names=FALSE)
Save.Pb.APAngle.Pec.Sum <- paste ("Pec_Pb_APAngleSum_", SaveDate, ".csv", sep="")
write.table(Pb.APAngle.Pec.Summary, file=Save.Pb.APAngle.Pec.Sum, sep =",", row.names=FALSE)
Save.Af.APAngle.Pel.Sum <- paste ("Pel_Af_APAngleSum_", SaveDate, ".csv", sep="")
write.table(Af.APAngle.Pel.Summary, file=Save.Af.APAngle.Pel.Sum, sep =",", row.names=FALSE)
Save.Pw.APAngle.Pec.Sum <- paste ("Pec_Pw_APAngleSum_", SaveDate, ".csv", sep="")
write.table(Pw.APAngle.Pec.Summary, file=Save.Pw.APAngle.Pec.Sum, sep =",", row.names=FALSE)
Save.Pw.APAngle.Pel.Sum <- paste ("Pel_Pw_APAngleSum_", SaveDate, ".csv", sep="")
write.table(Pw.APAngle.Pel.Summary, file=Save.Pw.APAngle.Pel.Sum, sep =",", row.names=FALSE)

setwd('/Users/SandyMKawano/Desktop/Research/Forelimbs on Force Plates/Force Data/R Analysis/Step 4 Results/Angle ML')
Save.Af.MLAngle.Pec.Sum <- paste ("Pec_Af_MLAngleSum_", SaveDate, ".csv", sep="")
write.table(Af.MLAngle.Pec.Summary, file=Save.Af.MLAngle.Pec.Sum, sep =",", row.names=FALSE)
Save.Pb.MLAngle.Pec.Sum <- paste ("Pec_Pb_MLAngleSum_", SaveDate, ".csv", sep="")
write.table(Pb.MLAngle.Pec.Summary, file=Save.Pb.MLAngle.Pec.Sum, sep =",", row.names=FALSE)
Save.Af.MLAngle.Pel.Sum <- paste ("Pel_Af_MLAngleSum_", SaveDate, ".csv", sep="")
write.table(Af.MLAngle.Pel.Summary, file=Save.Af.MLAngle.Pel.Sum, sep =",", row.names=FALSE)
Save.Pw.MLAngle.Pec.Sum <- paste ("Pec_Pw_MLAngleSum_", SaveDate, ".csv", sep="")
write.table(Pw.MLAngle.Pec.Summary, file=Save.Pw.MLAngle.Pec.Sum, sep =",", row.names=FALSE)
Save.Pw.MLAngle.Pel.Sum <- paste ("Pel_Pw_MLAngleSum_", SaveDate, ".csv", sep="")
write.table(Pw.MLAngle.Pel.Summary, file=Save.Pw.MLAngle.Pel.Sum, sep =",", row.names=FALSE)

setwd('/Users/SandyMKawano/Desktop/Research/Forelimbs on Force Plates/Force Data/R Analysis/Step 4 Results/GRF V')
Save.Af.InterpV.Pec.Sum <- paste ("Pec_Af_InterpVSum_", SaveDate, ".csv", sep="")
write.table(Af.InterpV.Pec.Summary, file=Save.Af.InterpV.Pec.Sum, sep =",", row.names=FALSE)
Save.Pb.InterpV.Pec.Sum <- paste ("Pec_Pb_InterpVSum_", SaveDate, ".csv", sep="")
write.table(Pb.InterpV.Pec.Summary, file=Save.Pb.InterpV.Pec.Sum, sep =",", row.names=FALSE)
Save.Af.InterpV.Pel.Sum <- paste ("Pel_Af_InterpVSum_", SaveDate, ".csv", sep="")
write.table(Af.InterpV.Pel.Summary, file=Save.Af.InterpV.Pel.Sum, sep =",", row.names=FALSE)
Save.Pw.InterpV.Pec.Sum <- paste ("Pec_Pw_InterpVSum_", SaveDate, ".csv", sep="")
write.table(Pw.InterpV.Pec.Summary, file=Save.Pw.InterpV.Pec.Sum, sep =",", row.names=FALSE)
Save.Pw.InterpV.Pel.Sum <- paste ("Pel_Pw_InterpVSum_", SaveDate, ".csv", sep="")
write.table(Pw.InterpV.Pel.Summary, file=Save.Pw.InterpV.Pel.Sum, sep =",", row.names=FALSE)

setwd('/Users/SandyMKawano/Desktop/Research/Forelimbs on Force Plates/Force Data/R Analysis/Step 4 Results/GRF ML')
Save.Af.InterpML.Pec.Sum <- paste ("Pec_Af_InterpMLSum_", SaveDate, ".csv", sep="")
write.table(Af.InterpML.Pec.Summary, file=Save.Af.InterpML.Pec.Sum, sep =",", row.names=FALSE)
Save.Pb.InterpML.Pec.Sum <- paste ("Pec_Pb_InterpMLSum_", SaveDate, ".csv", sep="")
write.table(Pb.InterpML.Pec.Summary, file=Save.Pb.InterpML.Pec.Sum, sep =",", row.names=FALSE)
Save.Af.InterpML.Pel.Sum <- paste ("Pel_Af_InterpMLSum_", SaveDate, ".csv", sep="")
write.table(Af.InterpML.Pel.Summary, file=Save.Af.InterpML.Pel.Sum, sep =",", row.names=FALSE)
Save.Pw.InterpML.Pec.Sum <- paste ("Pec_Pw_InterpMLSum_", SaveDate, ".csv", sep="")
write.table(Pw.InterpML.Pec.Summary, file=Save.Pw.InterpML.Pec.Sum, sep =",", row.names=FALSE)
Save.Pw.InterpML.Pel.Sum <- paste ("Pel_Pw_InterpMLSum_", SaveDate, ".csv", sep="")
write.table(Pw.InterpML.Pel.Summary, file=Save.Pw.InterpML.Pel.Sum, sep =",", row.names=FALSE)

setwd('/Users/SandyMKawano/Desktop/Research/Forelimbs on Force Plates/Force Data/R Analysis/Step 4 Results/GRF Hz')
Save.Af.InterpHz.Pec.Sum <- paste ("Pec_Af_InterpHzSum_", SaveDate, ".csv", sep="")
write.table(Af.InterpHz.Pec.Summary, file=Save.Af.InterpHz.Pec.Sum, sep =",", row.names=FALSE)
Save.Pb.InterpHz.Pec.Sum <- paste ("Pec_Pb_InterpHzSum_", SaveDate, ".csv", sep="")
write.table(Pb.InterpHz.Pec.Summary, file=Save.Pb.InterpHz.Pec.Sum, sep =",", row.names=FALSE)
Save.Af.InterpHz.Pel.Sum <- paste ("Pel_Af_InterpHzSum_", SaveDate, ".csv", sep="")
write.table(Af.InterpHz.Pel.Summary, file=Save.Af.InterpHz.Pel.Sum, sep =",", row.names=FALSE)
Save.Pw.InterpHz.Pec.Sum <- paste ("Pec_Pw_InterpHzSum_", SaveDate, ".csv", sep="")
write.table(Pw.InterpHz.Pec.Summary, file=Save.Pw.InterpHz.Pec.Sum, sep =",", row.names=FALSE)
Save.Pw.InterpHz.Pel.Sum <- paste ("Pel_Pw_InterpHzSum_", SaveDate, ".csv", sep="")
write.table(Pw.InterpHz.Pel.Summary, file=Save.Pw.InterpHz.Pel.Sum, sep =",", row.names=FALSE)

setwd('/Users/SandyMKawano/Desktop/Research/Forelimbs on Force Plates/Force Data/R Analysis/Step 4 Results/Net GRF')
Save.Af.NetGRF.Pec.Sum <- paste ("Pec_Af_NetGRFSum_", SaveDate, ".csv", sep="")
write.table(Af.NetGRF.Pec.Summary, file=Save.Af.NetGRF.Pec.Sum, sep =",", row.names=FALSE)
Save.Pb.NetGRF.Pec.Sum <- paste ("Pec_Pb_NetGRFSum_", SaveDate, ".csv", sep="")
write.table(Pb.NetGRF.Pec.Summary, file=Save.Pb.NetGRF.Pec.Sum, sep =",", row.names=FALSE)
Save.Af.NetGRF.Pel.Sum <- paste ("Pel_Af_NetGRFSum_", SaveDate, ".csv", sep="")
write.table(Af.NetGRF.Pel.Summary, file=Save.Af.NetGRF.Pel.Sum, sep =",", row.names=FALSE)
Save.Pw.NetGRF.Pec.Sum <- paste ("Pec_Pw_NetGRFSum_", SaveDate, ".csv", sep="")
write.table(Pw.NetGRF.Pec.Summary, file=Save.Pw.NetGRF.Pec.Sum, sep =",", row.names=FALSE)
Save.Pw.NetGRF.Pel.Sum <- paste ("Pel_Pw_NetGRFSum_", SaveDate, ".csv", sep="")
write.table(Pw.NetGRF.Pel.Summary, file=Save.Pw.NetGRF.Pel.Sum, sep =",", row.names=FALSE)



################  UNUSED CODE ##################

# a <- ggplot(data=okay, aes(x=Stance, y=Mean, ymin=Mean-SE, ymax=Mean+SE))
# b <- a+geom_line()
# c <- b+geom_ribbon(alpha=0.5) #alpha changes the transparency of the lines
# e <- b+geom_smooth(se=F)+stat_summary(fun.data=mean_se, geom="ribbon", alpha=0.25)
# 
# Af.APAngle.Pec.MaxMin <- aes(ymax=Af.APAngle.Pec.Mean.SE$Mean+Af.APAngle.Pec.Mean.SE$SE, ymin=Af.APAngle.Pec.Mean.SE$Mean-Af.APAngle.Pec.Mean.SE$SE)
# 
# ggplot(data=Af.APAngle.Pec.Mean.SE, aes(x=Af.APAngle.Pec.Mean.SE$Stance, y=Af.APAngle.Pec.Mean.SE$Mean))+
#   scale_y_continuous("Anteroposterior Angle (degrees)")+
#   scale_x_continuous("Stance")+
#   geom_line(colour="green", size=2, alpha=0.75)+
#   geom_ribbon(Af.APAngle.Pec.MaxMin, alpha=0.5, colour="green", fill="green")

# 



## Pectoral files
#GRF.APAngle.Convert.Pec <- rep(NA, length(GRF.Pec))   # Creating empty matrix to input data from loop
#for (i in 1:length(GRF.Pec)) {GRF.APAngle.Convert.Pec[i] <- GRF.Pec[[i]][2]}
#GRF.APAngle.Convert.Combined.Pec <- as.data.frame(do.call("rbind", GRF.APAngle.Convert.Pec))
#names(GRF.APAngle.Convert.Combined.Pec) <- PercentStance
#GRF.APAngle.Convert.Combined.Pec$File.name <- names(GRF.Pec)
#
#GRF.MLAngle.Convert.Pec <- rep(NA, length(GRF.Pec))   # Creating empty matrix to input data from loop
#for (i in 1:length(GRF.Pec)) {GRF.MLAngle.Convert.Pec[i] <- GRF.Pec[[i]][3]}
#GRF.MLAngle.Convert.Combined.Pec <- as.data.frame(do.call("rbind", GRF.MLAngle.Convert.Pec))
#names(GRF.MLAngle.Convert.Combined.Pec) <- PercentStance
#GRF.MLAngle.Convert.Combined.Pec$File.name <- names(GRF.Pec)
#
#GRF.InterpV.BW.Pec <- rep(NA, length(GRF.Pec))   # Creating empty matrix to input data from loop
#for (i in 1:length(GRF.Pec)) {GRF.InterpV.BW.Pec[i] <- GRF.Pec[[i]][4]}
#GRF.InterpV.BW.Combined.Pec <- as.data.frame(do.call("rbind", GRF.InterpV.BW.Pec))
#names(GRF.InterpV.BW.Combined.Pec) <- PercentStance
#GRF.InterpV.BW.Combined.Pec$File.name <- names(GRF.Pec)
#
#GRF.InterpML.BW.Pec <- rep(NA, length(GRF.Pec))   # Creating empty matrix to input data from loop
#for (i in 1:length(GRF.Pec)) {GRF.InterpML.BW.Pec[i] <- GRF.Pec[[i]][5]}
#GRF.InterpML.BW.Combined.Pec <- as.data.frame(do.call("rbind", GRF.InterpML.BW.Pec))
#names(GRF.InterpML.BW.Combined.Pec) <- PercentStance
#GRF.InterpML.BW.Combined.Pec$File.name <- names(GRF.Pec)
#
#GRF.InterpHz.BW.Pec <- rep(NA, length(GRF.Pec))   # Creating empty matrix to input data from loop
#for (i in 1:length(GRF.Pec)) {GRF.InterpHz.BW.Pec[i] <- GRF.Pec[[i]][6]}
#GRF.InterpHz.BW.Combined.Pec <- as.data.frame(do.call("rbind", GRF.InterpHz.BW.Pec))
#names(GRF.InterpHz.BW.Combined.Pec) <- PercentStance
#GRF.InterpHz.BW.Combined.Pec$File.name <- names(GRF.Pec)
#
#
## Pelvic files
#GRF.APAngle.Convert.Pel <- rep(NA, length(GRF.Pel))   # Creating empty matrix to input data from loop
#for (i in 1:length(GRF.Pel)) {GRF.APAngle.Convert.Pel[i] <- GRF.Pel[[i]][2]}
#GRF.APAngle.Convert.Combined.Pel <- as.data.frame(do.call("rbind", GRF.APAngle.Convert.Pel))
#names(GRF.APAngle.Convert.Combined.Pel) <- PercentStance
#GRF.APAngle.Convert.Combined.Pel$File.name <- names(GRF.Pel)
#
#GRF.MLAngle.Convert.Pel <- rep(NA, length(GRF.Pel))   # Creating empty matrix to input data from loop
#for (i in 1:length(GRF.Pel)) {GRF.MLAngle.Convert.Pel[i] <- GRF.Pel[[i]][3]}
#GRF.MLAngle.Convert.Combined.Pel <- as.data.frame(do.call("rbind", GRF.MLAngle.Convert.Pel))
#names(GRF.MLAngle.Convert.Combined.Pel) <- PercentStance
#GRF.MLAngle.Convert.Combined.Pel$File.name <- names(GRF.Pel)
#
#GRF.InterpV.BW.Pel <- rep(NA, length(GRF.Pel))   # Creating empty matrix to input data from loop
#for (i in 1:length(GRF.Pel)) {GRF.InterpV.BW.Pel[i] <- GRF.Pel[[i]][4]}
#GRF.InterpV.BW.Combined.Pel <- as.data.frame(do.call("rbind", GRF.InterpV.BW.Pel))
#names(GRF.InterpV.BW.Combined.Pel) <- PercentStance
#GRF.InterpV.BW.Combined.Pel$File.name <- names(GRF.Pel)
#
#GRF.InterpML.BW.Pel <- rep(NA, length(GRF.Pel))   # Creating empty matrix to input data from loop
#for (i in 1:length(GRF.Pel)) {GRF.InterpML.BW.Pel[i] <- GRF.Pel[[i]][5]}
#GRF.InterpML.BW.Combined.Pel <- as.data.frame(do.call("rbind", GRF.InterpML.BW.Pel))
#names(GRF.InterpML.BW.Combined.Pel) <- PercentStance
#GRF.InterpML.BW.Combined.Pel$File.name <- names(GRF.Pel)
#
#GRF.InterpHz.BW.Pel <- rep(NA, length(GRF.Pel))   # Creating empty matrix to input data from loop
#for (i in 1:length(GRF.Pel)) {GRF.InterpHz.BW.Pel[i] <- GRF.Pel[[i]][6]}
#GRF.InterpHz.BW.Combined.Pel <- as.data.frame(do.call("rbind", GRF.InterpHz.BW.Pel))
#names(GRF.InterpHz.BW.Combined.Pel) <- PercentStance
#GRF.InterpHz.BW.Combined.Pel$File.name <- names(GRF.Pel)

#
## Dividing Ambystoma files into pectoral and pelvic files
#Af.GRF.APAngle.Convert.Combined.Pec <- Af.GRF.APAngle.Convert.Combined[which(substring(Af.GRF.APAngle.Convert.Combined$File.name, 9,11)=="Pec"),]
#Af.GRF.MLAngle.Convert.Combined.Pec <- Af.GRF.MLAngle.Convert.Combined[which(substring(Af.GRF.MLAngle.Convert.Combined$File.name, 9,11)=="Pec"),]
#Af.GRF.InterpV.BW.Combined.Pec <- Af.GRF.InterpV.BW.Combined[which(substring(Af.GRF.InterpV.BW.Combined$File.name, 9,11)=="Pec"),]
#Af.GRF.InterpML.BW.Combined.Pec <- Af.GRF.InterpML.BW.Combined[which(substring(Af.GRF.InterpML.BW.Combined$File.name, 9,11)=="Pec"),]
#Af.GRF.InterpHz.BW.Combined.Pec <- Af.GRF.InterpHz.BW.Combined[which(substring(Af.GRF.InterpHz.BW.Combined$File.name, 9,11)=="Pec"),]
#
#Af.GRF.APAngle.Convert.Combined.Pel <- Af.GRF.APAngle.Convert.Combined[which(substring(Af.GRF.APAngle.Convert.Combined$File.name, 9,11)=="Pel"),]
#Af.GRF.MLAngle.Convert.Combined.Pel <- Af.GRF.MLAngle.Convert.Combined[which(substring(Af.GRF.MLAngle.Convert.Combined$File.name, 9,11)=="Pel"),]
#Af.GRF.InterpV.BW.Combined.Pel <- Af.GRF.InterpV.BW.Combined[which(substring(Af.GRF.InterpV.BW.Combined$File.name, 9,11)=="Pel"),]
#Af.GRF.InterpML.BW.Combined.Pel <- Af.GRF.InterpML.BW.Combined[which(substring(Af.GRF.InterpML.BW.Combined$File.name, 9,11)=="Pel"),]
#Af.GRF.InterpHz.BW.Combined.Pel <- Af.GRF.InterpHz.BW.Combined[which(substring(Af.GRF.InterpHz.BW.Combined$File.name, 9,11)=="Pel"),]