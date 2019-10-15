################  GRF of Fish and Salamander Appendages ########################
## Code for evaluating force data collected from the Blob lab force plate #####
# Abbreviations: Vert = Vertical, ML = mediolateral, Hz = Horizontal (refers to anteroposterior direction)

######################  OVERVIEW OF TRIALS TO KEEP #############################
####  Criteria for excluding trials:
### - if peak net GRF occurs at the point of overlap
### - if peak net GRF occurs during after after a significant body bump onto the plate
### - if there is any abnormal behaviors in video (e.g., high stepping, swinging down onto arm, etc.)

# Clear everything in the R workspace, so nothing gets mixed up between different trials
rm(list=ls(all=TRUE))

today <- Sys.Date()
SaveDate <- format(today, format="%y%m%d")

# Read the most recent Video Info file
setwd('E:/Forelimbs on Force Plates_150202/Trials Kept')
TrialsPath <- setwd('E:/Forelimbs on Force Plates_150202/Trials Kept')
Trials.Ind <- list.files(TrialsPath, pattern=".csv", full=TRUE)
Trials.Dates <-substring(Trials.Ind, 61, 66)
Trials.Newest <- max(as.numeric(Trials.Dates))
Trials.Use <- Trials.Ind[substring(Trials.Ind, 61, 66) %in% Trials.Newest]
myData<- data.frame(read.csv(Trials.Use, header=TRUE))

# Identifying individuals to file
myData$Ind <- substring(myData$File.name, 1,4)

# Subsetting pectoral vs. pelvic files
PectoralFiles <- subset(myData, Appendage =='Pectoral')
PelvicFiles <- subset(myData, Appendage == 'Pelvic')

# Producing table that organizes results
Pec.EndResult <- table(PectoralFiles$Status, as.factor(PectoralFiles$Ind))
Pel.EndResult <- table(PelvicFiles$Status, as.factor(PelvicFiles$Ind))

# Outputting table with which trials to use for analysis
Pec.Use <- table(PectoralFiles$Use, as.factor(PectoralFiles$Ind))
Pel.Use <- table(PelvicFiles$Use, as.factor(PelvicFiles$Ind))

# Compiling list of the files that are being included in the analysis
Pec.Save <- PectoralFiles[which(PectoralFiles$Use=="Yes"),]
Pel.Save <- PelvicFiles[which(PelvicFiles$Use=="Yes"),]

# Saving the list of files to analyze
library(xlsx)
Trials.SaveName <- paste("GRF Trials for Analysis_", SaveDate, ".xlsx", sep="")
write.xlsx(data.frame(Pec.Save), file=Trials.SaveName, sheetName="Pectoral", row.names=F)
write.xlsx(data.frame(Pel.Save), file=Trials.SaveName, sheetName="Pelvic", row.names=F, append=T)
