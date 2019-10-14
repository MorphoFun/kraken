################# Low pass zero phase Butterworth Filter ######################
##  Used to digitally filter GRF output for force platform analyses          ##
##                                                                           ##
##  Written by Sandy Kawano (sandy.kawano@gmail.com)                         ##
##                                                                           ##
###############################################################################

## Notes about pectoral appendage force plate analyses
# Filter specifications for digital low pass Butterworth filter
# Wp = passband corner frequency; Need to normalize to Nyquist frequency
# (divide by 1/2 your sampling frequency); this avoids aliasing
# Ws = stopband coner frequency; normalized by Nyquist frequency as well
# Rp = passband ripple in dB; max permissible passband loss
# Rs = stopband attenuation in dB; # dB the stopband is down from the passband

# For my salamander forceplate data
# Sampling rate (Fs) = 5000 Hz; 0.5Fs = 2500 Hz
# Limb Cycle Frequency ranges from about 0.7-2.78 Hz
# Natural, unloaded frequencies of the force platform were 190 Hz in all
# three directions
# Thus, stopband should be 190 Hz
# Passband should thus be from 0 to at least 5 Hz to make sure that it
# includes any limb cycles
# Stopband should start at 190 Hz to get rid of noise from force plate
# and then go to the Nyquist frequency for my data (5000/2 = 2500 Hz)
# Want less than 3dB of ripple in the passband (pretty standard)
# For now, setting at least 40dB of attenuation in the stopband for 99%
# attenuation

# Clear everything in the R workspace, so nothing gets mixed up between different trials
rm(list=ls(all=TRUE))

today <- Sys.Date()
SaveDate <- format(today, format="%y%m%d")

## Loading the library for the signal processing functions
## Need to install signal package if using for the first time
library(signal)


## Loading in file with filter specifications
# Read the most recent Filter Spec file
setwd('/Users/SandyMKawano/Desktop/Research/Forelimbs on Force Plates/Force Data/R Analysis/Step 2 Filter Data/Butterworth Filter Specs')
FilterFilePath <- setwd('/Users/SandyMKawano/Desktop/Research/Forelimbs on Force Plates/Force Data/R Analysis/Step 2 Filter Data/Butterworth Filter Specs')
FilterFile.Ind <- list.files(FilterFilePath, pattern=".csv", full=TRUE)
FilterFile.Dates <-substring(FilterFile.Ind, 144, 149)
FilterFile.Newest <- max(as.numeric(FilterFile.Dates))
FilterFile.Use <- FilterFile.Ind[substring(FilterFile.Ind, 144, 149) %in% FilterFile.Newest]
FilterFile <- data.frame(read.csv(FilterFile.Use, header=TRUE))

#R Assigning the filter specification variables
Freq <- FilterFile$Fs # Frequency of data
FreqN <- Freq/2 # Frequency normalized to Nyquist frequency
PassbandFreqN <- FilterFile$PbF/FreqN  # Passband frequency AKA Wp; normalized to Nyquist frequency
StopbandFreqN <- FilterFile$SbF/FreqN  # Stopband frequency AKA Ws; normalized to Nyquist frequency
PassbandRip <- FilterFile$Rp # Passband Ripple (dB) AKA Rp
StopbandAtt <- FilterFile$Rs # Stopband Attenuation (dB) AKA Rs

# Determing the order and cut-off frequency based on the filter specifications
ButtOrderCut <- buttord(PassbandFreqN, StopbandFreqN, PassbandRip, StopbandAtt)

# Creating low pass Butterworth filter of order n
ButtFiltLP <- butter(ButtOrderCut)

# Normalized cut-off frequency
CutFreqN <- ButtOrderCut$Wc/FreqN

## Reading in the file to be filtered
# Set the directory, to facilitate choosing file
setwd('/Users/SandyMKawano/Desktop/Research/Forelimbs on Force Plates/Force Data/R Analysis/Step 1 Calibrate and Organize Data/Prep Data for Filtering')

# Select the interpolated data files
myFile <- file.choose()
myData <- read.csv(myFile, header=TRUE)
Trial <- substring(myFile, 145, 151)

# Determining length of the jump in units of sweeps
# Remember that extra sweeps were added to the front and back of the trace to reduce edge effects (see Step 1)
# These extra sweeps were 25% the length of the jump.
# Thus, length(Limb) = length(FilterPrepData)/1.5
NumStanceSweeps <- round(nrow(myData)/1.5)

# Determining how many sweeps were added to the Limb Cycle for padding, and for that extra 25% of the trace to reduce edge effects
# The 100 is there because you pad the front and end of the trace with 100 sweeps of the padding value (see below for the calc's)
ExcessToRemove <- 100+c(round(0.25*NumStanceSweeps))

## Padding the front and end of data sets to make sure that
# Using the average of the trace from the first and last 30% of the FilterPrep data,
# and then repeating that average 100 times at the front and end of the trace
# Using 30% for taking the average in order to account for some overlap (5%) within the jump
FrontVAvg <- mean(myData$GRF0SumVN[1:c(round(0.3*length(myData$GRF0SumVN)))])
FrontVPad <- rep(FrontVAvg, 100)
EndVAvg <- mean(myData$GRF0SumVN[(length(myData$GRF0SumVN)-c(round(0.3*length(myData$GRF0SumVN))-1)):length(myData$GRF0SumVN)])
EndVPad <- rep(EndVAvg, 100)
# Adding the front and end padding to the vertical GRF data
SumVNPadded <- c(FrontVPad, myData$GRF0SumVN, EndVPad)
# Using forward and reverse filtering to prevent phase shifts (i.e., making this a zero phase filter)
FilterVPad <- filtfilt(ButtFiltLP$b, ButtFiltLP$a, SumVNPadded)
# Now removing the padding and 25% extra data on the front and back of trace, so that you're only saving the data within your jump
FilterVN <- FilterVPad[-c(1:ExcessToRemove, c(length(FilterVPad)-ExcessToRemove+1):length(FilterVPad))]


# Mediolateral
# See notes for vertical GRF for what the code is doing
FrontMLAvg <- mean(myData$GRF0SumMLN[1:c(round(0.3*length(myData$GRF0SumMLN)))])
FrontMLPad <- rep(FrontMLAvg, 100)
EndMLAvg <- mean(myData$GRF0SumMLN[(length(myData$GRF0SumMLN)-c(round(0.3*length(myData$GRF0SumMLN))-1)):length(myData$GRF0SumMLN)])
EndMLPad <- rep(EndMLAvg, 100)
SumMLNPadded <- c(FrontMLPad, myData$GRF0SumMLN, EndMLPad)
FilterMLPad <- filtfilt(ButtFiltLP$b, ButtFiltLP$a, SumMLNPadded)
FilterMLN <- FilterMLPad[-c(1:ExcessToRemove, c(length(FilterMLPad)-ExcessToRemove+1):length(FilterMLPad))]

# Horizontal
# See notes for the vertical GRF for what the code is doing
FrontHzAvg <- mean(myData$GRF0SumHzN[1:c(round(0.3*length(myData$GRF0SumHzN)))])
FrontHzPad <- rep(FrontHzAvg, 100)
EndHzAvg <- mean(myData$GRF0SumHzN[(length(myData$GRF0SumHzN)-c(round(0.3*length(myData$GRF0SumHzN))-1)):length(myData$GRF0SumHzN)])
EndHzPad <- rep(EndHzAvg, 100)
SumHzNPadded <- c(FrontHzPad, myData$GRF0SumHzN, EndHzPad)
FilterHzPad <- filtfilt(ButtFiltLP$b, ButtFiltLP$a, SumHzNPadded)
FilterHzN <- FilterHzPad[-c(1:ExcessToRemove, c(length(FilterHzPad)-ExcessToRemove+1):length(FilterHzPad))]

Appendage <- substring(myFile, 163, 165)

# Plotting the data
# Opening up a new window that you can send figure/graphs to
quartz(width=11)
# par() allows you to customize your window (in this case, saying you want 1 row of 3 graphs arranged in columns)
par(mfrow=c(1,3), oma = c(3, 0, 2, 0))  # oma = outer margin with 2 lines above the top of the graphs
      ## Vertical component of GRF graph
      # Raw
      plot(1:length(SumVNPadded), SumVNPadded, xlab='Sweep', ylab='GRF - Vertical (N)', main='Zeroed GRF - Vertical', type="l", col="blue")
      # Drawing the Filtered data as a line over the Raw data plot (type='l' is for drawing a line, col="black" draws that line in black)
      lines(1:length(FilterVPad), FilterVPad, type='l' , col="black", lwd=2)
      abline(v=ExcessToRemove, col='magenta', lty=2, lwd=2)
      abline(v=c(length(FilterVPad)-ExcessToRemove), col='magenta', lty=2, lwd=2)

      ## Mediolateral component of GRF graph
      # Raw
      plot(1:length(SumMLNPadded), SumMLNPadded, xlab='Sweep', ylab='GRF - Mediolateral (N)', main='Zeroed GRF - Mediolateral', type="l", col="red")
      # Drawing the Filtered data as a line over the Raw data plot (type='l' is for drawing a line, col="black" draws that line in black)
      lines(1:length(FilterMLPad), FilterMLPad, type='l' , col="black", lwd=2)
      abline(v=ExcessToRemove, col='magenta', lty=2, lwd=2)
      abline(v=c(length(FilterMLPad)-ExcessToRemove), col='magenta', lty=2, lwd=2)

      ## Horizontal (Anteroposterior) component of GRF graph
      # Raw
      plot(1:length(SumHzNPadded), SumHzNPadded, xlab='Sweep', ylab='GRF - Horizontal (N)', main='Zeroed GRF - Horizontal', type="l", col="forestgreen")
      # Drawing the Filtered data as a line over the Raw data plot (type='l' is for drawing a line, col="black" draws that line in black)
      lines(1:length(FilterHzPad), FilterHzPad, type='l' , col="black", lwd=2)
      abline(v=ExcessToRemove, col='magenta', lty=2, lwd=2)
      abline(v=c(length(FilterHzPad)-ExcessToRemove), col='magenta', lty=2, lwd=2)

      # Creating a new variable Graph Title from the Trial name
      GraphTitle <- paste(Trial, "_", Appendage, sep="")
      # writes an overall title over the graphs
      mtext(GraphTitle, line=0.5, outer=TRUE)
      mtext('Area within dashed pink lines = Stance', side=1, outer=TRUE, col = 'magenta')
      mtext('Colored traces = Raw Data; Black trace = Filtered Data', side=1, line=1.5, outer=TRUE, col = "slategrey")

# Saving the graph
# Creating the file path for the pdf file that you'd like to save now

PdfSave <- paste(substring(myFile, 1, 84), "/Step 2 Filter Data/Graphs/", Trial, "_", Appendage, "_Filter_", SaveDate, ".pdf", sep="")
pdf(PdfSave, width=11)
# Now preparing to save the graphs as a pdf
# Everything after the pdf() will be what you want to save as a pdf, and is the same code that generated the graphs above
# par() allows you to customize your window (in this case, saying you want 1 row of 3 graphs arranged in columns)
par(mfrow=c(1,3), oma = c(3, 0, 2, 0))  # oma = outer margin with 2 lines above the top of the graphs
      ## Vertical component of GRF graph
      # Raw
      plot(1:length(SumVNPadded), SumVNPadded, xlab='Sweep', ylab='GRF - Vertical (N)', main='Zeroed GRF - Vertical', type="l", col="blue")
      # Drawing the Filtered data as a line over the Raw data plot (type='l' is for drawing a line, col="black" draws that line in black)
      lines(1:length(FilterVPad), FilterVPad, type='l' , col="black", lwd=2)
      abline(v=ExcessToRemove, col='magenta', lty=2, lwd=2)
      abline(v=c(length(FilterVPad)-ExcessToRemove), col='magenta', lty=2, lwd=2)

      ## Mediolateral component of GRF graph
      # Raw
      plot(1:length(SumMLNPadded), SumMLNPadded, xlab='Sweep', ylab='GRF - Mediolateral (N)', main='Zeroed GRF - Mediolateral', type="l", col="red")
      # Drawing the Filtered data as a line over the Raw data plot (type='l' is for drawing a line, col="black" draws that line in black)
      lines(1:length(FilterMLPad), FilterMLPad, type='l' , col="black", lwd=2)
      abline(v=ExcessToRemove, col='magenta', lty=2, lwd=2)
      abline(v=c(length(FilterMLPad)-ExcessToRemove), col='magenta', lty=2, lwd=2)

      ## Horizontal (Anteroposterior) component of GRF graph
      # Raw
      plot(1:length(SumHzNPadded), SumHzNPadded, xlab='Sweep', ylab='GRF - Horizontal (N)', main='Zeroed GRF - Horizontal', type="l", col="forestgreen")
      # Drawing the Filtered data as a line over the Raw data plot (type='l' is for drawing a line, col="black" draws that line in black)
      lines(1:length(FilterHzPad), FilterHzPad, type='l' , col="black", lwd=2)
      abline(v=ExcessToRemove, col='magenta', lty=2, lwd=2)
      abline(v=c(length(FilterHzPad)-ExcessToRemove), col='magenta', lty=2, lwd=2)

      # Creating a new variable Graph Title from the Trial name
      GraphTitle <- paste(Trial, "_", Appendage, sep="")
      # writes an overall title over the graphs
      mtext(GraphTitle, line=0.5, outer=TRUE)
      mtext('Area within dashed pink lines = Stance', side=1, outer=TRUE, col = 'magenta')
      mtext('Colored traces = Raw Data; Black trace = Filtered Data', side=1, line=1.5, outer=TRUE, col = "slategrey")

dev.off()    # tells you that this is the end of everything that you wanted to save as a pdf



## Interpolating the data to 101 points (e.g., 0-100% at 1% increments)

# Remember to use just the Stance, and not the length of myData, which includes excess sweeps

# Figure out the interval that would give you 101 equidistant points within your data set
# Want 101 points because want 0% -> 100% at 1% intervals
InterpN <- 101 # Establishing how many points you want to interpolate the points to
N <- NumStanceSweeps/InterpN  # Determining the increments of data that would create 101 equidistant points

# Remove the excess 25% front the front and back of myData
StanceData <- myData[-c(1:c(round(0.25*NumStanceSweeps)), (length(myData$GRF0SumVN)-(c(round(0.25*NumStanceSweeps))-1)):length(myData$GRF0SumVN)),]

# Set your new interval that you would like your data interpolated to using your
# new sample size (e.g., 101)
X <- seq(StanceData$Sweep[1], StanceData$Sweep[length(StanceData$Sweep)], N)

# Interpolate dataset to a new sample size (InterpN)
# Using spline b/c I have a polynomial function and spline has a better capability of capturing it
# Some have said that spline is also less prone to error as well
# while still maintaining efficiency (i.e., does not take forever to compute; at least in MATLAB)
InterpV.N <- interp1(StanceData$Sweep, FilterVN, X, 'spline')
InterpML.N <- interp1(StanceData$Sweep, FilterMLN, X, 'spline')
InterpHz.N <- interp1(StanceData$Sweep, FilterHzN, X, 'spline')

# Organizing data
# Creating the increments of the stance in 1% increments
PercentStance <- seq(0,100,1)
# Creating a dataframe of all the Filtered data within stance
FilteredAll <- data.frame(StanceData$Sweep, FilterVN, FilterMLN, FilterHzN)
# Creating a dataframe of the Filtered data interpolated to 101 points for your jump
FilteredInterp <- data.frame(PercentStance, InterpV.N, InterpML.N, InterpHz.N)
# Creating a dataframe that includes both the raw and filtered data
RawAndFilterAll <- data.frame(StanceData, FilterVN, FilterMLN, FilterHzN)


# Saving the filtered data
# First, code the name of the file
#Appendage <- substring(myFile, 163, 165)
# All Data
if (Appendage == 'Pec') SaveFileNameAll.Pec <- paste(Trial,"_FilteredPecAll_", SaveDate, ".csv", sep="")
if (Appendage == 'Pel') SaveFileNameAll.Pel <- paste(Trial,"_FilteredPelAll_", SaveDate, ".csv", sep="")
# Interpolated data
if (Appendage == 'Pec') SaveFileNameInterp.Pec <- paste(Trial,"_FilteredPecInterp_", SaveDate, ".csv", sep="")
if (Appendage == 'Pel') SaveFileNameInterp.Pel <- paste(Trial,"_FilteredPelInterp_", SaveDate, ".csv", sep="")


# Changing directory
setwd('/Users/SandyMKawano/Desktop/Research/Forelimbs on Force Plates/Force Data/R Analysis/Step 2 Filter Data/Filtered All')
if (Appendage == 'Pec') write.table(RawAndFilterAll, file=SaveFileNameAll.Pec, sep =",", row.names=FALSE)
if (Appendage == 'Pel') write.table(RawAndFilterAll, file=SaveFileNameAll.Pel, sep =",", row.names=FALSE)

setwd('/Users/SandyMKawano/Desktop/Research/Forelimbs on Force Plates/Force Data/R Analysis/Step 2 Filter Data/Filtered Interpolated')
if (Appendage == 'Pec') write.table(FilteredInterp, file=SaveFileNameInterp.Pec, sep =",", row.names=FALSE)
if (Appendage == 'Pel') write.table(FilteredInterp, file=SaveFileNameInterp.Pel, sep =",", row.names=FALSE)
