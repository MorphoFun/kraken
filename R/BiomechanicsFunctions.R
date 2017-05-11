############ Biomechanics Functions ########

#### PCSA ####
#' @title Physiological cross-sectional area of muscles (PCSA)
#'
#' @name pcsa
#'
#' @description \code{pcsa} allows one to estimate the physiological cross-sectional area of a muscle
#'
#' @usage pcsa(mass, pennationAngle, fiberLength, density, stringAsFactors = FALSE, ...)
#'
#' @param \code{mass} Numerical value of the muscle mass, in units of kilograms.
#' @param \code{pennationAngle} Numerical value of the pennational angle of the muscle, in units of degrees.
#' @param \code{fasicleLength} Numerical value of the muscle fiber length, in units of meters. Often assumed to be average muscle fiber length.
#' @param \code{density} Numerical value of the muscle density. Defaults to 1060 kg/m^3, a typical value for striated muscles (Biewener 2003)
#'
#' @details See description in Hutchinson et al. (2015) and Sacks and Roy (1982) for more details regarding the calculation of PCSA.
#' @references Biewener AA. 2003. Animal locomotion. Oxford, UK: Oxford University Press.
#' @references Hutchinson JR, Rankin JW, Rubenson J, Rosenbluth KH, Siston RA, Delp SL. 2015. Musculoskeletal modelling of an ostrich (Struthio camelus) pelvic limb: influence of limb orientation on muscular capacity during locomotion. \url{http://dx.doi.org/10.7717/peerj.1001}
#' @references Sacks RD, Roy RR. 1982. Architecture of The Hind Limb Muscles of Cats: Functional Significance. Journal of Morphology, 185-195. \url{http://onlinelibrary.wiley.com/doi/10.1002/jmor.1051730206/abstract}
#'
#' @examples
#'
#' pcsa(0.3788, 0, 0.174)
#'
#' @export

pcsa <- function(mass, pennationAngle, fascicleLength, density = 1060, ...) {
  d <- (mass * cos(pennationAngle))/(fascicleLength*density)
  return(d)
}


##### PROFILE PLOTS ####
#' @title Generating profile plots for longitudinal data, with repeated measures
#'
#' @name profilePlotR
#'
#' @description Generating profile plot of multiple trials
#'
#' @usage profilePlotR(d = d, xname = xname, yname = yname, groupname = groupname, subgroupname = subgroupname, rowname = rowname, colors = c("red", "blue"), title = "plot", xlab = "x", ylab = "y", highlight = NULL, ...)
#'
#' @param \code{d} data (currently only accepts input for one variable at a time).
#' @param \code{xname} x-axis variable name.
#' @param \code{yname} y-axis variable name.
#' @param \code{groupname} variable name for the overall group that is being evaluated (e.g., species).
#' @param \code{subgroupname} variable name for a subgroup of the overall group (e.g., individual within species)
#' @param \code{rowname} variable name for the rows (e.g., ID number).
#' @param \code{title} character string for the title of the plot.
#' @param \code{xlab} character string for the x-axis label.
#' @param \code{ylab} character string for the y-axis label.
#' @param \code{colorlinesby} variable name for the grouping by which to color the individual lines by.
#' @param \code{highlight} Optional feature to highlight certain data points.
#'
#' @details Function to quickly generate profile plots for data. For instance, kinematic plots over time for multiple individuals that have multiple trials of data collected.

#'
#' @examples
#'
#' profilePlotR(subset(AT_Kine2, Variable == "AbductAdductAngle"), "PercentStance", "value", groupname = "Appendage", subgroupname = "Ind", rowname = "Filename", highlight = AT_Kine2_AAA_subset, title = "Abduction versus Adduction", xlab = "PercentStance", ylab = "Degrees")
#'
#' @import ggplot2
#' @import vegan
#' @export

profilePlotR <- function(d = d, xname = xname, yname = yname, groupname = groupname, subgroupname = subgroupname, rowname = rowname, title = "plot", xlab = "x", ylab = "y", colorlinesby = subgroupname, highlight = NULL, ...) {
  interact <- c(groupname, subgroupname, rowname)
  ggplot(d, aes_string(x = xname, y = yname)) +
    geom_line(aes_string(group = paste0('interaction(', paste0(interact, collapse = ', ' ),')'), color = colorlinesby, linetype = groupname), alpha = 0.3) +
    geom_smooth(aes_string(fill = groupname, linetype = groupname, color = groupname), color = "black",  alpha = 0.6) + # include means for each ind with 95% CI shading
    theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank()) + # get rid of gridlines
    theme(panel.background=element_blank()) + # make background white
    theme(axis.line.x =element_line(colour="black", linetype="solid"),
      axis.line.y =element_line(colour="black", linetype="solid")) + # put black lines for axes
    ggtitle(title) + theme(plot.title=element_text(hjust=0.5, size=15, face="bold")) +
    labs(x = xlab, y = ylab) +
    if(is.null(highlight) == FALSE) {
      geom_line(data = highlight, aes_string(group = paste0('interaction(', paste0(interact, collapse = ', ' ),')'), color = subgroupname, linetype = groupname),  size = 1.5, alpha = 0.5)
    }
}


##### IMPULSE ####
#' @title Calculate impulse from ground reaction force data
#'
#' @name impulse
#'
#' @description Estimates impulse from data of force over time by calculating the area under the curve.
#'
#' @usage impulse(time, GRF)
#'
#' @param \code{time} a vector of numerical data on the time sequence
#' @param \code{GRF} an array of columns for the force data (assumed that force data are already synchronized to the time data)
#'
#' @details Impulse is a measure of the force applied over a specific time period. The time and force data should already be ordered so that the first row is the beginning of the trial and the last row is the end of the trial.
#'
#' @examples
#' time <- seq(1:10)
#' set.seed(123)
#' GRF <- data.frame(x = rnorm(10), y = rnorm(10), z = rnorm(10))
#'
#' impulse(time, GRF)
#'
#'
#' @import zoo
#' @export

impulse <- function(time, GRF) {
  if (ncol(GRF) == 1) {
    totalImpulse <- sum(diff(time)*rollmean(GRF,2))
    rollImpulse <- diff(time)*rollmean(GRF,2)
  } else {
    totalImpulse <- sapply(GRF, FUN = function(x) sum(diff(time)*rollmean(x,2)))
    rollImpulse <- sapply(GRF, FUN = function(x) diff(time)*rollmean(x,2))
  }
  output <- list(
    totalImpulse = totalImpulse,
    rollImpulse = rollImpulse
  )
  return(output)
}



#### voltToForce ####
#' @title Convert voltage data to ground reaction force data
#'
#' @name voltToForce
#'
#' @description Takes voltage data from a multi-axis force platform and then converts those data to ground reation forces.
#'
#' @usage voltToForce(calib, lightStartFrame, startFrame, endFrame, videoHz, forceHz, zeroStart)
#'
#' @param \code{df} A data.frame of data recorded from 12 channels of a force platform in the following order: 1 column trigger, 4 columns for the verticals, 1 column for the sum of the verticals, 2 columns for the mediolaterals, 1 column for the sum of the mediolaterals, 2 anteroposteriors, and 1 column for the sum of the anteroposteriors.
#' @param \code{calib} A vector of numeric data containing the calibrations in the vertical, mediolateral, and anteroposterior directions, respectively.
#' @param \code{lightStartFrame} A numeric / integer value depicting the frame in which the light is triggered on the high-speed video.
#' @param \code{startFrame} A numeric / integer value depicting the frame number in which the behavior started (e.g., first contact of the foot onto the force plate)
#' @param \code{endFrame} A numeric / integer value depicting the frame number in which the behavior ended (e.g., penultimate frame to the foot lifting off from the force plate.
#' @param \code{zeroStart} A numeric / integer value depicting the sweep number in which to start zeroing the data (i.e., estimate the baseline). Should be at least 1000 sweeps away from any activity. 
#' @param \code{videoHz} A numeric / integer value depicting the frame rate for the high-speed videos. 100 Hz set as a default.
#' @param \code{forceHz} A numeric / integer value depicting the recording rate for the force plate data. 5000 Hz set as a default.
#' @param \code{filename} Option to identify the name of the file for indexing later.
#' @param \code{saveData} Option to save all of the data from this function as .csv files. Default is set to "no".
#' @details These procedures follow the methodology used in Kawano and Blob (2013) and Kawano et al. 2016. It is assumed that the output from the force platform contain 12 channels in the following order: trigger, four verticals, sum of the verticals, two mediolateral, sum of the mediolaterals, two anteroposterior, and the sum of the anteroposteriors.
#' @references Kawano SM, Blob RW. 2013. Propulsive forces of mudskipper fins and salamander limbs during terrestrial locomotion: implications for the invasion of land. Integrative and Comparative Biology 53(2): 283-294. \url{https://academic.oup.com/icb/article/53/2/283/806410/Propulsive-Forces-of-Mudskipper-Fins-and}
#' @references Kawano SM, Economy DR, Kennedy MS, Dean D, Blob RW. 2016. Comparative limb bone loading in the humerus and femur of the tiger salamander Ambystoma tigrinum: testing the "mixed-chain" hypothesis for skeletal safety factors. Journal of Experimental Biology 219: 341-353. \url{http://jeb.biologists.org/content/219/3/341}
#'
#' @examples
#' 
#' voltToForce(df, calib = c(-0.710, 1.337, 1.563), lightStartFrame = 248, startFrame = 20, endFrame = 196, zeroStart = 22000)
#'
#' @export
#' 

#### TO DO FOR voltToForce: add option to insert pectoral vs. pelvic vs. tail start and end points on the plot

voltToForce <- function(df, calib, lightStartFrame, startFrame, endFrame, zeroStart, videoHz = 100, forceHz= 5000, filename = NULL, saveData = "no", ...) {
  myData <- df[,c(1:12)] # ensuring only the desired columns are used for analysis
  names(myData) <- c("light_Volts", "vert1_Volts", "vert2_Volts", "vert3_Volts", "vert4_Volts", "vertSum_Volts", "ML1_Volts", "ML2_Volts", "MLSum_Volts", "AP1_Volts", "AP2_Volts", "APSum_Volts")
  myData$sweep <- 1:nrow(myData)
  
  # Calibrating data
  myData$vertSumCalib_N <- myData$vertSum_Volts*as.numeric(calib[1])
  myData$MLSumCalib_N <- myData$MLSum_Volts*as.numeric(calib[2])
  myData$APSumCalib_N <- myData$APSum_Volts*as.numeric(calib[3])
  
  # Putting forces in terms of GRF (which is opposite in direction to the force produced by the limb onto the force plate)
  myData$GRF_vertSumCalib_N <- myData$vertSumCalib_N # Already made negative based on the calibration calculations conducted earlier in the excel calibration files
  myData$GRF_MLSumCalib_N <- -myData$MLSumCalib_N
  myData$GRF_APSumCalib_N <- -myData$APSumCalib_N
  
  # Determining what sweep number the light is turned on, so we can sync with video frames
  lightSwitch <- myData[which(myData$light_Volts<0),]
  lightOnset <- lightSwitch[which(lightSwitch$sweep == min(lightSwitch$sweep)),]
  
  # Synching video frames with force sweep numbers
  startSweep <- lightOnset$sweep-((lightStartFrame-startFrame)*(forceHz/videoHz))
  endSweep <- lightOnset$sweep-((lightStartFrame-endFrame)*(forceHz/videoHz))
  
  # Correcting for the offset from the baseline (zero)
  offsetCalcStart <- zeroStart-2000
  offsetCalcEnd <- zeroStart-1000
  GRF_vertSumCalib_N_Offset <- mean(myData$GRF_vertSumCalib_N[offsetCalcStart:offsetCalcEnd])
  GRF_MLSumCalib_N_Offset <- mean(myData$GRF_MLSumCalib_N[offsetCalcStart:offsetCalcEnd])
  GRF_APSumCalib_N_Offset <- mean(myData$GRF_APSumCalib_N[offsetCalcStart:offsetCalcEnd])
  
  # Zeroing the force trace data using the offset value
  myData$GRF_vertSumCalib_N_Zero <- myData$GRF_vertSumCalib_N - GRF_vertSumCalib_N_Offset
  myData$GRF_MLSumCalib_N_Zero <- myData$GRF_MLSumCalib_N - GRF_MLSumCalib_N_Offset
  myData$GRF_APSumCalib_N_Zero <- myData$GRF_APSumCalib_N - GRF_APSumCalib_N_Offset
  
  # Plotting the force data that has been converted to Newtons and zero'd
  plotStart <- as.numeric(startSweep-1000)
  plotEnd <- as.numeric(endSweep+1000)
  windows(width=10)
  par(mfrow=c(1,3), oma = c(0, 0, 2, 0))  # oma = outer margin with 2 lines above the top of the graphs
  # Vertical component of GRF graph
  plot(myData$sweep[plotStart:plotEnd], myData$GRF_vertSumCalib_N_Zero[plotStart:plotEnd], xlab='Sweep', ylab='GRF - Vertical (N)', main='Zeroed GRF (Vertical) Force', type="l", col="blue")
  #points(ImpPointsX[1,], ImpPoints_GRFVert[1,], type='p', pch='O', col='cyan')
  #text(ImpPointsX[1,], ImpPoints_GRFVert[1,], labels=names(ImpPointsX), pos=3, font=2) # pos: 1 = below, 2 = left, 3 = above, 4 = right
  
  # Mediolateral component of GRF graph
  plot(myData$sweep[plotStart:plotEnd], myData$GRF_MLSumCalib_N_Zero[plotStart:plotEnd], xlab='Sweep', ylab='GRF - Mediolateral (N)', main='Zeroed GRF (Mediolateral) Force', type="l", col="red")
  #points(ImpPointsX[1,], ImpPoints_GRFML[1,], type='p', pch='O', col='cyan')
  #text(ImpPointsX[1,], ImpPoints_GRFML[1,], labels=names(ImpPointsX), pos=3, font=2) # pos: 1 = below, 2 = left, 3 = above, 4 = right
  
  # Anteroposterior component of GRF graph
  plot(myData$sweep[plotStart:plotEnd], myData$GRF_APSumCalib_N_Zero[plotStart:plotEnd], xlab='Sweep', ylab='GRF - Anteroposterior (N)', main='Zeroed GRF (Anteroposterior) Force', type="l", col="forestgreen")
  #points(ImpPointsX[1,], ImpPoints_GRFHz[1,], type='p', pch='O', col='cyan')
  #text(ImpPointsX[1,], ImpPoints_GRFHz[1,], labels=names(ImpPointsX), pos=3, font=2) # pos: 1 = below, 2 = left, 3 = above, 4 = right
  
  mtext(filename, line=0.5, outer=TRUE)  # writes an overall title over the graphs
  
  # Prepping the data to be filtered
  cycleSweeps <- myData$sweep[startSweep:endSweep]
  cycleGRFVert <- myData$GRF_vertSumCalib_N_Zero[startSweep:endSweep]
  cycleGRFML <- myData$GRF_MLSumCalib_N_Zero[startSweep:endSweep]
  cycleGRFAP <- myData$GRF_APSumCalib_N_Zero[startSweep:endSweep]
  filterPrep <- data.frame(cycleSweeps,cycleGRFVert, cycleGRFML, cycleGRFAP)
  names(filterPrep) <- c('sweep', 'GRF0SumVN', 'GRF0SumMLN', 'GRF0SumAPN')
  
  # Saving all of the data
  if (!saveData == "no") {
    today <- Sys.Date()
    saveDate <- format(today, format="%y%m%d")
    saveAllDataName <- paste(filename,"_allPrep_",saveDate, ".csv", sep="")
    write.table(myData, file=saveAllDataName, sep =",", row.names=FALSE)
    
    # Saving the filter prep data
    saveFileName <- paste(filename,"_filterPrep_", saveDate, ".csv", sep="")
    write.table(filterPrep, file=saveFileName, sep =",", row.names=FALSE)
  }
  
  # output the data
  output <- list(
    filterPrep = filterPrep,
    filename = filename,
    allData = myData,
    startSweep = startSweep,
    endSweep = endSweep
  )

}

#### buttFilteR ####
#' @title Apply a butterworth filter to data
#'
#' @name buttFilteR
#'
#' @description Applies a butterworth filter to data, using information about the data to determine what polynomial to use.
#'
#' @usage 
#'
#' @param \code{df} A list containing the file name and an array of data containing the independent varabile (time) as the first column that is followed by the dependent variables. Can take voltToForce output as an input.
#' @param \code{Fs} A numeric value indicating the sampling frequency. 5000 Hz is set as a default.
#' @param \code{PbF} A numeric value indicating the pass-band frequency. 6 is set as a default.
#' @param \code{SbF} A numeric value indicating the stop-band frequency. 190 is set as a default.
#' @param \code{Rp} A numeric value indicating passband ripple in dB; represents the max permissible passband loss. 2 dB is set as a default.
#' @param \code{Rs} A numeric value indicating stopband attenuation in dB; respresents the dB the stopband is down from the passband. 40 dB is set as a default.

#' @details These procedures follow the methodology used in Kawano and Blob (2013) and Kawano et al. 2016. It is assumed that the output from the force platform contain 12 channels in the following order: trigger, four verticals, sum of the verticals, two mediolateral, sum of the mediolaterals, two anteroposterior, and the sum of the anteroposteriors.
#' @references Kawano SM, Blob RW. 2013. Propulsive forces of mudskipper fins and salamander limbs during terrestrial locomotion: implications for the invasion of land. Integrative and Comparative Biology 53(2): 283-294. \url{https://academic.oup.com/icb/article/53/2/283/806410/Propulsive-Forces-of-Mudskipper-Fins-and}
#' @references Kawano SM, Economy DR, Kennedy MS, Dean D, Blob RW. 2016. Comparative limb bone loading in the humerus and femur of the tiger salamander Ambystoma tigrinum: testing the "mixed-chain" hypothesis for skeletal safety factors. Journal of Experimental Biology 219: 341-353. \url{http://jeb.biologists.org/content/219/3/341}
#'
#' @examples
#' 
#' GRF <- voltToForce(df, calib = c(-0.710, 1.337, 1.563), lightStartFrame = 248, startFrame = 20, endFrame = 196, zeroStart = 22000)
#' GRF_filtered <- buttFilteR(GRF)
#'
#' @import signal
#' @export

buttFilteR <- function(df, Fs = 5000, PbF = 6, SbF = 190, Rp = 2, Rs = 40, ...) {
  # Assigning the filter specification variables
  freq <- Fs # Frequency of data
  freqN <- freq/2 # Frequency normalized to Nyquist frequency
  passbandFreqN <- PbF/freqN  # Passband frequency AKA Wp; normalized to Nyquist frequency
  stopbandFreqN <- SbF/freqN  # Stopband frequency AKA Ws; normalized to Nyquist frequency
  passbandRip <- Rp # Passband Ripple (dB) AKA Rp
  stopbandAtt <- Rs # Stopband Attenuation (dB) AKA Rs
  
  # Determing the order and cut-off frequency based on the filter specifications
  buttOrderCut <- buttord(passbandFreqN, stopbandFreqN, passbandRip, stopbandAtt)
  
  # Creating low pass Butterworth filter of order n
  buttFiltLP <- butter(buttOrderCut)
  
  # Normalized cut-off frequency
  cutFreqN <- buttOrderCut$Wc/freqN

  filename <- df$filename
  
  # Padding the edges to remove edge effects
  # Doing this by mirroring the data on each side
  pad <- data.frame(df[[1]][,1], rev(df[[1]][,2]), rev(df[[1]][,3]), rev(df[[1]][,4]))
  names(pad) <- names(df[[1]]) 
  paddedData <- rbind(pad, df[[1]], pad)
  
  # Using forward and reverse filtering to prevent phase shifts (i.e., making this a zero phase filter)
  filterVPad <- filtfilt(buttFiltLP$b, buttFiltLP$a, paddedData[,2])
  filterMLPad <- filtfilt(buttFiltLP$b, buttFiltLP$a, paddedData[,3])
  filterAPPad <- filtfilt(buttFiltLP$b, buttFiltLP$a, paddedData[,4])
  
  # Now removing the padding
  filterVN <- filterVPad[(nrow(pad)+1):(nrow(pad)*2)]
  filterMLN <- filterMLPad[(nrow(pad)+1):(nrow(pad)*2)]
  filterAPN <- filterAPPad[(nrow(pad)+1):(nrow(pad)*2)]
  
  # Plotting the data
  # Opening up a new window that you can send figure/graphs to
  dev.new(width = 11)

  # par() allows you to customize your window (in this case, saying you want 1 row of 3 graphs arranged in columns)
  par(mfrow=c(1,3), oma = c(3, 0, 2, 0))  # oma = outer margin with 2 lines above the top of the graphs
  ## Vertical component of GRF graph
  # Raw
  plot(1:length(df[[1]][,2]), df[[1]][,2], xlab='Sweep', ylab='GRF - Vertical (N)', main='Zeroed GRF - Vertical', type="l", col="blue")
  # Drawing the Filtered data as a line over the Raw data plot (type='l' is for drawing a line, col="black" draws that line in black)
  lines(1:length(filterVN), filterVN, type='l' , col="black", lwd=2)

  
  ## Mediolateral component of GRF graph
  # Raw
  plot(1:length(df[[1]][,3]), df[[1]][,3], xlab='Sweep', ylab='GRF - Mediolateral (N)', main='Zeroed GRF - Mediolateral', type="l", col="red")
  # Drawing the Filtered data as a line over the Raw data plot (type='l' is for drawing a line, col="black" draws that line in black)
  lines(1:length(filterMLN), filterMLN, type='l' , col="black", lwd=2)
  
  ## Horizontal (Anteroposterior) component of GRF graph
  plot(1:length(df[[1]][,4]), df[[1]][,4], xlab='Sweep', ylab='GRF - Anteroposterior (N)', main='Zeroed GRF - Anteroposterior', type="l", col="forestgreen")
  
  # Drawing the Filtered data as a line over the Raw data plot (type='l' is for drawing a line, col="black" draws that line in black)
  lines(1:length(filterAPN), filterAPN, type='l' , col="black", lwd=2)
  
  # Creating a new variable Graph Title from the Trial name
  GraphTitle <- filename
  # writes an overall title over the graphs
  mtext(GraphTitle, line=0.5, outer=TRUE)
  mtext('Colored traces = Raw Data; Black trace = Filtered Data', side=1, line=1.5, outer=TRUE, col = "slategrey")
  
  ##### Interpotating the data to 101 points
  # Want 101 points because want 0% -> 100% at 1% intervals
  InterpN <- 100 # Establishing how many spaces between the points you want to interpolate to
  N <- (df$endSweep - df$startSweep)/InterpN  # Determining the increments of data that would create 101 equidistant points
  
  # Set your new interval that you would like your data interpolated to using your
  # new sample size (e.g., 101)
  X <- seq(df$filterPrep[1,1], df$filterPrep[nrow(df$filterPrep),1], N)
  
  # Interpolate dataset to a new sample size (InterpN)
  # Using spline b/c I have a polynomial function and spline has a better capability of capturing it
  # Some have said that spline is also less prone to error as well
  # while still maintaining efficiency (i.e., does not take forever to compute; at least in MATLAB)
  InterpV_N <- interp1(df$filterPrep[,1], filterVN, X, 'spline')
  InterpML_N <- interp1(df$filterPrep[,1], filterMLN, X, 'spline')
  InterpAP_N <- interp1(df$filterPrep[,1], filterAPN, X, 'spline')
  
  # Organizing data
  # Creating the increments of the stance in 1% increments
  PercentStance <- seq(0,100,1)
  # Creating a dataframe of all the Filtered data within stance
  FilteredAll <- data.frame(sweep = df$filterPrep[,1], filterVN, filterMLN, filterAPN)
  # Creating a dataframe of the Filtered data interpolated to 101 points for your jump
  FilteredInterp <- data.frame(PercentStance, InterpV_N, InterpML_N, InterpAP_N)
  # Creating a dataframe that includes both the raw and filtered data
  RawAndFilterAll <- data.frame(df$filterPrep, filterVN, filterMLN, filterAPN)
  
  
  output <- list(
    GRF0Sum_filter_sweeps = FilteredAll,
    GRF0Sum_filter_interp = FilteredInterp,
    GRF0Sum_filter_all = RawAndFilterAll
  )
  
  return(output)
}
