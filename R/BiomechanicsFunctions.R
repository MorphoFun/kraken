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
#' 
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
#' Note: this code was originally written with a force plate that was oriented with the following configuration: Front left = Yellow (Y), Front Right = Pink (P), Back Left = Orange (O), and Back Right = Green (G).
#' Deviations from this configuration may require adjustments in the signage of the GRF data.
#'
#' @usage voltToForce(calib, lightStartFrame, startFrame, endFrame, videoHz, forceHz, zeroStart)
#'
#' @param \code{df} A data.frame of data recorded from 12 channels of a force platform in the following order: 1 column trigger, 4 columns for the verticals, 1 column for the sum of the verticals, 2 columns for the mediolaterals, 1 column for the sum of the mediolaterals, 2 anteroposteriors, and 1 column for the sum of the anteroposteriors.
#' @param \code{calib} A vector of numeric data containing the calibrations in the vertical, mediolateral, and anteroposterior directions, respectively.
#' @param \code{zeroStart} A numeric / integer value depicting the frame number in which the data move away from the baseline. The data will be zero'd about 1000 sweeps prior to this frame. 
#' @param \code{lightStartFrame} A numeric / integer value depicting the frame in which the light is triggered on the high-speed video.
#' @param \code{startFrame} A numeric / integer value depicting the frame number in which the behavior started (e.g., first contact of the foot onto the force plate)
#' @param \code{endFrame} A numeric / integer value depicting the frame number in which the behavior ended (e.g., penultimate frame to the foot lifting off from the force plate.
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
#' voltToForce(df, calib = c(-0.710, 1.337, 1.563), lightStartFrame = 248, startFrame = 20, endFrame = 196, filename = "af01f01")
#'
#' @export
#' 

voltToForce <- function(df, calib, zeroStart, lightStartFrame, startFrame, endFrame, videoHz = 100, forceHz= 5000, filename = NULL, BW = NULL, saveData = "no", ...) {
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
  myData$GRF_APSumCalib_N <- myData$APSumCalib_N
  
  # Determining what sweep number the light is turned on, so we can sync with video frames
  lightSwitch <- myData[which(myData$light_Volts<0),]
  lightOnset <- lightSwitch[which(lightSwitch$sweep == min(lightSwitch$sweep)),]
  
  # Synching video frames with force sweep numbers
  startSweep <- lightOnset$sweep-((lightStartFrame-startFrame)*(forceHz/videoHz))
  endSweep <- lightOnset$sweep-((lightStartFrame-endFrame)*(forceHz/videoHz))
  zeroSweep <- lightOnset$sweep-((lightStartFrame-zeroStart)*(forceHz/videoHz))
  
  # Correcting for the offset from the baseline (zero)
  offsetCalcStart <- zeroSweep-2000
  offsetCalcEnd <- zeroSweep-1000
  GRF_vertSumCalib_N_Offset <- mean(myData$GRF_vertSumCalib_N[offsetCalcStart:offsetCalcEnd])
  GRF_MLSumCalib_N_Offset <- mean(myData$GRF_MLSumCalib_N[offsetCalcStart:offsetCalcEnd])
  GRF_APSumCalib_N_Offset <- mean(myData$GRF_APSumCalib_N[offsetCalcStart:offsetCalcEnd])
  
  # Zeroing the force trace data using the offset value
  myData$GRF_vertSumCalib_N_Zero <- myData$GRF_vertSumCalib_N - GRF_vertSumCalib_N_Offset
  myData$GRF_MLSumCalib_N_Zero <- myData$GRF_MLSumCalib_N - GRF_MLSumCalib_N_Offset
  myData$GRF_APSumCalib_N_Zero <- myData$GRF_APSumCalib_N - GRF_APSumCalib_N_Offset
  
  # Adding a time column
  myData$Time_s <- c(0, cumsum(rep(1/forceHz, (nrow(myData)-1))))
  
  # Plotting the force data that has been converted to Newtons and zero'd
  plotStart <- as.numeric(startSweep)
  plotEnd <- as.numeric(endSweep)
  
  ImpPoints_X <- data.frame(startSweep, endSweep)
  ImpPoints_GRFVert <- data.frame(myData$GRF_vertSumCalib_N_Zero[startSweep], myData$GRF_vertSumCalib_N_Zero[endSweep])
  ImpPoints_GRFML <- data.frame(myData$GRF_MLSumCalib_N_Zero[startSweep], myData$GRF_MLSumCalib_N_Zero[endSweep])
  ImpPoints_GRFAP <- data.frame(myData$GRF_APSumCalib_N_Zero[startSweep], myData$GRF_APSumCalib_N_Zero[endSweep])

  par(mfrow=c(1,3), oma = c(0, 0, 2, 0))  # oma = outer margin with 2 lines above the top of the graphs
  # Vertical component of GRF graph
  plot(myData$sweep[plotStart:plotEnd], myData$GRF_vertSumCalib_N_Zero[plotStart:plotEnd], xlab='Sweep', ylab='GRF - Vertical (N)', main='Zeroed GRF (Vertical) Force', type="l", col="blue")
  points(ImpPoints_X[1,], ImpPoints_GRFVert[1,], type='p', pch='O', col='cyan')
  text(ImpPoints_X[1,], ImpPoints_GRFVert[1,], labels=names(ImpPoints_X), pos=3, font=2) # pos: 1 = below, 2 = left, 3 = above, 4 = right
  
  # Mediolateral component of GRF graph
  plot(myData$sweep[plotStart:plotEnd], myData$GRF_MLSumCalib_N_Zero[plotStart:plotEnd], xlab='Sweep', ylab='GRF - Mediolateral (N)', main='Zeroed GRF (Mediolateral) Force', type="l", col="red")
  points(ImpPoints_X[1,], ImpPoints_GRFML[1,], type='p', pch='O', col='cyan')
  text(ImpPoints_X[1,], ImpPoints_GRFML[1,], labels=names(ImpPoints_X), pos=3, font=2) # pos: 1 = below, 2 = left, 3 = above, 4 = right
  
  # Anteroposterior component of GRF graph
  plot(myData$sweep[plotStart:plotEnd], myData$GRF_APSumCalib_N_Zero[plotStart:plotEnd], xlab='Sweep', ylab='GRF - Anteroposterior (N)', main='Zeroed GRF (Anteroposterior) Force', type="l", col="forestgreen")
  points(ImpPoints_X[1,], ImpPoints_GRFAP[1,], type='p', pch='O', col='cyan')
  text(ImpPoints_X[1,], ImpPoints_GRFAP[1,], labels=names(ImpPoints_X), pos=3, font=2) # pos: 1 = below, 2 = left, 3 = above, 4 = right
  
  mtext(filename, line=0.5, outer=TRUE)  # writes an overall title over the graphs
  
  # Prepping the data to be filtered
  cycleSweeps <- myData$sweep[startSweep:endSweep]
  cycleGRFVert <- myData$GRF_vertSumCalib_N_Zero[startSweep:endSweep]
  cycleGRFML <- myData$GRF_MLSumCalib_N_Zero[startSweep:endSweep]
  cycleGRFAP <- myData$GRF_APSumCalib_N_Zero[startSweep:endSweep]
  cycleGRFTime <- myData$Time_s[startSweep:endSweep]
  filterPrep <- data.frame(cycleSweeps,cycleGRFVert, cycleGRFML, cycleGRFAP, cycleGRFTime)
  names(filterPrep) <- c('sweep', 'GRF0SumVN', 'GRF0SumMLN', 'GRF0SumAPN', 'Time_s')
  
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
  if (is.null(BW) == TRUE) {
    output <- list(
      filterPrep = filterPrep,
      filename = filename,
      allData = myData,
      startSweep = startSweep,
      endSweep = endSweep
    )
  } else {
    BW_N <- BW*9.8 
    filterPrep_BW <- data.frame(sweep = filterPrep[,1], (filterPrep[,2:4]/BW_N), Time_s = filterPrep[,5])
    output <- list(
      filterPrep = filterPrep,
      filterPrep_BW = filterPrep_BW,
      filename = filename,
      allData = myData,
      startSweep = startSweep,
      endSweep = endSweep
    )
  }

}

#### butterFilteR ####
#' @title Apply a butterworth filter to data
#'
#' @name butterFilteR
#'
#' @description Applies a butterworth filter to data, using information about the data to determine what polynomial to use.
#'
#' @usage 
#'
#' @param \code{df} A list containing the file name and an array of data containing the independent variable (time) as the first column that is followed by the dependent variables. Can take voltToForce output as an input.
#' @param \code{Fs} A numeric value indicating the sampling frequency. 5000 Hz is set as a default.
#' @param \code{PbF} A numeric value indicating the pass-band frequency. 6 is set as a default.
#' @param \code{SbF} A numeric value indicating the stop-band frequency. 190 is set as a default.
#' @param \code{Rp} A numeric value indicating passband ripple in dB; represents the max permissible passband loss. 2 dB is set as a default.
#' @param \code{Rs} A numeric value indicating stopband attenuation in dB; respresents the dB the stopband is down from the passband. 40 dB is set as a default.
#' @param \code{saveAs} A character string containing the name for the resulting graph.
#' @param \code{saveGraph} A character string (options = "yes" or "no") to indicate whether you want to save the graph as a PDF or not.
#' @details These procedures follow the methodology used in Kawano and Blob (2013) and Kawano et al. 2016. It is assumed that the output from the force platform contain 12 channels in the following order: trigger, four verticals, sum of the verticals, two mediolateral, sum of the mediolaterals, two anteroposterior, and the sum of the anteroposteriors.
#' @references Kawano SM, Blob RW. 2013. Propulsive forces of mudskipper fins and salamander limbs during terrestrial locomotion: implications for the invasion of land. Integrative and Comparative Biology 53(2): 283-294. \url{https://academic.oup.com/icb/article/53/2/283/806410/Propulsive-Forces-of-Mudskipper-Fins-and}
#' @references Kawano SM, Economy DR, Kennedy MS, Dean D, Blob RW. 2016. Comparative limb bone loading in the humerus and femur of the tiger salamander Ambystoma tigrinum: testing the "mixed-chain" hypothesis for skeletal safety factors. Journal of Experimental Biology 219: 341-353. \url{http://jeb.biologists.org/content/219/3/341}
#'
#' @examples
#' 
#' GRF <- voltToForce(df, calib = c(-0.710, 1.337, 1.563), lightStartFrame = 248, startFrame = 20, endFrame = 196, zeroStart = 22000)
#' GRF_filtered <- butterFilteR(GRF, saveAs = "af01f18_Pec_Filter.pdf", saveGraph = "yes")
#'
#' @import signal
#' @export
#' 

### Include description of output?

butterFilteR <- function(df, Fs = 5000, PbF = 6, SbF = 190, Rp = 2, Rs = 40, saveAs = NULL, saveGraph = c("yes", "no"), ...) {
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
  pad <- data.frame(df[[1]][,1], rev(df[[1]][,2]), rev(df[[1]][,3]), rev(df[[1]][,4]), df[[1]][,5])
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
  
  if (saveGraph == "yes") {
  # Plotting the data
  pdf(saveAs)

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
  dev.off()
  }
  
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
  FilteredAll <- data.frame(sweep = df$filterPrep[,1], filterVN, filterMLN, filterAPN, Time_s = df$filterPrep[,5])
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



#### jointAngle ####
#' @title Calculate joint angles (in degrees) from XYZ coordinate data
#'
#' @name jointAngle
#'
#' @description Calculates the angle of a joint (in degrees), formed by three points with XYZ coordinates.
#'
#' @usage  jointAngle(P1, P2, P3)
#'
#' @param \code{P1} A data.frame of numeric values containing the X, Y, and Z coordinate data, respectively, for the first point (e.g, the shoulder)
#' @param \code{P2} A data.frame of numeric values containing the X, Y, and Z coordinate data, respectively, for the second point, which is assumed to be the vertex of the angle (e.g, the elbow)
#' @param \code{P3} A data.frame of numeric values containing the X, Y, and Z coordinate data, respectively, for the third point (e.g, the wrist)

#' @details These procedures follow the methodology used in Kawano and Blob (2013) and Kawano et al. 2016 to calculate angles formed about the limb joints in animals. 
#' @references Kawano SM, Blob RW. 2013. Propulsive forces of mudskipper fins and salamander limbs during terrestrial locomotion: implications for the invasion of land. Integrative and Comparative Biology 53(2): 283-294. \url{https://academic.oup.com/icb/article/53/2/283/806410/Propulsive-Forces-of-Mudskipper-Fins-and}
#' @references Kawano SM, Economy DR, Kennedy MS, Dean D, Blob RW. 2016. Comparative limb bone loading in the humerus and femur of the tiger salamander Ambystoma tigrinum: testing the "mixed-chain" hypothesis for skeletal safety factors. Journal of Experimental Biology 219: 341-353. \url{http://jeb.biologists.org/content/219/3/341}
#' @references "R - Comute Cross Product of Vectors (Physics)" answer posted by user Kevin on April 22, 2016. \url{https://stackoverflow.com/questions/36798301/r-compute-cross-product-of-vectors-physics}
#' @references "how to calculate the Euclidean norm of a vector in R?" answer posted by user joran on June 7, 2012. \url{https://stackoverflow.com/questions/10933945/how-to-calculate-the-euclidean-norm-of-a-vector-in-r}
#' 
#' @examples
#' 
#' shoulder <- matrix(c(0.006305306, 0.006526961, 0.006747555, -0.08206114, -0.08207707, -0.08207049, 0.006997669, 0.006980824, 0.006975157), 3, 3)
#' elbow <- matrix(c(0.007826633, 0.007959096, 0.008068560, -0.07183020, -0.07185459, -0.07186337, 0.005754819, 0.005764666, 0.005774707), 3, 3)
#' wrist <- matrix(c(0.01164592, 0.01160690, 0.01157642, -0.07348876, -0.07345559, -0.07342105, -0.000631402, -0.000671288, -0.000709513), 3, 3)
#' 
#' elbowAngle <- jointAngle(P1 = shoulder, P2 = elbow, P3 = wrist)
#'
#' @export


jointAngle <- function(P1, P2, P3, ...) {
  #### Create vectors of the XYZ coordinates for each anatomical landmark ####
  # Assume that P2 is the vertex and P1 and P3 are the other two points that form the angle
  # data are storted in df
  # Most proximal joint (e.g., Shoulder, S) XYZ are columns 1 - 3, 
  # Middle joint (e.g., elbow, E) XYZ are 4-6, and 
  # most distal joint (e.g., wrist, W) XYZ are 7 - 9
  # Only focusing on first row for now
  
  # norm is the Euclidean norm of a vector, which can be calculated in R using
  Euc_norm <- function(x) sqrt(sum(x^2))
  
  # To produce a vector cross-product: 
  vec_cross <- function(ab,ac){
    abci = ab[2] * ac[3] - ac[2] * ab[3];
    abcj = ac[1] * ab[3] - ab[1] * ac[3];
    abck = ab[1] * ac[2] - ac[1] * ab[2];
    return (c(abci, abcj, abck))
  }
  
  # P32 <- matrix(0, 3, 3)
  # P12 <- matrix(0, 3, 3)
  
  P32 <- matrix(0, nrow(P2), 3)
  P12 <- matrix(0, nrow(P2), 3)
  angle_degrees <- numeric(length = nrow(P1))
  
  for (i in 1:nrow(P1)) {
    if (isTRUE(is.vector(P1))) {
      # Create vectors of the P3-P2 and P1-P2 segments
      P32 <- as.numeric(P3 - P2)
      P12 <- as.numeric(P1 - P2)
      
      angle_degrees <- atan2(Euc_norm(vec_cross(P32, P12)), sum(P32 * P12))*(180/pi)
      
    } else {
      # Create vectors of the P3-P2 and P1-P2 segments
      P32[i,] <- as.numeric(P3[i,] - P2[i,])
      P12[i,] <- as.numeric(P1[i,] - P2[i,])
      
      angle_degrees[i] <- atan2(Euc_norm(vec_cross(P32[i,], P12[i,])), sum(P32[i,] * P12[i,]))*(180/pi)
    }
  }
  return(angle_degrees)
}



#### yaw ####
#' @title Calculate yaw (in degrees) from XY coordinate data
#'
#' @name yawAngle
#'
#' @description Calculates the yaw angle (in degrees), formed by two segments produced by the X and Y coordiantes of two points.
#'
#' @usage yawAngle(P1, P2)
#'
#' @param \code{P1} A data.frame of numeric values containing the X and Y coordinate data, respectively, for the first point (e.g, point along the midline)
#' @param \code{P2} A data.frame of numeric values containing the X and Y coordinate data, respectively, for the second point (e.g., pelvic girdle)

#' @details These procedures follow the methodology used in Kawano and Blob (2013) and Kawano et al. 2016 to calculate angles formed about the limb joints in animals. 
#' @references Kawano SM, Blob RW. 2013. Propulsive forces of mudskipper fins and salamander limbs during terrestrial locomotion: implications for the invasion of land. Integrative and Comparative Biology 53(2): 283-294. \url{https://academic.oup.com/icb/article/53/2/283/806410/Propulsive-Forces-of-Mudskipper-Fins-and}
#' @references Kawano SM, Economy DR, Kennedy MS, Dean D, Blob RW. 2016. Comparative limb bone loading in the humerus and femur of the tiger salamander Ambystoma tigrinum: testing the "mixed-chain" hypothesis for skeletal safety factors. Journal of Experimental Biology 219: 341-353. \url{http://jeb.biologists.org/content/219/3/341}
 
#' @examples
#' 

#' P1 <- data.frame(X = c(0.004978444,
#'                     0.005226819,
#'                     0.005483069,
#'                     0.005747358,
#'                     0.006019848,
#'                     0.006300701,
#'                     0.006590101,
#'                     0.00688825,
#'                     0.007195354, 
#'                     0.007511617),
#'               Y = c(-0.09374413,
#'                     -0.09380451,
#'                     -0.09382974,
#'                     -0.09382203,
#'                     -0.09378361,
#'                     -0.09371668,
#'                     -0.09362362,
#'                     -0.09350694,
#'                     -0.09336917, 
#'                     -0.09321286
#'               ))
#'
#'P2 <- data.frame(X = c(0.001005251,
#'                     0.001195392,
#'                     0.001406054,
#'                     0.001636047,
#'                     0.001884185,
#'                     0.002149271,
#'                     0.002429912,
#'                     0.002724486,
#'                     0.003031363,
#'                     0.003348911),
#'               Y = c(-0.09343679,
#'                     -0.09347614,
#'                     -0.09350069,
#'                     -0.09350983,
#'                     -0.09350293,
#'                     -0.09347937,
#'                     -0.09343846,
#'                     -0.09337947,
#'                     -0.09330165,
#'                     -0.09320426
#'               ))
#'
#' Yaw <- yawAngle(P1, P2)
#' 
#' @export

yawAngle <- function(P1, P2, ...) {
  SegVectorBackY <- P1[,2] - P2[,2]
  SegVectorBackX <- P1[,1] - P2[,1]
  TanYaw <- SegVectorBackY/SegVectorBackX
  yaw <- (atan(TanYaw))*(180/pi)
  return(yaw)
}
  


#### Protraction versus retraction ####
#' @title Calculate protraction versus retraction (in degrees) from XYZ coordinate data
#'
#' @name protraction
#'
#' @description Calculates the protraction angle (in degrees), of a limb based on the XYZ coordinates of three points along the limb. 
#'
#' @usage protraction(P1, P2, P3, yaw)
#'
#' @param \code{P1} A data.frame of numeric values containing the X, Y, and Z coordinate data, respectively, for the first point (e.g, the hip)
#' @param \code{P2} A data.frame of numeric values containing the X, Y, and Z coordinate data, respectively, for the second point, which is assumed to be the vertex of the angle (e.g, the knee)
#' @param \code{P3} A data.frame of numeric values containing the X, Y, and Z coordinate data, respectively, for the third point (e.g, the ankle)
#' @param \code{Yaw} A vector of numeric values containing the yaw angle values

#' @details These procedures follow the methodology used in Kawano and Blob (2013) and Kawano et al. 2016 to calculate angles formed about the limb joints in animals. 
#' @references Kawano SM, Blob RW. 2013. Propulsive forces of mudskipper fins and salamander limbs during terrestrial locomotion: implications for the invasion of land. Integrative and Comparative Biology 53(2): 283-294. \url{https://academic.oup.com/icb/article/53/2/283/806410/Propulsive-Forces-of-Mudskipper-Fins-and}
#' @references Kawano SM, Economy DR, Kennedy MS, Dean D, Blob RW. 2016. Comparative limb bone loading in the humerus and femur of the tiger salamander Ambystoma tigrinum: testing the "mixed-chain" hypothesis for skeletal safety factors. Journal of Experimental Biology 219: 341-353. \url{http://jeb.biologists.org/content/219/3/341}

#' @examples
#' 

#' shoulder <- matrix(c(0.006305306, 0.006526961, 0.006747555, -0.08206114, -0.08207707, -0.08207049, 0.006997669, 0.006980824, 0.006975157), 3, 3)
#' elbow <- matrix(c(0.007826633, 0.007959096, 0.008068560, -0.07183020, -0.07185459, -0.07186337, 0.005754819, 0.005764666, 0.005774707), 3, 3)
#' wrist <- matrix(c(0.01164592, 0.01160690, 0.01157642, -0.07348876, -0.07345559, -0.07342105, -0.000631402, -0.000671288, -0.000709513), 3, 3)
#' yaw <- c(-4.4232170, -4.6566089, -4.6142598)
#'
#' protractionAngle <- protraction(shoulder, elbow, wrist, yaw)
#' 
#' @importFrom pracma cross
#' @export

protraction <- function(P1,P2,P3, Yaw, ...) {
  HipPoint <- P1
  KneePoint <- P2
  AnklePoint <- P3
  
  FemurVector <- HipPoint - KneePoint
  HipPointTrans <- HipPoint - HipPoint
  KneePointTrans1A <- KneePoint - HipPoint 
  KneePointTrans1B <- KneePoint - HipPoint
  FemurVectorTrans1A <- HipPointTrans - KneePointTrans1A
  FemurVectorTrans1B <- HipPointTrans - KneePointTrans1B
  
  TibiaVector <- KneePoint - AnklePoint
  KneePointTrans2A <- KneePoint - KneePoint
  AnklePointTransA <- AnklePoint - KneePoint
  KneePointTrans2B <- KneePoint - KneePoint
  AnklePointTransB <- AnklePoint-KneePoint
  KneePointTrans2C <- KneePoint-KneePoint
  AnklePointTransC <- AnklePoint-KneePoint
  TibiaVectorTransA <- KneePointTrans2A - AnklePointTransA
  TibiaVectorTransB <- KneePointTrans2B-AnklePointTransB
  TibiaVectorTransC <- KneePointTrans2C-AnklePointTransC
  
  ## dot product between two row vectors
  wdot <- function(a, b) {
    y <- a*b
    y <- t(sum(t(y)))
    return(y)
  }
  
  ## Cross rpdouct between two row vectors 
  library(pracma) # for the cross() function that returns a vector
  wcross <- function(a, b) {
    c <- t(pracma::cross(t(a),t(b)))
    return(c)
  }
  
  ## vlength
  vlength <- function(x) {
    v <- (wdot(x, x))^0.5
    return(v)
  }
  
  ## Setting up Transverse Plane
  TVPlane1 <- c(0, 0, 0)
  TVPlane2 <- c(0, .1, 0)
  TVPlane3 <- c(0, 0, .1)
  TVVector1 <- TVPlane2 - TVPlane1
  TVVector2 <- TVPlane3 - TVPlane1
  # had to take the transpose to get it into row-form instead of column form
  TVNorm <- t(wcross(TVVector1,TVVector2))
  
  FemTVAngInit <- matrix()
  for (i in 1:nrow(KneePointTrans1A)) {
    FemurVectorTrans1AA <- FemurVectorTrans1A[i,]
    dotFemurVectorTV <- wdot(FemurVectorTrans1AA,TVNorm)
    MagFemurVectorTV <- vlength(FemurVectorTrans1AA)
    MagTVNorm <- vlength(TVNorm)
    MagFemurVectorTVNorm <- MagFemurVectorTV*MagTVNorm
    CosdotFemurVectorTV <- dotFemurVectorTV/MagFemurVectorTVNorm
    FemTVAngA <- (acos(CosdotFemurVectorTV))*RAD2DEG
    # %makes a one column matrix of angles of the protraction/retraction angle
    FemTVAngInit[i] <- FemTVAngA
  }
  
    FemTVAng <- abs(FemTVAngInit)
    #setting the zero of the angles to be perpendicular to x axis
    FemTVAng <- (90-FemTVAng)
    FemTVAng <- (90-FemTVAng)+Yaw -90
    
    return(FemTVAng)
  }




#### Pitch ####
#' @title Calculate pitch angle (in degrees) from XYZ coordinate data of two points
#'
#' @name pitchAngle
#'
#' @description Calculates the pitch angle (in degrees), from the XYZ coordinates of two points 
#'
#' @usage pitchAngle(P1, P2)
#'
#' @param \code{P1} A data.frame of numeric values containing the X, Y, and Z coordinate data, respectively, for the first point (e.g, the pelvic girdle)
#' @param \code{P2} A data.frame of numeric values containing the X, Y, and Z coordinate data, respectively, for the second point (e.g., the hip)

#' @details These procedures follow the methodology used in Kawano and Blob (2013) and Kawano et al. 2016 to calculate angles formed about the limb joints in animals. 
#' @references Kawano SM, Blob RW. 2013. Propulsive forces of mudskipper fins and salamander limbs during terrestrial locomotion: implications for the invasion of land. Integrative and Comparative Biology 53(2): 283-294. \url{https://academic.oup.com/icb/article/53/2/283/806410/Propulsive-Forces-of-Mudskipper-Fins-and}
#' @references Kawano SM, Economy DR, Kennedy MS, Dean D, Blob RW. 2016. Comparative limb bone loading in the humerus and femur of the tiger salamander Ambystoma tigrinum: testing the "mixed-chain" hypothesis for skeletal safety factors. Journal of Experimental Biology 219: 341-353. \url{http://jeb.biologists.org/content/219/3/341}

#' @examples
#' 

#' P1 <- matrix(c(0.001005251, 0.001195392, 0.001406054, -0.09343679, -0.09347614, -0.09350069, 0.01603717, 0.0161097, 0.01616693), 3, 3)
#' P2 <- matrix(c(0.006305306, 0.006526961, 0.006747555, -0.08206114, -0.08207707, -0.08207049, 0.006997669, 0.006980824, 0.006975157), 3, 3)
#' 
#' pitch <- pitchAngle(P1, P2)
#' 
#' @importFrom pracma cross
#' @export


pitchAngle <- function(P1, P2, ...) {
  SacPointX <- P1[,1]
  SacPointZ <- P1[,3]
  HipPointX <- P2[,1]
  HipPointZ <- P2[,3]
  ShiftPelvisPlaneX <- SacPointX - HipPointX
  ShiftPelvisPlaneZ <- SacPointZ - HipPointZ
  
  ## dot product between two row vectors
  wdot <- function(a, b) {
    y <- a*b
    y <- t(sum(t(y)))
    return(y)
  }
  
  ## Cross rpdouct between two row vectors 
  library(pracma) # for the cross() function that returns a vector
  wcross <- function(a, b) {
    c <- t(pracma::cross(t(a),t(b)))
    return(c)
  }
  
  ## vlength
  vlength <- function(x) {
    v <- (wdot(x, x))^0.5
    return(v)
  }
  
  ## ABDUCTION/ADDUCTION OF FEMUR: MOVEMENT IN THE DORSOVENTRAL DIRECTION ABSOLUTE FRAME OF REFERENCE
  
  # Roll calculations
  
  # Definitions of normal vectors to planes for calculation of segment angles to those vectors &, thereby, those planes
  
  # define horizontal plane and normal to it:  (wcross equals a cross product function written by Will Corning for convenient orientation)
  HZPlane1 <- c(0, 0, 0)
  HZPlane2 <- c(.1, 0, 0)
  HZPlane3 <- c(0, .1, 0)
  HZVector1 <- HZPlane2 - HZPlane1
  HZVector2 <- HZPlane3 - HZPlane1
  # had to take the tranpose to get this to work with wdot
  HZNorm <- t(wcross(HZVector1,HZVector2))

  # define pitch-adjusted horizontal plane and normal to it:
  # define plane of body based on sacral point and hip point, defining based on ilium rather than trunk
  
  HZPlane2Pitch <- matrix(NA, length(ShiftPelvisPlaneX), 3)
  HZNormPitch <- matrix(NA, length(ShiftPelvisPlaneX), 3)
  HZNormPitchMatrix <- matrix(NA, length(ShiftPelvisPlaneX), 3)
  MagHZNormPitch <- matrix(NA, length(ShiftPelvisPlaneX), 1)
  MagPitchVector <- matrix(NA, length(ShiftPelvisPlaneX), 1)
  CosdotPitch <- matrix(NA, length(ShiftPelvisPlaneX), 1)
  PitchAng <- matrix(NA, length(ShiftPelvisPlaneX), 1)
  
  for (i in 1:length(ShiftPelvisPlaneX)) {
    HZPlane1Pitch <- c(0, 0, 0)
    HZPlane2Pitch[i,] <- c(ShiftPelvisPlaneX[i], 0, ShiftPelvisPlaneZ[i])
    HZPlane3Pitch <- c(0, .1, 0)
    HZVector1Pitch <- HZPlane2Pitch - HZPlane1Pitch
    HZVector2Pitch <- HZPlane3Pitch - HZPlane1Pitch
    HZNormPitch[i,] <- wcross(HZVector1Pitch[i,], HZVector2Pitch)
    HZNormPitchMatrix[i,] <- HZNormPitch[i,]
    
    dotPitch <- wdot(HZNormPitch[i,],HZNorm)
    MagHZNormPitch[i] <- vlength(HZNormPitch[i,])
    MagHZNorm <- vlength(HZNorm)
    MagPitchVector[i] <- MagHZNormPitch[i,]*MagHZNorm
    CosdotPitch[i] <- dotPitch/MagPitchVector[i]
    PitchAng[i] <- (acos(CosdotPitch[i]))*(180/pi)
    Pitch <- PitchAng
  }
  return(Pitch)
}





############# GRFangles ##########
## assumes that your time component is in the first column and is followed by three columns containing the components of the GRF
## assumes that the order is Vert, ML, and AP
## must be data.frame

#' @title Calculate angles of orientation for ground reaction forces
#'
#' @name GRFAngles
#'
#' @description Calculates the net ground reaction force (GRF) and angles of the GRF orientation using information about the vertical, mediolaterl, and anteroposterior components of the GRF.
#'
#' @usage GRFAngles <- function(myData, ...)
#'
#' @param \code{myData} A data.frame of four columns of numeric data in the following order: Measure of time (e.g., seconds, percent of stance), vertical component of the GRF, mediolateral component of the GRF, and anteroposterior component of the GRF.
#' @details These procedures follow the methodology used in Kawano and Blob (2013) and Kawano et al. 2016.
#' @references Kawano SM, Blob RW. 2013. Propulsive forces of mudskipper fins and salamander limbs during terrestrial locomotion: implications for the invasion of land. Integrative and Comparative Biology 53(2): 283-294. \url{https://academic.oup.com/icb/article/53/2/283/806410/Propulsive-Forces-of-Mudskipper-Fins-and}
#' @references Kawano SM, Economy DR, Kennedy MS, Dean D, Blob RW. 2016. Comparative limb bone loading in the humerus and femur of the tiger salamander Ambystoma tigrinum: testing the "mixed-chain" hypothesis for skeletal safety factors. Journal of Experimental Biology 219: 341-353. \url{http://jeb.biologists.org/content/219/3/341}
#'
#' @examples
#' 
#' GRFAngles(GRFs_Filtered$GRF0Sum_filter_interp)
#'
#' @export
#' 

GRFAngles <- function(myData, ...) {
  # Calculating the angles of orientation
  myData$GRF_Vert_Sq <- myData[,2]^2
  myData$GRF_ML_Sq <- myData[,3]^2
  myData$GRF_AP_Sq <- myData[,4]^2
  myData$NetGRF_N <- sqrt(myData$GRF_Vert_Sq + myData$GRF_ML_Sq + myData$GRF_AP_Sq)
  myData$MLAngle <- (acos(myData[,3]/(sqrt(myData$GRF_ML_Sq + myData$GRF_Vert_Sq))))*(180/pi)
  myData$MLAngle_Convert <- 90 - myData$MLAngle
  myData$APAngle <- (acos(myData[,4]/(sqrt(myData$GRF_AP_Sq + myData$GRF_Vert_Sq))))*(180/pi)
  myData$APAngle_Convert <- 90 - myData$APAngle
  
  output <- data.frame(myData[,1:5], "NetGRF_N" = myData$NetGRF_N, "MLAngle_deg" = myData$MLAngle, "MLAngle_Convert_deg" = myData$MLAngle_Convert, "APAngle_deg" = myData$APAngle, "APAngle_Convert_deg" = myData$APAngle_Convert)
  
}


############# removeOverlaps ##########


#' @title Remove sections of the data where there is overlap
#'
#' @name removeOverlaps
#'
#' @description Identifies portions of the dataset where two occurrences overlap and replaces with NAs.
#'
#' @usage removeOverlaps <- function(df, primary, secondary, filmRate, ...) 
#'
#' @param \code{df} A data.frame containing variables of interest.
#' @param \code{primary} A data.frame containing two numeric values that correspond to the beginning and end of an event (e.g., stance phase of forelimb).
#' @param \code{secondary} A data.frame containing two numeric values that correspond to the beginning and end of an event (e.g., stance phase of hind limb) that would affect interpretations of the primary event.
#' @param \code{filmRaate} A numeric value of the filming rate. 
#' @details These procedures follow the methodology used in Kawano and Blob (2013) and Kawano et al. 2016.
#' @references Kawano SM, Blob RW. 2013. Propulsive forces of mudskipper fins and salamander limbs during terrestrial locomotion: implications for the invasion of land. Integrative and Comparative Biology 53(2): 283-294. \url{https://academic.oup.com/icb/article/53/2/283/806410/Propulsive-Forces-of-Mudskipper-Fins-and}
#' @references Kawano SM, Economy DR, Kennedy MS, Dean D, Blob RW. 2016. Comparative limb bone loading in the humerus and femur of the tiger salamander Ambystoma tigrinum: testing the "mixed-chain" hypothesis for skeletal safety factors. Journal of Experimental Biology 219: 341-353. \url{http://jeb.biologists.org/content/219/3/341}
#'
#' @examples
#' 
#' primary <- data.frame(35, 148)
#' secondary <- data.frame(128, 244)
#' removeOverlaps(Pel_GRFs_Filtered_dataset, primary, secondary, 100)
#'
#' @export
#' 


removeOverlaps <- function(df, primary, secondary, filmRate, ...) {
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



############# Yank ##########


#' @title Calculating yank for force data
#'
#' @name yank
#'
#' @description Calculates yank as a measure of the change in force per unit of time.
#'
#' @usage yank <- function(myData, ...)
#'
#' @param \code{time} Numeric values representing the time elapsed in units of seconds.
#' @param \code{force} Numeric values representing the force measured at each time point.
#' @details These procedures follow the methodology used in Lin et al. 2019.
#' @references Lin DC, McGowan CP, Blum KP, Ting LH. Yank: the time derivative of force is an important biomechanical variable in sensorimotor systems. Journal of Experimental Biology 222: jeb180414. \url{https://jeb.biologists.org/content/222/18/jeb180414}
#'
#' @examples
#' 
#' yank(seq(1:10), rnorm(10, 1, 1), )
#'
#' @export
#' 

yank <- function(time, force, ...) {
  dt <- c(NA, diff(time))
  df <- c(NA, diff(force))
  dy <- df/dt
  dataset <- data.frame(time = time, force = force, delta_force = df, yank = dy)
  return(dataset)
}


############################################################################################################################
######################################## CALCULATING CRUDE SPEED FOR A SINGLE FILE #########################################
#####                                                                                                                  #####
##### CrudeSpd is used to calculate a crude estimate of speed, whereby speed is calculated based on the                #####
#####   distance traveled over a given time.  The assumption is that a matrix of X,Y coordinates are entered and       #####
#####   distance is calculated based on the first and last row of these matrices.                                      #####
#####                                                                                                                  #####
##### - x = 2-column matrix that contain X-coordinates in column 1 and Y-coordinates in column 2                       #####
##### - rate = filming rate in units of frames per second                                                              #####
##### - calib = value to calibrate the X,Y coordinate data from units of pixels to another unit                        #####
#####   (e.g., cm, mm); if no vector is provided, output will not be calibrated and will remain in the                 #####
#####   original units used                                                                                            #####
##### - distance.scale = character string of the units of measure related to distance traveled                         #####
##### - time.scale = character string of the units of time. Default is seconds ("s")                                   #####
#####                                                                                                                  #####
#####                                                                                                                  #####
##### Note: If you are evaluating more than one trial, use CrudeSpd.Multi() instead                                    #####    
#####                                                                                                                  #####
##### This function will output:                                                                                       #####
##### 1. time elapsed                                                                                                  #####
##### 2. distance traveled (uncalibrated)                                                                              #####
##### 3. distance traveled (calibrated)                                                                                #####
##### 4. tangential velocity (calibrated)                                                                              #####  
##### 5. tangential speed (calibrated)                                                                                 ##### 
##### 6. names of your observations/rows (optional; only outputted if data entered for rownames)                       #####  
#####                                                                                                                  #####
#####                                                                                                                  #####
############################################################################################################################

CrudeSpd <- function(x, rate, calib=1, distance.scale="units", time.scale="s") {
  xy <- matrix(NA, 1, 5)
  xy[,1] <- nrow(x)*(1/rate)
  xy[,2] <- sqrt(((x[nrow(x),1]-x[1,1])^2)+((x[nrow(x),2]-x[1,2])^2))		
  xy[,3] <- xy[,2]/calib
  xy[,4] <- xy[,3]/xy[,1]
  xy[,5] <- abs(xy[,4])
  xy <- data.frame(xy)
  names(xy) <- c(paste("Time",time.scale, sep="."), "Distance.Uncalib", paste("Distance",distance.scale,sep="."), 
                 paste("Velocity",distance.scale,time.scale, sep="."), paste("Speed",distance.scale,time.scale, sep="."))
  return(xy)
}