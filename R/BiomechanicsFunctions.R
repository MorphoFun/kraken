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
#' @description The following function calculates the summary statistics for groups that have multiple observations and then plots the mean and confidence interval for each group. The confidence interval is based on the lower and upper Gaussian confidence limits based on the t-distribution.
#'
#' @usage profilePlotR(d = d, xvar = xvar, yvar = yvar, groupname = groupname, xlab = "x", ylab = "y", colorPalette = c("#D55E00", "#0072B2", "#56B4E9"), linetypes = NULL, grouplevels = NULL, title = NULL, ...)
#'
#' @param \code{d} data (currently only accepts input for one variable at a time).
#' @param \code{xvar} x-axis variable name.
#' @param \code{yvar} y-axis variable name.
#' @param \code{groupname} variable name for the overall group that is being evaluated (e.g., species).
#' @param \code{xlab} character string for the x-axis label.
#' @param \code{ylab} character string for the y-axis label.
#' @param \code{colorPalette} character string of colors to use for plot.
#' @param \code{linetypes} character string of lines to use for each group. 
#' @param \code{grouplevels} character string of the text to print for each group in the legend.
#' @param \code{title} character string for the title of the plot.
#'
#' @details Function to quickly generate profile plots for data. For instance, kinematic plots over time for multiple individuals that have multiple trials of data collected.

#'
#' @examples
#' percentStance = rep(seq(1,5),4)
#' yank = c(0.00062, 0.00172, 0.00269, 0.00346, 0.00412, 0.0022, 0.00072, 0.00169, 0.00246, 0.00312, 0.00028, 0.00084, 0.00151, 0.00239, 0.00340, 0.00041, 0.00122, 0.00202, 0.00277, 0.00341)
#' species = c("ab" , "ab", "ab", "ab", "ab", "ab", "ab", "ab", "ab", "ab", "cd", "cd", "cd", "cd", "cd", "cd", "cd", "cd", "cd", "cd")
#' df <- data.frame(percentStance, yank, species)
#' cbPalette <- c("#D55E00", "#0072B2")
#' grouplevels = c("Aa bb", "Cc dd")
#' profilePlotR(df, "percentstance", "yank", "species", "Percent Stance", "Yank (BW per sec)", colorpalette = cbPalette, grouplevels = grouplevels)
#'
#' @import ggplot2
#' 
#' @export

profilePlotR <- function(d = d, xvar = xvar, yvar = yvar, groupname = groupname, xlab = "x", ylab = "y", colorPalette = c("#D55E00", "#0072B2", "#56B4E9"), linetypes = NULL, grouplevels = NULL, title = NULL, xrange = NULL, yrange = NULL, ...) {
  if(is.null(grouplevels)) {grouplevels = unique(d[,groupname])}
  if(is.null(linetypes)) {linetypes = 1:(length(unique(d[,groupname])))}
  ggplot(d, aes_string(x= xvar, y = yvar)) + 
    scale_y_continuous(limits = c(yrange[1], yrange[2]),  paste(ylab, "\n")) +
    scale_x_continuous(limits = c(xrange[1], xrange[2]), paste("\n ", xlab)) +
    geom_rect(aes(xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=0), fill="gray95") +
    stat_summary(aes_string(linetype = groupname), fun = mean, geom = 'line', size=1, alpha=0.9) +
    stat_summary(aes_string(fill = groupname), fun.data = mean_cl_normal, geom = "ribbon", fun.args = list(mult = 1), alpha = 0.5) + 
    scale_color_manual(name = groupname, # changing legend title
                       labels = grouplevels, # Changing legend labels
                       values = colorPalette) +
    scale_fill_manual(name = groupname,
                      labels = grouplevels,
                      values = colorPalette) +
    scale_linetype_manual(name = groupname,
                          labels = grouplevels,
                          values = linetypes) +
    theme(axis.title.x=element_text(colour="black", size = 25))+ 
    theme(axis.title.y=element_text(colour='black', size = 25))+
    theme(axis.text.x=element_text(colour='black', size = 20))+
    theme(axis.text.y=element_text(colour='black', size = 20))+
    theme_classic() + 
    theme(legend.position="bottom", legend.direction="horizontal", legend.text = element_text(size = 12), legend.title = element_text(size = 15))
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
#' @usage GRFAngles(myData, ...)
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
#' @usage removeOverlaps(df, primary, secondary, filmRate, ...) 
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
#' @usage yank(time, force ...)
#'
#' @param \code{time} Numeric values representing the time elapsed in units of seconds.
#' @param \code{force} Numeric values representing the force measured at each time point.
#' @details These procedures follow the methodology described in Lin et al. 2019.
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


############# boneload_extinct ##########


#' @title Estimating peak stresses to the limb bones in extinct tetrapods
#'
#' @name boneload_extinct
#'
#' @description Biomechanical model to estimate the peak stresses experienced by the femur of a quadrupedal tetrapod during terrestrial locomotion. This computational model was developed to estimate peak limb bone stresses in extinct nonmammalian therapsids but could be applied to living tetrapods as well. GRF represents the ground reaction force. 
#'
#' @usage boneload_extinct(Anat, BW, minalpha, maxalpha, stepalpha, species = NULL, ...) 
#'
#' @param \code{Anat} Numeric values representing the anatomy of the study species. Anatomical variables include: mass_kg (in kilograms), rQuads_m (moment arm of the knee extensors, in units of meters), rAextAnk_m (moment arm of the ankle extensors about the ankle, in units of meters), rAextKnee_m (moment arm of the ankle extensors about the knee, in units of meters), thetaQuad_degrees (angle between extensor force and the knee extensor moment arm about the midshaft centroid, in units of degrees), LFemur_m (length of the femur, in units of meters), AFem_m2 (cross-sectional area of the femur, in units of meters squared), yDV_m (distance from neutral axis to cortex in the dorsoventral direction, in units of meters), IDV_m4 (second moment of area in the dorsoventral direction, in units of meters to the fourth power), RCFemDV_m (moment arm due to bone curvature in the dorsoventral direction, in units of meters), yAP_m (distance from neutral axis to cortex in the anteroposterior direction, in units of meters), IAP_m4 IDV_m4 (second moment of area in the anteroposterior direction, in units of meters to the fourth power), RCFemAP_m (moment arm due to bone curvature in the anteroposterior direction, in units of meters), RGRFKnee_m (moment arm of the GRF to the joint, in units of meters), LFoot_m (length of the foot, in units of meters), rQuadShaft_m (knee extensor moment arm about midshaft centroid, in units of meters)
#' @param \code{BW} Numeric values representing the magnitude of the net GRF in units of proportion of body mass (e.g., 0.5 for 50 percent of the animal's body mass).
#' @param \code{minalpha} Single numerical value for the minimum alpha value (angle between the GRF and limb bone, in units of degrees).
#' @param \code{maxalpha} Single numerical value for the maximum alpha value (angle between the GRF and limb bone, in units of degrees).
#' @param \code{stepalpha} Single numerical value for the increments of alpha (in units of degrees).
#' @param \code{species} Character string for the name of the species being tested. 
#' @details These procedures follow the methodology described in Blob 2001.
#' @references Blob, RW. 2001. The evolution of hindlimb posture in nonmammalian therapsids: biomechanical tests of paleontological hypotheses. Paleobiology 27(1): 14-38. \url{https://doi.org/10.1666/0094-8373(2001)027<0014:EOHPIN>2.0.CO;2}
#'
#' @examples
#' 
#' Anat <- data.frame(mass_kg = 3, rQuads_m = 0.0054315, rAextAnk_m = 0.00288353, rAextKnee_m = 0.008213, thetaQuad_degrees = 0, LFemur_m = 0.039931, AFem_m2 = 0.0000375, yDV_m = 0.00367, IDV_m4 = 0.000000000104, RCFemDV_m = 0.000501, yAP_m = 0.003377, IAP_m4 = 0.000000000127, RCFemAP_m = -0.001721, RGRFKnee_m = 0.012, LFoot_m = 0.062503, rQuadShaft_m = 0.021297, thetaQuadShaft_degrees = 14.495)
#' boneload_extinct(Anat, 0.5, 10, 70, 5, "Fakeosaurus") 
#'
#' @export
#' 

boneload_extinct <- function(Anat, BW, minalpha, maxalpha, stepalpha, species = NULL, ...) {
  DEG2RAD <- 180/pi
  
  ## Net theta calculations
  # Since zero divided by a numeric value is undefined, going to force it so that if theta Quad = zero in radians, it will equal zero in degrees as well
  # so it matches MATLAB output
  Anat$thetaQuad_radians <- Anat$thetaQuad_degrees/DEG2RAD
  Anat$thetaQuadShaft_radians <- Anat$thetaQuadShaft_degrees/DEG2RAD
  
  Anat$thetaQuadNet_radians <- Anat$thetaQuad_radians + Anat$thetaQuadShaft_radians
  
  ## GRF=Mass*9.807;note in line below, adjusted GRF to 0.4 of mass like salamanders for 
  ## Greererpeton; makes snese since the values would otherwise be >100MPa
  ## BW is the magnitude of the body weight of the GRF (e.g., 1 for gator or 0.4 for salamanders)
  Anat$BW = BW
  Anat$GRF <- Anat$BW*Anat$mass_kg*9.807
  Anat$RSTQuads <- Anat$rQuadShaft*sin(Anat$thetaQuadNet_radians)
  
  #### SELECTING ALPHA VALUES ####
  ## values are in units of degrees
  Anat$minalpha <- minalpha
  Anat$maxalpha <- maxalpha
  Anat$stepalpha <- stepalpha
  
  alpha <- seq(minalpha/DEG2RAD, maxalpha/DEG2RAD, stepalpha/DEG2RAD)
  nalpha <- length(alpha)
  
  #### CALCULATIONS ####
  
  ## Muscle force calculations
  # see Appendix 3 in Blob (2001) in Paleobiology
  # check around line 900 of the Ambystoma MATLAB bone loading code to see if that produces similar data for the regression
  # alternate lines to double RGRFKnee
  # RGRFKnee2=RGRFKnee*2;
  Anat$RGRFKnee2 <- Anat$RGRFKnee
  Anat$RrQuadsKnee <- Anat$RGRFKnee2/Anat$rQuads_m
  # RGRFAnkle=0.5*LFoot;
  # model RGRFAnkle based on regression increase of RGRFANkle/LFoot in lizards
  temp1 <- .242*DEG2RAD*alpha
  temp2 <- (-1*temp1)+40.361
  temp3 <- sin(temp2/DEG2RAD)
  temp4 <- temp3^2
  
  angles <- data.frame(alpha = alpha*DEG2RAD, RGRFAnkle = Anat$LFoot_m*temp4)
  
  # RGRFAnkle=LFoot*((sin(40.361-(.242*(RAD2DEG*alpha)))).^2); # I'm not sure what this extra line of code is for
  
  
  angles$RrAextAnk <- angles$RGRFAnkle/Anat$rAextAnk_m
  angles$FAext <- Anat$GRF*angles$RrAextAnk
  angles$FQuads <- (Anat$GRF*Anat$RGRFKnee2 + (angles$FAext*Anat$rAextKnee_m))/Anat$rQuads_m
  
  ## 
  angles$GRFaxFem <- Anat$GRF*cos(alpha) 
  angles$GRFtvFem <- Anat$GRF*sin(alpha); 
  angles$FQuadAx <- angles$FQuads*cos(Anat$thetaQuad_radians)
  angles$StressAxFem <- (-1)*(angles$GRFaxFem + angles$FQuadAx)/Anat$AFem_m2
  
  # Femoral bending moments due to GRF
  angles$BMomGRFtvDV <- (-1)*(Anat$LFemur_m*0.5*angles$GRFtvFem)
  angles$BMomGRFtvAP <- (Anat$LFemur_m*0.5*angles$GRFtvFem)
  angles$BMomGRFaxDV <- Anat$RCFemDV_m*angles$GRFaxFem
  angles$BMomGRFaxAP <- Anat$RCFemAP_m*angles$GRFaxFem
  
  
  # Femoral bending moments due to muscles calculated using cross products
  angles$QuadsBMom <- (-1)*angles$FQuads*Anat$RSTQuads
  angles$MuscStressDV <- angles$QuadsBMom*(Anat$yDV_m/Anat$IDV_m4)
  angles$NetGRFBMomDV <- angles$BMomGRFtvDV + angles$BMomGRFaxDV # if curv negative, augments GRF; if curv +, mitigates GRF
  angles$NetGRFBMomAP <- angles$BMomGRFtvAP + angles$BMomGRFaxAP # if curv negative, mitigates GRF; if curv +, augments GRF
  angles$NetGRFBStressDV <- angles$NetGRFBMomDV*(Anat$yDV_m/Anat$IDV_m4) # for dorsal surface
  angles$NetGRFBStressAP <- angles$NetGRFBMomAP*(Anat$yAP_m/Anat$IAP_m4) # for anterior surface
  angles$CrvGRFBStressDV <- angles$BMomGRFaxDV*(Anat$yDV_m/Anat$IDV_m4) # for dorsal surface
  angles$CrvGRFBStressAP <- angles$BMomGRFaxAP*(Anat$yAP_m/Anat$IAP_m4) # for anterior surface
  # Based on IG/gator data, GRF BMom at peak stress is mostly AP,
  # but still slight DV component in same direction as Quads;
  # therefore, set minimum as though Bmoms were perpendicular,
  # max as though they summed completely
  angles$NetBendFemMin <- ((angles$NetGRFBStressAP)^2 + (angles$MuscStressDV + angles$CrvGRFBStressDV)^2)^0.5
  angles$NetBendFemMax <- (((angles$MuscStressDV + angles$NetGRFBStressDV)^2) + (angles$CrvGRFBStressAP^2))^0.5
  angles$NetBendMinTens <- angles$NetBendFemMin + angles$StressAxFem
  angles$NetBendMinComp <- (-1*angles$NetBendFemMin) + angles$StressAxFem
  angles$NetBendMaxTens <- angles$NetBendFemMax + angles$StressAxFem
  angles$NetBendMaxComp <- (-1*angles$NetBendFemMax) + angles$StressAxFem
  # angle from DV, opposite willizrf
  angles$NetBendFemAngMin <- DEG2RAD*atan(angles$NetGRFBStressAP/(angles$MuscStressDV + angles$CrvGRFBStressDV))
  angles$NetBendFemAngMax <- DEG2RAD*atan(angles$CrvGRFBStressAP/(angles$MuscStressDV + angles$NetGRFBStressDV))
  
  angles$species <- species
  angles$mass_kg <- Anat$mass_kg
  angles$BW <- BW
  
  
  ## Plotting RGRFAnkle versus alpha
  plot(alpha*DEG2RAD, angles$RGRFAnkle, type = "l", col = "magenta", xlab = "alpha (degrees)", ylab = "RGRFAnkle (m)")
  
  # Need to fix these plots so min and max ranges of axes are large enough to include all points
  
  ## plotting Bending moments vs. alpha
  plot(alpha*DEG2RAD, angles$BMomGRFtvDV, col = "green")
  points(alpha*DEG2RAD, angles$QuadsBMom, col = "magenta")
  points(alpha*DEG2RAD, angles$BMomGRFaxDV, col = "blue")
  points(alpha*DEG2RAD, angles$BMomGRFaxAP, col = "yellow")
  
  ## Plotting Stresses vs. alpha
  plot(alpha*DEG2RAD, angles$StressAxFem, col = "red")   
  points(alpha*DEG2RAD, angles$MuscStressDV, col = "magenta")  
  points(alpha*DEG2RAD, angles$NetGRFBStressDV, col = "green")
  points(alpha*DEG2RAD, angles$NetGRFBStressAP, col = "blue")
  points(alpha*DEG2RAD, angles$CrvGRFBStressDV, col = "black", pch = 1)
  points(alpha*DEG2RAD, angles$CrvGRFBStressAP,col = "black", pch = 4)
  points(alpha*DEG2RAD, angles$NetBendFemMin, col = "yellow", pch = 8)
  points(alpha*DEG2RAD, angles$NetBendFemMax, col = "cyan", pch = 3)
  
  ## Plotting Net bending stresses (Pa) in compression vs. tension
  plot(alpha*DEG2RAD, angles$NetBendMinTens, col = "red")   
  points(alpha*DEG2RAD, angles$NetBendMinComp, col = "magenta")  
  points(alpha*DEG2RAD, angles$NetBendMaxTens, col = "green")
  points(alpha*DEG2RAD, angles$NetBendMaxComp, col = "blue")
  
  FMX_variables <- c("species", "mass_kg", "BW", "alpha", "GRFaxFem", "GRFtvFem", "FQuads", "FAext", "BMomGRFtvDV", "BMomGRFaxDV", "BMomGRFaxAP", "QuadsBMom")
  STX_variables <- c("species", "mass_kg", "BW", "alpha", "StressAxFem", "MuscStressDV", "CrvGRFBStressDV", "CrvGRFBStressAP", "NetGRFBStressDV", "NetGRFBStressAP", "NetBendFemMin", "NetBendFemMax", "NetBendFemAngMin", "NetBendFemAngMax", "NetBendMinTens", "NetBendMinComp", "NetBendMaxTens", "NetBendMaxComp")
  
  Anat$species <- species 
  
  output <- list(
    anatomy = Anat,
    angles = angles,
    FMX = angles[,FMX_variables],
    STX = angles[,STX_variables]
  )
  
  return(output)
}



############# smootheR ##########


#' @title Smooth numerical data using a polynomial smoothing spline
#'
#' @name smootheR
#'
#' @description This is a basic procedure to smooth numerica data using splines, which is partly useful for kinematic analyses.
#'
#' @usage smootheR(df, method = 1, norder = 3, spar = NULL, ...)
#'
#' @param \code{df} Numerical values that represent points to smooth. This could either be a vector or a data.frame.
#' @param \code{method} An integer numver that codes for the method is based on the different options available to smooth the data, following the description in the signal::smooth.Pspline() function. 
#' @param \code{norder} An integer number that codes the order for the polynomial used for the spline, following the description in the signal::smooth.Pspline() function.
#' @param \code{spar} An optional vector that can be provided to use custom parameters for smoothing the data, following the description in the signal::smooth.Pspline() function.
#' @details These procedures follow the methodology described in Kawano et al. 2016. 
#' @references Kawano SM, Economy DR, Kennedy MS, Dean D, Blob RW. 2016. Comparative limb bone loading in the humerus and femur of the tiger salamander Ambystoma tigrinum: testing the "mixed-chain" hypothesis for skeletal safety factors. Journal of Experimental Biology 219: 341-353. \url{http://jeb.biologists.org/content/219/3/341}
#'
#' @examples
#' 
#' 
#'
#' @export
#' @importFrom pspline smooth.Pspline
#' 

smootheR <- function(df, method = 1, norder = 3, spar = NULL, ...) {
  
  if (length(dim(df)) ==2) {
    smooth <- apply(df, 2, FUN = function(x) smooth.Pspline(1:nrow(df), x, norder= norder, method = method, spar = spar)$ysmth)
  } else {
    smooth <- smooth.Pspline(1:length(df), df, norder= norder, method = method, spar = spar)$ysmth
  }
  return(smooth)
}



############# interpolateR ##########


#' @title Interpolate numerical data
#'
#' @name interpolateR
#'
#' @description This is a basic procedure to interpolate numerica data into a different number of points, which is partly useful to standardize data in kinematic analyses.
#'
#' @usage interpolateR(df, n, method = 'spline', ...) 
#'
#' @param \code{df} Numerical values that represent points to interpolate. This could either be a vector or a data.frame.
#' @param \code{n} An integer numver for the number of points that you want to divide your data into. 
#' @param \code{method} A character string to indicate the interpolation method, following the description in the signal::interp1() function.The 'spline' is chosen as a default. 
#' @details These procedures follow the methodology described in Kawano et al. 2016. 
#' @references Kawano SM, Economy DR, Kennedy MS, Dean D, Blob RW. 2016. Comparative limb bone loading in the humerus and femur of the tiger salamander Ambystoma tigrinum: testing the "mixed-chain" hypothesis for skeletal safety factors. Journal of Experimental Biology 219: 341-353. \url{http://jeb.biologists.org/content/219/3/341}
#'
#' @examples
#' 
#' x <- seq(1, 10, 1)
#' y <- rnorm(10, 5, 0.1)
#' xy <- data.frame(x,y)
#' interpolateR(xy, 20, method = "cubic")
#' 
#'
#' @export
#' @importFrom signal interp1
#' 


interpolateR <- function(df, n, method = 'spline', ...) {
  
  # Determining 101 equally spaced points between the total number of frames / rows
  N <- nrow(df)/n
  
  # determining the points to interpolated depending on whether it's a vector vs. data.frame / matrix
  if (length(dim(df))==2) {
    X <- seq(1, nrow(df), length.out = n)
    apply(df, 2, FUN=function(x) interp1(1:nrow(df), x, X, method = method))
  } else {
    X <- seq(1, length(df), length.out = n)
    interp1(1:length(df), df, X, method = method)
  }
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