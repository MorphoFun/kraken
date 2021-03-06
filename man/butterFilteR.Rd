% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BiomechanicsFunctions.R
\name{butterFilteR}
\alias{butterFilteR}
\title{Apply a butterworth filter to data}
\usage{

}
\arguments{
\item{\code{df}}{A list containing the file name and an array of data containing the independent variable (time) as the first column that is followed by the dependent variables. Can take voltToForce output as an input.}

\item{\code{Fs}}{A numeric value indicating the sampling frequency. 5000 Hz is set as a default.}

\item{\code{PbF}}{A numeric value indicating the pass-band frequency. 6 is set as a default.}

\item{\code{SbF}}{A numeric value indicating the stop-band frequency. 190 is set as a default.}

\item{\code{Rp}}{A numeric value indicating passband ripple in dB; represents the max permissible passband loss. 2 dB is set as a default.}

\item{\code{Rs}}{A numeric value indicating stopband attenuation in dB; respresents the dB the stopband is down from the passband. 40 dB is set as a default.}

\item{\code{saveAs}}{A character string containing the name for the resulting graph.}

\item{\code{saveGraph}}{A character string (options = "yes" or "no") to indicate whether you want to save the graph as a PDF or not.}
}
\description{
Applies a butterworth filter to data, using information about the data to determine what polynomial to use.
}
\details{
These procedures follow the methodology used in Kawano and Blob (2013) and Kawano et al. 2016. It is assumed that the output from the force platform contain 12 channels in the following order: trigger, four verticals, sum of the verticals, two mediolateral, sum of the mediolaterals, two anteroposterior, and the sum of the anteroposteriors.
}
\examples{

GRF <- voltToForce(af01f18, calib = c(-0.7126, 1.3243, 1.3476), lightStartFrame = 276, startFrame = 35, endFrame = 148, zeroStart = 35, filename = "af01f18")
GRF_filtered <- butterFilteR(GRF, saveAs = "af01f18_Pec_Filter.pdf", saveGraph = "yes")

}
\references{
Kawano SM, Blob RW. 2013. Propulsive forces of mudskipper fins and salamander limbs during terrestrial locomotion: implications for the invasion of land. Integrative and Comparative Biology 53(2): 283-294. \url{https://academic.oup.com/icb/article/53/2/283/806410/Propulsive-Forces-of-Mudskipper-Fins-and}

Kawano SM, Economy DR, Kennedy MS, Dean D, Blob RW. 2016. Comparative limb bone loading in the humerus and femur of the tiger salamander Ambystoma tigrinum: testing the "mixed-chain" hypothesis for skeletal safety factors. Journal of Experimental Biology 219: 341-353. \url{http://jeb.biologists.org/content/219/3/341}
}
