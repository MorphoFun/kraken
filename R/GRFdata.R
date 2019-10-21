######### Kinetic data for fish and salamanders ############

#' @title Peak net GRFs of fins and limbs
#'
#' @name FinLimbGRFs_Peak
#'
#' @description Ground reaction force (GRF) data from isolated appendages of mudskipper fish (\emph{Periophthalmus barbarus}), tiger salamanders (\emph{Ambystoma tigrinum}), and Iberian ribbed newts (\emph{Pleurodeles waltl}).
#'
#'  \describe{
#'			\item{\code{PercentStance}}{The percentage of the stance phase when the peak net GRF occurs for that trial.}
#'      \item{\code{APAngle.Convert}}{Angle of the GRF in the anteroposterior direction, in the -90 (posterior) to +90 (posterior) degree range so that 0 degrees represent the transition between anterorior and posterior.}
#'      \item{\code{MLAngle.Convert}}{Angle of the GRF in the mediolateral direction, in the -90 (medial) to +90 (lateral) degree range so that 0 degrees represent the transition between anterorior and posterior.}
#'      \item{\code{InterpV.BW}}{Magnitude of the vertical component of the GRF, in units of body weight.}
#'      \item{\code{InterpML.BW}}{Magnitude of the mediolateral component of the GRF, in units of body weight.}
#'      \item{\code{InterpAP.BW}}{Magnitude of the anteroposterior component of the GRF, in units of body weight.}
#'      \item{\code{NetGRF.BW}}{Magnitude of the net GRF (summation of the vertical, mediolateral, and anteroposterior components), in units of body weight.}
#'      \item{\code{Group}}{The group for comparison. Species codes include: "at" (\emph{Ambystoma tigrinum}), "pb" (\emph{Periophthalmus barbarus}), and "pw" (\emph{Pleurodeles waltl}). Appendage codes include: "Pec" (Pectoral), and "Pel" (Pelvic).}
#'      \item{\code{Species}}{Two-letter abbreviation for the genus and species epithet. Species codes include: "at" (\emph{Ambystoma tigrinum}), "pb" (\emph{Periophthalmus barbarus}), and "pw" (\emph{Pleurodeles waltl}).}
#'      \item{\code{Ind}}{Indentification code for individual animal, with the first two letters as the Species.}
#'      \item{\code{Appendage}}{Appendage codes include: "Pec" (Pectoral), and "Pel" (Pelvic).}
#'      \item{\code{Trial}}{Specific identification code for the Trial. Note: data from for the pectoral and pelvic appendages of an individual can be obtained from the same trial.}
#'      \item{\code{InterpV.N}}{Magnitude of the vertical component of the GRF, in units of Newtons.}
#'      \item{\code{InterpML.N}}{Magnitude of the mediolateral component of the GRF, in units of Newtons.}
#'      \item{\code{InterpAP.N}}{Magnitude of the anteroposterior component of the GRF, in units of Newtons.}
#'      \item{\code{SquaredV.N}}{Squared value of the vertical component of the GRF, in units of Newtons.}
#'      \item{\code{SquaredML.N}}{Squared value of the mediolateral component of the GRF, in units of Newtons.}
#'      \item{\code{SquaredAP.N}}{Squared value of the anteroposterior component of the GRF, in units of Newtons.}
#'      \item{\code{NetGRF.N}}{Magnitude of the Net GRF, in units of Newtons.}
#'      \item{\code{APAngle}}{Angle of the GRF in the anteroposterior direction, in the 0-180 degree range}
#'      \item{\code{MLAngle}}{Angle of the GRF in the mediolateral direction, in the 0-180 degree range}
#'  }
#'
#' @references Kawano SM, Blob RW. 2013. Propulsive forces of mudskipper fins and salamander limbs during terrestrial locomotion: implications for the invasion of land.  \emph{Integrative and Comparative Biology} 53(2): 283-294.
#' @usage FinLimbGRFs_Peak
#' @keywords datasets
#' @export FinLimbGRFs_Peak
#'
#' @examples 
#' data(FinLimbGRFs_Peak)

FinLimbGRFs_Peak <- read.csv("./dataraw/PeakNetGRFData.csv", sep=",")

setwd("./data")
save(FinLimbGRFs_Peak, file="FinLimbGRFs_Peak.rda")


#' @title Example of raw data file for salaamander ground reaction force analysis 
#'
#' @name af01f18
#'
#' @description Example LabView output collected while a tiger salamander (\emph{Ambystoma tigrinum}) was walking on a force plate. The name of this trial is af01f18. Includes data from the forelimb and hind limb. 
#'
#'  \describe{
#'			\item{\code{PercentStance}}{The percentage of the stance phase when the peak net GRF occurs for that trial.}
#'      \item{\code{APAngle.Convert}}{Angle of the GRF in the anteroposterior direction, in the -90 (posterior) to +90 (posterior) degree range so that 0 degrees represent the transition between anterorior and posterior.}
#'      \item{\code{MLAngle.Convert}}{Angle of the GRF in the mediolateral direction, in the -90 (medial) to +90 (lateral) degree range so that 0 degrees represent the transition between anterorior and posterior.}
#'      \item{\code{InterpV.BW}}{Magnitude of the vertical component of the GRF, in units of body weight.}
#'      \item{\code{InterpML.BW}}{Magnitude of the mediolateral component of the GRF, in units of body weight.}
#'      \item{\code{InterpAP.BW}}{Magnitude of the anteroposterior component of the GRF, in units of body weight.}
#'      \item{\code{NetGRF.BW}}{Magnitude of the net GRF (summation of the vertical, mediolateral, and anteroposterior components), in units of body weight.}
#'      \item{\code{Group}}{The group for comparison. Species codes include: "at" (\emph{Ambystoma tigrinum}), "pb" (\emph{Periophthalmus barbarus}), and "pw" (\emph{Pleurodeles waltl}). Appendage codes include: "Pec" (Pectoral), and "Pel" (Pelvic).}
#'      \item{\code{Species}}{Two-letter abbreviation for the genus and species epithet. Species codes include: "at" (\emph{Ambystoma tigrinum}), "pb" (\emph{Periophthalmus barbarus}), and "pw" (\emph{Pleurodeles waltl}).}
#'      \item{\code{Ind}}{Indentification code for individual animal, with the first two letters as the Species.}
#'      \item{\code{Appendage}}{Appendage codes include: "Pec" (Pectoral), and "Pel" (Pelvic).}
#'      \item{\code{Trial}}{Specific identification code for the Trial. Note: data from for the pectoral and pelvic appendages of an individual can be obtained from the same trial.}
#'      \item{\code{InterpV.N}}{Magnitude of the vertical component of the GRF, in units of Newtons.}
#'      \item{\code{InterpML.N}}{Magnitude of the mediolateral component of the GRF, in units of Newtons.}
#'      \item{\code{InterpAP.N}}{Magnitude of the anteroposterior component of the GRF, in units of Newtons.}
#'      \item{\code{SquaredV.N}}{Squared value of the vertical component of the GRF, in units of Newtons.}
#'      \item{\code{SquaredML.N}}{Squared value of the mediolateral component of the GRF, in units of Newtons.}
#'      \item{\code{SquaredAP.N}}{Squared value of the anteroposterior component of the GRF, in units of Newtons.}
#'      \item{\code{NetGRF.N}}{Magnitude of the Net GRF, in units of Newtons.}
#'      \item{\code{APAngle}}{Angle of the GRF in the anteroposterior direction, in the 0-180 degree range}
#'      \item{\code{MLAngle}}{Angle of the GRF in the mediolateral direction, in the 0-180 degree range}
#'  }
#'
#' @references Kawano SM, Blob RW. 2013. Propulsive forces of mudskipper fins and salamander limbs during terrestrial locomotion: implications for the invasion of land.  \emph{Integrative and Comparative Biology} 53(2): 283-294.
#' @usage af01f18
#' @keywords datasets
#' @export af01f18
#'
#' @examples 
#' data(af01f18)

af01f18 <- read.table("./dataraw/af01f18.txt")
names(af01f18) <- c("light_Volts", "Vert1.Volts", "Vert2.Volts", "Vert3.Volts", "Vert4.Volts", "VertSum.Volts", "ML1.Volts", "ML2.Volts", "MLSum.Volts", "Hz1.Volts", "Hz2.Volts", "HzSum.Volts", "BLANK", "BLANK", "BLANK", "BLANK", "BLANK")
af01f18$Sweep <- 1:nrow(af01f18)


setwd("./data")
save(af01f18, file="af01f18.rda")
