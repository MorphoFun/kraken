######### Kinetic data at the timing of the peak net GRF for fish and salamanders ############

#' @title GRFs of fins and limbs
#'
#' @name FinLimbGRFs
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
#' @usage FinLimbGRFs
#' @keywords datasets
#' @export devtools
#'
#' @examples 
#' data(FinLimbGRFs)

# setwd(paste(getwd(),"/dataraw",sep=""))
# FinLimbGRFs <- read.csv("PeakNetGRFData_150518.csv", sep=",")
# 
#  setwd("../data")
#  save(FinLimbGRFs, file="FinLimbGRFs.rda")



#' @title GRFs of fins and limbs
#'
#' @name FinLimb_Forces
#'
#' @description Ground reaction force (GRF) data from isolated appendages of mudskipper fish (\emph{Periophthalmus barbarus}), tiger salamanders (\emph{Ambystoma tigrinum}), Iberian ribbed newts (\emph{Pleurodeles waltl}, and greater siren (\emph{Siren lacertina})).
#'
#'  \describe{
#'			\item{\code{Light.Volts}}{Activity of the LED light used to synchronize the force and video data, in units of volts.}
#'			\item{\code{Vert1.Volts}}{First of four vertical channels, in units of volts.}
#'			\item{\code{Vert2.Volts}}{Second of four vertical channels, in units of volts.}
#'			\item{\code{Vert3.Volts}}{Third of four vertical channels, in units of volts.}
#'			\item{\code{Vert4.Volts}}{Fourth of four vertical channels, in units of volts.}
#'			\item{\code{VertSum.Volts}}{Sum of the four vertical channels, in units of volts.}
#'			\item{\code{ML1.Volts}}{One of two mediolateral channels, in units of volts.}
#'			\item{\code{ML2.Volts}}{Second of two mediolateral channels, in units of volts.}
#'			\item{\code{MLSum.Volts}}{Sum of the two mediolateral channels, in units of volts.}
#'			\item{\code{Hz1.Volts}}{One of two horizontal (anteroposterior) channels, in units of volts.}
#'			\item{\code{Hz2.Volts}}{Second of two horizontal (anteroposterior) channels, in units of volts.}
#'			\item{\code{HzSum.Volts}}{Sum of two horizontal (anteroposterior) channels, in units of volts.}
#'			\item{\code{EMPTY}}{Empty channel (no data collected), in units of volts.}
#'			
#'  }
#'
#' @references Kawano SM, Blob RW. 2013. Propulsive forces of mudskipper fins and salamander limbs during terrestrial locomotion: implications for the invasion of land.  \emph{Integrative and Comparative Biology} 53(2): 283-294.
#' @usage FinLimb_Forces
#' @keywords datasets
#' @export devtools
#'
#' @examples 
#' data(FinLimb_Forces)

# setwd(paste(getwd(),"/dataraw/Force_LabView_Output",sep=""))
# FinLimb_Forces <- lapply(list.files(), FUN = read.table, header = FALSE)
# for (i in 1:length(FinLimb_Forces)) {
#   names(FinLimb_Forces[[i]]) <- c("Light.Volts", "Vert1.Volts", "Vert2.Volts", "Vert3.Volts", "Vert4.Volts", "VertSum.Volts", "ML1.Volts", "ML2.Volts", "MLSum.Volts", "Hz1.Volts", "Hz2.Volts", "HzSum.Volts", "EMPTY", "EMPTY", "EMPTY", "EMPTY", "EMPTY")
# }
# names(FinLimb_Forces) <- substring(list.files(), nchar(list.files())-10, nchar(list.files())-4)
# 
# setwd("../../data")
# save(FinLimb_Forces, file="FinLimb_Forces.rda")
