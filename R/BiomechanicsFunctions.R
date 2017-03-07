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
#' @references Sacks RD, Roy RR. 1982. Architecture of The Hind Limb Muscles of Cats: Functional Significance. Journal of Morphology, 185–195. \url{http://onlinelibrary.wiley.com/doi/10.1002/jmor.1051730206/abstract}
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
  ifelse(ncol(cbind(time,GRF) == 2), totalImpulse <- sum(diff(time)*rollmean(GRF,2)), totalImpulse <- sapply(GRF, FUN = function(x) sum(diff(time)*rollmean(x,2))))
  ifelse(ncol(cbind(time,GRF) == 2), rollImpulse <- diff(time)*rollmean(GRF,2), rollImpulse <- sapply(GRF, FUN = function(x) sum(diff(time)*rollmean(x,2))))
    output <- list(totalImpulse = totalImpulse, rollImpulse = rollImpulse)
  return(output)
}

