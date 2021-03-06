% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BiomechanicsFunctions.R
\name{pcsa}
\alias{pcsa}
\title{Physiological cross-sectional area of muscles (PCSA)}
\usage{
pcsa(mass, pennationAngle, fiberLength, density, stringAsFactors = FALSE, ...)
}
\arguments{
\item{\code{mass}}{Numerical value of the muscle mass, in units of kilograms.}

\item{\code{pennationAngle}}{Numerical value of the pennational angle of the muscle, in units of degrees.}

\item{\code{fasicleLength}}{Numerical value of the muscle fiber length, in units of meters. Often assumed to be average muscle fiber length.}

\item{\code{density}}{Numerical value of the muscle density. Defaults to 1060 kg/m^3, a typical value for striated muscles (Biewener 2003)}
}
\description{
\code{pcsa} allows one to estimate the physiological cross-sectional area of a muscle
}
\details{
See description in Hutchinson et al. (2015) and Sacks and Roy (1982) for more details regarding the calculation of PCSA.
}
\examples{

pcsa(0.3788, 0, 0.174)

}
\references{
Biewener AA. 2003. Animal locomotion. Oxford, UK: Oxford University Press.

Hutchinson JR, Rankin JW, Rubenson J, Rosenbluth KH, Siston RA, Delp SL. 2015. Musculoskeletal modelling of an ostrich (Struthio camelus) pelvic limb: influence of limb orientation on muscular capacity during locomotion. \url{http://dx.doi.org/10.7717/peerj.1001}

Sacks RD, Roy RR. 1982. Architecture of The Hind Limb Muscles of Cats: Functional Significance. Journal of Morphology, 185-195. \url{http://onlinelibrary.wiley.com/doi/10.1002/jmor.1051730206/abstract}
}
