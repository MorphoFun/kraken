% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BiomechanicsFunctions.R
\name{jointAngle}
\alias{jointAngle}
\title{Calculate joint angles (in degrees) from XYZ coordinate data}
\usage{
jointAngle(P1, P2, P3)
}
\arguments{
\item{\code{P1}}{A data.frame of numeric values containing the X, Y, and Z coordinate data, respectively, for the first point (e.g, the shoulder)}

\item{\code{P2}}{A data.frame of numeric values containing the X, Y, and Z coordinate data, respectively, for the second point, which is assumed to be the vertex of the angle (e.g, the elbow)}

\item{\code{P3}}{A data.frame of numeric values containing the X, Y, and Z coordinate data, respectively, for the third point (e.g, the wrist)}
}
\description{
Calculates the angle of a joint (in degrees), formed by three points with XYZ coordinates.
}
\details{
These procedures follow the methodology used in Kawano and Blob (2013) and Kawano et al. 2016 to calculate angles formed about the limb joints in animals.
}
\examples{

shoulder <- matrix(c(0.006305306, 0.006526961, 0.006747555, -0.08206114, -0.08207707, -0.08207049, 0.006997669, 0.006980824, 0.006975157), 3, 3)
elbow <- matrix(c(0.007826633, 0.007959096, 0.008068560, -0.07183020, -0.07185459, -0.07186337, 0.005754819, 0.005764666, 0.005774707), 3, 3)
wrist <- matrix(c(0.01164592, 0.01160690, 0.01157642, -0.07348876, -0.07345559, -0.07342105, -0.000631402, -0.000671288, -0.000709513), 3, 3)

elbowAngle <- jointAngle(P1 = shoulder, P2 = elbow, P3 = wrist)

}
\references{
Kawano SM, Blob RW. 2013. Propulsive forces of mudskipper fins and salamander limbs during terrestrial locomotion: implications for the invasion of land. Integrative and Comparative Biology 53(2): 283-294. \url{https://academic.oup.com/icb/article/53/2/283/806410/Propulsive-Forces-of-Mudskipper-Fins-and}

Kawano SM, Economy DR, Kennedy MS, Dean D, Blob RW. 2016. Comparative limb bone loading in the humerus and femur of the tiger salamander Ambystoma tigrinum: testing the "mixed-chain" hypothesis for skeletal safety factors. Journal of Experimental Biology 219: 341-353. \url{http://jeb.biologists.org/content/219/3/341}

"R - Comute Cross Product of Vectors (Physics)" answer posted by user Kevin on April 22, 2016. \url{https://stackoverflow.com/questions/36798301/r-compute-cross-product-of-vectors-physics}

"how to calculate the Euclidean norm of a vector in R?" answer posted by user joran on June 7, 2012. \url{https://stackoverflow.com/questions/10933945/how-to-calculate-the-euclidean-norm-of-a-vector-in-r}
}
