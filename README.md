[![Lifecycle:Experimental](https://img.shields.io/badge/Lifecycle-Experimental-339999)](<Redirect-URL>)

# kraken

## Description
```kraken``` is a package with biomechanics function related to bone loading and fracture.

## Features
Some of the main features include:  
* Physiological cross-sectional area: pcsa()
* Impulse from GRF data: impulse()
* Digigal zero-phase butterworth filter with custom polynomial for GRF data that also interpolates to 101 points using 1st and 2nd derivatives: butterFilteR()
* Yaw from XY coordinate data: yaw()
* Protraction / Retraction from XYZ coordinate data: protraction()
* Pitch from XYZ coordinate data: pitch()
* Angles of GRF orientation: GRFAngles()
* Identifying regions of a video where two structures overlap and then replacing those observations with NA values: removeOverlaps()
* Rate of change in force production: yank()
* Estimation of peak bone stresses in the limb bones of fossil specimens: boneload_extinct()


## Vignettes
Stay tuned!

## Installation
```kraken``` is currently in a development version, and can be access via:

```
library(devtools)
install_github("MorphoFun/kraken", dependencies = TRUE)
```
