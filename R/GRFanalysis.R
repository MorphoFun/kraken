######### Analysis of FinLimb GRF data ##########
## All data were taken at the peak net GRF for each individual trial

library(rstan)
library(devtools)
library(psa)

# loading from remote folder
url <- ("https://raw.githubusercontent.com/MorphoFun/kraken/master/dataraw/PeakNetGRFData_150518.csv")
FinLimbGRFs <- read.csv(file = url)
  
  
########## GRF DATA ############

# stan can't handle factor types (e.g., Species in character format), so need to convert factor levels to dummy variables using as.integer()
FinLimbGRFs$SpeciesInt <- as.integer(FinLimbGRFs$Species)
d <- FinLimbGRFs
# converting percentages as proportions to make them closer in magnitude to the .BW data
d[,1] <- FinLimbGRFs[,1]/100 
# converting degrees to radians to make them closer in magnitude to the .BW data
d[,2:3] <- FinLimbGRFs[,2:3]*(pi/180)


# decided against standardizing data b/c
# 1. I want to be able to interpret the results in their original units
# 2. I standardized the data to units of BW (body weight), and then converted the PercentStance and
#    angle data (from degrees to radians), so that the variables of interest are within the range of 
#    about -1 to 1. Not much of a scaling difference, so maybe not necessary to scale() the data.


######## TESTING FOR OMNIBUS ##########
library(fBasics)
# omnibus null = data are normal
fullDago <- sapply(d[,1:7], function(x) dagoTest(x))
pecDago <- sapply(subset(d, Appendage=="Pec", select = c(1:7)), function(x) dagoTest(x))
HLDago <- sapply(subset(d, Appendage=="Pel", select = c(1:7)), function(x) dagoTest(x))
FLDago <- sapply(subset(d, Appendage=="Pec" & !Species=="pb", select = c(1:7)), function(x) dagoTest(x))


######## PLOTTING HISTOGRAMS OF THE DATA #######
library(ggplot2) 
 
ggplot(d) + geom_histogram(aes(d$PercentStance, fill = interaction(d$Species, d$Appendage)), binwidth = 0.05) + facet_grid(Species~Appendage) 
ggplot(d) + geom_histogram(aes(d$APAngle.Convert, fill = interaction(d$Species, d$Appendage)), binwidth = 0.05) + facet_grid(Species~Appendage) 
ggplot(d) + geom_histogram(aes(d$MLAngle.Convert, fill = interaction(d$Species, d$Appendage)), binwidth = 0.05) + facet_grid(Species~Appendage) 
ggplot(d) + geom_histogram(aes(d$InterpV.BW, fill = interaction(d$Species, d$Appendage)), binwidth = 0.05) + facet_grid(Species~Appendage) 
ggplot(d) + geom_histogram(aes(d$InterpML.BW, fill = interaction(d$Species, d$Appendage)), binwidth = 0.05) + facet_grid(Species~Appendage) 
ggplot(d) + geom_histogram(aes(d$InterpAP.BW, fill = interaction(d$Species, d$Appendage)), binwidth = 0.05) + facet_grid(Species~Appendage) 
ggplot(d) + geom_histogram(aes(d$NetGRF.BW, fill = interaction(d$Species, d$Appendage)), binwidth = 0.05) + facet_grid(Species~Appendage) 


####### TESTING COEFFICIENT OF VARIATION TO IND VARIABILITY #####
## Looking to pool the data for each individual, but first need to determine how variable the data are 
## for each individual. 

indSummaryData <- aggregate(d[,1:7], by = list(d$Ind), function(x) {c(MEAN = mean(x), SD = sd(x), CoV = sd(x)/mean(x))})
indCoV <- aggregate(d[,1:7], by = list(d$Ind), function(x) {c(CoV = (sd(x)/mean(x))*100)})

## Ugh, a number of the variables seem to have a lot of variable within an individual.... 


######### REDUNDANCY ANALYSIS ###############
# Example: https://rstudio-pubs-static.s3.amazonaws.com/64619_2f93b223a318410bbf999d092ecf05ec.html
# For these data, the "sites" in vegan are actually individual trials,
# and the "Species" are the variables that are being compared amongst the species 

library(vegan)

###### RDA: pectoral appendages of the three species ####
pecrda <- rda(subset(d, Appendage=="Pec")[,1:7] ~ Species, data = subset(d, Appendage=="Pec"))
RsquareAdj(pecrda) # ~ 0.35

# pecindrda <- rda(subset(d[,1:7], d$Appendage=="Pec") ~ Species + Ind, data = subset(d, d$Appendage=="Pec"))
# RsquareAdj(pecindrda) # ~ 0.48

# eigenvalues for both RD1 and RD2 are <1

colvec <- c("brown", "blue", "lightblue")
pecscors <- scores(pecrda, display = c("sites", "species"), scaling = 3)
pecxlimit <- with(pecscors, range(species[,1], sites[,1]))
pecylimit <- with(pecscors, range(species[,2], sites[,2]))

plot(pecrda, xlab = paste("RD1 (", summary(pecrda)$cont[[1]][2]*100, "%)"), 
     ylab = paste("RD2 (", summary(pecrda)$cont[[1]][5]*100, "%)"),
     type = "n",
     scaling = 3, 
     xlim = c(pecxlimit[1]-0.5, pecxlimit[2]+0.5),
     ylim = c(pecylimit[1]-0.5, pecylimit[2]+0.5)
)
with(subset(d, Appendage=="Pec"), points(pecrda, display = "sites", col = colvec[Species], pch = 19),
     scaling = 3)
with(subset(d, Appendage=="Pec"), legend("topleft", legend = levels(Species), bty = "n", col = colvec, pch = 19))
ellipses <- ordiellipse(pecrda, subset(d, Appendage=="Pec")$Species, kind = "sd", conf = 0.95, lwd=2, draw = "polygon", col = "brown", alpha = 80, show.groups = "at") 
ellipses <- ordiellipse(pecrda, subset(d, Appendage=="Pec")$Species, kind = "sd", conf = 0.95, lwd=2, draw = "polygon", col = "lightblue", alpha = 80, show.groups = "pw") 
ellipses <- ordiellipse(pecrda, subset(d, Appendage=="Pec")$Species, kind = "sd", conf = 0.95, lwd=2, draw = "polygon", col = "blue", alpha = 80, show.groups = "pb") 


####### RDA: pelvic appendages of the two tetrapod species ####
pelrda <- rda(subset(d, Appendage=="Pel")[,1:7] ~ Species, data = subset(d, Appendage=="Pel"))
RsquareAdj(pelrda) # ~ 0.04

# pelindrda <- rda(subset(d[,1:7], d$Appendage=="Pel") ~ Species + Ind, data = subset(d, d$Appendage=="Pel"))
# RsquareAdj(pelindrda) # ~ 0.23

# eigenvalues for RD <1
# no bivariate plot since there is only one axis
boxplot(summary(pelrda)$sites ~ subset(d, Appendage=="Pel")$Species, xlab = "species", ylab = paste("RD1 (", summary(pelrda)$cont[[1]][2]*100, "%)"))

####### RDA: 5 species:appendage combos (stored under "Group") ####
grouprda <- rda(d[,1:7] ~ Group, data = d)
RsquareAdj(grouprda) # ~ 0.53

# groupindrda <- rda(d[,1:7] ~ Group + Ind, data = d)
# RsquareAdj(groupindrda) # ~ 0.58

# eigenvalues for all RD axes are <1

colvec <- c("brown", "blue", "lightblue")
groupscors <- scores(grouprda, display = c("sites", "species"), scaling = 3)
groupxlimit <- with(groupscors, range(species[,1], sites[,1]))
groupylimit <- with(groupscors, range(species[,2], sites[,2]))

colgroup <- c("tan", "brown", "blue", "lightblue", "cyan")
plot(grouprda, xlab = paste("RD1 (", summary(grouprda)$cont[[1]][2]*100, "%)"), 
     ylab = paste("RD2 (", summary(grouprda)$cont[[1]][5]*100, "%)"),
     type = "n",
     scaling = 3, 
     xlim = c(groupxlimit[1], groupxlimit[2]),
     ylim = c(groupylimit[1], groupylimit[2])
)
with(d, points(grouprda, display = "sites", col = colgroup[Group], pch = 19),
     scaling = 3)
with(d, legend("topleft", legend = levels(Group), bty = "n", col = colgroup, pch = 19))
ellipses <- ordiellipse(grouprda, d$Group, kind = "sd", conf = 0.95, lwd=2, draw = "polygon", col = "tan", alpha = 80, show.groups = "at_Pec") 
ellipses <- ordiellipse(grouprda, d$Group, kind = "sd", conf = 0.95, lwd=2, draw = "polygon", col = "brown", alpha = 80, show.groups = "at_Pel") 
ellipses <- ordiellipse(grouprda, d$Group, kind = "sd", conf = 0.95, lwd=2, draw = "polygon", col = "blue", alpha = 80, show.groups = "pb_Pec")
ellipses <- ordiellipse(grouprda, d$Group, kind = "sd", conf = 0.95, lwd=2, draw = "polygon", col = "lightblue", alpha = 80, show.groups = "pw_Pec") 
ellipses <- ordiellipse(grouprda, d$Group, kind = "sd", conf = 0.95, lwd=2, draw = "polygon", col = "cyan", alpha = 80, show.groups = "pw_Pel") 

####### RDA: Pooled data across the three species (just for shits and giggles) ####
speciesrda <- rda(d[,1:7] ~ Species, data = d)
# constrained variable only explains about 10% of the variance? yikes...
RsquareAdj(speciesrda) # pretty weak

# speciesindrda <- rda(d[,1:7] ~ Species + Ind, data = d)
# # Species with inds explains about 20% of the variance
# RsquareAdj(speciesindrda)

# eigenvalues for both RD1 and RD2 are <1, so the rda tells me that species classification doesn't explain a lot of the variation in the GRF variables?

colvec <- c("brown", "blue", "lightblue")
scors <- scores(speciesrda, display = c("sites", "species"), scaling = 3)
xlimit <- with(scors, range(species[,1], sites[,1]))
ylimit <- with(scors, range(species[,2], sites[,2]))

plot(speciesrda, xlab = paste("RD1 (", summary(speciesrda)$cont[[1]][2]*100, "%)"), 
     ylab = paste("RD2 (", summary(speciesrda)$cont[[1]][5]*100, "%)"),
     type = "n",
     scaling = 3, 
     xlim = c(xlimit[1]-0.5, xlimit[2]+0.5),
     ylim = c(ylimit[1]-0.5, ylimit[2]+0.5)
)
with(d, points(speciesrda, display = "sites", col = colvec[Species], pch = 19),
     scaling = 3)
with(d, legend("topleft", legend = levels(Species), bty = "n", col = colvec, pch = 19))
ellipses <- ordiellipse(speciesrda, d$Species, kind = "sd", conf = 0.95, lwd=2, draw = "polygon", col = "brown", alpha = 80, show.groups = "at") #ellipses with se are way small
ellipses <- ordiellipse(speciesrda, d$Species, kind = "sd", conf = 0.95, lwd=2, draw = "polygon", col = "lightblue", alpha = 80, show.groups = "pw") #ellipses with se are way small
ellipses <- ordiellipse(speciesrda, d$Species, kind = "sd", conf = 0.95, lwd=2, draw = "polygon", col = "blue", alpha = 80, show.groups = "pb") #ellipses with se are way small

### Things to try:
## Try to calculate area of overlap between the ellipses? 
## Calculate significance of the axes? http://onlinelibrary.wiley.com/doi/10.1111/j.2041-210X.2010.00078.x/pdf


# defining response for the stan model later
response <- d$NetGRF.BW


# ########## SELECTING VARIABLES FOR COMPARISON ##############
# species_dat <- list(N = nrow(FinLimbGRFs),
#                     S = as.integer(FinLimbGRFs$Species),
#                     I = as.integer(FinLimbGRFs$Ind),
#                     y = as.numeric(unlist(aggregate(response, list(predictor), mean)[2])),
#                     sigma = as.numeric(unlist(aggregate(response, list(predictor), sd)[2]))
# )
# 
# # For comparison of the pooled data for each species
# species_dat <- list(S = length(levels(FinLimbGRFs$Species)), 
#                     y = as.numeric(unlist(aggregate(response, list(predictor), mean)[2])),
#                     sigma = as.numeric(unlist(aggregate(response, list(predictor), sd)[2]))
# )
# 
# dat <-list(spp = as.integer(d$Species),
#             app = as.integer(d$Appendage),
#             ind = as.integer(d$Ind),
#             y = response,
#             N = nrow(d),
#             J = length(unique(d$Ind))
#           )
# 
# # see page 11 of the Rstan: the R interface to Stan manual; do I need to convert y to an array? S should never equal 1, so I'm guessing the answer is "no". 
# 
# ########## STAN MODEL FOR COMPARING SPECIES MEANS ###########
# # Should the data and parameters be defined as vectors instead? https://groups.google.com/forum/#!msg/stan-users/4PgOF38Mnwk/hgUPCA768w0J
# # Should vector types always be used for linear (algebra) problems?
# 
# stanMod <- '
# data {
#   int<lower=1> N;                  //number of data points
#   real y[N];                       //response
#   int<lower=1,upper=3> spp [N];    //predictor
#   int<lower=1> J;                  //number of individuals
#   int<lower=1, upper=J> ind[N];    //ind id
# }
# 
# parameters {
#   vector[2] beta;            //fixed intercept and slope
#   vector[J] u;               //ind intercepts
#   real<lower=0> sigma_e;     //error sd
#   real<lower=0> sigma_u;     //ind sd
# }
# 
# model {
#   real mu;
#   //priors
#   u ~ normal(0,sigma_u);    //ind random effects
#   // likelihood
#   for (i in 1:N){
#   mu <- beta[1] + u[ind[i]] + beta[2]*spp[i];
#   y[i] ~ lognormal(mu,sigma_e);
#   }
# }'
# 
# ## Sample from posterior distribution:
# stanFit <- stan(model_code = stanMod, data=dat, iter=2000, chains=4)
# print(stanFit)
# 
# ## Summarize results:
# print(stanFit,pars=c("beta","sigma_e","sigma_u"),probs=c(0.025,0.5,0.975))
# 
# beta1 <- extract(stanFit,pars=c("beta[2]"))$beta
# print(signif(quantile(beta1,probs=c(0.025,0.5,0.975)),2))
# 
# ## Posterior probability of beta1 being less than 0:
# mean(beta1<0)
