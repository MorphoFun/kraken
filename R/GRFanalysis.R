######### Analysis of FinLimb GRF data ##########
## All data were taken at the peak net GRF for each individual trial

library(rstan)
library(devtools)

readGitHub <- function(url, header = TRUE, fill = TRUE, stringAsFactors = FALSE, ...) {
  store <- tempfile()
  download.file(url, destfile=store, method="curl")
  dd <- read.table(store, header = header, ...)
  return(dd)
}


# loading from remote folder
FinLimbGRFs <- readGitHub("https://raw.githubusercontent.com/MorphoFun/kraken/master/dataraw/PeakNetGRFData_150518.csv", sep = ",")

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
 
ggplot(d) + geom_histogram(aes(d$PercentStance, fill = interaction(d$Species, d$Appendage))) + facet_grid(Species~Appendage, binwidth = 40) 
ggplot(d) + geom_histogram(aes(d$APAngle.Convert, fill = interaction(d$Species, d$Appendage))) + facet_grid(Species~Appendage) 
ggplot(d) + geom_histogram(aes(d$MLAngle.Convert, fill = interaction(d$Species, d$Appendage))) + facet_grid(Species~Appendage) 
ggplot(d) + geom_histogram(aes(d$InterpV.BW, fill = interaction(d$Species, d$Appendage))) + facet_grid(Species~Appendage) 
ggplot(d) + geom_histogram(aes(d$InterpML.BW, fill = interaction(d$Species, d$Appendage))) + facet_grid(Species~Appendage) 
ggplot(d) + geom_histogram(aes(d$InterpAP.BW, fill = interaction(d$Species, d$Appendage))) + facet_grid(Species~Appendage) 
ggplot(d) + geom_histogram(aes(d$NetGRF.BW, fill = interaction(d$Species, d$Appendage))) + facet_grid(Species~Appendage) 


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


# defining response vs. predictor for the stan model later
response <- d$NetGRF.BW
predictor <- d$SpeciesInt

########## SELECTING VARIABLES FOR COMPARISON ##############
species_dat <- list(N = nrow(FinLimbGRFs),
                    S = as.integer(FinLimbGRFs$Species),
                    I = as.integer(FinLimbGRFs$Ind),
                    y = as.numeric(unlist(aggregate(response, list(predictor), mean)[2])),
                    sigma = as.numeric(unlist(aggregate(response, list(predictor), sd)[2]))
)

# For comparison of the pooled data for each species
species_dat <- list(S = length(levels(FinLimbGRFs$Species)), 
                    y = as.numeric(unlist(aggregate(response, list(predictor), mean)[2])),
                    sigma = as.numeric(unlist(aggregate(response, list(predictor), sd)[2]))
)

# see page 11 of the Rstan: the R interface to Stan manual; do I need to convert y to an array? S should never equal 1, so I'm guessing the answer is "no". 

########## STAN MODEL FOR COMPARING SPECIES MEANS ###########
# Should the data and parameters be defined as vectors instead? https://groups.google.com/forum/#!msg/stan-users/4PgOF38Mnwk/hgUPCA768w0J
# Should vector types always be used for linear (algebra) problems?

### stan mcode for hierarchical model with multivariate priors on teh group-level coefficients and group-level
## prior means (p. 66 on stan manual):
hierstanmod <- '
data {
  int<lower=0> N;             // num individuals
  int<lower=1> K;             // num ind predictors
  int<lower=1> J;             // num groups
  int<lower=1> L;             // num group predictors
  int<lower=1,upper=J> jj[N]; // group for individual
  matrix[N,K] x;              // individual predictors
  row_vector[L] u[J];         // group predictors
  vector[N] y;                // outcomes
}
parameters {
  corr_matrix[K] Omega;       // prior correlation
  vector<lower=0>[K] tau;     // prior scale
  matrix[L,K] gamma;          // group coeffs
  vector[K] beta[J];          // ind coeffs by group
  real<lower=0> sigma;        // prediction error scale
}
model {
  tau ~ cauchy(0,2.5);
  Omega ~ lkj_corr(2);
  to_vector(gamma) ~ normal(0,
{
  row_vector[K] u_gamma[J];
  for (j in 1:J)
    u_gamma[j] <- u[j] * gamma;
    beta ~ multi_normal(u_gamma, quad_form_diag(Omega, tau));
}
{
  vector[N] x_beta_jj;
  for (n in 1:N)
    x_beta_jj[n] <- x[n] * beta[jj[n]];
  y ~ normal(x_beta_jj, sigma);
}
}'

mod <- '
data {
  int<lower=1> S;        // dummy coding for species
  int<lower=1> App;      // dummy coding for appendage type
  int<lower=1> N;        // total number of observations
  vector[N] y;           // response variable
}
parameters {
  real alpha_spp;
  real beta_spp;
}
transformed parameters {
  real mu_spp;
  
}
model {
  
  y ~ normal(mu_spp, sigma_spp);
  
}
'
#  y ~ normal(beta[1] + beta[2]*App + beta[3]*App + beta[4]*App, sigma);

### This code currently seems to be working to reproduce species means (only doing for practice; not for the paper)
speciesMeansModel <-'
data {
  int<lower=1> S; // number of species
  real y[S]; // estimated treatment effects
  real<lower=0> sigma[S]; // SE of the effect estimates
}
parameters {
  real mu; // overall means of all data
  real<lower=0> tau; // overall sd of all data
  vector[S] eta; // species-level errors
}
transformed parameters {
  vector[S] theta; // species effects
    theta <- mu + tau * eta;
}
model {
  eta ~ normal(0, 1);
  y ~ normal(theta, sigma);
}
'

speciesMeansFit <- stan(model_code = speciesMeansModel, data = species_dat, iter = 1000, chains = 4)
# Warning messages:
# 1: There were 75 divergent transitions after warmup. Increasing adapt_delta may help. 
# 2: Examine the pairs() plot to diagnose sampling problems
## NOTE: number of divergent transitions after warmup is not consistent; will change if you re-run the code; same is true if you use the file option below
pairs(speciesMeansFit)

## Increasing adapt delta to see if it helps with the divergent transitions
speciesMeansfitAD9 <- stan(model_code = speciesMeansModel, data = species_dat, iter = 1000, chains = 4, control = list(adapt_delta = 0.9))
# Nope, didn't help. 
# Now it also spits out:
# The following numerical problems occured the indicated number of times on chain  4 count
# Exception thrown at line 19: normal_log: Location parameter[1] is inf, but must be finite!     1
# If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.

#speciesMeansfit <- stan(file = 'speciesMeans.stan', data = species_dat, iter = 1000, chains = 4)
### Get warning messages
# Warning messages:
# 1: There were 21 divergent transitions after warmup. Increasing adapt_delta may help. 
# 2: Examine the pairs() plot to diagnose sampling problems

print(speciesMeansFit)
# theta means are very similar to the means of y by species (as we would expect)

## Evaluating whether there were problems in the model
# yellow points = transitions where the max treedepth__ was hit
# red points = transitions where n_divergent__ = 1 
# the below-diagonal intersection (draws with below-median accept_stat__) and the above diagonal intersection (draws with above-media accept_stat__) of the same two 
# variables should have distributions that are mirror images of each other. 
pairs(speciesMeansFit, pars = c("mu", "tau", "lp__"))
# below- and above-diagonal plots do not seem to mirror each other, so the data may be skewed?

# looking at the sampler directly, by each individual chain
lapply(get_sampler_params(speciesMeansFit, inc_warmup = TRUE), summary, digits = 2)
# max of n_divergent__ reaches one for each chain but the mean values are relatively small, so this suggests there were only a small number of divergent transitions? 
# the mean accept_stat__ is about 0.75 but the median is about 0.95; this suggests that it is pretty skewed


###### Attempting to account for individual variation ####
# Based on: http://rstudio-pubs-static.s3.amazonaws.com/64315_bc3a395edd104095a8384db8d9952f43.html
# So far I haven't gotten a working model figured out yet... 
speciesLookupVec <- unique(d[c("Ind", "Species")])[,"Species"]

dat <- with(d, 
            list(Ni = nrow(d), 
                 Nj = length(unique(Ind)),
                 Nk = length(unique(Species)),
                 indID = as.integer(Ind),
                 speciesID = as.integer(Species), 
                 speciesLookup = as.integer(speciesLookupVec),
                 y = response
            ))

speciesMeansIndModel <-'
data {
int<lower=1> Ni; // number of observations
int<lower=1> Nj; // number of individuals
int<lower=1> Nk; // number of species

int<lower=1> indID[Ni];
int<lower=1> speciesID[Ni];
int<lower=1> speciesLookup[Nj];
real y[Ni];
}
parameters {
real beta_0; // sample intercept
real beta_1; // sample slope
real<lower=0> sigma_e0; // sample error

real u_0jk[Nj];
real<lower=0> sigma_u0jk;

real u_0k[Nk];
real<lower=0> sigma_u0k; 
}
transformed parameters {
// Varying intercepts
real beta_0jk[Nj];
real beta_0k[Nk];

// Individual mean
real mu[Ni];

// Varying intercepts definition
for (k in 1:Nk) {
beta_0k[k] <- beta_0 + u_0k[k]; // species level
}
for (j in 1:Nj) {
beta_0jk[j] <- beta_0k[speciesLookup[j]] + u_0jk[j]; // individual level
}

// sample mean
for (i in 1:Ni) {
mu[i] <- beta_0jk[indID[i]]; 
}
}

model {
// assuming flat prior for now, so do not need to code anything to specify that

// Random effects distribution
u_0k ~ normal(0, sigma_u0k);
u_0jk ~ normal(0, sigma_u0jk);

// Likelihood 
// outcome model is N(mu, sigma^2) whereby sigma represents SD rather than variance
for (i in 1:Ni) {
y[i] ~ normal(mu[i], sigma_e0);
}
}
'

speciesIndStanFit <- stan(model_code = speciesMeansIndModel, data = dat,
                          chains = 4, iter = 1000, warmup = 1000, thin = 10)


