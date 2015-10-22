######## Analysis of FinLimb GRF daa ######


# Currently having issues with installing straight from GitHub
# install_github("MorphoFun/kraken")
# 
# data("FinLimbGRFs")

######### Analysis of FinLimb GRF data ##########
## All data were taken at the peak net GRF for each individual trial

library(rstan)
library(devtools)

# loading from local folder
load("FinLimbGRFs.rda")

########## GRF DATA ############
load("FinLimbGRFs.rda")

# stan can't handle factor types (e.g., Species in character format), so need to convert factor levels to dummy variables using as.integer()
FinLimbGRFs$SpeciesInt <- as.integer(FinLimbGRFs$Species)
d <- FinLimbGRFs
# converting percentages as proportions to make them closer in magnitude to the .BW data
d[,1] <- FinLimbGRFs[,1]/100 
# converting degrees to radians to make them closer in magnitude to the .BW data
d[,2:3] <- FinLimbGRFs[,2:3]*(pi/180)

# defining response vs. predictor for the stan model later
response <- d$NetGRF.BW
predictor <- d$SpeciesInt

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

######### REDUNDANCY ANALYSIS ###############
# Example: https://rstudio-pubs-static.s3.amazonaws.com/64619_2f93b223a318410bbf999d092ecf05ec.html
library(vegan)

speciesrda <- rda(FinLimbGRFs[,1:7] ~ Species, data = FinLimbGRFs)
# constrained variable only explains about 10% of the variance? yikes...
RsquareAdj(speciesrda) # pretty weak

speciesindrda <- rda(FinLimbGRFs[,1:7] ~ Species + Ind, data = FinLimbGRFs)
# Species with inds explains about 20% of the variance
RsquareAdj(speciesindrda)

plot(speciesrda, xlab = paste("RDA1 (", summary(speciesrda)$cont[[1]][2]*100, "%)"), 
     ylab = paste("RDA1 (", summary(speciesrda)$cont[[1]][5]*100, "%)"),
     type = "n"
     )
    #trialscores <- scores(speciesrda, choices=1:2, scaling = 2)

########## SELECTING VARIABLES FOR COMPARISON ##############
species_dat <- list(S = length(levels(FinLimbGRFs$Species)), 
                    y = as.numeric(unlist(aggregate(response, list(predictor), mean)[2])),
                    sigma = as.numeric(unlist(aggregate(response, list(predictor), sd)[2]))
)
# see page 11 of the Rstan: the R interface to Stan manual; do I need to convert y to an array? S should never equal 1, so I'm guessing the answer is "no". 

########## STAN MODEL FOR COMPARING SPECIES MEANS ###########
# Should the data and parameters be defined as vectors instead? https://groups.google.com/forum/#!msg/stan-users/4PgOF38Mnwk/hgUPCA768w0J
# Should vector types always be used for linear (algebra) problems?

speciesMeansModel <-'
data {
  int<lower=0> S; // number of species
  real y[S]; // estimated treatment effects
  real<lower=0> sigma[S]; // SE of the effect estimates
}
parameters {
  real mu;
  real<lower=0> tau;
  real eta[S];
}
transformed parameters {
  real theta[S];
  for (s in 1:S)
    theta[s] <- mu + tau * eta[s];
}
model {
  eta ~ normal(0, 1);
  y ~ normal(theta, sigma);
}
'
speciesMeansfitString <- stan(model_code = speciesMeansModel, data = species_dat, iter = 1000, chains = 4)
# Warning messages:
# 1: There were 75 divergent transitions after warmup. Increasing adapt_delta may help. 
# 2: Examine the pairs() plot to diagnose sampling problems
## NOTE: number of divergent transitions after warmup is not consistent; will change if you re-run the code; same is true if you use the file option below

## Increasing adapt delta to see if it helps with the divergent transitions
speciesMeansfitStringAD9 <- stan(model_code = speciesMeansModel, data = species_dat, iter = 1000, chains = 4, control = list(adapt_delta = 0.9))
# Nope, didn't help. 
# Now it also spits out:
# The following numerical problems occured the indicated number of times on chain  4 count
# Exception thrown at line 19: normal_log: Location parameter[1] is inf, but must be finite!     1
# If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.

speciesMeansfit <- stan(file = 'speciesMeans.stan', data = species_dat, iter = 1000, chains = 4)
### Get warning messages
# Warning messages:
# 1: There were 21 divergent transitions after warmup. Increasing adapt_delta may help. 
# 2: Examine the pairs() plot to diagnose sampling problems

print(speciesMeansfit)
# theta means are very similar to the means of y by species (as we would expect)

## Evaluating whether there were problems in the model
# yellow points = transitions where the max treedepth__ was hit
# red points = transitions where n_divergent__ = 1 
# the below-diagonal intersection (draws with below-median accept_stat__) and the above diagonal intersection (draws with above-media accept_stat__) of the same two 
# variables should have distributions that are mirror images of each other. 
pairs(speciesMeansfit, pars = c("mu", "tau", "lp__"))
# below- and above-diagonal plots do not seem to mirror each other, so the data may be skewed?

# looking at the sampler directly, by each individual chain
lapply(get_sampler_params(speciesMeansfit, inc_warmup = TRUE), summary, digits = 2)
# max of n_divergent__ reaches one for each chain but the mean values are relatively small, so this suggests there were only a small number of divergent transitions? 
# the mean accept_stat__ is about 0.75 but the median is about 0.95; this suggests that it is pretty skewed









##### SIMPLE EXAMPLE 
## based on code from https://github.com/mclark--/Miscellaneous-R-Code/blob/master/ModelFitting/Bayesian/rstan_linregwithprior.R

##########################
### UNVECTORIZED MODEL ###
##########################


### Stan related stuff ###


# Create the data list object
dat = list(N = nrow(d), y=y, X1=d[,1], X2=d[,2], X3=d[,3])

# Create the stan model object.
stanmodelcode <-'
data {                      // all of data noted here must be in the list that is imported
int<lower=0> N;           // Sample size
vector[N] X1;             // Predictor X1
vector[N] X2;
vector[N] X3;
vector[N] y;              // Response
}
parameters {                // which parameters will be estimated?
real alpha;               
real beta1;
real beta2;
real beta3;
real<lower=0> sigma;      
}
model {                     // Model setup of priors and likelihood
//priors
alpha ~ normal(0, 10);    // note that in Stan, normal(0,2) means distributed as mean 0 and standard deviation 2           
beta1 ~ normal(0, 10);
beta2 ~ normal(0, 10);
beta3 ~ normal(0, 10);
sigma ~ cauchy(0, 2.5);   // see Gelman 2006 or 2013 for example

//likelihood
for (n in 1:N)
y[n] ~ normal(alpha + beta1 * X1[n] + beta2 * X2[n] + beta3 * X3[n], sigma);
}
'

library(rstan)

### Run the model and examine results ###
# fit
fit <- stan(model_code = stanmodelcode, model_name = "example", 
            data = dat, iter = 12000, warmup=2000, thin=10, chains = 3, # sample_file = 'norm.csv', if you want to save
            verbose = F) 

# summary
print(fit, digits_summary=4)

# compare
summary(modlm)
