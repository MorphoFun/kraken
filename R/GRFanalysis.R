######## Analysis of FinLimb GRF daa ######

library(devtools)
# I'm having issues with installing kraken from github... 
# install_github("MorphoFun/kraken")

load("FinLimbGRFs.rda")

library(rstan)

# stan can't handle factor types (e.g., Species in character format), so need to convert factor levels to dummy variables using as.integer()
FinLimbGRFs$SpeciesInt <- as.integer(FinLimbGRFs$Species)
d <- FinLimbGRFs
response <- FinLimbGRFs$NetGRF.BW
predictor <- FinLimbGRFs$SpeciesInt

species_dat <- list(S = length(levels(FinLimbGRFs$Species)), 
                    y = as.numeric(unlist(aggregate(response, list(predictor), mean)[2])),
                    sigma = as.numeric(unlist(aggregate(response, list(predictor), sd)[2]))
)
# see page 11 of the Rstan: the R interface to Stan manual; do I need to convert y to an array? S should never equal 1, so I'm guessing the answer is "no". 

### Stan code
speciesModel <-'
data {
  int<lower=0> S; // number of species
  real y[S]; // estimated treatment effects
  real<lower=0> sigma[S]; SE of the effect estimates
}
parameters {
  real mu;
  real<lower=0> tau;
  vector[S] eta;
}
transformed parameters {
  vector[S] theta;
  theta <- mu + tau * eta;
}
model {
  eta ~ normal(0,1);
  y ~ normal(theta, sigma);
}
'

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
