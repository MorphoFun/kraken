# Based primarily from: https://github.com/mclark--/Miscellaneous-R-Code/blob/master/ModelFitting/Bayesian/rstan_linregwithprior.R
# Right now I don't think the varying slope contribution from "Ind" is being incorporated correctly into the model; so far this is simple linear regression
# I need to edit this so that a different beta value is created for each species
### This code is currently not working; I need to mess with it some more.
# speciesMeansModel <- '
# data {
#   int<lower=0> N; // sample size
#   vector[N] S; // Predictor - species ID
#   vector[N] I; // "random effect" - individual ID
#   vector[N] y; // response
# }
# parameters {
#   real alpha; 
#   real beta1;
#   real<lower=0> sigma;
# }
# model {
# //priors
#   alpha ~ normal(0,10); // 0,10 thrown in arbitrarily as a broad prior
#   beta1 ~ normal(0,10);
#   sigma ~ cauchy(0, 2.5); // Gelman 2006 and 2013 have an example using this
# 
# //likelihood
# for (n in 1:N)
#   y[n] ~ normal(alpha + beta1*S[n], sigma);
# }
# '


data {
  int<lower=1> N; // # total number of observations
  int<lower=1> S; // # total number of species
  int<lower=1> I; // # total number of individuals
  int<lower=1> FE; // # total number of fixed effects
    
  int<lower=1, upper=S> spp[N]; // # species ID
  int<lower=1, upper=I> ind[N]; // # individual ID

  vector y[N]; // # response
  vector spp[N]; // # predictor
}

parameters {
  vector beta[FE];
  
  real<lower=0> sigma_FE; // # error from fixed effect
  real RE[I];   
  real<lower=0> sigma_RE; // # error from random effect
  
}

transformed parameters {
  real beta_RE[I]; // # Varying intercepts model with individual variation incorporated
  real mu[N]; // # individual mean
  
}


#### EXAMPLE 2 ####
SpeciesIndModel <- '
data {
  int<lower=1> N; // Sample size
  int<lower=1> K; // Dimension of model matrix (one column is a vector of ones for the intercept)   
  matrix[N,K] X; // Model Matrix
  vector[N] y; // response variable
}

transformed data {
// nothing yet
}

parameters {
vector[S] 
}
'

#### EXAMPLE 3 #####
modelString = "
data {
  int<lower=1> Ntotal ;
  int x[Ntotal] ;
  real y[Ntotal] ;
  real meanY ;
  real sdY ;
}
transformed data {
  real unifLo ;
  real unifHi ;
  real normalSigma ;
  real expLambda ;
  unifLo <- sdY/1000 ;
  unifHi <- sdY*1000 ;
  normalSigma <- sdY*100 ;
  expLambda <- 1/29.0 ;
}
parameters {
real<lower=0> nuMinusOne ; 
real mu[3] ;               // 3 groups
real<lower=0> sigma[3] ;   // 3 groups
}
transformed parameters {
real<lower=0> nu ;         // actually lower=1
nu <- nuMinusOne + 1 ;
}
model {
sigma ~ uniform( unifLo , unifHi ) ; // vectorized
mu ~ normal( meanY , normalSigma ) ; // vectorized
nuMinusOne ~ exponential( expLambda ) ;
for ( i in 1:Ntotal ) {
y[i] ~ normal( nu , mu[x[i]] , sigma[x[i]] ) ;
}
}
" 



#### SOURCE CODE FOR genMCMC ######

genMCMC2 <- function( datFrm, yName="y" , xName="x" , numSavedSteps=50000 , 
          saveName=NULL ) { 
  require(rstan)
  #-----------------------------------------------------------------------------
  # THE DATA.
  y = as.numeric(datFrm[,yName])
  x = as.numeric(as.factor(datFrm[,xName]))
  xLevels = levels(as.factor(datFrm[,xName]))
  # Do some checking that data make sense:
  if ( any( x!=1 & x!=2 ) ) { stop("All x values must be 1 or 2.") }
  if ( any( !is.finite(y) ) ) { stop("All y values must be finite.") }
  if ( length(x) != length(y) ) { stop("x and y must be same length.") }
  Ntotal = length(y)
  # Specify the data in a list, for later shipment to JAGS:
  dataList = list(
    y = y ,
    x = x ,
    Ntotal = Ntotal ,
    meanY = mean(y) ,
    sdY = sd(y)
  )
  #-----------------------------------------------------------------------------
  # THE MODEL.
  modelString = "
  data {
  int<lower=1> Ntotal ;
  int x[Ntotal] ;
  real y[Ntotal] ;
  real meanY ;
  real sdY ;
}
  transformed data {
  real unifLo ;
  real unifHi ;
  real normalSigma ;
  real expLambda ;
  unifLo <- sdY/1000 ;
  unifHi <- sdY*1000 ;
  normalSigma <- sdY*100 ;
  expLambda <- 1/29.0 ;
  }
  parameters {
  real<lower=0> nuMinusOne ; 
  real mu[2] ;               // 2 groups
  real<lower=0> sigma[2] ;   // 2 groups
  }
  transformed parameters {
  real<lower=0> nu ;         // actually lower=1
  nu <- nuMinusOne + 1 ;
  }
  model {
  sigma ~ uniform( unifLo , unifHi ) ; // vectorized
  mu ~ normal( meanY , normalSigma ) ; // vectorized
  nuMinusOne ~ exponential( expLambda ) ;
  for ( i in 1:Ntotal ) {
  y[i] ~ student_t( nu , mu[x[i]] , sigma[x[i]] ) ;
  }
  }
  " # close quote for modelString
  # Write out modelString to a text file
  writeLines( modelString , con="TEMPmodel.txt" )
  #-----------------------------------------------------------------------------
  # INTIALIZE THE CHAINS.
  # Initial values of MCMC chains based on data:
  mu = c( mean(y[x==1]) , mean(y[x==2]) )
  sigma = c( sd(y[x==1]) , sd(y[x==2]) )
  # Regarding initial values in next line: (1) sigma will tend to be too big if 
  # the data have outliers, and (2) nu starts at 5 as a moderate value. These
  # initial values keep the burn-in period moderate.
  initsList = list( mu = mu , sigma = sigma , nuMinusOne = 4 )
  #-----------------------------------------------------------------------------
  # RUN THE CHAINS
  parameters = c( "mu" , "sigma" , "nu" )     # The parameters to be monitored
  adaptSteps = 500               # Number of steps to "tune" the samplers
  burnInSteps = 1000
  nChains = 4 
  thinSteps = 1
  # Translate to C++ and compile to DSO:
  stanDso <- stan_model( model_code=modelString ) 
  # Get MC sample of posterior:
  stanFit <- sampling( object=stanDso , 
                       data = dataList , 
                       #pars = parameters , # optional
                       chains = nChains ,
                       iter = ( ceiling(numSavedSteps/nChains)*thinSteps
                                +burnInSteps ) , 
                       warmup = burnInSteps , 
                       #init = initsList , # optional
                       thin = thinSteps )
  # For consistency with JAGS-oriented functions in DBDA2E collection, 
  # convert stan format to coda format:
  codaSamples = mcmc.list( lapply( 1:ncol(stanFit) , 
                                   function(x) { mcmc(as.array(stanFit)[,x,]) } ) )
  # resulting codaSamples object has these indices: 
  #   codaSamples[[ chainIdx ]][ stepIdx , paramIdx ]
  if ( !is.null(saveName) ) {
    save( codaSamples , file=paste(saveName,"Mcmc.Rdata",sep="") )
    save( stanFit , file=paste(saveName,"StanFit.Rdata",sep="") )
    save( stanDso , file=paste(saveName,"StanDso.Rdata",sep="") )
  }  
  return( codaSamples )
}




stan_code <- '
data {
  // Define variables in data
  // Number of level-1 observations (an integer)
  int<lower=0> Ni;
  // Number of level-2 clusters
  int<lower=0> Nj;
  // Number of level-3 clusters
  int<lower=0> Nk;

  // Cluster IDs
  int<lower=1> classid[Ni];
  int<lower=1> schoolid[Ni];

  // Level 3 look up vector for level 2
  int<lower=1> schoolLookup[Nj];

  // Continuous outcome
  real mathgain[Ni];
  
  // Continuous predictor
  // real X_1ijk[Ni];
}

parameters {
  // Define parameters to estimate
  // Population intercept (a real number)
  real beta_0;
  // Population slope
  // real beta_1;

  // Level-1 errors
  real<lower=0> sigma_e0;

  // Level-2 random effect
  real u_0jk[Nj];
  real<lower=0> sigma_u0jk;

  // Level-3 random effect
  real u_0k[Nk];
  real<lower=0> sigma_u0k;
}

transformed parameters  {
  // Varying intercepts
  real beta_0jk[Nj];
  real beta_0k[Nk];

  // Individual mean
  real mu[Ni];

  // Varying intercepts definition
  // Level-3 (10 level-3 random intercepts)
  for (k in 1:Nk) {
    beta_0k[k] <- beta_0 + u_0k[k];
  }
  // Level-2 (100 level-2 random intercepts)
  for (j in 1:Nj) {
    beta_0jk[j] <- beta_0k[schoolLookup[j]] + u_0jk[j];
  }
  // Individual mean
  for (i in 1:Ni) {
    mu[i] <- beta_0jk[classid[i]];
  }
}

model {
  // Prior part of Bayesian inference
  // Flat prior for mu (no need to specify if non-informative)

  // Random effects distribution
  u_0k  ~ normal(0, sigma_u0k);
  u_0jk ~ normal(0, sigma_u0jk);

  // Likelihood part of Bayesian inference
  // Outcome model N(mu, sigma^2) (use SD rather than Var)
  for (i in 1:Ni) {
    mathgain[i] ~ normal(mu[i], sigma_e0);
  }
}'


#### VECTORIZED 8 SCHOOL EXAMPLE ####
schoolsind_data <- 
  list(J=8, 
       n = 80, 
       score = seq(1, 80),
       school = rep(1:8, 10)
  )

schoolindStan <- '
data {
  int<lower=0> J; // number of schools
  int<lower=0> n; //number of students
  real score[n]; // test score for each student
  int school[n]; //school for each student
}
parameters {
  real mu;
  real<lower=0> sigma;
}

model {
     vector[n] mu_school; 
     for (i in 1:n) 
       mu_school[i] <- mu[school[i]]; 
     y ~ normal(mu_school, sigma); 
}'

fit2 <- stan(model_code = schoolindStan, data = schoolsind_data, iter = 2000, chains =4)

#### UNVECTORIZED 8 SCHOOL EXAMPLE ###
schools_data <- 
  list(J=8, 
       y = c(28, 8, -3, 7, -1, 1, 18, 12),
       sigma = c(15, 10, 16, 11, 9, 11, 10, 18))

schoolStan <- ' data {
  int<lower=0> J; // number of schools
  real y[J]; // estimated treatment effects
  real<lower=0> sigma[J]; // s.e.’s of effect estimates
}
parameters {
  real mu; // population mean
  real<lower=0> tau; // population sd
  vector[J] eta; // school-level errors
}
transformed parameters {
  vector[J] theta; // school effects
  theta <- mu + tau*eta;
}
model {
  eta ~ normal(0, 1);
  y ~ normal(theta, sigma);
}'

fit1 <- stan(model_code = schoolStan, data = schools_data, iter = 2000, chains =4)