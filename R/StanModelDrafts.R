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

######### STAN CODE FOR SINGLE RANDOM EFFECT #########
stancode <- '
data{
int<lower=0> N;
int<lower=0> K;
matrix[N,K] X;
vector[N] price;
int J;
int<lower=1,upper=J> re[N];
}
parameters{
vector[J] a;
real mu_a;
real tau;
real<lower=0> sigma_a;
real<lower=0> sigma;
vector[K] beta;
}
transformed parameters{
vector[N] mu_hat;
for(i in 1:N)
mu_hat[i] <- a[re[i]];
}
model {
mu_a ~ normal(0,10);
tau ~ cauchy(0,5);
a ~ normal(mu_a,sigma_a);
for(i in 1:N)
price[i] ~ normal(X[i]*beta + mu_hat[i], sigma);
}'


#### simple linear regression ####
data(ChickWeight)
dat <- list(
    N <- nrow(ChickWeight),
    x <- ChickWeight$Time,
    y <- ChickWeight$weight
  )

linreg <- '
data {
  int<lower=0> N;
  vector[N] x; 
  vector[N] y;
}
parameters {
  real alpha;
  real beta;
  real<lower=0> sigma;
}
model {
  y ~ normal(alpha + beta * x, sigma);
}
'

ChickWeightFit <- stan(model_code = linreg, data = dat, iter = 1000, chains = 4)
ChickWeightLM <- lm(Weight ~ Time, data = ChickWeight)


#####   multivariate mixture model ########
#### based on: https://kylenpayne.wordpress.com/2014/01/10/bayesian-normal-mixture-models-with-rstan/

#Simulate Normal Distributions
library(mvtnorm)

### first cluster
mu1<-c(0,0,0,0)
sigma1<-matrix(c(1,0,0,0,
                 0,1,0,0,
                 0,0,1,0,
                 0,0,0,1), ncol=4, nrow=4, byrow=TRUE)

norm1<-rmvnorm(30, mean = mu1, sigma=sigma1)
### second cluster

mu2<-c(2,2,2,2)
sigma2<-matrix(c(1,0,0,0,
                 0,1,0,0,
                 0,0,1,0,
                 0,0,0,1
) + rnorm(1,0,.1), ncol=4, nrow=4, byrow=TRUE)

norm2<-rmvnorm(30, mean=mu2, sigma=sigma2)
### third cluster

mu3 <-c(5,3.4,2,0)
sigma3<-matrix(c(5,0,0,0,
                 0,2.5,0,0,
                 0,0,1.25,0,
                 0,0,0,3) + rnorm(1,0,.2), ncol=4, nrow=4, byrow=TRUE)

norm3<-rmvnorm(30, mean=mu3, sigma=sigma3)
norms<-rbind(norm1,norm2,norm3)

N = 90
Dim = 4
y<-array(as.vector(norms), dim=c(N,Dim))

#Fit Mixture Models

mixture_dat <-list(
  N = N,
  D = 4, K = 3, y=y )

mixture_model <-'
data {
  int<lower=1> D; // number of variables
  int<lower=1> K; // number of groups
  int<lower=1> N; // number of rows
  vector[D] y[N]; //observations
}

parameters {
  simplex[K] theta;
  vector[D] mu[K];
  real<lower=0,upper=10> sigma[K];  // scales of mixture components
}

model {
  real ps[K];
  //for (k in 1:K) {
    //  mu[k] ~ normal(0,10);
    //}
  mu[1] ~ normal(0,1);
  mu[2] ~ normal(2,1);
  mu[3] ~ normal(4,1);
  
  for (n in 1:N){
    for(k in 1:K){
      
      ps[k] <- log(theta[k])
      
      + normal_log(y[n],mu[k],sigma[k]);
    }
    
    increment_log_prob(log_sum_exp(ps));
  }
}'

fit.1 <- stan(model_code = mixture_model, model_name ="Mixture", data = mixture_dat)



########## READING TIME - SIMPLE LINEAR REGRESSION ########
#### http://arxiv.org/abs/1506.06201
## read in data:
rDat <- readGitHub("https://raw.githubusercontent.com/vasishth/BayesLMMTutorial/master/data/gibsonwu2012data.txt")
## subset critical region:
rDat <- subset( rDat , region == "headnoun" )

## create data as list for Stan, and fit model:
stanDat <- list( rt = rDat$rt, so = ifelse(rDat$type=="subj-ext", -1, 1), N = nrow(rDat) )

stanMod <- '
data {
  int<lower=1> N;                //number of data points
  real rt[N];                    //reading time
  real<lower=-1,upper=1> so[N];  //predictor
}

parameters {
  vector[2] beta;            //intercept and slope
  real<lower=0> sigma_e;     //error sd
}

model {
  real mu;
  for (i in 1:N){                   // likelihood
    mu <- beta[1] + beta[2]*so[i];
    rt[i] ~ lognormal(mu,sigma_e);
  }
}
'
library(rstan)
fixEfFit <- stan ( model_code = stanMod, data = stanDat, iter = 2000 , chains = 4 )

## plot traceplot, excluding warm-up:
traceplot( fixEfFit , pars = c("beta","sigma_e"), inc_warmup = FALSE)

## examine quantiles of posterior distributions:
print( fixEfFit , pars = c("beta","sigma_e") ,
probs = c(0.025,0.5,0.975))

## examine quantiles of parameter of interest:
beta1 <- extract ( fixEfFit , pars=c("beta[2]"))$beta
print ( signif ( quantile ( beta1,probs = c(0.025,0.5,0.975)), 2))

##### READING TIME - RANDOM INTERCEPTS MODEL ####

## format data for Stan:
stanDat<-list(subj=as.integer(factor(rDat$subj)),
                item=as.integer(factor(rDat$item)),
                rt=rDat$rt,
                so = ifelse(rDat$type=="subj-ext", -1, 1),
                N=nrow(rDat),
                J=length(unique(rDat$subj)),
                K=length(unique(rDat$item)))

ranInt <- '
data {
int<lower=1> N;                  //number of data points
real rt[N];                      //reading time
real<lower=-1,upper=1> so[N];    //predictor
int<lower=1> J;                  //number of subjects
int<lower=1> K;                  //number of items
int<lower=1, upper=J> subj[N];   //subject id
int<lower=1, upper=K> item[N];   //item id
}

parameters {
vector[2] beta;            //fixed intercept and slope
vector[J] u;               //subject intercepts
vector[K] w;               //item intercepts
real<lower=0> sigma_e;     //error sd
real<lower=0> sigma_u;     //subj sd
real<lower=0> sigma_w;     //item sd
}

model {
real mu;
//priors
u ~ normal(0,sigma_u);    //subj random effects
w ~ normal(0,sigma_w);    //item random effects
// likelihood
for (i in 1:N){
mu <- beta[1] + u[subj[i]] + w[item[i]] + beta[2]*so[i];
rt[i] ~ lognormal(mu,sigma_e);
}
}'

## Sample from posterior distribution:
ranIntFit <- stan(model_code = ranInt, data=stanDat,
                     iter=2000, chains=4)
## Summarize results:
print(ranIntFit,pars=c("beta","sigma_e","sigma_u","sigma_w"),
         probs=c(0.025,0.5,0.975))

beta1 <- extract(ranIntFit,pars=c("beta[2]"))$beta
print(signif(quantile(beta1,probs=c(0.025,0.5,0.975)),2))

## Posterior probability of beta1 being less than 0:
mean(beta1<0)


##### READING TIME - RANDOM INTERCEPTS MODEL WITHOUT ITEM ####

## format data for Stan:
stanDatRINoItem<-list(subj=as.integer(factor(rDat$subj)),
              rt=rDat$rt,
              so = ifelse(rDat$type=="subj-ext", -1, 1),
              N=nrow(rDat),
              J=length(unique(rDat$subj)))

ranIntNoItem <- '
data {
int<lower=1> N;                  //number of data points
real rt[N];                      //reading time
real<lower=-1,upper=1> so[N];    //predictor
int<lower=1> J;                  //number of subjects
int<lower=1, upper=J> subj[N];   //subject id
}

parameters {
vector[2] beta;            //fixed intercept and slope
vector[J] u;               //subject intercepts
real<lower=0> sigma_e;     //error sd
real<lower=0> sigma_u;     //subj sd
}

model {
real mu;
//priors
u ~ normal(0,sigma_u);    //subj random effects
// likelihood
for (i in 1:N){
mu <- beta[1] + u[subj[i]] + beta[2]*so[i];
rt[i] ~ lognormal(mu,sigma_e);
}
}'

## Sample from posterior distribution:
ranIntFit_noItem <- stan(model_code = ranIntNoItem, data=stanDatRINoItem,
                     iter=2000, chains=4)
print(ranIntFit_noItem)

## Summarize results:
print(ranIntFit_noItem,pars=c("beta","sigma_e","sigma_u"),
probs=c(0.025,0.5,0.975))

beta1 <- extract(ranIntFit_noItem,pars=c("beta[2]"))$beta
print(signif(quantile(beta1,probs=c(0.025,0.5,0.975)),2))

## Posterior probability of beta1 being less than 0:
mean(beta1<0)
