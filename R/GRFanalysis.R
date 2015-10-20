######## Analysis of FinLimb GRF daa ######

library(devtools)
install_github("MorphoFun/kraken")

data("FinLimbGRFs")

library(rstan)

##### SIMPLE EXAMPLE 
## from https://github.com/mclark--/Miscellaneous-R-Code/blob/master/ModelFitting/Bayesian/rstan_linregwithprior.R

cormat = matrix(c(1, .2, -.1, .3,
  .2, 1, .1, .2,
  -.1, .1, 1, .1,
  .3, .2, .1, 1),
  ncol=4, byrow=T)

cormat

# ensure pos def
library(Matrix)
cormat = nearPD(cormat, corr=T)$mat

### generate data ###

library(MASS)
means = rep(0, ncol(cormat))
n = 1000
d = mvrnorm(n, means, cormat, empirical=T)
colnames(d) = c('X1', 'X2', 'X3', 'y')
d[,'y'] = d[,'y'] -.1 # unnecessary, just to model a non-zero intercept
str(d)
cor(d)

### prepare for later processing ###

# strip X (add intercept column) and y for vectorized version later
X = cbind(1, d[,1:3]); colnames(X) = c('Intercept', 'X1', 'X2', 'X3')
y = d[,4]

# for comparison
modlm = lm(y~., data.frame(d))


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
