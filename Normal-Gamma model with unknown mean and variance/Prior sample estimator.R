#Sample from prior using STAN
setwd("~/Desktop/Project IV")   # set working directory to Project IV folder
library(rstan)
library(durhamSLR)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# Compile the prior model once
prior_model <- stan_model("Prior Normal Normal with unknown mean and precision.stan")

start=proc.time()
#Normal-Gamma model with unknown mean and variance

#Set seed for reproducibility
set.seed(123)

#Data
x<-rnorm(10,mean = 10,sd=2)

#Inputs
m0<-10
k0<-0.1
alpha0<-1
beta0<-0.1
n<- length(x)

#Tractable posterior inputs
k1<- k0 +n
m1<- (k0*m0 + n*mean(x))/k1
alpha1<- alpha0 + n/2
beta1<- beta0 + 1/2 * (n-1) * var(x) + (k0 *n * (mean(x) - m0)^2)/(2 * k1)


#Mean of likelihood using prior samples:



prior_fit <- sampling(
  prior_model,
  data = list(
    k0 = k0,
    m0 = m0,
    alpha0 = alpha0,
    beta0 = beta0
  ),
  iter = 1000,
  chains = 1,
  algorithm = "Fixed_param",
  refresh = 0
)

#Diagnostic checks
#print(prior_fit)
#output = as.array(prior_fit)
#diagnostics(output)

#Extract Prior samples
prior_sample_mu<- extract(prior_fit, pars = 'mu')$'mu'  
prior_sample_tau<- extract(prior_fit, pars = 'tau')$'tau'

#Calculate log likelihood at each prior sample (mu,tau)
likelihood_prior_log<- rep(0, length(prior_sample_mu))
for(i in 1:length(prior_sample_mu)){
  likelihood_prior_log[i]<- sum(dnorm(x,prior_sample_mu[i], sqrt(1/prior_sample_tau[i]), log = TRUE))
}

#Apply log sum exp trick to log(mean(likelihood))
m <- max(likelihood_prior_log)
stan_prior_le <- m + log(mean(exp(likelihood_prior_log - m)))
stan_prior_le


end=proc.time()

timer<- end-start
timer[3]
