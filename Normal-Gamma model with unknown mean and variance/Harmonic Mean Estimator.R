
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

#(Log) Harmonic mean estimator


#Sample from posterior using STAN

#Apply MCMC to model- get posterior mu samples
library(rstan)
library(durhamSLR)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())


setwd("~/Desktop/Project IV")   # set working directory to Project IV folder

#Stan Data
stan_data<- list(x=x, N=length(x),k0=k0, m0=m0,alpha0=alpha0,beta0=beta0)

posterior_sample <- stan(
  file = "Posterior Normal Normal model with unknown mean and precision.stan",
  data = stan_data, iter = 10000)

#Diagnostic checks
#print(posterior_sample)
#output = as.array(posterior_sample)
#diagnostics(output)

#Extract Posterior samples
posterior_sample_mu<- extract(posterior_sample, pars = 'mu')$'mu' 
posterior_sample_tau<- extract(posterior_sample, pars = 'tau')$'tau'

#Calculate likelihood at each posterior sample of (mu,tau)
likelihood_log<- rep(0, length(posterior_sample_mu))
for(i in 1:length(posterior_sample_mu)){
  likelihood_log[i]<- sum(dnorm(x,posterior_sample_mu[i], sqrt(1/posterior_sample_tau[i]), log = TRUE))
}

#Calculate log harmonic mean using log sum exp trick for numerical stability
m<-max(likelihood_log)
hme_log<- log(length(posterior_sample_mu)) - log(exp(-m) * sum(exp(-likelihood_log + m)))
#Tends to be larger as HME overestimates as it is statistical unstable
hme_log

