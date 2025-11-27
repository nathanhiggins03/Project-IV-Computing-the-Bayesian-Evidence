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

#SMC

#Now create log power posterior function in R

power_post_log<- function(t, theta,X){
  out<-t*sum(dnorm(X, theta[1], sqrt(1/theta[2]), log = TRUE)) +dnorm(theta[1], m0, sqrt(1/(k0*theta[2])), log = TRUE) + dgamma(theta[2], alpha0,beta0, log = TRUE)
  return(out)
}

#T is number of tempered distributions
T<- 10

#N is number of (gibbs) samples for each tempered distribution
N<-2000

#Sample from joint prior

setwd("~/Desktop/Project IV")   # set working directory to Project IV folder

library(rstan)
library(durhamSLR)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

prior_fit = stan('Prior Normal Normal with unknown mean and precision.stan', 
                 data = list(k0=k0, m0=m0,alpha0=alpha0,beta0=beta0), 
                 iter =4000, chains = 1, algorithm = "Fixed_param")

#Extract Prior samples
prior_sample_mu<- extract(prior_fit, pars = 'mu')$'mu'  
prior_sample_tau<- extract(prior_fit, pars = 'tau')$'tau'


#Need to store outputs- samples of mu,tau in seperate T+1 x N matrices

mu_samples<-matrix(0,nrow = T+1, ncol = N)
tau_samples<- matrix(0,nrow = T+1, ncol = N)

#Store initial samples in matrix
mu_samples[1,]<-prior_sample_mu
tau_samples[1,]<-prior_sample_tau

#Need to store log of unnormalised weights
log_w<-matrix(0,nrow = T, ncol = N)

#Tempering function- want uniform on log scale 
Beta_t<-function(t, c=2){
  return((t/T)^c)
}

#Gives T+1 tempered posteriors(inlcuding prior(t=0) and posterior(t=1))
t_list<-Beta_t(0:T)




#SMC loop- weight, normalise, resample
for(t in 2:(T+1)){
  power<-t_list[t]
  for (i in 1:N) {
    #Weight
    log_w[t-1,i]<- power_post_log(power, c(mu_samples[t-1,i], tau_samples[t-1,i]), x) - power_post_log(t_list[t-1], c(mu_samples[t-1,i], tau_samples[t-1,i]), x)
  }
  # Log sum exp trick
  max_logw <- max(log_w[t-1, ])
  w <- exp(log_w[t-1, ] - max_logw)
  #Explicit normalisation not needed- R can handle unnormalised weights
  #norm_w <- w / sum(w)
  #Would swap w for norm_w in resample below
  
  # Joint resampling
  resample_index <- sample(1:N, size = N, replace = TRUE, prob = w)
  mu_samples[t, ]  <- mu_samples[t-1, resample_index]
  tau_samples[t, ] <- tau_samples[t-1, resample_index]
}

#Calculating log evidence evidence

smc_log <- 0
for (t in 1:T) {
  max_logw <- max(log_w[t, ])
  smc_log <- smc_log + log(mean(exp(log_w[t, ] - max_logw))) + max_logw
}
smc_log
