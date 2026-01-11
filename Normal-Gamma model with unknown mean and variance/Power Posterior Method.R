setwd("~/Desktop/Project IV")   # set working directory to Project IV folder

library(rstan)
library(durhamSLR)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
prior_model <- stan_model('Prior Normal Normal with unknown mean and precision.stan')
power_model <- stan_model('Power Posterior Normal Normal model with unknown mean and precision.stan')


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

#Power Posterior Method

#We set our temperature scale

#Then at each temperature we : 
#Sample from power posterior
#Estimate Ej = expectation of log likelhood using samples from power posterior at temp tj

#At end can use trapezoidal rule formula to estimate log evidence(need log sum exp trick for Ejs and log evidence calculations)

#Now create power log likelihood function  for the Ej's

power_like_log<- function(t, theta,X){
  out<-t*sum(dnorm(X, theta[1], sqrt(1/theta[2]), log = TRUE)) 
  return(out)
}

#N is number of (gibbs) samples for each tempered distribution
N<-4000

#Sample from joint prior


prior_fit <- sampling(
  prior_model,
  data = list(k0=k0, m0=m0, alpha0=alpha0, beta0=beta0),
  iter = N,
  chains = 1,
  algorithm = "Fixed_param",
  refresh = 0
)
#Extract Prior samples
prior_sample_mu<- extract(prior_fit, pars = 'mu')$'mu'  
prior_sample_tau<- extract(prior_fit, pars = 'tau')$'tau'







#T is number of tempered distributions
T<- 10


#Temperature scale:

#Tempering function- want uniform on log scale 
Beta_t<-function(t, c=2){
  return((t/T)^c)
}

#Gives T+1 tempered posteriors(inlcuding prior(t=0) and posterior(t=1))
t_list<-Beta_t(0:T)

#Need to store outputs- samples of mu,tau in seperate T+1 x N matrices

mu_samples<-matrix(0,nrow = T+1, ncol = N)
tau_samples<- matrix(0,nrow = T+1, ncol = N)

#Store initial samples in matrix
mu_samples[1,]<-prior_sample_mu
tau_samples[1,]<-prior_sample_tau

#Set E0 to expectation of log likelihood under prior
E<- matrix(0,nrow = T+1, ncol = N)
for (i in 1:N) {
  E[1,i]<- sum(dnorm(x, mu_samples[1,i], sqrt(1/tau_samples[1,i]), log = TRUE)) 
  
}



N_samples_warmup <- 0
N_samples_kept   <- N
#Now need to structure calculating weights then getting new samples in loops
for(t in 2:(T+1)){
  power<-t_list[t]
  #run<- stan('Power Posterior Normal Normal model with unknown mean and precision.stan', 
  #           data =  list(x=x, N=length(x),k0=k0, m0=m0,alpha0=alpha0,beta0=beta0, t=power), iter =N/2)
  run <- sampling(
    power_model,
    data = list(x=x, N=length(x), k0=k0, m0=m0, alpha0=alpha0, beta0=beta0, t=power),
    iter = N_samples_kept,
    warmup = N_samples_warmup,
    chains = 1,
    refresh = 0
  )
  
  mu_samples[t,]<- extract(run, pars = 'mu')$'mu'  
  tau_samples[t,]<- extract(run, pars = 'tau')$'tau'  
  for (i in 1:N) {
    E[t,i]<- sum(dnorm(x, mu_samples[t,i], sqrt(1/tau_samples[t,i]), log = TRUE)) 
  }
}


#Now work out power posterior estimate of log evidence
#Need to apply log sum exp trick to this to get it more accurate
power_post_est<-0
for(t in 2:(T+1)){
  power_post_est<- power_post_est + (t_list[t] - t_list[t-1]) *(0.5*(rowMeans(E)[t] + rowMeans(E)[t-1]))
}
power_post_est

end=proc.time()

timer<- end-start
timer[3]
