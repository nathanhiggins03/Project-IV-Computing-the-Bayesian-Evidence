start=proc.time()

#Bayesian Linear Regression with unknown beta and precision

#Set seed for reproducibility
set.seed(123)
#Data
#y=2+3x+ϵ,ϵ∼N(0,0.5^2)

N <- 30
x <- runif(N, 0, 5)
X <- cbind(1, x)   # include intercept
beta_true <- c(2, 3)
sigma_true <- 0.5
y <- as.vector(X %*% beta_true + rnorm(N, 0, sigma_true))

#Data we have is
X
y

#Prior inputs
d <- ncol(X)
m0 <- rep(1, d)
Lambda0 <- diag(5, d)   # vague prior
alpha0 <- 3
beta0  <- 36


#Posterior distribution inputs
LambdaN <- Lambda0 + t(X) %*% X
mN <- solve(LambdaN, Lambda0 %*% m0 + t(X) %*% y)
alphaN <- alpha0 + N / 2
betaN <- beta0 + 0.5 * (t(y) %*% y + t(m0) %*% Lambda0 %*% m0 - t(mN) %*% LambdaN %*% mN)



#SMC

#Now create log power posterior function in R
library(mvtnorm)
library(extraDistr)

power_post_log<- function(t, theta,X,y){
  Beta<- theta[1:2]
  sigma_sq<-theta[3]
  mean<-X %*% Beta
  sd <- sqrt(sigma_sq)
  out <- t*sum(dnorm(y, mean = mean, sd = sd, log = TRUE)) +
    dmvnorm(Beta, mean = m0, sigma = sigma_sq * solve(Lambda0), log = TRUE) +
    dinvgamma(sigma_sq, alpha0, beta0, log = TRUE)
  return(out)
}

#T is number of tempered distributions
T<- 50

#N is number of (gibbs) samples for each tempered distribution
Nsim<-20000

#Sample from joint prior

setwd("~/Desktop/Project IV")   # set working directory to Project IV folder

library(rstan)
library(durhamSLR)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

#Need lots of draws for it to work
prior_fit = stan('Prior Bayesian Linear Regression.stan', 
                 data = list(N=length(y),d=d, X=X, m0=m0,alpha0=alpha0,beta0=beta0, Lambda0 = Lambda0), 
                 iter =200000, chains = 1, algorithm = "Fixed_param")

#Extract Prior samples
prior_sample_sigma_sq_all<- extract(prior_fit, pars = 'sigma_sq')$'sigma_sq'  
prior_sample_Beta1_all<- extract(prior_fit, pars = 'Beta[1]')$'Beta[1]'
prior_sample_Beta2_all<- extract(prior_fit, pars = 'Beta[2]')$'Beta[2]'
prior_sample_Beta_all<- cbind(prior_sample_Beta1_all,prior_sample_Beta2_all)

prior_sample_sigma_sq<- tail(prior_sample_sigma_sq_all, Nsim)
prior_sample_Beta1<- tail(prior_sample_Beta1_all, Nsim)
prior_sample_Beta2<- tail(prior_sample_Beta2_all, Nsim)
prior_sample_Beta<- tail(prior_sample_Beta_all, Nsim)

#Need to store outputs- samples of mu,tau in seperate T+1 x N matrices

beta_samples<-array(0, dim = c(2,T + 1, Nsim))
sigma_sq_samples<- matrix(0,nrow = T+1, ncol = Nsim)

#Store initial samples in matrix
beta_samples[1,1,]<-prior_sample_Beta1
beta_samples[2,1,]<-prior_sample_Beta2
sigma_sq_samples[1,]<-prior_sample_sigma_sq



#Need to store log of unnormalised weights
log_w<-matrix(0,nrow = T, ncol = Nsim)

#Tempering function- want uniform on log scale 
Beta_t<-function(t, c=2){
  return((t/T)^c)
}

#Gives T+1 tempered posteriors(inlcuding prior(t=0) and posterior(t=1))
t_list<-Beta_t(0:T, c=4)
#t_list <- exp(seq(log(1e-6), log(1), length.out = T+1))




#SMC loop- weight, normalise, resample
for(t in 2:(T+1)){
  power<-t_list[t]
  for (i in 1:Nsim) {
    #Weight
    log_w[t-1,i]<- power_post_log(power, c(beta_samples[1,t-1,i], beta_samples[2,t-1,i],sigma_sq_samples[t-1,i]), X,y) - power_post_log(t_list[t-1], c(beta_samples[1,t-1,i], beta_samples[2,t-1,i],sigma_sq_samples[t-1,i]), X,y)
  }
  # Log sum exp trick
  max_logw <- max(log_w[t-1, ])
  w <- exp(log_w[t-1, ] - max_logw)
  #Explicit normalisation not needed- R can handle unnormalised weights
  #norm_w <- w / sum(w)
  #Would swap w for norm_w in resample below
  
  # Joint resampling
  resample_index <- sample(1:Nsim, size = Nsim, replace = TRUE, prob = w)
  beta_samples[,t, ]  <- beta_samples[,t-1, resample_index]
  sigma_sq_samples[t, ] <- sigma_sq_samples[t-1, resample_index]
}

#Calculating log evidence evidence

smc_log <- 0
for (t in 1:T) {
  max_logw <- max(log_w[t, ])
  smc_log <- smc_log + log(mean(exp(log_w[t, ] - max_logw))) + max_logw
}
smc_log



end=proc.time()

timer<- end-start
timer[3]
