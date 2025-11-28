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
beta0  <- 2 * var(y)


#Posterior distribution inputs
LambdaN <- Lambda0 + t(X) %*% X
mN <- solve(LambdaN, Lambda0 %*% m0 + t(X) %*% y)
alphaN <- alpha0 + N / 2
betaN <- beta0 + 0.5 * (t(y) %*% y + t(m0) %*% Lambda0 %*% m0 - t(mN) %*% LambdaN %*% mN)


#Power Posterior Method

#Following the 'review' paper

#We set our temperature scale

#Then at each temperature we : 
#Sample from power posterior
#Estimate Ej = expectation of log likelhood using samples from power posterior at temp tj

#At end can use trapezoidal rule formula to estimate log evidence(need log sum exp trick for Ejs and log evidence calculations)

#Now create power log likelihood function  for the Ej's

power_like_log<- function(t, theta,X,y){
  out<-t*sum(dnorm(y, mean = X %*% theta[1:2], sd = sqrt(theta[3]), log = TRUE))
  return(out)
}

#N is number of (gibbs) samples for each tempered distribution
Nsim<-10000
#Sample from prior using STAN

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


#T is number of tempered distributions
T<- 10
#N is number of (gibbs) samples for each tempered distribution
Nsim<-10000

#Temperature scale:

#Tempering function- want uniform on log scale 
Beta_t<-function(t, c=2){
  return((t/T)^c)
}

#Gives T+1 tempered posteriors(inlcuding prior(t=0) and posterior(t=1))
t_list<-Beta_t(0:T, c=2)

#Need to store outputs- samples of mu,tau in seperate T+1 x N matrices

beta_samples<-array(0, dim = c(2,T + 1, Nsim))
sigma_sq_samples<- matrix(0,nrow = T+1, ncol = Nsim)

#Store initial samples in matrix
beta_samples[1,1,]<-prior_sample_Beta1
beta_samples[2,1,]<-prior_sample_Beta2
sigma_sq_samples[1,]<-prior_sample_sigma_sq


#Set E0 to expectation of log likelihood under prior

E<- matrix(0,nrow = T+1, ncol = Nsim)
for (i in 1:Nsim) {
  E[1,i]<- sum(dnorm(y, mean = X %*% beta_samples[,1,i], sd = sqrt(sigma_sq_samples[1,i]), log = TRUE))
  
}



#Now need to structure calculating weights then getting new samples in loops
for(t in 2:(T+1)){
  power<-t_list[t]
  init_list <- lapply(1:4, function(j) {
    idx <- sample(1:Nsim, 1)
    list(Beta = beta_samples[,t-1, idx], sigma_sq = sigma_sq_samples[t-1, idx])
  })
  run<- stan('Power Posterior Bayes Linear Regression.stan', 
             data =  list(N=length(y),d=d, X=X, m0=m0,alpha0=alpha0,beta0=beta0, Lambda0 = Lambda0,y=y, t=power), 
             iter =Nsim/2, 
             init = init_list)
  beta_samples[,t,]<- t(extract(run, pars = 'Beta')$'Beta')  
  sigma_sq_samples[t,]<- extract(run, pars = 'sigma_sq')$'sigma_sq'  
  for (i in 1:Nsim) {
    E[t,i]<- sum(dnorm(y, mean = X %*% beta_samples[,t,i], sd = sqrt(sigma_sq_samples[t,i]), log = TRUE))
  }
}



#Now work out power posterior estimate of log evidence
#Need to apply log sum exp trick to this to get it more accurate
power_post_est<-0
for(t in 2:(T+1)){
  power_post_est<- power_post_est + (t_list[t] - t_list[t-1]) *(0.5*(rowMeans(E)[t] + rowMeans(E)[t-1]))
}
#power_post_est<- N*power_post_est
power_post_est



end=proc.time()

timer<- end-start
timer[3]
