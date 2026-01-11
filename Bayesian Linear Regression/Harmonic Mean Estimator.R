library(rstan)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

setwd("~/Desktop/Project IV")

posterior_model <- stan_model(
  file = "Posterior Bayesian Linear Regression.stan"
)

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


#(Log) Harmonic mean estimator


#Sample from posterior using STAN

#Apply MCMC to model- get posterior mu samples
library(rstan)
library(durhamSLR)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())


setwd("~/Desktop/Project IV")   # set working directory to Project IV folder

#Stan Data
stan_data<- list(N=N,d=d, X=X, m0=m0,alpha0=alpha0,beta0=beta0, Lambda0 = Lambda0, y=y)

posterior_sample <- sampling(
  posterior_model,
  data = stan_data,
  iter = 50000,
  chains = 1,
  refresh = 0
)

#Diagnostic checks
#print(posterior_sample)
#output = as.array(posterior_sample)
#diagnostics(output)

#Extract Prior samples
post_sample_sigma_sq<- extract(posterior_sample, pars = 'sigma_sq')$'sigma_sq'  
post_sample_Beta1<- extract(posterior_sample, pars = 'Beta[1]')$'Beta[1]'
post_sample_Beta2<- extract(posterior_sample, pars = 'Beta[2]')$'Beta[2]'
post_sample_Beta<- cbind(post_sample_Beta1,post_sample_Beta2)

#Calculate likelihood at each posterior sample of (mu,tau)
likelihood_log<- rep(0, length(post_sample_Beta1))
for(i in 1:length(post_sample_Beta1)){
  mean<-as.vector(X %*% post_sample_Beta[i,])
  likelihood_log[i]<- sum(dnorm(y,mean = mean, sd = sqrt(post_sample_sigma_sq[i]) , log = TRUE))
}

#Calculate log harmonic mean using log sum exp trick for numerical stability
m<-max(likelihood_log)
hme_log<- log(length(post_sample_Beta1)) - log(exp(-m) * sum(exp(-likelihood_log + m)))
hme_log

end=proc.time()

timer<- end-start
timer[3]
