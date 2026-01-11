library(rstan)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

setwd("~/Desktop/Project IV")

prior_model <- stan_model(
  file = "Prior Bayesian Linear Regression.stan"
)

start=proc.time()
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



#Mean of likelihood using prior samples:

#Sample from prior using STAN

setwd("~/Desktop/Project IV")   # set working directory to Project IV folder

library(rstan)
library(durhamSLR)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

#Need lots of draws for it to work
prior_fit <- sampling(
  prior_model,
  data = list(
    N = N, d = d, X = X,
    m0 = m0, alpha0 = alpha0,
    beta0 = beta0, Lambda0 = Lambda0
  ),
  iter = 100000,
  chains = 1,
  algorithm = "Fixed_param",
  refresh = 0
)

#Diagnostic checks
print(prior_fit)
#output = as.array(prior_fit)
#diagnostics(output)


#Extract Prior samples
prior_sample_sigma_sq<- extract(prior_fit, pars = 'sigma_sq')$'sigma_sq'  
prior_sample_Beta1<- extract(prior_fit, pars = 'Beta[1]')$'Beta[1]'
prior_sample_Beta2<- extract(prior_fit, pars = 'Beta[2]')$'Beta[2]'
prior_sample_Beta<- cbind(prior_sample_Beta1,prior_sample_Beta2)

#Calculate log likelihood at each prior sample (mu,tau)

likelihood_prior_log<- rep(0, length(prior_sample_Beta1))
for(i in 1:length(prior_sample_Beta1)){
  mean<-as.vector(X %*% prior_sample_Beta[i,])
  likelihood_prior_log[i]<- sum(dnorm(y,mean = mean, sd = sqrt(prior_sample_sigma_sq[i]) , log = TRUE))
}

#Apply log sum exp trick to log(mean(likelihood))
m <- max(likelihood_prior_log)
stan_prior_le <- m + log(mean(exp(likelihood_prior_log - m)))
stan_prior_le

end=proc.time()

timer<- end-start
timer[3]
