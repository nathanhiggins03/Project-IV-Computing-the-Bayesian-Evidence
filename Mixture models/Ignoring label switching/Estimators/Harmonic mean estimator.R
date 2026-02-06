library(rstan)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

setwd("~/Desktop/Project IV")

posterior_model <- stan_model(
  file = "Posterior Mixture model.stan"
)

start=proc.time()

#Bayesian Linear Regression with unknown beta and precision

#Set seed for reproducibility
set.seed(123)

#Data
data(galaxies, package = "MASS")
K <- 3   # or any K you want to test
#Note paper divides units by 1000
y<- galaxies/1000



#Prior inputs
N = length(y)
K = K
y = y

alpha = rep(1, K)        # uniform Dirichlet
mu0 = mean(y)     # data-centered prior
lambda0 = 2.6/(max(y)-min(y))           # weak prior on means
a0 = 1.28                   # weak Inv-Gamma
b0 = 0.36*(mean(y^2) - (mean(y)^2))

#(Log) Harmonic mean estimator


#Sample from posterior using STAN

#Apply MCMC to model- get posterior mu samples
library(rstan)
library(durhamSLR)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())


setwd("~/Desktop/Project IV")   # set working directory to Project IV folder

#Stan Data
stan_data<- list(
  N = N,
  K = K,
  y = y,
  
  alpha = alpha,        # uniform Dirichlet
  mu0 = mu0,     # data-centered prior
  lambda0 = lambda0,         # weak prior on means
  a0 = a0,                # weak Inv-Gamma
  b0 = b0
)

posterior_sample <- sampling(
  posterior_model,
  data = stan_data,
  iter = 4000,
  chains = 4,
  refresh = 0
)

#Diagnostic checks
#print(posterior_sample)
#output = as.array(posterior_sample)
#diagnostics(output)

#Extract Prior samples
post_sample_sigma_sq<- extract(posterior_sample, pars = 'sigma2')$'sigma2'  
post_sample_omega<- extract(posterior_sample, pars = 'omega')$'omega'
post_sample_mu<- extract(posterior_sample, pars = 'mu')$'mu'


#Calculate log likelihood at each prior sample (mu,tau)

S <- dim(post_sample_mu)[1]
N <- length(y)
K <- ncol(post_sample_mu)

likelihood_log <- rep(0, S)

for (s in 1:S) {
  
  likelihood_log[s] <- sum(
    sapply(y, function(yi) {
      
      # mixture log density at yi
      m <- max(
        log(post_sample_omega[s, ]) +
          dnorm(
            yi,
            mean = post_sample_mu[s, ],
            sd   = sqrt(post_sample_sigma_sq[s, ]),
            log  = TRUE
          )
      )
      
      m + log(sum(exp(
        log(post_sample_omega[s, ]) +
          dnorm(
            yi,
            mean = post_sample_mu[s, ],
            sd   = sqrt(post_sample_sigma_sq[s, ]),
            log  = TRUE
          ) - m
      )))
    }))
}


#Calculate log harmonic mean using log sum exp trick for numerical stability
m<-max(likelihood_log)
hme_log<- log(S) - log(exp(-m) * sum(exp(-likelihood_log + m)))
hme_log

end=proc.time()

timer<- end-start
timer[3]



