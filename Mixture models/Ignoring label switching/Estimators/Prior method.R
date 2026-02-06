library(rstan)
library(durhamSLR)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

setwd("~/Desktop/Project IV")

prior_model <- stan_model(
  file = "Prior Mixture model.stan"
)

start=proc.time()
#Set seed for reproducibility
set.seed(123)

#Data
data(galaxies, package = "MASS")
K <- 3   # or any K you want to test
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



#Mean of likelihood using prior samples:

#Sample from prior using STAN

setwd("~/Desktop/Project IV")   # set working directory to Project IV folder



#Need lots of draws for it to work
prior_fit <- sampling(
  prior_model,
  data = list(
    N = N,
    K = K,
    y = y,
    
    alpha = alpha,        # uniform Dirichlet
    mu0 = mu0,     # data-centered prior
    lambda0 = lambda0,         # weak prior on means
    a0 = a0,                # weak Inv-Gamma
    b0 = b0
  ),
  iter = 3000000,
  chains = 1,
  algorithm = "Fixed_param",
  refresh = 0
)

#Diagnostic checks
#print(prior_fit)
#output = as.array(prior_fit)
#diagnostics(output)


#Extract Prior samples
prior_sample_sigma_sq<- extract(prior_fit, pars = 'sigma2')$'sigma2'  
prior_sample_omega<- extract(prior_fit, pars = 'omega')$'omega'
prior_sample_mu<- extract(prior_fit, pars = 'mu')$'mu'
prior_sample_z<- extract(prior_fit, pars = 'z')$'z'

#Calculate log likelihood at each prior sample (mu,tau)

S <- dim(prior_sample_mu)[1]
N <- length(y)
K <- ncol(prior_sample_mu)

likelihood_prior_log <- rep(0, S)

for (s in 1:S) {
  
  likelihood_prior_log[s] <- sum(
    sapply(y, function(yi) {
      
      # mixture log density at yi
      m <- max(
        log(prior_sample_omega[s, ]) +
          dnorm(
            yi,
            mean = prior_sample_mu[s, ],
            sd   = sqrt(prior_sample_sigma_sq[s, ]),
            log  = TRUE
          )
      )
      
      m + log(sum(exp(
        log(prior_sample_omega[s, ]) +
          dnorm(
            yi,
            mean = prior_sample_mu[s, ],
            sd   = sqrt(prior_sample_sigma_sq[s, ]),
            log  = TRUE
          ) - m
      )))
    }))
}


#Apply log sum exp trick to log(mean(likelihood))
m <- max(likelihood_prior_log)
stan_prior_le <- m + log(mean(exp(likelihood_prior_log - m)))
stan_prior_le

end=proc.time()

timer<- end-start
timer[3]

#Can compare prior method to "arithmetic mean" in paper for each K, know we are
#reproducing correct results
