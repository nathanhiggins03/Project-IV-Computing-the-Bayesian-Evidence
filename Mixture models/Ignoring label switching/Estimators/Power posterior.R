setwd("~/Desktop/Project IV")

library(rstan)


rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())


# Compile Stan models

prior_model <- stan_model("Prior Mixture model.stan")
power_model <- stan_model("Power Posterior Mixture model.stan")

start <- proc.time()
set.seed(123)

# Data

data(galaxies, package = "MASS")
y <- galaxies / 1000
N <- length(y)
K <- 3

# Priors

alpha   <- rep(1, K)
mu0     <- mean(y)
lambda0 <- 2.6 / (max(y) - min(y))
a0      <- 1.28
b0      <- 0.36 * (mean(y^2) - mean(y)^2)

# Power posterior settings

T     <- 30
Nsim  <- 2000
c_pow <- 2

Beta_t <- function(t) (t / T)^c_pow
t_list <- Beta_t(0:T)

# Mixture log-likelihood

mix_loglik <- function(y, omega, mu, sigma2) {
  ll <- 0
  for (n in seq_along(y)) {
    log_terms <- log(omega) +
      dnorm(y[n], mu, sqrt(sigma2), log = TRUE)
    m <- max(log_terms)
    ll <- ll + m + log(sum(exp(log_terms - m)))
  }
  ll
}


#Sample from prior
prior_fit <- sampling(
  prior_model,
  data = list(
    K = K,
    alpha = alpha,
    mu0 = mu0,
    lambda0 = lambda0,
    a0 = a0,
    b0 = b0
  ),
  iter = Nsim,
  chains = 1,
  algorithm = "Fixed_param",
  refresh = 0
)

omega_samples  <- array(NA, c(T + 1, Nsim, K))
mu_samples     <- array(NA, c(T + 1, Nsim, K))
sigma2_samples <- array(NA, c(T + 1, Nsim, K))

omega_samples[1,,]  <- extract(prior_fit, "omega")$omega
mu_samples[1,,]     <- extract(prior_fit, "mu")$mu
sigma2_samples[1,,] <- extract(prior_fit, "sigma2")$sigma2


# Expectation storage

E <- matrix(0, nrow = T + 1, ncol = Nsim)

for (i in 1:Nsim) {
  E[1, i] <- mix_loglik(
    y,
    omega_samples[1, i, ],
    mu_samples[1, i, ],
    sigma2_samples[1, i, ]
  )
}


# Power posterior loop

for (t in 2:(T + 1)) {
  
  power <- t_list[t]
  cat(sprintf("Power posterior %d / %d   (t = %.3f)\n",
              t - 1, T, power))
  
  init_list <- list(list(
    omega  = omega_samples[t - 1, Nsim, ],
    mu     = mu_samples[t - 1, Nsim, ],
    sigma2 = sigma2_samples[t - 1, Nsim, ]
  ))
  
  fit <- sampling(
    power_model,
    data = list(
      N = N,
      K = K,
      y = y,
      alpha = alpha,
      mu0 = mu0,
      lambda0 = lambda0,
      a0 = a0,
      b0 = b0,
      t = power
    ),
    iter = Nsim,
    warmup = 0,
    chains = 1,
    init = init_list,
    refresh = 0
  )
  
  omega_samples[t,,]  <- extract(fit, "omega")$omega
  mu_samples[t,,]     <- extract(fit, "mu")$mu
  sigma2_samples[t,,] <- extract(fit, "sigma2")$sigma2
  
  for (i in 1:Nsim) {
    E[t, i] <- mix_loglik(
      y,
      omega_samples[t, i, ],
      mu_samples[t, i, ],
      sigma2_samples[t, i, ]
    )
  }
}


# Power posterior estimator

logZ <- 0
for (t in 2:(T + 1)) {
  logZ <- logZ +
    (t_list[t] - t_list[t - 1]) *
    0.5 * (mean(E[t, ]) + mean(E[t - 1, ]))
}

cat("\nEstimated log marginal likelihood:", logZ, "\n")

end <- proc.time()
(end - start)[3]
