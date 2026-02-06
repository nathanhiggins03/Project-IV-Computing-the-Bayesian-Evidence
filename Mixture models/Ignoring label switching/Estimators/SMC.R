setwd("~/Desktop/Project IV")

library(rstan)
library(MASS)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())


# Compile prior model

prior_model <- stan_model("Prior Mixture model.stan")

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


# SMC parameters

T    <- 10
Nsim <- 10000

Beta_t <- function(t, c = 4) (t / T)^c
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


# Sample from prior

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


# Log incremental weights

log_w <- matrix(0, nrow = T, ncol = Nsim)


# SMC loop

for (t in 2:(T + 1)) {
  
  t_prev <- t_list[t - 1]
  t_curr <- t_list[t]
  
  cat(sprintf("SMC step %d / %d   (%.4f → %.4f)\n",
              t - 1, T, t_prev, t_curr))
  
  # Weight update
  for (i in 1:Nsim) {
    ll <- mix_loglik(
      y,
      omega_samples[t - 1, i, ],
      mu_samples[t - 1, i, ],
      sigma2_samples[t - 1, i, ]
    )
    log_w[t - 1, i] <- (t_curr - t_prev) * ll
  }
  
  # Normalise weights (log-sum-exp)
  max_logw <- max(log_w[t - 1, ])
  w <- exp(log_w[t - 1, ] - max_logw)
  
  # Resample
  idx <- sample(1:Nsim, size = Nsim, replace = TRUE, prob = w)
  
  omega_samples[t,,]  <- omega_samples[t - 1, idx, ]
  mu_samples[t,,]     <- mu_samples[t - 1, idx, ]
  sigma2_samples[t,,] <- sigma2_samples[t - 1, idx, ]
}

# SMC estimate
logZ <- 0
for (t in 1:T) {
  max_logw <- max(log_w[t, ])
  logZ <- logZ + log(mean(exp(log_w[t, ] - max_logw))) + max_logw
}

cat("\nEstimated log marginal likelihood:", logZ, "\n")

end <- proc.time()
(end - start)[3]
