start <- proc.time()

set.seed(123)


library(extraDistr)

data(galaxies, package = "MASS")

K <- 3
y <- galaxies / 1000
N <- length(y)

# Priors
alpha <- rep(1, K)
mu0 <- mean(y)
lambda0 <- 2.6 / (max(y) - min(y))
a0 <- 1.28
b0 <- 0.36 * (mean(y^2) - mean(y)^2)


# Log posterior (joint)

l_theta <- function(theta, y, K,
                    alpha, mu0, lambda0, a0, b0) {
  
  eta <- c(theta[1:(K-1)], 0)
  omega <- exp(eta) / sum(exp(eta))
  
  mu <- theta[K:(2*K-1)]
  log_sigma2 <- theta[(2*K):(3*K-1)]
  sigma2 <- exp(log_sigma2)
  
  # log likelihood
  ll <- sum(sapply(y, function(yi) {
    m <- max(log(omega) + dnorm(yi, mu, sqrt(sigma2), log = TRUE))
    m + log(sum(exp(
      log(omega) + dnorm(yi, mu, sqrt(sigma2), log = TRUE) - m
    )))
  }))
  
  # log prior
  lp <- ddirichlet(omega, alpha, log = TRUE) +
    sum(dnorm(mu, mu0, sqrt(1 / lambda0), log = TRUE)) +
    sum(dinvgamma(sigma2, a0, b0, log = TRUE) + log_sigma2)
  
  ll + lp
}

make_theta_init_prior <- function(K, alpha, mu0, lambda0, a0, b0) {
  
  # mixture weights
  omega <- as.numeric(rdirichlet(1, alpha))
  eta <- log(omega[-K] / omega[K])
  
  # component means
  mu <- rnorm(K, mu0, sqrt(1 / lambda0))
  
  # component variances
  sigma2 <- rinvgamma(K, a0, b0)
  
  c(eta, mu, log(sigma2))
}


# Optimisation (θ*)
#Options showing estimate depend on which mode found

#theta_init <- c(
#  rep(0, K - 1),
#  quantile(y, probs = seq(0.2, 0.8, length.out = K)),
#  rep(log(var(y)), K)
#)

theta_init <- make_theta_init_prior(
  K = K,
  alpha = alpha,
  mu0 = mu0,
  lambda0 = lambda0,
  a0 = a0,
  b0 = b0
)


optimiser <- optim(
  par = theta_init,
  fn = l_theta,
  y = y,
  K = K,
  alpha = alpha,
  mu0 = mu0,
  lambda0 = lambda0,
  a0 = a0,
  b0 = b0,
  control = list(fnscale = -1)
)

theta_star <- optimiser$par


# Gibbs sampler (stores z)

gibbs_mix <- function(S, y, K, alpha, mu0, lambda0, a0, b0) {
  
  n <- length(y)
  z <- sample(1:K, n, replace = TRUE)
  mu <- rep(mean(y), K)
  sigma2 <- rep(var(y), K)
  omega <- rep(1 / K, K)
  
  out <- list(
    theta = matrix(NA, S, 3*K - 1),
    z = matrix(NA, S, n)
  )
  
  for (s in 1:S) {
    
    # z | rest
    for (i in 1:n) {
      p <- omega * dnorm(y[i], mu, sqrt(sigma2))
      z[i] <- sample(1:K, 1, prob = p)
    }
    
    # mu_k | rest
    for (k in 1:K) {
      yk <- y[z == k]
      nk <- length(yk)
      lambda_n <- lambda0 + nk
      mu_n <- (lambda0 * mu0 + sum(yk)) / lambda_n
      mu[k] <- rnorm(1, mu_n, sqrt(sigma2[k] / lambda_n))
    }
    
    # sigma2_k | rest
    for (k in 1:K) {
      yk <- y[z == k]
      nk <- length(yk)
      a_n <- a0 + nk / 2
      b_n <- b0 + 0.5 * sum((yk - mu[k])^2)
      sigma2[k] <- rinvgamma(1, a_n, b_n)
    }
    
    # omega | rest
    counts <- tabulate(z, K)
    omega <- as.numeric(rdirichlet(1, alpha + counts))
    
    eta <- log(omega[-K] / omega[K])
    
    out$theta[s, ] <- c(eta, mu, log(sigma2))
    out$z[s, ] <- z
  }
  
  out
}

set.seed(456)
gibbs_out <- gibbs_mix(
  S = 5000,
  y = y,
  K = K,
  alpha = alpha,
  mu0 = mu0,
  lambda0 = lambda0,
  a0 = a0,
  b0 = b0
)


# Posterior ordinate at θ* (Chib)

eta_star <- c(theta_star[1:(K-1)], 0)
omega_star <- exp(eta_star) / sum(exp(eta_star))
mu_star <- theta_star[K:(2*K-1)]
sigma2_star <- exp(theta_star[(2*K):(3*K-1)])

log_post_vals <- numeric(nrow(gibbs_out$z))

for (s in 1:nrow(gibbs_out$z)) {
  
  z <- gibbs_out$z[s, ]
  
  # omega | z
  log_p_omega <- ddirichlet(
    omega_star,
    alpha + tabulate(z, K),
    log = TRUE
  )
  
  # mu_k | sigma2_k, z
  log_p_mu <- 0
  for (k in 1:K) {
    yk <- y[z == k]
    nk <- length(yk)
    lambda_n <- lambda0 + nk
    mu_n <- (lambda0 * mu0 + sum(yk)) / lambda_n
    log_p_mu <- log_p_mu +
      dnorm(mu_star[k], mu_n, sqrt(sigma2_star[k] / lambda_n), log = TRUE)
  }
  
  # sigma2_k | mu_k, z
  log_p_sigma2 <- 0
  for (k in 1:K) {
    yk <- y[z == k]
    nk <- length(yk)
    a_n <- a0 + nk / 2
    b_n <- b0 + 0.5 * sum((yk - mu_star[k])^2)
    log_p_sigma2 <- log_p_sigma2 +
      dinvgamma(sigma2_star[k], a_n, b_n, log = TRUE)
  }
  
  log_post_vals[s] <- log_p_omega + log_p_mu + log_p_sigma2
}

# log-sum-exp
m <- max(log_post_vals)
log_post_mc <- m + log(mean(exp(log_post_vals - m)))


# Prior + likelihood at θ*

lprior <- ddirichlet(omega_star, alpha, log = TRUE) +
  sum(dnorm(mu_star, mu0, sqrt(1 / lambda0), log = TRUE)) +
  sum(dinvgamma(sigma2_star, a0, b0, log = TRUE) +
        log(sigma2_star))

llike <- sum(sapply(y, function(yi) {
  m <- max(log(omega_star) + dnorm(yi, mu_star, sqrt(sigma2_star), log = TRUE))
  m + log(sum(exp(
    log(omega_star) +
      dnorm(yi, mu_star, sqrt(sigma2_star), log = TRUE) - m
  )))
}))


# Chib's estimator

chibs <- lprior + llike - log_post_mc
chibs

end <- proc.time()
(end - start)[3]


#theta_init changes the estimation as finds different posterior modes 
#due to label switching

#Theta_init1 gives -2 and theta_init2 gives -190
#Effect of label swithcing/multiple modes


#Chib’s estimator assumes:
#  “There is a single dominant posterior mode, and I can evaluate the posterior density at that point.”
#Finite mixtures violate this assumption by symmetry.
