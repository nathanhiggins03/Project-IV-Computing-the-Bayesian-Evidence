#Prior method

library(ggplot2)

Sim <- 30
MC_values <- c(3000000)   # <-- choose MC sizes

df_all <- data.frame()

library(rstan)
library(durhamSLR)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

setwd("~/Desktop/Project IV")

prior_model <- stan_model(
  file = "Prior Mixture model.stan"
)

for (MC_sample in MC_values) {
  
  time_simulation <- rep(0, Sim)
  est_simulation  <- rep(0, Sim)
  
  for (k in 1:Sim) {
    
    start=proc.time()
    #Set seed for reproducibility
    set.seed(123+k)
    
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
      iter = 2*MC_sample,
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
    est_simulation[k]<-stan_prior_le
    
    end=proc.time()
    
    timer<- end-start
    time_simulation[k]<- timer[3]
    
    #Can compare prior method to "arithmetic mean" in paper for each K, know we are
    #reproducing correct results
    
    
  }
  
  #Store results
  df_all <- rbind(
    df_all,
    data.frame(
      Estimate = est_simulation,
      Time     = time_simulation,   # ← ADD THIS LINE
      MC = factor(paste0("N = ", MC_sample),
                  levels = paste0("N = ", MC_values))
    )
  )
  
}

df_prior <- df_all   # for Prior MC



#HME

#Need to change to incorporate more efficient code for HME
library(ggplot2)

Sim <- 30
MC_values <- c(100000)    # <-- choose MC sizes

df_all <- data.frame()

library(rstan)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

setwd("~/Desktop/Project IV")

posterior_model <- stan_model(
  file = "Posterior Mixture model.stan"
)


for (MC_sample in MC_values) {
  
  time_simulation <- rep(0, Sim)
  est_simulation  <- rep(0, Sim)
  
  for (k in 1:Sim) {
    
    start=proc.time()
    
    #Bayesian Linear Regression with unknown beta and precision
    
    #Set seed for reproducibility
    set.seed(123+k)
    
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
      iter = 2*MC_sample,
      chains = 1,
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
    est_simulation[k]<-hme_log
    
    time_simulation[k] <- (proc.time() - start)[3]
    
  }
  
  #Store results
  df_all <- rbind(
    df_all,
    data.frame(
      Estimate = est_simulation,
      Time     = time_simulation,   # ← ADD THIS LINE
      MC = factor(paste0("N = ", MC_sample),
                  levels = paste0("N = ", MC_values))
    )
  )
  
}

df_hme <- df_all     # for HME



#Chib

library(extraDistr)
library(ggplot2)
library(dplyr)

data(galaxies, package = "MASS")
y <- galaxies / 1000
K <- 3
alpha <- rep(1, K)
mu0 <- mean(y)
lambda0 <- 2.6 / (max(y) - min(y))
a0 <- 1.28
b0 <- 0.36 * (mean(y^2) - mean(y)^2)
MC_sample <- 100000  # Gibbs sample size

# Functions
logsumexp <- function(x) { m <- max(x); m + log(sum(exp(x - m))) }

l_theta <- function(theta, y, K, alpha, mu0, lambda0, a0, b0) {
  eta <- c(theta[1:(K-1)], 0)
  omega <- exp(eta) / sum(exp(eta))
  mu <- theta[K:(2*K-1)]
  sigma2 <- exp(theta[(2*K):(3*K-1)])
  
  ll <- sum(sapply(y, function(yi) {
    m <- max(log(omega) + dnorm(yi, mu, sqrt(sigma2), log=TRUE))
    m + log(sum(exp(log(omega) + dnorm(yi, mu, sqrt(sigma2), log=TRUE) - m)))
  }))
  
  lp <- ddirichlet(omega, alpha, log=TRUE) +
    sum(dnorm(mu, mu0, sqrt(1/lambda0), log=TRUE)) +
    sum(dinvgamma(sigma2, a0, b0, log=TRUE) + log(sigma2))
  
  ll + lp
}

make_theta_init_prior <- function(K, alpha, mu0, lambda0, a0, b0) {
  omega <- as.numeric(rdirichlet(1, alpha))
  eta <- log(omega[-K] / omega[K])
  mu <- rnorm(K, mu0, sqrt(1/lambda0))
  sigma2 <- rinvgamma(K, a0, b0)
  c(eta, mu, log(sigma2))
}

gibbs_mix <- function(S, y, K, alpha, mu0, lambda0, a0, b0) {
  n <- length(y)
  z <- sample(1:K, n, replace=TRUE)
  mu <- rep(mean(y), K)
  sigma2 <- rep(var(y), K)
  omega <- rep(1/K, K)
  
  theta_out <- matrix(NA, S, 3*K-1)
  z_out <- matrix(NA, S, n)
  
  for (s in 1:S) {
    for (i in 1:n) {
      p <- omega * dnorm(y[i], mu, sqrt(sigma2))
      z[i] <- sample(1:K, 1, prob=p)
    }
    for (k in 1:K) {
      yk <- y[z==k]; nk <- length(yk)
      lambda_n <- lambda0 + nk
      mu_n <- (lambda0*mu0 + sum(yk))/lambda_n
      mu[k] <- rnorm(1, mu_n, sqrt(sigma2[k]/lambda_n))
    }
    for (k in 1:K) {
      yk <- y[z==k]; nk <- length(yk)
      a_n <- a0 + nk/2
      b_n <- b0 + 0.5*sum((yk - mu[k])^2)
      sigma2[k] <- rinvgamma(1, a_n, b_n)
    }
    counts <- tabulate(z, K)
    omega <- as.numeric(rdirichlet(1, alpha + counts))
    theta_out[s,] <- c(log(omega[-K]/omega[K]), mu, log(sigma2))
    z_out[s,] <- z
  }
  
  list(theta=theta_out, z=z_out)
}

# Wrap the full Chib procedure in a function
run_chib <- function() {
  theta_init <- make_theta_init_prior(K, alpha, mu0, lambda0, a0, b0)
  optimiser <- optim(theta_init, l_theta, y=y, K=K, alpha=alpha,
                     mu0=mu0, lambda0=lambda0, a0=a0, b0=b0,
                     control=list(fnscale=-1))
  theta_star <- optimiser$par
  
  gibbs_out <- gibbs_mix(MC_sample, y, K, alpha, mu0, lambda0, a0, b0)
  
  eta_star <- c(theta_star[1:(K-1)], 0)
  omega_star <- exp(eta_star) / sum(exp(eta_star))
  mu_star <- theta_star[K:(2*K-1)]
  sigma2_star <- exp(theta_star[(2*K):(3*K-1)])
  
  log_post_vals <- numeric(nrow(gibbs_out$z))
  for (s in 1:nrow(gibbs_out$z)) {
    z <- gibbs_out$z[s,]
    log_post_vals[s] <- ddirichlet(omega_star, alpha + tabulate(z,K), log=TRUE) +
      sum(sapply(1:K, function(k) {
        yk <- y[z==k]; nk <- length(yk)
        lambda_n <- lambda0 + nk; mu_n <- (lambda0*mu0 + sum(yk))/lambda_n
        log_p_mu <- dnorm(mu_star[k], mu_n, sqrt(sigma2_star[k]/lambda_n), log=TRUE)
        a_n <- a0 + nk/2; b_n <- b0 + 0.5*sum((yk - mu_star[k])^2)
        log_p_sigma <- dinvgamma(sigma2_star[k], a_n, b_n, log=TRUE)
        log_p_mu + log_p_sigma
      }))
  }
  
  log_post_mc <- logsumexp(log_post_vals) - log(length(log_post_vals))
  
  lprior <- ddirichlet(omega_star, alpha, log=TRUE) +
    sum(dnorm(mu_star, mu0, sqrt(1/lambda0), log=TRUE)) +
    sum(dinvgamma(sigma2_star, a0, b0, log=TRUE) + log(sigma2_star))
  
  llike <- sum(sapply(y, function(yi) {
    m <- max(log(omega_star) + dnorm(yi, mu_star, sqrt(sigma2_star), log=TRUE))
    m + log(sum(exp(log(omega_star) + dnorm(yi, mu_star, sqrt(sigma2_star), log=TRUE) - m)))
  }))
  
  chibs <- lprior + llike - log_post_mc
  chibs
}

# Run 30 simulations

pt_start <- proc.time()

chib_estimates <- replicate(30, run_chib())

pt_end <- proc.time()

# Elapsed time in seconds
elapsed_time <- (pt_end - pt_start)["elapsed"]

# Mean runtime per Chib run
mean_runtime <- elapsed_time / 30

df_chib <- data.frame(Estimate = chib_estimates,
                      Time     = rep(mean_runtime,30),   # ← ADD THIS LINE
                      MC = factor(paste0("N = ", 0),
                                  levels = paste0("N = ", 0)))









#Laplace

############################
## Laplace approximation ##
## Multiple initialisations
############################

library(extraDistr)
library(mvtnorm)
library(ggplot2)
library(dplyr)
library(ggbeeswarm)

data(galaxies, package = "MASS")
y <- galaxies / 1000
K <- 3
N <- length(y)

# Priors
alpha <- rep(1, K)
mu0 <- mean(y)
lambda0 <- 2.6 / (max(y) - min(y))
a0 <- 1.28
b0 <- 0.36 * (mean(y^2) - mean(y)^2)

############################
# Log posterior
############################
l_theta <- function(theta, y, K, alpha, mu0, lambda0, a0, b0) {
  
  eta <- c(theta[1:(K-1)], 0)
  omega <- exp(eta) / sum(exp(eta))
  mu <- theta[K:(2*K-1)]
  log_sigma2 <- theta[(2*K):(3*K-1)]
  sigma2 <- exp(log_sigma2)
  
  ll <- sum(sapply(y, function(yi) {
    m <- max(log(omega) + dnorm(yi, mu, sqrt(sigma2), log = TRUE))
    m + log(sum(exp(
      log(omega) + dnorm(yi, mu, sqrt(sigma2), log = TRUE) - m
    )))
  }))
  
  lp <- ddirichlet(omega, alpha, log = TRUE) +
    sum(dnorm(mu, mu0, sqrt(1 / lambda0), log = TRUE)) +
    sum(dinvgamma(sigma2, a0, b0, log = TRUE) + log_sigma2)
  
  ll + lp
}

############################
# Random initialisation - two options
############################
make_theta_init_prior <- function() {
  omega <- as.numeric(rdirichlet(1, alpha))
  eta <- log(omega[-K] / omega[K])
  mu <- rnorm(K, mu0, sqrt(1 / lambda0))
  sigma2 <- rinvgamma(K, a0, b0)
  c(eta, mu, log(sigma2))
}

make_theta_init_wild <- function(K, y) {
  c(
    rnorm(K - 1, 0, 3),                 # eta (weights)
    runif(K, min(y) - 5, max(y) + 5),   # mu (very wide)
    rnorm(K, 0, 3)                      # log sigma2
  )
}

############################
# Single Laplace run
############################
run_laplace <- function() {
  
  start <- proc.time()
  
  #theta_init <- make_theta_init_prior()
  theta_init <- make_theta_init_wild(K, y)
  
  
  opt <- optim(
    par = theta_init,
    fn = l_theta,
    y = y,
    K = K,
    alpha = alpha,
    mu0 = mu0,
    lambda0 = lambda0,
    a0 = a0,
    b0 = b0,
    hessian = TRUE,
    control = list(fnscale = -1)
  )
  
  H <- opt$hessian
  if (any(!is.finite(H)) || det(-H) <= 0) return(c(NA, NA))
  
  Sigma <- solve(-H)
  d <- length(opt$par)
  
  log_laplace <- (d / 2) * log(2 * pi) +
    0.5 * log(det(Sigma)) +
    l_theta(opt$par, y, K, alpha, mu0, lambda0, a0, b0)
  
  time <- (proc.time() - start)[3]
  
  c(log_laplace, time)
}

############################
# Run 30 times
############################
res <- replicate(30, run_laplace())
res <- t(res)
colnames(res) <- c("Estimate", "Time")
df <- as.data.frame(res)
df <- df %>% filter(is.finite(Estimate))


df_laplace <- data.frame(
  Estimate = df$Estimate,
  Time     = df$Time,
  MC       = factor(
    paste0("N = ", 0),
    levels = paste0("N = ", 0)
  ))




#AIS

#Need to change to incorporate more efficient code for HME
library(ggplot2)

Sim <- 30
MC_values <- c(12000)   # <-- choose MC sizes
T_value<-20
df_all <- data.frame()

library(rstan)


rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

setwd("~/Desktop/Project IV")

# Compile Stan models

prior_model <- stan_model("Prior Mixture model.stan")
pp_model <- stan_model("Power Posterior Mixture model.stan")


for (MC_sample in MC_values) {
  
  time_simulation <- rep(0, Sim)
  est_simulation  <- rep(0, Sim)
  
  for (k in 1:Sim) {
    
    start=proc.time()
    set.seed(123+k)
    
    # Data
    
    data(galaxies, package = "MASS")
    y <- galaxies / 1000
    N <- length(y)
    K <- 3
    
    # Priors
    
    alpha <- rep(1, K)
    mu0 <- mean(y)
    lambda0 <- 2.6 / (max(y) - min(y))
    a0 <- 1.28
    b0 <- 0.36 * (mean(y^2) - mean(y)^2)
    
    # AIS parameters
    
    T <- T_value
    Nsim <- MC_sample
    c_power <- 2
    
    Beta_t <- function(t) (t / T)^c_power
    t_list <- Beta_t(0:T)
    
    
    # 1. Sample from prior
    
    prior_fit <- sampling(
      prior_model,
      data = list(K = K, alpha = alpha, mu0 = mu0,
                  lambda0 = lambda0, a0 = a0, b0 = b0),
      iter = Nsim,
      chains = 1,
      algorithm = "Fixed_param",
      refresh = 0
    )
    
    omega_store  <- array(NA, c(T + 1, Nsim, K))
    mu_store     <- array(NA, c(T + 1, Nsim, K))
    sigma2_store <- array(NA, c(T + 1, Nsim, K))
    
    omega_store[1,,]  <- extract(prior_fit, "omega")$omega
    mu_store[1,,]     <- extract(prior_fit, "mu")$mu
    sigma2_store[1,,] <- extract(prior_fit, "sigma2")$sigma2
    
    log_w <- matrix(0, nrow = T + 1, ncol = Nsim)
    
    # Mixture log-likelihood
    
    mix_loglik <- function(y, omega, mu, sigma2) {
      ll <- 0
      for (n in seq_along(y)) {
        log_terms <- log(omega) +
          dnorm(y[n], mu, sqrt(sigma2), log = TRUE)
        ll <- ll + log(sum(exp(log_terms - max(log_terms)))) +
          max(log_terms)
      }
      ll
    }
    
    
    
    #AIS loop
    
    warmup <- 3000
    
    for (t_idx in 2:(T + 1)) {
      
      t_prev <- t_list[t_idx - 1]
      t_curr <- t_list[t_idx]
      
      cat(sprintf("AIS step %d / %d  (%.3f → %.3f)\n",
                  t_idx - 1, T, t_prev, t_curr))
      
      init <- list(list(
        omega  = omega_store[t_idx - 1, Nsim, ],
        mu     = mu_store[t_idx - 1, Nsim, ],
        sigma2 = sigma2_store[t_idx - 1, Nsim, ]
      ))
      
      fit <- sampling(
        pp_model,
        data = list(
          N = N, K = K, y = y,
          alpha = alpha,
          mu0 = mu0,
          lambda0 = lambda0,
          a0 = a0,
          b0 = b0,
          t = t_curr
        ),
        iter = warmup + Nsim,
        warmup = warmup,
        chains = 1,
        init = init,
        refresh = 0
      )
      
      omega_store[t_idx,,]  <- extract(fit, "omega")$omega
      mu_store[t_idx,,]     <- extract(fit, "mu")$mu
      sigma2_store[t_idx,,] <- extract(fit, "sigma2")$sigma2
      
      for (i in 1:Nsim) {
        ll <- mix_loglik(
          y,
          omega_store[t_idx - 1, i, ],
          mu_store[t_idx - 1, i, ],
          sigma2_store[t_idx - 1, i, ]
        )
        log_w[t_idx, i] <- log_w[t_idx - 1, i] +
          (t_curr - t_prev) * ll
      }
    }
    
    
    #Log marginal likelihood
    
    final_log_w <- log_w[T + 1, ]
    m <- max(final_log_w)
    logZ <- m + log(mean(exp(final_log_w - m)))
    
    cat("\nEstimated log marginal likelihood:", logZ, "\n")
    
    est_simulation[k]<-logZ
    
    time_simulation[k] <- (proc.time() - start)[3]
    
    
    #Important to note in report that we can't do trick in calculating weights
    #seen on page 18 as priors don't cancel in a mixture as it is a sum
    #hence have to caclculate weights using general formula in AIS algorithm
    #This is main calculation change for AIS
  }
  
  #Store results
  df_all <- rbind(
    df_all,
    data.frame(
      Estimate = est_simulation,
      Time     = time_simulation,   # ← ADD THIS LINE
      MC = factor(paste0("N = ", MC_sample),
                  levels = paste0("N = ", MC_values))
    )
  )
  
}



df_ais <- df_all     # for AIS


#Power posterior

#Need to change to incorporate more efficient code for HME
library(ggplot2)

Sim <- 30
MC_values <- c(12000)   # <-- choose MC sizes
T_value<-20

df_all <- data.frame()

setwd("~/Desktop/Project IV")

library(rstan)


rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())


# Compile Stan models

prior_model <- stan_model("Prior Mixture model.stan")
power_model <- stan_model("Power Posterior Mixture model.stan")


for (MC_sample in MC_values) {
  
  time_simulation <- rep(0, Sim)
  est_simulation  <- rep(0, Sim)
  
  for (k in 1:Sim) {
    
    start <- proc.time()
    set.seed(123+k)
    
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
    
    T     <- T_value
    Nsim  <- MC_sample
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
    
    est_simulation[k]<-logZ
    
    time_simulation[k] <- (proc.time() - start)[3]
    
  }
  
  #Store results
  df_all <- rbind(
    df_all,
    data.frame(
      Estimate = est_simulation,
      Time     = time_simulation,   # ← ADD THIS LINE
      MC = factor(paste0("N = ", MC_sample),
                  levels = paste0("N = ", MC_values))
    )
  )
  
}


df_pp <- df_all      # for TI



#SMC


#Need to change to incorporate more efficient code for HME
library(ggplot2)

Sim <- 30
MC_values <- c(12000)  # <-- choose MC sizes
T_value<-20
df_all <- data.frame()

setwd("~/Desktop/Project IV")

library(rstan)
library(MASS)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())


# Compile prior model

prior_model <- stan_model("Prior Mixture model.stan")

for (MC_sample in MC_values) {
  
  time_simulation <- rep(0, Sim)
  est_simulation  <- rep(0, Sim)
  
  for (k in 1:Sim) {
    
    start <- proc.time()
    set.seed(123+k)
    
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
    
    T    <- T_value
    Nsim <- MC_sample
    
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
    
    est_simulation[k]<-logZ
    
    time_simulation[k] <- (proc.time() - start)[3]
    
    
    
  }
  
  #Store results
  df_all <- rbind(
    df_all,
    data.frame(
      Estimate = est_simulation,
      Time     = time_simulation,   # ← ADD THIS LINE
      MC = factor(paste0("N = ", MC_sample),
                  levels = paste0("N = ", MC_values))
    )
  )
  
}


df_smc <- df_all     # for SMC


library(dplyr)

df_prior <- df_prior %>% mutate(Estimator = "Prior")
df_hme   <- df_hme   %>% mutate(Estimator = "HME")
df_chib  <- df_chib  %>% mutate(Estimator = "Chib")
df_laplace <- df_laplace  %>% mutate(Estimator = "Laplace")
df_ais   <- df_ais   %>% mutate(Estimator = "AIS")
df_pp    <- df_pp    %>% mutate(Estimator = "Power Posterior")
df_smc   <- df_smc   %>% mutate(Estimator = "SMC")


df_compare <- bind_rows(
  df_prior,
  df_hme,
  df_ais,
  df_pp,
  df_smc,
  df_chib,
  df_laplace
)

df_compare$Estimator <- factor(
  df_compare$Estimator,
  levels = c("Prior", "HME", "AIS", "Power Posterior", "SMC","Chib","Laplace")
)

library(ggplot2)

ggplot(df_compare, aes(x = Estimator, y = Estimate, fill = Estimator)) +
  geom_violin(trim = FALSE, alpha = 0.7) +
  geom_boxplot(width = 0.12,
               fill = "white",
               outlier.shape = NA) +
  labs(
    title = "Comparing Evidence estimators",
    y = "Log Evidence",
    x = ""
  ) +
  theme_minimal() +
  theme(legend.position = "none")




df_efficiency <- df_compare %>%
  group_by(Estimator) %>%
  summarise(
    mean_time = mean(Time),
    mc_sd     = sd(Estimate),
    .groups = "drop"
  )

ggplot(df_efficiency,
       aes(x = mean_time, y = mc_sd, colour = Estimator)) +
  geom_point(size = 4) +
  labs(
    title = "Efficiency comparison",
    x = "Mean runtime (seconds)",
    y = "Monte Carlo SD",
    colour = "Estimator"
  ) +
  theme_minimal()
