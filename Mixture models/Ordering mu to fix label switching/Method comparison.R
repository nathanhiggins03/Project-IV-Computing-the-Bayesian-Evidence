#Prior method

library(ggplot2)

Sim <- 30
MC_values <- c(100000)   # <-- choose MC sizes

df_all <- data.frame()

library(rstan)
library(durhamSLR)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

setwd("~/Desktop/Project IV")

prior_model <- stan_model(
  file = "Ordered Prior Mixture model.stan"
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
MC_values <- c(10000)    # <-- choose MC sizes

df_all <- data.frame()

library(rstan)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

setwd("~/Desktop/Project IV")

posterior_model <- stan_model(
  file = "Ordered Posterior Mixture model.stan"
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
MC_sample <- 10000  # Gibbs sample size

# Functions
logsumexp <- function(x) { m <- max(x); m + log(sum(exp(x - m))) }

relabel_by_mu <- function(z, mu, sigma2, omega) {
  ord <- order(mu)
  
  mu_new     <- mu[ord]
  sigma2_new <- sigma2[ord]
  omega_new  <- omega[ord]
  
  # relabel allocations: old label -> new ordered label
  z_new <- match(z, ord)
  
  list(
    z      = z_new,
    mu     = mu_new,
    sigma2 = sigma2_new,
    omega  = omega_new
  )
}

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
    
    ## 5. ***RELABEL HERE***
    tmp <- relabel_by_mu(z, mu, sigma2, omega)
    z      <- tmp$z
    mu     <- tmp$mu
    sigma2 <- tmp$sigma2
    omega  <- tmp$omega
    
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
  
  # RELABEL theta_star to match Gibbs ordering
  ord <- order(mu_star)
  
  mu_star     <- mu_star[ord]
  sigma2_star <- sigma2_star[ord]
  omega_star  <- omega_star[ord]
  
  
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
  
  ## --- mixture weights ---
  eta <- c(theta[1:(K-1)], 0)
  omega <- exp(eta) / sum(exp(eta))
  
  ## --- ordered means ---
  mu_raw <- theta[K:(2*K-1)]
  
  mu <- numeric(K)
  mu[1] <- mu_raw[1]
  for (k in 2:K) {
    mu[k] <- mu[k-1] + exp(mu_raw[k])
  }
  
  ## --- variances ---
  log_sigma2 <- theta[(2*K):(3*K-1)]
  sigma2 <- exp(log_sigma2)
  
  ## --- log-likelihood ---
  ll <- sum(sapply(y, function(yi) {
    m <- max(log(omega) + dnorm(yi, mu, sqrt(sigma2), log = TRUE))
    m + log(sum(exp(
      log(omega) + dnorm(yi, mu, sqrt(sigma2), log = TRUE) - m
    )))
  }))
  
  ## --- log-prior ---
  lp <- ddirichlet(omega, alpha, log = TRUE) +
    sum(dnorm(mu, mu0, sqrt(1 / lambda0), log = TRUE)) +
    sum(dinvgamma(sigma2, a0, b0, log = TRUE) + log_sigma2) +
    sum(mu_raw[2:K])  # Jacobian of exp()
  
  ll + lp
}


############################
# Random initialisation - two options
############################
make_theta_init_prior <- function() {
  
  ## weights
  omega <- as.numeric(rdirichlet(1, alpha))
  eta <- log(omega[-K] / omega[K])
  
  ## ordered means
  mu <- sort(rnorm(K, mu0, sqrt(1 / lambda0)))
  
  mu_raw <- numeric(K)
  mu_raw[1] <- mu[1]
  for (k in 2:K) {
    mu_raw[k] <- log(mu[k] - mu[k-1])
  }
  
  ## variances
  sigma2 <- rinvgamma(K, a0, b0)
  
  c(eta, mu_raw, log(sigma2))
}

make_theta_init_wild <- function(K, y) {
  
  mu <- sort(runif(K, min(y) - 5, max(y) + 5))
  
  mu_raw <- numeric(K)
  mu_raw[1] <- mu[1]
  for (k in 2:K) {
    mu_raw[k] <- log(mu[k] - mu[k-1])
  }
  
  c(
    rnorm(K - 1, 0, 3),   # eta
    mu_raw,               # ordered means
    rnorm(K, 0, 3)        # log sigma2
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
MC_values <- c(2000)   # <-- choose MC sizes
T_value<-10
df_all <- data.frame()

library(rstan)


rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

setwd("~/Desktop/Project IV")

# Compile Stan models

prior_model <- stan_model("Ordered Prior Mixture model.stan")
pp_model <- stan_model("Ordered Power Posterior Mixture model.stan")


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
      iter = 2*Nsim,
      chains = 1,
      algorithm = "Fixed_param",
      refresh = 0
    )
    
    #omega_store  <- array(NA, c(T + 1, Nsim, K))
    #mu_store     <- array(NA, c(T + 1, Nsim, K))
    #sigma2_store <- array(NA, c(T + 1, Nsim, K))
    
    #omega_store[1,,]  <- extract(prior_fit, "omega")$omega
    #mu_store[1,,]     <- extract(prior_fit, "mu")$mu
    #sigma2_store[1,,] <- extract(prior_fit, "sigma2")$sigma2
    
    
    omega_raw  <- extract(prior_fit, "omega")$omega   # [Nsim, K]
    mu_raw     <- extract(prior_fit, "mu")$mu        # [Nsim, K]
    sigma2_raw <- extract(prior_fit, "sigma2")$sigma2 # [Nsim, K]
    
    omega_store  <- array(NA, c(T + 1, Nsim, K))
    mu_store     <- array(NA, c(T + 1, Nsim, K))
    sigma2_store <- array(NA, c(T + 1, Nsim, K))
    
    for (i in 1:Nsim) {
      # Get permutation that sorts mu ascending
      idx <- order(mu_raw[i, ])
      
      mu_store[1, i, ]     <- mu_raw[i, idx]
      sigma2_store[1, i, ] <- sigma2_raw[i, idx]
      omega_store[1, i, ]  <- omega_raw[i, idx]
    }
    
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
MC_values <- c(2000)   # <-- choose MC sizes
T_value<-10

df_all <- data.frame()

setwd("~/Desktop/Project IV")

library(rstan)


rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())


# Compile Stan models

prior_model <- stan_model("Ordered Prior Mixture model.stan")
power_model <- stan_model("Ordered Power Posterior Mixture model.stan")


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
      iter = 2*Nsim,
      chains = 1,
      algorithm = "Fixed_param",
      refresh = 0
    )
    
    #omega_samples  <- array(NA, c(T + 1, Nsim, K))
    #mu_samples     <- array(NA, c(T + 1, Nsim, K))
    #sigma2_samples <- array(NA, c(T + 1, Nsim, K))
    
    #omega_samples[1,,]  <- extract(prior_fit, "omega")$omega
    #mu_samples[1,,]     <- extract(prior_fit, "mu")$mu
    #sigma2_samples[1,,] <- extract(prior_fit, "sigma2")$sigma2
    
    
    
    omega_raw  <- extract(prior_fit, "omega")$omega   # [Nsim, K]
    mu_raw     <- extract(prior_fit, "mu")$mu        # [Nsim, K]
    sigma2_raw <- extract(prior_fit, "sigma2")$sigma2 # [Nsim, K]
    
    
    omega_samples  <- array(NA, c(T + 1, Nsim, K))
    mu_samples     <- array(NA, c(T + 1, Nsim, K))
    sigma2_samples <- array(NA, c(T + 1, Nsim, K))
    
    for (i in 1:Nsim) {
      # Get permutation that sorts mu ascending
      idx <- order(mu_raw[i, ])
      
      mu_samples[1, i, ]     <- mu_raw[i, idx]
      sigma2_samples[1, i, ] <- sigma2_raw[i, idx]
      omega_samples[1, i, ]  <- omega_raw[i, idx]
    }
    
    
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
MC_values <- c(20000)  # <-- choose MC sizes
T_value<-20
df_all <- data.frame()

setwd("~/Desktop/Project IV")

library(rstan)
library(MASS)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())


# Compile prior model

prior_model <- stan_model("Ordered Prior Mixture model.stan")

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

#library(ggplot2)
#Violin plots
#ggplot(df_compare, aes(x = Estimator, y = Estimate, fill = Estimator)) +
#  geom_violin(trim = FALSE, alpha = 0.7) +
#  geom_boxplot(width = 0.12,
#               fill = "white",
#               outlier.shape = NA) +
#  labs(
#    title = "Comparing Evidence estimators",
#    y = "Log Evidence",
#    x = ""
#  ) +
#  theme_minimal() +
#  theme(legend.position = "none")




df_efficiency <- df_compare %>%
  group_by(Estimator) %>%
  summarise(
    mean_time = mean(Time),
    mc_sd     = sd(Estimate),
    mc_mean = mean(Estimate),
    .groups = "drop"
  )

#Efficiency : SD vs mean run time
ggplot(df_efficiency,
       aes(x = mean_time, y = mc_sd, colour = Estimator)) +
  geom_point(size = 4) +
  labs(
    title = "Efficiency comparison",
    x = "Mean runtime (seconds)",
    y = "Monte Carlo SD",
    colour = "Estimator"
  ) +
  theme_minimal()+
  theme(
    legend.title = element_text(size = 14),
    legend.text  = element_text(size = 12)
  )

ggplot(df_compare,
       aes(x = Estimate, y = Estimator, fill = Estimator)) +
  geom_boxplot(alpha = 0.7) +
  labs(
    title = "Comparing evidence estimators",
    x = "Log Evidence",
    y = ""
  ) +
  theme_minimal() +
  theme(legend.position = "none")


ggplot(df_compare,
       aes(x = Estimate, y = Estimator, fill = Estimator)) +
  geom_boxplot(alpha = 0.7) +
  labs(
    title = "Comparing evidence estimators",
    x = "Log Evidence",
    y = ""
  ) +
  theme_minimal() +
  coord_cartesian(xlim = c(NA, -200))+
  theme(legend.position = "none")





# Print mean and SD for each estimator
df_efficiency %>%
  select(Estimator, mc_mean, mc_sd) %>%
  arrange(factor(Estimator, levels = c("Prior", "HME", "AIS", "Power Posterior", "SMC","Chib","Laplace")))

#Plots without Chib
df_compare_nochib <- df_compare %>%
  filter(Estimator != "Chib") %>%
  mutate(
    Estimator = factor(
      Estimator,
      levels = c("Prior", "HME", "AIS", "Power Posterior", "SMC", "Laplace")
    )
  )

#ggplot(df_compare_nochib,
#       aes(x = Estimator, y = Estimate, fill = Estimator)) +
#  geom_violin(trim = FALSE, alpha = 0.7) +
#  geom_boxplot(
#    width = 0.12,
#    fill = "white",
#    outlier.shape = NA
#  ) +
#  labs(
#    title = "Comparing Evidence Estimators",
#    y = "Log Evidence",
#    x = ""
#  ) +
#  theme_minimal() +
#  theme(legend.position = "none")

ggplot(df_compare_nochib,
       aes(x = Estimate, y = Estimator, fill = Estimator)) +
  geom_boxplot(alpha = 0.7) +
  labs(
    title = "Comparing evidence estimators",
    x = "Log Evidence",
    y = ""
  ) +
  theme_minimal() +
  theme(legend.position = "none")

df_efficiency_nochib <- df_compare_nochib %>%
  group_by(Estimator) %>%
  summarise(
    mean_time = mean(Time),
    mc_sd     = sd(Estimate),
    mc_mean   = mean(Estimate),
    .groups   = "drop"
  )

ggplot(df_efficiency_nochib,
       aes(x = mean_time, y = mc_sd, colour = Estimator)) +
  geom_point(size = 4) +
  labs(
    title = "Efficiency comparison",
    x = "Mean runtime (seconds)",
    y = "Monte Carlo SD",
    colour = "Estimator"
  ) +
  theme_minimal()+
  theme(
    legend.title = element_text(size = 14),
    legend.text  = element_text(size = 12)
  )


#Plots without Chib and Laplace

estimators_keep <- c("Prior", "HME", "AIS", "Power Posterior", "SMC")

df_compare_nolaplace <- df_compare %>%
  filter(Estimator %in% estimators_keep) %>%
  mutate(
    Estimator = factor(Estimator, levels = estimators_keep)
  )

ggplot(df_compare_nolaplace,
       aes(x = Estimator, y = Estimate, fill = Estimator)) +
  geom_violin(trim = FALSE, alpha = 0.7) +
  geom_boxplot(
    width = 0.12,
    fill = "white",
    outlier.shape = NA
  ) +
  labs(
    title = "Comparing evidence estimators",
    y = "Log Evidence",
    x = ""
  ) +
  theme_minimal() +
  theme(legend.position = "none")

df_efficiency_nolaplace <- df_compare_nolaplace %>%
  group_by(Estimator) %>%
  summarise(
    mean_time = mean(Time),
    mc_sd     = sd(Estimate),
    mc_mean   = mean(Estimate),
    .groups   = "drop"
  )

ggplot(df_efficiency_nolaplace,
       aes(x = mean_time, y = mc_sd, colour = Estimator)) +
  geom_point(size = 4) +
  labs(
    title = "Efficiency comparison",
    x = "Mean runtime (seconds)",
    y = "Monte Carlo SD",
    colour = "Estimator"
  ) +
  theme_minimal()


#Boxplots
#ggplot(df_compare_nolaplace,
#       aes(x = Estimator, y = Estimate, fill = Estimator)) +
#  geom_boxplot(
#    alpha = 0.7,
#    width = 0.6,
#    outlier.shape = 16
#  ) +
#  labs(
#    title = "Comparing Evidence Estimators",
#    y = "Log Evidence",
#    x = ""
#  ) +
#  theme_minimal() +
#  theme(legend.position = "none")


ggplot(df_compare_nolaplace,
       aes(x = Estimate, y = Estimator, fill = Estimator)) +
  geom_boxplot(alpha = 0.7) +
  labs(
    title = "Comparing evidence estimators",
    x = "Log Evidence",
    y = ""
  ) +
  theme_minimal() +
  theme(legend.position = "none")





#Plots with correct colour for without laplace and chib

all_levels <- c(
  "Prior", "HME", "AIS",
  "Power Posterior", "SMC",
  "Chib", "Laplace"
)

full_pal <- scales::hue_pal()(length(all_levels))
names(full_pal) <- all_levels


ggplot(df_compare_nolaplace,
       aes(x = Estimator, y = Estimate, fill = Estimator)) +
  geom_violin(trim = FALSE, alpha = 0.7) +
  geom_boxplot(
    width = 0.12,
    fill = "white",
    outlier.shape = NA
  ) +
  scale_fill_manual(values = full_pal[estimators_keep]) +
  labs(
    title = "Comparing evidence estimators",
    y = "Log Evidence",
    x = ""
  ) +
  theme_minimal() +
  theme(legend.position = "none")


ggplot(df_efficiency_nolaplace,
       aes(x = mean_time, y = mc_sd, colour = Estimator)) +
  geom_point(size = 4) +
  scale_colour_manual(values = full_pal[estimators_keep]) +
  labs(
    title = "Efficiency comparison",
    x = "Mean runtime (seconds)",
    y = "Monte Carlo SD",
    colour = "Estimator"
  ) +
  theme_minimal()+
  theme(
    legend.title = element_text(size = 14),
    legend.text  = element_text(size = 12)
  )

ggplot(df_compare_nolaplace,
       aes(x = Estimate, y = Estimator, fill = Estimator)) +
  geom_boxplot(alpha = 0.7) +
  scale_fill_manual(values = full_pal[estimators_keep]) +
  labs(
    title = "Comparing evidence estimators",
    x = "Log Evidence",
    y = ""
  ) +
  theme_minimal() +
  theme(legend.position = "none")




#Plots without Laplace

estimators_keep2 <- c("Prior", "HME", "AIS", "Power Posterior", "SMC", "Chib")

df_compare_nolaplace2 <- df_compare %>%
  filter(Estimator %in% estimators_keep2) %>%
  mutate(
    Estimator = factor(Estimator, levels = estimators_keep2)
  )

ggplot(df_compare_nolaplace2,
       aes(x = Estimator, y = Estimate, fill = Estimator)) +
  geom_violin(trim = FALSE, alpha = 0.7) +
  geom_boxplot(
    width = 0.12,
    fill = "white",
    outlier.shape = NA
  ) +
  labs(
    title = "Comparing evidence estimators",
    y = "Log Evidence",
    x = ""
  ) +
  theme_minimal() +
  theme(legend.position = "none")

df_efficiency_nolaplace2 <- df_compare_nolaplace2 %>%
  group_by(Estimator) %>%
  summarise(
    mean_time = mean(Time),
    mc_sd     = sd(Estimate),
    mc_mean   = mean(Estimate),
    .groups   = "drop"
  )

ggplot(df_efficiency_nolaplace2,
       aes(x = mean_time, y = mc_sd, colour = Estimator)) +
  geom_point(size = 4) +
  labs(
    title = "Efficiency comparison",
    x = "Mean runtime (seconds)",
    y = "Monte Carlo SD",
    colour = "Estimator"
  ) +
  theme_minimal()





# get default ggplot discrete palette
levs <- levels(factor(df_efficiency_nolaplace2$Estimator))
pal  <- scales::hue_pal()(length(levs))
names(pal) <- levs

# override only chib
pal["Chib"] <- "#C77CFF"
  # ggplot default purple

ggplot(df_efficiency_nolaplace2,
       aes(x = mean_time, y = mc_sd, colour = Estimator)) +
  geom_point(size = 4) +
  scale_colour_manual(values = pal) +
  labs(
    title = "Efficiency comparison",
    x = "Mean runtime (seconds)",
    y = "Monte Carlo SD",
    colour = "Estimator"
  ) +
  theme_minimal()







#Boxplots
#ggplot(df_compare_nolaplace2,
#       aes(x = Estimator, y = Estimate, fill = Estimator)) +
#  geom_boxplot(
#    alpha = 0.7,
#    width = 0.6,
#    outlier.shape = 16
#  ) +
#  labs(
#    title = "Comparing Evidence Estimators",
#    y = "Log Evidence",
#    x = ""
#  ) +
#  theme_minimal() +
#  theme(legend.position = "none")


ggplot(df_compare_nolaplace2,
       aes(x = Estimate, y = Estimator, fill = Estimator)) +
  geom_boxplot(alpha = 0.7) +
  labs(
    title = "Comparing evidence estimators",
    x = "Log Evidence",
    y = ""
  ) +
  theme_minimal() +
  coord_cartesian(xlim = c(NA, -200))+
  theme(legend.position = "none")



df_laplace <- df_compare %>% 
  filter(Estimator == "Laplace")

df_chib <- df_compare %>% 
  filter(Estimator == "Chib")

ggplot(df_laplace,
       aes(x = Estimator, y = Estimate)) +
  geom_boxplot(
    width = 0.3,
    fill = "#FF61C3",   # choose any colour you like
    outlier.shape = NA
  ) +
  labs(
    title = "Laplace evidence estimates",
    y = "Log Evidence",
    x = ""
  ) +
  coord_flip() +
  theme_minimal()


ggplot(df_chib,
       aes(x = Estimator, y = Estimate)) +
  geom_boxplot(
    width = 0.3,
    fill = "#C77CFF",
    outlier.shape = NA
  ) +
  labs(
    title = "Chib evidence estimates",
    y = "Log Evidence",
    x = ""
  ) +
  #coord_cartesian(xlim = c(NA, -200))+
  coord_flip() +
  theme_minimal()

#df_laplace_chib <- df_compare %>%
#  filter(Estimator %in% c("Laplace", "Chib")) %>%
#  mutate(
#    Estimator = factor(Estimator, levels = c("Laplace", "Chib"))
#  )

#ggplot(df_laplace_chib,
#       aes(x = Estimator, y = Estimate, fill = Estimator)) +
#  geom_boxplot(
#    width = 0.4,
#    outlier.shape = NA
#  ) +
#  labs(
#    title = "Laplace vs Chib Evidence Estimates",
#    y = "Log Evidence",
#    x = ""
#  ) +
#  theme_minimal() +
#  theme(legend.position = "none")


#Decide what plots to save and put in project and put stats in df_efficiency 
#in table in project


