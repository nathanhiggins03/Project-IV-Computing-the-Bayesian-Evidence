#AIS

#Need to change to incorporate more efficient code for HME
library(ggplot2)

Sim <- 30
K_values <- c(3,4,5,6,8,11)

N_lookup <- c(
  `3`  = 2000,
  `4`  = 2000,
  `5`  = 4000,
  `6`  = 4000,
  `8`  = 5000,
  `11` = 5000
)

T_lookup <- c(
  `3`  = 10,
  `4`  = 10,
  `5`  = 15,
  `6`  = 15,
  `8`  = 20,
  `11` = 20
)
df_all <- data.frame()

library(rstan)


rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

setwd("~/Desktop/Project IV")

# Compile Stan models

prior_model <- stan_model("Prior Mixture model.stan")
pp_model <- stan_model("Power Posterior Mixture model.stan")

for (K in K_values) {
  Nsim <- N_lookup[as.character(K)]
  T    <- T_lookup[as.character(K)]
  time_simulation <- rep(0, Sim)
  est_simulation  <- rep(0, Sim)
  
  for (k in 1:Sim) {
    
    start=proc.time()
    set.seed(123+k)
    
    # Data
    
    data(galaxies, package = "MASS")
    y <- galaxies / 1000
    N <- length(y)
    
    # Priors
    
    alpha <- rep(1, K)
    mu0 <- mean(y)
    lambda0 <- 2.6 / (max(y) - min(y))
    a0 <- 1.28
    b0 <- 0.36 * (mean(y^2) - mean(y)^2)
    
    # AIS parameters
    
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
      
      omega_init <- as.numeric(omega_store[t_idx - 1, Nsim, ])
      omega_init <- omega_init / sum(omega_init)
      init <- list(list(
        omega  = omega_init,
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
      Time     = time_simulation,
      MC       = Nsim,
      T        = T,
      K = factor(paste0("K = ", K),
                 levels = paste0("K = ", K_values))
    )
    )
  
}



df_ais <- df_all     # for AIS


#Power posterior

#Need to change to incorporate more efficient code for HME
library(ggplot2)

Sim <- 30
K_values <- c(3,4,5,6,8,11)

N_lookup <- c(
  `3`  = 2000,
  `4`  = 2000,
  `5`  = 4000,
  `6`  = 4000,
  `8`  = 5000,
  `11` = 5000
)

T_lookup <- c(
  `3`  = 10,
  `4`  = 10,
  `5`  = 15,
  `6`  = 15,
  `8`  = 20,
  `11` = 20
)

df_all <- data.frame()

setwd("~/Desktop/Project IV")

library(rstan)


rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())


# Compile Stan models

prior_model <- stan_model("Prior Mixture model.stan")
power_model <- stan_model("Power Posterior Mixture model.stan")

for (K in K_values) {
  Nsim <- N_lookup[as.character(K)]
  T    <- T_lookup[as.character(K)]
  time_simulation <- rep(0, Sim)
  est_simulation  <- rep(0, Sim)
  
  for (k in 1:Sim) {
    
    start <- proc.time()
    set.seed(123+k)
    
    # Data
    
    data(galaxies, package = "MASS")
    y <- galaxies / 1000
    N <- length(y)
    
    # Priors
    
    alpha   <- rep(1, K)
    mu0     <- mean(y)
    lambda0 <- 2.6 / (max(y) - min(y))
    a0      <- 1.28
    b0      <- 0.36 * (mean(y^2) - mean(y)^2)
    
    # Power posterior settings
    
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
      
      #init_list <- list(list(
      #  omega  = omega_samples[t - 1, Nsim, ],
      #  mu     = mu_samples[t - 1, Nsim, ],
      #  sigma2 = sigma2_samples[t - 1, Nsim, ]
      #))
      
      omega_init <- as.numeric(omega_samples[t - 1, Nsim, ])
      omega_init <- omega_init / sum(omega_init)
      
      init_list <- list(list(
        omega  = omega_init,
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
      Time     = time_simulation,
      MC       = Nsim,
      T        = T,
      K = factor(paste0("K = ", K),
                 levels = paste0("K = ", K_values))
    )
    )
  
  
}


df_pp <- df_all      # for TI




library(dplyr)

df_ais   <- df_ais   %>% mutate(Estimator = "AIS")
df_pp    <- df_pp    %>% mutate(Estimator = "Power Posterior")

df_compare <- bind_rows(
  df_ais,
  df_pp
)

df_compare$Estimator <- factor(
  df_compare$Estimator,
  levels = c("AIS", "Power Posterior")
)

ggplot(
  subset(df_compare, Estimator == "AIS"),
  aes(x = K, y = Estimate)
) +
  geom_boxplot() +
  labs(
    title = "AIS estimates for different K",
    x = "K",
    y = "Log evidence"
  ) +
  theme_minimal()

ggplot(
  subset(df_compare, Estimator == "Power Posterior"),
  aes(x = K, y = Estimate)
) +
  geom_boxplot() +
  labs(
    title = "Power posterior estimates for different K",
    x = "K",
    y = "Log evidence"
  ) +
  theme_minimal()

df_efficiency <- df_compare %>%
  group_by(K, Estimator) %>%
  summarise(
    mean_time = mean(Time),
    mc_sd     = sd(Estimate),
    .groups = "drop"
  )

ggplot(df_efficiency,
       aes(x = mean_time,
           y = mc_sd,
           colour = Estimator,
           shape = K)) +
  geom_point(size = 4) +
  labs(
    title = "Efficiency comparison",
    x = "Mean runtime (seconds)",
    y = "Monte Carlo SD"
  ) +
  theme_minimal()




