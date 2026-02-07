#Need to change to incorporate more efficient code for HME
library(ggplot2)

Sim <- 30
MC_values <- c(1000)   # <-- choose MC sizes
T_value<-10
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


#Plotting evidence violin plots against MC sample size
ggplot(df_all, aes(x = MC, y = Estimate)) +
  geom_violin(aes(fill = MC),
              trim = FALSE,
              alpha = 0.7) +
  geom_boxplot(width = 0.12,
               fill = "white",
               outlier.shape = NA) +
  scale_color_manual(
    name = ""
  ) +
  labs(
    title = "Monte Carlo convergence of the AIS estimator",
    x = "N",
    y = "Log Evidence"
  ) +
  theme_minimal()+
  theme(
    legend.position = "right",
    legend.title = element_blank()
  )


#Mean runtime vs MC sample size
library(dplyr)

df_time <- df_all %>%
  group_by(MC) %>%
  summarise(
    mean_time = mean(Time),
    sd_time   = sd(Time),
    .groups = "drop"
  )

#ggplot(df_time, aes(x = MC, y = mean_time, group = 1)) +
#  geom_point(size = 3) +
#  geom_line() +
#  labs(
#    title = "Mean Runtime vs Monte Carlo Sample Size",
#    x = "Monte Carlo sample size",
#    y = "Mean runtime (seconds)"
#  ) +
#  theme_minimal()



#Plotting MC standard deviation against run time
df_error <- df_all %>%
  group_by(MC) %>%
  summarise(
    mc_sd = sd(Estimate)
  )
df_efficiency <- left_join(df_time, df_error, by = "MC")
ggplot(df_efficiency, aes(x = mean_time, y = mc_sd)) +
  geom_point(size = 3) +
  geom_line()+
  geom_text(aes(label = MC), vjust = -0.5, hjust=-0.1) +
  labs(
    title = "Monte Carlo efficiency - AIS",
    x = "Mean runtime (seconds)",
    y = "Monte Carlo SD"
  ) +
#  coord_cartesian(xlim = c(NA, 13.5))+
  theme_minimal()
