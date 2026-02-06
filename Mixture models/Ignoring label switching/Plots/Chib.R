#NEED TO FIX

#Need to change to incorporate more efficient code for HME
library(ggplot2)
library(extraDistr)

Sim <- 30
MC_values <- c(100)   # <-- choose MC sizes

df_all <- data.frame()


for (MC_sample in MC_values) {
  
  time_simulation <- rep(0, Sim)
  est_simulation  <- rep(0, Sim)
  
  for (k in 1:Sim) {
    
    start <- proc.time()
    
    set.seed(123+10*k)
    
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
      S = MC_sample,
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
    
    est_simulation[k]<-chibs
    
    time_simulation[k] <- (proc.time() - start)[3]
    
    
    #theta_init changes the estimation as finds different posterior modes 
    #due to label switching
    
    #Theta_init1 gives -2 and theta_init2 gives -190
    #Effect of label swithcing/multiple modes
    
    
    #Chib’s estimator assumes:
    #  “There is a single dominant posterior mode, and I can evaluate the posterior density at that point.”
    #Finite mixtures violate this assumption by symmetry.
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
    title = "Monte Carlo convergence of Chib's estimator",
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
  geom_text(aes(label = MC), vjust = -0.4, hjust=-0.2) +
  labs(
    title = "Monte Carlo efficiency - Chib's estimator",
    x = "Mean runtime (seconds)",
    y = "Monte Carlo SD"
  ) +
#  coord_cartesian(xlim = c(NA, 20))+
  theme_minimal()


