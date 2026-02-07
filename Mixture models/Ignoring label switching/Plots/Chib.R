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

chib_estimates <- replicate(30, run_chib())

# Plot results
df_plot <- data.frame(Estimate = chib_estimates)

#Regular
ggplot(df_plot, aes(x = "Chib Estimates", y = Estimate)) +
  geom_violin(fill = "skyblue", alpha = 0.6) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  labs(title = "Distribution of Chib's estimator over 30 runs",
       x = "", y = "Log Evidence") +
  theme_minimal()

#Log of log evidence
ggplot(df_plot, aes(x = "Chib Estimates", y = Estimate)) +
  geom_violin(fill = "skyblue", alpha = 0.6) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  labs(title = "Distribution of Chib's estimator over 30 runs",
       x = "", y = "Log Evidence") +
  theme_minimal() +
  scale_y_continuous(trans = "log10")  # log scale handles large differences

library(ggbeeswarm)  # install.packages("ggbeeswarm")
#Raincloud
ggplot(df_plot, aes(x = "Chib Estimates", y = Estimate)) +
  geom_violin(fill = "skyblue", alpha = 0.4) +
  geom_quasirandom(size = 2, color = "darkblue", alpha = 0.7) +
  labs(title = "Chib estimates - raincloud view", x = "", y = "Log Evidence") +
  theme_minimal()


library(ggplot2)
#Boxplot
# df_plot contains 30 Chib estimates
df_plot <- data.frame(
  Estimate = chib_estimates,
  Run = factor(1:30)
)

ggplot(df_plot, aes(x = "", y = Estimate)) +
  geom_boxplot(fill = "skyblue", alpha = 0.7, outlier.color = "red", outlier.shape = 16) +
  geom_jitter(width = 0.2, height = 0, color = "darkblue", alpha = 0.6, size = 2) +
  labs(
    title = "Chib's estimator over 30 runs",
    x = "",
    y = "Log Evidence"
  ) +
  theme_minimal()

library(ggplot2)
#Scatterplot
df_plot <- data.frame(
  Estimate = chib_estimates,
  Run = 1:30
)

ggplot(df_plot, aes(x = Run, y = Estimate)) +
  geom_point(color = "blue", size = 3, alpha = 0.7) +
  geom_line(color = "gray70", linetype = "dashed") +  # optional, connects points
  labs(
    title = "Chib's estimator over 30 runs",
    x = "Run",
    y = "Log Evidence"
  ) +
  theme_minimal()


