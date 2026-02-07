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

############################
# 1. Violin + boxplot
############################
ggplot(df, aes(x = "", y = Estimate)) +
  geom_violin(fill = "skyblue", alpha = 0.6) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  labs(
    title = "Laplace approximation across 30 initialisations",
    x = "",
    y = "Log Evidence"
  ) +
  theme_minimal()

############################
# 2. Raincloud plot
############################
ggplot(df, aes(x = "", y = Estimate)) +
  geom_violin(fill = "skyblue", alpha = 0.4) +
  geom_quasirandom(size = 2, alpha = 0.8) +
  labs(
    title = "Laplace estimates (raincloud plot)",
    x = "",
    y = "Log Evidence"
  ) +
  theme_minimal()

############################
# 3. Efficiency plot
############################
ggplot(df, aes(x = Time, y = Estimate)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_smooth(se = FALSE, linetype = "dashed") +
  labs(
    title = "Laplace efficiency: runtime vs estimate",
    x = "Runtime (seconds)",
    y = "Log Evidence"
  ) +
  theme_minimal()


############################
# Efficiency summary
############################

df_eff <- df %>%
  summarise(
    mean_time = mean(Time),
    sd_est    = sd(Estimate)
  )

############################
# Efficiency plot
############################
ggplot(df_eff, aes(x = mean_time, y = sd_est)) +
  geom_point(size = 4) +
  geom_text(
    label = "Laplace",
    vjust = -0.8
  ) +
  labs(
    title = "Monte Carlo efficiency of Laplace approximation",
    x = "Mean runtime (seconds)",
    y = "SD of log evidence"
  ) +
  theme_minimal()
