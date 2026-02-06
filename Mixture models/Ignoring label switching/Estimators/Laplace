start=proc.time()

#Bayesian Linear Regression with unknown beta and precision

#Set seed for reproducibility
set.seed(123)

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

#Laplace version 

#l(theta)= log(prior x likelihood)
library(extraDistr)
library(mvtnorm)


library(extraDistr)

l_theta <- function(theta, y, K,
                    alpha, mu0, lambda0, a0, b0) {
  
  # unpack parameters
  eta <- c(theta[1:(K-1)], 0)
  omega <- exp(eta) / sum(exp(eta))
  
  mu <- theta[K:(2*K-1)]
  
  log_sigma2 <- theta[(2*K):(3*K-1)]
  sigma2 <- exp(log_sigma2)
  
  # log-likelihood (marginal mixture)
  ll <- sum(sapply(y, function(yi) {
    m <- max(log(omega) + dnorm(yi, mu, sqrt(sigma2), log = TRUE))
    m + log(sum(exp(
      log(omega) + dnorm(yi, mu, sqrt(sigma2), log = TRUE) - m
    )))
  }))
  
  # log-prior
  lp <- ddirichlet(omega, alpha, log = TRUE) +
    sum(dnorm(mu, mu0, sqrt(1 / lambda0), log = TRUE)) +
    sum(dinvgamma(sigma2, a0, b0, log = TRUE) + log_sigma2)
  
  ll + lp
}


#Hessian using optim function- more general as only need l_theta function

theta_init <- c(
  rep(0, K - 1),
  c(0,11,5),
  rep(log(var(y)), K)
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
  hessian = TRUE,
  control = list(fnscale = -1)
)

#Note we can extract optimal theta(mode) from this
theta_mode<-optimiser$par
#and hessian
Hessian<-optimiser$hessian

#Laplace estimate- on log scale for numerical stability

# Dimension of parameter vector
d <- length(theta_mode)

sigma<- solve(-Hessian)

LogLaplace <- (d / 2) * log(2*pi) +
  0.5 * log(det(sigma)) +
  l_theta(theta_mode, y, K, alpha, mu0, lambda0, a0, b0)
LogLaplace


end=proc.time()

timer<- end-start
timer[3]

#Varies with initial conditions dramatically as several modes
#Could illustrate with visualising data and changing starting points
#and see which mode you get to
