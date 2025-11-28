start=proc.time()

#Bayesian Linear Regression with unknown beta and precision

#Set seed for reproducibility
set.seed(123)
#Data
#y=2+3x+ϵ,ϵ∼N(0,0.5^2)

N <- 30
x <- runif(N, 0, 5)
X <- cbind(1, x)   # include intercept
beta_true <- c(2, 3)
sigma_true <- 0.5
y <- as.vector(X %*% beta_true + rnorm(N, 0, sigma_true))

#Data we have is
X
y

#Prior inputs
d <- ncol(X)
m0 <- rep(1, d)
Lambda0 <- diag(5, d)   # vague prior
alpha0 <- 3
beta0  <- 2 * var(y)


#Posterior distribution inputs
LambdaN <- Lambda0 + t(X) %*% X
mN <- solve(LambdaN, Lambda0 %*% m0 + t(X) %*% y)
alphaN <- alpha0 + N / 2
betaN <- beta0 + 0.5 * (t(y) %*% y + t(m0) %*% Lambda0 %*% m0 - t(mN) %*% LambdaN %*% mN)


#Laplace version 

#l(theta)= log(prior x likelihood)
library(extraDistr)
library(mvtnorm)


l_theta<- function(theta, X, y){
  Beta<- theta[1:2]
  sigma_sq<-theta[3]
  mean<-X %*% Beta
  sd <- sqrt(sigma_sq)
  out <- sum(dnorm(y, mean = mean, sd = sd, log = TRUE)) +
    dmvnorm(Beta, mean = m0, sigma = sigma_sq * solve(Lambda0), log = TRUE) +
    dinvgamma(sigma_sq, alpha0, beta0, log = TRUE)
  return(out)
}

#Hessian using optim function- more general as only need l_theta function

optimiser<- optim(par= c(2,3,0.5), fn=l_theta,X=X, y=y,hessian = TRUE,control = list(fnscale = -1))
#Note we can extract optimal theta(mode) from this
theta_mode<-optimiser$par
#and hessian
Hessian<-optimiser$hessian

#Laplace estimate- on log scale for numerical stability

# Dimension of parameter vector
d <- length(theta_mode)

sigma<- solve(-Hessian)

LogLaplace<- (d / 2)*log(2*pi) +0.5*log(det(sigma)) + l_theta(theta_mode, X, y)



end=proc.time()

timer<- end-start
timer[3]
