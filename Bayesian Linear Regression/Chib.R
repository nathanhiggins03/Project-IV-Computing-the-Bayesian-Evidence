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


#Chibs Method

#First need to set up a gibbs sampler for the model. Need FCDs found by getting joint posterior up to 
#proportionality and getting rid of any terms that don't have beta/sigma_sq for FCD of beta/sigma_sq respectively
library(mvtnorm)
library(extraDistr)

gibbs <- function(N,d, X,y, m0, Lambda0, alpha0, beta0) {
  #Initialise with prior means
  mat <- matrix(0, ncol=3, nrow=N)
  Beta <- m0
  sigma_sq <- beta0 / (alpha0-1)
  mat[1,]<-c(Beta,sigma_sq)
  
  #FCD paramters
  n <- length(y)
  for (i in 2:N) {
    LambdaN<- Lambda0 + t(X)%*%X
    mN <- solve(LambdaN, Lambda0 %*% m0 + t(X) %*% y)
    Beta <- as.numeric(rmvnorm(1, mean = mN, sigma = sigma_sq * solve(LambdaN)))
    
    alphanN<- alpha0 + (n+d)/2
    epsilon <- y - X %*% Beta
    betaN<- beta0 + 0.5 * (t(epsilon) %*% epsilon + t(Beta - m0) %*% Lambda0 %*% (Beta - m0))
    sigma_sq <- rinvgamma(1, alpha = alphanN, beta = betaN) 
    mat[i,] <- c(Beta, sigma_sq)
  }
  return(mat)
}

set.seed(3421)
out1=gibbs(N=5000,d=d, X=X,y=y, m0=m0, Lambda0=Lambda0, alpha0=alpha0, beta0=beta0)
#Column mean
colMeans(out1)
#Compare to STAN posterior samples
#print(posterior_sample)

#l(theta)
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

#Need to continue from here









#Now estimate log posterior at theta*
beta_gibbs<-out1[,1:2]
sigma_sq_gibbs<- out1[,3]

#We have access to this FCD(for Beta)
LambdaN<- Lambda0 + t(X)%*%X
mN <- solve(LambdaN, Lambda0 %*% m0 + t(X) %*% y)
log_posterior_p1<-dmvnorm(theta_mode[1:2], mean = mN, sigma = theta_mode[3] * solve(LambdaN), log = TRUE)

#Need to estimate this distribution using Monte Carlo using samples from FCD(gibbs samples)
counter<-rep(0,length(sigma_sq_gibbs))
for(i in 1:length(sigma_sq_gibbs)){
  n <- length(y)
  alphanN<- alpha0 + (n+d)/2
  epsilon <- y - X %*% beta_gibbs[i,]
  betaN<- beta0 + 0.5 * (t(epsilon) %*% epsilon + t(beta_gibbs[i,] - m0) %*% Lambda0 %*% (beta_gibbs[i,] - m0))
  counter[i] <- dinvgamma(theta_mode[3], alpha = alphanN, beta = betaN, log = TRUE) 
}
# log-sum-exp trick
max_log <- max(counter)
log_posterior_p2 <- max_log + log(mean(exp(counter - max_log)))

#Complete log posterior
log_posterior<- log_posterior_p1 + log_posterior_p2


#Chibs estimator

lprior<-dmvnorm(theta_mode[1:2], mean = m0, sigma = theta_mode[3] * solve(Lambda0), log = TRUE) + dinvgamma(theta_mode[3], alpha0, beta0, log = TRUE)
llike<- sum(dnorm(y, mean = X %*% theta_mode[1:2], sd = sqrt(theta_mode[3]), log = TRUE))


chibs<- lprior + llike - log_posterior
chibs


end=proc.time()

timer<- end-start
timer[3]
