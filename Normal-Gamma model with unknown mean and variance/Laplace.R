#Normal-Gamma model with unknown mean and variance

#Set seed for reproducibility
set.seed(123)

#Data
x<-rnorm(10,mean = 10,sd=2)

#Inputs
m0<-10
k0<-0.1
alpha0<-1
beta0<-0.1
n<- length(x)

#Tractable posterior inputs
k1<- k0 +n
m1<- (k0*m0 + n*mean(x))/k1
alpha1<- alpha0 + n/2
beta1<- beta0 + 1/2 * (n-1) * var(x) + (k0 *n * (mean(x) - m0)^2)/(2 * k1)

#Laplace Estimator

#l(theta)= log(prior x likelihood)

l_theta<- function(theta, X){
  out<- sum(dnorm(X, theta[1], sqrt(1/theta[2]), log = TRUE)) +dnorm(theta[1], m0, sqrt(1/(k0*theta[2])), log = TRUE) + dgamma(theta[2], alpha0,beta0, log = TRUE)
  return(out)
}

#Hessian using optim function- more general as only need l_theta function

optimiser<- optim(par= c(10,1), fn=l_theta,X=x,hessian = TRUE,control = list(fnscale = -1))
#Note we can extract optimal theta(mode) from this
theta_mode<-optimiser$par
#and hessian
Hessian<-optimiser$hessian

#Laplace estimate- on log scale for numerical stability

# Dimension of parameter vector
d <- length(theta_mode)

sigma<- solve(-Hessian)

LogLaplace<- (d / 2)*log(2*pi) +0.5*log(det(sigma)) + l_theta(theta_mode, x)
LogLaplace
