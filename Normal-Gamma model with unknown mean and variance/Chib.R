start=proc.time()
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

#Chibs Method

#First need to set up a gibbs sampler for the model. Need FCDs found by getting joint posterior up to 
#proportionality and getting rid of any terms that don't have mu/tau for FCD of mu/tau respectively

gibbs <- function(N, x, m0, k0, alpha0, beta0) {
  n <- length(x)
  xbar <- mean(x)
  mat <- matrix(0, ncol=2, nrow=N)
  mu <- m0
  tau <- alpha0 / beta0
  #Initialise with prior means
  mat[1,]<-c(mu,tau)
  for (i in 2:N) {
    mu <- rnorm(1, (k0*m0 + n*xbar)/(k0 + n), sqrt(1/(k0*tau + n*tau)))
    tau <- rgamma(1, shape = alpha0 + n/2,
                  rate= beta0 + 0.5*sum((x - mu)^2) + 0.5*k0*(mu - m0)^2)
    mat[i,] <- c(mu, tau)
  }
  return(mat)
}

#set.seed(3421)
out1=gibbs(N=5000,x,m0,k0,alpha0,beta0)
#Column mean
#colMeans(out1)

#Direct sampling from posterior(or STAN) to check Gibbs correct
#direct_sample <- function(N, m1, k1, alpha1, beta1) {
#  tau <- rgamma(N, alpha1, rate = beta1)
#  mu <- rnorm(N, m1, sqrt(1/(k1 * tau)))
#  cbind(mu, tau)
#}

#set.seed(3421)
#out2 <- direct_sample(5000, m1, k1, alpha1, beta1)
#Column means
#colMeans(out2)

#Checking gibbs sampling converges to STAN sample
#set.seed(42)
#G <- gibbs(50000, x, m0, k0, alpha0, beta0)
#D <- direct_sample(50000, m1, k1, alpha1, beta1)

#mean(G[,2]); mean(D[,2])   # compare E[tau]
#sd(G[,2]); sd(D[,2])       # compare SD[tau]




#Now find mu* and tau* points of highest posterior density
#theta_mode<- c(m1, (alpha1-1)/beta1)

#l(theta)= log(prior x likelihood)

l_theta<- function(theta, X){
  out<- sum(dnorm(X, theta[1], sqrt(1/theta[2]), log = TRUE)) +dnorm(theta[1], m0, sqrt(1/(k0*theta[2])), log = TRUE) + dgamma(theta[2], alpha0,beta0, log = TRUE)
  return(out)
}

#Hessian using optim function- more general as only need l_theta function
optimiser<- optim(par= c(10,1), fn=l_theta,X=x,hessian = TRUE,control = list(fnscale = -1))
#Note we can extract optimal theta(mode) from this
theta_mode<-optimiser$par

#Now estimate log posterior at theta*
mu_gibbs<-out1[,1]
#We have access to this FCD
log_posterior_p1 <- dnorm(theta_mode[1],(k0*m0 + sum(x))/(k0 + n), 
                          sqrt(1/(k0*theta_mode[2] + n*theta_mode[2])), 
                          log = TRUE)

#Need to estimate this distribution using Monte Carlo using samples from FCD(gibbs samples)
counter<-rep(0,length(mu_gibbs))
for(i in 1:length(mu_gibbs)){
  counter[i]<-dgamma(theta_mode[2], shape = alpha0 + n/2,
                     rate= beta0 + 0.5*sum((x - mu_gibbs[i])^2) + 0.5*k0*(mu_gibbs[i] - m0)^2, 
                     log = TRUE)
}
# log-sum-exp trick
max_log <- max(counter)
log_posterior_p2 <- max_log + log(mean(exp(counter - max_log)))

#Complete log posterior
log_posterior<- log_posterior_p1 + log_posterior_p2


#Chibs estimator

lprior<-dnorm(theta_mode[1], m0,sqrt(1/(k0*theta_mode[2])), log = TRUE) + dgamma(theta_mode[2], alpha0, beta0, log = TRUE)
llike<- sum(dnorm(x,theta_mode[1], sqrt(1/theta_mode[2]), log=TRUE))

chibs<- lprior + llike - log_posterior
chibs

end=proc.time()

timer<- end-start
timer[3]
