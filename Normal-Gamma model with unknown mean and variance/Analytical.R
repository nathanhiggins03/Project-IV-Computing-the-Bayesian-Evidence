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

#Analytical log evidence
true_le<- log(gamma(alpha1)) -log(gamma(alpha0)) + alpha0* log(beta0) - alpha1 * log(beta1) + 0.5*(log(k0)-log(k1)) - n/2 * log(2*pi)
true_le



