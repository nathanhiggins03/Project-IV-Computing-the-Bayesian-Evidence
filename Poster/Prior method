library(ggplot2)

Sim <- 30
MC_values <- c(1000, 10000, 100000,1000000)   # <-- choose MC sizes

df_all <- data.frame()

library(rstan)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

setwd("~/Desktop/Project IV")

prior_model <- stan_model(
  file = "Prior Bayesian Linear Regression.stan"
)

for (MC_sample in MC_values) {
  
  time_simulation <- rep(0, Sim)
  est_simulation  <- rep(0, Sim)
  
  for (k in 1:Sim) {
    
    start=proc.time()
    #Set seed for reproducibility
    set.seed(123+k)
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
    beta0  <- 36
    
    
    #Posterior distribution inputs
    LambdaN <- Lambda0 + t(X) %*% X
    mN <- solve(LambdaN, Lambda0 %*% m0 + t(X) %*% y)
    alphaN <- alpha0 + N / 2
    betaN <- beta0 + 0.5 * (t(y) %*% y + t(m0) %*% Lambda0 %*% m0 - t(mN) %*% LambdaN %*% mN)
    
    
    
    #Mean of likelihood using prior samples:
    
    #Sample from prior using STAN
    
    
    
    #Need lots of draws for it to work
    prior_fit <- sampling(
      prior_model,
      data = list(
        N = N, d = d, X = X,
        m0 = m0, alpha0 = alpha0,
        beta0 = beta0, Lambda0 = Lambda0
      ),
      iter = 2*MC_sample,
      chains = 1,
      algorithm = "Fixed_param",
      refresh = 0
    )
    
    #Diagnostic checks
    #print(prior_fit)
    #output = as.array(prior_fit)
    #diagnostics(output)
    
    
    #Extract Prior samples
    prior_sample_sigma_sq<- extract(prior_fit, pars = 'sigma_sq')$'sigma_sq'  
    prior_sample_Beta1<- extract(prior_fit, pars = 'Beta[1]')$'Beta[1]'
    prior_sample_Beta2<- extract(prior_fit, pars = 'Beta[2]')$'Beta[2]'
    prior_sample_Beta<- cbind(prior_sample_Beta1,prior_sample_Beta2)
    
    #Calculate log likelihood at each prior sample (mu,tau)
    
    likelihood_prior_log<- rep(0, length(prior_sample_Beta1))
    for(i in 1:length(prior_sample_Beta1)){
      mean<-as.vector(X %*% prior_sample_Beta[i,])
      likelihood_prior_log[i]<- sum(dnorm(y,mean = mean, sd = sqrt(prior_sample_sigma_sq[i]) , log = TRUE))
    }
    
    #Apply log sum exp trick to log(mean(likelihood))
    m <- max(likelihood_prior_log)
    est_simulation[k] <- m + log(mean(exp(likelihood_prior_log - m)))
    
    time_simulation[k] <- (proc.time() - start)[3]
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

#True value
#Data
#y=2+3x+ϵ,ϵ∼N(0,0.5^2)
set.seed(123)
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
beta0  <- 36


#Posterior distribution inputs
LambdaN <- Lambda0 + t(X) %*% X
mN <- solve(LambdaN, Lambda0 %*% m0 + t(X) %*% y)
alphaN <- alpha0 + N / 2
betaN <- beta0 + 0.5 * (t(y) %*% y + t(m0) %*% Lambda0 %*% m0 - t(mN) %*% LambdaN %*% mN)


#Analytical log evidence

term1 <- as.numeric(0.5 * (determinant(Lambda0, log = TRUE)$modulus - determinant(LambdaN, log = TRUE)$modulus))
term2 <- alpha0 * log(beta0) - alphaN * log(betaN)
term3 <- log(gamma(alphaN)) - log(gamma(alpha0))
term4 <- -(N/2) * log(2*pi)
true_le<- term1 + term2 + term3 + term4
true_le


#Plotting evidence violin plots against MC sample size
ggplot(df_all, aes(x = MC, y = Estimate)) +
  geom_violin(aes(fill = MC),
              trim = FALSE,
              alpha = 0.7) +
  geom_boxplot(width = 0.12,
               fill = "white",
               outlier.shape = NA) +
  geom_hline(aes(yintercept = true_le, color = "Analytical\nvalue"),
             linewidth = 1.2,
             na.rm = TRUE) +
  scale_color_manual(
    name = "",
    values = c("Analytical\nvalue" = "red")
  ) +
  labs(
    title = "Monte Carlo convergence of the Prior method",
    x = "N",
    y = "Log Evidence"
  ) +
  theme_minimal()+
  theme(
    legend.position = "right",
    axis.text  = element_text(size = 18, colour = "black"),
    axis.title = element_text(size = 20, colour = "black"),
    plot.title = element_text(size = 20, colour = "black"),
    legend.text  = element_text(size = 20, colour = "black"),
    legend.title = element_text(size = 20, colour = "black")
  )

