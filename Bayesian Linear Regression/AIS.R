setwd("~/Desktop/Project IV") 
library(rstan)
library(durhamSLR) 
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
# Compile the Stan model once
stan_model_obj <- stan_model('Power Posterior Bayes Linear Regression.stan')
prior_model <- stan_model(
  'Prior Bayesian Linear Regression.stan'
)
start=proc.time()
#Adapt for Bayes Linear Regression
#Set seed for reproducibility
set.seed(123)

library(mvtnorm)
library(extraDistr)


# --- 1. DATA, PRIORS, and HYPERPARAMETERS ---

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


# AIS Parameters
T <- 40             # Number of tempered distributions 
Nsim <- 4000           # Number of AIS particles/samples
c_power <- 2        # Power schedule exponent


# Function for the log-power posterior (used for calculating weights)
# log p_t(theta) = t * log L(x|theta) + log pi(theta)
power_post_log<- function(t, theta,X,y){
  Beta<- theta[1:2]
  sigma_sq<-theta[3]
  mean<-X %*% Beta
  sd <- sqrt(sigma_sq)
  out <- t*sum(dnorm(y, mean = mean, sd = sd, log = TRUE)) +
    dmvnorm(Beta, mean = m0, sigma = sigma_sq * solve(Lambda0), log = TRUE) +
    dinvgamma(sigma_sq, alpha0, beta0, log = TRUE)
  return(out)
}

# Tempering function
Beta_t <- function(t, c=c_power){
  return((t/T)^c)
}
t_list <- Beta_t(0:T) # T+1 tempered posteriors


mh_step <- function(theta, beta_temp) {
  Beta <- theta[1:2]
  sigma_sq <- theta[3]
  
  # --- Proposals ---
  Beta_prop <- Beta + rnorm(2, 0, 0.2)
  
  # log-scale proposal for variance (IMPORTANT)
  sigma_sq_prop <- sigma_sq * exp(rnorm(1, 0, 0.1))
  
  theta_prop <- c(Beta_prop, sigma_sq_prop)
  
  # --- Log densities ---
  log_curr <- power_post_log(beta_temp, theta, X, y)
  log_prop <- power_post_log(beta_temp, theta_prop, X, y)
  
  #M-H accept/reject
  log_accept_ratio <- (log_prop - log_curr) + 
    log(sigma_sq_prop) - log(sigma_sq)
  
  if (log(runif(1)) < log_accept_ratio) {
    return(theta_prop)
  } else {
    return(theta)
  }
}

# Sample from joint prior to initialise the arrays
prior_fit <- sampling(
  prior_model,
  data = list(N=length(y), d=d, X=X,
              m0=m0, alpha0=alpha0,
              beta0=beta0, Lambda0 = Lambda0),
  iter = 200000,
  chains = 1,
  algorithm = "Fixed_param",
  refresh = 0
)

#Extract Prior samples
prior_sample_sigma_sq_all<- extract(prior_fit, pars = 'sigma_sq')$'sigma_sq'  
prior_sample_Beta1_all<- extract(prior_fit, pars = 'Beta[1]')$'Beta[1]'
prior_sample_Beta2_all<- extract(prior_fit, pars = 'Beta[2]')$'Beta[2]'
prior_sample_Beta_all<- cbind(prior_sample_Beta1_all,prior_sample_Beta2_all)

prior_sample_sigma_sq<- tail(prior_sample_sigma_sq_all, Nsim)
prior_sample_Beta1<- tail(prior_sample_Beta1_all, Nsim)
prior_sample_Beta2<- tail(prior_sample_Beta2_all, Nsim)
prior_sample_Beta<- tail(prior_sample_Beta_all, Nsim)

#Need to store outputs- samples of mu,tau in seperate T+1 x N matrices

beta_samples<-array(0, dim = c(2,T + 1, Nsim))
sigma_sq_samples<- matrix(0,nrow = T+1, ncol = Nsim)

#Store initial samples in matrix
beta_samples[1,1,]<-prior_sample_Beta1
beta_samples[2,1,]<-prior_sample_Beta2
sigma_sq_samples[1,]<-prior_sample_sigma_sq

#Need to store log of unnormalised weights
log_w<-matrix(0,nrow = T+1, ncol = Nsim)

# Set initial log weights.
log_w[1,] <- 0

# AIS MAIN LOOP

# Define the extended mixing parameters
#N_samples_kept <- Nsim          # 2000 samples kept (for weight calculation)
#N_samples_warmup <- 4000     # 4000 iterations for warmup


for(t_index in 2:(T+1)){
  power_curr <- t_list[t_index]
  power_prev <- t_list[t_index - 1]
  
  cat(sprintf("Running step %d of %d: t_prev=%.4f -> t_curr=%.4f\n", 
              t_index - 1, T, power_prev, power_curr))
  
  for (i in 1:Nsim) {
    
    # --- Previous particle ---
    theta_prev <- c(beta_samples[,t_index - 1, i], 
                    sigma_sq_samples[t_index - 1, i])
    
    # ---- 1. Weight update (UNCHANGED) ----
    log_w_curr <- power_post_log(power_curr, theta_prev, X, y)
    log_w_prev <- power_post_log(power_prev, theta_prev, X, y)
    
    log_w[t_index, i] <- log_w[t_index - 1, i] + (log_w_curr - log_w_prev)
    
    # ---- 2. MH transition ----
    theta_new <- mh_step(theta_prev, power_curr)
    
    # Store updated particle
    beta_samples[,t_index, i] <- theta_new[1:2]
    sigma_sq_samples[t_index, i] <- theta_new[3]
  }
}
# The final estimate is log(E_i[w_i]) using Log-Sum-Exp for stability.


final_log_weights <- log_w[T+1, ]

# Find max log-weight (m)
m <- max(final_log_weights)

# Calculate log(sum(exp(log_w - m))) - log(N)
AIS_log_evidence <- m + log(mean(exp(final_log_weights - m)))

# Output Results
cat("\n--- AIS Results ---\n")
cat("Number of particles (N):", Nsim, "\n")
cat("Number of steps (T):", T, "\n")
cat("Estimated Log Evidence (log Z):", AIS_log_evidence, "\n")
# --- AIS SCRIPT END ---

end=proc.time()

timer<- end-start
timer[3]

