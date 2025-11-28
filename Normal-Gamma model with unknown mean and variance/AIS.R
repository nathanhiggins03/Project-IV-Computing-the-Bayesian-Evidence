start=proc.time()
#Set seed for reproducibility
set.seed(123)

setwd("~/Desktop/Project IV") 
library(rstan)
library(durhamSLR) 
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# --- 1. DATA, PRIORS, and HYPERPARAMETERS ---

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

# AIS Parameters
T <- 50             # Number of tempered distributions 
N <- 4000           # Number of AIS particles/samples
c_power <- 2        # Power schedule exponent


# Function for the log-power posterior (used for calculating weights)
# log p_t(theta) = t * log L(x|theta) + log pi(theta)
power_post_log <- function(t, theta, X, m0, k0, alpha0, beta0){
  mu <- theta[1]
  tau <- theta[2]
  
  # 1. Tempered Log-Likelihood: t * log L(x | mu, tau)
  ll_t <- t * sum(dnorm(X, mean = mu, sd = sqrt(1/tau), log = TRUE))
  
  # 2. Log-Priors: log pi(mu | tau) + log pi(tau)
  lp_mu <- dnorm(mu, mean = m0, sd = sqrt(1/(k0 * tau)), log = TRUE)
  lp_tau <- dgamma(tau, shape = alpha0, rate = beta0, log = TRUE)
  
  return(ll_t + lp_mu + lp_tau)
}

# Tempering function
Beta_t <- function(t, c=c_power){
  return((t/T)^c)
}
t_list <- Beta_t(0:T) # T+1 tempered posteriors


# Compile the Stan model once
stan_model_obj <- stan_model('Power Posterior Normal Normal model with unknown mean and precision.stan')

# Sample from joint prior to initialise the arrays
prior_fit = stan('Prior Normal Normal with unknown mean and precision.stan', 
                 data = list(k0=k0, m0=m0, alpha0=alpha0, beta0=beta0), 
                 iter = 2*N, chains = 1, algorithm = "Fixed_param")

# Extract Prior samples (taking the last N samples)
prior_sample_mu <- tail(extract(prior_fit, pars = 'mu')$'mu', N) 
prior_sample_tau <- tail(extract(prior_fit, pars = 'tau')$'tau', N)


# Need to store outputs: [T+1 x N matrices]
mu_samples <- matrix(0, nrow = T + 1, ncol = N)
tau_samples <- matrix(0, nrow = T + 1, ncol = N)
log_w <- matrix(0, nrow = T + 1, ncol = N)

# Store initial samples (t=0)
mu_samples[1,] <- prior_sample_mu
tau_samples[1,] <- prior_sample_tau

# Set initial log weights.
log_w[1,] <- 0

# AIS MAIN LOOP

# Define the extended mixing parameters
N_samples_kept <- N          # 2000 samples kept (for weight calculation)
N_samples_warmup <- 4000     # 4000 iterations for warmup

for(t_index in 2:(T+1)){
  power_curr <- t_list[t_index]
  power_prev <- t_list[t_index - 1]
  
  cat(sprintf("Running step %d of %d: t_prev=%.4f -> t_curr=%.4f\n", 
              t_index - 1, T, power_prev, power_curr))
  
  
  # 1. Setup warm-start initialisation using the last particle from the previous step.
  init_list <- list(list(mu = mu_samples[t_index - 1, N] , tau = tau_samples[t_index - 1, N]))
  
  # 2. Run MCMC
  run <- sampling(stan_model_obj, 
                  data = list(x=x, N=length(x), k0=k0, m0=m0, alpha0=alpha0, beta0=beta0, t=power_curr), 
                  iter = N_samples_warmup + N_samples_kept, # Total 6000 iterations
                  warmup = N_samples_warmup,                # 4000 warmup iterations
                  chains = 1,        
                  init = init_list,  
                  refresh = 0)
  
  # 3. Extract new samples theta_t
  # The 'extract' function automatically returns the N_samples_kept iterations
  mu_samples[t_index,] <- extract(run, pars = 'mu')$'mu' 	
  tau_samples[t_index,] <- extract(run, pars = 'tau')$'tau' 	
  
  for (i in 1:N) {
    theta_prev <- c(mu_samples[t_index - 1, i], tau_samples[t_index - 1, i])
    
    # log(p_t/p_{t-1}) = log p_t - log p_{t-1}
    log_w_curr <- power_post_log(power_curr, theta_prev, x, m0, k0, alpha0, beta0)
    log_w_prev <- power_post_log(power_prev, theta_prev, x, m0, k0, alpha0, beta0)
    
    log_w[t_index, i] <- log_w[t_index - 1, i] + (log_w_curr - log_w_prev)
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
cat("Number of particles (N):", N, "\n")
cat("Number of steps (T):", T, "\n")
cat("Estimated Log Evidence (log Z):", AIS_log_evidence, "\n")
# --- AIS SCRIPT END ---

end=proc.time()

timer<- end-start
timer[3]

