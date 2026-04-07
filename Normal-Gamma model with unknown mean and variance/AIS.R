setwd("~/Desktop/Project IV") 
library(rstan)
library(durhamSLR) 
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# Compile the Stan model once
stan_model_obj <- stan_model('Power Posterior Normal Normal model with unknown mean and precision.stan')
prior_model <- stan_model('Prior Normal Normal with unknown mean and precision.stan')

start=proc.time()
#Set seed for reproducibility
set.seed(123)


#Data/priors

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
c_power <- 2        # Power schedule


# Function for the log-power posterior 
# log p_t(theta) = t * log L(x|theta) + log pi(theta)
power_post_log <- function(t, theta, X, m0, k0, alpha0, beta0){
  mu <- theta[1]
  tau <- theta[2]
  
  ll_t <- t * sum(dnorm(X, mean = mu, sd = sqrt(1/tau), log = TRUE))

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

prior_fit <- sampling(
  prior_model,
  data = list(k0=k0, m0=m0, alpha0=alpha0, beta0=beta0),
  iter = 2*N,
  chains = 1,
  algorithm = "Fixed_param",
  refresh = 0
)
# Extract Prior samples (taking the last N samples)
prior_sample_mu <- tail(extract(prior_fit, pars = 'mu')$'mu', N) 
prior_sample_tau <- tail(extract(prior_fit, pars = 'tau')$'tau', N)


# Need to store outputs:
mu_samples <- matrix(0, nrow = T + 1, ncol = N)
tau_samples <- matrix(0, nrow = T + 1, ncol = N)
log_w <- matrix(0, nrow = T + 1, ncol = N)

# Store initial samples (t=0)
mu_samples[1,] <- prior_sample_mu
tau_samples[1,] <- prior_sample_tau

# Set initial log weights.
log_w[1,] <- 0

# AIS MAIN LOOP

for(t_index in 2:(T+1)){
  power_curr <- t_list[t_index]
  power_prev <- t_list[t_index - 1]
  
  cat(sprintf("Running step %d of %d: t_prev=%.4f -> t_curr=%.4f\n", 
              t_index - 1, T, power_prev, power_curr))
  
  for (i in 1:N) {
    
    # Previous particle (this is now a TRUE AIS trajectory)
    theta_prev <- c(mu_samples[t_index - 1, i], 
                    tau_samples[t_index - 1, i])
    
#Weight update
    log_w_curr <- power_post_log(power_curr, theta_prev, x, m0, k0, alpha0, beta0)
    log_w_prev <- power_post_log(power_prev, theta_prev, x, m0, k0, alpha0, beta0)
    
    log_w[t_index, i] <- log_w[t_index - 1, i] + (log_w_curr - log_w_prev)
    
# Metropolis–Hastings step 
    
    mu <- theta_prev[1]
    tau <- theta_prev[2]
    
    # Propose new state
    mu_prop <- rnorm(1, mu, 0.5)
    tau_prop <- abs(rnorm(1, tau, 0.1))  # keep tau > 0
    
    theta_prop <- c(mu_prop, tau_prop)
    
    # Compute log densities at current temperature
    log_curr <- power_post_log(power_curr, theta_prev, x, m0, k0, alpha0, beta0)
    log_prop <- power_post_log(power_curr, theta_prop, x, m0, k0, alpha0, beta0)
    
    # Accept / reject
    if (log(runif(1)) < (log_prop - log_curr)) {
      theta_new <- theta_prop
    } else {
      theta_new <- theta_prev
    }
    
    # Store updated particle
    mu_samples[t_index, i] <- theta_new[1]
    tau_samples[t_index, i] <- theta_new[2]
  }
}

# The final estimate is log(E_i[w_i]) using log-sum-exp for stability.


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


