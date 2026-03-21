#UPDATE SMC WITH NEW T AND GENERAL NEW CHANGES AS WELL
#Need to change to incorporate more efficient code for HME

library(ggplot2)

Sim <- 30
MC_values <- c(25,150,500,5000)   # <-- choose MC sizes

df_all <- data.frame()

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

#Data
set.seed(123)
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

#True value
#Analytical log evidence

term1 <- as.numeric(0.5 * (determinant(Lambda0, log = TRUE)$modulus - determinant(LambdaN, log = TRUE)$modulus))
term2 <- alpha0 * log(beta0) - alphaN * log(betaN)
term3 <- log(gamma(alphaN)) - log(gamma(alpha0))
term4 <- -(N/2) * log(2*pi)
true_le<- term1 + term2 + term3 + term4
true_le

for (MC_sample in MC_values) {
  
  time_simulation <- rep(0, Sim)
  est_simulation  <- rep(0, Sim)
  
  for (k in 1:Sim) {
    
    start=proc.time()
    #Adapt for Bayes Linear Regression
    #Set seed for reproducibility
    set.seed(123 + k)
    
    library(mvtnorm)
    library(extraDistr)
    
    
    
    
    # AIS Parameters
    T <- 25             # Number of tempered distributions 
    Nsim <- MC_sample           # Number of AIS particles/samples
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
    
    mh_step <- function(theta, beta_temp, X, y, m0, Lambda0, alpha0, beta0) {
      Beta <- theta[1:2]
      sigma_sq <- theta[3]
      
      # Proposals
      Beta_prop <- Beta + rnorm(2, 0, 0.2)
      sigma_sq_prop <- sigma_sq * exp(rnorm(1, 0, 0.1))
      
      theta_prop <- c(Beta_prop, sigma_sq_prop)
      
      # Log densities
      log_curr <- power_post_log(beta_temp, theta, X, y)
      log_prop <- power_post_log(beta_temp, theta_prop, X, y)
      
      # Correct MH ratio (log-scale correction)
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
      iter = 10*MC_sample,
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
    N_samples_kept <- Nsim          # 2000 samples kept (for weight calculation)
    N_samples_warmup <- 4000     # 4000 iterations for warmup
    
    for(t_index in 2:(T+1)){
      power_curr <- t_list[t_index]
      power_prev <- t_list[t_index - 1]
      
      for (i in 1:Nsim) {
        
        theta_prev <- c(beta_samples[,t_index - 1, i], 
                        sigma_sq_samples[t_index - 1, i])
        
        # ---- 1. Weight update ----
        log_w_curr <- power_post_log(power_curr, theta_prev, X, y)
        log_w_prev <- power_post_log(power_prev, theta_prev, X, y)
        
        log_w[t_index, i] <- log_w[t_index - 1, i] + (log_w_curr - log_w_prev)
        
        # ---- 2. MH transitions (IMPROVED: multiple steps) ----
        theta_tmp <- theta_prev
        for (s in 1:3) {
          theta_tmp <- mh_step(theta_tmp, power_curr, X, y, m0, Lambda0, alpha0, beta0)
        }
        
        # Store updated particle
        beta_samples[,t_index, i] <- theta_tmp[1:2]
        sigma_sq_samples[t_index, i] <- theta_tmp[3]
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
    
    est_simulation[k]<-AIS_log_evidence
    
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
    title = "Monte Carlo convergence of the AIS estimator",
    x = "N",
    y = "Log Evidence"
  ) +
  theme_minimal()+
  theme(
    legend.position = "right",
    legend.title = element_blank()
  )


#Mean runtime vs MC sample size
library(dplyr)

df_time <- df_all %>%
  group_by(MC) %>%
  summarise(
    mean_time = mean(Time),
    sd_time   = sd(Time),
    .groups = "drop"
  )

#ggplot(df_time, aes(x = MC, y = mean_time, group = 1)) +
#  geom_point(size = 3) +
#  geom_line() +
#  labs(
#    title = "Mean Runtime vs Monte Carlo Sample Size",
#    x = "Monte Carlo sample size",
#    y = "Mean runtime (seconds)"
#  ) +
#  theme_minimal()



#Plotting MC standard deviation against run time
df_error <- df_all %>%
  group_by(MC) %>%
  summarise(
    mc_sd = sd(Estimate),
    mc_bias = mean(Estimate) - true_le
  )
df_efficiency <- left_join(df_time, df_error, by = "MC")
ggplot(df_efficiency, aes(x = mean_time, y = mc_sd)) +
  geom_point(size = 3) +
  geom_line()+
  geom_text(aes(label = MC), vjust = -0.5, hjust=-0.1) +
  labs(
    title = "Monte Carlo efficiency - AIS",
    x = "Mean runtime (seconds)",
    y = "Monte Carlo SD"
  ) +
  coord_cartesian(xlim = c(NA, 13.5))+
  theme_minimal()



#RMSE vs mean run time

#RMSE vs mean runtime
df_rmse <- df_all %>%
  group_by(MC) %>%
  summarise(
    mean_time = mean(Time),
    rmse = sqrt(mean((Estimate - true_le)^2)),
    .groups = "drop"
  )

ggplot(df_rmse, aes(x = mean_time, y = rmse)) +
  geom_point(size = 3) +
  geom_line() +
  geom_text(
    aes(label = MC),
    vjust = 0.1,
    hjust = -0.5
  ) +
  labs(
    title = "Monte Carlo efficiency – AIS estimator",
    x = "Mean runtime (seconds)",
    y = "RMSE of log evidence"
  ) +
  coord_cartesian(xlim = c(NA, 100)) +
  theme_minimal() +
  theme(
    legend.position = "none"
  )
