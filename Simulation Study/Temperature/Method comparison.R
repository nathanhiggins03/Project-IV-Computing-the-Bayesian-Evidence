#Need to change to incorporate more efficient code for HME
library(ggplot2)

Sim <- 30
T_values <- c(50)   # <-- choose T

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

for (T_sample in T_values) {
  
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
    T <- T_sample             # Number of tempered distributions 
    Nsim <- 1000           # Number of AIS particles/samples
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
      
      # Log-scale correction
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
      iter = 100000,
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
        
        # ---- 1. Weight update (UNCHANGED) ----
        log_w_curr <- power_post_log(power_curr, theta_prev, X, y)
        log_w_prev <- power_post_log(power_prev, theta_prev, X, y)
        
        log_w[t_index, i] <- log_w[t_index - 1, i] + (log_w_curr - log_w_prev)
        
        # ---- 2. MH transitions (NEW, replaces Stan) ----
        theta_tmp <- theta_prev
        
        # Do multiple MH steps for better mixing
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
    #cat("\n--- AIS Results ---\n")
    #cat("Number of particles (N):", Nsim, "\n")
    # cat("Number of steps (T):", T, "\n")
    #cat("Estimated Log Evidence (log Z):", AIS_log_evidence, "\n")
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
      T = factor(paste0("T = ", T_sample),
                 levels = paste0("T = ", T_values))
    )
  )
  
}


df_ais <- df_all     # for AIS



#Need to change to incorporate more efficient code for HME
library(ggplot2)

Sim <- 30
T_values <- c(50)   # <-- choose T

df_all <- data.frame()

setwd("~/Desktop/Project IV")   # set working directory to Project IV folder

library(rstan)
library(durhamSLR)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())


prior_model <- stan_model("Prior Bayesian Linear Regression.stan")
power_model <- stan_model("Power Posterior Bayes Linear Regression.stan")


for (T_sample in T_values) {
  
  time_simulation <- rep(0, Sim)
  est_simulation  <- rep(0, Sim)
  
  for (k in 1:Sim) {
    
    start=proc.time()
    
    #Bayesian Linear Regression with unknown beta and precision
    
    #Set seed for reproducibility
    set.seed(123 + k)
    
    #Power Posterior Method
    
    #Following the 'review' paper
    
    #We set our temperature scale
    
    #Then at each temperature we : 
    #Sample from power posterior
    #Estimate Ej = expectation of log likelhood using samples from power posterior at temp tj
    
    #At end can use trapezoidal rule formula to estimate log evidence(need log sum exp trick for Ejs and log evidence calculations)
    
    #Now create power log likelihood function  for the Ej's
    
    power_like_log<- function(t, theta,X,y){
      out<-t*sum(dnorm(y, mean = X %*% theta[1:2], sd = sqrt(theta[3]), log = TRUE))
      return(out)
    }
    
    #N is number of (gibbs) samples for each tempered distribution
    Nsim<-1000
    #Sample from prior using STAN
    
    #Need lots of draws for it to work
    prior_fit <- sampling(
      prior_model,
      data = list(
        N = length(y),
        d = d,
        X = X,
        y = y,
        m0 = m0,
        Lambda0 = Lambda0,
        alpha0 = alpha0,
        beta0 = beta0
      ),
      iter = 100000,
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
    
    
    #T is number of tempered distributions
    T<- T_sample
    #N is number of (gibbs) samples for each tempered distribution
    Nsim<-1000
    
    #Temperature scale:
    
    #Tempering function- want uniform on log scale 
    Beta_t<-function(t, c=2){
      return((t/T)^c)
    }
    
    #Gives T+1 tempered posteriors(inlcuding prior(t=0) and posterior(t=1))
    t_list<-Beta_t(0:T, c=2)
    
    #Need to store outputs- samples of mu,tau in seperate T+1 x N matrices
    
    beta_samples<-array(0, dim = c(2,T + 1, Nsim))
    sigma_sq_samples<- matrix(0,nrow = T+1, ncol = Nsim)
    
    #Store initial samples in matrix
    beta_samples[1,1,]<-prior_sample_Beta1
    beta_samples[2,1,]<-prior_sample_Beta2
    sigma_sq_samples[1,]<-prior_sample_sigma_sq
    
    
    #Set E0 to expectation of log likelihood under prior
    
    E<- matrix(0,nrow = T+1, ncol = Nsim)
    for (i in 1:Nsim) {
      E[1,i]<- sum(dnorm(y, mean = X %*% beta_samples[,1,i], sd = sqrt(sigma_sq_samples[1,i]), log = TRUE))
      
    }
    
    
    N_samples_warmup <- 0
    N_samples_kept   <- Nsim
    
    #Now need to structure calculating weights then getting new samples in loops
    for(t in 2:(T+1)){
      power <- t_list[t]
      
      # Init: pick Nsim indices from previous temperature
      init_list <- lapply(1:1, function(j) {  # only 1 chain here
        idx <- sample(1:Nsim, 1)
        list(Beta = beta_samples[,t-1, idx],
             sigma_sq = sigma_sq_samples[t-1, idx])
      })
      
      run <- sampling(
        power_model,
        data = list(
          N = length(y),
          d = d,
          X = X,
          y = y,
          m0 = m0,
          Lambda0 = Lambda0,
          alpha0 = alpha0,
          beta0 = beta0,
          t = power
        ),
        iter = N_samples_kept,
        warmup = N_samples_warmup,
        chains = 1,          # make sure chains = 1
        init = init_list,
        refresh = 0
      )
      
      # Extract Beta and sigma_sq safely
      beta_post <- extract(run, pars = "Beta")$Beta
      if(is.null(dim(beta_post))) beta_post <- matrix(beta_post, nrow=2, ncol=1)
      beta_samples[,t,] <- t(beta_post)
      
      sigma_post <- extract(run, pars = "sigma_sq")$sigma_sq
      if(is.null(dim(sigma_post))) sigma_post <- matrix(sigma_post, nrow=1, ncol=1)
      sigma_sq_samples[t,] <- sigma_post
      
      # Compute E[t,]
      for(i in 1:Nsim){
        E[t,i] <- sum(dnorm(y, mean = X %*% beta_samples[,t,i], sd = sqrt(sigma_sq_samples[t,i]), log = TRUE))
      }
    }
    
    
    #Now work out power posterior estimate of log evidence
    #Need to apply log sum exp trick to this to get it more accurate
    power_post_est<-0
    for(t in 2:(T+1)){
      power_post_est<- power_post_est + (t_list[t] - t_list[t-1]) *(0.5*(rowMeans(E)[t] + rowMeans(E)[t-1]))
    }
    est_simulation[k]<-power_post_est
    
    time_simulation[k] <- (proc.time() - start)[3]
  }
  
  #Store results
  df_all <- rbind(
    df_all,
    data.frame(
      Estimate = est_simulation,
      Time     = time_simulation,   # ← ADD THIS LINE
      T = factor(paste0("T = ", T_sample),
                 levels = paste0("T = ", T_values))
    )
  )
  
}


df_pp <- df_all      # for TI




#Need to change to incorporate more efficient code for HME
library(ggplot2)

Sim <- 30
T_values <- c(5)   # <-- choose T

df_all <- data.frame()

setwd("~/Desktop/Project IV")   # set working directory to Project IV folder

library(rstan)
library(durhamSLR)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
prior_model <- stan_model("Prior Bayesian Linear Regression.stan")

for (T_sample in T_values) {
  
  time_simulation <- rep(0, Sim)
  est_simulation  <- rep(0, Sim)
  
  for (k in 1:Sim) {
    
    start=proc.time()
    
    #Bayesian Linear Regression with unknown beta and precision
    
    #Set seed for reproducibility
    set.seed(123 + k)
    
    #SMC
    
    #Now create log power posterior function in R
    library(mvtnorm)
    library(extraDistr)
    
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
    
    #T is number of tempered distributions
    T<- T_sample
    
    #N is number of (gibbs) samples for each tempered distribution
    Nsim<-1000
    
    #Sample from joint prior
    
    
    #Need lots of draws for it to work
    # Prior sampling
    prior_fit <- sampling(
      prior_model,
      data = list(
        N = length(y),
        d = d,
        X = X,
        m0 = m0,
        Lambda0 = Lambda0,
        alpha0 = alpha0,
        beta0 = beta0
      ),
      iter = 100000,          # number of draws you want
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
    log_w<-matrix(0,nrow = T, ncol = Nsim)
    
    #Tempering function- want uniform on log scale 
    Beta_t<-function(t, c=2){
      return((t/T)^c)
    }
    
    #Gives T+1 tempered posteriors(inlcuding prior(t=0) and posterior(t=1))
    t_list<-Beta_t(0:T, c=4)
    #t_list <- exp(seq(log(1e-6), log(1), length.out = T+1))
    
    
    
    
    #SMC loop- weight, normalise, resample
    for(t in 2:(T+1)){
      power<-t_list[t]
      for (i in 1:Nsim) {
        #Weight
        log_w[t-1,i]<- power_post_log(power, c(beta_samples[1,t-1,i], beta_samples[2,t-1,i],sigma_sq_samples[t-1,i]), X,y) - power_post_log(t_list[t-1], c(beta_samples[1,t-1,i], beta_samples[2,t-1,i],sigma_sq_samples[t-1,i]), X,y)
      }
      # Log sum exp trick
      max_logw <- max(log_w[t-1, ])
      w <- exp(log_w[t-1, ] - max_logw)
      #Explicit normalisation not needed- R can handle unnormalised weights
      #norm_w <- w / sum(w)
      #Would swap w for norm_w in resample below
      
      # Joint resampling
      resample_index <- sample(1:Nsim, size = Nsim, replace = TRUE, prob = w)
      beta_samples[,t, ]  <- beta_samples[,t-1, resample_index]
      sigma_sq_samples[t, ] <- sigma_sq_samples[t-1, resample_index]
      
    }
    
    #Calculating log evidence evidence
    
    smc_log <- 0
    for (t in 1:T) {
      max_logw <- max(log_w[t, ])
      smc_log <- smc_log + log(mean(exp(log_w[t, ] - max_logw))) + max_logw
    }
    
    est_simulation[k]<-smc_log
    
    time_simulation[k] <- (proc.time() - start)[3]
  }
  
  #Store results
  df_all <- rbind(
    df_all,
    data.frame(
      Estimate = est_simulation,
      Time     = time_simulation,   # ← ADD THIS LINE
      T = factor(paste0("T = ", T_sample),
                 levels = paste0("T = ", T_values))
    )
  )
  
}


df_smc <- df_all     # for SMC


library(dplyr)


df_ais   <- df_ais   %>% mutate(Estimator = "AIS")
df_pp    <- df_pp    %>% mutate(Estimator = "Power Posterior")
df_smc   <- df_smc   %>% mutate(Estimator = "SMC")


df_compare <- bind_rows(
  df_ais,
  df_pp,
  df_smc
)

df_compare$Estimator <- factor(
  df_compare$Estimator,
  levels = c("AIS", "Power Posterior", "SMC")
)

library(ggplot2)

ggplot(df_compare, aes(x = Estimator, y = Estimate, fill = Estimator)) +
  geom_violin(trim = FALSE, alpha = 0.7, scale="width") +
  geom_boxplot(width = 0.12,
               fill = "white",
               outlier.shape = NA) +
  geom_hline(yintercept = true_le,
             colour = "red",
             linewidth = 0.8) +
  labs(
    title = "Evidence estimators using optimal T",
    y = "Log Evidence",
    x = ""
  ) +
  theme_minimal() +
  coord_cartesian(ylim = c(-70, NA))+
  theme(legend.position = "none")




df_efficiency <- df_compare %>%
  group_by(Estimator) %>%
  summarise(
    mean_time = mean(Time),
    mc_sd     = sd(Estimate),
    mc_bias   = mean(Estimate) - true_le,
    .groups = "drop"
  )

ggplot(df_efficiency,
       aes(x = mean_time, y = mc_sd, colour = Estimator)) +
  geom_point(size = 4) +
  labs(
    title = "Monte Carlo efficiency comparison",
    x = "Mean runtime (seconds)",
    y = "Monte Carlo SD",
    colour = "Estimator"
  ) +
  theme_minimal()

# RMSE vs mean runtime comparison

df_rmse <- df_compare %>%
  group_by(Estimator) %>%
  summarise(
    mean_time = mean(Time),
    rmse = sqrt(mean((Estimate - true_le)^2)),
    .groups = "drop"
  )

ggplot(df_rmse,
       aes(x = mean_time, y = rmse, colour = Estimator)) +
  geom_point(size = 5) +
  labs(
    title = "Efficiency comparison",
    x = "Mean runtime (seconds)",
    y = "RMSE of log evidence",
    colour = "Estimator"
  ) +
  theme_minimal()
