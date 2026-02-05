#Need to change to incorporate more efficient code for HME
library(ggplot2)

Sim <- 30
T_values <- c(5,10,50,100)   # <-- choose T

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
    Nsim<-500
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
    Nsim<-500
    
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


#Plotting evidence violin plots against MC sample size
ggplot(df_all, aes(x = T, y = Estimate)) +
  geom_violin(aes(fill = T),
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
    title = "Effect of the number of temperatures on the Power posterior estimator",
    x = "T",
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
  group_by(T) %>%
  summarise(
    mean_time = mean(Time),
    sd_time   = sd(Time),
    .groups = "drop"
  )

#ggplot(df_time, aes(x = T, y = mean_time, group = 1)) +
#  geom_point(size = 3) +
#  geom_line() +
#  labs(
#    title = "Mean Runtime vs number of temperatures",
#    x = "T",
#    y = "Mean runtime (seconds)"
#  ) +
#  theme_minimal()



#Plotting MC standard deviation against run time
df_error <- df_all %>%
  group_by(T) %>%
  summarise(
    mc_sd = sd(Estimate),
    mc_bias = mean(Estimate) - true_le
  )
df_efficiency <- left_join(df_time, df_error, by = "T")
ggplot(df_efficiency, aes(x = mean_time, y = mc_sd)) +
  geom_point(size = 3) +
  geom_line()+
  geom_text(aes(label = T), vjust = -0.7, hjust = -0.2) +
  labs(
    title = "Effect of the number of temperatures on efficiency - Power posterior",
    x = "Mean runtime (seconds)",
    y = "Monte Carlo SD"
  ) +
  coord_cartesian(xlim = c(NA, 2))+
  theme_minimal()

#RMSE vs mean runtime (Power posterior – varying number of temperatures)

df_rmse <- df_all %>%
  group_by(T) %>%
  summarise(
    mean_time = mean(Time),
    rmse = sqrt(mean((Estimate - true_le)^2)),
    .groups = "drop"
  )

ggplot(df_rmse, aes(x = mean_time, y = rmse)) +
  geom_point(size = 3) +
  geom_line() +
  geom_text(aes(label = T), vjust = -0.6, hjust = -0.2) +
  labs(
    title = "Effect of the number of temperatures on efficiency - Power posterior",
    x = "Mean runtime (seconds)",
    y = "RMSE of log evidence"
  ) +
  coord_cartesian(xlim = c(NA, 2)) +
  theme_minimal()
