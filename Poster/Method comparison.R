library(ggplot2)

Sim <- 30
MC_values <- c(1000000)   # <-- choose MC sizes

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

df_prior <- df_all   # for Prior MC


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





#Need to change to incorporate more efficient code for HME
library(ggplot2)

Sim <- 30
MC_values <- c(100000)    # <-- choose MC sizes

df_all <- data.frame()

library(rstan)
setwd("~/Desktop/Project IV")   # set working directory to Project IV folder

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

posterior_model <- stan_model(
  file = "Posterior Bayesian Linear Regression.stan"
)


for (MC_sample in MC_values) {
  
  time_simulation <- rep(0, Sim)
  est_simulation  <- rep(0, Sim)
  
  for (k in 1:Sim) {
    start <- proc.time() 
    #inlcuding STAN setup in time measurement
    
    set.seed(123 + k)
    
    #Bayesian Linear Regression with unknown beta and precision
    
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
    
    
    #(Log) Harmonic mean estimator
    
    
    #Sample from posterior using STAN
    
    #Apply MCMC to model- get posterior mu samples
    library(rstan)
    library(durhamSLR)
    rstan_options(auto_write = TRUE)
    options(mc.cores = parallel::detectCores())
    
    
    # setwd("~/Desktop/Project IV")   # set working directory to Project IV folder
    
    #Stan Data
    stan_data<- list(N=N,d=d, X=X, m0=m0,alpha0=alpha0,beta0=beta0, Lambda0 = Lambda0, y=y)
    
    #posterior_sample <- stan(
    #  file = "Posterior Bayesian Linear Regression.stan",
    #  data = stan_data, iter = 2*MC_sample)
    posterior_sample <- sampling(
      posterior_model,
      data = stan_data,
      iter = 2 * MC_sample,
      chains = 1,
      refresh = 0
    )
    
    
    #Diagnostic checks
    #print(posterior_sample)
    #output = as.array(posterior_sample)
    #diagnostics(output)
    
    #Extract Prior samples
    post_sample_sigma_sq<- extract(posterior_sample, pars = 'sigma_sq')$'sigma_sq'  
    post_sample_Beta1<- extract(posterior_sample, pars = 'Beta[1]')$'Beta[1]'
    post_sample_Beta2<- extract(posterior_sample, pars = 'Beta[2]')$'Beta[2]'
    post_sample_Beta<- cbind(post_sample_Beta1,post_sample_Beta2)
    
    #Calculate likelihood at each posterior sample of (mu,tau)
    likelihood_log<- rep(0, length(post_sample_Beta1))
    for(i in 1:length(post_sample_Beta1)){
      mean<-as.vector(X %*% post_sample_Beta[i,])
      likelihood_log[i]<- sum(dnorm(y,mean = mean, sd = sqrt(post_sample_sigma_sq[i]) , log = TRUE))
    }
    
    #Calculate log harmonic mean using log sum exp trick for numerical stability
    m<-max(likelihood_log)
    hme_log<- log(length(post_sample_Beta1)) - log(exp(-m) * sum(exp(-likelihood_log + m)))
    
    est_simulation[k]<-hme_log
    
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

df_hme <- df_all     # for HME








#Need to change to incorporate more efficient code for HME
library(ggplot2)

Sim <- 30
MC_values <- c(10000)   # <-- choose MC sizes

df_all <- data.frame()

setwd("~/Desktop/Project IV")   # set working directory to Project IV folder

library(rstan)
library(durhamSLR)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())


prior_model <- stan_model("Prior Bayesian Linear Regression.stan")
power_model <- stan_model("Power Posterior Bayes Linear Regression.stan")


for (MC_sample in MC_values) {
  
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
    Nsim<-MC_sample
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
    
    
    #T is number of tempered distributions
    T<- 10
    #N is number of (gibbs) samples for each tempered distribution
    Nsim<-MC_sample
    
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
      MC = factor(paste0("N = ", MC_sample),
                  levels = paste0("N = ", MC_values))
    )
  )
  
}


df_pp <- df_all      # for TI


library(dplyr)

df_prior <- df_prior %>% mutate(Estimator = "Prior")
df_hme   <- df_hme   %>% mutate(Estimator = "HME")
df_pp    <- df_pp    %>% mutate(Estimator = "Power Posterior")


df_compare <- bind_rows(
  df_prior,
  df_hme,
  df_pp
)

df_compare$Estimator <- factor(
  df_compare$Estimator,
  levels = c("Prior", "HME","Power Posterior")
)

library(ggplot2)

ggplot(df_compare, aes(x = Estimator, y = Estimate, fill = Estimator)) +
  geom_violin(trim = FALSE, alpha = 0.7) +
  geom_boxplot(width = 0.12,
               fill = "white",
               outlier.shape = NA) +
  geom_hline(yintercept = true_le,
             colour = "red",
             linewidth = 1.2) +
  labs(
    title = "Comparing evidence estimators",
    y = "Log Evidence",
    x = ""
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text  = element_text(size = 18, colour = "black"),
    axis.title = element_text(size = 20, colour = "black"),
    plot.title = element_text(size = 20, colour = "black")
  )





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
  scale_x_continuous(limits = c(0, 40))+
  theme_minimal()+
  theme(
    legend.position = "right",
    axis.text  = element_text(size = 18, colour = "black"),
    axis.title = element_text(size = 20, colour = "black"),
    plot.title = element_text(size = 20, colour = "black"),
    legend.text  = element_text(size = 20, colour = "black"),
    legend.title = element_text(size = 20, colour = "black")
  )



#RMSE vs mean run time

df_efficiency <- df_compare %>%
  group_by(Estimator) %>%
  summarise(
    mean_time = mean(Time),
    bias      = mean(Estimate) - true_le,
    variance  = var(Estimate),
    rmse      = sqrt(mean((Estimate - true_le)^2)),
    .groups = "drop"
  )

ggplot(df_efficiency,
       aes(x = mean_time, y = rmse, colour = Estimator)) +
  geom_point(size = 4) +
  labs(
    title = "Monte Carlo efficiency comparison",
    x = "Mean runtime (seconds)",
    y = "RMSE of log evidence",
    colour = "Estimator"
  ) +
  scale_x_continuous(limits = c(0, 25)) +
  theme_minimal() +
  theme(
    legend.position = "right",
    axis.text  = element_text(size = 18, colour = "black"),
    axis.title = element_text(size = 20, colour = "black"),
    plot.title = element_text(size = 20, colour = "black"),
    legend.text  = element_text(size = 20, colour = "black"),
    legend.title = element_text(size = 20, colour = "black")
  )
