library(ggplot2)

Sim <- 30
MC_values <- c(100000)   # <-- choose MC sizes

df_all <- data.frame()

library(rstan)
library(durhamSLR)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

setwd("~/Desktop/Project IV")

prior_model <- stan_model(
  file = "Prior Mixture model.stan"
)

for (MC_sample in MC_values) {
  
  time_simulation <- rep(0, Sim)
  est_simulation  <- rep(0, Sim)
  
  for (k in 1:Sim) {

    start=proc.time()
    #Set seed for reproducibility
    set.seed(123+k)
    
    #Data
    data(galaxies, package = "MASS")
    K <- 3   # or any K you want to test
    y<- galaxies/1000
    
    
    
    #Prior inputs
    N = length(y)
    K = K
    y = y
    
    alpha = rep(1, K)        # uniform Dirichlet
    mu0 = mean(y)     # data-centered prior
    lambda0 = 2.6/(max(y)-min(y))           # weak prior on means
    a0 = 1.28                   # weak Inv-Gamma
    b0 = 0.36*(mean(y^2) - (mean(y)^2))
    
    
    
    #Mean of likelihood using prior samples:
    
    #Sample from prior using STAN
    
    setwd("~/Desktop/Project IV")   # set working directory to Project IV folder
    
    
    
    #Need lots of draws for it to work
    prior_fit <- sampling(
      prior_model,
      data = list(
        N = N,
        K = K,
        y = y,
        
        alpha = alpha,        # uniform Dirichlet
        mu0 = mu0,     # data-centered prior
        lambda0 = lambda0,         # weak prior on means
        a0 = a0,                # weak Inv-Gamma
        b0 = b0
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
    prior_sample_sigma_sq<- extract(prior_fit, pars = 'sigma2')$'sigma2'  
    prior_sample_omega<- extract(prior_fit, pars = 'omega')$'omega'
    prior_sample_mu<- extract(prior_fit, pars = 'mu')$'mu'
    prior_sample_z<- extract(prior_fit, pars = 'z')$'z'
    
    #Calculate log likelihood at each prior sample (mu,tau)
    
    S <- dim(prior_sample_mu)[1]
    N <- length(y)
    K <- ncol(prior_sample_mu)
    
    likelihood_prior_log <- rep(0, S)
    
    for (s in 1:S) {
      
      likelihood_prior_log[s] <- sum(
        sapply(y, function(yi) {
          
          # mixture log density at yi
          m <- max(
            log(prior_sample_omega[s, ]) +
              dnorm(
                yi,
                mean = prior_sample_mu[s, ],
                sd   = sqrt(prior_sample_sigma_sq[s, ]),
                log  = TRUE
              )
          )
          
          m + log(sum(exp(
            log(prior_sample_omega[s, ]) +
              dnorm(
                yi,
                mean = prior_sample_mu[s, ],
                sd   = sqrt(prior_sample_sigma_sq[s, ]),
                log  = TRUE
              ) - m
          )))
        }))
    }
    
    
    #Apply log sum exp trick to log(mean(likelihood))
    m <- max(likelihood_prior_log)
    stan_prior_le <- m + log(mean(exp(likelihood_prior_log - m)))
    est_simulation[k]<-stan_prior_le
    
    end=proc.time()
    
    timer<- end-start
    time_simulation[k]<- timer[3]
    
    #Can compare prior method to "arithmetic mean" in paper for each K, know we are
    #reproducing correct results
    
    
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




#Plotting evidence violin plots against MC sample size
ggplot(df_all, aes(x = MC, y = Estimate)) +
  geom_violin(aes(fill = MC),
              trim = FALSE,
              alpha = 0.7) +
  geom_boxplot(width = 0.12,
               fill = "white",
               outlier.shape = NA) +
  scale_color_manual(
    name = ""
  ) +
  labs(
    title = "Monte Carlo convergence of the Prior method",
    x = "N",
    y = "Log Evidence"
  ) +
  theme_minimal()+
  theme(
    legend.position = "right",
    legend.title = element_blank(),
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
    mc_sd = sd(Estimate)
  )
df_efficiency <- left_join(df_time, df_error, by = "MC")
ggplot(df_efficiency, aes(x = mean_time, y = mc_sd)) +
  geom_point(size = 3) +
  geom_line()+
  geom_text(aes(label = MC), vjust = -0.4, hjust=-0.1) +
  labs(
    title = "Monte Carlo efficiency - Prior method",
    x = "Mean runtime (seconds)",
    y = "Monte Carlo SD"
  ) +
#  coord_cartesian(xlim = c(NA, 25))+
  theme_minimal()
