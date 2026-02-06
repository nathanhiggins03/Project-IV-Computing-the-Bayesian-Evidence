#Visualising data

library(ggplot2)

data(galaxies, package = "MASS")

y <- galaxies / 1000

ggplot(data.frame(y), aes(x = y)) +
  geom_histogram(aes(y = after_stat(density)),
                 bins = 30,
                 fill = "steelblue",
                 color = "black") +
  labs(
    title = "Galaxy velocities",
    x = "Velocity (km/s)",
    y = "Density"
  ) +
  theme_minimal()


#Inference ignoring label switching

#Posterior samples
set.seed(123)
K <- 3   # or any K you want to test
y<- galaxies/1000
stan_data <- list(
  N = length(y),
  K = K,
  y = y,
  
  alpha = rep(1, K),        # uniform Dirichlet
  mu0 = mean(y),     # data-centered prior
  lambda0 = 2.6/(max(y)-min(y)),           # weak prior on means
  a0 = 1.28,                   # weak Inv-Gamma
  b0 = 0.36*(mean(y^2) - (mean(y)^2))
)


library(rstan)

setwd("~/Desktop/Project IV")   # set working directory to Project IV folder

library(rstan)
library(durhamSLR)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
set.seed(123)
fit <- stan(
  file = "Posterior Mixture model.stan",
  data = stan_data,
  chains = 4,
  iter = 4000,
  warmup = 2000
)

#Diagnostic checks
print(fit)
output = as.array(fit,pars = "mu", include=TRUE)
diagnostics(output)



#Plots show label switching with k=3 for mu1 and mu3
library(bayesplot)

posterior <- as.array(fit, pars = "mu")

# Trace plots
mcmc_trace(posterior)

# Density plots
mcmc_dens(posterior)

mean(posterior[,,1])
sd(posterior[,,1])

mean(posterior[,,2])
sd(posterior[,,2])

mean(posterior[,,2])
sd(posterior[,,2])

library(bayesplot)
library(ggplot2)

posterior <- as.array(fit, pars = "mu")

trace_plot <- mcmc_trace(posterior)

trace_plot +
  facet_wrap(
    ~ parameter,
    labeller = as_labeller(c(
      "mu[1]" = expression(mu[1]),
      "mu[2]" = expression(mu[2]),
      "mu[3]" = expression(mu[3])
    ), label_parsed)
  ) +
  labs(
    title = expression(paste("Trace plots for ", mu)),
    x = "Iteration",
    y = "Value"
  )

library(bayesplot)

posterior <- as.array(fit, pars = "mu")

mcmc_dens(
  posterior,
  facet_args = list(
    labeller = as_labeller(
      c(
        "mu[1]" = expression(mu[1]),
        "mu[2]" = expression(mu[2]),
        "mu[3]" = expression(mu[3])
      ),
      label_parsed
    )
  ))+
  labs(
    title = expression(paste("Posterior Densities for ", mu)),
    x = "Value",
    y = "Density"
  )
  



#Inference with preprocessing model

#Fix label switching: Method 1 - relabelling

library(dplyr)
library(posterior)
library(bayesplot)
library(ggplot2)

draws <- as_draws_df(fit) %>%
  arrange(.chain, .iteration)

draws_mu_relabelled <- draws %>%
  group_by(.chain) %>%                 # <<< critical
  rowwise() %>%
  mutate(
    mu_ord = list(sort(c(`mu[1]`, `mu[2]`, `mu[3]`))),
    `mu[1]` = mu_ord[[1]],
    `mu[2]` = mu_ord[[2]],
    `mu[3]` = mu_ord[[3]]
  ) %>%
  select(.chain, .iteration, starts_with("mu[")) %>%  # <<< DROP EVERYTHING ELSE
  ungroup()



mu_array <- as_draws_array(draws_mu_relabelled)


variables(mu_array)
# should return: "mu[1]" "mu[2]" "mu[3]"

trace_plot <- mcmc_trace(mu_array)

trace_plot +
  facet_wrap(
    ~ parameter,
    labeller = as_labeller(
      c(
        "mu[1]" = expression(mu[1]),
        "mu[2]" = expression(mu[2]),
        "mu[3]" = expression(mu[3])
      ),
      label_parsed
    )
  ) +
  labs(
    title = expression(paste("Trace plots for ", mu)),
    x = "Iteration",
    y = "Value"
  )

mcmc_dens(
  mu_array,
  facet_args = list(
    labeller = as_labeller(
      c(
        "mu[1]" = expression(mu[1]),
        "mu[2]" = expression(mu[2]),
        "mu[3]" = expression(mu[3])
      ),
      label_parsed
    )
  )
) +
  labs(
    title = expression(paste("Posterior Densities for ", mu)),
    x = "Value",
    y = "Density"
  )

library(posterior)

mu_mat <- as_draws_matrix(mu_array)

apply(mu_mat, 2, mean)
apply(mu_mat, 2, sd)



#Inference with post processing data


