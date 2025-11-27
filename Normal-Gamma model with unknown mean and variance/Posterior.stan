data {
  int<lower=0> N;
  vector[N] x;
  real<lower=0> k0;
  real m0;
  real<lower=0> alpha0;
  real<lower=0> beta0;
}

parameters {
  real mu;
  real<lower=0> tau;
}


model {
  tau ~ gamma(alpha0,beta0);
  mu ~ normal(m0,sqrt(1/(k0*tau)));
  x ~ normal(mu, sqrt(1/tau));
}
