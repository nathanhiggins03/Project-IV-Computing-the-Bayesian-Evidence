data {
  real<lower=0> k0;
  real m0;
  real<lower=0> alpha0;
  real<lower=0> beta0;
}
generated quantities {
  real mu;
  real x;
  real<lower=0> tau;
  tau= gamma_rng(alpha0,beta0);
  mu = normal_rng(m0,sqrt(1/(k0*tau)));
  x = normal_rng(mu, sqrt(1/tau));
}
