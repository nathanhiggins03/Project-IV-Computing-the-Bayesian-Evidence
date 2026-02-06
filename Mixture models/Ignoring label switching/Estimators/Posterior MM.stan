data {
  int<lower=1> N;            // number of observations
  int<lower=1> K;            // number of mixture components
  vector[N] y;               // data

  // Hyperparameters
  vector[K] alpha;           // Dirichlet prior parameters
  real mu0;                  // prior mean for component means
  real<lower=0> lambda0;     // prior precision scaling
  real<lower=0> a0;          // Inv-Gamma shape
  real<lower=0> b0;          // Inv-Gamma scale
}

parameters {
  simplex[K] omega;          // mixture weights
  vector[K] mu;              // component means
  vector<lower=0>[K] sigma2; // component variances
}

model {
  // ---- Priors
  omega ~ dirichlet(alpha);
  sigma2 ~ inv_gamma(a0, b0);
  mu ~ normal(mu0, sqrt(sigma2 / lambda0));

  // ---- Marginalized mixture likelihood
  for (n in 1:N) {
    vector[K] lps;
    for (k in 1:K) {
      lps[k] =
        log(omega[k]) +
        normal_lpdf(y[n] | mu[k], sqrt(sigma2[k]));
    }
    target += log_sum_exp(lps);
  }
}

