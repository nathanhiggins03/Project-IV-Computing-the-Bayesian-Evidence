data {
  int<lower=1> K;            // number of mixture components

  // Hyperparameters
  vector[K] alpha;           // Dirichlet parameters
  real mu0;                  // prior mean
  real<lower=0> lambda0;     // mean precision scaling
  real<lower=0> a0;          // Inv-Gamma shape
  real<lower=0> b0;          // Inv-Gamma scale
}

generated quantities {
  simplex[K] omega;          // mixture weights
  vector[K] mu;              // component means
  vector<lower=0>[K] sigma2; // component variances

  int z;                     // latent component
  real y;                    // prior predictive draw

  // ---- sample mixture weights
  omega = dirichlet_rng(alpha);

  // ---- sample component parameters
  for (k in 1:K) {
    sigma2[k] = inv_gamma_rng(a0, b0);
    mu[k] = normal_rng(mu0, sqrt(sigma2[k] / lambda0));
  }

  // ---- sample latent component
  z = categorical_rng(omega);

  // ---- sample observation
  y = normal_rng(mu[z], sqrt(sigma2[z]));
}
