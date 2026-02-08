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

  simplex[K] omega_sorted;    // <- moved here, top-level
  vector[K] mu_sorted;        // temporary sorted vectors
  vector[K] sigma2_sorted;
  array[K] int idx;           // sorting indices

  int z;                     // latent component
  real y;                     // prior predictive draw

  // ---- sample mixture weights
  omega = dirichlet_rng(alpha);

  // ---- sample component parameters
  for (k in 1:K) {
    sigma2[k] = inv_gamma_rng(a0, b0);
    mu[k] = normal_rng(mu0, sqrt(sigma2[k] / lambda0));
  }

  // ---- sort components to fix label switching
  idx = sort_indices_asc(mu);
  for (k in 1:K) {
    mu_sorted[k] = mu[idx[k]];
    sigma2_sorted[k] = sigma2[idx[k]];
    omega_sorted[k] = omega[idx[k]];
  }

  mu = mu_sorted;
  sigma2 = sigma2_sorted;
  omega = omega_sorted;

  // ---- sample latent component
  z = categorical_rng(omega);

  // ---- sample observation
  y = normal_rng(mu[z], sqrt(sigma2[z]));
}


