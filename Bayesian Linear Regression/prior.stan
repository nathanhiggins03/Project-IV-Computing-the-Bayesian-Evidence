data {
  int<lower=0> N;
    int<lower=1> d;
  matrix[N,d] X;
  vector[d] m0;
  real<lower=0> alpha0;
  real<lower=0> beta0;
  cov_matrix[d] Lambda0;
}
generated quantities {
  real<lower=0> sigma_sq;
  vector[d] Beta;
  vector[N] y;
  sigma_sq = inv_gamma_rng(alpha0,beta0);
  Beta = multi_normal_rng(m0,sigma_sq * inverse(Lambda0));
  y = multi_normal_rng(X * Beta, sigma_sq * diag_matrix(rep_vector(1, N)));
}
