data {
  int<lower=0> N;
    int<lower=1> d;
  matrix[N,d] X;
  vector[d] m0;
  real<lower=0> alpha0;
  real<lower=0> beta0;
  cov_matrix[d] Lambda0;
  vector[N] y;
  real<lower=0,upper=1> t;
}

parameters {
 real<lower=0> sigma_sq;
 vector[d] Beta;
}


model {
  sigma_sq ~ inv_gamma(alpha0,beta0);
  Beta ~ multi_normal(m0,sigma_sq * inverse(Lambda0));
  target += t * normal_lpdf(y | X * Beta, sqrt(sigma_sq) );
}
