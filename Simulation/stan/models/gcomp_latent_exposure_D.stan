// Latent exposure A, error model Z = beta_Z1 * A + beta_Z2 * D + noise with known error parameters. Initialise A = Z in R.
data {
  int<lower=1> N;
  int<lower=1> G;
  array[N] int<lower=0> Y;
  vector<lower=1>[N] pop;
  matrix[N, 4] X;
  vector[N] Z;
  vector[N] D;                   // covariate that shifts the surrogate (mean 0, variance 1)
  real beta_Z1;                  // known error parameters
  real beta_Z2;
  real<lower=0> sigma_Z;
  vector[G] a_grid;
}

transformed data {
  vector[N] log_pop = log(pop);
  real mean_Z = mean(Z);
}

parameters {
  real beta0;
  real beta_A;
  vector[4] beta_X;
  real gamma0;
  vector[4] gamma_A;
  real<lower=0> sigma_A;
  vector[N] A;
}

model {
  beta0 ~ normal(0, 10);
  beta_A ~ normal(0, 10);
  beta_X ~ normal(0, 10);
  gamma0 ~ normal(mean_Z, 10);
  gamma_A ~ normal(0, 10);
  sigma_A ~ cauchy(0, 1);

  A ~ normal(gamma0 + X * gamma_A, sigma_A);
  Z ~ normal(beta_Z1 * A + beta_Z2 * D, sigma_Z);
  Y ~ poisson_log(log_pop + beta0 + beta_A * A + X * beta_X);
}

generated quantities {
  vector[G] theta;
  {
    vector[N] w = dirichlet_rng(rep_vector(1, N));
    vector[N] xb = X * beta_X;
    for (g in 1:G) {
      theta[g] = dot_product(w, pop .* exp(beta0 + beta_A * a_grid[g] + xb));
    }
  }
}
