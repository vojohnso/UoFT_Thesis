
// Latent exposure with CAR effect phi on the outcome; error model Z = A + beta_Z2 * D + noise (slope (beta_Z1) fixed at 1),
// sigma_Z ~ normal(1, 0.3), sigma_A >= 0.5. Initialise A = Z in R.
// Comments in the priors give the true values from data.generate() and whether the data identify the parameter
// (see Additional/identified_combinations.tex).

data {
  int<lower=1> N;
  int<lower=1> G;
  array[N] int<lower=0> Y;
  vector<lower=1>[N] pop;
  matrix[N, 4] X;
  vector[N] Z;
  vector[N] D;                   // covariate that shifts the surrogate (mean 0, variance 1)
  vector[G] a_grid;
  matrix[N, N] W;                // adjacency matrix
}

transformed data {
  vector[N] log_pop = log(pop);
  real mean_Z = mean(Z);
  vector[N] zeros = rep_vector(0, N);
  matrix[N, N] D_w = diag_matrix(to_vector(W * rep_vector(1, N)));   // number of neighbours per area
}

parameters {
  real beta0;
  real beta_A;
  vector[4] beta_X;
  real gamma0;
  vector[4] gamma_A;
  real<lower=0.5> sigma_A;
  vector[N] phi;
  real<lower=0> sigma_phi;
  real<lower=0, upper=1> rho;
  real beta_Z2;
  real<lower=0> sigma_Z;
  vector[N] A;
}

model {
  beta0 ~ normal(0, 10);                 // true -4, depends on sigma_Z
  beta_A ~ normal(0, 10);                // true 0.1, depends on sigma_Z
  beta_X ~ normal(0, 10);                // true (-0.5, -0.25, 0.25, 0.5), depends on sigma_Z
  gamma0 ~ normal(mean_Z, 10);           // true 10, identified
  gamma_A ~ normal(0, 10);               // true (-0.75, 0.75, -0.75, 0.75), identified
  sigma_A ~ cauchy(0, 1);                // true 1.5, not identified apart from sigma_Z
  sigma_phi ~ cauchy(0, 1);              // true 0.446 at 100 areas, weakly identified
  rho ~ beta(2, 1);                      // true 0.9, weakly identified
  beta_Z2 ~ normal(0, 10);               // true 0.5, identified
  sigma_Z ~ normal(1, 0.3);              // true 1, not identified apart from sigma_A

  phi ~ multi_normal_prec(zeros, (D_w - rho * W) / square(sigma_phi));   // CAR
  A ~ normal(gamma0 + X * gamma_A, sigma_A);
  Z ~ normal(A + beta_Z2 * D, sigma_Z);
  Y ~ poisson_log(log_pop + beta0 + beta_A * A + X * beta_X + phi);
}

generated quantities {
  vector[G] theta;
  {
    vector[N] w = dirichlet_rng(rep_vector(1, N));
    vector[N] xb = X * beta_X;
    for (g in 1:G) {
      theta[g] = dot_product(w, pop .* exp(beta0 + beta_A * a_grid[g] + xb + phi));
    }
  }
}
