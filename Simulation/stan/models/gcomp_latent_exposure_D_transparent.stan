// Latent exposure model with the Option A error model, written in a "transparent" parameterization
// (Gustafson): the latent exposure is T = beta_Z1 * A, the exposure on the scale of Z. Everything the data can
// identify is a parameter of T; the one unidentified quantity, the scale beta_Z1, appears only in generated quantities.
//   T_i = m0 + X_i m + noise(s_A),   Z_i = T_i + beta_Z2 * D_i + noise(sigma_Z),   Y_i ~ Poisson(pop_i exp(beta0 + kappa T_i + X_i beta_X))
// The original parameters are recovered as A = T / beta_Z1, beta_A = kappa * beta_Z1, gamma = m / beta_Z1, sigma_A = s_A / beta_Z1.
// No spatial term here (a CAR effect on T would be added to the mean of T exactly as before).
data {
  int<lower=1> N;                // number of areas
  int<lower=1> G;                // number of exposure grid values
  array[N] int<lower=0> Y;       // event counts
  vector<lower=1>[N] pop;        // area populations (the offset)
  matrix[N, 4] X;                // confounders
  vector[N] Z;                   // error-prone exposure (surrogate)
  vector[N] D;                   // covariate that shifts the surrogate (mean 0, variance 1)
  vector[G] a_grid;              // exposure values (on the scale of A) at which theta(a) is evaluated
}

transformed data {
  vector[N] log_pop = log(pop);
  real mean_Z = mean(Z);
}

parameters {
  // identified by the data
  real beta0;
  real kappa;                    // log-rate change per unit of T (= beta_A / beta_Z1)
  vector[4] beta_X;
  real m0;                       // mean of T (= beta_Z1 * gamma0)
  vector[4] m;                   // effect of X on T (= beta_Z1 * gamma)
  real<lower=0> s_A;             // residual SD of T (= beta_Z1 * sigma_A)
  real beta_Z2;                  // shift of Z per unit of D
  real<lower=0> sigma_Z;         // SD of the measurement error
  vector[N] T;                   // latent exposure on the Z scale
  // not identified by the data: only its prior matters
  real<lower=0> beta_Z1;         // calibration slope of Z on A, i.e. the scale of A
}

model {
  // priors for the identified parameters
  beta0   ~ normal(0, 1);
  kappa   ~ normal(0, 1);
  beta_X  ~ normal(0, 1);
  m0      ~ normal(mean_Z, 10);
  m       ~ normal(0, 1);
  s_A     ~ cauchy(0, 1);
  beta_Z2 ~ normal(0, 1);
  sigma_Z ~ normal(0, 2);        // half-normal through the lower bound

  // the sensitivity assumption: nothing in the likelihood depends on beta_Z1
  beta_Z1 ~ normal(1, 0.5);

  T ~ normal(m0 + X * m, s_A);                          // exposure model (on the Z scale)
  Z ~ normal(T + beta_Z2 * D, sigma_Z);                 // measurement model
  Y ~ poisson_log(log_pop + beta0 + kappa * T + X * beta_X);
}

generated quantities {
  vector[G] theta;               // theta(a) at each grid value of A: depends on the unidentified beta_Z1
  real beta_A = kappa * beta_Z1;
  real sigma_A = s_A / beta_Z1;
  {
    vector[N] w = dirichlet_rng(rep_vector(1, N));   // Bayesian bootstrap weights over areas
    vector[N] xb = X * beta_X;
    for (g in 1:G) {
      theta[g] = dot_product(w, exp(beta0 + kappa * beta_Z1 * a_grid[g] + xb));
    }
  }
}
