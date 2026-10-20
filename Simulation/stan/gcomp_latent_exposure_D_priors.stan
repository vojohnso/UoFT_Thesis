// Latent exposure model with the Option A error model, estimated error parameters, and a CAR spatial effect on the exposure.

data {
  int<lower=1> N;                // number of areas
  int<lower=1> G;                // number of exposure grid values
  array[N] int<lower=0> Y;       // event counts
  vector<lower=1>[N] pop;        // area populations (the offset)
  matrix[N, 4] X;                // confounders
  vector[N] Z;                   // error-prone exposure (surrogate)
  vector[N] D;                   // covariate that shifts the surrogate (mean 0, variance 1)
  vector[G] a_grid;              // exposure values at which theta(a) is evaluated
  // spatial structure
  matrix[N, N] W;                // adjacency matrix (1 if two areas are neighbours)
  real<lower=0, upper=1> rho;    // CAR dependence, known
}

transformed data {
  vector[N] log_pop = log(pop);
  real mean_Z = mean(Z);
  // CAR prior has: D being diagonal with each area's number of neighbours
  vector[N] zeros = rep_vector(0, N);
  matrix[N, N] D_mat = diag_matrix(to_vector(W * rep_vector(1, N)));
}

parameters {
  // outcome model
  real beta0;
  real beta_A;
  vector[4] beta_X;
  // exposure model
  real gamma0;
  vector[4] gamma;
  real<lower=0> sigma_A;
  // spatial effect on the exposure (CAR)
  vector[N] phi;
  real<lower=0> sigma_phi;
  // error model Z = beta_Z1 * A + beta_Z2 * D + noise(sigma_Z)
  real beta_Z1;
  real beta_Z2;
  real<lower=0> sigma_Z;
  // latent true exposures
  vector[N] A;
}

// Chains must be initialised with A near (Z - beta_Z2*D)/beta_Z1
model {
  // outcome parameters
  beta0   ~ normal(0, 1);
  beta_A  ~ normal(0, 1);
  beta_X  ~ normal(0, 1);
  // exposure-model parameters
  gamma0  ~ normal(mean_Z, 10);  // exposure mean centred on the surrogate mean
  gamma   ~ normal(0, 1);
  sigma_A ~ cauchy(0, 1);        // error exposure variance
  sigma_phi ~ cauchy(0, 1);
  // error-model parameters
  beta_Z1 ~ normal(1, 1);        // calibration slope of Z on A; not identified by the data (sets the scale of A)
  beta_Z2 ~ normal(0, 1);        // shift of Z per unit of D
  sigma_Z ~ cauchy(0, 1);       

  // CAR prior on the spatial effect: phi ~ Normal(0, sigma_phi^2 (D - rho W)^-1)
  phi ~ multi_normal_prec(zeros, (D_mat - rho * W) / square(sigma_phi));

  A ~ normal(gamma0 + X * gamma + phi, sigma_A);   // exposure model with the spatial effect
  Z ~ normal(beta_Z1 * A + beta_Z2 * D, sigma_Z);    // measurement model (Option A)
  Y ~ poisson_log(log_pop + beta0 + beta_A * A + X * beta_X);
}

generated quantities {
  vector[G] theta;               // theta(a) at each grid value: expected events per person if everyone had exposure a
  {
    vector[N] w = dirichlet_rng(rep_vector(1, N));   // Bayesian bootstrap weights over areas
    vector[N] xb = X * beta_X;
    for (g in 1:G) {
      // the latent A is not used here: the exposure is set to a_grid[g] for everyone
      theta[g] = dot_product(w, exp(beta0 + beta_A * a_grid[g] + xb));
    }
  }
}
