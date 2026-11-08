// Anchored version of gcomp_latent_exposure_D_priors.stan: the calibration slope is fixed at 1 (A is on the scale of Z),
// which removes the ridge between the slope and the scale of A. Estimated: beta_Z2, sigma_Z, plus the CAR spatial effect on the outcome

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
  matrix[N, N] W;                // adjacency matrix (1 if two areas are neighbours) - numbers of neighbours per area
}

transformed data {
  vector[N] log_pop = log(pop);
  real mean_Z = mean(Z);
  // CAR prior has: D being diagonal with each area's number of neighbours
  vector[N] zeros = rep_vector(0, N);
  matrix[N, N] D_w = diag_matrix(to_vector(W * rep_vector(1, N))); // diagonal degree matrix for num neighbours per area
}

parameters {
  // outcome model
  real beta0;
  real beta_A;
  vector[4] beta_X;
  // exposure model
  real gamma0;
  vector[4] gamma_A;
  real<lower=0> sigma_A;
  // spatial effect on the outcome log-rate (CAR)
  vector[N] phi;
  real<lower=0> sigma_phi;
  real<lower=0, upper=1> rho;    // CAR spatial dependence (0 = independent areas, near 1 = strongly smooth)
  // error model Z = A + beta_Z2 * D + noise(sigma_Z)   (slope on A anchored at 1)
  real beta_Z2;
  real<lower=0> sigma_Z;
  // latent true exposures
  vector[N] A;
}

// Initialize the chains with A near Z in the R script
model {
  // outcome parameters
  beta0 ~ normal(0, 1);
  beta_A ~ normal(0, 1);
  beta_X ~ normal(0, 1);
  
  // exposure-model parameters
  gamma0 ~ normal(mean_Z, 10);  // exposure mean centred on the surrogate mean
  gamma_A ~ normal(0, 1);
  sigma_A ~ cauchy(0, 1);        // error exposure variance
  
  // car prior parameters
  sigma_phi ~ cauchy(0, 1);
  rho ~ beta(2, 1);              // favours strong spatial dependence (median 0.71) instead of Uniform(0, 1), which implies weak correlation: see sim_rho_detectability.R
  
  // error-model parameters
  beta_Z2 ~ normal(0, 1);        // shift of Z per unit of D
  sigma_Z ~ cauchy(0, 1);       

  // CAR prior on the spatial effect in the outcome (log-rate): phi ~ Normal(0, sigma_phi^2 (D - rho W)^-1)
  // multi_normal_prec uses precision instead of variance
  phi ~ multi_normal_prec(zeros, (D_w - rho * W) / square(sigma_phi));

  A ~ normal(gamma0 + X * gamma_A, sigma_A);   // exposure model \\ Justify X as a confounder
  Z ~ normal(A + beta_Z2 * D, sigma_Z);    // measurement model (Option A, slope anchored at 1)
  Y ~ poisson_log(log_pop + beta0 + beta_A * A + X * beta_X + phi);  // non-linear specification? \\ interaction term A * X
}

generated quantities {
  vector[G] theta;               // theta(a) at each grid value: expected number of cases per area if every area had exposure a
  {
    vector[N] w = dirichlet_rng(rep_vector(1, N));   // Bayesian bootstrap weights over areas
    vector[N] xb = X * beta_X;
    for (g in 1:G) {
      // the latent A is not used here: the exposure is set to a_grid[g] for everyone
      theta[g] = dot_product(w, pop .* exp(beta0 + beta_A * a_grid[g] + xb + phi)); // do on counts or log scale (look at it compared to response)
    }
  }
}
// Yi/N = linear predictor 