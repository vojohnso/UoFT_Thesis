// Latent exposure model with the Option A error model, estimated error parameters, and a CAR spatial effect on the exposure

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
  // spatial effect on the exposure (CAR)
  vector[N] phi;
  real<lower=0> sigma_phi;
  real<lower=0, upper=1> rho;    // CAR spatial dependence (0 = independent areas, near 1 = strongly smooth)
  // error model Z = beta_Z1 * A + beta_Z2 * D + noise(sigma_Z)
  real beta_Z1;
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
  rho ~ uniform(0, 1);           // non-informative prior for rho (spatial correlation in CAR)
  
  // error-model parameters
  beta_Z1 ~ normal(1, 1);        // calibration slope of Z on A 
  beta_Z2 ~ normal(0, 1);        // shift of Z per unit of D
  sigma_Z ~ cauchy(0, 1);       

  // CAR prior on the spatial effect: phi ~ Normal(0, sigma_phi^2 (D - rho W)^-1)
  // multi_normal_prec uses precision instead of variance
  phi ~ multi_normal_prec(zeros, (D_w - rho * W) / square(sigma_phi));

  A ~ normal(gamma0 + X * gamma_A + phi, sigma_A);   // exposure model with the spatial effect
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
