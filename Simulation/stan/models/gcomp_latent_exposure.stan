data {
  int<lower=1> N;                // number of areas
  int<lower=1> G;                // number of exposure grid values
  array[N] int<lower=0> Y;       // event counts
  vector<lower=1>[N] pop;        // area populations (the offset)
  matrix[N, 4] X;                // confounders
  vector[N] Z;                   // error-prone exposure (surrogate)
  real<lower=0> sigma_Z;         // known SD of the error in Z
  vector[G] a_grid;              // exposure values at which theta(a) is evaluated
}

transformed data {
  vector[N] log_pop = log(pop);
  real mean_Z = mean(Z);
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
  // latent true exposures
  vector[N] A;
}

// Chains must be initialised with A near Z
model {
  // outcome parameters
  beta0   ~ normal(0, 1);
  beta_A  ~ normal(0, 1);
  beta_X  ~ normal(0, 1);
  // exposure-model parameters
  // Sample A close to what Z is 
  gamma0  ~ normal(mean_Z, 10);
  gamma   ~ normal(0, 1);
  sigma_A ~ cauchy(0, 1); 

  A ~ normal(gamma0 + X * gamma, sigma_A);   // exposure model
  Z ~ normal(A, sigma_Z);                    // measurement model 
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
