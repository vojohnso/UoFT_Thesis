data {
  int<lower=1> N;                // number of areas
  int<lower=1> G;                // number of exposure grid values
  array[N] int<lower=0> Y;       // event counts
  vector<lower=1>[N] pop;        // area populations (the offset, and the g-computation weights)
  matrix[N, 4] X;                // confounders
  vector[N] Z;                   // error-prone exposure (surrogate)
  real<lower=0> sigma_Z;         // known SD of the error in Z
  vector[G] a_grid;              // exposure values at which theta(a) is evaluated
}

transformed data {
  vector[N] log_pop = log(pop);
}

parameters {
  // outcome model (beta0 is the log rate at A = 0, i.e. no centering of the exposure)
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

// Chains must be initialised with A near Z (see run_latent_exposure.R); with Stan's default random starts
// (A in about (-2, 2) when exposure is near 10) the chains settle in different places and disagree.
model {
  // outcome parameters
  beta0   ~ normal(-3, 3);
  beta_A  ~ normal(0, 1);
  beta_X  ~ normal(0, 1);
  // exposure-model parameters
  gamma0  ~ normal(10, 10);
  gamma   ~ normal(0, 2);
  sigma_A ~ normal(0, 3);        // half-normal through the lower bound

  A ~ normal(gamma0 + X * gamma, sigma_A);   // exposure model
  Z ~ normal(A, sigma_Z);                    // measurement model
  Y ~ poisson_log(log_pop + beta0 + beta_A * A + X * beta_X);
}

generated quantities {
  vector[G] theta;               // theta(a) at each grid value
  {
    vector[N] w = dirichlet_rng(rep_vector(1, N)) .* pop;
    w = w / sum(w);
    vector[N] xb = X * beta_X;
    for (g in 1:G) {
      // the latent A is not used here: the exposure is set to a_grid[g] for everyone
      theta[g] = dot_product(w, exp(beta0 + beta_A * a_grid[g] + xb));
    }
  }
}
