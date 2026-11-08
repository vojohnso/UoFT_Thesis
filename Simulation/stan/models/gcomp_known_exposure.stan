// Outcome model with the exposure A observed (also the naive fit, with Z as the exposure); spline in A
data {
  int<lower=1> N;                // areas
  int<lower=1> K;                // spline degrees of freedom
  int<lower=1> G;                // exposure grid values
  array[N] int<lower=0> Y;       // event counts
  vector<lower=1>[N] pop;        // populations (offset)
  matrix[N, 4] X;                // confounders
  matrix[N, K] B;                // spline basis at each area's exposure
  matrix[G, K] Bg;               // spline basis at the grid values
}

transformed data {
  vector[N] log_pop = log(pop);
}

parameters {
  real beta0;
  vector[K] beta_A;
  vector[4] beta_X;
}

model {
  beta0 ~ normal(0, 10);
  beta_A ~ normal(0, 10);
  beta_X ~ normal(0, 10);

  Y ~ poisson_log(log_pop + beta0 + B * beta_A + X * beta_X);
}

generated quantities {
  vector[G] theta;               // expected cases per area if every area had exposure a
  {
    vector[N] w = dirichlet_rng(rep_vector(1, N));   // Bayesian bootstrap weights
    vector[N] xb = X * beta_X;
    for (g in 1:G) {
      theta[g] = dot_product(w, pop .* exp(beta0 + Bg[g] * beta_A + xb));
    }
  }
}
