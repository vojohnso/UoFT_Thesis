data {
  int<lower=1> N;                // number of areas
  int<lower=1> K;                // spline degrees of freedom
  int<lower=1> G;                // number of exposure grid values
  array[N] int<lower=0> Y;       // event counts
  vector<lower=1>[N] pop;        // area populations (the offset)
  matrix[N, 4] X;                // confounders
  matrix[N, K] B;                // spline basis at each area's exposure
  matrix[G, K] Bg;               // spline basis at the grid values
}

transformed data {
  vector[N] log_pop = log(pop);
}

parameters {
  real beta0;                    // intercept
  vector[K] beta_A;              // coefficents of A
  vector[4] beta_X;              // coefficients of X
}

model {
  beta0 ~ normal(0, 1);
  beta_A ~ normal(0, 10);
  beta_X ~ normal(0, 1);
  
  Y ~ poisson_log(log_pop + beta0 + B * beta_A + X * beta_X);
}

generated quantities {
  vector[G] theta;               // theta(a) at each grid value: expected events per person if everyone had exposure a
  {
    vector[N] dirichlet_w = dirichlet_rng(rep_vector(1, N));
    vector[N] xb = X * beta_X;
    for (g in 1:G) {
      vector[N] rate = exp(beta0 + Bg[g] * beta_A + xb); 
      theta[g] = dot_product(dirichlet_w, rate);
    }
  }
}
