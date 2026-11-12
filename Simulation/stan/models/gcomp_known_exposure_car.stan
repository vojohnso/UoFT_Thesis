// Outcome model with the exposure observed (A for the oracle fit, Z for the naive fit) and a 
// CAR spatial effect phi on the outcome
// (sparse form, Joseph 2016); spline in the exposure
// Same CAR prior as gcomp_latent_exposure_D_reduced_sparse.stan

functions {
  // log density of phi ~ CAR(precision tau * (D - alpha * W)), up to a constant (the sum of log D is dropped)
  real sparse_car_lpdf(vector phi, real tau, real alpha, array[,] int W_sparse, vector D_sparse, vector lambda, int n, int W_n) {
    row_vector[n] phit_D = (phi .* D_sparse)';
    row_vector[n] phit_W = rep_row_vector(0, n);
    vector[n] ldet_terms;
    for (i in 1:W_n) {
      phit_W[W_sparse[i, 1]] += phi[W_sparse[i, 2]];
      phit_W[W_sparse[i, 2]] += phi[W_sparse[i, 1]];
    }
    for (i in 1:n) ldet_terms[i] = log1m(alpha * lambda[i]);
    return 0.5 * (n * log(tau) + sum(ldet_terms) - tau * (phit_D * phi - alpha * (phit_W * phi)));
  }
}

data {
  int<lower=1> N;                // areas
  int<lower=1> K;                // spline degrees of freedom
  int<lower=1> G;                // exposure grid values
  array[N] int<lower=0> Y;       // event counts
  vector<lower=1>[N] pop;        // populations (offset)
  matrix[N, 4] X;                // confounders
  matrix[N, K] B;                // spline basis at each area's exposure
  matrix[G, K] Bg;               // spline basis at the grid values
  matrix[N, N] W;                // adjacency matrix
}

transformed data {
  vector[N] log_pop = log(pop);
  int W_n = to_int(sum(W) / 2);                 // number of neighbour pairs
  array[W_n, 2] int W_sparse;                   // neighbour pairs (i < j)
  vector[N] D_sparse = W * rep_vector(1, N);    // number of neighbours per area
  vector[N] lambda;                             // eigenvalues of D^-1/2 W D^-1/2
  {
    int counter = 1;
    for (i in 1:(N - 1)) {
      for (j in (i + 1):N) {
        if (W[i, j] == 1) {
          W_sparse[counter, 1] = i;
          W_sparse[counter, 2] = j;
          counter += 1;
        }
      }
    }
    lambda = eigenvalues_sym(quad_form(W, diag_matrix(inv_sqrt(D_sparse))));
  }
}

parameters {
  real beta0;
  vector[K] beta_A;
  vector[4] beta_X;
  vector[N] phi;
  real<lower=0> sigma_phi;
  real<lower=0, upper=1> rho;
}

model {
  beta0 ~ normal(0, 10);
  beta_A ~ normal(0, 10);
  beta_X ~ normal(0, 10);
  sigma_phi ~ cauchy(0, 1);
  rho ~ beta(2, 1);

  phi ~ sparse_car(inv_square(sigma_phi), rho, W_sparse, D_sparse, lambda, N, W_n);
  Y ~ poisson_log(log_pop + beta0 + B * beta_A + X * beta_X + phi);
}

generated quantities {
  vector[G] theta;               // expected cases per area if every area had exposure a
  {
    vector[N] w = dirichlet_rng(rep_vector(1, N));   // Bayesian bootstrap weights
    vector[N] xb = X * beta_X;
    for (g in 1:G) {
      theta[g] = dot_product(w, pop .* exp(beta0 + Bg[g] * beta_A + xb + phi));
    }
  }
}
