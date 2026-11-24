// TIGHT-PRIOR diagnostic version of the reduced latent exposure model
// every prior is tight around the true value (oracle for all parameters).
// Reduced latent exposure model (beta_Z1 = 1 => fewer parameters to estimate) with the CAR prior in the sparse form 
// Latent exposure with CAR effect phi on the outcome; error model Z = A + beta_Z2 * D + noise (slope FIXED at 1)
// sigma_Z prior set by the data (sigma_Z_mean, sigma_Z_sd): wild guess normal(0, 5) or oracle normal(1, 0.3); sigma_A >= 0.5
// Initialise A = Z 

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
  beta0 ~ normal(-4, 0.1);   // true -4, depends on sigma_Z
  beta_A ~ normal(0.1, 0.1);   // true 0.1, depends on sigma_Z
  beta_X ~ normal([-0.5, -0.25, 0.25, 0.5]', 0.1); // true (-0.5, -0.25, 0.25, 0.5), depends on sigma_Z
  gamma0 ~ normal(10, 0.1);           // true 10, identified
  gamma_A ~ normal([-0.75, 0.75, -0.75, 0.75]', 0.1);  // true (-0.75, 0.75, -0.75, 0.75), identified
  sigma_A ~ normal(1.5, 0.1);    // true 1.5, not identified apart from sigma_Z
  sigma_phi ~ normal(0.446, 0.05);      // true 0.446 at 100 areas, weakly identified
  rho ~ beta(90, 10);       // true 0.9, weakly identified
  beta_Z2 ~ normal(0.5, 0.1);   // true 0.5, identified
  sigma_Z ~ normal(1, 0.05);  // true 1, not identified apart from sigma_A

  phi ~ sparse_car(inv_square(sigma_phi), rho, W_sparse, D_sparse, lambda, N, W_n);
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
