// Berkson latent exposure model: the true exposure is the observed value plus independent noise, A = Z + N,
// N ~ normal(0, sigma_N), with the CAR effect phi on the outcome (sparse form, as in the classical reduced model).
// There is no exposure model for A and no slope to fix: A | Z is normal around Z, so the only error parameter is sigma_N.
// sigma_N is weakly identified (only through curvature of the exposure effect), so its prior is set by the data
// (sigma_N_mean, sigma_N_sd), like sigma_Z in the classical reduced model. A is non-centred: A = Z + sigma_N * A_raw.

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
  vector[G] a_grid;
  matrix[N, N] W;                // adjacency matrix
  real sigma_N_mean;             // prior for sigma_N (truncated at 0)
  real<lower=0> sigma_N_sd;
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
  real beta_A;
  vector[4] beta_X;
  vector[N] phi;
  real<lower=0> sigma_phi;
  real<lower=0, upper=1> rho;
  real<lower=0> sigma_N;
  vector[N] A_raw;
}

transformed parameters {
  vector[N] A = Z + sigma_N * A_raw;
}

model {
  beta0 ~ normal(0, 10);
  beta_A ~ normal(0, 10);
  beta_X ~ normal(0, 10);
  sigma_phi ~ cauchy(0, 1);
  rho ~ beta(2, 1);
  sigma_N ~ normal(sigma_N_mean, sigma_N_sd);  // weakly identified: the prior drives it

  phi ~ sparse_car(inv_square(sigma_phi), rho, W_sparse, D_sparse, lambda, N, W_n);
  A_raw ~ std_normal();                        // A | Z ~ normal(Z, sigma_N)
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
