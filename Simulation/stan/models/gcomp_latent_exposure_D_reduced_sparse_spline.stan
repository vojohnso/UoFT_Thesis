// Reduced latent exposure model (beta_Z1 (SLOPE) = 1) with the sparse CAR prior and a natural cubic spline in the latent A
// The basis is built in Stan from A; knots are data on the standardised scale u = (A - centre) / scale.
// Latent exposure with CAR effect phi on the outcome; error model Z = A + beta_Z2 * D + noise,
// sigma_Z prior set by the data (sigma_Z_mean, sigma_Z_sd): wild guess normal(0, 5) or oracle normal(1, 0.3)
// Initialise A = Z in R.

functions {
  // natural cubic spline basis, K knots give K - 1 columns, linear beyond the boundary knots
  matrix ns_basis(vector x, vector knots) {
    int n = rows(x);
    int K = num_elements(knots);
    matrix[n, K - 1] B;
    vector[n] d_last;
    B[, 1] = x;
    d_last = (pow(fmax(x - knots[K - 1], 0), 3) - pow(fmax(x - knots[K], 0), 3)) / (knots[K] - knots[K - 1]);
    for (k in 1:(K - 2)) {
      B[, k + 1] = (pow(fmax(x - knots[k], 0), 3) - pow(fmax(x - knots[K], 0), 3)) / (knots[K] - knots[k]) - d_last;
    }
    return B;
  }
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
  int<lower=3> K;                // number of knots (K - 1 basis columns)
  vector[K] knots;               // on the standardised scale
  real centre;
  real<lower=0> scale;
  matrix[N, N] W;                // adjacency matrix
  real sigma_Z_mean;             // prior for sigma_Z (truncated at 0)
  real<lower=0> sigma_Z_sd;
}

transformed data {
  vector[N] log_pop = log(pop);
  real mean_Z = mean(Z);
  int W_n = to_int(sum(W) / 2);                 // number of neighbour pairs
  array[W_n, 2] int W_sparse;                   // neighbour pairs (i < j)
  vector[N] D_sparse = W * rep_vector(1, N);    // number of neighbours per area
  row_vector[K - 1] B0 = ns_basis([0]', knots)[1];   // basis at A = centre, subtracted so beta0 is the level there
  matrix[G, K - 1] Bg = ns_basis((a_grid - centre) / scale, knots);
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
  for (g in 1:G) Bg[g] -= B0;
}

parameters {
  real beta0;
  vector[K - 1] beta_A;
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
  matrix[N, K - 1] B = ns_basis((A - centre) / scale, knots);
  for (i in 1:N) B[i] -= B0;

  beta0 ~ normal(0, 10);                 // true -4, depends on sigma_Z
  beta_A ~ normal(0, 10);                // true 0.1, depends on sigma_Z
  beta_X ~ normal(0, 10);                // true (-0.5, -0.25, 0.25, 0.5), depends on sigma_Z
  gamma0 ~ normal(mean_Z, 10);           // true 10, identified
  gamma_A ~ normal(0, 10);               // true (-0.75, 0.75, -0.75, 0.75), identified
  sigma_A ~ cauchy(0, 1);                // true 1.5, not identified apart from sigma_Z
  sigma_phi ~ cauchy(0, 1);              // true 0.446 at 100 areas, weakly identified
  rho ~ beta(2, 1);                      // true 0.9, weakly identified
  beta_Z2 ~ normal(0, 10);               // true 0.5, identified
  sigma_Z ~ normal(sigma_Z_mean, sigma_Z_sd);  // true 1, not identified apart from sigma_A

  phi ~ sparse_car(inv_square(sigma_phi), rho, W_sparse, D_sparse, lambda, N, W_n);
  A ~ normal(gamma0 + X * gamma_A, sigma_A);
  Z ~ normal(A + beta_Z2 * D, sigma_Z);
  Y ~ poisson_log(log_pop + beta0 + B * beta_A + X * beta_X + phi);
}

generated quantities {
  vector[G] theta;
  {
    vector[N] w = dirichlet_rng(rep_vector(1, N));
    vector[N] xb = X * beta_X;
    for (g in 1:G) {
      theta[g] = dot_product(w, pop .* exp(beta0 + Bg[g] * beta_A + xb + phi));
    }
  }
}
