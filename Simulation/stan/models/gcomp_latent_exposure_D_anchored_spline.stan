// Anchored latent exposure model (see gcomp_latent_exposure_D_anchored.stan) with a natural cubic spline in A.
// The basis is computed in Stan from the latent A; knots are fixed data on the standardised scale u = (A - centre) / scale.
functions {
  // natural cubic spline basis (truncated-power form), K knots give K - 1 columns, linear beyond the boundary knots
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
}

data {
  int<lower=1> N;
  int<lower=1> G;
  int<lower=3> K;                // number of knots
  array[N] int<lower=0> Y;
  vector<lower=1>[N] pop;
  matrix[N, 4] X;
  vector[N] Z;
  vector[N] D;
  vector[G] a_grid;
  vector[K] knots;               // on the standardised scale
  real centre;
  real<lower=0> scale;
  matrix[N, N] W;
}

transformed data {
  vector[N] log_pop = log(pop);
  real mean_Z = mean(Z);
  vector[N] zeros = rep_vector(0, N);
  matrix[N, N] D_w = diag_matrix(to_vector(W * rep_vector(1, N)));
  row_vector[K - 1] B0 = ns_basis([0]', knots)[1];
  matrix[G, K - 1] Bg = ns_basis((a_grid - centre) / scale, knots);
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

  beta0 ~ normal(0, 1);
  beta_A ~ normal(0, 1);
  beta_X ~ normal(0, 1);
  gamma0 ~ normal(mean_Z, 10);
  gamma_A ~ normal(0, 1);
  sigma_A ~ cauchy(0, 1);
  sigma_phi ~ cauchy(0, 1);
  rho ~ beta(2, 1);
  beta_Z2 ~ normal(0, 1);
  sigma_Z ~ normal(1, 0.3);

  phi ~ multi_normal_prec(zeros, (D_w - rho * W) / square(sigma_phi));
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
