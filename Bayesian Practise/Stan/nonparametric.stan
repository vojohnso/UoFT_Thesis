
data {
  int N; // number of observations/subjects
  int P; // number of columns in design matrix (including A)
  matrix[N, P] X; 
  array[N] int<lower=0, upper=1> Y; // outcome
  int<lower=1> a_col ; // which column of X is treatment A
}

parameters {
  vector[P] beta; // parameters of all coefficients
}

model {
  Y ~ bernoulli_logit(X * beta);
  
  // specify priors
  beta ~ normal(0, 1);
}

generated quantities {
  vector[N] p1;
  vector[N] p0;
  real S1_marg;
  real S0_marg; 
  real psi; // risk difference: P(death | A=1) - P(death | A=0)
  vector[N] dirichlet_w;
  
  {
    matrix[N, P] X1 = X;
    matrix[N, P] X0 = X;
    for (i in 1:N) {
      X1[i, a_col] = 1;
      X0[i, a_col] = 0;
    }
    p1 = inv_logit(X1 * beta); 
    p0 = inv_logit(X0 * beta);
  }

  // Bayesian bootstrap: Dirichlet(1,...,1) weights drawn fresh every posterior draw
  dirichlet_w = dirichlet_rng(rep_vector(1.0, N));

  S1_marg = dot_product(dirichlet_w, p1);
  S0_marg = dot_product(dirichlet_w, p0);
  psi = S1_marg - S0_marg;
}
