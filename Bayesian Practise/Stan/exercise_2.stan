
data {
  int N; // number of observations
  int K; // number of strata
  vector[N] W;
  array[N] int<lower=1, upper=K> V;
  vector[N] A;
  array[N] int<lower=0, upper=1> Y;
  
  int n_v[K]; // length of sample size in each stratum
}

parameters {
  // outcome model parameters
  real beta_W; // confounder W effect
  vector[K] beta_V; // vector of confounder stratum effects
  
  real mu_A; // population averaged treatment effect
  real<lower=0> tau; // SD of how each strata differ from each other
  vector[K] beta_A; // each strata's own slope
}

model {
  vector[N] eta; // outcome model 
  beta_A ~ normal(mu_A, tau);
  for (i in 1:N) {
    eta[i] = beta_W*W[i] + beta_V[V[i]] + beta_A[V[i]]*A[i]; // V is categorical, W is cont
    Y[i] ~ bernoulli_logit(eta[i]);
  }
  
  // specify priors
  beta_V ~ normal(0, 1);
  beta_W ~ normal(0, 1);
  tau ~ normal(0, 1);
}

generated quantities {
  real mean_y1;
  real mean_y0;
  vector[K] odds_ratio;
  {
    for (v in 1:K) { // for each strata 
      int nv = n_v[v];
      vector[nv] bb_weights = dirichlet_rng(rep_vector(1.0, nv)); 
      mean_y1 = dot_product(bb_weights, inv_logit(beta_W*W[v] + beta_V[v] + beta_A[v]));
      mean_y0 = dot_product(bb_weights, inv_logit(beta_W*W[v] + beta_V[v]));
      odds_ratio[v] = (mean_y1 / (1- mean_y1)) / (mean_y0 / (1 - mean_y0));
    }
  }
}