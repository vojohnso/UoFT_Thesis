//
// This Stan program defines a simple model, with a
// vector of values 'y' modeled as normally distributed
// with mean 'mu' and standard deviation 'sigma'.
//
// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//

data {
  int N; // number of observations
  vector[N] L;
  vector[N] A;
  vector[N] Y;
}

parameters {
  real beta0;
  real beta1;
  real beta2;
  real<lower=0> sigma;
}

model {
  vector[N] mu_y = beta0 + beta1 * A + beta2 * L;
  Y ~ normal(mu_y, sigma);
  // specify priors
  beta0 ~ normal(0, 1);
  beta1 ~ normal(0, 1);
  beta2 ~ normal(0, 1);
  sigma ~ normal(0, 1);
}

generated quantities {
  real mean_ATE;
  {
    vector[N] mu_1 = beta0 + beta1 * 1 + beta2 * L;
    vector[N] mu_0 = beta0 + beta1 * 0 + beta2 * L;
    vector[N] Y1 = to_vector(normal_rng(mu_1, sigma));
    vector[N] Y0 = to_vector(normal_rng(mu_0, sigma));
    mean_ATE = mean(Y1 - Y0);
  }
}
