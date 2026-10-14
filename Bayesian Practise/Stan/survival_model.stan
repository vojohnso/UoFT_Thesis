data {
  int<lower=0> N;                                 // number of subjects
  vector[N] X1;                                    // confounders
  vector[N] X2;
  array[N] int<lower=0, upper=1> is_cens;           // censoring indicator
  vector[N] T;                                      // event/censoring time (use T[i] when observed)
  vector[N] C;                                      // censoring time (use C[i] when censored)
  array[N] int<lower=0, upper=1> A;                 // treatment (binary)
}

parameters {
  // outcome model: logit... no, Weibull AFT location, depends on X1, X2, U, A
  real eta0;
  real eta_x1;
  real eta_x2;
  real eta_u;
  real eta_a;
  real<lower=0> k_D;                                // Weibull shape (event)

  // treatment model: logit(P(A=1)) depends on X1, X2, U
  real alpha0;
  real alpha1;
  real alpha2;
  real alpha3;

  // latent U model: logit(P(U=1)) depends on X1, X2
  real gamma0;
  real gamma1;
  real gamma2;
}

model {
  // priors
  eta0   ~ normal(0, sqrt(1/10.0));   // JAGS precision=10  -> sd = sqrt(1/10)
  eta_x1 ~ normal(0, sqrt(1/4.0));    // precision=4  -> sd = sqrt(1/4)
  eta_x2 ~ normal(0, sqrt(1/4.0));
  eta_a  ~ normal(0, sqrt(1/4.0));
  eta_u  ~ normal(0, sqrt(1/1.0));

  gamma0 ~ normal(0, sqrt(1/5.0));
  gamma1 ~ normal(0, sqrt(1/1.0));
  gamma2 ~ normal(0, sqrt(1/1.0));

  alpha0 ~ normal(0, sqrt(1/10.0));
  alpha1 ~ normal(0, sqrt(1/4.0));
  alpha2 ~ normal(0, sqrt(1/4.0));
  alpha3 ~ normal(0, sqrt(1/1.0));

  k_D ~ gamma(9, 6);   // event shape

  // likelihood: marginalize the discrete latent U out of each subject's contribution
  for (i in 1:N) {

    // --- branch U_i = 1 ---
    real logit_pU1 = gamma0 + gamma1 * X1[i] + gamma2 * X2[i];
    real lp_U1 = bernoulli_logit_lpmf(1 | logit_pU1);

    real logit_pA_given_U1 = alpha0 + alpha1 * X1[i] + alpha2 * X2[i] + alpha3 * 1;
    real lp_A_given_U1 = bernoulli_logit_lpmf(A[i] | logit_pA_given_U1);

    real mu_i_U1 = eta0 + eta_x1 * X1[i] + eta_x2 * X2[i] + eta_u * 1 + eta_a * A[i];
    real sigma_i_U1 = exp(mu_i_U1);
    real lp_T_given_U1 = is_cens[i] == 1
                          ? weibull_lccdf(C[i] | k_D, sigma_i_U1)
                          : weibull_lpdf(T[i] | k_D, sigma_i_U1);

    real branch_U1 = lp_U1 + lp_A_given_U1 + lp_T_given_U1;

    // --- branch U_i = 0 ---
    real logit_pU0 = gamma0 + gamma1 * X1[i] + gamma2 * X2[i];
    real lp_U0 = bernoulli_logit_lpmf(0 | logit_pU0);

    real logit_pA_given_U0 = alpha0 + alpha1 * X1[i] + alpha2 * X2[i] + alpha3 * 0;
    real lp_A_given_U0 = bernoulli_logit_lpmf(A[i] | logit_pA_given_U0);

    real mu_i_U0 = eta0 + eta_x1 * X1[i] + eta_x2 * X2[i] + eta_u * 0 + eta_a * A[i];
    real sigma_i_U0 = exp(mu_i_U0);
    real lp_T_given_U0 = is_cens[i] == 1
                          ? weibull_lccdf(C[i] | k_D, sigma_i_U0)
                          : weibull_lpdf(T[i] | k_D, sigma_i_U0);

    real branch_U0 = lp_U0 + lp_A_given_U0 + lp_T_given_U0;

    // marginalize: sum over the two possible values of U_i
    target += log_sum_exp(branch_U1, branch_U0);
  }
}

generated quantities {
  // posterior P(U_i = 1 | data) for each subject, recovered via Bayes' rule
  // using the same two branch log-probabilities (not sampled directly, since
  // U was never a parameter -- just algebra on the marginalized pieces).
  vector[N] p_U1;
  for (i in 1:N) {
    real logit_pU1 = gamma0 + gamma1 * X1[i] + gamma2 * X2[i];
    real lp_U1 = bernoulli_logit_lpmf(1 | logit_pU1);
    real logit_pA_given_U1 = alpha0 + alpha1 * X1[i] + alpha2 * X2[i] + alpha3 * 1;
    real lp_A_given_U1 = bernoulli_logit_lpmf(A[i] | logit_pA_given_U1);
    real mu_i_U1 = eta0 + eta_x1 * X1[i] + eta_x2 * X2[i] + eta_u * 1 + eta_a * A[i];
    real sigma_i_U1 = exp(mu_i_U1);
    real lp_T_given_U1 = is_cens[i] == 1
                          ? weibull_lccdf(C[i] | k_D, sigma_i_U1)
                          : weibull_lpdf(T[i] | k_D, sigma_i_U1);
    real branch_U1 = lp_U1 + lp_A_given_U1 + lp_T_given_U1;

    real logit_pU0 = gamma0 + gamma1 * X1[i] + gamma2 * X2[i];
    real lp_U0 = bernoulli_logit_lpmf(0 | logit_pU0);
    real logit_pA_given_U0 = alpha0 + alpha1 * X1[i] + alpha2 * X2[i] + alpha3 * 0;
    real lp_A_given_U0 = bernoulli_logit_lpmf(A[i] | logit_pA_given_U0);
    real mu_i_U0 = eta0 + eta_x1 * X1[i] + eta_x2 * X2[i] + eta_u * 0 + eta_a * A[i];
    real sigma_i_U0 = exp(mu_i_U0);
    real lp_T_given_U0 = is_cens[i] == 1
                          ? weibull_lccdf(C[i] | k_D, sigma_i_U0)
                          : weibull_lpdf(T[i] | k_D, sigma_i_U0);
    real branch_U0 = lp_U0 + lp_A_given_U0 + lp_T_given_U0;

    p_U1[i] = exp(branch_U1 - log_sum_exp(branch_U1, branch_U0));
  }
}
