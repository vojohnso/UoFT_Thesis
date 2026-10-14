library(MASS)
library(spdep)

# ============================================================
# Spatial helpers
# ============================================================

# Build I x I rook contiguity adjacency matrix for a sqrt(I) x sqrt(I) grid
make_adj_matrix <- function(I) {
  g <- sqrt(I)
  if (g != floor(g)) stop("I must be a perfect square (e.g. 100, 196, 400)")
  nb <- cell2nb(nrow = g, ncol = g, type = "rook")
  W <- nb2mat(nb, style = "B", zero.policy = TRUE)
  W
}

# Draw one realisation of the CAR prior: phi ~ N(0, sigma2 * (D_W - rho*W)^-1)
sim_car <- function(W, rho, sigma2) {
  D <- diag(rowSums(W))
  Q <- D - rho * W
  Sigma <- sigma2 * solve(Q)
  as.numeric(mvrnorm(1, mu = rep(0, nrow(W)), Sigma = Sigma))
}

# Solve for sigma2_phi that achieves a target R2_phi = Var(phi) / Var(Ai)
# Uses the exact average marginal variance of the CAR draw, sigma2 * mean(diag(Q^-1))
compute_sigma2_phi <- function(R2_phi, var_X, sigma_A2, W, rho) {
  Q <- diag(rowSums(W)) - rho * W
  target_var_phi <- R2_phi * (var_X + sigma_A2) / (1 - R2_phi)
  target_var_phi / mean(diag(solve(Q)))
}

# Area-level confounders Xi1-Xi4: three Normal and one right-skewed Gamma (Xi2, shape 2),
# standardised to mean 0 and variance 1 so var_X = sum(gamma^2) still holds. The Gamma one sits on
# the smallest outcome coefficient so its skew adds little spread to the sample-vs-population truth gap.
draw_confounders <- function(n) {
  std_gamma <- function(n, shape) (rgamma(n, shape = shape, rate = 1) - shape) / sqrt(shape)
  list(Xi1 = rnorm(n, 0, 1), Xi2 = std_gamma(n, 2),
       Xi3 = rnorm(n, 0, 1), Xi4 = rnorm(n, 0, 1))
}

# Linear-predictor contribution of the confounders
cov_effects <- function(Xi1, Xi2, Xi3, Xi4) {
  (-0.5)*Xi1 + (-0.25)*Xi2 + (0.25)*Xi3 + (0.5)*Xi4
}

# Exposure-response term f(a). Only the nonlinear form uses Xi1, through a bounded exposure x Xi1 interaction
erf_terms <- function(a, Xi1, out_form) {
  if (out_form == "linear") {
    0.10 * (a - 10)
  } else if (out_form == "quadratic") {
    # Vertex at a = 13 (10 + 0.30/(2*0.05)), inside the [5,15] evaluation window
    0.30*(a - 10) - 0.05*(a - 10)^2
  } else {
    # Up, local max (a ~ 9), down, local min (a ~ 12.5), up: a linear trend plus one sine wave
    # The tanh interaction is bounded (<= 0.25 per SD of Xi1), unlike Josey's (a - 10)*x1
    0.12*(a - 10) + 0.55*(sin(2*pi*((a - 10) + 3.75)/9) - 0.5) - 0.25*tanh((a - 10)/4)*Xi1
  }
}

# Area populations: log-normal, median about 4,900 (log-scale mean 8.5, SD 1.2)
draw_populations <- function(n) pmax(round(exp(rnorm(n, mean = 8.5, sd = 1.2))), 1)

# ============================================================
# Data generating mechanism
# ============================================================

# DAG: Xi, phi -> Ai;  Ai -> Zi <- Di;  Xi, Ai -> Yi  (Ni is the outcome offset)
# Exposure scale: cumulative wildfire-smoke PM2.5 (ug/m3), mean 10, bulk range about 5-15
#
# seed:     random seed for reproducibility; NULL skips set.seed() and continues the current random stream
#           (used by the replicate loop in stan/sim_latent_replicates.R, which seeds once)
# n_areas:  number of FSA postal code areas I (must be a perfect square; 400 = 20 x 20 grid)
# gamma:    exposure model coefficients, gamma0 = 10 is the exposure mean; the signs of gamma1-4 against the outcome
#           coefficients in cov_effects set the confounding (sum of products 0.375, an unadjusted log-slope bias of about 0.07)
# sigma_A:  residual SD of true exposure; with gamma = +-0.75 the total SD is about 2.4, so A and Z are never negative in practice (no clipping)
# rho_A:    spatial correlation in exposure CAR prior, 0 = spatially independent
# R2_phi:   target proportion of Var(Ai) explained by spatial random effect phi
# eta0:     baseline log error variance, exp(eta0) is error at Di = 0
# eta1:     Di-slope, eta1 = 0 recovers Josey homoscedastic Assumption 2
# beta1:    attenuation slope in Di: E[Zi | Ai, Di] = exp(beta1*Di) * Ai; beta1 = 0 means no attenuation
# out_form: ERF functional form: linear / quadratic / nonlinear
# xi:       NegBin dispersion, xi = Inf gives Poisson (base model)

data.generate <- function(seed = 1006092577,
                          n_areas = 400,
                          gamma = c(gamma0 = 10, gamma1 = -0.75, gamma2 = 0.75,
                                    gamma3 = -0.75, gamma4 = 0.75),
                          sigma_A = 1.5,
                          rho_A = 0.9,
                          R2_phi = 0.20,
                          eta0 = -0.5,
                          eta1 = 0.5,
                          beta1 = -0.1,
                          out_form = c("linear", "quadratic", "nonlinear"),
                          xi = Inf) {
  if (!is.null(seed)) set.seed(seed)
  out_form <- match.arg(out_form)
  I <- n_areas

  # Step 1: Area-level covariates
  # Xi1-Xi4: independent confounders affecting both exposure and outcome (see draw_confounders)
  # Di: standardised Gamma(2,1) distance proxy, mean 0 var 1 after standardisation
  # Ni: log-normal population sizes anchored to Ontario FSA Census distribution
  Xc <- draw_confounders(I)
  Xi1 <- Xc$Xi1; Xi2 <- Xc$Xi2; Xi3 <- Xc$Xi3; Xi4 <- Xc$Xi4
  Di_raw <- rgamma(I, shape = 2, rate = 1)
  Di <- (Di_raw - 2) / sqrt(2)
  Ni <- draw_populations(I)

  # Step 2: Adjacency matrix and average neighbourhood size
  W <- make_adj_matrix(I)
  d_bar <- mean(rowSums(W))

  # Step 3: Spatial random effect on exposure
  # var_X = sum(gamma^2) since each Xi has variance 1
  # sigma2_phi chosen so phi explains R2_phi of Var(Ai)
  var_X <- sum(gamma[c("gamma1","gamma2","gamma3","gamma4")]^2)
  sigma2_phi <- compute_sigma2_phi(R2_phi, var_X, sigma_A^2, W, rho_A)
  phi <- if (rho_A > 0) sim_car(W, rho_A, sigma2_phi) else rep(0, I)

  # Step 4: Latent true exposure
  # Ai | Xi, phi ~ N(10 - 0.75*Xi1 + 0.75*Xi2 - 0.75*Xi3 + 0.75*Xi4 + phi, 1.5^2); no clipping
  # Signs give real confounding: sum(gamma_k * beta_k) = 0.375 (the outcome coefficients are in cov_effects)
  mu_A <- gamma["gamma0"] + gamma["gamma1"]*Xi1 + gamma["gamma2"]*Xi2 +
    gamma["gamma3"]*Xi3 + gamma["gamma4"]*Xi4
  Ai <- rnorm(I, mean = mu_A + phi, sd = sigma_A)

  # Step 5: Error-prone surrogate -- Version B (primary model)
  # Zi | Ai, Di ~ N(exp(beta1*Di) * Ai, exp(eta0 + eta1*Di))
  # Di standardised so eta1 = 1 means 1 SD increase in distance multiplies variance by e
  # exp(beta1*Di) is the calibration slope: about 1 on average (beta1 = -0.1) but attenuated far from monitors
  atten_i <- exp(beta1 * Di)
  sigma_Zi <- sqrt(exp(eta0 + eta1 * Di))
  Zi <- rnorm(I, mean = atten_i * Ai, sd = sigma_Zi)

  # Step 6: Poisson outcome with log-linear predictor and population offset
  # Yi ~ Poisson(Ni * exp(eta_i)), intercept -3 gives ~5% baseline event rate
  cov_terms <- cov_effects(Xi1, Xi2, Xi3, Xi4)
  exp_terms <- erf_terms(Ai, Xi1, out_form)
  log_mu_i <- log(Ni) + (-3) + exp_terms + cov_terms
  mu_i <- exp(log_mu_i)
  Yi <- if (is.infinite(xi)) {
    rpois(I, lambda = mu_i)
  } else {
    rnbinom(I, size = xi, prob = xi / (xi + mu_i))
  }

  data.frame(
    area_id = 1:I,
    Yi = Yi,
    Ni = Ni,
    Ai = Ai,
    Zi = Zi,
    Xi1 = Xi1, Xi2 = Xi2, Xi3 = Xi3, Xi4 = Xi4,
    Di = Di,
    Di_raw = Di_raw,
    atten_i = atten_i,
    phi = phi,
    sigma2_phi = sigma2_phi,
    out_form = out_form
  )
}

# ============================================================
# True ERF
# ============================================================

# Two truths for the average potential outcome theta(a) = E[Yi(a)/Ni], the event rate per person if every
# area's exposure were set to a:
#     theta(a) = sum_i Ni * r_i(a) / sum_i Ni,   r_i(a) = exp(-3 + f(a, Xi1) + cov(Xi))
# Neither uses the realised Ai, Zi, phi, Di or the Poisson noise: the intervention sets A = a, and phi and Di
# do not enter the outcome. Both truths use the same calculation (avg_rate) and differ only in the areas
# it is applied to:
#   true_erf:     this dataset's own areas (sample / finite-population truth, as in Josey et al. 2023)
#   true_erf_pop: a very large simulated population of areas (the proposal's population estimand)

avg_rate <- function(a, Xi1, Xi2, Xi3, Xi4, Ni, out_form) {
  weighted.mean(exp(-3 + erf_terms(a, Xi1, out_form) + cov_effects(Xi1, Xi2, Xi3, Xi4)), Ni)
}

true_erf <- function(a.vals, data) {
  out_form <- data$out_form[1]
  sapply(a.vals, function(a) avg_rate(a, data$Xi1, data$Xi2, data$Xi3, data$Xi4, data$Ni, out_form))
}

# Population version: draws n_mc areas (confounders and populations) and applies the same calculation
# but from many samples to reduce Monte Carlo error 
# so the truth is identical across runs.
true_erf_pop <- function(a.vals, out_form, n_mc = 2e6, seed = 5155) {
  set.seed(seed)
  Xc <- draw_confounders(n_mc)
  Ni <- draw_populations(n_mc)
  sapply(a.vals, function(a) avg_rate(a, Xc$Xi1, Xc$Xi2, Xc$Xi3, Xc$Xi4, Ni, out_form))
}
