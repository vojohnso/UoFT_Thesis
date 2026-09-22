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

# Solve for sigma2_phi that achieves a target R2_phi
# R2_phi = Var(phi) / Var(Ai) where Var(phi) ~= sigma2_phi / d_bar
compute_sigma2_phi <- function(R2_phi, var_X, sigma_A2, d_bar) {
  d_bar * R2_phi * (var_X + sigma_A2) / (1 - R2_phi)
}

# Piecewise-linear rescaling of (a - 7) onto Josey et al. (2023)'s own
# u = a - 10 in [-4, 4] (their window [6,14], mean 10). Our window [2,14] is
# asymmetric around our mean (7): 5 units below, 7 above. Scaling each side
# separately preserves Josey's exact curve shape/coefficients while keeping
# the mapped range comparable to theirs (u in [-4,4]) instead of overshooting
# on the longer right side.
make_z_nonlinear <- function(a) {
  raw <- a - 7
  scale <- ifelse(raw < 0, 4/5, 4/7)
  raw * scale
}

# ============================================================
# Data generating mechanism
# ============================================================

# DAG: Xi, Di -> Ai -> Zi;  Xi, Ai -> Yi
# Exposure scale: Ontario annual average PM2.5, bulk range 2-14 ug/m3
#
# seed:     random seed for reproducibility
# n_areas:  number of FSA postal code areas I (must be a perfect square)
# gamma:    exposure model coefficients, gamma0 = 7 anchors to Ontario PM2.5 mean
# sigma_A:  SD of true exposure, sigma_A = 3 gives realistic Ontario spatial variation
# rho_A:    spatial correlation in exposure CAR prior, 0 = spatially independent
# R2_phi:   target proportion of Var(Ai) explained by spatial random effect phi
# eta0:     baseline log error variance, exp(eta0) is error at Di = 0
# eta1:     Di-slope, eta1 = 0 recovers Josey homoscedastic Assumption 2
# out_form: ERF functional form: linear / quadratic / nonlinear
# xi:       NegBin dispersion, xi = Inf gives Poisson (base model)

data.generate <- function(seed = 1006092577,
                          n_areas = 196,
                          gamma = c(gamma0 = 7, gamma1 = 1.5, gamma2 = -1.5,
                                    gamma3 = -1.5, gamma4 = 1.5),
                          sigma_A = 3,
                          rho_A = 0,
                          R2_phi = 0.30,
                          eta0 = 0,
                          eta1 = 1.0,
                          out_form = c("linear", "quadratic", "nonlinear"),
                          xi = Inf) {
  set.seed(seed)
  out_form <- match.arg(out_form)
  I <- n_areas

  # Step 1: Area-level covariates
  # Xi1-Xi4: iid N(0,1) confounders affecting both exposure and outcome
  # Di: standardised Gamma(2,1) distance proxy, mean 0 var 1 after standardisation
  # Ni: log-normal population sizes anchored to Ontario FSA Census distribution
  Xi1 <- rnorm(I, 0, 1)
  Xi2 <- rnorm(I, 0, 1)
  Xi3 <- rnorm(I, 0, 1)
  Xi4 <- rnorm(I, 0, 1)
  Di_raw <- rgamma(I, shape = 2, rate = 1)
  Di <- (Di_raw - 2) / sqrt(2)
  Ni <- pmax(round(exp(rnorm(I, mean = 8.5, sd = 1.2))), 1)

  # Step 2: Adjacency matrix and average neighbourhood size
  W <- make_adj_matrix(I)
  d_bar <- mean(rowSums(W))

  # Step 3: Spatial random effect on exposure
  # var_X = sum(gamma^2) since Xi ~ N(0,1)
  # sigma2_phi chosen so phi explains R2_phi of Var(Ai)
  var_X <- sum(gamma[c("gamma1","gamma2","gamma3","gamma4")]^2)
  sigma2_phi <- compute_sigma2_phi(R2_phi, var_X, sigma_A^2, d_bar)
  phi <- if (rho_A > 0) sim_car(W, rho_A, sigma2_phi) else rep(0, I)

  # Step 4: Latent true exposure
  # Ai | Xi, phi ~ N(7 + 1.5*Xi1 - 1.5*Xi2 - 1.5*Xi3 + 1.5*Xi4 + phi, 9)
  mu_A <- gamma["gamma0"] + gamma["gamma1"]*Xi1 + gamma["gamma2"]*Xi2 +
    gamma["gamma3"]*Xi3 + gamma["gamma4"]*Xi4
  Ai <- pmax(rnorm(I, mean = mu_A + phi, sd = sigma_A), 0)

  # Step 5: Error-prone surrogate -- Version B (primary model)
  # Zi | Ai, Di ~ N(Ai, exp(eta0 + eta1*Di))
  # Di standardised so eta1 = 1 means 1 SD increase in distance multiplies variance by e
  sigma_Zi <- sqrt(exp(eta0 + eta1 * Di))
  Zi <- rnorm(I, mean = Ai, sd = sigma_Zi)

  # Step 6: Poisson outcome with log-linear predictor and population offset
  # Yi ~ Poisson(Ni * exp(eta_i)), intercept -3 gives ~5% baseline event rate
  cov_terms <- (-0.5)*Xi1 + (-0.25)*Xi2 + (0.25)*Xi3 + (0.5)*Xi4
  if (out_form == "linear") {
    exp_terms <- 0.10 * (Ai - 7)
  } else if (out_form == "quadratic") {
    # Vertex at Ai = 12 (7 + 0.35/(2*0.035)), inside the [2,14] evaluation window,
    # so theta(a) visibly rises then declines rather than staying on the rising
    # branch the whole time (the previous 0.15/-0.008 pair had its vertex at
    # Ai = 16.4, outside the window -- always rising, indistinguishable from linear)
    exp_terms <- 0.35*(Ai - 7) - 0.035*(Ai - 7)^2
  } else if (out_form == "nonlinear") {
    # Josey et al. (2023) nonlinear form (linear + cosine + exposure x confounder
    # interaction, coefficients 0.25/0.75/0.25), rescaled from their symmetric
    # window a-10 in [-4,4] onto our asymmetric window Ai-7 in [-5,7] so the
    # curve stays bounded (comparable dynamic range) instead of exploding on
    # the long right tail. See make_z_nonlinear() below.
    z <- make_z_nonlinear(Ai)
    exp_terms <- 0.25*(z + 2) - 0.75*cos(pi*(z + 4)/4) - 0.25*z*Xi1
  }
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
    phi = phi,
    sigma2_phi = sigma2_phi,
    out_form = out_form
  )
}

# ============================================================
# True ERF
# ============================================================

# Compute true APO theta(a) = E[Yi(a)/Ni] -- population-weighted average event rate
# following Josey et al. (2023) predict_example: weighted.mean(mu_out, Ni)
# Rate scale removes Ni dependence and is not affected by the population offset
# Uses THIS replicate's own simulated covariates (data$Xi1, ...), not a fixed external
# population, following Josey et al. (2023) predict_example -- avoids Jensen's inequality
# explosion from averaging exp(interaction term) over a freshly drawn covariate set

true_erf <- function(a.vals, data) {
  out_form <- data$out_form[1]
  cov_terms <- (-0.5)*data$Xi1 + (-0.25)*data$Xi2 +
    (0.25)*data$Xi3 + (0.5)*data$Xi4
  out <- numeric(length(a.vals))
  for (i in seq_along(a.vals)) {
    a <- a.vals[i]
    if (out_form == "linear") {
      exp_terms <- 0.10 * (a - 7)
    } else if (out_form == "quadratic") {
      exp_terms <- 0.35*(a - 7) - 0.035*(a - 7)^2
    } else if (out_form == "nonlinear") {
      z <- make_z_nonlinear(a)
      exp_terms <- 0.25*(z + 2) - 0.75*cos(pi*(z + 4)/4) - 0.25*z*data$Xi1
    }
    # Rate per person: exp(-3 + f(a) + X*beta), no Ni offset
    # Weighted mean by Ni following Josey predict_example
    mu_out <- exp(-3 + exp_terms + cov_terms)
    out[i] <- weighted.mean(mu_out, data$Ni)
  }
  out
}
