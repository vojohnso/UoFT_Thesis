library(MASS)
library(spdep)
library(BART)
library(gtools)
library(ggplot2)

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
    exp_terms <- 0.15*(Ai - 7) - 0.008*(Ai - 7)^2
  } else if (out_form == "nonlinear") {
    exp_terms <- 0.10*(Ai - 7) - 0.75*cos(pi*(Ai - 5)/8) - 0.25*(Ai - 7)*Xi1
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

# Reference population: 100k covariate draws fixed once so theta(a) is stable
# across all replicates and does not depend on the simulated dataset
make_reference_pop <- function(n_ref = 100000) {
  set.seed(1006092577)
  data.frame(
    Xi1 = rnorm(n_ref), Xi2 = rnorm(n_ref),
    Xi3 = rnorm(n_ref), Xi4 = rnorm(n_ref)
  )
}
REF_POP <- make_reference_pop()

# Compute true APO theta(a) = E[Yi(a)] following Josey et al. (2023)
# Sets A = a for all reference units, averages Ni * exp(eta_i(a)) over the reference population
true_erf <- function(a.vals, out_form, ref_pop = REF_POP) {
  cov_terms <- (-0.5)*ref_pop$Xi1 + (-0.25)*ref_pop$Xi2 +
    (0.25)*ref_pop$Xi3 + (0.5)*ref_pop$Xi4
  set.seed(42)
  Ni_ref <- pmax(round(exp(rnorm(nrow(ref_pop), mean = 8.5, sd = 1.2))), 1)
  out <- numeric(length(a.vals))
  for (i in seq_along(a.vals)) {
    a <- a.vals[i]
    if (out_form == "linear") {
      exp_terms <- 0.10 * (a - 7)
    } else if (out_form == "quadratic") {
      exp_terms <- 0.15*(a - 7) - 0.008*(a - 7)^2
    } else if (out_form == "nonlinear") {
      exp_terms <- 0.10*(a - 7) - 0.75*cos(pi*(a - 5)/8) - 0.25*(a - 7)*ref_pop$Xi1
    }
    out[i] <- mean(Ni_ref * exp(-3 + exp_terms + cov_terms))
  }
  out
}

# ============================================================
# BART g-computation
# ============================================================

# Fits nonparametric Bayesian g-computation using BART as the outcome model.
# Outcome modelled on log(Yi/Ni), predictions exponentiated and multiplied by Ni
# to recover count scale: theta(a) = E[Yi(a)].
# Bayesian bootstrap applies Dirichlet weights over areas per posterior draw.
#
# data:       output of data.generate()
# a_grid:     exposure values at which to evaluate the ERF
# ndpost:     number of posterior draws to keep
# nskip:      burn-in draws to discard
# use_true_A: TRUE = oracle (uses Ai), FALSE = naive (uses Zi)

bart_gcomp <- function(data,
                       a_grid = seq(2, 14, by = 1),
                       ndpost = 1000,
                       nskip = 500,
                       use_true_A = FALSE) {
  data$exposure <- if (use_true_A) data$Ai else data$Zi
  model_label <- if (use_true_A) "Oracle (true A)" else "Naive (surrogate Z)"
  
  # Centre exposure for numerical stability
  exp_mean <- mean(data$exposure)
  data$exp_c <- data$exposure - exp_mean
  a_grid_c <- a_grid - exp_mean
  
  cat("\nFitting", model_label, "...\n")
  
  # Outcome on log rate scale -- gbart uses wbart internally for continuous outcome
  Y_train <- log(data$Yi / data$Ni)
  Y_train[!is.finite(Y_train)] <- log(0.5 / data$Ni[!is.finite(Y_train)])
  X_train <- as.matrix(data[, c("exp_c", "Xi1", "Xi2", "Xi3", "Xi4")])
  n <- nrow(data)
  n_grid <- length(a_grid)
  
  # Counterfactual test matrix: n_grid copies of X_train with exp_c set to each grid value
  X_test <- do.call(rbind, lapply(a_grid_c, function(a) {
    X <- X_train
    X[, "exp_c"] <- a
    X
  }))
  
  # Fit BART -- yhat.test is ndpost x (n_grid * n) on log rate scale
  bart_fit <- gbart(X_train, Y_train, x.test = X_test,
                    type = "wbart", ndpost = ndpost, nskip = nskip, printevery = 0)
  
  # G-computation: Bayesian bootstrap over areas for each draw and grid point
  M <- ndpost
  theta_draws <- matrix(NA, nrow = M, ncol = n_grid)
  colnames(theta_draws) <- paste0("a=", a_grid)
  for (g in seq_len(n_grid)) {
    col_idx <- ((g - 1) * n + 1):(g * n)
    count_pred <- data$Ni * exp(bart_fit$yhat.test[, col_idx])
    for (m in seq_len(M)) {
      bb_w <- rdirichlet(1, rep(1, n))
      theta_draws[m, g] <- sum(bb_w * count_pred[m, ])
    }
  }
  
  theta_summary <- data.frame(
    a = a_grid,
    mean = apply(theta_draws, 2, mean),
    sd = apply(theta_draws, 2, sd),
    q025 = apply(theta_draws, 2, quantile, 0.025),
    q50 = apply(theta_draws, 2, quantile, 0.50),
    q975 = apply(theta_draws, 2, quantile, 0.975),
    row.names = NULL
  )
  
  cat("\n", model_label, ": theta(a)\n")
  print(theta_summary, digits = 4)
  
  list(theta_draws = theta_draws, theta_summary = theta_summary,
       bart_fit = bart_fit, model_label = model_label, a_grid = a_grid)
}

# ============================================================
# Plotting and comparison
# ============================================================

# Three-curve ERF plot: oracle, naive, and true ERF with 95% credible bands
make_erf_plot <- function(oracle, naive, true_theta, a_grid, title) {
  plot_df <- rbind(
    data.frame(a = a_grid, mean = oracle$theta_summary$mean,
               q025 = oracle$theta_summary$q025, q975 = oracle$theta_summary$q975,
               method = "Oracle (true A)"),
    data.frame(a = a_grid, mean = naive$theta_summary$mean,
               q025 = naive$theta_summary$q025, q975 = naive$theta_summary$q975,
               method = "Naive (surrogate Z)"),
    data.frame(a = a_grid, mean = true_theta, q025 = true_theta, q975 = true_theta,
               method = "True ERF")
  )
  ggplot(plot_df, aes(x = a, y = mean, colour = method, fill = method)) +
    geom_ribbon(data = subset(plot_df, method != "True ERF"),
                aes(ymin = q025, ymax = q975), alpha = 0.15, colour = NA) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    scale_colour_manual(values = c("True ERF" = "black",
                                   "Oracle (true A)" = "#2166ac",
                                   "Naive (surrogate Z)" = "#d6604d")) +
    scale_fill_manual(values = c("True ERF" = "black",
                                 "Oracle (true A)" = "#2166ac",
                                 "Naive (surrogate Z)" = "#d6604d")) +
    labs(x = expression(paste("PM"[2.5], " (", mu, "g/m"^3, ")")),
         y = expression(theta(a) == E*"["*Y[i](a)*"]"),
         colour = NULL, fill = NULL, title = title,
         subtitle = "Ribbons show 95% posterior credible intervals") +
    theme_bw(base_size = 13) +
    theme(legend.position = "bottom", panel.grid.minor = element_blank())
}

# Comparison table: true theta, oracle, naive, bias, percent bias at each grid point
make_comparison <- function(oracle, naive, a_grid, true_theta) {
  data.frame(
    a = a_grid,
    true = round(true_theta, 4),
    oracle = round(oracle$theta_summary$mean, 4),
    naive = round(naive$theta_summary$mean, 4),
    bias = round(naive$theta_summary$mean - true_theta, 4),
    pct_bias = round(100 * (naive$theta_summary$mean - true_theta) / true_theta, 2)
  )
}

# ============================================================
# Sanity checks
# ============================================================

cat("Baseline check: no spatial, linear ERF\n")
sim1 <- data.generate(seed = 401, n_areas = 196, out_form = "linear",
                      rho_A = 0, R2_phi = 0.30, eta1 = 1.0)
cat("Ai:   mean =", round(mean(sim1$Ai), 2), "sd =", round(sd(sim1$Ai), 2),
    "min =", round(min(sim1$Ai), 2), "\n")
cat("Zi:   mean =", round(mean(sim1$Zi), 2), "sd =", round(sd(sim1$Zi), 2), "\n")
cat("Di:   mean =", round(mean(sim1$Di), 2), "sd =", round(sd(sim1$Di), 2),
    "(should be ~0, ~1)\n")
cat("Ni:   median =", as.integer(median(sim1$Ni)),
    "5th =", as.integer(quantile(sim1$Ni, 0.05)),
    "95th =", as.integer(quantile(sim1$Ni, 0.95)), "\n")
cat("Rate: mean =", round(mean(sim1$Yi / sim1$Ni), 4), "\n")

cat("\nSpatial check: rho_A = 0.9, R2_phi = 0.30\n")
sim2 <- data.generate(seed = 401, n_areas = 196, out_form = "linear",
                      rho_A = 0.9, R2_phi = 0.30)
cat("sigma2_phi used:", round(sim2$sigma2_phi[1], 2), "\n")
cat("phi sd:", round(sd(sim2$phi), 3), "\n")
cat("R2_phi actual:", round(var(sim2$phi) / var(sim2$Ai), 3), "\n")

cat("\nTrue ERF theta(a) at selected grid points\n")
a_grid_check <- seq(2, 14, by = 2)
print(data.frame(a = a_grid_check, theta = round(true_erf(a_grid_check, "linear"), 3)))

# ============================================================
# Run all three ERF scenarios
# ============================================================

a_grid <- seq(2, 14, by = 1)

for (form in c("linear", "quadratic", "nonlinear")) {
  cat("\n", toupper(form), "\n")
  cat(strrep("-", 40), "\n")
  sim <- data.generate(seed = 401, n_areas = 196, out_form = form)
  cat("Ai:   mean =", round(mean(sim$Ai), 1), "sd =", round(sd(sim$Ai), 1),
      "range = [", round(min(sim$Ai), 1), ",", round(max(sim$Ai), 1), "]\n")
  cat("Zi:   mean =", round(mean(sim$Zi), 1), "sd =", round(sd(sim$Zi), 1),
      "cor(A,Z) =", round(cor(sim$Ai, sim$Zi), 3), "\n")
  cat("Rate: mean =", round(mean(sim$Yi / sim$Ni), 4), "\n")
  true_theta <- true_erf(a_grid, form)
  naive <- bart_gcomp(sim, a_grid = a_grid, use_true_A = FALSE)
  oracle <- bart_gcomp(sim, a_grid = a_grid, use_true_A = TRUE)
  cat("\nComparison table\n")
  print(make_comparison(oracle, naive, a_grid, true_theta), digits = 4)
  print(make_erf_plot(oracle, naive, true_theta, a_grid,
                      title = paste("ERF:", form, "outcome -- Ontario PM2.5")))
}