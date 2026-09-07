library(MASS)
library(gtools)
library(BART)
library(ggplot2)

# DAG: Xi, Di -> Ai -> Zi; Xi, Ai -> Yi
# Exposure scale: short-term wildfire smoke PM2.5, range 0-100 ug/m3

make_adj_matrix <- function(I) {
  g <- sqrt(I)
  if (g != floor(g)) stop("I must be a perfect square")
  W <- matrix(0, I, I)
  for (i in 1:I) {
    row_i <- ceiling(i / g)
    col_i <- i - (row_i - 1) * g
    if (col_i < g) W[i, i + 1] <- 1
    if (col_i > 1) W[i, i - 1] <- 1
    if (row_i < g) W[i, i + g] <- 1
    if (row_i > 1) W[i, i - g] <- 1
  }
  W
}

sim_car <- function(W, rho, sigma2) {
  D <- diag(rowSums(W))
  Sigma <- sigma2 * solve(D - rho * W)
  as.numeric(mvrnorm(1, mu = rep(0, nrow(W)), Sigma = Sigma))
}


# data.generate: 
# (1) Generates four independent standard normal confounders Xi1-Xi4 which affect true exposure and outcome
# Di ~ Uniform(0, 1) as the calibration between true and observed exposure (Zi) and Ni is population offset ~ Uniform(10, 1000)
# (2) Builds an adjacency matrix and if rho_A > 0, draws spatial correlation random effects from CAR (phi_i and nu_i)
# (3) Generate true latent Ai | X_i, phi_i ~ N(gamma0 + gamma1*Xi1 + ... + gamma4*Xi4 + phi_i, sigma_A^2)
# (4) Generate error prone surrogate Z_i | A_i, D_i ~ N(beta_Z1*A_i + beta_Z2*D_i, sigma_Z^2)
# (5) Generate outcome count data using spatial correlation Y_i ~ Poisson(N_i * exp(eta_i) where eta_i has
# different functional forms dependent on the confounders Xi and nu_i

# data.generate arguments:
# seed: random seed
# n_areas: number of postal code areas I
# gamma: exposure model coefficients
#   Ai ~ Normal(gamma0 + gamma1*Xi1 + ... + gamma4*Xi4 + phi_i, sigma_A^2)
#   gamma0 = 30 centres exposure at 30 ug/m3 (realistic wildfire smoke mean)
# beta_Z1: scaling of Ai in surrogate, 1 = unbiased, <1 = attenuation
# beta_Z2: directional bias from Di on surrogate mean, 0 = no directional bias
# sigma_Z: measurement error SD
# sigma_A: SD of true exposure, 15 gives realistic spread across 0-100 ug/m3
# out_form: functional form of exposure-response, linear/quadratic/nonlinear
# rho_A: spatial correlation in exposure CAR prior, 0 = independent
# rho_Y: spatial correlation in outcome CAR prior, 0 = independent
# sigma2_phi: marginal variance of exposure spatial random effect
# sigma2_nu: marginal variance of outcome spatial random effect

data.generate <- function(seed = 401,
                          n_areas = 100,
                          gamma = c(gamma0 = 30, gamma1 = 3, gamma2 = -3, gamma3 = -3, gamma4 = 3),
                          beta_Z1 = 1,
                          beta_Z2 = 0,
                          sigma_Z = 5,
                          sigma_A = 15,
                          out_form = c("linear", "quadratic", "nonlinear"),
                          rho_A = 0,
                          rho_Y = 0,
                          sigma2_phi = 1,
                          sigma2_nu = 1) {
  
  options(digits = 4)
  set.seed(seed)
  out_form <- match.arg(out_form)
  I <- n_areas
  
  # Step 1: Area-level covariates
  # Xi1-Xi4 standard normal confounders affecting both exposure and outcome
  # Di calibration covariate explaining measurement error direction/magnitude
  # Ni population size offset, Uniform(10, 1000)
  Xi1 <- rnorm(I, 0, 1)
  Xi2 <- rnorm(I, 0, 1)
  Xi3 <- rnorm(I, 0, 1)
  Xi4 <- rnorm(I, 0, 1)
  Di <- runif(I, 0, 1)
  Ni <- round(runif(I, 10, 1000))
  
  # Step 2: Adjacency matrix -- sqrt(I) x sqrt(I) grid, rook contiguity
  W <- make_adj_matrix(I)
  
  # Step 3: Spatial random effect on exposure
  # phi_i captures unmeasured spatial drivers of true exposure
  # phi = 0 when rho_A = 0
  if (rho_A > 0) {
    phi <- sim_car(W, rho_A, sigma2_phi)
  } else {
    phi <- rep(0, I)
  }
  
  # Step 4: Latent true exposure
  # Ai | Xi, phi ~ Normal(mu_A + phi, sigma_A^2)
  # gamma0 = 30 centres distribution at 30 ug/m3 wildfire smoke PM2.5
  # sigma_A = 15 gives realistic spread, P(Ai < 0) is negligible
  mu_A <- gamma["gamma0"] + gamma["gamma1"]*Xi1 + gamma["gamma2"]*Xi2 +
    gamma["gamma3"]*Xi3 + gamma["gamma4"]*Xi4
  Ai <- rnorm(I, mean = mu_A + phi, sd = sigma_A)
  Ai <- pmax(Ai, 0)
  
  # Step 5: Error-prone surrogate
  # Zi | Ai, Di ~ Normal(beta_Z1*Ai + beta_Z2*Di, sigma_Z^2)
  # sigma_Z = 5 gives realistic CanOSSEM-like prediction error
  Zi <- rnorm(I, mean = beta_Z1 * Ai + beta_Z2 * Di, sd = sigma_Z)
  Zi <- pmax(Zi, 0)
  
  # Step 6: Spatial random effect on outcome
  # nu_i captures unmeasured spatial clustering in health outcomes
  # nu = 0 when rho_Y = 0
  if (rho_Y > 0) {
    nu <- sim_car(W, rho_Y, sigma2_nu)
  } else {
    nu <- rep(0, I)
  }
  
  # Step 7: Area-level count outcome Yi ~ Poisson(Ni * exp(eta_i))
  # intercept -3 gives baseline event rate exp(-3) ~ 5% at mean exposure
  cov_terms <- (-0.5)*Xi1 + (-0.25)*Xi2 + (0.25)*Xi3 + (0.5)*Xi4
  
  if (out_form == "linear") {
    # Monotone increasing
    exp_terms <- 0.015 * (Ai - 30)
  } else if (out_form == "quadratic") {
    # Concave - rises then flattens around 80-90 ug/m3
    exp_terms <- 0.025*(Ai - 30) - 0.0002*(Ai - 30)^2
  } else if (out_form == "nonlinear") {
    # S-shaped - slow rise, steep increase around 30-50, flattens then rises again
    # Rescaled Josey et al. cosine form for 0-100 ug/m3 range
    exp_terms <- 0.015*(Ai - 20) - 0.75*cos(pi*(Ai - 10)/40)
  }
  
  log_mu_i <- log(Ni) + (-3) + exp_terms + cov_terms + nu
  Yi <- rpois(I, lambda = exp(log_mu_i))
  
  data.frame(
    area_id = 1:I,
    Yi = Yi,
    Ni = Ni,
    Ai = Ai,
    Zi = Zi,
    Xi1 = Xi1, Xi2 = Xi2, Xi3 = Xi3, Xi4 = Xi4,
    Di = Di,
    phi = phi,
    nu = nu,
    out_form = out_form
  )
}

# Computes the true average potential outcome µ(a) at each exposure level in a.vals.
# It applies the exact same outcome formula as data.generate but sets exposure = a for every area 
# while keeping all covariates at their observed values. It then takes a population-weighted average
# using Ni as weights. 

true_erf <- function(a.vals, data) {
  out_form <- data$out_form[1]
  out <- rep(NA, length(a.vals))
  for (i in seq_along(a.vals)) {
    a.vec <- rep(a.vals[i], nrow(data))
    cov_terms <- (-0.5)*data$Xi1 + (-0.25)*data$Xi2 +
      (0.25)*data$Xi3 + (0.5)*data$Xi4
    if (out_form == "linear") {
      exp_terms <- 0.015 * (a.vec - 30)
    } else if (out_form == "quadratic") {
      exp_terms <- 0.025*(a.vec - 30) - 0.0002*(a.vec - 30)^2
    } else if (out_form == "nonlinear") {
      exp_terms <- 0.015*(a.vec - 20) - 0.75*cos(pi*(a.vec - 10)/40)
    }
    mu_out <- exp(-3 + exp_terms + cov_terms + data$nu)
    out[i] <- weighted.mean(mu_out, data$Ni)
  }
  return(out)
}

# Fits a nonparametric Bayesian g-computation using BART as the outcome model. 
# Argument use_true_A is a flag that switches between true exposure Ai and naive Zi and exposure 
# is centered for stability. Outcome is modelled on log(Yi/Ni) since gbart on type = "wbart" assumes
# continuous Gaussian outcome and then predictions are exponentiated back to rate scale
# Bayesian bootstrap applies Dirichlet weights over areas for each posterior draw

bart_gcomp <- function(data, a_grid = seq(0, 100, by = 10),
                       ndpost = 1000,
                       nskip = 500,
                       use_true_A = FALSE) {
  
  data$exposure <- if (use_true_A) data$Ai else data$Zi
  model_label <- if (use_true_A) "Oracle BART (true A)" else "Naive BART (surrogate Z)"
  
  exp_mean <- mean(data$exposure)
  data$exposure_c <- data$exposure - exp_mean
  a_grid_c <- a_grid - exp_mean
  
  cat(sprintf("\nFitting %s model...\n", model_label))
  
  Y_train <- log(data$Yi / data$Ni)
  # Non-negative
  Y_train[!is.finite(Y_train)] <- log(0.5 / data$Ni[!is.finite(Y_train)])
  X_train <- as.matrix(data[, c("exposure_c", "Xi1", "Xi2", "Xi3", "Xi4")])
  
  n <- nrow(data)
  n_grid <- length(a_grid)
  
  X_test <- do.call(rbind, lapply(a_grid_c, function(a) {
    X <- X_train
    X[, "exposure_c"] <- a
    X
  }))
  
  bart_fit <- gbart(X_train, Y_train, x.test = X_test,
                    type = "wbart", ndpost = ndpost, nskip = nskip, printevery = 0)
  
  M <- ndpost
  mu_draws <- matrix(NA, nrow = M, ncol = n_grid)
  colnames(mu_draws) <- paste0("a=", a_grid)
  
  for (g in seq_len(n_grid)) {
    col_idx <- ((g - 1) * n + 1):(g * n)
    rate_preds <- exp(bart_fit$yhat.test[, col_idx])
    for (m in seq_len(M)) {
      bb_w <- rdirichlet(1, rep(1, n))
      mu_draws[m, g] <- sum(bb_w * rate_preds[m, ])
    }
  }
  
  delta_draws <- matrix(NA, nrow = M, ncol = n_grid - 1)
  colnames(delta_draws) <- paste0("a=", a_grid[-1], " vs a=", a_grid[-n_grid])
  for (g in seq_len(n_grid - 1)) {
    delta_draws[, g] <- mu_draws[, g + 1] - mu_draws[, g]
  }
  
  summarise_draws <- function(mat) {
    data.frame(
      label = colnames(mat),
      mean = apply(mat, 2, mean),
      sd = apply(mat, 2, sd),
      q025 = apply(mat, 2, quantile, 0.025),
      q50 = apply(mat, 2, quantile, 0.50),
      q975 = apply(mat, 2, quantile, 0.975),
      row.names = NULL
    )
  }
  
  mu_summary <- summarise_draws(mu_draws)
  delta_summary <- summarise_draws(delta_draws)
  
  cat(sprintf("\n%s: mu(a)\n", model_label))
  print(mu_summary, digits = 4)
  cat(sprintf("\n%s: delta(a, a')\n", model_label))
  print(delta_summary, digits = 4)
  
  list(mu_draws = mu_draws, delta_draws = delta_draws,
       mu_summary = mu_summary, delta_summary = delta_summary,
       bart_fit = bart_fit, model_label = model_label, a_grid = a_grid)
}


# Builds a three-curve plot showcasing the ERF for oracle, naive and true 
make_erf_plot <- function(oracle, naive, true_mu, a_grid, title) {
  plot_df <- rbind(
    data.frame(a = a_grid, mean = oracle$mu_summary$mean,
               q025 = oracle$mu_summary$q025, q975 = oracle$mu_summary$q975,
               method = "Oracle (true A)"),
    data.frame(a = a_grid, mean = naive$mu_summary$mean,
               q025 = naive$mu_summary$q025, q975 = naive$mu_summary$q975,
               method = "Naive (surrogate Z)"),
    data.frame(a = a_grid, mean = true_mu, q025 = true_mu, q975 = true_mu,
               method = "True ERF")
  )
  ggplot(plot_df, aes(x = a, y = mean, colour = method, fill = method)) +
    geom_ribbon(data = subset(plot_df, method != "True ERF"),
                aes(ymin = q025, ymax = q975), alpha = 0.15, colour = NA) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    scale_colour_manual(values = c(
      "True ERF" = "black",
      "Oracle (true A)" = "#2166ac",
      "Naive (surrogate Z)" = "#d6604d")) +
    scale_fill_manual(values = c(
      "True ERF" = "black",
      "Oracle (true A)" = "#2166ac",
      "Naive (surrogate Z)" = "#d6604d")) +
    scale_x_continuous(breaks = a_grid) +
    labs(
      x = expression(paste("Wildfire Smoke PM"[2.5], " (", mu, "g/m"^3, ")")),
      y = "Average Event Rate",
      colour = NULL, fill = NULL,
      title = title,
      subtitle = "Ribbons show 95% posterior credible intervals") +
    theme_bw(base_size = 13) +
    theme(legend.position = "bottom", panel.grid.minor = element_blank())
}

# Builds a table with true_mu, oracle_mean, naive_mean, bias and percent bias for each grid point
make_comparison <- function(oracle, naive, a_grid, true_mu) {
  data.frame(
    a = a_grid,
    true_mu = round(true_mu, 4),
    oracle_mean = oracle$mu_summary$mean,
    naive_mean = naive$mu_summary$mean,
    bias = naive$mu_summary$mean - oracle$mu_summary$mean,
    pct_bias = 100 * (naive$mu_summary$mean - oracle$mu_summary$mean) /
      oracle$mu_summary$mean
  )
}

# Run all three scenarios
# Grid search from 0 to 100 for "normal" PM2.5 values 
a_grid <- seq(0, 100, by = 10)

# Loops between the three different functional forms of outcome
# Generates a dataset and computed true ERF, naive and oracle BART g-comp, prints the table and plots
for (form in c("linear", "quadratic", "nonlinear")) {
  cat(sprintf("\n========== %s ==========\n", toupper(form)))
  sim <- data.generate(seed = 401, n_areas = 100, out_form = form)
  cat(sprintf("Exposure -- mean: %.1f  sd: %.1f  range: [%.1f, %.1f]\n",
              mean(sim$Ai), sd(sim$Ai), min(sim$Ai), max(sim$Ai)))
  cat(sprintf("Surrogate -- mean: %.1f  sd: %.1f  cor(A,Z): %.3f\n",
              mean(sim$Zi), sd(sim$Zi), cor(sim$Ai, sim$Zi)))
  cat(sprintf("Observed event rate: %.4f\n", mean(sim$Yi / sim$Ni)))
  true_mu <- true_erf(a_grid, sim)
  naive <- bart_gcomp(sim, a_grid = a_grid, use_true_A = FALSE)
  oracle <- bart_gcomp(sim, a_grid = a_grid, use_true_A = TRUE)
  cat("\nComparison:\n")
  print(make_comparison(oracle, naive, a_grid, true_mu), digits = 4)
  p <- make_erf_plot(oracle, naive, true_mu, a_grid,
                     title = sprintf("ERF: %s outcome, wildfire smoke PM2.5", form))
  print(p)
}