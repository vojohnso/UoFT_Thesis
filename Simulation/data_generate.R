library(BART)
library(gtools)
library(ggplot2)

source("data_generate.R")

# ============================================================
# BART g-computation
# ============================================================

# Fits nonparametric Bayesian g-computation using BART as the outcome model.
# Outcome modelled on log(Yi/Ni) -- wbart assumes continuous Gaussian outcome.
# Predictions are exponentiated back to the rate scale then multiplied by Ni
# to recover the count scale: theta(a) = E[Yi(a)].
# Bayesian bootstrap applies Dirichlet weights over areas for each posterior draw
# to propagate covariate distribution uncertainty into the posterior of theta(a).

# Arguments:
#   data       -- output of data.generate()
#   a_grid     -- exposure values at which to evaluate the ERF
#   ndpost     -- number of posterior draws to keep
#   nskip      -- burn-in draws to discard
#   use_true_A -- TRUE = oracle (uses Ai), FALSE = naive (uses Zi)

bart_gcomp <- function(data,
                       a_grid     = seq(2, 14, by = 1),
                       ndpost     = 1000,
                       nskip      = 500,
                       use_true_A = FALSE) {
  
  data$exposure  <- if (use_true_A) data$Ai else data$Zi
  model_label    <- if (use_true_A) "Oracle (true A)" else "Naive (surrogate Z)"
  
  # Centre exposure for numerical stability
  exp_mean       <- mean(data$exposure)
  data$exp_c     <- data$exposure - exp_mean
  a_grid_c       <- a_grid - exp_mean
  
  cat("\nFitting", model_label, "...\n")
  
  # Outcome on log rate scale
  Y_train <- log(data$Yi / data$Ni)
  Y_train[!is.finite(Y_train)] <- log(0.5 / data$Ni[!is.finite(Y_train)])
  X_train <- as.matrix(data[, c("exp_c", "Xi1", "Xi2", "Xi3", "Xi4")])
  
  n      <- nrow(data)
  n_grid <- length(a_grid)
  
  # Build counterfactual test matrix: n_grid copies of X_train, each with exp_c = a_grid_c[g]
  X_test <- do.call(rbind, lapply(a_grid_c, function(a) {
    X           <- X_train
    X[, "exp_c"] <- a
    X
  }))
  
  # Fit BART -- yhat.test is ndpost x (n_grid * n) on log rate scale
  bart_fit <- gbart(X_train, Y_train, x.test = X_test,
                    type      = "wbart",
                    ndpost    = ndpost,
                    nskip     = nskip,
                    printevery = 0)
  
  # G-computation: for each draw m and grid point g, Bayesian bootstrap over areas
  M          <- ndpost
  theta_draws <- matrix(NA, nrow = M, ncol = n_grid)
  colnames(theta_draws) <- paste0("a=", a_grid)
  
  for (g in seq_len(n_grid)) {
    col_idx    <- ((g - 1) * n + 1):(g * n)
    count_pred <- data$Ni * exp(bart_fit$yhat.test[, col_idx])  # back to count scale
    for (m in seq_len(M)) {
      bb_w              <- rdirichlet(1, rep(1, n))
      theta_draws[m, g] <- sum(bb_w * count_pred[m, ])
    }
  }
  
  # Summarise posterior
  summarise_draws <- function(mat) {
    data.frame(
      a     = as.numeric(sub("a=", "", colnames(mat))),
      mean  = apply(mat, 2, mean),
      sd    = apply(mat, 2, sd),
      q025  = apply(mat, 2, quantile, 0.025),
      q50   = apply(mat, 2, quantile, 0.50),
      q975  = apply(mat, 2, quantile, 0.975),
      row.names = NULL
    )
  }
  
  theta_summary <- summarise_draws(theta_draws)
  
  cat("\n", model_label, ": theta(a)\n")
  print(theta_summary, digits = 4)
  
  list(
    theta_draws   = theta_draws,
    theta_summary = theta_summary,
    bart_fit      = bart_fit,
    model_label   = model_label,
    a_grid        = a_grid
  )
}

# ============================================================
# Plotting
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
    data.frame(a = a_grid, mean = true_theta,
               q025 = true_theta, q975 = true_theta,
               method = "True ERF")
  )
  ggplot(plot_df, aes(x = a, y = mean, colour = method, fill = method)) +
    geom_ribbon(data = subset(plot_df, method != "True ERF"),
                aes(ymin = q025, ymax = q975), alpha = 0.15, colour = NA) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    scale_colour_manual(values = c(
      "True ERF"          = "black",
      "Oracle (true A)"   = "#2166ac",
      "Naive (surrogate Z)" = "#d6604d")) +
    scale_fill_manual(values = c(
      "True ERF"          = "black",
      "Oracle (true A)"   = "#2166ac",
      "Naive (surrogate Z)" = "#d6604d")) +
    labs(
      x        = expression(paste("PM"[2.5], " (", mu, "g/m"^3, ")")),
      y        = expression(theta(a) == E*"["*Y[i](a)*"]"),
      colour   = NULL,
      fill     = NULL,
      title    = title,
      subtitle = "Ribbons show 95% posterior credible intervals"
    ) +
    theme_bw(base_size = 13) +
    theme(legend.position = "bottom", panel.grid.minor = element_blank())
}

# Comparison table: true theta, oracle, naive, bias, percent bias
make_comparison <- function(oracle, naive, a_grid, true_theta) {
  data.frame(
    a         = a_grid,
    true      = round(true_theta, 4),
    oracle    = round(oracle$theta_summary$mean, 4),
    naive     = round(naive$theta_summary$mean, 4),
    bias      = round(naive$theta_summary$mean - true_theta, 4),
    pct_bias  = round(100 * (naive$theta_summary$mean - true_theta) / true_theta, 2)
  )
}

# ============================================================
# Run all three ERF scenarios
# ============================================================

a_grid <- seq(2, 14, by = 1)

for (form in c("linear", "quadratic", "nonlinear")) {
  
  cat("\n", toupper(form), "\n")
  cat(strrep("-", 40), "\n")
  
  sim <- data.generate(seed = 401, n_areas = 196, out_form = form)
  cat("Ai:   mean =", round(mean(sim$Ai), 1), " sd =", round(sd(sim$Ai), 1),
      " range = [", round(min(sim$Ai), 1), ",", round(max(sim$Ai), 1), "]\n")
  cat("Zi:   mean =", round(mean(sim$Zi), 1), " sd =", round(sd(sim$Zi), 1),
      " cor(A,Z) =", round(cor(sim$Ai, sim$Zi), 3), "\n")
  cat("Rate: mean =", round(mean(sim$Yi / sim$Ni), 4), "\n")
  
  true_theta <- true_erf(a_grid, form)
  naive      <- bart_gcomp(sim, a_grid = a_grid, use_true_A = FALSE)
  oracle     <- bart_gcomp(sim, a_grid = a_grid, use_true_A = TRUE)
  
  cat("\nComparison table\n")
  print(make_comparison(oracle, naive, a_grid, true_theta), digits = 4)
  
  p <- make_erf_plot(oracle, naive, true_theta, a_grid,
                     title = paste("ERF:", form, "outcome -- Ontario PM2.5"))
  print(p)
}