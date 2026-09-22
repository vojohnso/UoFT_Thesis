library(BART)
library(gtools)

# ============================================================
# BART g-computation
# ============================================================

# Fits nonparametric Bayesian g-computation using BART as the outcome model.
# Outcome modelled on log(Yi/Ni) -- the log event rate per person.
# Predictions are exponentiated to recover the rate scale: theta(a) = E[Yi(a)/Ni].
# Bayesian bootstrap applies Dirichlet weights over areas per posterior draw,
# with Ni as importance weights to match the weighted mean in true_erf.
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

  # Outcome on log rate scale: log(Yi/Ni) = eta_i without the offset
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

  # G-computation: for each draw m and grid point g, compute Ni-weighted mean rate
  # Bayesian bootstrap draws Dirichlet weights then computes weighted mean of rates
  M <- ndpost
  theta_draws <- matrix(NA, nrow = M, ncol = n_grid)
  colnames(theta_draws) <- paste0("a=", a_grid)
  for (g in seq_len(n_grid)) {
    col_idx <- ((g - 1) * n + 1):(g * n)
    rate_pred <- exp(bart_fit$yhat.test[, col_idx])  # rate scale, no Ni
    for (m in seq_len(M)) {
      bb_w <- as.numeric(rdirichlet(1, rep(1, n))) * data$Ni  # Ni-weighted Dirichlet
      bb_w <- bb_w / sum(bb_w)
      theta_draws[m, g] <- sum(bb_w * rate_pred[m, ])
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
         y = expression(theta(a) == E*"["*Y[i](a)/N[i]*"]"),
         colour = NULL, fill = NULL, title = title,
         subtitle = "Ribbons show 95% posterior credible intervals") +
    theme_bw(base_size = 13) +
    theme(legend.position = "bottom", panel.grid.minor = element_blank())
}

# Comparison table: true theta, oracle, naive, bias, percent bias at each grid point
make_comparison <- function(oracle, naive, a_grid, true_theta) {
  data.frame(
    a = a_grid,
    true = round(true_theta, 5),
    oracle = round(oracle$theta_summary$mean, 5),
    naive = round(naive$theta_summary$mean, 5),
    bias = round(naive$theta_summary$mean - true_theta, 5),
    pct_bias = round(100 * (naive$theta_summary$mean - true_theta) / true_theta, 2)
  )
}
