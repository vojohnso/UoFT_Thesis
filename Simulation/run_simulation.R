source("data_generate.R")

# Simulation study: evaluate naive BART g-computation under measurement error
# For each scenario we run n_rep replicates and compute:
# bias: mean(mu_hat(a) - true_mu(a)) across replicates
# rmse: sqrt(mean((mu_hat(a) - true_mu(a))^2))
# coverage: fraction of 95% CIs containing true_mu(a)
# width: mean width of 95% CIs


# Global parameters for the simulation, n_rep = number of simulation generations
# n_areas = number of spatial areas
n_rep <- 100
n_areas <- 100
a_grid <- seq(0, 100, by = 10)
forms <- c("linear", "quadratic", "nonlinear")

# Take a specification (outcome functional form for example) and returns metrics across replications
run_sim_scenario <- function(out_form, 
                             n_rep, 
                             n_areas, 
                             a_grid,
                             beta_Z1 = 1, 
                             beta_Z2 = 0, 
                             sigma_Z = 5, 
                             use_true_A = FALSE) {
  
  n_grid <- length(a_grid)
  
  # Storage: one row per replicate per grid point
  results <- data.frame(
    rep = rep(1:n_rep, each = n_grid),
    a = rep(a_grid, times = n_rep),
    true_mu = NA,
    mu_hat = NA,
    q025 = NA,
    q975 = NA
  )
  
  # Reuses r as seed for reproducibility purposes 
  # Generates a dataset, computes true erf and BART g-comp
  for (r in 1:n_rep) {
    cat(sprintf("\r  Replicate %d / %d", r, n_rep))
    sim <- data.generate(
      seed = r,
      n_areas = n_areas,
      out_form = out_form,
      beta_Z1 = beta_Z1,
      beta_Z2 = beta_Z2,
      sigma_Z = sigma_Z
    )
    
    true_mu <- true_erf(a_grid, sim)
    
    fit <- bart_gcomp(sim, a_grid = a_grid, use_true_A = use_true_A,
                      ndpost = 500, nskip = 250)
    
    idx <- results$rep == r
    results$true_mu[idx] <- true_mu
    results$mu_hat[idx] <- fit$mu_summary$mean
    results$q025[idx] <- fit$mu_summary$q025
    results$q975[idx] <- fit$mu_summary$q975
  }
  cat("\n")
  
  # Compute performance metrics per grid point
  metrics <- do.call(rbind, lapply(a_grid, function(a) {
    sub <- results[results$a == a, ]
    bias <- mean(sub$mu_hat - sub$true_mu)
    rmse <- sqrt(mean((sub$mu_hat - sub$true_mu)^2))
    coverage <- mean(sub$true_mu >= sub$q025 & sub$true_mu <= sub$q975)
    width <- mean(sub$q975 - sub$q025)
    data.frame(a = a, bias = bias, rmse = rmse, coverage = coverage, width = width)
  }))
  
  list(results = results, metrics = metrics)
}

# Run naive and oracle for each outcome form 
all_metrics <- list()

for (form in forms) {
  cat(sprintf("\n %s: NAIVE \n", toupper(form)))
  naive_sim <- run_sim_scenario(form, n_rep, n_areas, a_grid, use_true_A = FALSE)
  
  cat(sprintf("\n%s: ORACLE \n", toupper(form)))
  oracle_sim <- run_sim_scenario(form, n_rep, n_areas, a_grid, use_true_A = TRUE)
  
  all_metrics[[form]] <- list(naive = naive_sim$metrics, oracle = oracle_sim$metrics)
  
  cat(sprintf("\n %s: Naive performance \n", form))
  print(naive_sim$metrics, digits = 3)
  cat(sprintf("\n %s: Oracle performance \n", form))
  print(oracle_sim$metrics, digits = 3)
}

# Plot bias across exposure levels for each scenario
plot_metrics <- function(naive_metrics, oracle_metrics, form, save_dir = ".") {
  df <- rbind(
    cbind(naive_metrics, method = "Naive (surrogate Z)"),
    cbind(oracle_metrics, method = "Oracle (true A)")
  )
  
  p_bias <- ggplot(df, aes(x = a, y = bias, colour = method)) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    scale_colour_manual(values = c(
      "Naive (surrogate Z)" = "#d6604d",
      "Oracle (true A)" = "#2166ac")) +
    labs(
      x = expression(paste("PM"[2.5], " (", mu, "g/m"^3, ")")),
      y = "Bias",
      colour = NULL,
      title = sprintf("Bias across exposure levels: %s outcome", form),
      subtitle = "Dashed line at zero indicates no bias"
    ) +
    theme_bw(base_size = 13) +
    theme(legend.position = "bottom", panel.grid.minor = element_blank())
  
  p_coverage <- ggplot(df, aes(x = a, y = coverage, colour = method)) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    geom_hline(yintercept = 0.95, linetype = "dashed", colour = "grey50") +
    scale_colour_manual(values = c(
      "Naive (surrogate Z)" = "#d6604d",
      "Oracle (true A)" = "#2166ac")) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
    labs(
      x = expression(paste("PM"[2.5], " (", mu, "g/m"^3, ")")),
      y = "Coverage",
      colour = NULL,
      title = sprintf("95%% CI coverage: %s outcome", form),
      subtitle = "Dashed line at nominal 0.95"
    ) +
    theme_bw(base_size = 13) +
    theme(legend.position = "bottom", panel.grid.minor = element_blank())
  
  p_rmse <- ggplot(df, aes(x = a, y = rmse, colour = method)) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    scale_colour_manual(values = c(
      "Naive (surrogate Z)" = "#d6604d",
      "Oracle (true A)" = "#2166ac")) +
    labs(
      x = expression(paste("PM"[2.5], " (", mu, "g/m"^3, ")")),
      y = "RMSE",
      colour = NULL,
      title = sprintf("RMSE across exposure levels: %s outcome", form)
    ) +
    theme_bw(base_size = 13) +
    theme(legend.position = "bottom", panel.grid.minor = element_blank())
  
  p_width <- ggplot(df, aes(x = a, y = width, colour = method)) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    scale_colour_manual(values = c(
      "Naive (surrogate Z)" = "#d6604d",
      "Oracle (true A)" = "#2166ac")) +
    labs(
      x = expression(paste("PM"[2.5], " (", mu, "g/m"^3, ")")),
      y = "Mean CI width",
      colour = NULL,
      title = sprintf("95%% CI width: %s outcome", form),
      subtitle = "Wider intervals indicate more propagated uncertainty"
    ) +
    theme_bw(base_size = 13) +
    theme(legend.position = "bottom", panel.grid.minor = element_blank())
  
  ggsave(file.path(save_dir, sprintf("bias_%s.png", form)),
         p_bias, width = 7, height = 5, dpi = 150)
  ggsave(file.path(save_dir, sprintf("coverage_%s.png", form)),
         p_coverage, width = 7, height = 5, dpi = 150)
  ggsave(file.path(save_dir, sprintf("rmse_%s.png", form)),
         p_rmse, width = 7, height = 5, dpi = 150)
  ggsave(file.path(save_dir, sprintf("width_%s.png", form)),
         p_width, width = 7, height = 5, dpi = 150)
  
  cat(sprintf("Plots saved for %s scenario\n", form))
  list(bias = p_bias, coverage = p_coverage, rmse = p_rmse, width = p_width)
}

for (form in forms) {
  plots <- plot_metrics(all_metrics[[form]]$naive,
                        all_metrics[[form]]$oracle, form)
  print(plots$bias)
  print(plots$coverage)
}

# Save results
saveRDS(all_metrics, "sim_results.rds")
cat("\nResults saved to sim_results.rds\n")