library(ggplot2)

source("gen_data.R")
source("estimators.R")

# ============================================================
# Output location
# ============================================================

output_dir <- "output"
dir.create(output_dir, showWarnings = FALSE)

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

cat("\nTrue ERF theta(a) = E[Yi(a)/Ni] at selected grid points\n")
a_grid_check <- seq(2, 14, by = 2)
print(data.frame(a = a_grid_check,
                 theta = round(true_erf(a_grid_check, sim1), 5)))

# ============================================================
# Run all three ERF scenarios: naive (Zi) vs. oracle (Ai) vs. true theta(a)
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
  true_theta <- true_erf(a_grid, sim)
  naive <- bart_gcomp(sim, a_grid = a_grid, use_true_A = FALSE)
  oracle <- bart_gcomp(sim, a_grid = a_grid, use_true_A = TRUE)

  comparison <- make_comparison(oracle, naive, a_grid, true_theta)
  cat("\nComparison table\n")
  print(comparison, digits = 5)

  erf_plot <- make_erf_plot(oracle, naive, true_theta, a_grid,
                            title = paste("ERF:", form, "outcome -- Ontario PM2.5"))
  print(erf_plot)

  # Save outputs: plot (PNG), comparison table (CSV), trimmed result objects (RDS)
  # Drop bart_fit (full posterior tree object, ~45MB each) -- theta_draws/theta_summary
  # already contain everything needed for downstream analysis or re-plotting
  naive_slim <- naive[c("theta_draws", "theta_summary", "model_label", "a_grid")]
  oracle_slim <- oracle[c("theta_draws", "theta_summary", "model_label", "a_grid")]

  ggsave(file.path(output_dir, sprintf("erf_plot_%s.png", form)),
         erf_plot, width = 8, height = 6, dpi = 150)
  write.csv(comparison, file.path(output_dir, sprintf("comparison_%s.csv", form)),
            row.names = FALSE)
  saveRDS(list(sim = sim, true_theta = true_theta, naive = naive_slim,
              oracle = oracle_slim, comparison = comparison, a_grid = a_grid),
          file.path(output_dir, sprintf("results_%s.rds", form)))
  cat("Saved plot, comparison table, and results to", output_dir, "\n")
}

cat("\nAll outputs saved to:", normalizePath(output_dir), "\n")
