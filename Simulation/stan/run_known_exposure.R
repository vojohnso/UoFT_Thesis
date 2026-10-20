# Known exposure: Poisson outcome model with the true exposure A observed (no measurement error), g-computation inside Stan.
# Run from the Simulation/ folder:  source("stan/run_known_exposure.R")
library(rstan)
rstan_options(auto_write = TRUE)

source("gen_data.R")
source("estimators.R")
source("stan/prep_data.R")
source("stan/plot_helpers.R")

stan_file <- "stan/gcomp_known_exposure.stan"
a_grid <- seq(5, 15, by = 1)
# Spline degrees of freedom per scenario: 1 is a straight line, 5 allows curves 
df_by_form <- c(linear = 1, quadratic = 5, nonlinear = 5)

mod <- stan_model(stan_file)
results <- list()

for (form in names(df_by_form)) {
  sim <- data.generate(seed = 1006092577, out_form = form)   # the base seed (replicate r of sim_replicates.R uses base_seed + r)
  dat <- make_stan_data(sim, exposure = "oracle", df = df_by_form[[form]], a_grid = a_grid)
  fit <- sampling(mod, data = dat, chains = 4, iter = 2000, cores = 2, seed = 1006092577, refresh = 0)

  # Posterior draws of theta(a): draws x grid
  theta <- rstan::extract(fit, "theta")$theta

  # Compare with the truth and with the BART oracle fit to the same dataset
  invisible(capture.output(bart <- bart_gcomp(sim, a_grid = a_grid, use_true_A = TRUE)))
  comp <- data.frame(
    a = a_grid,
    truth_pop = true_erf_pop(a_grid, form),
    truth_sample = true_erf(a_grid, sim),
    stan = colMeans(theta),
    stan_lo = apply(theta, 2, quantile, 0.025),
    stan_hi = apply(theta, 2, quantile, 0.975),
    bart = bart$theta_summary$mean
  )
  results[[form]] <- list(comp = comp, fit = fit)

  cat("\n", toupper(form), " (df =", df_by_form[[form]], ")\n")
  print(round(comp, 4), row.names = FALSE)
  cat("max relative difference from the dataset's truth: Stan",
      round(100 * max(abs(comp$stan / comp$truth_sample - 1)), 1), "%, BART",
      round(100 * max(abs(comp$bart / comp$truth_sample - 1)), 1), "%\n")
}

# Plot: posterior mean and 95% band from Stan, the BART oracle, and the two truths
est <- do.call(rbind, lapply(names(results), function(f) {
  cp <- results[[f]]$comp
  rbind(data.frame(scenario = f, a = cp$a, method = "Stan", mean = cp$stan, lo = cp$stan_lo, hi = cp$stan_hi),
        data.frame(scenario = f, a = cp$a, method = "BART", mean = cp$bart, lo = NA, hi = NA))
}))
truth <- do.call(rbind, lapply(names(results), function(f) {
  cp <- results[[f]]$comp
  rbind(data.frame(scenario = f, a = cp$a, truth = "population truth", value = cp$truth_pop),
        data.frame(scenario = f, a = cp$a, truth = "this dataset's truth", value = cp$truth_sample))
}))
est$scenario <- factor(est$scenario, levels = names(df_by_form))
truth$scenario <- factor(truth$scenario, levels = names(df_by_form))
p <- plot_erf(est, truth, title = "Known exposure: oracle Poisson model in Stan",
              subtitle = "Posterior mean and 95% band (Stan), BART oracle, and the true curves")
dir.create("output/stan", recursive = TRUE, showWarnings = FALSE)
ggsave("output/stan/known_exposure_erf.png", p, width = 11, height = 4.2, dpi = 130)
print(p)
