# Latent exposure: the true exposure A is unknown; only the surrogate Z is observed, with a known error model.
# Test data: no spatial effect, no slope, constant error variance, so the fitted model is exactly right.
# Run from the Simulation/ folder:  source("stan/scripts/run_latent_exposure.R")
library(rstan)
rstan_options(auto_write = TRUE)

source("gen_data.R")
source("stan/R/prep_data.R"); source("stan/R/compile_model.R")
source("stan/R/plot_helpers.R")

latent_model <- "gcomp_latent_exposure"      # Stan models live in stan/models/
known_model <- "gcomp_known_exposure"
a_grid <- seq(5, 15, by = 1)
sigma_Z <- 1                   # error SD (variance 1, reliability about 0.82), also the known value passed to Stan

sim <- data.generate(seed = 1006092577, out_form = "linear",
                     rho_A = 0, bZ1 = 1, bZ2 = 0, sigma_Z = sigma_Z)

fit_draws <- function(model_name, dat, par, init = "random") {
  fit <- sampling(compile_model(model_name), data = dat, chains = 4, iter = 2000, cores = 2, seed = 1006092577, refresh = 0, init = init)
  list(fit = fit, theta = rstan::extract(fit, par)$theta)
}
oracle <- fit_draws(known_model, make_stan_data(sim, "oracle", df = 1, a_grid = a_grid), "theta")  # true A
naive <- fit_draws(known_model, make_stan_data(sim, "naive", df = 1, a_grid = a_grid), "theta")    # Z as if it were A
# Latent model: the latent exposures start at the surrogate Z (random starts around 0 give chains that disagree)
dat_latent <- make_stan_data_latent(sim, sigma_Z, a_grid)
init_latent <- function() list(beta0 = -3, beta_A = 0, beta_X = rep(0, 4), gamma0 = mean(sim$Zi),
                               gamma = rep(0, 4), sigma_A = 1, A = dat_latent$Z)
latent <- fit_draws(latent_model, dat_latent, "theta", init = init_latent)

summ <- function(th) data.frame(mean = colMeans(th), lo = apply(th, 2, quantile, 0.025), hi = apply(th, 2, quantile, 0.975))
s_or <- summ(oracle$theta); s_nv <- summ(naive$theta); s_lt <- summ(latent$theta)
comp <- data.frame(a = a_grid, truth_pop = true_erf_pop(a_grid, "linear"), truth_sample = true_erf(a_grid, sim),
                   oracle = s_or$mean, naive = s_nv$mean, latent = s_lt$mean,
                   latent_lo = s_lt$lo, latent_hi = s_lt$hi)
print(round(comp, 4), row.names = FALSE)
cat("max relative difference from the dataset's truth: oracle",
    round(100 * max(abs(comp$oracle / comp$truth_sample - 1)), 1), "%, naive",
    round(100 * max(abs(comp$naive / comp$truth_sample - 1)), 1), "%, latent",
    round(100 * max(abs(comp$latent / comp$truth_sample - 1)), 1), "%\n")

# Plot 1: the exposure-response curves
# Only the latent fit gets a band; oracle and naive are lines (lo/hi = NA) to keep the plot readable
est <- rbind(data.frame(a = a_grid, method = "Oracle", s_or["mean"], lo = NA, hi = NA),
             data.frame(a = a_grid, method = "Naive", s_nv["mean"], lo = NA, hi = NA),
             data.frame(a = a_grid, method = "Latent", s_lt))
truth <- rbind(data.frame(a = a_grid, truth = "population truth", value = comp$truth_pop),
               data.frame(a = a_grid, truth = "this dataset's truth", value = comp$truth_sample))
p1 <- plot_erf(est, truth, title = "Latent exposure with a known error model",
               subtitle = "Linear scenario, no spatial effect, constant error variance (reliability about 0.82)")
dir.create("output/stan", recursive = TRUE, showWarnings = FALSE)
ggsave("output/stan/latent_exposure_erf.png", p1, width = 7, height = 5, dpi = 130)
print(p1)

# Plot 2: how well the latent exposures are recovered (needs the model to save A as a transformed parameter)
if ("A" %in% latent$fit@model_pars) {
  A_post <- colMeans(rstan::extract(latent$fit, "A")$A)
  d <- rbind(data.frame(panel = "Surrogate Z", true_A = sim$Ai, value = sim$Zi),
             data.frame(panel = "Posterior mean of A (latent model)", true_A = sim$Ai, value = A_post))
  p2 <- ggplot(d, aes(true_A, value)) + geom_point(alpha = 0.4) + geom_abline(colour = "red") +
    facet_wrap(~panel) + labs(x = "true exposure A", y = NULL, title = "Recovery of the latent exposures",
                              subtitle = "Red line: estimate = truth") + theme_bw(base_size = 12)
  ggsave("output/stan/latent_exposure_recovery.png", p2, width = 9, height = 4.2, dpi = 130)
  print(p2)
}
