# Latent exposure model with the Option A error model (D shifts the surrogate), on the default spatial data.
# Fit 1: error parameters known (no spatial term in the model).
# Fit 2: beta_Z1, beta_Z2 and sigma_Z all estimated (priors are in the Stan file), with the CAR spatial effect.
# beta_Z1 sets the scale of A and is not identified, so expect chains that disagree for it.
# Run from the Simulation/ folder:  source("stan/scripts/run_latent_exposure_D.R")
library(rstan)
rstan_options(auto_write = TRUE)
source("gen_data.R"); source("stan/R/prep_data.R"); source("stan/R/compile_model.R")

a_grid <- seq(5, 15, by = 1)
n_areas <- 100      # 10 x 10 grid: the plain CAR prior is too slow at 400 areas
truth_pars <- c(bZ1 = 1, bZ2 = 0.5, sigma_Z = 1)
sim <- data.generate(seed = 1006092577, n_areas = n_areas, out_form = "linear",                # default spatial exposure (rho_A = 0.9, R2_phi = 0.2)
                     bZ1 = truth_pars["bZ1"], bZ2 = truth_pars["bZ2"], sigma_Z = truth_pars["sigma_Z"])
ts <- true_erf(a_grid, sim)

init_fun <- function(dat, free_error = FALSE) function() {
  Z_adj <- dat$Z - 0.5 * dat$D                       # undo a typical D shift; A starts near this
  out <- list(beta0 = -3, beta_A = 0, beta_X = rep(0, 4), gamma0 = mean(dat$Z), gamma = rep(0, 4), sigma_A = 1,
              A = Z_adj + rnorm(dat$N, 0, 0.1))
  if (free_error) out <- c(out, list(beta_Z1 = 1, beta_Z2 = 0, sigma_Z = 1, phi = rep(0, dat$N), sigma_phi = 1))
  out
}
fit_it <- function(model_name, dat, free_error) {
  sampling(compile_model(model_name), data = dat, chains = 4, iter = 2000, cores = 4, seed = 1006092577, refresh = 0,
       init = init_fun(dat, free_error), control = list(adapt_delta = 0.9))
}
report <- function(label, fit, pars) {
  sp <- get_sampler_params(fit, inc_warmup = FALSE)
  th <- colMeans(rstan::extract(fit, "theta")$theta)
  cat("\n==", label, "==  divergences:", sum(sapply(sp, function(m) sum(m[, "divergent__"]))),
      " max Rhat:", round(max(summary(fit)$summary[, "Rhat"], na.rm = TRUE), 3),
      " max error in theta:", round(100 * max(abs(th / ts - 1)), 1), "%\n")
  s <- summary(fit, pars = pars)$summary
  print(round(s[, c("mean", "sd", "2.5%", "97.5%", "n_eff", "Rhat")], 3))
}

# Fit 1: known error model
f_known <- fit_it("gcomp_latent_exposure_D",
                  make_stan_data_latent_D(sim, bZ1 = 1, bZ2 = 0.5, sigma_Z = 1, a_grid = a_grid), FALSE)
report("Known error model (beta_Z1 = 1, beta_Z2 = 0.5, sigma_Z = 1)", f_known, c("beta_A", "sigma_A"))

# Fit 2: error parameters estimated, with CAR
# plain CAR file: with 400 areas this is very slow (hours per chain); use a smaller n_areas in data.generate() to test
f_est <- fit_it("gcomp_latent_exposure_D_priors", make_stan_data_latent_D_priors(sim, a_grid), TRUE)
report("Error parameters estimated, with CAR", f_est,
       c("beta_Z1", "beta_Z2", "sigma_Z", "sigma_A", "sigma_phi", "beta_A"))
