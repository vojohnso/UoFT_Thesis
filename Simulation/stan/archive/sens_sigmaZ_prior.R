# Prior sensitivity for sigma_Z in the anchored latent model (gcomp_latent_exposure_D_anchored_sens.stan):
# same dataset as run_latent_exposure_D_anchored.R, sigma_Z ~ normal(m, 0.3) for several centres m (truth = 1)

library(rstan)
source("gen_data.R"); source("stan/R/prep_data.R"); source("stan/R/compile_model.R")

a_grid <- seq(5, 15, by = 1)
sim <- data.generate(seed = 1006092577, n_areas = 100, out_form = "linear", spatial_in = "outcome", bZ1 = 1, bZ2 = 0.5, sigma_Z = 1)
dat <- make_stan_data_latent_D_priors(sim, a_grid)
init <- function() list(beta0 = -3, beta_A = 0, beta_X = rep(0, 4), gamma0 = mean(dat$Z), gamma_A = rep(0, 4), sigma_A = 1,
                        A = dat$Z, beta_Z2 = 0, sigma_Z = 1, phi = rep(0, dat$N), sigma_phi = 1, rho = 0.5)
model <- compile_model("gcomp_latent_exposure_D_anchored_sens")
pars <- c("sigma_A", "sigma_phi", "rho", "beta_A", "beta_Z2", "sigma_Z")

res <- do.call(rbind, lapply(c(0.5, 0.7, 1, 1.5, 2), function(m) {
  fit <- sampling(model, data = c(dat, list(sigmaZ_prior_mean = m, sigmaZ_prior_sd = 0.3)), chains = 4, iter = 3000, cores = 4,
                  seed = 1006092577, refresh = 0, init = init, control = list(adapt_delta = 0.9))
  s <- summary(fit, pars = pars)$summary
  sp <- get_sampler_params(fit, inc_warmup = FALSE)
  data.frame(prior_mean = m, sigma_Z = s["sigma_Z", "mean"], sigma_Z_sd = s["sigma_Z", "sd"], sigma_A = s["sigma_A", "mean"],
             beta_A = s["beta_A", "mean"], max_Rhat = max(s[, "Rhat"]), min_neff = min(s[, "n_eff"]),
             divergences = sum(sapply(sp, function(x) sum(x[, "divergent__"]))))
}))
cat("\nTruth: sigma_Z = 1, sigma_A = 1.5, beta_A = 0.1; prior sd = 0.3\n")
print(round(res, 3), row.names = FALSE)
dir.create("output/stan/latent_D_car", recursive = TRUE, showWarnings = FALSE)
write.csv(res, "output/stan/latent_D_car/sens_sigmaZ_prior.csv", row.names = FALSE)
