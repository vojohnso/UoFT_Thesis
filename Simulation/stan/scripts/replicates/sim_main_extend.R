# Main study scale-up: replicates first..last (arguments, default 21..100) of the core methods, with the same master seed so replicates 1..20 are the
# existing ones. Linear truth, 100 areas, spatial effect in the outcome. Six fits per replicate:
#   oracle, naive, known errors (all with CAR); reduced latent model with sigma_Z fixed at 1, with sigma_Z ~ normal(1, 0.3), with the wild guess normal(0, 5).
# Saves after every replicate to output/stan/latent_D_car/main_reps_<first>_<last>.rds so partial results survive an interruption.
library(rstan)
source("gen_data.R"); source("stan/R/prep_data.R"); source("stan/R/compile_model.R")
args <- commandArgs(trailingOnly = TRUE)
first <- as.integer(args[1]); last <- as.integer(args[2]); n_iter <- 2000
base_seed <- 1006092577; n_areas <- 100; a_grid <- seq(5, 15, by = 1)
set.seed(base_seed); rep_seeds <- sample.int(.Machine$integer.max, 100)
exposure_mod <- compile_model("gcomp_known_exposure_car"); known_mod <- compile_model("gcomp_latent_exposure_D_known_car")
red_mod <- compile_model("gcomp_latent_exposure_D_reduced_sparse")
pop_truth <- true_erf_pop(a_grid, "linear", spatial_in = "outcome", n_areas = n_areas)
out_file <- sprintf("output/stan/latent_D_car/main_reps_%d_%d.rds", first, last)

fit_it <- function(mod, dat, seed, init = "random") {
  fit <- sampling(mod, data = dat, chains = 4, iter = n_iter, cores = 4, seed = seed, refresh = 0, init = init, control = list(adapt_delta = 0.9))
  sp <- get_sampler_params(fit, inc_warmup = FALSE)
  list(theta = rstan::extract(fit, "theta")$theta, div = sum(sapply(sp, function(m) sum(m[, "divergent__"]))),
       rhat = max(summary(fit)$summary[, "Rhat"], na.rm = TRUE))
}
run_one <- function(r) {
  s <- rep_seeds[r]
  sim <- data.generate(seed = s, n_areas = n_areas, out_form = "linear", spatial_in = "outcome", bZ1 = 1, bZ2 = 0.5, sigma_Z = 1)
  car <- car_structure(nrow(sim))
  d_known <- c(make_stan_data_latent_D(sim, bZ1 = 1, bZ2 = 0.5, sigma_Z = 1, a_grid = a_grid), car)
  d_red <- make_stan_data_latent_D_priors(sim, a_grid)
  init_known <- function() list(beta0 = -3, beta_A = 0, beta_X = rep(0, 4), gamma0 = mean(sim$Zi), gamma_A = rep(0, 4), sigma_A = 1, A = sim$Zi,
                                phi = rep(0, nrow(sim)), sigma_phi = 1, rho = 0.5)
  init_red <- function(sz) function() c(init_known(), list(beta_Z2 = 0, sigma_Z = sz))
  fits <- list(
    "oracle" = fit_it(exposure_mod, c(make_stan_data(sim, "oracle", df = 1, a_grid = a_grid), car), s),
    "naive" = fit_it(exposure_mod, c(make_stan_data(sim, "naive", df = 1, a_grid = a_grid), car), s),
    "latent, known errors" = fit_it(known_mod, d_known, s, init_known),
    "reduced, sigma_Z fixed at 1" = fit_it(red_mod, c(d_red, list(sigma_Z_mean = 1, sigma_Z_sd = 0.001)), s, init_red(1)),
    "reduced, sigma_Z ~ normal(1, 0.3)" = fit_it(red_mod, c(d_red, list(sigma_Z_mean = 1, sigma_Z_sd = 0.3)), s, init_red(1)),
    "reduced, sigma_Z wild guess" = fit_it(red_mod, c(d_red, list(sigma_Z_mean = 0, sigma_Z_sd = 5)), s, init_red(1)))
  list(summary = do.call(rbind, lapply(names(fits), function(m) { th <- fits[[m]]$theta
         data.frame(rep = r, method = m, a = a_grid, true = pop_truth, true_sample = true_erf(a_grid, sim), est = colMeans(th),
                    lo = apply(th, 2, quantile, 0.025), hi = apply(th, 2, quantile, 0.975)) })),
       diag = data.frame(rep = r, method = names(fits), div = sapply(fits, `[[`, "div"), max_rhat = sapply(fits, `[[`, "rhat")))
}
res <- list(); t0 <- Sys.time()
for (r in first:last) {
  res[[as.character(r)]] <- run_one(r)
  saveRDS(list(summary = do.call(rbind, lapply(res, `[[`, "summary")), diag = do.call(rbind, lapply(res, `[[`, "diag")), n_iter = n_iter), out_file)
  cat("replicate", r, "(", r - first + 1, "of", last - first + 1, ") done,", round(difftime(Sys.time(), t0, units = "mins"), 1), "min elapsed\n")
}
