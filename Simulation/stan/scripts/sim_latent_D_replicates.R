# Replicates for the latent exposure model with the Option A error model (D shifts the surrogate) and the CAR
# spatial effect, on the default data 
# Fits per replicate: oracle (known-exposure model), naive (Z instead of A), latent, known error parameters known, 
# latent, estimated   beta_Z1, beta_Z2, sigma_Z estimated, with CAR
library(rstan)
rstan_options(auto_write = TRUE)
source("gen_data.R"); source("stan/R/prep_data.R"); source("stan/R/compile_model.R")

args <- commandArgs(trailingOnly = TRUE)
n_reps <- if (length(args) >= 1) as.integer(args[1]) else 20
n_iter <- if (length(args) >= 2) as.integer(args[2]) else 2000
base_seed <- 1006092577
n_areas <- 100      # 10 x 10 grid
a_grid <- seq(5, 15, by = 1)
set.seed(base_seed)
rep_seeds <- sample.int(.Machine$integer.max, n_reps)

known_mod <- compile_model("gcomp_latent_exposure_D")
est_mod <- compile_model("gcomp_latent_exposure_D_priors")   # CAR takes a long time
exposure_mod <- compile_model("gcomp_known_exposure")
pop_truth <- true_erf_pop(a_grid, "linear")

fit_it <- function(mod, dat, seed, init = "random") {
  fit <- sampling(mod, data = dat, chains = 4, iter = n_iter, cores = 4, seed = seed, refresh = 0, init = init,
                  control = list(adapt_delta = 0.9))
  sp <- get_sampler_params(fit, inc_warmup = FALSE)
  list(fit = fit, theta = rstan::extract(fit, "theta")$theta,
       div = sum(sapply(sp, function(m) sum(m[, "divergent__"]))),
       rhat = max(summary(fit)$summary[, "Rhat"], na.rm = TRUE))
}

run_one <- function(r) {
  s <- rep_seeds[r]
  sim <- data.generate(seed = s, n_areas = n_areas, out_form = "linear", bZ1 = 1, bZ2 = 0.5, sigma_Z = 1)
  d_known <- make_stan_data_latent_D(sim, bZ1 = 1, bZ2 = 0.5, sigma_Z = 1, a_grid = a_grid)
  d_est <- make_stan_data_latent_D_priors(sim, a_grid)
  init_known <- function() list(beta0 = -3, beta_A = 0, beta_X = rep(0, 4), gamma0 = mean(sim$Zi), gamma = rep(0, 4),
                                sigma_A = 1, A = sim$Zi)
  init_est <- function() c(init_known(), list(beta_Z1 = 1, beta_Z2 = 0, sigma_Z = 1, phi = rep(0, nrow(sim)), sigma_phi = 1, rho = 0.5))
  fits <- list(
    "oracle" = fit_it(exposure_mod, make_stan_data(sim, "oracle", df = 1, a_grid = a_grid), s),
    "naive" = fit_it(exposure_mod, make_stan_data(sim, "naive", df = 1, a_grid = a_grid), s),
    "latent, known errors" = fit_it(known_mod, d_known, s, init_known),
    "latent, estimated errors + CAR" = fit_it(est_mod, d_est, s, init_est))
  pars <- rstan::extract(fits[["latent, estimated errors + CAR"]]$fit, c("beta_Z1", "beta_Z2", "sigma_Z", "rho"))
  list(summary = do.call(rbind, lapply(names(fits), function(m) {
         th <- fits[[m]]$theta
         data.frame(rep = r, method = m, a = a_grid, true = pop_truth, true_sample = true_erf(a_grid, sim),
                    est = colMeans(th), lo = apply(th, 2, quantile, 0.025), hi = apply(th, 2, quantile, 0.975))})),
       diag = data.frame(rep = r, method = names(fits), div = sapply(fits, `[[`, "div"), max_rhat = sapply(fits, `[[`, "rhat")),
       error_pars = data.frame(rep = r, beta_Z1 = mean(pars$beta_Z1), beta_Z2 = mean(pars$beta_Z2), sigma_Z = mean(pars$sigma_Z), rho = mean(pars$rho)))
}

t0 <- Sys.time()
res <- vector("list", n_reps)
for (r in seq_len(n_reps)) {
  res[[r]] <- run_one(r)
  cat("replicate", r, "of", n_reps, "done (", round(difftime(Sys.time(), t0, units = "mins"), 1), "min elapsed )\n")
}
out <- list(summary = do.call(rbind, lapply(res, `[[`, "summary")), diag = do.call(rbind, lapply(res, `[[`, "diag")),
            error_pars = do.call(rbind, lapply(res, `[[`, "error_pars")), n_iter = n_iter)
dir.create("output/stan", recursive = TRUE, showWarnings = FALSE)
saveRDS(out, "output/stan/latent_D_replicates.rds")

raw <- out$summary
summ <- do.call(rbind, lapply(split(raw, list(raw$method, raw$a), drop = TRUE), function(d)
  data.frame(method = d$method[1], a = d$a[1], rel_bias_pct = 100 * mean(d$est - d$true) / d$true[1],
             rmse = sqrt(mean((d$est - d$true)^2)), coverage = mean(d$lo <= d$true & d$true <= d$hi))))
tab <- aggregate(cbind(abs_rel_bias_pct = abs(rel_bias_pct), rmse, coverage) ~ method, summ, mean)
cat("\nAveraged over the grid,", n_reps, "replicates, coverage against the population truth:\n"); print(tab, digits = 3)
cat("\nSampler diagnostics: replicates with any divergence, and worst / median of each replicate's max R-hat:\n")
print(do.call(rbind, lapply(split(out$diag, out$diag$method), function(d) data.frame(method = d$method[1], replicates_with_divergence = sum(d$div > 0),
      worst_rhat = max(d$max_rhat), median_rhat = median(d$max_rhat)))), digits = 3, row.names = FALSE)
cat("\nError parameters, posterior means across replicates (truth: beta_Z1 = 1, beta_Z2 = 0.5, sigma_Z = 1, rho = 0.9):\n")
print(round(sapply(out$error_pars[, -1], function(x) c(mean = mean(x), sd = sd(x))), 3))
write.csv(summ, "output/stan/latent_D_replicates_summary.csv", row.names = FALSE)
