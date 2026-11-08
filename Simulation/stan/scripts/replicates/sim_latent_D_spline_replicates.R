# Replicates for the spline version of the anchored latent model (natural cubic spline in A, df = 3) on a quadratic (default) or
# linear exposure-response, with the CAR spatial effect on the outcome. Same replicate seeds as sim_latent_D_anchored_replicates.R
# Fits per replicate: oracle (true A), naive (Z), latent spline + CAR
library(rstan)
source("gen_data.R"); source("stan/R/prep_data.R"); source("stan/R/compile_model.R")

args <- commandArgs(trailingOnly = TRUE)
n_reps <- if (length(args) >= 1) as.integer(args[1]) else 5
n_iter <- if (length(args) >= 2) as.integer(args[2]) else 1500
out_form <- if (length(args) >= 3) args[3] else "quadratic"
sfx <- if (out_form == "quadratic") "" else paste0("_", out_form)   # keeps the quadratic output names
base_seed <- 1006092577
n_areas <- 100
a_grid <- seq(5, 15, by = 1)
set.seed(base_seed)
rep_seeds <- sample.int(.Machine$integer.max, n_reps)

est_mod <- compile_model("gcomp_latent_exposure_D_anchored_spline")
exposure_mod <- compile_model("gcomp_known_exposure")
pop_truth <- true_erf_pop(a_grid, out_form, spatial_in = "outcome", n_areas = n_areas)

fit_it <- function(mod, dat, seed, init = "random") {
  fit <- sampling(mod, data = dat, chains = 4, iter = n_iter, cores = 4, seed = seed, refresh = 0, init = init,
                  control = list(adapt_delta = 0.9))
  sp <- get_sampler_params(fit, inc_warmup = FALSE)
  list(theta = rstan::extract(fit, "theta")$theta,
       div = sum(sapply(sp, function(m) sum(m[, "divergent__"]))),
       treedepth = sum(sapply(sp, function(m) sum(m[, "treedepth__"] >= 10))),
       rhat = max(summary(fit)$summary[, "Rhat"], na.rm = TRUE))
}

run_one <- function(r) {
  s <- rep_seeds[r]
  sim <- data.generate(seed = s, n_areas = n_areas, out_form = out_form, spatial_in = "outcome", bZ1 = 1, bZ2 = 0.5, sigma_Z = 1)
  d_est <- make_stan_data_latent_D_spline(sim, a_grid)
  init_est <- function() list(beta0 = -3, beta_A = rep(0, 3), beta_X = rep(0, 4), gamma0 = mean(sim$Zi), gamma_A = rep(0, 4),
                              sigma_A = 1, A = sim$Zi, beta_Z2 = 0, sigma_Z = 1, phi = rep(0, nrow(sim)), sigma_phi = 1, rho = 0.5)
  fits <- list(
    "oracle" = fit_it(exposure_mod, make_stan_data(sim, "oracle", df = 3, a_grid = a_grid), s),
    "naive" = fit_it(exposure_mod, make_stan_data(sim, "naive", df = 3, a_grid = a_grid), s),
    "latent spline + CAR" = fit_it(est_mod, d_est, s, init_est))
  list(summary = do.call(rbind, lapply(names(fits), function(m) {
         th <- fits[[m]]$theta
         data.frame(rep = r, method = m, a = a_grid, true = pop_truth, est = colMeans(th),
                    lo = apply(th, 2, quantile, 0.025), hi = apply(th, 2, quantile, 0.975))})),
       diag = data.frame(rep = r, method = names(fits), div = sapply(fits, `[[`, "div"),
                         treedepth = sapply(fits, `[[`, "treedepth"), max_rhat = sapply(fits, `[[`, "rhat")))
}

t0 <- Sys.time()
res <- vector("list", n_reps)
dir.create("output/stan/latent_D_car", recursive = TRUE, showWarnings = FALSE)
for (r in seq_len(n_reps)) {
  res[[r]] <- run_one(r)
  saveRDS(res[seq_len(r)], paste0("output/stan/latent_D_car/latent_D_spline_replicates", sfx, "_partial.rds"))
  cat("replicate", r, "of", n_reps, "done (", round(difftime(Sys.time(), t0, units = "mins"), 1), "min elapsed )\n")
}
out <- list(summary = do.call(rbind, lapply(res, `[[`, "summary")), diag = do.call(rbind, lapply(res, `[[`, "diag")), n_iter = n_iter)
saveRDS(out, paste0("output/stan/latent_D_car/latent_D_spline_replicates", sfx, ".rds"))

raw <- out$summary
summ <- do.call(rbind, lapply(split(raw, list(raw$method, raw$a), drop = TRUE), function(d)
  data.frame(method = d$method[1], a = d$a[1], rel_bias_pct = 100 * mean(d$est - d$true) / d$true[1],
             rmse = sqrt(mean((d$est - d$true)^2)), coverage = mean(d$lo <= d$true & d$true <= d$hi))))
tab <- aggregate(cbind(abs_rel_bias_pct = abs(rel_bias_pct), rmse, coverage) ~ method, summ, mean)
cat("\nAveraged over the grid,", n_reps, "replicates,", out_form, "truth, coverage against the population truth:\n"); print(tab, digits = 3)
cat("\nDiagnostics per method: replicates with divergences, with tree-depth hits, worst / median max R-hat:\n")
print(do.call(rbind, lapply(split(out$diag, out$diag$method), function(d) data.frame(method = d$method[1],
      with_divergence = sum(d$div > 0), with_treedepth = sum(d$treedepth > 0), worst_rhat = max(d$max_rhat), median_rhat = median(d$max_rhat)))),
      digits = 3, row.names = FALSE)
write.csv(summ, paste0("output/stan/latent_D_car/latent_D_spline_replicates", sfx, "_summary.csv"), row.names = FALSE)
