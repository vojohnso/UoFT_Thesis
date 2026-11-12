# Diagnostic: the reduced latent model with sigma_Z fixed at its true value (1) through a very
# tight prior, normal(1, 0.001). Same 20 seeds as sim_latent_D_reduced_replicates.R.
#
# Output: output/stan/latent_D_car/latent_D_reduced_fixedsigmaZ_replicates.rds (+ a summary csv)
# Usage:  Rscript stan/scripts/replicates/sim_latent_D_reduced_fixedsigmaZ.R [n_reps] [n_iter]

library(rstan)
source("gen_data.R")
source("stan/R/prep_data.R")
source("stan/R/compile_model.R")


# Settings
args <- commandArgs(trailingOnly = TRUE)
n_reps <- if (length(args) >= 1) as.integer(args[1]) else 20
n_iter <- if (length(args) >= 2) as.integer(args[2]) else 2000

base_seed <- 1006092577
n_areas <- 100 # 10 x 10 grid
a_grid <- seq(5, 15, by = 1)

set.seed(base_seed)
rep_seeds <- sample.int(.Machine$integer.max, n_reps)

est_mod <- compile_model("gcomp_latent_exposure_D_reduced_sparse")

# Population truth, with the spatial effect in the outcome
pop_truth <- true_erf_pop(a_grid, "linear", spatial_in = "outcome", n_areas = n_areas)


# Helpers
# Fit one Stan model; keep the fit, the posterior draws of theta and two convergence summaries
fit_it <- function(mod, dat, seed, init = "random") {
  fit <- sampling(
    mod, data = dat, chains = 4, iter = n_iter, cores = 4, seed = seed,
    refresh = 0, init = init, control = list(adapt_delta = 0.9)
  )
  sampler_params <- get_sampler_params(fit, inc_warmup = FALSE)
  list(
    fit = fit,
    theta = rstan::extract(fit, "theta")$theta,
    div = sum(sapply(sampler_params, function(m) sum(m[, "divergent__"]))),
    rhat = max(summary(fit)$summary[, "Rhat"], na.rm = TRUE)
  )
}

# Fit the reduced model with sigma_Z fixed at 1 on replicate r
run_one <- function(r) {
  s <- rep_seeds[r]
  sim <- data.generate(
    seed = s, n_areas = n_areas, out_form = "linear", spatial_in = "outcome",
    bZ1 = 1, bZ2 = 0.5, sigma_Z = 1
  )
  d_est <- make_stan_data_latent_D_priors(sim, a_grid)

  # Starting values: A starts at the observed surrogate Z
  init_known <- function() {
    list(
      beta0 = -3, beta_A = 0, beta_X = rep(0, 4),
      gamma0 = mean(sim$Zi), gamma_A = rep(0, 4), sigma_A = 1,
      A = sim$Zi
    )
  }
  init_est <- function() {
    c(init_known(), list(
      beta_Z2 = 0, sigma_Z = 1,
      phi = rep(0, nrow(sim)), sigma_phi = 1, rho = 0.5
    ))
  }

  fits <- list(
    "reduced, sigma_Z fixed at 1" = fit_it(
      est_mod, c(d_est, list(sigma_Z_mean = 1, sigma_Z_sd = 0.001)), s, init_est
    )
  )

  # Posterior means of the error parameters
  ep <- function(m) {
    p <- rstan::extract(fits[[m]]$fit, c("beta_Z2", "sigma_Z", "rho"))
    data.frame(
      rep = r, method = m,
      beta_Z2 = mean(p$beta_Z2), sigma_Z = mean(p$sigma_Z), rho = mean(p$rho)
    )
  }

  list(
    # Posterior mean and 95% interval of theta(a) at each grid value
    summary = do.call(rbind, lapply(names(fits), function(m) {
      th <- fits[[m]]$theta
      data.frame(
        rep = r, method = m, a = a_grid,
        true = pop_truth, true_sample = true_erf(a_grid, sim),
        est = colMeans(th),
        lo = apply(th, 2, quantile, 0.025),
        hi = apply(th, 2, quantile, 0.975)
      )
    })),
    # Divergences and max R-hat
    diag = data.frame(
      rep = r, method = names(fits),
      div = sapply(fits, `[[`, "div"),
      max_rhat = sapply(fits, `[[`, "rhat")
    ),
    error_pars = ep(names(fits)[1])
  )
}


# Run and save
t0 <- Sys.time()
res <- vector("list", n_reps)
for (r in seq_len(n_reps)) {
  res[[r]] <- run_one(r)
  cat(
    "replicate", r, "of", n_reps, "done (",
    round(difftime(Sys.time(), t0, units = "mins"), 1), "min elapsed )\n"
  )
}

out <- list(
  summary = do.call(rbind, lapply(res, `[[`, "summary")),
  diag = do.call(rbind, lapply(res, `[[`, "diag")),
  error_pars = do.call(rbind, lapply(res, `[[`, "error_pars")),
  n_iter = n_iter
)
dir.create("output/stan/latent_D_car", recursive = TRUE, showWarnings = FALSE)
saveRDS(out, "output/stan/latent_D_car/latent_D_reduced_fixedsigmaZ_replicates.rds")


# Summary printed to the console, and saved
raw <- out$summary

# At each exposure value: relative bias (%), RMSE and coverage over the replicates
summ <- do.call(rbind, lapply(split(raw, list(raw$method, raw$a), drop = TRUE), function(d) {
  data.frame(
    method = d$method[1],
    a = d$a[1],
    rel_bias_pct = 100 * mean(d$est - d$true) / d$true[1],
    rmse = sqrt(mean((d$est - d$true)^2)),
    coverage = mean(d$lo <= d$true & d$true <= d$hi)
  )
}))

# Average over the grid, per method
tab <- aggregate(cbind(abs_rel_bias_pct = abs(rel_bias_pct), rmse, coverage) ~ method, summ, mean)

cat("\nAveraged over the grid,", n_reps, "replicates, coverage against the population truth:\n")
print(tab, digits = 3)

cat("\nSampler diagnostics: replicates with any divergence, and worst / median of each replicate's max R-hat:\n")
print(
  do.call(rbind, lapply(split(out$diag, out$diag$method), function(d) {
    data.frame(
      method = d$method[1],
      replicates_with_divergence = sum(d$div > 0),
      worst_rhat = max(d$max_rhat),
      median_rhat = median(d$max_rhat)
    )
  })),
  digits = 3, row.names = FALSE
)

cat("\nError parameters, posterior means across replicates (truth: beta_Z2 = 0.5, sigma_Z = 1, rho = 0.9):\n")
print(aggregate(cbind(beta_Z2, sigma_Z, rho) ~ method, out$error_pars, function(x) {
  round(c(mean = mean(x), sd = sd(x)), 3)
}))

write.csv(summ, "output/stan/latent_D_car/latent_D_reduced_fixedsigmaZ_replicates_summary.csv", row.names = FALSE)
