# Diagnostic: the reduced latent model with tight priors around the true value on EVERY
# parameter (Stan model gcomp_latent_exposure_D_reduced_sparse_tight). Same 20 seeds as the
# other reduced studies.
#
# Question: once every parameter is nearly known, does anything other than sigma_Z explain
# the bias? (It did not: this is not better than fixing sigma_Z alone.)

library(rstan)
source("gen_data.R")
source("stan/R/prep_data.R")
source("stan/R/compile_model.R")
source("stan/R/diagnostics.R")

# Settings
# Optional arguments: number of replicates and iterations
args <- commandArgs(trailingOnly = TRUE)
n_reps <- if (length(args) >= 1) as.integer(args[1]) else 20
n_iter <- if (length(args) >= 2) as.integer(args[2]) else 2000

base_seed <- 1006092577
n_areas <- 100
a_grid <- seq(5, 15, by = 1)

set.seed(base_seed)
rep_seeds <- sample.int(.Machine$integer.max, n_reps)

red_mod <- compile_model("gcomp_latent_exposure_D_reduced_sparse_tight")

# Population truth, with the spatial effect in the outcome
pop_truth <- true_erf_pop(a_grid, "linear", spatial_in = "outcome", n_areas = n_areas)


# Fit one Stan model; keep the fit, theta and the convergence summary
fit_it <- function(mod, dat, seed, init = "random") {
  fit <- sampling(
    mod, data = dat, chains = 4, iter = n_iter, cores = 4, seed = seed,
    refresh = 0, init = init, control = list(adapt_delta = 0.9)
  )
  list(
    fit = fit,
    theta = rstan::extract(fit, "theta")$theta,
    conv = convergence_summary(fit)
  )
}

# Fit the model on replicate r
run_one <- function(r) {
  s <- rep_seeds[r]
  sim <- data.generate(
    seed = s, n_areas = n_areas, out_form = "linear", spatial_in = "outcome",
    bZ1 = 1, bZ2 = 0.5, sigma_Z = 1
  )
  d_red <- make_stan_data_latent_D_priors(sim, a_grid)

  # Starting values: A starts at the observed surrogate Z
  init_known <- function() {
    list(
      beta0 = -3, beta_A = 0, beta_X = rep(0, 4),
      gamma0 = mean(sim$Zi), gamma_A = rep(0, 4), sigma_A = 1,
      A = sim$Zi,
      phi = rep(0, nrow(sim)), sigma_phi = 1, rho = 0.5
    )
  }
  init_red <- function() c(init_known(), list(beta_Z2 = 0, sigma_Z = 1))

  fits <- list("reduced, tight priors on everything" = fit_it(red_mod, d_red, s, init_red))

  list(
    # Posterior mean and 95% interval of theta(a) at each grid value, per method
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
    # Convergence of each fit
    diag = do.call(rbind, lapply(names(fits), function(m) {
      data.frame(rep = r, method = m, fits[[m]]$conv)
    })),
    # Posterior means of the error parameters
    error_pars = do.call(rbind, lapply(names(fits), function(m) {
      p <- rstan::extract(fits[[m]]$fit, c("beta_Z2", "sigma_Z", "rho"))
      data.frame(
        rep = r, method = m,
        beta_Z2 = mean(p$beta_Z2), sigma_Z = mean(p$sigma_Z), rho = mean(p$rho)
      )
    }))
  )
}


# Run and save
t0 <- Sys.time()
res <- list()
for (r in seq_len(n_reps)) {
  res[[r]] <- run_one(r)
  cat("replicate", r, "of", n_reps, "done,", round(difftime(Sys.time(), t0, units = "mins"), 1), "min elapsed\n")
}

out <- list(
  summary = do.call(rbind, lapply(res, `[[`, "summary")),
  diag = do.call(rbind, lapply(res, `[[`, "diag")),
  error_pars = do.call(rbind, lapply(res, `[[`, "error_pars"))
)
dir.create("output/stan/latent_D_car", recursive = TRUE, showWarnings = FALSE)
saveRDS(out, "output/stan/latent_D_car/latent_D_reduced_tight_replicates.rds")
