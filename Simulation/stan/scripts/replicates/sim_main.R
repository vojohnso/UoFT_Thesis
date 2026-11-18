# Main study: replicates first..last (arguments, default 1..100) of the core methods.
#
# Setting: linear truth, 100 areas, spatial (CAR) effect in the outcome, true sigma_Z = 1.
#
# Seven fits per replicate, all on the same simulated dataset:
#   1. oracle                              (true exposure A, with CAR)
#   2. naive                               (surrogate Z used as if it were A, with CAR)
#   3. latent, known errors                (error model given as data, with CAR)
#   4. reduced, sigma_Z fixed at 1         (sigma_Z ~ normal(1, 0.001))
#   5. reduced, sigma_Z ~ normal(1, 0.3)   (informative prior around the truth)
#   6. reduced, sigma_Z wild guess         (sigma_Z ~ normal(0, 5), truncated at 0)
#   7. reduced, tight priors on everything (diagnostic: every parameter almost known)

library(rstan)
source("gen_data.R")
source("stan/R/prep_data.R")
source("stan/R/compile_model.R")
source("stan/R/diagnostics.R")

# Settings
# Optional arguments: first and last replicate (default 1 100)
args <- commandArgs(trailingOnly = TRUE)
first <- if (length(args) >= 1) as.integer(args[1]) else 1L
last <- if (length(args) >= 2) as.integer(args[2]) else 100L

base_seed <- 1006092577
n_areas <- 100
a_grid <- seq(5, 15, by = 1)

# One seed per replicate, drawn once from the master seed
set.seed(base_seed)
rep_seeds <- sample.int(.Machine$integer.max, 100)

# Models
exposure_mod <- compile_model("gcomp_known_exposure_car")
known_mod <- compile_model("gcomp_latent_exposure_D_known_car")
red_mod <- compile_model("gcomp_latent_exposure_D_reduced_sparse")
tight_mod <- compile_model("gcomp_latent_exposure_D_reduced_sparse_tight")

# Population truth of the ERF, the same in every replicate
pop_truth <- true_erf_pop(a_grid, "linear", spatial_in = "outcome", n_areas = n_areas)

out_file <- sprintf("output/stan/latent_D_car/main_reps_%d_%d.rds", first, last)

# Fit one Stan model (warm-up is half the iterations); keep theta and the convergence summary
fit_it <- function(mod, dat, seed, init = "random", iter = 8000) {
  fit <- sampling(
    mod, data = dat, chains = 4, iter = iter, warmup = iter / 2, cores = 4,
    seed = seed, refresh = 0, init = init, control = list(adapt_delta = 0.9)
  )
  list(
    theta = rstan::extract(fit, "theta")$theta,
    conv = data.frame(iter = iter, convergence_summary(fit))
  )
}

# Run all seven methods on replicate r
run_one <- function(r) {
  s <- rep_seeds[r]
  sim <- data.generate(
    seed = s, n_areas = n_areas, out_form = "linear", spatial_in = "outcome",
    bZ1 = 1, bZ2 = 0.5, sigma_Z = 1
  )
  car <- car_structure(nrow(sim))

  # Stan data
  d_oracle <- c(make_stan_data(sim, "oracle", df = 1, a_grid = a_grid), car)
  d_naive <- c(make_stan_data(sim, "naive", df = 1, a_grid = a_grid), car)
  d_known <- c(make_stan_data_latent_D(sim, bZ1 = 1, bZ2 = 0.5, sigma_Z = 1, a_grid = a_grid), car)
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

  # The two methods with a loose sigma_Z prior converge slowly, so they get more iterations
  fits <- list(
    "oracle" = fit_it(exposure_mod, d_oracle, s),
    "naive" = fit_it(exposure_mod, d_naive, s),
    "latent, known errors" = fit_it(known_mod, d_known, s, init_known),
    "reduced, sigma_Z fixed at 1" = fit_it(
      red_mod, c(d_red, list(sigma_Z_mean = 1, sigma_Z_sd = 0.001)), s, init_red
    ),
    "reduced, sigma_Z ~ normal(1, 0.3)" = fit_it(
      red_mod, c(d_red, list(sigma_Z_mean = 1, sigma_Z_sd = 0.3)), s, init_red, iter = 16000
    ),
    "reduced, sigma_Z wild guess" = fit_it(
      red_mod, c(d_red, list(sigma_Z_mean = 0, sigma_Z_sd = 5)), s, init_red, iter = 16000
    ),
    "reduced, tight priors on everything" = fit_it(tight_mod, d_red, s, init_red)
  )

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
    }))
  )
}


# Run, saving after every replicate
res <- list()
t0 <- Sys.time()
for (r in first:last) {
  res[[length(res) + 1]] <- run_one(r)
  saveRDS(
    list(
      summary = do.call(rbind, lapply(res, `[[`, "summary")),
      diag = do.call(rbind, lapply(res, `[[`, "diag"))
    ),
    out_file
  )
  cat("replicate", r, "done,", round(difftime(Sys.time(), t0, units = "mins"), 1), "min elapsed\n")
}
