# Main study: replicates first..last (arguments, default 1..100) of the core methods.
#
# Setting: linear truth, 100 areas, spatial (CAR) effect in the outcome, true sigma_Z = 1.
# Replicate r always uses the r-th seed drawn from the master seed, so any sub-range
# reproduces the corresponding replicates of the full run.
#
# Seven fits per replicate, all on the same simulated dataset:
#   1. oracle                              (true exposure A, with CAR)
#   2. naive                               (surrogate Z used as if it were A, with CAR)
#   3. latent, known errors                (error model given as data, with CAR)
#   4. reduced, sigma_Z fixed at 1         (sigma_Z ~ normal(1, 0.001))
#   5. reduced, sigma_Z ~ normal(1, 0.3)   (informative prior around the truth)
#   6. reduced, sigma_Z wild guess         (sigma_Z ~ normal(0, 5), truncated at 0)
#   7. reduced, tight priors on everything (diagnostic: every parameter almost known)
#
# Each method was first tried on its own at 20 replicates (the other scripts in this folder);
# this script runs all of them on the same 100 replicates, and report_all_methods.R summarises them.
#
# Results are saved after every replicate to
#   output/stan/latent_D_car/main_reps_<first>_<last>.rds          (all methods)
#   output/stan/latent_D_car/main_reps_<first>_<last>_tight.rds    (only the tight-priors method)
# so partial results survive an interruption.
#
# Usage (from the Simulation folder):
#   Rscript stan/scripts/replicates/sim_main.R               # replicates 1-100, all methods
#   Rscript stan/scripts/replicates/sim_main.R 21 60         # replicates 21-60, all methods
#   Rscript stan/scripts/replicates/sim_main.R 21 100 tight  # replicates 21-100, tight-priors method only

library(rstan)
source("gen_data.R")
source("stan/R/prep_data.R")
source("stan/R/compile_model.R")


# Settings
args <- commandArgs(trailingOnly = TRUE)
first <- if (length(args) >= 1) as.integer(args[1]) else 1L
last <- if (length(args) >= 2) as.integer(args[2]) else 100L
tight_only <- length(args) >= 3 && args[3] == "tight" # run only the tight-priors method
n_iter <- 2000

base_seed <- 1006092577
n_areas <- 100
a_grid <- seq(5, 15, by = 1)

# One seed per replicate, drawn once from the master seed
set.seed(base_seed)
rep_seeds <- sample.int(.Machine$integer.max, 100)



# Models and truth
exposure_mod <- compile_model("gcomp_known_exposure_car")
known_mod <- compile_model("gcomp_latent_exposure_D_known_car")
red_mod <- compile_model("gcomp_latent_exposure_D_reduced_sparse")
tight_mod <- compile_model("gcomp_latent_exposure_D_reduced_sparse_tight")

# Population truth of the ERF, the same in every replicate
pop_truth <- true_erf_pop(a_grid, "linear", spatial_in = "outcome", n_areas = n_areas)

# Where this run saves its results
out_file <- sprintf("output/stan/latent_D_car/main_reps_%d_%d%s.rds", first, last, if (tight_only) "_tight" else "")


# Helpers
# Fit one Stan model; keep the posterior draws of theta and two convergence summaries
fit_it <- function(mod, dat, seed, init = "random") {
  fit <- sampling(
    mod, data = dat, chains = 4, iter = n_iter, cores = 4, seed = seed,
    refresh = 0, init = init, control = list(adapt_delta = 0.9)
  )
  sampler_params <- get_sampler_params(fit, inc_warmup = FALSE)
  list(
    theta = rstan::extract(fit, "theta")$theta,
    div = sum(sapply(sampler_params, function(m) sum(m[, "divergent__"]))),
    rhat = max(summary(fit)$summary[, "Rhat"], na.rm = TRUE)
  )
}

# Run all six methods on replicate r
run_one <- function(r) {
  s <- rep_seeds[r]

  # Simulate one dataset
  sim <- data.generate(
    seed = s, n_areas = n_areas, out_form = "linear", spatial_in = "outcome",
    bZ1 = 1, bZ2 = 0.5, sigma_Z = 1
  )
  car <- car_structure(nrow(sim))

  # Stan data for the latent models
  d_known <- c(
    make_stan_data_latent_D(sim, bZ1 = 1, bZ2 = 0.5, sigma_Z = 1, a_grid = a_grid),
    car
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
  init_red <- function(sz) function() c(init_known(), list(beta_Z2 = 0, sigma_Z = sz))

  # The fits. Each is wrapped in a function so that a run can fit only some of them.
  fit_specs <- list(
    "oracle" = function() {
      fit_it(exposure_mod, c(make_stan_data(sim, "oracle", df = 1, a_grid = a_grid), car), s)
    },
    "naive" = function() {
      fit_it(exposure_mod, c(make_stan_data(sim, "naive", df = 1, a_grid = a_grid), car), s)
    },
    "latent, known errors" = function() fit_it(known_mod, d_known, s, init_known),
    "reduced, sigma_Z fixed at 1" = function() {
      fit_it(red_mod, c(d_red, list(sigma_Z_mean = 1, sigma_Z_sd = 0.001)), s, init_red(1))
    },
    "reduced, sigma_Z ~ normal(1, 0.3)" = function() {
      fit_it(red_mod, c(d_red, list(sigma_Z_mean = 1, sigma_Z_sd = 0.3)), s, init_red(1))
    },
    "reduced, sigma_Z wild guess" = function() {
      fit_it(red_mod, c(d_red, list(sigma_Z_mean = 0, sigma_Z_sd = 5)), s, init_red(1))
    },
    "reduced, tight priors on everything" = function() fit_it(tight_mod, d_red, s, init_red(1))
  )
  to_run <- if (tight_only) "reduced, tight priors on everything" else names(fit_specs)
  fits <- lapply(fit_specs[to_run], function(f) f())

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
    # Divergences and max R-hat per method
    diag = data.frame(
      rep = r, method = names(fits),
      div = sapply(fits, `[[`, "div"),
      max_rhat = sapply(fits, `[[`, "rhat")
    )
  )
}


# Run, saving after every replicate
res <- list()
t0 <- Sys.time()
for (r in first:last) {
  res[[as.character(r)]] <- run_one(r)

  saveRDS(
    list(
      summary = do.call(rbind, lapply(res, `[[`, "summary")),
      diag = do.call(rbind, lapply(res, `[[`, "diag")),
      n_iter = n_iter
    ),
    out_file
  )

  cat(
    "replicate", r, "(", r - first + 1, "of", last - first + 1, ") done,",
    round(difftime(Sys.time(), t0, units = "mins"), 1), "min elapsed\n"
  )
}
