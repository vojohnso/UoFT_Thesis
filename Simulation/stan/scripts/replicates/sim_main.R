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
# Optional arguments: first and last replicate (default 1 100), and "tight" to run only the tight-priors method
args <- commandArgs(trailingOnly = TRUE)
first <- if (length(args) >= 1) as.integer(args[1]) else 1L
last <- if (length(args) >= 2) as.integer(args[2]) else 100L
tight_only <- length(args) >= 3 && args[3] == "tight" # run only the tight-priors method

# Iterations and warm-up per method, from the tuning run (stan/scripts/checks/tune_iterations.R).
# The two methods with a loose sigma_Z prior converge slowly and need the longer run.
iter_default <- c(iter = 8000, warmup = 4000)
iter_long <- c(iter = 16000, warmup = 8000)
iter_for <- function(method) {
  if (method %in% c("reduced, sigma_Z ~ normal(1, 0.3)", "reduced, sigma_Z wild guess")) iter_long else iter_default
}

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
# Fit one Stan model; keep the posterior draws of theta and the convergence summaries
fit_it <- function(mod, dat, seed, init = "random", method) {
  it <- iter_for(method)
  fit <- sampling(
    mod, data = dat, chains = 4, iter = it[["iter"]], warmup = it[["warmup"]], cores = 4,
    seed = seed, refresh = 0, init = init, control = list(adapt_delta = 0.9)
  )
  list(
    theta = rstan::extract(fit, "theta")$theta,
    conv = convergence_summary(fit),
    iter = it[["iter"]], warmup = it[["warmup"]]
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

  fit_specs <- list(
    "oracle" = function() {
      d <- c(make_stan_data(sim, "oracle", df = 1, a_grid = a_grid), car)
      fit_it(exposure_mod, d, s, method = "oracle")
    },
    "naive" = function() {
      d <- c(make_stan_data(sim, "naive", df = 1, a_grid = a_grid), car)
      fit_it(exposure_mod, d, s, method = "naive")
    },
    "latent, known errors" = function() {
      fit_it(known_mod, d_known, s, init_known, method = "latent, known errors")
    },
    "reduced, sigma_Z fixed at 1" = function() {
      d <- c(d_red, list(sigma_Z_mean = 1, sigma_Z_sd = 0.001))
      fit_it(red_mod, d, s, init_red(1), method = "reduced, sigma_Z fixed at 1")
    },
    "reduced, sigma_Z ~ normal(1, 0.3)" = function() {
      d <- c(d_red, list(sigma_Z_mean = 1, sigma_Z_sd = 0.3))
      fit_it(red_mod, d, s, init_red(1), method = "reduced, sigma_Z ~ normal(1, 0.3)")
    },
    "reduced, sigma_Z wild guess" = function() {
      d <- c(d_red, list(sigma_Z_mean = 0, sigma_Z_sd = 5))
      fit_it(red_mod, d, s, init_red(1), method = "reduced, sigma_Z wild guess")
    },
    "reduced, tight priors on everything" = function() {
      fit_it(tight_mod, d_red, s, init_red(1), method = "reduced, tight priors on everything")
    }
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
    # Convergence per method: divergences, R-hat (all parameters and the key ones), Geweke
    diag = do.call(rbind, lapply(names(fits), function(m) {
      cv <- fits[[m]]$conv
      data.frame(
        rep = r, method = m, iter = fits[[m]]$iter, warmup = fits[[m]]$warmup,
        div = cv$div, max_rhat = cv$rhat_all, rhat_key = cv$rhat_key, worst_param = cv$worst_param,
        geweke_frac_flag = cv$geweke_frac_flag, geweke_max_abs_z = cv$geweke_max_abs_z
      )
    }))
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
      n_iter = rbind(default = iter_default, long = iter_long)
    ),
    out_file
  )

  cat(
    "replicate", r, "(", r - first + 1, "of", last - first + 1, ") done,",
    round(difftime(Sys.time(), t0, units = "mins"), 1), "min elapsed\n"
  )
}
