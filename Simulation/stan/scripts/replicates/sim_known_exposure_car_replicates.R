# Oracle (true exposure A) and naive (surrogate Z) outcome models, both WITH the CAR spatial
# effect in the outcome (gcomp_known_exposure_car.stan).
#
# Same 20 seeds, data and truth as sim_latent_D_reduced_replicates.R, so the rows can be compared
# with the oracle/naive fits that have no spatial term.

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

mod <- compile_model("gcomp_known_exposure_car")

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

# Fit oracle and naive on replicate r
run_one <- function(r) {
  s <- rep_seeds[r]
  sim <- data.generate(
    seed = s, n_areas = n_areas, out_form = "linear", spatial_in = "outcome",
    bZ1 = 1, bZ2 = 0.5, sigma_Z = 1
  )
  car <- car_structure(nrow(sim))

  fits <- list(
    "oracle + CAR" = fit_it(mod, c(make_stan_data(sim, "oracle", df = 1, a_grid = a_grid), car), s),
    "naive + CAR" = fit_it(mod, c(make_stan_data(sim, "naive", df = 1, a_grid = a_grid), car), s)
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


# Run and save
t0 <- Sys.time()
res <- list()
for (r in seq_len(n_reps)) {
  res[[r]] <- run_one(r)
  cat("replicate", r, "of", n_reps, "done,", round(difftime(Sys.time(), t0, units = "mins"), 1), "min elapsed\n")
}

out <- list(
  summary = do.call(rbind, lapply(res, `[[`, "summary")),
  diag = do.call(rbind, lapply(res, `[[`, "diag"))
)
dir.create("output/stan/latent_D_car", recursive = TRUE, showWarnings = FALSE)
saveRDS(out, "output/stan/latent_D_car/known_exposure_car_replicates.rds")
