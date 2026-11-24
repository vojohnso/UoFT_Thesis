# Does the spline df explain the underestimated peak? Oracle fit (true exposure A) at df = 2, 3, 4, 5
# on the same datasets as the spline study (replicates first..last, default 1..10), quadratic truth.

library(rstan)
source("gen_data.R")
source("stan/R/prep_data.R")
source("stan/R/compile_model.R")
source("stan/R/diagnostics.R")

# Settings
# Optional arguments: first and last replicate (default 1 10)
args <- commandArgs(trailingOnly = TRUE)
first <- if (length(args) >= 1) as.integer(args[1]) else 1L
last <- if (length(args) >= 2) as.integer(args[2]) else 10L

base_seed <- 1006092577
n_areas <- 100
a_grid <- seq(5, 15, by = 1)

# One seed per replicate, drawn once from the master seed
set.seed(base_seed)
rep_seeds <- sample.int(.Machine$integer.max, 100)

# Models
exposure_mod <- compile_model("gcomp_known_exposure_car")
dfs <- c(2, 3, 4, 5)

# Population truth of the ERF, the same in every replicate
pop_truth <- true_erf_pop(a_grid, "quadratic", spatial_in = "outcome", n_areas = n_areas)

dir.create("output/stan/latent_D_car", recursive = TRUE, showWarnings = FALSE)
out_file <- sprintf("output/stan/latent_D_car/spline_df_oracle_reps_%d_%d.rds", first, last)

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

# Run all five methods on replicate r
run_one <- function(r) {
  s <- rep_seeds[r]
  sim <- data.generate(
    seed = s, n_areas = n_areas, out_form = "quadratic", spatial_in = "outcome",
    bZ1 = 1, bZ2 = 0.5, sigma_Z = 1
  )
  car <- car_structure(nrow(sim))

  # Oracle fit (true exposure A) at each spline df
  fits <- list()
  for (k in dfs) {
    d <- c(make_stan_data(sim, "oracle", df = k, a_grid = a_grid), car)
    fits[[paste0("oracle spline df=", k)]] <- fit_it(exposure_mod, d, s)
  }

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
