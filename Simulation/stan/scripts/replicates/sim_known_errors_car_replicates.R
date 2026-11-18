# Known-errors latent model (error parameters given as data) WITH the CAR spatial effect in the
# outcome (gcomp_latent_exposure_D_known_car.stan).
#
# Same 20 seeds, data and truth as sim_latent_D_reduced_replicates.R, so the rows can be compared
# with the known-errors fit that has no spatial term.

library(rstan)
source("gen_data.R")
source("stan/R/prep_data.R")
source("stan/R/compile_model.R")


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

mod <- compile_model("gcomp_latent_exposure_D_known_car")
pop_truth <- true_erf_pop(a_grid, "linear", spatial_in = "outcome", n_areas = n_areas)


# Helpers
# Fit the model; keep the posterior draws of theta and two convergence summaries
fit_it <- function(dat, seed, init) {
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

# Fit the known-errors model on replicate r
run_one <- function(r) {
  s <- rep_seeds[r]
  sim <- data.generate(
    seed = s, n_areas = n_areas, out_form = "linear", spatial_in = "outcome",
    bZ1 = 1, bZ2 = 0.5, sigma_Z = 1
  )

  # Stan data: the error parameters are given, plus the neighbour structure for the CAR term
  d <- c(
    make_stan_data_latent_D(sim, bZ1 = 1, bZ2 = 0.5, sigma_Z = 1, a_grid = a_grid),
    car_structure(nrow(sim))
  )

  # Starting values: A starts at the observed surrogate Z
  init <- function() {
    list(
      beta0 = -3, beta_A = 0, beta_X = rep(0, 4),
      gamma0 = mean(sim$Zi), gamma_A = rep(0, 4), sigma_A = 1,
      A = sim$Zi,
      phi = rep(0, nrow(sim)), sigma_phi = 1, rho = 0.5
    )
  }

  fits <- list("latent, known errors + CAR" = fit_it(d, s, init))

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
    )
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
  n_iter = n_iter
)
saveRDS(out, "output/stan/latent_D_car/known_errors_car_replicates.rds")


# Summary printed to the console
raw <- out$summary

# At each exposure value: relative bias (%), RMSE and coverage over the replicates
by_a <- do.call(rbind, lapply(split(raw, list(raw$method, raw$a), drop = TRUE), function(d) {
  data.frame(
    method = d$method[1],
    a = d$a[1],
    rel_bias_pct = 100 * mean(d$est - d$true) / d$true[1],
    rmse = sqrt(mean((d$est - d$true)^2)),
    coverage = mean(d$lo <= d$true & d$true <= d$hi)
  )
}))

cat("\nAveraged over the grid,", n_reps, "replicates, coverage against the population truth:\n")
print(aggregate(cbind(abs_rel_bias_pct = abs(rel_bias_pct), rmse, coverage) ~ method, by_a, mean), digits = 3)

cat("\nSampler diagnostics:\n")
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
