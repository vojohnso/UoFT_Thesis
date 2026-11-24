# Stage 3a: is naive (surrogate Z used as if it were A) biased under Berkson error as it is under classical error?
# Oracle (true A) and naive (Z) fits on replicates first..last (default 1..10), for four settings:
#   truth: linear (df = 1 basis) or quadratic (natural spline, df = 3), both with the CAR effect in the outcome
#   error: classical (the default of data.generate) or Berkson (A = Z + N, sigma_N = 1)
# Within a replicate the true exposure A and the outcomes Y are identical across the two error types
# (same seed); only the surrogate Z differs. The Berkson error SD matches the classical sigma_Z = 1
# (the classical Z also carries the 0.5*D shift, so its total error variance is 1.25).

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

# Model (known exposure with CAR; oracle and naive differ only in the exposure data)
exposure_mod <- compile_model("gcomp_known_exposure_car")

# The four settings: truth, spline df used for the exposure effect, error type
settings <- data.frame(
  truth = c("linear", "linear", "quadratic", "quadratic"),
  df = c(1, 1, 3, 3),
  error_type = c("classical", "berkson", "classical", "berkson")
)

# Population truth of the ERF, the same in every replicate of a truth
pop_truth <- lapply(c(linear = "linear", quadratic = "quadratic"), function(f) {
  true_erf_pop(a_grid, f, spatial_in = "outcome", n_areas = n_areas)
})

dir.create("output/stan/latent_D_car", recursive = TRUE, showWarnings = FALSE)
out_file <- sprintf("output/stan/latent_D_car/berkson_3a_reps_%d_%d.rds", first, last)

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

# Oracle and naive for every setting on replicate r
run_one <- function(r) {
  s <- rep_seeds[r]
  summ <- list()
  dgn <- list()
  for (i in seq_len(nrow(settings))) {
    st <- settings[i, ]
    sim <- data.generate(
      seed = s, n_areas = n_areas, out_form = st$truth, spatial_in = "outcome",
      bZ1 = 1, bZ2 = 0.5, sigma_Z = 1, error_type = st$error_type, sigma_N = 1
    )
    car <- car_structure(nrow(sim))
    for (m in c("oracle", "naive")) {
      d <- c(make_stan_data(sim, m, df = st$df, a_grid = a_grid), car)
      f <- fit_it(exposure_mod, d, s)
      th <- f$theta
      summ[[length(summ) + 1]] <- data.frame(
        rep = r, truth = st$truth, error_type = st$error_type, method = m, a = a_grid,
        true = pop_truth[[st$truth]], true_sample = true_erf(a_grid, sim),
        est = colMeans(th),
        lo = apply(th, 2, quantile, 0.025),
        hi = apply(th, 2, quantile, 0.975)
      )
      dgn[[length(dgn) + 1]] <- data.frame(
        rep = r, truth = st$truth, error_type = st$error_type, method = m, f$conv
      )
    }
  }
  list(summary = do.call(rbind, summ), diag = do.call(rbind, dgn))
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
