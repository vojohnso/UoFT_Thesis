# Berkson study, linear truth: replicates first..last of the Berkson latent method against oracle and naive.
#
# Setting: linear truth, 100 areas, spatial (CAR) effect in the outcome, Berkson error A = Z + N, N ~ normal(0, sigma_N^2).
#
# Five fits per replicate, all on the same simulated dataset:
#   1. oracle                              (true exposure A, with CAR)
#   2. naive                               (surrogate Z used as if it were A, with CAR)
#   3. Berkson, sigma_N fixed at truth     (sigma_N ~ normal(true, 0.001))
#   4. Berkson, sigma_N ~ normal(truth, 0.3)  (informative prior around the truth)
#   5. Berkson, sigma_N wild guess         (sigma_N ~ normal(0, 5), truncated at 0)
#
# Arguments: first, last replicate (default 1 10) and sigma_N, the true Berkson error SD (default 1; must be below
# sigma_A = 1.5). The sigma_N in the file name keeps different error sizes apart.

library(rstan)
source("gen_data.R")
source("stan/R/prep_data.R")
source("stan/R/compile_model.R")
source("stan/R/diagnostics.R")

# Settings
args <- commandArgs(trailingOnly = TRUE)
first <- if (length(args) >= 1) as.integer(args[1]) else 1L
last <- if (length(args) >= 2) as.integer(args[2]) else 10L
sigma_N <- if (length(args) >= 3) as.numeric(args[3]) else 1

base_seed <- 1006092577
n_areas <- 100
a_grid <- seq(5, 15, by = 1)

# One seed per replicate, drawn once from the master seed
set.seed(base_seed)
rep_seeds <- sample.int(.Machine$integer.max, 100)

# Models
exposure_mod <- compile_model("gcomp_known_exposure_car")
berkson_mod <- compile_model("gcomp_latent_exposure_berkson_sparse")

# Population truth of the ERF, the same in every replicate
pop_truth <- true_erf_pop(a_grid, "linear", spatial_in = "outcome", n_areas = n_areas)

dir.create("output/stan/latent_D_car", recursive = TRUE, showWarnings = FALSE)
out_file <- sprintf("output/stan/latent_D_car/berkson_latent_sN%s_reps_%d_%d.rds", sigma_N, first, last)

# Fit one Stan model (warm-up is half the iterations); keep theta and the convergence summary
fit_it <- function(mod, dat, seed, init = "random", iter = 8000) {
  fit <- sampling(
    mod, data = dat, chains = 4, iter = iter, warmup = iter / 2, cores = 4,
    seed = seed, refresh = 0, init = init, control = list(adapt_delta = 0.9)
  )
  # Posterior mean and sd of beta_A and (Berkson fits) sigma_N, to see whether the data inform sigma_N
  ps <- rstan::extract(fit, c("beta_A", if ("sigma_N" %in% fit@model_pars) "sigma_N"))
  post <- data.frame(beta_A_mean = mean(ps$beta_A), beta_A_sd = sd(ps$beta_A))
  if (!is.null(ps$sigma_N)) post <- cbind(post, sigma_N_mean = mean(ps$sigma_N), sigma_N_sd = sd(ps$sigma_N))
  list(
    theta = rstan::extract(fit, "theta")$theta,
    conv = data.frame(iter = iter, convergence_summary(fit), post)
  )
}

# Run all five methods on replicate r
run_one <- function(r) {
  s <- rep_seeds[r]
  sim <- data.generate(
    seed = s, n_areas = n_areas, out_form = "linear", spatial_in = "outcome",
    error_type = "berkson", sigma_N = sigma_N
  )
  car <- car_structure(nrow(sim))

  # Stan data
  d_oracle <- c(make_stan_data(sim, "oracle", df = 1, a_grid = a_grid), car)
  d_naive <- c(make_stan_data(sim, "naive", df = 1, a_grid = a_grid), car)
  d_berk <- make_stan_data_berkson(sim, a_grid)

  # Starting values: the noise starts at zero (A = Z)
  init_berk <- function() {
    list(
      beta0 = -3, beta_A = 0, beta_X = rep(0, 4), phi = rep(0, nrow(sim)),
      sigma_phi = 1, rho = 0.5, sigma_N = 1, A_raw = rep(0, nrow(sim))
    )
  }

  # The two methods with a loose sigma_N prior may converge slowly, so they get more iterations
  fits <- list(
    "oracle" = fit_it(exposure_mod, d_oracle, s),
    "naive" = fit_it(exposure_mod, d_naive, s),
    "Berkson, sigma_N fixed at truth" = fit_it(
      berkson_mod, c(d_berk, list(sigma_N_mean = sigma_N, sigma_N_sd = 0.001)), s, init_berk
    ),
    "Berkson, sigma_N ~ normal(truth, 0.3)" = fit_it(
      berkson_mod, c(d_berk, list(sigma_N_mean = sigma_N, sigma_N_sd = 0.3)), s, init_berk, iter = 16000
    ),
    "Berkson, sigma_N wild guess" = fit_it(
      berkson_mod, c(d_berk, list(sigma_N_mean = 0, sigma_N_sd = 5)), s, init_berk, iter = 16000
    )
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
      cv <- fits[[m]]$conv
      for (v in c("sigma_N_mean", "sigma_N_sd")) if (is.null(cv[[v]])) cv[[v]] <- NA
      data.frame(rep = r, method = m, cv)
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
