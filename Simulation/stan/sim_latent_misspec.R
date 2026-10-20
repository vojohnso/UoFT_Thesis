# Misspecified exposure model: what happens to the latent model when the exposure is NOT linear in the confounders
# but the Stan model assumes it is (A ~ normal(gamma0 + X * gamma, sigma_A)).
# Each replicate is a pair of datasets with the same confounders, populations and noise, differing only in the
# exposure mean: linear (the latent model is exactly right) and nonlinear (exposure_nonlin = 0.5, misspecified).
# Fits: latent model on both, plus oracle and naive on the nonlinear dataset (the benchmarks).
# One master seed; the replicate seeds are drawn from it up front. Linear outcome scenario, no spatial term,
# error Z = A + N(0, 1) with sigma_Z known.
# Run from the Simulation/ folder:  Rscript stan/sim_latent_misspec.R [n_reps=30]
library(rstan)
rstan_options(auto_write = TRUE)
source("gen_data.R"); source("stan/prep_data.R")

args <- commandArgs(trailingOnly = TRUE)
n_reps <- if (length(args) >= 1) as.integer(args[1]) else 30
base_seed <- 1006092577
a_grid <- seq(5, 15, by = 1)
sigma_Z <- 1
nl <- 0.5
set.seed(base_seed)
rep_seeds <- sample.int(.Machine$integer.max, n_reps)    # one data seed per replicate, from the master seed

known_mod <- stan_model("stan/gcomp_known_exposure.stan")
latent_mod <- stan_model("stan/gcomp_latent_exposure.stan")
pop_truth <- true_erf_pop(a_grid, "linear")

fit_theta <- function(mod, dat, seed, init = "random") {
  fit <- sampling(mod, data = dat, chains = 4, iter = 2000, cores = 4, seed = seed, refresh = 0, init = init)
  list(fit = fit, theta = rstan::extract(fit, "theta")$theta)
}
fit_latent <- function(sim, seed) {
  dat <- make_stan_data_latent(sim, sigma_Z, a_grid)
  init <- function() list(beta0 = -3, beta_A = 0, beta_X = rep(0, 4), gamma0 = mean(sim$Zi), gamma = rep(0, 4),
                          sigma_A = 1, A = dat$Z + rnorm(dat$N, 0, 0.1))
  f <- fit_theta(latent_mod, dat, seed, init)
  f$A_rmse <- sqrt(mean((colMeans(rstan::extract(f$fit, "A")$A) - sim$Ai)^2))
  f$sigma_A <- mean(rstan::extract(f$fit, "sigma_A")$sigma_A)
  f
}

run_one <- function(r) {
  s <- rep_seeds[r]
  gen <- function(nonlin) data.generate(seed = s, out_form = "linear", rho_A = 0, exposure_nonlin = nonlin,
                                        bZ1 = 1, bZ2 = 0, sigma_Z = sigma_Z)
  sim_lin <- gen(0); sim_nl <- gen(nl)
  fits <- list(
    "latent, linear exposure (correct)" = fit_latent(sim_lin, s),
    "latent, nonlinear exposure (misspecified)" = fit_latent(sim_nl, s),
    "oracle (true A)" = fit_theta(known_mod, make_stan_data(sim_nl, "oracle", df = 1, a_grid = a_grid), s),
    "naive (Z as A)" = fit_theta(known_mod, make_stan_data(sim_nl, "naive", df = 1, a_grid = a_grid), s))
  do.call(rbind, lapply(names(fits), function(m) {
    th <- fits[[m]]$theta
    data.frame(rep = r, method = m, a = a_grid, true = pop_truth, est = colMeans(th),
               lo = apply(th, 2, quantile, 0.025), hi = apply(th, 2, quantile, 0.975),
               A_rmse = if (is.null(fits[[m]]$A_rmse)) NA else fits[[m]]$A_rmse,
               sigma_A = if (is.null(fits[[m]]$sigma_A)) NA else fits[[m]]$sigma_A,
               Z_rmse = sqrt(mean((sim_nl$Zi - sim_nl$Ai)^2)))
  }))
}

res <- vector("list", n_reps)
for (r in seq_len(n_reps)) { res[[r]] <- run_one(r); cat("replicate", r, "of", n_reps, "done\n") }
res <- do.call(rbind, res)
dir.create("output/stan", recursive = TRUE, showWarnings = FALSE)
saveRDS(res, "output/stan/latent_misspec.rds")

summ <- do.call(rbind, lapply(split(res, list(res$method, res$a)), function(d)
  data.frame(method = d$method[1], a = d$a[1], rel_bias_pct = 100 * mean(d$est - d$true) / d$true[1],
             rmse = sqrt(mean((d$est - d$true)^2)), coverage = mean(d$lo <= d$true & d$true <= d$hi))))
tab <- aggregate(cbind(abs_rel_bias_pct = abs(rel_bias_pct), rmse, coverage) ~ method, summ, mean)
extra <- aggregate(cbind(A_rmse, sigma_A) ~ method, res[!is.na(res$A_rmse), ], mean)
cat("\nAveraged over the grid, ", n_reps, " replicates (truth sigma_A about 1.5 for the linear data; Z alone: RMSE of A about ",
    round(mean(res$Z_rmse), 2), "):\n", sep = "")
print(merge(tab, extra, all.x = TRUE), digits = 3)
write.csv(summ, "output/stan/latent_misspec_summary.csv", row.names = FALSE)
