# Replicates for the latent exposure model: Oracle, Naive and Latent (all Stan) on the linear scenario with
# no spatial effect and constant, known error variance (so the latent model is exactly right).
# Run from the Simulation/ folder:  Rscript stan/scripts/sim_latent_replicates.R [n_reps=50]
# One seed for the whole study: set.seed() once, then the replicates are a plain loop that keeps drawing from that
# stream (data, initial values, and the seed handed to each Stan fit). The 4 chains of each fit run in parallel.
# To reproduce replicate r on its own you have to rerun the loop from the start.
library(rstan)
rstan_options(auto_write = TRUE)

source("gen_data.R")
source("stan/R/prep_data.R"); source("stan/R/compile_model.R")

args <- commandArgs(trailingOnly = TRUE)
n_reps <- if (length(args) >= 1) as.integer(args[1]) else 50
base_seed <- 1006092577
a_grid <- seq(5, 15, by = 1)
sigma_Z <- 1                  # error SD, also the known value passed to Stan
methods <- c("oracle", "naive", "latent")

# Compile once, before forking
known_mod <- compile_model("gcomp_known_exposure")
latent_mod <- compile_model("gcomp_latent_exposure")
pop_truth <- true_erf_pop(a_grid, "linear")

fit_theta <- function(mod, dat, seed, init = "random") {
  fit <- sampling(mod, data = dat, chains = 4, iter = 2000, cores = 4, seed = seed, refresh = 0, init = init)
  sp <- get_sampler_params(fit, inc_warmup = FALSE)
  list(theta = rstan::extract(fit, "theta")$theta,
       div = sum(sapply(sp, function(m) sum(m[, "divergent__"]))),
       rhat = max(summary(fit)$summary[, "Rhat"], na.rm = TRUE),
       fit = fit)
}

run_one <- function(rep_id) {
  sim <- data.generate(seed = NULL, out_form = "linear", rho_A = 0, bZ1 = 1, bZ2 = 0, sigma_Z = sigma_Z)
  seed <- sample.int(.Machine$integer.max, 1)   # Stan seed for this replicate, drawn from the same stream
  dat_lat <- make_stan_data_latent(sim, sigma_Z, a_grid)
  init_latent <- function() list(beta0 = -3, beta_A = 0, beta_X = rep(0, 4), gamma0 = mean(sim$Zi),
                                 gamma = rep(0, 4), sigma_A = 1, A = dat_lat$Z)
  fits <- list(
    oracle = fit_theta(known_mod, make_stan_data(sim, "oracle", df = 1, a_grid = a_grid), seed),
    naive = fit_theta(known_mod, make_stan_data(sim, "naive", df = 1, a_grid = a_grid), seed),
    latent = fit_theta(latent_mod, dat_lat, seed, init = init_latent))
  true_sample <- true_erf(a_grid, sim)
  summary <- do.call(rbind, lapply(methods, function(m) {
    th <- fits[[m]]$theta
    data.frame(rep = rep_id, method = m, a = a_grid, true = pop_truth, true_sample = true_sample,
               est = colMeans(th), lo = apply(th, 2, quantile, 0.025), hi = apply(th, 2, quantile, 0.975))
  }))
  pars <- rstan::extract(fits$latent$fit, c("beta_A", "sigma_A"))
  list(summary = summary,
       diag = data.frame(rep = rep_id, method = methods, div = sapply(methods, function(m) fits[[m]]$div),
                         max_rhat = sapply(methods, function(m) fits[[m]]$rhat)),
       latent_pars = data.frame(rep = rep_id, beta_A = mean(pars$beta_A), beta_A_lo = quantile(pars$beta_A, .025),
                                beta_A_hi = quantile(pars$beta_A, .975), sigma_A = mean(pars$sigma_A),
                                sigma_A_lo = quantile(pars$sigma_A, .025), sigma_A_hi = quantile(pars$sigma_A, .975)))
}

t0 <- Sys.time()
cat(n_reps, "replicates\n")
set.seed(base_seed)
res <- vector("list", n_reps)
for (r in seq_len(n_reps)) {
  res[[r]] <- run_one(r)
  cat("replicate", r, "of", n_reps, "done\n")
}
ok <- rep(TRUE, n_reps)
out <- list(summary = do.call(rbind, lapply(res, `[[`, "summary")),
            diag = do.call(rbind, lapply(res, `[[`, "diag")),
            latent_pars = do.call(rbind, lapply(res, `[[`, "latent_pars")),
            rep = which(ok), truth_beta_A = 0.10, truth_sigma_A = 1.5)
dir.create("output/stan", recursive = TRUE, showWarnings = FALSE)
saveRDS(out, "output/stan/latent_replicates.rds")
cat("done in", round(difftime(Sys.time(), t0, units = "mins"), 1), "min\n")
