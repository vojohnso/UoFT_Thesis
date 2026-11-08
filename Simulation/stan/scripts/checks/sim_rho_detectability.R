# How well can the CAR dependence parameter rho be detected and estimated?
# Best case, no MCMC: simulate the exposure residual r = phi + noise (phi ~ CAR(rho, sigma_phi^2), the same construction and
# R2_phi = 0.2 as data.generate(), residual SD sigma_A = 1.5), treat the exposure mean as known, and compute the profile
# likelihood of rho (sigma_phi and sigma_A maximised out). Two settings:
#   "A observed"          the residual of the true exposure A is available
#   "A seen through Z"    only Z = A + N(0, 1) is available (sigma_Z known), as in our model
# Reports, over replicate datasets: the spread of the maximum-likelihood estimate of rho, the width of the 95% likelihood-ratio
# interval, how often that interval covers the truth, and how often it excludes rho = 0 (i.e. spatial dependence is detected).
# Run from the Simulation/ folder:  Rscript stan/scripts/checks/sim_rho_detectability.R [n_reps=100]
suppressMessages(source("gen_data.R"))

args <- commandArgs(trailingOnly = TRUE)
n_reps <- if (length(args) >= 1) as.integer(args[1]) else 100
rho_grid <- c(0, seq(0.1, 0.9, by = 0.1), 0.95, 0.99)
true_rhos <- c(0.3, 0.5, 0.7, 0.9, 0.99)
sigma_A <- 1.5; sigma_Z <- 1; var_X <- 2.25; R2_phi <- 0.2
set.seed(1006092577)

run_n <- function(n_areas) {
  W <- make_adj_matrix(n_areas); Dd <- rowSums(W)
  eig <- lapply(rho_grid, function(rho) eigen(diag(Dd) - rho * W, symmetric = TRUE))     # depends on rho, not on the data
  names(eig) <- as.character(rho_grid)
  prof_ll <- function(y, noise_var) {                                                    # profile log-likelihood over rho_grid
    sapply(seq_along(rho_grid), function(j) {
      e <- eig[[j]]; yy <- as.vector(crossprod(e$vectors, y))^2; lam <- e$values
      nll <- function(th) { v <- exp(th[1]) / lam + exp(th[2]) + noise_var; 0.5 * sum(log(v) + yy / v) }
      -optim(c(log(var(y) / 2), log(var(y) / 2)), nll, method = "BFGS")$value
    })
  }
  res <- list()
  for (rho0 in true_rhos) {
    s2 <- compute_sigma2_phi(R2_phi, var_X, sigma_A^2, W, rho0)
    e0 <- eigen(diag(Dd) - rho0 * W, symmetric = TRUE)
    for (setting in c("A observed", "A seen through Z")) {
      nv <- if (setting == "A observed") 0 else sigma_Z^2
      out <- t(replicate(n_reps, {
        phi <- as.vector(e0$vectors %*% (sqrt(s2 / e0$values) * rnorm(n_areas)))
        r <- phi + rnorm(n_areas, 0, sigma_A) + if (nv > 0) rnorm(n_areas, 0, sqrt(nv)) else 0
        ll <- prof_ll(r, nv); ok <- 2 * (max(ll) - ll) <= 3.84
        c(rho_hat = rho_grid[which.max(ll)], lo = min(rho_grid[ok]), hi = max(rho_grid[ok]),
          covers = as.numeric(rho0 >= min(rho_grid[ok]) & rho0 <= max(rho_grid[ok])),
          detects = as.numeric(2 * (max(ll) - ll[1]) > 3.84))
      }))
      res[[length(res) + 1]] <- data.frame(n_areas = n_areas, setting = setting, true_rho = rho0,
        median_rho_hat = median(out[, "rho_hat"]), q25 = quantile(out[, "rho_hat"], .25), q75 = quantile(out[, "rho_hat"], .75),
        median_interval = paste0("[", median(out[, "lo"]), ", ", median(out[, "hi"]), "]"),
        share_interval_covers_truth = mean(out[, "covers"]), share_excluding_rho0 = mean(out[, "detects"]))
    }
  }
  do.call(rbind, res)
}
tab <- rbind(run_n(100), run_n(400)); rownames(tab) <- NULL
dir.create("output/stan/checks", recursive = TRUE, showWarnings = FALSE)
write.csv(tab, "output/stan/checks/rho_detectability.csv", row.names = FALSE)
cat("\nProfile-likelihood detectability of rho (", n_reps, " datasets per row; rho grid ", paste(rho_grid, collapse = ", "), ")\n", sep = "")
print(tab, digits = 2, row.names = FALSE)
