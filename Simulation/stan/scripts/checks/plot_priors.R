# Prior distributions of the reduced latent model (gcomp_latent_exposure_D_reduced_sparse.stan)
# against the true values used in the simulation, to show where each prior supports the truth.
#
# Identified parameters have non-informative priors. sigma_Z is not identified, so it is shown
# under its two density cases (informative and wild guess); the third case, fixing it at the
# truth, is a point mass and is not drawn here.
#
# sigma_A has the lower bound 0.5, and the half-Cauchy priors are truncated at 0, so their
# densities are renormalised (doubled for the half-Cauchy).
#
# Writes output/stan/checks/prior_vs_truth.png and prior_vs_truth.csv

library(ggplot2)
source("gen_data.R")

# One simulated dataset, only to read off the true sigma_phi
sim <- data.generate(
  seed = 1006092577, n_areas = 100, out_form = "linear", spatial_in = "outcome",
  bZ1 = 1, bZ2 = 0.5, sigma_Z = 1
)
sigma_phi_true <- sqrt(sim$sigma2_phi[1])


# Priors
# Truncate a prior at `lo` (density d, cdf p) and renormalise it
trunc_pri <- function(d, p, lo) {
  list(
    d = function(x) ifelse(x >= lo, d(x) / (1 - p(lo)), 0),
    p = function(x) pmax(p(x) - p(lo), 0) / (1 - p(lo))
  )
}
sA <- trunc_pri(function(x) dcauchy(x), function(x) pcauchy(x), 0.5) # sigma_A >= 0.5
sP <- trunc_pri(function(x) dcauchy(x), function(x) pcauchy(x), 0) # sigma_phi >= 0
sZ <- trunc_pri(function(x) dnorm(x, 1, 0.3), function(x) pnorm(x, 1, 0.3), 0) # sigma_Z, informative
sZw <- trunc_pri(function(x) dnorm(x, 0, 5), function(x) pnorm(x, 0, 5), 0) # sigma_Z, wild guess

# Each prior: name, label, density d, cdf p, plotting range, and the true value(s)
priors <- list(
  list(
    name = "beta0", label = "normal(0, 10)",
    d = function(x) dnorm(x, 0, 10), p = function(x) pnorm(x, 0, 10),
    range = c(-30, 30), truth = -4
  ),
  list(
    name = "beta_A", label = "normal(0, 10)",
    d = function(x) dnorm(x, 0, 10), p = function(x) pnorm(x, 0, 10),
    range = c(-30, 30), truth = 0.1
  ),
  list(
    name = "beta_X", label = "normal(0, 10)",
    d = function(x) dnorm(x, 0, 10), p = function(x) pnorm(x, 0, 10),
    range = c(-30, 30), truth = c(-0.5, -0.25, 0.25, 0.5)
  ),
  list(
    name = "gamma0", label = "normal(mean(Z), 10)",
    d = function(x) dnorm(x, 10, 10), p = function(x) pnorm(x, 10, 10),
    range = c(-20, 40), truth = 10
  ),
  list(
    name = "gamma_A", label = "normal(0, 10)",
    d = function(x) dnorm(x, 0, 10), p = function(x) pnorm(x, 0, 10),
    range = c(-30, 30), truth = c(-0.75, 0.75, -0.75, 0.75)
  ),
  list(
    name = "sigma_A", label = "cauchy(0, 1), >= 0.5",
    d = sA$d, p = sA$p,
    range = c(0, 6), truth = 1.5
  ),
  list(
    name = "sigma_phi", label = "cauchy(0, 1), >= 0",
    d = sP$d, p = sP$p,
    range = c(0, 3), truth = sigma_phi_true
  ),
  list(
    name = "rho", label = "beta(2, 1)",
    d = function(x) dbeta(x, 2, 1), p = function(x) pbeta(x, 2, 1),
    range = c(0, 1), truth = 0.9
  ),
  list(
    name = "beta_Z2", label = "normal(0, 10)",
    d = function(x) dnorm(x, 0, 10), p = function(x) pnorm(x, 0, 10),
    range = c(-30, 30), truth = 0.5
  ),
  list(
    name = "sigma_Z", label = "informative: normal(1, 0.3), >= 0",
    d = sZ$d, p = sZ$p,
    range = c(0, 3), truth = 1
  ),
  list(
    name = "sigma_Z", label = "wild guess: normal(0, 5), >= 0",
    d = sZw$d, p = sZw$p,
    range = c(0, 10), truth = 1
  )
)


# Plot: prior density with the true value as a red dashed line
# Prior density on a grid, one facet per prior
curves <- do.call(rbind, lapply(priors, function(pr) {
  x <- seq(pr$range[1], pr$range[2], length.out = 400)
  data.frame(
    param = paste0(pr$name, "\n", pr$label),
    prior = pr$label,
    x = x,
    density = pr$d(x)
  )
}))

# True values, one or more per facet
truths <- unique(do.call(rbind, lapply(priors, function(pr) {
  data.frame(param = paste0(pr$name, "\n", pr$label), truth = pr$truth)
})))

curves$param <- factor(curves$param, levels = unique(curves$param))
truths$param <- factor(truths$param, levels = levels(curves$param))

pl <- ggplot(curves, aes(x, density)) +
  geom_line() +
  geom_vline(data = truths, aes(xintercept = truth), colour = "red", linetype = "dashed", inherit.aes = FALSE) +
  facet_wrap(~param, scales = "free", ncol = 3) +
  labs(
    x = NULL, y = "prior density",
    title = "Priors of the reduced latent model against the true values",
    subtitle = "Red dashed = true value used in the simulation"
  ) +
  theme_bw(base_size = 11)

dir.create("output/stan/checks", recursive = TRUE, showWarnings = FALSE)
ggsave("output/stan/checks/prior_vs_truth.png", pl, width = 11, height = 10, dpi = 120)


# Table: where the truth sits under each prior
# The prior cdf at the truth (near 0 or 1 means the truth is far out in a tail)
tab <- do.call(rbind, lapply(priors, function(pr) {
  data.frame(
    parameter = pr$name,
    prior = pr$label,
    truth = pr$truth,
    prior_cdf_at_truth = round(pr$p(pr$truth), 4)
  )
}))
tab$inside_central_95 <- tab$prior_cdf_at_truth >= 0.025 & tab$prior_cdf_at_truth <= 0.975

write.csv(tab, "output/stan/checks/prior_vs_truth.csv", row.names = FALSE)
print(tab, row.names = FALSE)
