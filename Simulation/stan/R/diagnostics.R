# Convergence diagnostics for one stanfit: R-hat (max and which parameter) and Geweke z-scores, per chain, on the scalar/low-dimensional
# parameters plus the estimand theta. Geweke compares the first 10% and last 50% of each chain (coda default); |z| > 1.96 flags non-stationarity.
library(coda)
key_pars <- c("beta0", "beta_A", "beta_X", "gamma0", "gamma_A", "sigma_A", "sigma_phi", "rho", "beta_Z2", "sigma_Z", "theta")
convergence_summary <- function(fit) {
  pars <- intersect(key_pars, fit@model_pars)
  mc <- As.mcmc.list(fit, pars = pars)
  z <- unlist(lapply(mc, function(ch) { g <- geweke.diag(ch, frac1 = 0.1, frac2 = 0.5)$z; g[is.finite(g)] }))
  sm <- summary(fit)$summary; rh <- sm[, "Rhat"]; rh_key <- rh[grepl(paste0("^(", paste(pars, collapse = "|"), ")(\\[|$)"), names(rh))]
  sp <- get_sampler_params(fit, inc_warmup = FALSE)
  list(div = sum(sapply(sp, function(m) sum(m[, "divergent__"]))),
       rhat_all = max(rh, na.rm = TRUE), rhat_key = max(rh_key, na.rm = TRUE), worst_param = names(rh_key)[which.max(rh_key)],
       geweke_max_abs_z = max(abs(z)), geweke_frac_flag = mean(abs(z) > 1.96), n_z = length(z))
}
