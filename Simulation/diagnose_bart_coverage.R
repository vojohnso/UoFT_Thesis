# Why does the BART oracle undercover in the quadratic scenario (mostly at low exposure)?
# Oracle fits only (true A), quadratic scenario, population truth. Variants of the fit:
#   base    : current bart_gcomp (log rate, homoscedastic Gaussian noise), 400 areas
#   weights : same, but gbart weights w = max(Y, 0.5), since Var(log(Y/N)) is about 1/Y for a Poisson count
#   n1600   : base fit on 1600 areas (40 x 40 grid), to see whether the problem is sparse data
# Run from Simulation/:  Rscript diagnose_bart_coverage.R [n_reps=100] [n_cores]
library(parallel)
source("gen_data.R")
source("estimators.R")

args <- commandArgs(trailingOnly = TRUE)
n_reps <- if (length(args) >= 1) as.integer(args[1]) else 100
n_cores <- if (length(args) >= 2) as.integer(args[2]) else max(1, detectCores() - 1)
base_seed <- 1006092577
a_grid <- seq(5, 15, by = 1)
pop_truth <- true_erf_pop(a_grid, "quadratic")

# bart_gcomp from estimators.R, oracle only, with an optional weight vector for the noise variance
fit_variant <- function(sim, w = NULL) {
  X <- as.matrix(cbind(exp_c = sim$Ai - mean(sim$Ai), sim[, c("Xi1", "Xi2", "Xi3", "Xi4")]))
  y <- log(sim$Yi / sim$Ni)
  y[!is.finite(y)] <- log(0.5 / sim$Ni[!is.finite(y)])
  n <- nrow(sim)
  Xt <- do.call(rbind, lapply(a_grid - mean(sim$Ai), function(a) { Z <- X; Z[, "exp_c"] <- a; Z }))
  fit <- if (is.null(w)) gbart(X, y, x.test = Xt, type = "wbart", ndpost = 1000, nskip = 500, printevery = 0)
         else gbart(X, y, x.test = Xt, type = "wbart", w = w, ndpost = 1000, nskip = 500, printevery = 0)
  bb <- rdirichlet(1000, rep(1, n))
  th <- sapply(seq_along(a_grid), function(g) rowSums(bb * exp(fit$yhat.test[, ((g - 1) * n + 1):(g * n)])))
  data.frame(a = a_grid, est = colMeans(th), lo = apply(th, 2, quantile, 0.025), hi = apply(th, 2, quantile, 0.975))
}

run_one <- function(r) {
  sim <- data.generate(seed = base_seed + r, out_form = "quadratic")
  big <- data.generate(seed = base_seed + r, n_areas = 1600, out_form = "quadratic")
  invisible(capture.output({
    out <- list(base = fit_variant(sim), weights = fit_variant(sim, w = pmax(sim$Yi, 0.5)), n1600 = fit_variant(big))
  }))
  do.call(rbind, lapply(names(out), function(v) data.frame(rep = r, variant = v, out[[v]], true = pop_truth)))
}

res <- mclapply(seq_len(n_reps), run_one, mc.cores = n_cores)
res <- do.call(rbind, res[!sapply(res, inherits, "try-error")])
dir.create("output/diagnose_bart", showWarnings = FALSE)
saveRDS(res, "output/diagnose_bart/coverage_variants.rds")
cat("replicates used:", length(unique(res$rep)), "\n")
tab <- do.call(rbind, lapply(split(res, list(res$variant, res$a)), function(d)
  data.frame(variant = d$variant[1], a = d$a[1], coverage = mean(d$lo <= d$true & d$true <= d$hi),
             rel_bias_pct = 100 * mean(d$est - d$true) / d$true[1], rmse = sqrt(mean((d$est - d$true)^2)),
             half_width = mean(d$hi - d$lo) / 2)))
tab$ratio_halfwidth_to_rmse <- tab$half_width / tab$rmse   # about 2 for a calibrated 95% interval
write.csv(tab, "output/diagnose_bart/coverage_variants.csv", row.names = FALSE)
print(reshape(tab[, c("variant", "a", "coverage")], idvar = "a", timevar = "variant", direction = "wide"), digits = 3)
print(aggregate(cbind(coverage, ratio_halfwidth_to_rmse) ~ variant, tab, mean), digits = 3)
