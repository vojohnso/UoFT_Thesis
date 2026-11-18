# One aggregate table and one set of plots for all methods of the main study (sim_main.R):
# linear truth, 100 areas, spatial effect in the outcome
# Methods, all with the CAR spatial term include (1) oracle, (2) naive, (3) latent with known errors,
# (4) reduced model with sigma_Z fixed at 1, (5) reduced model with sigma_Z ~ normal(1, 0.3),
# (6) reduced model with a wild guess on sigma_Z ~ normal(0, 5) and (7) tight priors on everything

library(ggplot2)
d <- "output/stan/latent_D_car"

# Optional arguments: which main_reps files to use (comma-separated, default all) and the folder to save to
args <- commandArgs(trailingOnly = TRUE)
files <- if (length(args) >= 1) {
  file.path(d, strsplit(args[1], ",")[[1]])
} else {
  list.files(d, pattern = "^main_reps_.*\\.rds$", full.names = TRUE)
}
out_dir <- if (length(args) >= 2) args[2] else d
dir.create(file.path(out_dir, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "plots"), recursive = TRUE, showWarnings = FALSE)
if (length(files) == 0) stop("no main_reps_*.rds files in ", d, " (run sim_main.R first)")
runs <- lapply(files, readRDS)

raw <- do.call(rbind, lapply(runs, function(r) r$summary[, c("rep", "method", "a", "true", "est", "lo", "hi")]))
dg_cols <- c(
  "rep", "method", "iter", "div", "max_rhat", "rhat_key", "worst_param",
  "geweke_max_abs_z", "geweke_frac_flag"
)
dg <- do.call(rbind, lapply(runs, function(r) {
  d <- r$diag
  for (v in setdiff(dg_cols, names(d))) d[[v]] <- NA # older runs did not store these
  d[, dg_cols]
}))

# Method order used in the table and plots
lv <- c(
  "oracle", "naive", "latent, known errors", "reduced, sigma_Z fixed at 1",
  "reduced, tight priors on everything",
  "reduced, sigma_Z ~ normal(1, 0.3)", "reduced, sigma_Z wild guess"
)
stopifnot(all(raw$method %in% lv))
raw$method <- factor(raw$method, levels = lv)


# Metrics

# At each exposure value a: relative bias (%), RMSE and coverage, over the replicates
by_a <- do.call(rbind, lapply(split(raw, list(raw$method, raw$a), drop = TRUE), function(x) {
  data.frame(
    method = x$method[1],
    a = x$a[1],
    truth = x$true[1],
    mean_est = mean(x$est),
    rel_bias_pct = 100 * mean(x$est - x$true) / x$true[1],
    rmse = sqrt(mean((x$est - x$true)^2)),
    coverage = mean(x$lo <= x$true & x$true <= x$hi)
  )
}))

# Average over the grid, per method: |relative bias|, RMSE, coverage
tab <- aggregate(cbind(abs_rel_bias_pct = abs(rel_bias_pct), rmse, coverage) ~ method, by_a, mean)

# Convergence summary per method: replicates with divergences, worst and median of the per-fit max R-hat
dgt <- do.call(rbind, lapply(split(dg, dg$method), function(x) {
  data.frame(
    method = x$method[1],
    reps_with_divergence = sum(x$div > 0),
    worst_rhat = max(x$max_rhat),
    median_rhat = median(x$max_rhat),
    fits_rhat_key_above_1.05 = sum(x$rhat_key > 1.05, na.rm = TRUE),
    mean_geweke_frac_flag = mean(x$geweke_frac_flag, na.rm = TRUE)
  )
}))

# Number of replicates per method
nrep <- aggregate(rep ~ method, raw, function(x) length(unique(x)))
names(nrep)[2] <- "n_reps"

# Combine and add the Monte Carlo standard error of the coverage
out <- merge(merge(nrep, tab, by = "method"), dgt, by = "method")
out <- out[order(match(out$method, lv)), ]
out$mc_se_coverage <- sqrt(out$coverage * (1 - out$coverage) / out$n_reps)

dir.create(file.path(out_dir, "tables"), showWarnings = FALSE)
dir.create(file.path(out_dir, "plots"), showWarnings = FALSE)
write.csv(out, file.path(out_dir, "tables", "aggregate_all_methods.csv"), row.names = FALSE)
print(out, digits = 3, row.names = FALSE)

# Every fit (replicate x method) with its convergence diagnostics
dg <- dg[order(match(dg$method, lv), dg$rep), ]
write.csv(dg, file.path(out_dir, "tables", "per_fit_diagnostics.csv"), row.names = FALSE)


# Plots
by_a$method <- factor(by_a$method, levels = lv)
n <- paste(unique(range(out$n_reps)), collapse = "-") # replicate count for the titles

# Estimated curves: one line per replicate, with the mean of the estimates and the truth
p1 <- ggplot(raw, aes(a, est, group = rep)) +
  geom_line(alpha = 0.3, linewidth = 0.35) +
  geom_line(data = by_a, aes(y = mean_est, group = 1), colour = "blue", linetype = "dashed", linewidth = 0.8) +
  geom_line(data = by_a, aes(y = truth, group = 1), colour = "red", linewidth = 0.9) +
  facet_wrap(~method, ncol = 2) +
  coord_cartesian(ylim = c(0, 3500)) +
  labs(
    x = "exposure a", y = "expected counts per area",
    title = paste0("Estimated curves, all methods (", n, " replicates)"),
    subtitle = "Red = truth, blue dashed = mean of estimates, grey = replicates (y-axis cut at 3500)"
  ) +
  theme_bw(base_size = 10)

# Relative bias across the exposure grid
p2 <- ggplot(by_a, aes(a, rel_bias_pct, colour = method)) +
  geom_hline(yintercept = 0, colour = "grey50") +
  geom_line() +
  geom_point(size = 1) +
  labs(
    x = "exposure a", y = "relative bias (%)",
    title = paste0("Relative bias, all methods (", n, " replicates)"),
    colour = NULL
  ) +
  theme_bw(base_size = 10)

# Coverage across the exposure grid
p3 <- ggplot(by_a, aes(a, coverage, colour = method)) +
  geom_hline(yintercept = 0.95, linetype = "dashed", colour = "grey50") +
  geom_line() +
  geom_point(size = 1) +
  ylim(0, 1) +
  labs(
    x = "exposure a", y = "coverage of the 95% interval",
    title = paste0("Coverage, all methods (", n, " replicates, dashed = 0.95)"),
    colour = NULL
  ) +
  theme_bw(base_size = 10)

ggsave(file.path(out_dir, "plots", "all_methods_curves.png"), p1, width = 10, height = 14, dpi = 110)
ggsave(file.path(out_dir, "plots", "all_methods_bias.png"), p2, width = 9.5, height = 4.8, dpi = 120)
ggsave(file.path(out_dir, "plots", "all_methods_coverage.png"), p3, width = 9.5, height = 4.8, dpi = 120)
