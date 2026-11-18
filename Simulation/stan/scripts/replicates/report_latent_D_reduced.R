# Aggregate table and plots for the reduced-model study with wide priors
# (sim_latent_D_reduced_replicates.R, 20 replicates): oracle, naive, known errors, and the
# reduced latent model with the oracle and wild-guess sigma_Z priors.

library(ggplot2)
dir_in <- "output/stan/latent_D_car"
res <- readRDS(file.path(dir_in, "latent_D_reduced_replicates.rds"))


# Method labels
short <- c(
  "latent, reduced + CAR (oracle sigma_Z prior)" = "reduced, oracle sigma_Z prior",
  "latent, reduced + CAR (wild-guess sigma_Z prior)" = "reduced, wild-guess sigma_Z prior"
)
rename <- function(x) ifelse(x %in% names(short), short[x], x)
res$summary$method <- rename(res$summary$method)
res$diag$method <- rename(res$diag$method)

lv <- c(
  "oracle", "naive", "latent, known errors",
  "reduced, oracle sigma_Z prior", "reduced, wild-guess sigma_Z prior"
)


# Table
raw <- res$summary

# At each exposure value: relative bias (%), RMSE and coverage over the replicates
by_a <- do.call(rbind, lapply(split(raw, list(raw$method, raw$a), drop = TRUE), function(d) {
  data.frame(
    method = d$method[1],
    a = d$a[1],
    truth = d$true[1],
    mean_est = mean(d$est),
    rel_bias_pct = 100 * mean(d$est - d$true) / d$true[1],
    rmse = sqrt(mean((d$est - d$true)^2)),
    coverage = mean(d$lo <= d$true & d$true <= d$hi)
  )
}))

# Average over the grid, per method
tab <- aggregate(cbind(abs_rel_bias_pct = abs(rel_bias_pct), rmse, coverage) ~ method, by_a, mean)

# Convergence summary per method
dg <- do.call(rbind, lapply(split(res$diag, res$diag$method), function(d) {
  data.frame(
    method = d$method[1],
    reps_with_divergence = sum(d$div > 0),
    worst_rhat = max(d$max_rhat),
    median_rhat = median(d$max_rhat)
  )
}))

out <- merge(tab, dg, by = "method")
out <- cbind(n_reps = length(unique(raw$rep)), out)
out <- out[order(match(out$method, lv)), ]

dir.create(file.path(dir_in, "tables"), showWarnings = FALSE)
write.csv(out, file.path(dir_in, "tables", "aggregate_reduced_20reps.csv"), row.names = FALSE)
print(out, digits = 3, row.names = FALSE)


# Plots
# Leave known errors out of the plots
raw <- raw[raw$method != "latent, known errors", ]
by_a <- by_a[by_a$method != "latent, known errors", ]
raw$method <- factor(raw$method, levels = lv)
by_a$method <- factor(by_a$method, levels = lv)

n <- length(unique(raw$rep))
dir_out <- file.path(dir_in, "plots")
dir.create(dir_out, showWarnings = FALSE)

# Estimated curves: one line per replicate, with the mean of the estimates and the truth
p1 <- ggplot(raw, aes(a, est, group = rep)) +
  geom_line(alpha = 0.35, linewidth = 0.4) +
  geom_line(data = by_a, aes(y = mean_est, group = 1), colour = "blue", linetype = "dashed", linewidth = 0.8) +
  geom_line(data = by_a, aes(y = truth, group = 1), colour = "red", linewidth = 0.9) +
  facet_wrap(~method, ncol = 2) +
  labs(
    x = "exposure a", y = "expected counts per area",
    title = paste0("Estimated curves, reduced model with wide priors (", n, " replicates)"),
    subtitle = "Red = truth, blue dashed = mean of estimates, grey = individual replicates"
  ) +
  theme_bw(base_size = 11)

# Relative bias across the exposure grid
p2 <- ggplot(by_a, aes(a, rel_bias_pct, colour = method)) +
  geom_hline(yintercept = 0, colour = "grey50") +
  geom_line() +
  geom_point() +
  labs(
    x = "exposure a", y = "relative bias (%)",
    title = paste0("Relative bias, reduced model with wide priors (", n, " replicates)"),
    colour = NULL
  ) +
  theme_bw(base_size = 11)

# Coverage across the exposure grid
p3 <- ggplot(by_a, aes(a, coverage, colour = method)) +
  geom_hline(yintercept = 0.95, linetype = "dashed", colour = "grey50") +
  geom_line() +
  geom_point() +
  ylim(0, 1) +
  labs(
    x = "exposure a", y = "coverage of the 95% interval",
    title = paste0("Coverage, reduced model with wide priors (", n, " replicates, dashed = 0.95)"),
    colour = NULL
  ) +
  theme_bw(base_size = 11)

ggsave(file.path(dir_out, "reduced_curves.png"), p1, width = 10, height = 12, dpi = 120)
ggsave(file.path(dir_out, "reduced_bias.png"), p2, width = 8.5, height = 4.5, dpi = 120)
ggsave(file.path(dir_out, "reduced_coverage.png"), p3, width = 8.5, height = 4.5, dpi = 120)
