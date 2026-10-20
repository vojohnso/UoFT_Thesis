# Summarises output/stan/latent_replicates.rds: bias, RMSE and 95% coverage (population truth and dataset truth)
# for Oracle, Naive and Latent, plus parameter recovery and sampler diagnostics.
# Run from the Simulation/ folder:  Rscript stan/scripts/summarize_latent_replicates.R
library(ggplot2)
out <- readRDS("output/stan/latent_replicates.rds")
raw <- out$summary
n <- length(unique(raw$rep))
cat("Replicates:", n, "\n")

summ <- do.call(rbind, lapply(split(raw, list(raw$method, raw$a), drop = TRUE), function(d) {
  cov_pop <- mean(d$lo <= d$true & d$true <= d$hi)
  data.frame(method = d$method[1], a = d$a[1], true = d$true[1], mean_est = mean(d$est),
             rel_bias_pct = 100 * mean(d$est - d$true) / d$true[1],
             rmse = sqrt(mean((d$est - d$true)^2)),
             coverage_pop = cov_pop, coverage_sample = mean(d$lo <= d$true_sample & d$true_sample <= d$hi),
             ci_width = mean(d$hi - d$lo))
}))
summ <- summ[order(summ$method, summ$a), ]; rownames(summ) <- NULL
write.csv(summ, "output/stan/latent_replicates_summary.csv", row.names = FALSE)

cat("\nAveraged over the grid (coverage SE about", round(sqrt(.95 * .05 / n), 3), "per grid value):\n")
print(aggregate(cbind(abs_rel_bias_pct = abs(rel_bias_pct), rmse, coverage_pop, coverage_sample, ci_width) ~ method,
                data = summ, FUN = mean), digits = 3)

cat("\nLatent model parameters (truth: beta_A = ", out$truth_beta_A, ", sigma_A = ", out$truth_sigma_A, "):\n", sep = "")
lp <- out$latent_pars
cat(sprintf("  beta_A : mean %.3f, SD across replicates %.3f, 95%% interval covers truth in %.2f\n", mean(lp$beta_A),
            sd(lp$beta_A), mean(lp$beta_A_lo <= out$truth_beta_A & out$truth_beta_A <= lp$beta_A_hi)))
cat(sprintf("  sigma_A: mean %.3f, SD across replicates %.3f, 95%% interval covers truth in %.2f\n", mean(lp$sigma_A),
            sd(lp$sigma_A), mean(lp$sigma_A_lo <= out$truth_sigma_A & out$truth_sigma_A <= lp$sigma_A_hi)))

cat("\nSampler diagnostics (replicates with any divergence, max R-hat):\n")
print(aggregate(cbind(any_div = div > 0, max_rhat) ~ method, data = out$diag,
                FUN = function(x) if (is.logical(x)) sum(x) else max(x)), digits = 4)

cols <- c(oracle = "#2166ac", naive = "#d6604d", latent = "#1b7837")
p_bias <- ggplot(summ, aes(a, rel_bias_pct, colour = method)) + geom_hline(yintercept = 0, linetype = "dashed") +
  geom_line() + geom_point() + scale_colour_manual(values = cols) +
  labs(x = "exposure a", y = "relative bias of theta(a) (%)", colour = NULL,
       title = sprintf("Latent model: bias over %d replicates (population truth)", n)) +
  theme_bw() + theme(legend.position = "bottom")
p_cov <- ggplot(summ, aes(a, coverage_pop, colour = method)) + geom_hline(yintercept = 0.95, linetype = "dashed") +
  geom_line() + geom_point() + ylim(0, 1) + scale_colour_manual(values = cols) +
  labs(x = "exposure a", y = "95% interval coverage", colour = NULL,
       title = sprintf("Latent model: coverage over %d replicates (population truth)", n)) +
  theme_bw() + theme(legend.position = "bottom")
ggsave("output/stan/latent_replicates_bias.png", p_bias, width = 7, height = 4.5, dpi = 130)
ggsave("output/stan/latent_replicates_coverage.png", p_cov, width = 7, height = 4.5, dpi = 130)
