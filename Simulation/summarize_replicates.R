library(ggplot2)

# Reads output/replicates_<scenario>.rds (from sim_replicates.R) and computes bias, RMSE and 95%
# credible-interval coverage.
# Truth = population ERF (true_erf_pop); coverage against the per-replicate sample truth is also reported.
# Optional upper limit on the reported nonlinear grid, report_max_nonlinear (15 = top of the grid, i.e. no limit).
output_dir <- "output"
report_max_nonlinear <- 15

forms <- c("linear", "quadratic", "nonlinear")
raw <- do.call(rbind, lapply(forms, function(f)
  readRDS(file.path(output_dir, sprintf("replicates_%s.rds", f)))$summary))
n_reps <- tapply(raw$rep, raw$form, function(r) length(unique(r)))
cat("Replicates:", paste(names(n_reps), n_reps, collapse = ", "), "\n")

raw <- raw[!(raw$form == "nonlinear" & raw$a > report_max_nonlinear), ]

summ <- do.call(rbind, lapply(split(raw, list(raw$form, raw$method, raw$a), drop = TRUE), function(d) {
  cover <- mean(d$lo <= d$true & d$true <= d$hi)
  data.frame(form = d$form[1], method = d$method[1], a = d$a[1],
             n_reps = nrow(d),
             true = mean(d$true),
             mean_est = mean(d$est),
             bias = mean(d$est - d$true),
             rel_bias_pct = 100 * mean(d$est - d$true) / mean(d$true),
             rmse = sqrt(mean((d$est - d$true)^2)),
             coverage = cover,
             mc_se_cov = sqrt(cover * (1 - cover) / nrow(d)),
             coverage_vs_sample_truth = mean(d$lo <= d$true_sample & d$true_sample <= d$hi),
             ci_width = mean(d$hi - d$lo))
}))
summ <- summ[order(summ$form, summ$method, summ$a), ]
rownames(summ) <- NULL
write.csv(summ, file.path(output_dir, "replicates_summary.csv"), row.names = FALSE)

for (form in forms) {
  cat("\n", toupper(form), "\n")
  print(format(summ[summ$form == form, -1], digits = 3), row.names = FALSE)
}
cat("\nMean coverage across grid (vs population truth / vs sample truth):\n")
print(aggregate(cbind(coverage, coverage_vs_sample_truth) ~ form + method, data = summ, FUN = mean),
      digits = 3)

summ$form <- factor(summ$form, levels = forms)
cols <- c(oracle = "#2166ac", naive = "#d6604d")
p_bias <- ggplot(summ, aes(a, rel_bias_pct, colour = method)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_line() + geom_point() + facet_wrap(~form, scales = "free") +
  scale_colour_manual(values = cols) +
  labs(x = "PM2.5 (a)", y = "Relative bias of theta(a) (%)", colour = NULL,
       title = sprintf("Bias over %d replicates (population truth)", max(n_reps))) +
  theme_bw() + theme(legend.position = "bottom")
p_cov <- ggplot(summ, aes(a, coverage, colour = method)) +
  geom_hline(yintercept = 0.95, linetype = "dashed") +
  geom_line() + geom_point() + facet_wrap(~form, scales = "free_x") + ylim(0, 1) +
  scale_colour_manual(values = cols) +
  labs(x = "PM2.5 (a)", y = "95% credible-interval coverage", colour = NULL,
       title = sprintf("Coverage over %d replicates (population truth)", max(n_reps))) +
  theme_bw() + theme(legend.position = "bottom")
ggsave(file.path(output_dir, "replicates_bias.png"), p_bias, width = 11, height = 4, dpi = 130)
ggsave(file.path(output_dir, "replicates_coverage.png"), p_cov, width = 11, height = 4, dpi = 130)
cat("\nSaved summary and plots to", output_dir, "\n")
