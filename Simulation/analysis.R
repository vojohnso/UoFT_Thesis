source("data_generate.R")
library(ggplot2)
library(gridExtra)

# Load simulation results saved by run_simulation.R
all_metrics <- readRDS("sim_results.rds")
forms <- names(all_metrics)

# Table 1: Performance metrics averaged across exposure grid
# For each scenario report mean bias, RMSE, coverage, width
# for both naive and oracle

make_summary_table <- function(all_metrics) {
  rows <- list()
  for (form in forms) {
    for (method in c("naive", "oracle")) {
      m <- all_metrics[[form]][[method]]
      rows[[length(rows) + 1]] <- data.frame(
        out_form = form,
        method = method,
        mean_bias = mean(m$bias),
        mean_abs_bias = mean(abs(m$bias)),
        mean_rmse = mean(m$rmse),
        mean_coverage = mean(m$coverage),
        mean_width = mean(m$width)
      )
    }
  }
  do.call(rbind, rows)
}

summary_table <- make_summary_table(all_metrics)
cat("=== Table 1: Average performance metrics across exposure grid ===\n")
print(summary_table, digits = 3, row.names = FALSE)

# Table 2: Performance at a single exposure level a = 50
# More interpretable than averages across the full grid

make_pointwise_table <- function(all_metrics, a_target = 50) {
  rows <- list()
  for (form in forms) {
    for (method in c("naive", "oracle")) {
      m <- all_metrics[[form]][[method]]
      sub <- m[m$a == a_target, ]
      rows[[length(rows) + 1]] <- data.frame(
        out_form = form,
        method = method,
        a = a_target,
        bias = sub$bias,
        rmse = sub$rmse,
        coverage = sub$coverage,
        width = sub$width
      )
    }
  }
  do.call(rbind, rows)
}

pointwise_table <- make_pointwise_table(all_metrics, a_target = 50)
cat("\n=== Table 2: Performance at a = 50 ug/m3 ===\n")
print(pointwise_table, digits = 3, row.names = FALSE)

# Figure 1: Bias across exposure levels
# One panel per outcome form, naive vs oracle

make_bias_df <- function(all_metrics) {
  rows <- list()
  for (form in forms) {
    for (method in c("naive", "oracle")) {
      m <- all_metrics[[form]][[method]]
      rows[[length(rows) + 1]] <- cbind(m, out_form = form, method = method)
    }
  }
  do.call(rbind, rows)
}

metric_df <- make_bias_df(all_metrics)
metric_df$method <- factor(metric_df$method,
                           levels = c("naive", "oracle"),
                           labels = c("Naive (surrogate Z)", "Oracle (true A)"))
metric_df$out_form <- factor(metric_df$out_form,
                             levels = c("linear", "quadratic", "nonlinear"),
                             labels = c("Linear", "Quadratic", "Nonlinear"))

fig1 <- ggplot(metric_df, aes(x = a, y = bias, colour = method)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  facet_wrap(~out_form, nrow = 1) +
  scale_colour_manual(values = c(
    "Naive (surrogate Z)" = "#d6604d",
    "Oracle (true A)" = "#2166ac")) +
  labs(
    x = expression(paste("PM"[2.5], " (", mu, "g/m"^3, ")")),
    y = "Bias",
    colour = NULL,
    title = "Figure 1: Bias in average potential outcome across exposure levels",
    subtitle = "Dashed line at zero indicates no bias"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

print(fig1)
ggsave("fig1_bias.png", fig1, width = 10, height = 4, dpi = 150)

# Figure 2: 95% CI coverage across exposure levels
# Dashed line at nominal 0.95

fig2 <- ggplot(metric_df, aes(x = a, y = coverage, colour = method)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  geom_hline(yintercept = 0.95, linetype = "dashed", colour = "grey50") +
  facet_wrap(~out_form, nrow = 1) +
  scale_colour_manual(values = c(
    "Naive (surrogate Z)" = "#d6604d",
    "Oracle (true A)" = "#2166ac")) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
  labs(
    x = expression(paste("PM"[2.5], " (", mu, "g/m"^3, ")")),
    y = "Coverage",
    colour = NULL,
    title = "Figure 2: 95% CI coverage across exposure levels",
    subtitle = "Dashed line at nominal 0.95"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

print(fig2)
ggsave("fig2_coverage.png", fig2, width = 10, height = 4, dpi = 150)

# Figure 3: RMSE across exposure levels

fig3 <- ggplot(metric_df, aes(x = a, y = rmse, colour = method)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  facet_wrap(~out_form, nrow = 1) +
  scale_colour_manual(values = c(
    "Naive (surrogate Z)" = "#d6604d",
    "Oracle (true A)" = "#2166ac")) +
  labs(
    x = expression(paste("PM"[2.5], " (", mu, "g/m"^3, ")")),
    y = "RMSE",
    colour = NULL,
    title = "Figure 3: RMSE across exposure levels"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

print(fig3)
ggsave("fig3_rmse.png", fig3, width = 10, height = 4, dpi = 150)

# Figure 4: Credible interval width across exposure levels
# Wider intervals = more uncertainty

fig4 <- ggplot(metric_df, aes(x = a, y = width, colour = method)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  facet_wrap(~out_form, nrow = 1) +
  scale_colour_manual(values = c(
    "Naive (surrogate Z)" = "#d6604d",
    "Oracle (true A)" = "#2166ac")) +
  labs(
    x = expression(paste("PM"[2.5], " (", mu, "g/m"^3, ")")),
    y = "Mean CI width",
    colour = NULL,
    title = "Figure 4: Mean 95% credible interval width across exposure levels",
    subtitle = "Wider intervals indicate more uncertainty from measurement error"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

print(fig4)
ggsave("fig4_width.png", fig4, width = 10, height = 4, dpi = 150)

# Figure 5: Single dataset ERF illustration
# Show one example dataset with true ERF, oracle, naive
# One panel per outcome form

a_grid <- seq(0, 100, by = 10)

erf_rows <- list()
for (form in c("linear", "quadratic", "nonlinear")) {
  sim <- data.generate(seed = 401, n_areas = 100, out_form = form)
  true_mu <- true_erf(a_grid, sim)
  naive <- bart_gcomp(sim, a_grid = a_grid, use_true_A = FALSE,
                      ndpost = 500, nskip = 250)
  oracle <- bart_gcomp(sim, a_grid = a_grid, use_true_A = TRUE,
                       ndpost = 500, nskip = 250)
  
  erf_rows[[length(erf_rows) + 1]] <- rbind(
    data.frame(a = a_grid, mean = oracle$mu_summary$mean,
               q025 = oracle$mu_summary$q025, q975 = oracle$mu_summary$q975,
               method = "Oracle (true A)", out_form = form),
    data.frame(a = a_grid, mean = naive$mu_summary$mean,
               q025 = naive$mu_summary$q025, q975 = naive$mu_summary$q975,
               method = "Naive (surrogate Z)", out_form = form),
    data.frame(a = a_grid, mean = true_mu,
               q025 = true_mu, q975 = true_mu,
               method = "True ERF", out_form = form)
  )
}

erf_df <- do.call(rbind, erf_rows)
erf_df$out_form <- factor(erf_df$out_form,
                          levels = c("linear", "quadratic", "nonlinear"),
                          labels = c("Linear", "Quadratic", "Nonlinear"))

fig5 <- ggplot(erf_df, aes(x = a, y = mean, colour = method, fill = method)) +
  geom_ribbon(data = subset(erf_df, method != "True ERF"),
              aes(ymin = q025, ymax = q975), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  facet_wrap(~out_form, nrow = 1) +
  scale_colour_manual(values = c(
    "True ERF" = "black",
    "Oracle (true A)" = "#2166ac",
    "Naive (surrogate Z)" = "#d6604d")) +
  scale_fill_manual(values = c(
    "True ERF" = "black",
    "Oracle (true A)" = "#2166ac",
    "Naive (surrogate Z)" = "#d6604d")) +
  scale_x_continuous(breaks = seq(0, 100, by = 20)) +
  labs(
    x = expression(paste("Wildfire Smoke PM"[2.5], " (", mu, "g/m"^3, ")")),
    y = "Average Event Rate",
    colour = NULL, fill = NULL,
    title = "Figure 5: Estimated ERF for one simulated dataset",
    subtitle = "Ribbons show 95% posterior credible intervals"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

print(fig5)
ggsave("fig5_erf_illustration.png", fig5, width = 12, height = 4.5, dpi = 150)

cat("\nAll figures saved.\n")