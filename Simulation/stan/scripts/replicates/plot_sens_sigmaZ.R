# Sensitivity plot: bias, RMSE and coverage of the reduced latent model against the
# assumed (fixed) sigma_Z. The true sigma_Z is 1.
#
# Reads  output/stan/latent_D_car/tables/sensitivity_sigmaZ.csv
#        (one row per assumed value, summarised from the sens_sigmaZ runs)
# Writes output/stan/latent_D_car/plots/sensitivity_sigmaZ.png

library(ggplot2)
d <- read.csv("output/stan/latent_D_car/tables/sensitivity_sigmaZ.csv")


# Panel 1: relative bias at a = 5, 10, 15 and averaged over the grid
bias <- rbind(
  data.frame(sigma_Z = d$sigma_Z, y = d$bias_a5, series = "a = 5"),
  data.frame(sigma_Z = d$sigma_Z, y = d$bias_a10, series = "a = 10"),
  data.frame(sigma_Z = d$sigma_Z, y = d$bias_a15, series = "a = 15"),
  data.frame(sigma_Z = d$sigma_Z, y = d$abs_bias, series = "average |bias|")
)
bias$series <- factor(bias$series, levels = c("a = 5", "a = 10", "a = 15", "average |bias|"))

p1 <- ggplot(bias, aes(sigma_Z, y, colour = series, linetype = series == "average |bias|")) +
  geom_hline(yintercept = 0, colour = "grey50") +
  geom_vline(xintercept = 1, linetype = "dotted") + # true sigma_Z
  geom_line() +
  geom_point() +
  scale_linetype_manual(values = c("solid", "dashed"), guide = "none") +
  scale_x_continuous(breaks = d$sigma_Z) +
  labs(
    x = "assumed sigma_Z (fixed)", y = "relative bias (%)", colour = NULL,
    title = "Bias against the assumed sigma_Z",
    subtitle = "Truth sigma_Z = 1 (dotted); 20 replicates, linear truth, 100 areas"
  ) +
  theme_bw(base_size = 11)


# Panel 2: RMSE
p2 <- ggplot(d, aes(sigma_Z, rmse)) +
  geom_vline(xintercept = 1, linetype = "dotted") +
  geom_line() +
  geom_point() +
  scale_x_continuous(breaks = d$sigma_Z) +
  labs(x = "assumed sigma_Z (fixed)", y = "RMSE", title = "RMSE") +
  theme_bw(base_size = 11)


# Panel 3: coverage
p3 <- ggplot(d, aes(sigma_Z, cover)) +
  geom_hline(yintercept = 0.95, linetype = "dashed", colour = "grey50") +
  geom_vline(xintercept = 1, linetype = "dotted") +
  geom_line() +
  geom_point() +
  ylim(0, 1) +
  scale_x_continuous(breaks = d$sigma_Z) +
  labs(x = "assumed sigma_Z (fixed)", y = "coverage", title = "Coverage (dashed = 0.95)") +
  theme_bw(base_size = 11)


# Save: one combined figure if patchwork is available, else three files
if (requireNamespace("patchwork", quietly = TRUE)) {
  library(patchwork)
  p <- p1 / (p2 | p3) + plot_layout(heights = c(1.5, 1))
  ggsave("output/stan/latent_D_car/plots/sensitivity_sigmaZ.png", p, width = 9, height = 9, dpi = 120)
} else {
  ggsave("output/stan/latent_D_car/plots/sensitivity_sigmaZ_bias.png", p1, width = 8, height = 4.5, dpi = 120)
  ggsave("output/stan/latent_D_car/plots/sensitivity_sigmaZ_rmse.png", p2, width = 6, height = 4, dpi = 120)
  ggsave("output/stan/latent_D_car/plots/sensitivity_sigmaZ_coverage.png", p3, width = 6, height = 4, dpi = 120)
}
