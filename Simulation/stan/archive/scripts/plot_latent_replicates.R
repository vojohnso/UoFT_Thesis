# Plots for the latent exposure replicate studies: estimated curves against the truth, relative bias and coverage by exposure value
# Reads the .rds files written by sim_latent_D_replicates.R (free slope), sim_latent_D_anchored_replicates.R and sim_latent_D_spline_replicates.R
library(ggplot2)

dir_in <- "output/stan/latent_D_car"
dir_out <- file.path(dir_in, "plots")
dir.create(dir_out, showWarnings = FALSE)

studies <- c(freeslope_linear = "latent_D_replicates.rds",
             anchored_linear = "latent_D_anchored_replicates.rds",
             spline_quadratic = "latent_D_spline_replicates.rds",
             spline_linear = "latent_D_spline_replicates_linear.rds")

summarise_by_a <- function(raw) do.call(rbind, lapply(split(raw, list(raw$method, raw$a), drop = TRUE), function(d)
  data.frame(method = d$method[1], a = d$a[1], truth = d$true[1], mean_est = mean(d$est),
             rel_bias_pct = 100 * mean(d$est - d$true) / d$true[1], coverage = mean(d$lo <= d$true & d$true <= d$hi))))

make_plots <- function(raw, s) {
  n_reps <- length(unique(raw$rep))
  by_a <- summarise_by_a(raw)
  p1 <- ggplot(raw, aes(a, est, group = rep)) + geom_line(alpha = 0.35, linewidth = 0.4) +
    geom_line(data = by_a, aes(y = mean_est, group = 1), colour = "blue", linetype = "dashed", linewidth = 0.8) +
    geom_line(data = by_a, aes(y = truth, group = 1), colour = "red", linewidth = 0.9) +
    facet_wrap(~method) + labs(x = "exposure a", y = "expected counts per area",
      title = paste0("Estimated curves, ", s, " (", n_reps, " replicates)"), subtitle = "Red = truth, blue dashed = mean of estimates, grey = individual replicates") +
    theme_bw(base_size = 11)
  p2 <- ggplot(by_a, aes(a, rel_bias_pct, colour = method)) + geom_hline(yintercept = 0, colour = "grey50") + geom_line() + geom_point() +
    labs(x = "exposure a", y = "relative bias (%)", title = paste0("Relative bias, ", s, " (", n_reps, " replicates)"), colour = NULL) + theme_bw(base_size = 11)
  p3 <- ggplot(by_a, aes(a, coverage, colour = method)) + geom_hline(yintercept = 0.95, linetype = "dashed", colour = "grey50") +
    geom_line() + geom_point() + ylim(0, 1) +
    labs(x = "exposure a", y = "coverage of the 95% interval", title = paste0("Coverage, ", s, " (", n_reps, " replicates, dashed = 0.95)"), colour = NULL) +
    theme_bw(base_size = 11)
  ggsave(file.path(dir_out, paste0(s, "_curves.png")), p1, width = 10, height = 4.2 * ceiling(length(unique(raw$method)) / 3), dpi = 120)
  ggsave(file.path(dir_out, paste0(s, "_bias.png")), p2, width = 7.5, height = 4.2, dpi = 120)
  ggsave(file.path(dir_out, paste0(s, "_coverage.png")), p3, width = 7.5, height = 4.2, dpi = 120)
  cat("Saved plots for", s, "\n")
}

for (s in names(studies)) {
  f <- file.path(dir_in, studies[[s]])
  if (!file.exists(f)) next
  make_plots(readRDS(f)$summary, s)
}

# Free-slope vs anchored latent model on the same datasets: the first 10 replicates of the anchored study use the same seeds
# as the free-slope study, so the oracle and naive fits are identical in both. The known-errors model is left out.
free <- readRDS(file.path(dir_in, "latent_D_replicates.rds"))$summary
anch <- readRDS(file.path(dir_in, "latent_D_anchored_replicates.rds"))$summary
anch <- anch[anch$rep <= max(free$rep), ]
free$method[free$method == "latent, estimated errors + CAR"] <- "latent, free slope + CAR"
anch$method[anch$method == "latent, anchored + CAR"] <- "latent, anchored + CAR"
cols <- c("rep", "method", "a", "true", "est", "lo", "hi")
cmp <- rbind(free[free$method == "latent, free slope + CAR", cols], anch[anch$method != "latent, known errors", cols])
cmp$method <- factor(cmp$method, levels = c("oracle", "naive", "latent, free slope + CAR", "latent, anchored + CAR"))
make_plots(cmp, "freeslope_vs_anchored_linear")
