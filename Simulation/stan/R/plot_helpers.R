library(ggplot2)

# Plots theta(a) curves with 95% posterior bands, for any of the models.
# est:   data frame with columns a, method, mean, lo, hi (lo and hi are NA for point estimates only),
#        and optionally scenario (one panel per scenario)
# truth: data frame with columns a, truth ("population truth" or "this dataset's truth"), value,
#        and optionally scenario
plot_erf <- function(est, truth, title = NULL, subtitle = NULL) {
  cols <- c(Oracle = "#2166ac", Naive = "#d6604d", Latent = "#1b7837", Stan = "#2166ac", BART = "#e08214")
  p <- ggplot() +
    geom_ribbon(data = est[!is.na(est$lo), ], aes(a, ymin = lo, ymax = hi, fill = method), alpha = 0.15) +
    geom_line(data = est, aes(a, mean, colour = method), linewidth = 0.9) +
    geom_line(data = truth, aes(a, value, linetype = truth), colour = "black", linewidth = 0.8) +
    scale_colour_manual(values = cols) +
    scale_fill_manual(values = cols) +
    scale_linetype_manual(values = c("population truth" = "solid", "this dataset's truth" = "dashed")) +
    labs(x = expression(paste("exposure ", a, " (", mu, "g/m"^3, ")")),
         y = expression(theta(a)), colour = NULL, fill = NULL, linetype = NULL,
         title = title, subtitle = subtitle) +
    theme_bw(base_size = 12) +
    theme(legend.position = "bottom", panel.grid.minor = element_blank())
  if ("scenario" %in% names(est)) p <- p + facet_wrap(~scenario, scales = "free_y")
  p
}
