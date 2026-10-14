library(ggplot2)

# Extra plots from the saved replicate results (output/replicates_<scenario>.rds):
#  1. fan plot: for each method, the mean estimate and the middle 95% of estimates across replicates
#  2. contrasts: bias and coverage of two contrasts computed from the posterior draws
# Contrasts use the nearest integer grid values: extreme vs typical = theta(15) - theta(10)
# (the 99th and 50th percentiles of the exposure are about 15.4 and 10) and interquartile = theta(12) - theta(8)
# (the 75th and 25th percentiles are about 11.6 and 8.4).
output_dir <- "output"
forms <- c("linear", "quadratic", "nonlinear")
cols <- c(naive = "#d6604d", oracle = "#2166ac")
res <- lapply(setNames(forms, forms), function(f) readRDS(file.path(output_dir, sprintf("replicates_%s.rds", f))))
summ <- do.call(rbind, lapply(res, `[[`, "summary"))
summ$form <- factor(summ$form, levels = forms)

# 1. Fan plot
key <- interaction(summ$form, summ$method, summ$a, drop = TRUE)
fan <- do.call(rbind, lapply(split(summ, key), function(d)
  data.frame(form = d$form[1], method = d$method[1], a = d$a[1], truth = d$true[1], mean = mean(d$est),
             lo = quantile(d$est, 0.025), hi = quantile(d$est, 0.975))))
p1 <- ggplot(fan, aes(a, colour = method, fill = method)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.18, colour = NA) +
  geom_line(aes(y = mean), linewidth = 0.9) +
  geom_line(aes(y = truth), colour = "black", linewidth = 0.8, linetype = "solid") +
  facet_wrap(~form, scales = "free_y") +
  scale_colour_manual(values = cols) + scale_fill_manual(values = cols) +
  labs(x = expression(paste("exposure ", a, " (", mu, "g/m"^3, ")")), y = expression(theta(a)), colour = NULL, fill = NULL,
       title = "Estimates across 200 replicates",
       subtitle = "Line: mean estimate; band: middle 95% of the estimates across replicates; black: population truth") +
  theme_bw(base_size = 12) + theme(legend.position = "bottom", panel.grid.minor = element_blank())
ggsave(file.path(output_dir, "replicates_fan.png"), p1, width = 11, height = 4.2, dpi = 130)

# 2. Contrasts from the posterior draws
contrasts <- list("extreme vs typical: theta(15) - theta(10)" = c("a=15", "a=10"),
                  "interquartile: theta(12) - theta(8)" = c("a=12", "a=8"))
rows <- list()
for (f in forms) {
  o <- res[[f]]; tr <- o$summary[o$summary$rep == o$rep[1] & o$summary$method == "oracle", c("a", "true")]
  for (cn in names(contrasts)) {
    ab <- contrasts[[cn]]; a1 <- as.numeric(sub("a=", "", ab[1])); a2 <- as.numeric(sub("a=", "", ab[2]))
    truth <- tr$true[tr$a == a1] - tr$true[tr$a == a2]
    for (m in c("naive", "oracle")) {
      for (k in seq_along(o$rep)) {
        d <- o$draws[, ab[1], m, k] - o$draws[, ab[2], m, k]
        rows[[length(rows) + 1]] <- data.frame(form = f, contrast = cn, method = m, rep = o$rep[k], truth = truth,
                                               est = mean(d), lo = quantile(d, 0.025), hi = quantile(d, 0.975))
      }
    }
  }
}
con <- do.call(rbind, rows)
con$form <- factor(con$form, levels = forms)
con$err <- con$est - con$truth
con$covered <- con$lo <= con$truth & con$truth <= con$hi
tab <- do.call(rbind, lapply(split(con, interaction(con$form, con$contrast, con$method, drop = TRUE)), function(d)
  data.frame(form = d$form[1], contrast = d$contrast[1], method = d$method[1], truth = d$truth[1],
             bias = mean(d$err), rel_bias_pct = 100 * mean(d$err) / d$truth[1], rmse = sqrt(mean(d$err^2)),
             coverage = mean(d$covered))))
rownames(tab) <- NULL
write.csv(tab, file.path(output_dir, "replicates_contrasts.csv"), row.names = FALSE)
print(format(tab, digits = 3), row.names = FALSE)

lab <- unique(tab[, c("form", "contrast", "method", "coverage")])
lab$text <- sprintf("coverage %.2f", lab$coverage)
p2 <- ggplot(con, aes(method, err / abs(truth) * 100, fill = method)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_boxplot(outlier.size = 0.7, alpha = 0.6) +
  geom_text(data = lab, aes(x = method, y = Inf, label = text), vjust = 1.4, size = 3, inherit.aes = FALSE) +
  facet_grid(contrast ~ form, scales = "free_y") +
  scale_fill_manual(values = cols) +
  labs(x = NULL, y = "estimate minus truth, as % of the true contrast", fill = NULL,
       title = "Contrast estimates across 200 replicates",
       subtitle = "Box: middle 50% of replicates; dashed line: no bias; text: coverage of the 95% credible interval") +
  theme_bw(base_size = 12) + theme(legend.position = "none")
ggsave(file.path(output_dir, "replicates_contrasts.png"), p2, width = 10, height = 6, dpi = 130)
