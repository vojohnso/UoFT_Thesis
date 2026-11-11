# One aggregate table and one set of plots for all methods on the same 20 datasets (linear truth, 100 areas, spatial effect in the outcome):
#   with the CAR spatial term: oracle, naive, known errors (error parameters given), reduced latent model with sigma_Z fixed at 1,
#   with the oracle prior normal(1, 0.3), with the wild-guess prior normal(0, 5), and with tight priors on every parameter.
# Table: output/stan/latent_D_car/tables/aggregate_all_methods.csv; plots: plots/all_methods_{curves,bias,coverage}.png
library(ggplot2)
d <- "output/stan/latent_D_car"
get <- function(file, keep = NULL, rename = NULL) {
  r <- readRDS(file.path(d, file)); s <- r$summary; g <- r$diag
  if (!is.null(keep)) { s <- s[s$method %in% keep, ]; g <- g[g$method %in% keep, ] }
  if (!is.null(rename)) for (k in names(rename)) { s$method[s$method == k] <- rename[[k]]; g$method[g$method == k] <- rename[[k]] }
  list(s = s[, c("rep", "method", "a", "true", "est", "lo", "hi")], g = g[, c("rep", "method", "div", "max_rhat")])
}
parts <- list(
  get("known_exposure_car_replicates.rds", rename = c("oracle + CAR" = "oracle", "naive + CAR" = "naive")),
  get("known_errors_car_replicates.rds", rename = c("latent, known errors + CAR" = "latent, known errors")),
  get("latent_D_reduced_fixedsigmaZ_replicates.rds", rename = c("reduced, sigma_Z fixed at 1" = "reduced, sigma_Z fixed at 1")),
  get("latent_D_reduced_replicates.rds", keep = c("latent, reduced + CAR (oracle sigma_Z prior)", "latent, reduced + CAR (wild-guess sigma_Z prior)"),
      rename = c("latent, reduced + CAR (oracle sigma_Z prior)" = "reduced, sigma_Z ~ normal(1, 0.3)", "latent, reduced + CAR (wild-guess sigma_Z prior)" = "reduced, sigma_Z wild guess")),
  get("latent_D_reduced_tight_replicates.rds", rename = c("reduced, tight priors on everything" = "reduced, tight priors on everything")))
# Main study extension (replicates 21-100, sim_main_extend.R), same method names; the tight-priors row stays at the original replicates
for (f in c("main_reps_21_60.rds", "main_reps_61_100.rds")) if (file.exists(file.path(d, f))) parts[[length(parts) + 1]] <- get(f)
raw <- do.call(rbind, lapply(parts, `[[`, "s")); dg <- do.call(rbind, lapply(parts, `[[`, "g"))
lv <- c("oracle", "naive", "latent, known errors", "reduced, sigma_Z fixed at 1", "reduced, tight priors on everything",
        "reduced, sigma_Z ~ normal(1, 0.3)", "reduced, sigma_Z wild guess")
stopifnot(all(raw$method %in% lv)); raw$method <- factor(raw$method, levels = lv)

by_a <- do.call(rbind, lapply(split(raw, list(raw$method, raw$a), drop = TRUE), function(x)
  data.frame(method = x$method[1], a = x$a[1], truth = x$true[1], mean_est = mean(x$est), rel_bias_pct = 100 * mean(x$est - x$true) / x$true[1],
             rmse = sqrt(mean((x$est - x$true)^2)), coverage = mean(x$lo <= x$true & x$true <= x$hi))))
tab <- aggregate(cbind(abs_rel_bias_pct = abs(rel_bias_pct), rmse, coverage) ~ method, by_a, mean)
dgt <- do.call(rbind, lapply(split(dg, dg$method), function(x) data.frame(method = x$method[1], reps_with_divergence = sum(x$div > 0),
        worst_rhat = max(x$max_rhat), median_rhat = median(x$max_rhat))))
nrep <- aggregate(rep ~ method, raw, function(x) length(unique(x))); names(nrep)[2] <- "n_reps"
out <- merge(merge(nrep, tab, by = "method"), dgt, by = "method"); out <- out[order(match(out$method, lv)), ]
out$mc_se_coverage <- sqrt(out$coverage * (1 - out$coverage) / out$n_reps)
dir.create(file.path(d, "tables"), showWarnings = FALSE); dir.create(file.path(d, "plots"), showWarnings = FALSE)
write.csv(out, file.path(d, "tables", "aggregate_all_methods.csv"), row.names = FALSE); print(out, digits = 3, row.names = FALSE)

by_a$method <- factor(by_a$method, levels = lv); n <- paste(range(out$n_reps), collapse = "-")
p1 <- ggplot(raw, aes(a, est, group = rep)) + geom_line(alpha = 0.3, linewidth = 0.35) +
  geom_line(data = by_a, aes(y = mean_est, group = 1), colour = "blue", linetype = "dashed", linewidth = 0.8) +
  geom_line(data = by_a, aes(y = truth, group = 1), colour = "red", linewidth = 0.9) + facet_wrap(~method, ncol = 2) + coord_cartesian(ylim = c(0, 3500)) +
  labs(x = "exposure a", y = "expected counts per area", title = paste0("Estimated curves, all methods (", n, " replicates)"),
       subtitle = "Red = truth, blue dashed = mean of estimates, grey = replicates (y-axis cut at 3500)") + theme_bw(base_size = 10)
p2 <- ggplot(by_a, aes(a, rel_bias_pct, colour = method)) + geom_hline(yintercept = 0, colour = "grey50") + geom_line() + geom_point(size = 1) +
  labs(x = "exposure a", y = "relative bias (%)", title = paste0("Relative bias, all methods (", n, " replicates)"), colour = NULL) + theme_bw(base_size = 10)
p3 <- ggplot(by_a, aes(a, coverage, colour = method)) + geom_hline(yintercept = 0.95, linetype = "dashed", colour = "grey50") + geom_line() + geom_point(size = 1) +
  ylim(0, 1) + labs(x = "exposure a", y = "coverage of the 95% interval", title = paste0("Coverage, all methods (", n, " replicates, dashed = 0.95)"), colour = NULL) +
  theme_bw(base_size = 10)
ggsave(file.path(d, "plots", "all_methods_curves.png"), p1, width = 10, height = 14, dpi = 110)
ggsave(file.path(d, "plots", "all_methods_bias.png"), p2, width = 9.5, height = 4.8, dpi = 120)
ggsave(file.path(d, "plots", "all_methods_coverage.png"), p3, width = 9.5, height = 4.8, dpi = 120)
