# Aggregate tables (bias, RMSE, coverage averaged over the exposure grid, plus sampler diagnostics) for the linear-truth latent studies:
# the anchored study (20 replicates), the free-slope study (10), and the two latent models on their 10 shared datasets
dir_in <- "output/stan/latent_D_car"
dir_out <- file.path(dir_in, "tables")
dir.create(dir_out, showWarnings = FALSE)

aggregate_study <- function(res, label) {
  raw <- res$summary
  by_a <- do.call(rbind, lapply(split(raw, list(raw$method, raw$a), drop = TRUE), function(d)
    data.frame(method = d$method[1], rel_bias_pct = 100 * mean(d$est - d$true) / d$true[1],
               rmse = sqrt(mean((d$est - d$true)^2)), coverage = mean(d$lo <= d$true & d$true <= d$hi))))
  tab <- aggregate(cbind(abs_rel_bias_pct = abs(rel_bias_pct), rmse, coverage) ~ method, by_a, mean)
  dg <- do.call(rbind, lapply(split(res$diag, res$diag$method), function(d)
    data.frame(method = d$method[1], reps_with_divergence = sum(d$div > 0), worst_rhat = max(d$max_rhat), median_rhat = median(d$max_rhat))))
  out <- merge(tab, dg, by = "method")
  cbind(study = label, n_reps = length(unique(raw$rep)), out)
}

free <- readRDS(file.path(dir_in, "latent_D_replicates.rds"))
anch <- readRDS(file.path(dir_in, "latent_D_anchored_replicates.rds"))
free$summary$method[free$summary$method == "latent, estimated errors + CAR"] <- "latent, free slope + CAR"
free$diag$method[free$diag$method == "latent, estimated errors + CAR"] <- "latent, free slope + CAR"

# the first 10 anchored replicates use the same seeds, hence the same datasets, as the free-slope study
anch10 <- anch
anch10$summary <- anch$summary[anch$summary$rep <= 10, ]
anch10$diag <- anch$diag[anch$diag$rep <= 10, ]

tabs <- list(anchored_20reps = aggregate_study(anch, "anchored, 20 replicates"),
             freeslope_10reps = aggregate_study(free, "free slope, 10 replicates"),
             anchored_first10reps = aggregate_study(anch10, "anchored, first 10 replicates"))
for (n in names(tabs)) {
  write.csv(tabs[[n]], file.path(dir_out, paste0("aggregate_", n, ".csv")), row.names = FALSE)
  cat("\n", tabs[[n]]$study[1], "\n"); print(tabs[[n]][, -1], digits = 3, row.names = FALSE)
}

# Comparison table: oracle, naive, known-errors and anchored rows come from the 20-replicate anchored study; the free-slope study only
# has 10 replicates, so that row keeps its own n_reps
cmp <- rbind(tabs$anchored_20reps, tabs$freeslope_10reps[tabs$freeslope_10reps$method == "latent, free slope + CAR", ])[, -1]
cmp$method <- factor(cmp$method, levels = c("oracle", "naive", "latent, known errors", "latent, free slope + CAR", "latent, anchored + CAR"))
cmp <- cmp[order(cmp$method), ]
write.csv(cmp, file.path(dir_out, "aggregate_freeslope_vs_anchored.csv"), row.names = FALSE)
cat("\nfree slope vs anchored (free slope has 10 replicates, the others 20)\n"); print(cmp, digits = 3, row.names = FALSE)
