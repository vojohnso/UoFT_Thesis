# Combine the sigma_Z sensitivity runs into one table, one row per assumed sigma_Z.
# Bias is the relative bias (%) of the posterior mean of theta(a); the true sigma_Z is 1.

# Assumed sigma_Z values that have been run
values <- c(0.25, 0.5, 0.75, 1, 1.25, 1.5, 2)

rows <- lapply(values, function(sz) {
  f <- sprintf("output/stan/latent_D_car/latent_D_reduced_sens_sigmaZ_%s.rds", sz)
  if (!file.exists(f)) return(NULL)
  res <- readRDS(f)
  raw <- res$summary

  # Relative bias, RMSE and coverage at each exposure value
  by_a <- do.call(rbind, lapply(split(raw, raw$a), function(d) {
    data.frame(
      a = d$a[1],
      bias = 100 * mean(d$est - d$true) / d$true[1],
      rmse = sqrt(mean((d$est - d$true)^2)),
      cover = mean(d$lo <= d$true & d$true <= d$hi)
    )
  }))

  data.frame(
    sigma_Z = sz,
    n_reps = length(unique(raw$rep)),
    abs_bias = mean(abs(by_a$bias)),
    bias_a5 = by_a$bias[by_a$a == 5],
    bias_a10 = by_a$bias[by_a$a == 10],
    bias_a15 = by_a$bias[by_a$a == 15],
    rmse = mean(by_a$rmse),
    cover = mean(by_a$cover),
    div_reps = sum(res$diag$div > 0),
    worst_rhat = max(res$diag$max_rhat),
    med_rhat = median(res$diag$max_rhat)
  )
})

tab <- do.call(rbind, rows)
dir.create("output/stan/latent_D_car/tables", recursive = TRUE, showWarnings = FALSE)
write.csv(tab, "output/stan/latent_D_car/tables/sensitivity_sigmaZ.csv", row.names = FALSE)
print(tab, digits = 3, row.names = FALSE)
