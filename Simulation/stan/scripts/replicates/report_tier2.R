# Tier 2 table: per setting and method, |bias| %, RMSE, coverage (+MC SE), convergence (R-hat, Geweke, divergences). Writes tables/aggregate_tier2.csv
d <- "output/stan/latent_D_car"; set <- c("sz0.5", "sz2", "rho0.5")
read_set <- function(s) { fs <- list.files(d, sprintf("^tier2_%s_reps_.*rds$", s), full.names = TRUE); l <- lapply(fs, readRDS)
  list(summary = do.call(rbind, lapply(l, `[[`, "summary")), diag = do.call(rbind, lapply(l, `[[`, "diag"))) }
out <- do.call(rbind, lapply(set, function(s) { x <- read_set(s); sm <- x$summary; dg <- x$diag
  sm$err <- sm$est - sm$true; sm$cov <- sm$lo <= sm$true & sm$true <= sm$hi
  per <- aggregate(cbind(err, cov) ~ rep + method, sm, mean)
  # |relative bias| of the mean estimate over reps, then averaged over the grid (as in the other tables)
  m <- aggregate(est ~ method + a, sm, mean); tr <- aggregate(true ~ a, sm, mean); m <- merge(m, tr, by = "a")
  bias <- aggregate(I(abs(est - true) / abs(true) * 100) ~ method, m, mean); names(bias)[2] <- "abs_rel_bias_pct"
  rm <- aggregate(I((est - true)^2) ~ method, sm, mean); rm[[2]] <- sqrt(rm[[2]]); names(rm)[2] <- "rmse"
  cv <- aggregate(cov ~ method, sm, mean)
  n <- aggregate(rep ~ method, sm, function(v) length(unique(v))); names(n)[2] <- "n_reps"
  g <- aggregate(cbind(rhat_key, geweke_frac_flag, geweke_max_abs_z) ~ method, dg, median)
  c1 <- aggregate(cbind(div > 0, rhat_key > 1.05, geweke_frac_flag > 0.1) ~ method, dg, mean); names(c1)[2:4] <- c("frac_reps_divergent", "frac_reps_rhat_gt_1.05", "frac_reps_geweke_flag_gt_10pct")
  r <- Reduce(function(a, b) merge(a, b, by = "method"), list(n, bias, rm, cv, g, c1)); r$mc_se_cov <- sqrt(r$cov * (1 - r$cov) / r$n_reps); cbind(setting = s, r) }))
out <- out[order(match(out$setting, set)), ]; dir.create(file.path(d, "tables"), FALSE); write.csv(out, file.path(d, "tables/aggregate_tier2.csv"), row.names = FALSE)
print(format(out, digits = 3), row.names = FALSE)
