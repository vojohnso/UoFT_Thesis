library(parallel)

source("gen_data.R")
source("estimators.R")

# Replicate r uses seed base_seed + r, so each replicate is reproducible on its own and the
# result does not depend on the number of cores. The three scenarios use the same seeds, so they
# share covariates, exposure, error and populations and differ only in the outcome shape.
# Data settings (areas, spatial effect, error model) are the defaults in data.generate().
args <- commandArgs(trailingOnly = TRUE)
n_reps <- if (length(args) >= 1) as.integer(args[1]) else 200
n_cores <- if (length(args) >= 2) as.integer(args[2]) else max(1, detectCores() - 1)
base_seed <- 1006092577

output_dir <- "output"
dir.create(output_dir, showWarnings = FALSE)

a_grid <- seq(5, 15, by = 1)
forms <- c("linear", "quadratic", "nonlinear")
methods <- c("naive", "oracle")

# Population truth, fixed for each scenario 
pop_truth <- lapply(setNames(forms, forms), function(f) true_erf_pop(a_grid, f))

# One replicate: summary table (mean and 95% interval at each grid value) and the posterior draws 
run_one <- function(rep_id, form) {
  sim <- data.generate(seed = base_seed + rep_id, out_form = form)
  true_sample <- true_erf(a_grid, sim)
  fits <- list()
  invisible(capture.output({
    fits$naive <- bart_gcomp(sim, a_grid = a_grid, use_true_A = FALSE)
    fits$oracle <- bart_gcomp(sim, a_grid = a_grid, use_true_A = TRUE)
  }))
  summary <- do.call(rbind, lapply(methods, function(m) {
    s <- fits[[m]]$theta_summary
    data.frame(form = form, rep = rep_id, method = m, a = a_grid, true = pop_truth[[form]],
               true_sample = true_sample, est = s$mean, lo = s$q025, hi = s$q975)
  }))
  draws <- sapply(methods, function(m) fits[[m]]$theta_draws, simplify = "array")
  list(summary = summary, draws = draws)   # draws: posterior draw x grid value x method
}

for (form in forms) {
  t0 <- Sys.time()
  cat(form, ":", n_reps, "replicates on", n_cores, "cores\n")
  res <- mclapply(seq_len(n_reps), run_one, form = form, mc.cores = n_cores)
  ok <- !sapply(res, inherits, "try-error")
  if (!all(ok)) cat("  WARNING:", sum(!ok), "replicates failed\n")
  res <- res[ok]
  out <- list(
    summary = do.call(rbind, lapply(res, `[[`, "summary")),
    rep = which(ok),
    draws = sapply(res, `[[`, "draws", simplify = "array")   # draw x grid x method x replicate
  )
  saveRDS(out, file.path(output_dir, sprintf("replicates_%s.rds", form)))
  cat("  done in", round(difftime(Sys.time(), t0, units = "mins"), 1), "min\n")
}
