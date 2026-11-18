# Find the number of iterations and the burn-in (warm-up) needed for convergence.
#
# For ONE simulated replicate, every method is refitted over a ladder of (iterations, warm-up)
# settings, with the same data, seed and chains. For each fit it records the convergence
# diagnostics (R-hat, Geweke, divergences), the run time, and the posterior mean of theta(a) at
# a = 5, 10, 15 (to see when the estimate stops moving). The smallest setting at which the
# diagnostics are acceptable for every method is the one to use in the main study.

library(rstan)
library(ggplot2)
source("gen_data.R")
source("stan/R/prep_data.R")
source("stan/R/compile_model.R")
source("stan/R/diagnostics.R")


# Settings
# Arguments: replicate number, then optionally which methods (e.g. wild,0.3) and the iter:warmup ladder
args <- commandArgs(trailingOnly = TRUE)
r <- if (length(args) >= 1) as.integer(args[1]) else 1L
method_filter <- if (length(args) >= 2) args[2] else "all"
ladder_txt <- if (length(args) >= 3) args[3] else "2000:1000,4000:2000,8000:4000,16000:8000"

ladder <- do.call(rbind, lapply(strsplit(strsplit(ladder_txt, ",")[[1]], ":"), function(x) {
  data.frame(iter = as.integer(x[1]), warmup = as.integer(x[2]))
}))

base_seed <- 1006092577
n_areas <- 100
a_grid <- seq(5, 15, by = 1)
out_dir <- "output/stan/tuning"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Same seed scheme as sim_main.R, so replicate r is the same dataset
set.seed(base_seed)
s <- sample.int(.Machine$integer.max, 100)[r]


# Data and models
sim <- data.generate(
  seed = s, n_areas = n_areas, out_form = "linear", spatial_in = "outcome",
  bZ1 = 1, bZ2 = 0.5, sigma_Z = 1
)
car <- car_structure(nrow(sim))
d_known <- c(make_stan_data_latent_D(sim, bZ1 = 1, bZ2 = 0.5, sigma_Z = 1, a_grid = a_grid), car)
d_red <- make_stan_data_latent_D_priors(sim, a_grid)

exposure_mod <- compile_model("gcomp_known_exposure_car")
known_mod <- compile_model("gcomp_latent_exposure_D_known_car")
red_mod <- compile_model("gcomp_latent_exposure_D_reduced_sparse")
tight_mod <- compile_model("gcomp_latent_exposure_D_reduced_sparse_tight")

# Starting values: A starts at the observed surrogate Z
init_known <- function() {
  list(
    beta0 = -3, beta_A = 0, beta_X = rep(0, 4),
    gamma0 = mean(sim$Zi), gamma_A = rep(0, 4), sigma_A = 1,
    A = sim$Zi,
    phi = rep(0, nrow(sim)), sigma_phi = 1, rho = 0.5
  )
}
init_red <- function() c(init_known(), list(beta_Z2 = 0, sigma_Z = 1))

# Each method: model, data and starting values
methods <- list(
  "oracle" = list(exposure_mod, c(make_stan_data(sim, "oracle", df = 1, a_grid = a_grid), car), "random"),
  "naive" = list(exposure_mod, c(make_stan_data(sim, "naive", df = 1, a_grid = a_grid), car), "random"),
  "latent, known errors" = list(known_mod, d_known, init_known),
  "reduced, sigma_Z fixed at 1" = list(red_mod, c(d_red, list(sigma_Z_mean = 1, sigma_Z_sd = 0.001)), init_red),
  "reduced, sigma_Z ~ normal(1, 0.3)" = list(red_mod, c(d_red, list(sigma_Z_mean = 1, sigma_Z_sd = 0.3)), init_red),
  "reduced, sigma_Z wild guess" = list(red_mod, c(d_red, list(sigma_Z_mean = 0, sigma_Z_sd = 5)), init_red),
  "reduced, tight priors on everything" = list(tight_mod, d_red, init_red)
)
if (method_filter != "all") {
  keep <- Reduce(`|`, lapply(strsplit(method_filter, ",")[[1]], function(p) grepl(p, names(methods), fixed = TRUE)))
  methods <- methods[keep]
}


# Fit every method at every setting
rows <- list()
for (m in names(methods)) {
  for (k in seq_len(nrow(ladder))) {
    iter <- ladder$iter[k]
    warmup <- ladder$warmup[k]
    spec <- methods[[m]]

    t0 <- Sys.time()
    fit <- sampling(
      spec[[1]], data = spec[[2]], chains = 4, iter = iter, warmup = warmup, cores = 4,
      seed = s, refresh = 0, init = spec[[3]], control = list(adapt_delta = 0.9)
    )
    minutes <- as.numeric(difftime(Sys.time(), t0, units = "mins"))

    conv <- convergence_summary(fit)
    theta <- colMeans(rstan::extract(fit, "theta")$theta)
    rows[[length(rows) + 1]] <- data.frame(
      rep = r, method = m, iter = iter, warmup = warmup, minutes = round(minutes, 1),
      divergences = conv$div,
      rhat_key = round(conv$rhat_key, 3), worst_param = conv$worst_param,
      rhat_all = round(conv$max_rhat, 3),
      geweke_frac_flag = round(conv$geweke_frac_flag, 3),
      geweke_max_abs_z = round(conv$geweke_max_abs_z, 2),
      theta_a5 = round(theta[a_grid == 5]), theta_a10 = round(theta[a_grid == 10]),
      theta_a15 = round(theta[a_grid == 15])
    )

    # Trace plots of the worst parameters, at the largest setting only
    if (k == nrow(ladder)) {
      sm <- summary(fit)$summary
      key <- grep(paste0("^(", paste(intersect(key_pars, fit@model_pars), collapse = "|"), ")(\\[|$)"), rownames(sm), value = TRUE)
      worst <- head(key[order(-sm[key, "Rhat"])], 4)
      tag <- gsub("[^A-Za-z0-9]+", "_", m)
      ggsave(
        file.path(out_dir, sprintf("trace_rep%d_%s_%d.png", r, tag, iter)),
        rstan::traceplot(fit, pars = worst, inc_warmup = FALSE) + theme(legend.position = "bottom"),
        width = 9, height = 6, dpi = 110
      )
    }

    res <- do.call(rbind, rows)
    write.csv(res, file.path(out_dir, sprintf("tuning_rep%d.csv", r)), row.names = FALSE)
    saveRDS(res, file.path(out_dir, sprintf("tuning_rep%d.rds", r)))
    cat(sprintf("%-40s iter %5d warmup %5d  R-hat(key) %.3f  Geweke flag %.2f  div %d  %.1f min\n",
                m, iter, warmup, conv$rhat_key, conv$geweke_frac_flag, conv$div, minutes))
  }
}
