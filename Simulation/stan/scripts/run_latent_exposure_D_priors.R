# Latent exposure model with estimated error parameters (beta_Z1, beta_Z2, sigma_Z) and the CAR spatial effect
# (gcomp_latent_exposure_D_priors.stan), one dataset, with diagnostics: summary table, trace plots, joint-posterior
# plots, and per-chain means 

library(rstan)
library(ggplot2)
source("gen_data.R"); source("stan/R/prep_data.R"); source("stan/R/compile_model.R")

n_areas <- 100
a_grid <- seq(5, 15, by = 1)
sim <- data.generate(seed = 1006092577, n_areas = n_areas, out_form = "linear", bZ1 = 1, bZ2 = 0.5, sigma_Z = 1)
dat <- make_stan_data_latent_D_priors(sim, a_grid)
truth <- c(beta_Z1 = 1, beta_Z2 = 0.5, sigma_Z = 1, sigma_A = 1.5, sigma_phi = sqrt(sim$sigma2_phi[1]), rho = 0.9, beta_A = 0.10)
cat("True values:", paste(names(truth), round(truth, 2), collapse = ", "), "\n")

init <- function() list(beta0 = -3, beta_A = 0, beta_X = rep(0, 4), gamma0 = mean(dat$Z), gamma_A = rep(0, 4), sigma_A = 1,
                        A = dat$Z, beta_Z1 = 1, beta_Z2 = 0, sigma_Z = 1,
                        phi = rep(0, dat$N), sigma_phi = 1, rho = 0.5)
fit <- sampling(compile_model("gcomp_latent_exposure_D_priors"), data = dat, chains = 4, iter = 2000, cores = 4,
                seed = 1006092577, refresh = 0, init = init, control = list(adapt_delta = 0.9))

pars <- c("beta_Z1", "sigma_A", "sigma_phi", "rho", "beta_A", "beta_Z2", "sigma_Z")
sp <- get_sampler_params(fit, inc_warmup = FALSE)
cat("\nDivergences:", sum(sapply(sp, function(m) sum(m[, "divergent__"]))),
    "  hit max tree depth:", sum(sapply(sp, function(m) sum(m[, "treedepth__"] >= 10))), "\n")
print(round(summary(fit, pars = pars)$summary[, c("mean", "sd", "n_eff", "Rhat")], 3))

# draws: iteration x chain x parameter
draws <- rstan::extract(fit, pars = pars, permuted = FALSE)
long <- do.call(rbind, lapply(pars, function(p) data.frame(iter = rep(seq_len(dim(draws)[1]), dim(draws)[2]),
        chain = factor(rep(seq_len(dim(draws)[2]), each = dim(draws)[1])), param = p, value = as.vector(draws[, , p]))))
long$param <- factor(long$param, levels = pars)

# per-chain means: if chains disagree, the posterior is not being explored consistently
wide <- data.frame(chain = factor(rep(seq_len(dim(draws)[2]), each = dim(draws)[1])), apply(draws, 3, as.vector))
chain_means <- aggregate(. ~ chain, wide, mean)
chain_means$sigmaA_x_betaZ1 <- with(wide, tapply(sigma_A * beta_Z1, chain, mean))
chain_means$betaA_div_betaZ1 <- with(wide, tapply(beta_A / beta_Z1, chain, mean))
cat("\nPosterior means by chain (truth: beta_Z1 = 1, sigma_A = 1.5, beta_A = 0.1: \n")
print(round(chain_means[, -1], 3), row.names = chain_means$chain)

dir.create("output/stan", recursive = TRUE, showWarnings = FALSE)
tr <- ggplot(long, aes(iter, value, colour = chain)) + geom_line(alpha = 0.7, linewidth = 0.3) +
  geom_hline(data = data.frame(param = factor(names(truth), levels = pars), v = truth), aes(yintercept = v),
             linetype = "dashed", colour = "black") +
  facet_grid(param ~ ., scales = "free_y") + labs(x = "iteration (post-warmup)", y = NULL, colour = "chain",
  title = "Trace Plots: latent model with estimated error parameters + CAR", subtitle = "Dashed line = true value") +
  theme_bw(base_size = 11)
ggsave("output/stan/latent_exposure_D_priors_trace.png", tr, width = 9, height = 11.5, dpi = 120)

pr <- rbind(data.frame(pair = "beta_Z1 vs sigma_A", x = wide$beta_Z1, y = wide$sigma_A, chain = wide$chain),
            data.frame(pair = "beta_Z1 vs beta_A", x = wide$beta_Z1, y = wide$beta_A, chain = wide$chain),
            data.frame(pair = "sigma_A vs sigma_phi", x = wide$sigma_A, y = wide$sigma_phi, chain = wide$chain),
            data.frame(pair = "sigma_phi vs rho", x = wide$sigma_phi, y = wide$rho, chain = wide$chain))
pa <- ggplot(pr, aes(x, y, colour = chain)) + geom_point(alpha = 0.35, size = 0.7) + facet_wrap(~pair, scales = "free") +
  labs(x = NULL, y = NULL, title = "Joint posterior draws by chain") + theme_bw(base_size = 11)
ggsave("output/stan/latent_exposure_D_priors_pairs.png", pa, width = 15, height = 4.2, dpi = 120)
cat("\nSaved output/stan/latent_exposure_D_priors_trace.png and latent_exposure_D_priors_pairs.png\n")
