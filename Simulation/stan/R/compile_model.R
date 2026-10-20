# Compiles stan/models/<name>.stan and caches the compiled model in stan/compiled/<name>.rds, so the compiled
# files stay out of the models folder. The model is recompiled only when the .stan file is newer than the cache.
# Usage (from the Simulation/ folder):  mod <- compile_model("gcomp_known_exposure");  fit <- rstan::sampling(mod, data = ...)
compile_model <- function(name) {
  src <- file.path("stan", "models", paste0(name, ".stan"))
  rds <- file.path("stan", "compiled", paste0(name, ".rds"))
  if (!file.exists(src)) stop("no such Stan model: ", src)
  dir.create(dirname(rds), showWarnings = FALSE, recursive = TRUE)
  if (file.exists(rds) && file.mtime(rds) > file.mtime(src)) return(readRDS(rds))
  mod <- rstan::stan_model(src)
  saveRDS(mod, rds)
  mod
}
