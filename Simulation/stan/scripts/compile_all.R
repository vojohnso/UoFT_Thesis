# Check the packages and compile every Stan model once, so parallel jobs don't all compile at the same time.
# Run it once from the Simulation folder before submitting jobs: Rscript stan/scripts/compile_all.R

for (pkg in c("rstan", "MASS", "spdep", "coda", "splines", "ggplot2")) {
  if (!requireNamespace(pkg, quietly = TRUE)) stop("missing R package: ", pkg)
}
source("stan/R/compile_model.R")

for (f in list.files("stan/models", pattern = "\\.stan$")) {
  name <- sub("\\.stan$", "", f)
  compile_model(name)
  cat("compiled", name, "\n")
}
