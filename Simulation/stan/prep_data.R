library(splines)

# Builds the data list for the Stan models.
# sim:      one dataset from data.generate()
# exposure: "oracle" uses the true exposure Ai, "naive" uses the surrogate Zi
# df:       degrees of freedom of the natural cubic spline for the exposure effect
# a_grid:   exposure values at which theta(a) is evaluated
#
# The spline basis is built here (not in Stan): knots sit at quantiles of the observed exposure, and the
# same basis is evaluated at the grid values, so "set the exposure to a for everyone" is a row of Bg 
make_stan_data <- function(sim, exposure = c("oracle", "naive"), df = 5,
                           a_grid = seq(5, 15, by = 1)) {
  exposure <- match.arg(exposure)
  a_obs <- if (exposure == "oracle") sim$Ai else sim$Zi
  basis <- ns(a_obs, df = df)
  B <- matrix(as.numeric(basis), ncol = df)
  Bg <- matrix(as.numeric(predict(basis, a_grid)), ncol = df)
  centre <- colMeans(B)
  list(
    N = nrow(sim),
    K = df,
    G = length(a_grid),
    Y = as.integer(sim$Yi),
    pop = as.numeric(sim$Ni),
    X = as.matrix(sim[, c("Xi1", "Xi2", "Xi3", "Xi4")]),
    B = sweep(B, 2, centre),
    Bg = sweep(Bg, 2, centre)
  )
}

# Latent exposure model: the true exposure is not given to Stan, only the surrogate Z.
# sigma_Z is the known error SD (the fitted model is told the error model)
# The exposure effect is a straight line in A, so no spline basis is needed (a spline basis would
# have to be rebuilt inside Stan from the latent A)
make_stan_data_latent <- function(sim, sigma_Z, a_grid = seq(5, 15, by = 1)) {
  list(
    N = nrow(sim),
    G = length(a_grid),
    Y = as.integer(sim$Yi),
    pop = as.numeric(sim$Ni),
    X = as.matrix(sim[, c("Xi1", "Xi2", "Xi3", "Xi4")]),
    Z = as.numeric(sim$Zi),
    sigma_Z = sigma_Z,
    a_grid = as.numeric(a_grid)
  )
}
