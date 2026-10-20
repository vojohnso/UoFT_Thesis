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

# Adjacency matrix for the CAR spatial effect on the exposure. Areas are rows of the 20 x 20 rook grid used by
# data.generate(), in the same order (area_id). Stan gets the dense matrix (and the number of neighbouring pairs, which
# only a sparse CAR implementation uses).
# rho (CAR dependence) is treated as known; the simulation uses rho_A = 0.9.
car_structure <- function(n_areas, rho = 0.9) {
  W <- make_adj_matrix(n_areas)
  list(W = W, W_n = as.integer(sum(W) / 2), rho = rho)
}

# Latent exposure model with the Option A error model Z = bZ1 * A + bZ2 * D + noise.
# Known version: bZ1, bZ2 and sigma_Z are data (the fitted model is told the error model); no spatial term.
make_stan_data_latent_D <- function(sim, bZ1 = 1, bZ2 = 0.5, sigma_Z = 1, a_grid = seq(5, 15, by = 1)) {
  c(make_stan_data_latent(sim, sigma_Z, a_grid), list(D = as.numeric(sim$Di), beta_Z1 = bZ1, beta_Z2 = bZ2))
}

# Estimated version: beta_Z1, beta_Z2 and sigma_Z are parameters with priors written in the Stan file
# (beta_Z1 sets the scale of A and is not identified by the data; see the identifiability discussion in sim_versions.tex).
# The model also has the CAR spatial effect on the exposure, so the neighbour structure is included (rho known).
make_stan_data_latent_D_priors <- function(sim, a_grid = seq(5, 15, by = 1), rho = 0.9) {
  d <- make_stan_data_latent(sim, sigma_Z = 1, a_grid = a_grid)
  d$sigma_Z <- NULL
  c(d, car_structure(nrow(sim), rho), list(D = as.numeric(sim$Di)))
}
