############################
# A Practical Introduction to Bayesian Estimation of Causal
# Effects: Parametric and Nonparametric Approaches
#############################

set.seed(1006092577)
library(BART)

############################
# Exercise 1: Standardization Practise
############################

# Based on Appendix E

# L ~ N(0, 1)
# A | L ~ Bern(sigma(1 - L/2))
# Y | A, L ~ N((L + 1/2L^2)A, 1/4)
# ATE = E(E(Y | A = 1, L) - E(Y | A = 0, L))
# = E(L+ 1/2L^2)
# = 1/2

# 1. Simulate data with true ATE (0.5)
# 2. Fit a regression
# 3. Predict everyone twice under A = 1 and A = 0
# 4. The ATE is the averaged differences 

n <- 10000

L <- rnorm(n, 0, 1)
A <- rbinom(n, 1, plogis(1-0.5*L))
Y <- rnorm(n, (L + 0.5*L^2)*A, 1/2)
TRUE_ATE = 1/2 

d <- data.frame(Y = Y, A = A, L = L)

fit <- lm(Y ~ A + L, d)
X1 <- predict(fit, newdata = transform(d, A = 1))
X0 <- predict(fit, newdata = transform(d, A = 0))
ATE = mean(X1 - X0)


############################
# Exercise 2: Parametric Bayesian Bootstrap
############################

# 1. Confounder W ~ N(0, 1)
# 2 . Stratum membership V with P(Vi = v) = p_v for v in 1, 2, ..., 5
# p_v = {0.3, 0.3, 0.2, 0.1, 0.1}
# 3. Treat assignment as Bernoulli P(A | W, V = v) = sigma(W + gamma), gamma = {0, -0.5, 0.5, 0.5, -0.5}
# 4. Bernoulli outcome P(Y | A, W, V) = logit[-1 + W + (1 + sum(n_v I(Vi = v)))A] where n_v = (-.5, 0, .5, .6)
library(rstanarm)
library(gtools)
n <- 500

gamma_v <- c(0, -.5, .5, .5, -.5) # stratum effect on treatment assignment
eta_v <- c(0, -.5, 0, .5, .6) # stratum effect on outcome (V = 1 is the reference group)g

W <- rnorm(n, 0, 1)
V <- sample(1:5, n, replace = TRUE, prob = c(0.3, 0.3, 0.2, 0.1, 0.1))
A <- rbinom(n, 1, plogis(W + gamma_v[V]))
Y <- rbinom(n, 1, plogis(-1 + W + (1  + eta_v[V]) * A))

d <- data.frame(Y, A, W, V = factor(V))

d <- d[order(d$V),]

fit <- stan_glmer(Y ~ W + V + A + (0 + A | V), # random slope on treatment by stratum with no random intercept 
                  data = d,
                  family = binomial, 
                  prior_covariance = decov(scale = 0.5))

M <- nrow(as.matrix(fit))

p1 <- posterior_epred(fit, newdata = transform(d, A = 1))
p0 <- posterior_epred(fit, newdata = transform(d, A = 0))

n_v <- as.vector(table(d$V))
psi_post <- matrix(NA, M, 5)
for(v in 1:5) {
  idx <- which(d$V == v) # all participants in stratum V
  for(m in 1:M){
    bb_w <- as.vector(rdirichlet(1, rep(1, length(idx)))) # Direchlet with length of the stratum
    m1 <- sum(bb_w * p1[m, idx])
    m0 <- sum(bb_w * p0[m, idx])
    psi_post[m, v] <- (m1/(1-m1) / (m0/(1-m0)))
  }
  
}

apply(psi_post, 2, quantile, c(0.025, 0.5, 0.975))

############################
# Exercise 3: Time-varying G-Computation
############################

# L0 -> A0 -> L1 -> A1 -> Y

alpha0 <- 0;    alpha1 <- 0.8     # A0 depends on L0
gamma0 <- 0;    gamma1 <- 0.6;  gamma2 <- -0.8    # treatment damages L1
beta0  <- 0;    beta1  <- 0.8;  beta2  <- 0.3     # L1 governs A1
delta0 <- 0;    delta1 <- 0.5
delta2 <- 0.5                                      # direct effect of A0
delta3 <- 1.0                                      # L1 helps outcome
delta4 <- 0.5                                      # effect of A1
sigma_L <- 1;   sigma_Y <- 1
n <- 2000

L0 <- rnorm(n, 0, 1)
A0 <- rbinom(n, 1, plogis(alpha0 + alpha1 * L0))
L1 <- rnorm(n, gamma0 + gamma1 * L0 + gamma2 * A0, sigma_L)
A1 <- rbinom(n, 1, plogis(beta0 + beta1 * L1 + beta2 * A0))
Y <- rnorm(n, delta0 + delta1 * L0 + delta2 * A0 + delta3 * L1 + delta4 * A1, sigma_Y)

d <- data.frame(Y, L0, A0, L1, A1)
TRUE_PSI <- delta2 + delta3*gamma2 + delta4

# Step 1: G-Formula requires every factor in the product estimated - outcome model, confounder models
# outcome model
fit_Y <- stan_glm(Y ~ L0 + A0 + L1 + A1, 
                  data = d,
                  family = gaussian(),
                  prior = normal(0, 0.25),
                  chains = 4,
                  iter = 2000)
# Step 2: confounder model
fit_L1 <- stan_glm(L1 ~ L0 + A0,
                     data = d,
                     family = gaussian(),
                     prior = normal(0, 0.25),
                     chains = 4,
                     iter = 2000)
# P(L0) is pre-treatment and hence the bayesian bootstrap should cover it.

# We want to compute mu(1, 1) for the causal estimand. But we can't just plug in A0 = A1 = 1
# Since we still have L1 values at observed values
# Instead, simulate L0 and L1.
M <- nrow(as.matrix(fit_Y))
post_Y <- as.matrix(fit_Y)
post_L1 <- as.matrix(fit_L1)

# Step 3: Simulate confounders sequentially for t in {0, ..., T}
# Repeat this B times to obtain the average confounders L0, L1, ...
B <- 5000
gcomp_one_draw <- function(m, a0, a1, B = 1000) {
  bb_w <- rdirichlet(1, rep(1, n))
  L0_sim <- sample(d$L0, B, replace = TRUE, prob = bb_w)
  
  # Simulate L1
  mu_L1 <- post_L1[m, "(Intercept)"] +
    post_L1[m, "L0"] * L0_sim +
    post_L1[m, "A0"] * a0
  L1_sim <- rnorm(B, mu_L1, post_L1[m, "sigma"])
  
  mu_Y <- post_Y[m, "(Intercept)"] +
    post_Y[m, "L0"] * L0_sim +
    post_Y[m, "A0"] * a0 +
    post_Y[m, "L1"] * L1_sim +
    post_Y[m, "A1"] * a1
  
  mean(mu_Y)        
}

# Step 4: Integrate the outcome model conditional on the current set of draws under both interventions 
psi <- numeric(M)
for (m in 1:M) {
  psi[m] <- gcomp_one_draw(m, a0 = 1, a1 = 1, B) - gcomp_one_draw(m, a0 = 0, a1 = 0, B) 
}
quantile(psi, c(0.025, 0.5, 0.975))
sd(psi)

m_fix <- 1
sapply(c(50, 100, 500, 1000, 5000, 20000), function(B) {
  replicate(20, gcomp_one_draw(m_fix, 1, 1, B) - gcomp_one_draw(m_fix, 0, 0, B))
}) |> apply(2, sd)
# B = 5000 seems enough

# Validate the covariate simulation
check_L1 <- function(m, a0, B = 20000) {
  w <- rdirichlet(1, rep(1, n))
  L0s <- sample(d$L0, B, TRUE, prob = w)
  mu  <- post_L1[m,"(Intercept)"] + post_L1[m,"L0"]*L0s + post_L1[m,"A0"]*a0
  mean(rnorm(B, mu, post_L1[m,"sigma"]))
}

c(treated = check_L1(1, 1), truth_t = gamma0 + gamma2,
  untreated = check_L1(1, 0), truth_u = gamma0)

data.frame(
  method = c("truth", "g-computation", "naive adjusted", "naive omitting L1"),
  value  = c(TRUE_PSI,
             median(psi),
             sum(coef(lm(Y ~ L0 + A0 + L1 + A1, d))[c("A0","A1")]),
             sum(coef(lm(Y ~ L0 + A0 + A1, d))[c("A0","A1")]))
)

############################
# Exercise 4: Dynamic Treatment Regimes
############################

# Similar to Exercise 3 in terms of time-varying confounders.
# Main difference is that confounders and interventions are dynamic and typically depend on a decision rule

# Rule: r(Lt) = I(L_t > kappa)

# Use the same simulation data as Exercise 3

alpha0 <- 0;    alpha1 <- 0.8     # A0 depends on L0
gamma0 <- 0;    gamma1 <- 0.6;  gamma2 <- -0.8    # treatment damages L1
beta0  <- 0;    beta1  <- 0.8;  beta2  <- 0.3     # L1 governs A1
delta0 <- 0;    delta1 <- 0.5
delta2 <- 0.5                                      # direct effect of A0
delta3 <- 1.0                                      # L1 helps outcome
delta4 <- 0.5                                      # effect of A1
sigma_L <- 1;   sigma_Y <- 1
n <- 2000

L0 <- rnorm(n, 0, 1)
A0 <- rbinom(n, 1, plogis(alpha0 + alpha1 * L0))
L1 <- rnorm(n, gamma0 + gamma1 * L0 + gamma2 * A0, sigma_L)
A1 <- rbinom(n, 1, plogis(beta0 + beta1 * L1 + beta2 * A0))
Y <- rnorm(n, delta0 + delta1 * L0 + delta2 * A0 + delta3 * L1 + delta4 * A1, sigma_Y)

d <- data.frame(Y, L0, A0, L1, A1)
TRUE_PSI <- delta2 + delta3*gamma2 + delta4

# Step 1: Confounder model: remains the same
# outcome model
fit_Y <- stan_glm(Y ~ L0 + A0 + L1 + A1, 
                  data = d,
                  family = gaussian(),
                  prior = normal(0, 0.25),
                  chains = 4,
                  iter = 2000)
# Step 2: confounder model
fit_L1 <- stan_glm(L1 ~ L0 + A0,
                   data = d,
                   family = gaussian(),
                   prior = normal(0, 0.25),
                   chains = 4,
                   iter = 2000)

M <- nrow(as.matrix(fit_Y))
post_Y <- as.matrix(fit_Y)
post_L1 <- as.matrix(fit_L1)

gcomp_dynamic <- function(m, kappa, B = 1000) {
  bb_w <- rdirichlet(1, rep(1, n))
  L0_sim <- sample(d$L0, B, replace = TRUE, prob = bb_w)
  
  # Apply decision rule
  a0 <- as.numeric(L0_sim > kappa)
  # Simulate L1
  mu_L1 <- post_L1[m, "(Intercept)"] +
    post_L1[m, "L0"] * L0_sim +
    post_L1[m, "A0"] * a0
  L1_sim <- rnorm(B, mu_L1, post_L1[m, "sigma"])
  
  # Apply decision rule 
  a1 <- as.numeric(L1_sim > kappa)
  mu_Y <- post_Y[m, "(Intercept)"] +
    post_Y[m, "L0"] * L0_sim +
    post_Y[m, "A0"] * a0 +
    post_Y[m, "L1"] * L1_sim +
    post_Y[m, "A1"] * a1
  
  mean(mu_Y)        
}
kappas <- seq(-2, 2, by = 0.25)

mu_k <- sapply(kappas, function(k)
  vapply(1:M, function(m) gcomp_dynamic(m, k), numeric(1)))

curve <- t(apply(mu_k, 2, quantile, c(0.025, 0.5, 0.975)))

plot(kappas, curve[,2], type = "o", pch = 19, ylim = range(curve),
     xlab = expression(kappa), ylab = expression(mu(kappa)),
     main = "Mean outcome under threshold rules")
polygon(c(kappas, rev(kappas)), c(curve[,1], rev(curve[,3])),
        col = adjustcolor("steelblue", 0.25), border = NA)
lines(kappas, curve[,2], type = "o", pch = 19, col = "steelblue")

############################
# Exercise 5: Extend to T Points
############################

Nt <- 10 # Number of time points
N <- 300 # Number of subjects
beta_L_true <- c(0.8, 0, 0, 0.25, -0.25, 0.5, -0.5, 1.0, -1.0) # True confounder effects
theta_true  <- -0.3 # Effect of A on current L 

L <- A <- matrix(NA, N, Nt) # Matrix of N x Nt to hold all simulated confounder effects

L[, 1] <- rnorm(N, 0, 1) # Simulate L0 for all N
A[, 1] <- rbinom(N, 1, plogis(L[, 1])) # Simulate A0

# Simulate all confounders, treatments and outcome
for (t in 2:Nt) {
  h <- 1:(t-1)
  lag <- (t-1):1
  mu <- L[, h, drop=FALSE] %*% beta_L_true[lag] + A[, h, drop=FALSE] %*% rep(theta_true, t-1)
  L[, t] <- rnorm(N, mu, sd = 2)
  A[, t] <- rbinom(N, 1, plogis(L[, t]))
}

Y <- rnorm(N, A %*% rep(0.5, Nt) + L %*% rep(0.4, Nt), 1)

# Compute the true causal estimand
true_mu <- function(a_vec, B = 200000) {
  Ls <- matrix(NA, B, Nt)
  Ls[,1] <- rnorm(B, 0, 1)
  for (t in 2:Nt) {
    h <- 1:(t-1)
    lag <- (t-1):1
    mu <- Ls[, h, drop=FALSE] %*% beta_L_true[lag] + sum(a_vec[h]) * theta_true
    Ls[,t] <- rnorm(B, mu, sd = 2)
  }
  mean(Ls %*% rep(0.4, Nt) + sum(a_vec) * 0.5)
}

TRUE_PSI <- true_mu(rep(1, Nt)) - true_mu(rep(0, Nt))

# Causal Step

fit_Y <- fit_L <- matrix(NA, N, Nt)
make_df <- function(t) {
  h <- 1:(t-1)
  df <- data.frame(Lt = L[,t], L[, h, drop=FALSE], A[, h, drop=FALSE])
  names(df) <- c("Lt", paste0("L", h), paste0("A", h))
  df
}

fits_L <- lapply(2:Nt, function(t)
  stan_glm(Lt ~ ., data = make_df(t), family = gaussian(),
           prior = normal(0, 2), chains = 2, iter = 1000, seed = 1, refresh = 0))
names(fits_L) <- paste0("t", 2:Nt)

df_Y <- data.frame(Y = Y, L, A)
names(df_Y) <- c("Y", paste0("L", 1:Nt), paste0("A", 1:Nt))
fit_Y <- stan_glm(Y ~ ., data = df_Y, family = gaussian(),
                  prior = normal(0, 2), chains = 2, iter = 1000, seed = 1, refresh = 0)

post_L <- lapply(fits_L, as.matrix)
post_Y <- as.matrix(fit_Y)
M <- nrow(post_Y)

gcomp_T <- function(m, a_vec, B = 500) {
  
  Ls <- matrix(NA, B, Nt)
  
  ## t = 1: BB-resample the observed baseline (pre-treatment)
  w <- as.vector(rdirichlet(1, rep(1, N)))
  Ls[,1] <- sample(L[,1], B, replace = TRUE, prob = w)
  
  ## t = 2..Nt: simulate forward under the regime
  for (t in 2:Nt) {
    p <- post_L[[t-1]][m, ]                      # draws for THIS time point's model
    h <- 1:(t-1)
    
    mu <- p["(Intercept)"] +
      Ls[, h, drop=FALSE] %*% p[paste0("L", h)] +
      sum(a_vec[h] * p[paste0("A", h)])       # a_vec is fixed, so this is a scalar
    
    Ls[,t] <- rnorm(B, mu, p["sigma"])
  }
  
  ## outcome model on the simulated histories
  q  <- post_Y[m, ]
  mu_Y <- q["(Intercept)"] +
    Ls %*% q[paste0("L", 1:Nt)] +
    sum(a_vec * q[paste0("A", 1:Nt)])
  
  mean(mu_Y)
}
psi <- vapply(1:M, function(m)
  gcomp_T(m, rep(1, Nt)) - gcomp_T(m, rep(0, Nt)), numeric(1))

quantile(psi, c(0.025, 0.5, 0.975))
TRUE_PSI

# Exercise 6: BART

n <- 2000

L <- rnorm(n, 0, 1)
A <- rbinom(n, 1, plogis(1-0.5*L))
Y <- rnorm(n, (L + 0.5 * L^2) * A, 1/2)
TRUE_ATE = 1/2 

d <- data.frame(Y = Y, A = A, L = L)
X_train <- d[, c("A", "L")]
Y_train <- Y

X_a1 <- X_a0 <- X_train
X_a1[, "A"] <- 1
X_a0[, "A"] <- 0
X_test <- rbind(X_a1, X_a0)

fit_bart <- gbart(X_train, Y_train, x.test = X_test, type = "wbart", ndpost = 1000)
mu_a1 <- fit_bart$yhat.test[, 1:n]          # 1000 x 2000 matrix (M x n) 
mu_a0 <- fit_bart$yhat.test[, (n+1):(2*n)] 


bayes_boot <- function(m1, m0) {
  M <- nrow(m1); n <- ncol(m1)
  psi <- rep(NA, M)
  for (m in 1:M) {
    bb_weights <- rdirichlet(1, rep(1, n))
    psi[m] <- sum(bb_weights * (m1[m, ] - m0[m, ]))
  }
  psi
}

psi_bart <- bayes_boot(mu_a1, mu_a0)
quantile(psi_bart, c(0.025, 0.5, 0.975))

############################
# Kuan Liu's Causal Workshop
############################

library(tidyverse)
sas_origin <- as.Date("1960-01-01")

rhc <- read_csv("https://hbiostat.org/data/repo/rhc.csv") |>
  mutate(
    sadmdte = as.Date(sadmdte, origin = sas_origin),
    dschdte = as.Date(dschdte, origin = sas_origin),
    dthdte  = as.Date(dthdte,  origin = sas_origin),
    lstctdte = as.Date(lstctdte, origin = sas_origin),
    A       = as.integer(swang1 == "RHC"),
    Y_death = as.integer(death == "Yes"),
    Y_los   = as.numeric(coalesce(dschdte, lstctdte) - sadmdte)
  )

rhc |> count(A, swang1)
############################
# Parametric Bayesian Causal
############################

# Step 1: Specify a Bayesian outcome model
# Step 2: Obtain posterior draws via MCMC (posterior draws of the model parameters)
# Step 3: For each draw, compute counterfactual predictions
# Step 4: Integrate over the Bayesian Bootstrap using Dirichlet(1_n) weights 

# Outcome = death (Y_death)
# Treatment = RHC (A)
# Confounders: age + sex + race + cat1 + meanbp1 + hrt1 + resp1 + temp1 + wtkilo1
# Note the outcome is binary - specify the outcome model as a binomial model
# Step 1/2:
outcome_fit <- stan_glm(Y_death ~ A + age + sex + race + cat1 + meanbp1 + hrt1 + resp1 + temp1 + wtkilo1,
                        data = rhc,
                        family = binomial(link="logit"),
                        prior = normal(0, 0.25),
                        prior_intercept = normal(0, 5),
                        chains = 4,
                        iter = 1000,
                        seed = 1006092577)

# Step 3:
rhc_m1 <- rhc_m0 <- rhc
rhc_m1$A <- 1
rhc_m0$A <- 0

m1 <- posterior_epred(outcome_fit, newdata = rhc_m1)
m0 <- posterior_epred(outcome_fit, newdata = rhc_m0)

# Step 4:
M <- nrow(as.matrix(outcome_fit))
N <- nrow(rhc)
psi <- numeric(M)
for (i in 1:M) {
  w <- rdirichlet(1, rep(1, N))
  psi[i] <- sum(w * (m1[i, ] - m0[i, ]))
}
quantile(psi, c(0.025, 0.5, 0.975))


ggplot(data.frame(psi = psi), aes(x = psi)) +
  geom_density(fill = "steelblue", alpha = 0.4) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(x = expression(Psi ~ "(Risk Difference)"),
       title = "Posterior ATE: Effect of RHC on 30-day Mortality") +
  theme_minimal()

############################
# Non-parametric G-Computation 
############################

# Instead of specifying an outcome model we can use a BART approach which avoids misspecification
# Use gbart function which requires train and test sets for the model
# Use type="pbart" for probit

covars  <- c("age","sex","race","cat1","meanbp1",
             "hrt1","resp1","temp1","wtkilo1")
X_train <- model.matrix(~ A + ., data = rhc[, c("A", covars)])[, -1]
Y_train <- rhc$Y_death
X_a1 <- X_a0 <- X_train
X_a1[, "A"] <- 1
X_a0[, "A"] <- 0
# Test set is them the counterfactuals
X_test <- rbind(X_a1, X_a0)

fit_bart <- gbart(X_train, Y_train, x.test = X_test, type = "pbart", ndpost = 1000) # M = 1000 number of draws

# M x 2N matrix, the first N are A = 1 and the second N are A = 0
# 
mu_a1 <- pnorm(fit_bart$yhat.test[, 1:N])      
mu_a0 <- pnorm(fit_bart$yhat.test[, (N+1):(2*N)])
# Once again we run the Bayesian bootstrap to get our causal estimand
M <- nrow(fit_bart$yhat.test)

bayes_boot <- function(mu_a1, mu_a0) {
  M <- nrow(mu_a1)
  N <- ncol(mu_a1)
  sapply(1:M, function(m) {
    w <- rdirichlet(1, rep(1, N))
    sum(w * (mu_a1[m, ] - mu_a0[m, ]))
  })
}

psi_bart <- bayes_boot(mu_a1, mu_a0)
quantile(psi_bart, c(0.025, 0.5, 0.975))
# Looks pretty similar

psi_comparison <- data.frame(
  psi = c(psi, psi_bart),
  method = rep(c("Parametric", "BART"), times = c(length(psi), length(psi_bart)))
)

ggplot(psi_comparison, aes(x = psi, fill = method)) +
  geom_density(alpha = 0.4) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
  scale_fill_manual(values = c("Parametric" = "tomato", "BART" = "steelblue")) +
  labs(
    x = expression(Psi ~ "(Risk Difference)"),
    title = "Posterior ATE: Effect of RHC on 30-day Mortality",
    subtitle = "Parametric vs BART g-computation",
    fill = "Method"
  ) +
  theme_minimal()

############################
# Bayesian Propensity Score Weighting 
############################

# Propensity score models the treatment given covariates
# Propogate uncertainty through Bayesian decision theory
# "We want to draw inference from an ideal experimental framework population (unfounded)
# but we observe data from a population where A depends on L (confounded)
# Utility function as the log-likelihood of the marginal outcome model log P(Y_i | A_i)
# Maximize the utlity function under the experimental distribution but integrate over the observed one



