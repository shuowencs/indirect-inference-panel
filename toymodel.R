# 1/26/2021 Shuowen Chen
# This script runs indirect inference bias correction for a toy model


# Packages
library(plm)
library(doParallel)
library(doRNG)
library(alpaca)
library(pracma) # for Jacobian

########### 1.  Set Parameters #############
N <- 200 # cross section
T <- 5 # time series
mu <- 0.5
sigma <- 1.5

H <- 5 # number of simulated samples
S <- 10 # number of Monte Carlo simulations

########### 2.  Define Functions #############
# Simulated panel only individual fixed effects
# y_it = 1(alpha_i +  e_it > 0)
sim_dat <- function(shocks, unitfe){
  T <- dim(shocks)[1] # time series dimension
  N <- dim(shocks)[2] # cross section dimension
  
  # simulate AR(1) panel sequence
  y <- matrix(0, nrow = T, ncol = N)
  for (t in 1:T) y[t, ] <- unitfe > shocks[t, ] 
  
  # construct panel data set
  y <- matrix(y, ncol = 1)
  dat <- data.frame(id = kronecker(c(1:N), rep(1, T)), 
                    time = kronecker(rep(1, N), c(1:T)), y = y)
  dat <- pdata.frame(dat, index = c("id", "time"))
  return(dat)
}

# objective function for indirect inference
# inputs:
#  coef: fixed effect coefficients
objective <- function(phi, shocks_sim, alpha_i, coef) {
  # compute simulated beta, returns |beta-beta_sim|
  H <- dim(shocks_sim)[3]
  alpha_i <- phi[1] + phi[2]*alpha_i
  print(phi)
  # loop over sample paths
  coef_sim <- rep(0, 2)
  for (h in 1:H) {
    dats <- sim_dat(shocks_sim[, , h], alpha_i[, h])
    fit_sim <- speedglm(y ~ factor(id) - 1, data = dats, family = binomial(link = "probit"))
    coef_sim <- coef_sim + c(mean(coef(fit_sim)), sd(coef(fit_sim)))/H
  }
  # Use identity weighting matrix
  return( sum( (coef - coef_sim)^2 ) )
}


# a function to obtain indirect inference estimator for one simulation
# inputs:
# N: number of cross section units
# T: number of time series 

estimate_ii <- function(N, T, H, mu0, sigma0) {
  # generate panel data
  alpha <- rnorm(N, 0, 1) # individual fixed effect
  e <- matrix(rnorm(N*T, 0, 1), nrow = T, ncol = N) # panel error
  
  y <- matrix(0, nrow = T, ncol = N)
  for (t in 1:T) y[t, ] <- mu0 + sigma0*alpha > e[t, ]
  # construct panel data set
  y <- matrix(y, ncol = 1)
  dat <- data.frame(id = kronecker(c(1:N), rep(1, T)), 
                    time = kronecker(rep(1, N), c(1:T)), y = y)
  dat <- pdata.frame(dat, index = c("id", "time"))
  # fixed effect estimate
  fit <- speedglm(y ~ factor(id) - 1, data = dat, family = binomial(link = "probit"))
  # FE estimate (auxiliary estimator)
  mu_hat <- mean(coef(fit))
  sigma_hat <- sd(coef(fit))
  coef_fe <- c(mu_hat, sigma_hat)
  
  # shocks and unit fe for sim data in indirect inference
  shocks_sim <- array(rnorm(N*T*H, 0, 1), c(T, N, H))
  al_sim <- matrix(rnorm(N*H, 0, 1), nrow = N, ncol = H)
  
  # estimation
  est <- optim(c(mu0, sigma0), objective, method = "Nelder-Mead",
               shocks_sim = shocks_sim, alpha_i = al_sim,
               coef = coef_fe)
  results <- est$par
  #est <- optimize(objective, c(-1, 1), shocks_sim = shocks_sim, 
  #                alpha_i = alpha_sim, coef = coef_fe)
  
  # return results
  return(list(ii = results, auxiliary = coef_fe))
}

########### 3. Simulation #############
set.seed(88)

nCores <- 1   # number of CPUs for parallelization
registerDoParallel(cores = nCores)

# loop over number of simulations
results_par <- foreach(s = 1:S) %dorng% {
  results <- estimate_ii(N, T, H, mu, sigma)
}

# convert results to a matrix
ii_matrix <- sapply(results_par, function(x) return(x$ii))
fe_matrix <- sapply(results_par, function(x) return(x$auxiliary))







