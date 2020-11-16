# 11/4/2020 Shuowen Chen
# R Code for Indirect inference for nonlinear panel models

# This script replicates the approach in Gourieoux, Phillips and Yu (2010)

# Housekeeping
rm(list = ls())

# Packages
library(plm)
library(doParallel)
library(doRNG)

########### 1.  Set Parameters #############
N <- 100 # cross section
T <- 5 # time series
phi0 <- 0.95 # AR 1 coefficient
H <- 5 # number of simulated samples
S <- 1 # number of Monte Carlo simulations

# placeholder for indirect inference estimator
estimates <- rep(NA, S)

########### 2.  Define Functions #############
# Simulated panel AR(1) with only individual fixed effects
# y_it = alpha_i + phi*y_it-1 + e_it
sim_dat <- function(phi, shocks, unitfe){
  T <- dim(shocks)[1] # time series dimension
  N <- dim(shocks)[2] # cross section dimension
  
  # initial values
  y0 <- rep(0, N)
  for (i in 1:N) y0[i] <- rnorm(1, unitfe[i]/(1 - phi), 1/sqrt(1 - phi^2))
  # simulate AR(1) panel sequence
  y <- matrix(0, nrow = T, ncol = N)
  y[1, ] <- unitfe + phi*y0 + shocks[1, ]
  for (t in 2:T) y[t, ] <- unitfe + phi*y[t - 1, ] + shocks[t, ]
  
  # construct panel data set
  y <- matrix(y, ncol = 1)
  dat <- data.frame(id = kronecker(c(1:N), rep(1, T)), 
                    time = kronecker(rep(1, N), c(1:T)), y = y)
  dat <- pdata.frame(dat, index = c("id", "time"))
  # add lagged variable
  dat$l1y <- lag(dat$y, 1)
  return(dat)
}

# objective function for indirect inference
# inputs:
#  coef: fixed effect coefficients
objective <- function(phi, shocks_sim, alpha_i, coef) {
  # compute simulated beta, returns |beta-beta_sim|
  H <- dim(shocks_sim)[3]
  # loop over sample paths
  coef_sim <- 0
  for (h in 1:H) {
    dats <- sim_dat(phi, shocks_sim[, , h], alpha_i[, h])
    fit_sim <- lm(y ~ l1y + factor(id), dats)
    coef_sim <- coef_sim + coef(fit_sim)[2]/H
  }
  # Use identity weighting matrix
  return( sum( (coef - coef_sim)^2 ) )
}


# a function to obtain indirect inference estimator for one simulation
# inputs:
# N: number of cross section units
# T: number of time series 

estimate_ii <- function(N, T, H, phi0) {
  # generate panel data
  alpha <- rnorm(N, 0, 1) # individual fixed effect
  e <- matrix(rnorm(N*T, 0, 1), nrow = T, ncol = N) # panel error
  # initial values
  y0 <- rep(0, N)
  for (i in 1:N) y0[i] <- rnorm(1, alpha[i]/(1 - phi0), 1/sqrt(1 - phi0^2))
  y <- matrix(0, nrow = T, ncol = N)
  y[1, ] <- alpha + phi0*y0 + e[1, ]
  for (t in 2:T) y[t, ] <- alpha + phi0*y[t - 1, ] + e[t, ]
  # construct panel data set
  y <- matrix(y, ncol = 1)
  dat <- data.frame(id = kronecker(c(1:N), rep(1, T)), 
                    time = kronecker(rep(1, N), c(1:T)), y = y)
  dat <- pdata.frame(dat, index = c("id", "time"))
  dat$l1y <- lag(dat$y, 1)
  # fixed effect estimate
  fit <- lm(y ~ l1y + factor(id), dat)
  coef_fe <- coef(fit)[2]
  
  # shocks and unit fe for sim data in indirect inference
  shocks_sim <- array(rnorm(N*T*H, 0, 1), c(T, N, H))
  alpha_sim <- matrix(rnorm(N*H, 0, 1), nrow = N, ncol = H)
  
  # estimation
  est <- optimize(objective, c(-1, 1), shocks_sim = shocks_sim, 
                  alpha_i = alpha_sim, coef = coef_fe)
  
  # return results
  return(est$minimum)
}

########### 3. Simulation #############
set.seed(88)

nCores <- 1   # number of CPUs for parallelization
#registerDoParallel(cores = nCores)

# loop over number of simulations
#results_par <- foreach(s = 1:S) %dorng% {
#  results <- estimate_ii(N, T, H, phi0)
#}

# convert results to a matrix
#sim_results <- t(matrix(unlist(results_par), 1, S)) 

# report the result
#print(c(mean(sim_results), sd(sim_results)), digits = 2)

results <- estimate_ii(N, T, H, phi0)





