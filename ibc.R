# 12/9/2020 Shuowen Chen
# Code for indirect inference bias correction for nonlinear panel

library(optimr)
library(alpaca) # for command feglm
library(plm)
library(foreign)
library(doParallel)
library(doRNG)

lfpdata <- read.dta("LFPmovers.dta")
lfpdata <- pdata.frame(lfpdata, index = c("id", "year"))

H <- 10 # number of simulation path
S <- 10 # number of simulations

# Regression specification
spec_fe <- lfp ~ kids0_2 + age + kids3_5 | id

# Check number of missing LFP values
nomissing <- !is.na(lfpdata$lfp)
pop <- length(lfpdata$lfp[nomissing])

# Estimate fitted value from the true dataset to create calibrated dataset
fit <- feglm(spec_fe, data = lfpdata, family = binomial(link = "probit"))
# true coefficient of the model
coef_fe <- coef(fit) 
# beta*X_{i,t}+\alpha_{i} and stack up by time series dimension
index <- matrix(predict(fit), nrow = length(unique(lfpdata$year)), byrow = FALSE) 

################# Create a calibrated panel dataset
cali_dat <- function(data, index) {
  N <- length(unique(data$id))  # number of individuals
  T <- length(unique(data$year)) # number of time periods
  lfp_s <- matrix(0, nrow = T, ncol = N)
  # Interpretation: beta*X_{i,t}+\alpha_{i} > rnorm(N)
  for (t in 1:T) lfp_s[t, ] <- (index[t, ] > rnorm(N))
  lfp_s <- matrix(lfp_s, ncol = 1, byrow = FALSE)
  # All other covariates are taken from the original data
  data_s <- data.frame(lfp = lfp_s, data[, c("kids0_2", "kids3_5", "kids6_17",
                                             "loghusbandincome", "age", "age2")], 
                       id = kronecker(sample(N), rep(1, T)), 
                       year = kronecker(rep(1, N), c(1:T)))
  # eliminate obs with all 0's or all 1's
  D <- matrix(data_s$id, ncol = 1)
  td <- kronecker(tapply(lfp_s, D, sum), matrix(1, T, 1)) 
  insample <- ((td > 0) & (td < T))
  # sort the data by id
  data_s <- data_s[order(data_s$id), ]
  data_s <- data_s[insample, ]
  return(data_s)
}

# A function to simulate data in indirect inference
sim_dat_multi <- function(phi, shocks, unitfe, covariate) {
  T <- dim(shocks)[1] # time series dimension
  N <- dim(shocks)[2] # cross section dimension
  
  lfp_s <- matrix(0, nrow = T, ncol = N)
  # Interpretation: beta*X_{i,t}+\alpha_{i} > rnorm(N)
  xbeta <- matrix(as.matrix(covariate) %*% matrix(phi), nrow = T, ncol = N)
  for (t in 1:T) lfp_s[t, ] <- (xbeta[t, ] + unitfe > shocks[t, ])
  lfp_s <- matrix(lfp_s, ncol = 1, byrow = FALSE)
  # All other covariates are taken from the original data
  data_s <- data.frame(lfp = lfp_s, covariate,
                       id = kronecker(c(1:N), rep(1, T)), 
                       year = kronecker(rep(1, N), c(1:T)))
  # eliminate obs with all 0's or all 1's
  D <- matrix(data_s$id, ncol = 1)
  td <- kronecker(tapply(lfp_s, D, sum),  matrix(1, T, 1)) 
  insample <- ((td > 0) & (td < T))
  data_s <- data_s[order(data_s$id), ] # sort data by id
  data_s <- data_s[insample, ]
  return(data_s)
}

# objective function for indirect inference
# inputs:
#  phi:      parameters to be optimized
#  coef:     fixed effect coefficients of the auxiliary model
#  reg_form: regression formula
#  X:        exogenous covariate
# output:
#  objective function for estimation
objective_multi <- function(phi, shocks_sim, alpha_i, coef, reg_form, X) {
  # compute simulated beta, returns |beta-beta_sim|
  print(phi)
  H <- dim(shocks_sim)[3]
  coef_sim <- rep(0, length(phi))
  # loop over sample paths
  for (h in 1:H) {
    dats <- sim_dat_multi(phi, shocks_sim[, , h], alpha_i[, h], X)
    fit_sim <- feglm(reg_form, dats, family = binomial(link = "probit"))
    coef_sim <- coef_sim + coef(fit_sim)/H
  }
  # Use identity weighting matrix
  return( sum( (coef - coef_sim)^2 ) )
}


# a function to obtain indirect inference estimator for one simulation
# inputs:
# data:      true dataset at hand
# index:     fitted value for constructing the data for simulation
# H:         number of simulation paths
# spec_fe:   regression specification for feglm
# xnames:    names of RHS exogenous covariate
# initials:  starting points for the optimization routine, should be a matrix
estimate_ii_multi <- function(data, index, H, spec_fe, xnames, initials) {
  # generate a calibrated panel data
  dat_syn <- cali_dat(data, index)
  fit_cali <- feglm(spec_fe, data = dat_syn, family = binomial(link = "probit"))
  # coefficients from the auxiliary model
  coef_cali <- coef(fit_cali)
  # number of pars of interest
  dim_par <- length(coef_cali)
  
  N <- length(unique(dat_syn$id))  # number of individuals
  T <- length(unique(dat_syn$year)) # number of time periods
  
  # shocks and unit fe for sim data in indirect inference
  shocks_sim <- array(rnorm(N*T*H, 0, 1), c(T, N, H))
  # Try different ways of simulating the individual fixed effect
  # 1. From normal dist.
  #al_sim <- matrix(rnorm(N*H, 0, 1), nrow = N, ncol = H) 
  # 2. Use the fixed effect estimate from fit_cali
  alf <- getFEs(fit_cali)$id 
  al_sim <- matrix(rep(alf, H), nrow = N, ncol = H)
  # 3. Use the sd of estimated individual fe in the synthetic data
  #sd_alpha <- sd(getFEs(fit_cali)$id)
  #al_sim <- matrix(rnorm(N*H, 0, sd_alpha), nrow = N, ncol = H)
  # 4. From logistic distribution
  #al_sim <- matrix(rlogis(N*H), nrow = N, ncol = H)
  # estimation procedure for indirect inference estimator
  if (dim_par == 1) {
    # one-dimensional optimization
    est <- optimize(objective_multi, c(-1, 1), shocks_sim = shocks_sim, 
                    alpha_i = al_sim, coef = coef_cali, reg_form = spec_fe, 
                    X = dat_syn[, ..xnames])
    results <- est$minimum
  } else {
    # results are sensitive to number of H
    est <- optim(initials, objective_multi, method = "Nelder-Mead",
                 shocks_sim = shocks_sim, alpha_i = al_sim,
                 coef = coef_cali, reg_form = spec_fe, 
                 X = dat_syn[, ..xnames])
    results <- est$par
  }
  return(list(ii = results, auxiliary = coef_cali))
}


# A function to generate DGP in which individual fixed effect is user-specified
# Inputs:
#  mu:    mean of individual fixed effects
#  sigma: standard error of individual fixed effects
dgp_mc <- function(coef0, covariate, N, T, mu, sigma) {
  # generate dependent variable
  lfp_s <- matrix(0, nrow = T, ncol = N) # placeholder
  alpha <- rnorm(N, mu, sigma) # individual fixed effects
  e <- matrix(rnorm(N*T, 0, 1), nrow = T, ncol = N) # panel error
  xbeta <- matrix(as.matrix(covariate) %*% matrix(coef0), nrow = T, ncol = N)
  for (t in 1:T) lfp_s[t, ] <- (xbeta[t, ] + alpha > e[t, ])
  # construct panel data
  lfp_s <- matrix(lfp_s, ncol = 1, byrow = FALSE)
  data_s <- data.frame(lfp = lfp_s, covariate, id = kronecker(sample(N), rep(1, T)), 
                       year = kronecker(rep(1, N), c(1:T)))
  # eliminate obs with all 0's or all 1's
  D <- matrix(data_s$id, ncol = 1)
  td <- kronecker(tapply(lfp_s, D, sum), matrix(1, T, 1)) 
  insample <- ((td > 0) & (td < T))
  # sort the data by id
  data_s <- data_s[order(data_s$id), ]
  data_s <- data_s[insample, ]
  return(data_s)
}

# a function to obtain indirect inference estimator for one simulation
# with user specified individual fixed effects
# inputs:
# data:      true dataset at hand
# index:     fitted value for constructing the data for simulation
# H:         number of simulation paths
# spec_fe:   regression specification for feglm
# xnames:    names of RHS exogenous covariate
# initials:  starting points for the optimization routine, should be a matrix
# mu:        mean of user specified individual fixed effects
# sigma:     standard deviation of user specified individual fixed effects
estimate_ii_mc <- function(data, H, spec_fe, xnames, initials, 
                           coef0, mu = 0, sigma = 1) {
  dat_syn <- dgp_mc(coef0, data[, ..xnames], length(unique(data$id)),
                    length(unique(data$year)), mu, sigma)
  fit <- feglm(spec_fe, data = dat_syn, family = binomial(link = "probit"))
  # coefficients from the auxiliary model
  coef_aux <- coef(fit)
  alf_est <- mean(getFEs(fit)$id)
  alf_sd <- sd(getFEs(fit)$id)
  aux <- c(coef_aux, alf_est, alf_sd)
  # number of pars of interest
  dim_par <- length(aux)
  
  N <- length(unique(dat_syn$id))  # number of individuals
  T <- length(unique(dat_syn$year)) # number of time periods
  
  # shocks and unit fe for sim data in indirect inference
  shocks_sim <- array(rnorm(N*T*H, 0, 1), c(T, N, H))
  # from normal dist.
  al_sim <- matrix(rnorm(N*H, mu, sigma), nrow = N, ncol = H) 
  
  # estimation procedure for indirect inference estimator
  if (dim_par == 1) {
    # one-dimensional optimization
    est <- optimize(objective_mc, c(-1, 1), shocks_sim = shocks_sim, 
                    alpha_i = al_sim, coef = aux, reg_form = spec_fe, 
                    X = dat_syn[, ..xnames])
    results <- est$minimum
  } else {
    # results are sensitive to number of H
    est <- optim(initials, objective_mc, method = "Nelder-Mead",
                 shocks_sim = shocks_sim, alpha_i = al_sim,
                 coef = aux, reg_form = spec_fe, 
                 X = dat_syn[, ..xnames])
    results <- est$par
  }
  return(list(ii = results, auxiliary = aux))
}

objective_mc <- function(phi, shocks_sim, alpha_i, coef, reg_form, X) {
  # compute simulated beta, returns |beta-beta_sim|
  print(phi)
  H <- dim(shocks_sim)[3]
  coef_sim <- rep(0, length(phi))
  # loop over sample paths
  for (h in 1:H) {
    dats <- sim_dat_mc(phi, shocks_sim[, , h], alpha_i[, h], X)
    fit_sim <- feglm(reg_form, dats, family = binomial(link = "probit"))
    alf_est <- mean(getFEs(fit)$id)
    alf_sd <- sd(getFEs(fit)$id)
    coef_tot <- c(coef(fit_sim), alf_est, alf_sd)
    coef_sim <- coef_sim + coef_tot/H
  }
  # Use identity weighting matrix
  return( sum( (coef - coef_sim)^2 ) )
}

# A function to simulate data in indirect inference
sim_dat_mc <- function(phi, shocks, unitfe, covariate) {
  T <- dim(shocks)[1] # time series dimension
  N <- dim(shocks)[2] # cross section dimension
  
  lfp_s <- matrix(0, nrow = T, ncol = N)
  # Interpretation: beta*X_{i,t}+\alpha_{i} > rnorm(N)
  xbeta <- matrix(as.matrix(covariate) %*% matrix(head(phi, -2)), 
                  nrow = T, ncol = N)
  for (t in 1:T) lfp_s[t, ] <- (xbeta[t, ] + unitfe > shocks[t, ])
  lfp_s <- matrix(lfp_s, ncol = 1, byrow = FALSE)
  # All other covariates are taken from the original data
  data_s <- data.frame(lfp = lfp_s, covariate,
                       id = kronecker(c(1:N), rep(1, T)), 
                       year = kronecker(rep(1, N), c(1:T)))
  # eliminate obs with all 0's or all 1's
  D <- matrix(data_s$id, ncol = 1)
  td <- kronecker(tapply(lfp_s, D, sum),  matrix(1, T, 1)) 
  insample <- ((td > 0) & (td < T))
  data_s <- data_s[order(data_s$id), ] # sort data by id
  data_s <- data_s[insample, ]
  return(data_s)
}

########### 3. Simulation #############
set.seed(88, kind = "L'Ecuyer-CMRG") # for reproduction
nCores <- 1 # number of CPUs for parallelization
registerDoParallel(cores = nCores)
mu <- 0
sigma <- 1
pm <- c(coef_fe, mu, sigma) # initial values

#results_mc <- estimate_ii_mc(lfpdata, H, spec_fe, c("kids0_2", "age", "kids3_5"), 
#                          pm, coef_fe, mu, sigma)

#results <- estimate_ii_multi(lfpdata, index, H, spec_fe, c("kids0_2", "age", "kids3_5"), coef_fe)
#results_par <- foreach(s = 1:S) %dorng% {
#  results <- estimate_ii_multi(lfpdata, index, H, spec_fe,
#                               c("kids0_2", "age", "kids3_5"), pm)
#}

# convert results to a matrix
#ii_matrix <- sapply(results_par, function(x) return(x$ii))
#rowMeans(ii_matrix)
#apply(ii_matrix, 1, sd)
results_par <- foreach(s = 1:S) %dorng% {
  results_mc <- estimate_ii_mc(lfpdata, H, spec_fe, c("kids0_2", "age", "kids3_5"), 
                               pm, coef_fe, mu, sigma)
}


