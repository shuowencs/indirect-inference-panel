# 12/11/2020 Shuowen Chen
# Code for indirect inference bias correction for nonlinear panel
# The script allows for static and dynamic nonlinear panel

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
spec_s <- lfp ~ kids0_2 + kids3_5 + kids6_17 + loghusbandincome + age + age2 | id
spec_d <- lfp ~ laglfp + kids0_2 + kids3_5 + kids6_17 + loghusbandincome + age + age2 | id

# Estimate fitted value from the true dataset to create calibrated dataset
fit_s <- feglm(spec_s, data = lfpdata, family = binomial(link = "probit"))
fit_d <- feglm(spec_d, data = lfpdata, family = binomial(link = "probit"))
coef0_s <- coef(fit_s)
coef0_d <- coef(fit_d)

################# Create a calibrated panel dataset
cali_dat <- function(data, fit, model = c("static", "dynamic")) {
  model <- match.arg(model)
  N <- length(unique(data$id))  # number of individuals
  T <- length(unique(data$year)) # number of time periods
  # true coefficient of the model
  coef0 <- coef(fit) 
  lfp_s <- matrix(0, nrow = T, ncol = N)
  if (model == "static") {
    # beta*X_{i,t}+\alpha_{i} and stack up by time series dimension
    index <- matrix(predict(fit), nrow = T, byrow = FALSE) 
    # Interpretation: beta*X_{i,t}+\alpha_{i} > rnorm(N)
    for (t in 1:T) lfp_s[t, ] <- (index[t, ] > rnorm(N))
    lfp_s <- matrix(lfp_s, ncol = 1, byrow = FALSE)
    # All other covariates are taken from the original data
    data_s <- data.frame(lfp = lfp_s, data[, c("kids0_2", "kids3_5", "kids6_17",
                                               "loghusbandincome", "age", "age2")], 
                         id = kronecker(sample(N), rep(1, T)), 
                         year = kronecker(rep(1, N), c(1:T)))
  } else {
    index <- matrix(predict(fit) - coef0[1]*data$laglfp, nrow = T, byrow = FALSE)
    lfp0 <- data$laglfp[data$year == 1] #initial value
    for (t in 1:T) {
      # Interpretation: lfp0*coef0[1] + X'*coef0[2:end] > rnorm(N0),
      # Determines the lfp in the next period
      lfp_s[t, ] <- (lfp0*coef0[1] + index[t, ] > rnorm(N))
      # update so as to simulate the next period lfp
      lfp0 <- lfp_s[t, ] 
    }
    # simulate lagged dependent variable
    laglfp_s <- matrix(rbind(data$laglfp[data$year == 1], lfp_s[-T, ]), 
                       ncol = 1, byrow = FALSE)
    lfp_s <- matrix(lfp_s, ncol = 1, byrow = FALSE)
    # All other covariates are taken from the original data
    data_s <- data.frame(lfp = lfp_s, laglfp = laglfp_s, 
                         data[, c("kids0_2", "kids3_5", "kids6_17",
                                  "loghusbandincome", "age", "age2")],
                         id = kronecker(sample(N), rep(1, T)), 
                         year = kronecker(rep(1, N), c(1:T)))
  }
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
sim_dat_multi <- function(phi, shocks, unitfe, covariate,
                          model = c("static", "dynamic")) {
  model <- match.arg(model)
  T <- dim(shocks)[1] # time series dimension
  N <- dim(shocks)[2] # cross section dimension
  lfp_s <- matrix(0, nrow = T, ncol = N)
  if (model == "static") {
    # Interpretation: beta*X_{i,t}+\alpha_{i} > rnorm(N)
    xbeta <- matrix(as.matrix(covariate) %*% matrix(phi), nrow = T, ncol = N)
    for (t in 1:T) lfp_s[t, ] <- (xbeta[t, ] + unitfe > shocks[t, ])
    lfp_s <- matrix(lfp_s, ncol = 1, byrow = FALSE)
    # All other covariates are taken from the original data
    data_s <- data.frame(lfp = lfp_s, covariate,
                         id = kronecker(c(1:N), rep(1, T)), 
                         year = kronecker(rep(1, N), c(1:T)))
  } else {
    # x'beta without lagged dependent variable
    xbeta <- matrix(as.matrix(covariate[, -1]) %*% matrix(tail(phi, -1)),
                    nrow = T, ncol = N)
    # initial value of lfp: 1st laglfp obs for each unit
    lfp0 <- covariate[, 1][seq(1, N*T, T)]
    for (t in 1:T) {
      # Interpretation: lfp0*phi[1] + X'*phi[2:end] > rnorm(N0),
      # Determines the lfp in the next period
      lfp_s[t, ] <- (lfp0*phi[1] + xbeta[t, ] + unitfe > shocks[t, ])
      # update so as to simulate the next period lfp
      lfp0 <- lfp_s[t, ] 
    }
    # simulate lagged dependent variable
    laglfp_s <- matrix(rbind(unlist(covariate[, 1][seq(1, N*T, T)]), 
                             lfp_s[-T, ]), ncol = 1, byrow = FALSE)
    lfp_s <- matrix(lfp_s, ncol = 1, byrow = FALSE)
    # All other covariates are taken from the original data
    data_s <- data.frame(lfp = lfp_s, laglfp = laglfp_s, covariate[, -1],
                         id = kronecker(sample(N), rep(1, T)), 
                         year = kronecker(rep(1, N), c(1:T)))
  }
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
objective_multi <- function(phi, shocks_sim, alpha_i, coef, reg_form, X,
                            model = c("static", "dynamic")) {
  model <- match.arg(model)
  # compute simulated beta, returns |beta-beta_sim|
  print(phi)
  H <- dim(shocks_sim)[3]
  coef_sim <- rep(0, length(phi))
  # loop over sample paths
  for (h in 1:H) {
    dats <- sim_dat_multi(phi, shocks_sim[, , h], alpha_i[, h], X, model)
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
estimate_ii_multi <- function(data, fit0, H, spec_fe, xnames, initials,
                              model = c("static", "dynamic")) {
  model <- match.arg(model)
  # generate a calibrated panel data
  dat_syn <- cali_dat(data, fit0, model)
  fit_aux <- feglm(spec_fe, data = dat_syn, family = binomial(link = "probit"))
  # coefficients from the auxiliary model
  coef_aux <- coef(fit_aux)
  # number of pars of interest
  dim_par <- length(coef_aux)
  
  N <- length(unique(dat_syn$id))  # number of individuals
  T <- length(unique(dat_syn$year)) # number of time periods
  
  # shocks and unit fe for sim data in indirect inference
  shocks_sim <- array(rnorm(N*T*H, 0, 1), c(T, N, H))
  # Try different ways of simulating the individual fixed effect
  # 1. From normal dist.
  #al_sim <- matrix(rnorm(N*H, 0, 1), nrow = N, ncol = H) 
  # 2. Use the fixed effect estimate from fit_cali
  alf <- getFEs(fit_aux)$id 
  al_sim <- matrix(rep(alf, H), nrow = N, ncol = H)
  # 3. Use the sd of estimated individual fe in the synthetic data
  #sd_alpha <- sd(getFEs(fit_cali)$id)
  #al_sim <- matrix(rnorm(N*H, 0, sd_alpha), nrow = N, ncol = H)
  # 4. From logistic distribution
  #al_sim <- matrix(rlogis(N*H), nrow = N, ncol = H)
  # estimation procedure for indirect inference estimator
  est <- optim(initials, objective_multi, method = "Nelder-Mead",
               shocks_sim = shocks_sim, alpha_i = al_sim,
               coef = coef_aux, reg_form = spec_fe, model = model,
               X = dat_syn[, ..xnames])
  results <- est$par
  return(list(ii = results, auxiliary = coef_aux))
}

########### 3. Simulation #############
set.seed(88, kind = "L'Ecuyer-CMRG") # for reproduction
nCores <- 1 # number of CPUs for parallelization
registerDoParallel(cores = nCores)

#results_static <- estimate_ii_multi(lfpdata, fit_s, H, spec_s, 
#                                    c("kids0_2", "kids3_5", "kids6_17", 
#                                      "loghusbandincome", "age", "age2"), coef0_s,
#                                    model = "static")

results_dynamic <- estimate_ii_multi(lfpdata, fit_d, H, spec_d, 
                                     c("laglfp", "kids0_2", "kids3_5", "kids6_17", 
                                       "loghusbandincome", "age", "age2"), coef0_d,
                                     model = "dynamic")
#results_par <- foreach(s = 1:S) %dorng% {
#  results <- estimate_ii_multi(lfpdata, index, H, spec_fe,
#                               c("kids0_2", "age", "kids3_5"), pm)
#}

# convert results to a matrix
#ii_matrix <- sapply(results_par, function(x) return(x$ii))
#rowMeans(ii_matrix)
#apply(ii_matrix, 1, sd)



