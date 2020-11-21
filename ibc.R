# 11/20/2020 Shuowen Chen
# Code for indirect inference bias correction for nonlinear panel

library(optimr)
library(alpaca)
# goal: rewrite the code tomorrow to allow for arbitrary number of coefficients for estimation

lfpdata <- read.dta("LFPmovers.dta")
lfpdata <- pdata.frame(lfpdata, index = c("id", "year"))

N <- length(unique(lfpdata$id))  # number of individuals
T <- length(unique(lfpdata$year)) # number of time periods
H <- 10
S <- 1

# Regression specification
spec_fe <- lfp ~ kids0_2 + age | id
spec2 <- lfp ~ kids0_2 + age + factor(id)

# Check number of missing LFP values
nomissing <- !is.na(lfpdata$lfp)
pop <- length(lfpdata$lfp[nomissing])

# Estimate fitted value from the true dataset to create calibrated dataset
fit <- feglm(spec_fe, data = lfpdata, family = binomial(link = "probit"))
# true coefficient of the model
coef_fe <- coef(fit) 
# beta*X_{i,t}+\alpha_{i} and stack up by time series dimension
index <- matrix(predict(fit), nrow = T, byrow = FALSE) 

################# Create a calibrated panel dataset
cali_dat <- function(data, T, N, index) {
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
  td <- kronecker(tapply(lfp_s, D, sum),  matrix(1, T, 1)) 
  insample <- ((td > 0) & (td < T))
  data_s <- data_s[insample, ]
  return(data_s)
}

# A function to simulate data in indirect inference
sim_dat_multi <- function(phi, shocks, unitfe, covariate) {
  T <- dim(shocks)[1] # time series dimension
  N <- dim(shocks)[2] # cross section dimension
  
  lfp_s <- matrix(0, nrow = T, ncol = N)
  # Interpretation: beta*X_{i,t}+\alpha_{i} > rnorm(N)
  #mat <- as.matrix(covariate)
  #temp <- matrix(phi*mat[, 1] + coef_true[2]*mat[, 2], nrow = T, ncol = N)
  temp <- matrix(as.matrix(covariate) %*% matrix(phi), nrow = T, ncol = N)
  for (t in 1:T) lfp_s[t, ] <- (temp[t, ] + unitfe > shocks[t, ])
  lfp_s <- matrix(lfp_s, ncol = 1, byrow = FALSE)
  # All other covariates are taken from the original data
  data_s <- data.frame(lfp = lfp_s, covariate,
                       id = kronecker(c(1:N), rep(1, T)), 
                       year = kronecker(rep(1, N), c(1:T)))
  #data_s <- data.frame(lfp = lfp_s, covariate,
  #                     id = kronecker(sample(N), rep(1, T)), 
  #                     year = kronecker(rep(1, N), c(1:T)))
  # eliminate obs with all 0's or all 1's
  D <- matrix(data_s$id, ncol = 1)
  td <- kronecker(tapply(lfp_s, D, sum),  matrix(1, T, 1)) 
  insample <- ((td > 0) & (td < T))
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
# N:         number of cross section units
# T:         number of time series 
# H:         number of simulation paths
# spec_fe:   regression specification for feglm
# spec_2:    specification for speedglm (to get sd of estimated individual fixed effects)
# covariate: RHS exogenous covariate
# initials:  starting points for the optimization routine, should be a matrix
estimate_ii_multi <- function(data, index, N, T, H, spec_fe, spec_2, covariate,
                              initials) {
  # generate a calibrated panel data
  dat_syn <- cali_dat(data, T, N, index)
  
  fit_cali <- feglm(spec_fe, data = dat_syn, family = binomial(link = "probit"))
  # coefficients from the auxiliary model
  coef_cali <- coef(fit_cali)
  # number of pars of interest
  dim_par <- length(coef_cali)
  
  # shocks and unit fe for sim data in indirect inference
  shocks_sim <- array(rnorm(N*T*H, 0, 1), c(T, N, H))
  # Try different ways of simulating the individual fixed effect
  # 1. From normal dist.
  #alpha_sim <- matrix(rnorm(N*H, 0, 1), nrow = N, ncol = H) 
  # 2. Use the true DGP's fixed effect
  #alf <- matrix(predict(fit) - coef_fe*lfpdata$kids0_2, nrow = T, byrow = FALSE)
  #alpha <- matrix(rep(alf[1, ], H), ncol = H)
  # 3. Use the sd of estimated individual fe in the synthetic data
  # Note: this doesn't work well because individual fe is biased estimate
  f <- speedglm(spec_2, dat_syn, family = binomial(link = "probit"))
  sd_alpha <- sd(coef(f)[-c(1:(dim_par + 1))])
  al_sim <- matrix(rnorm(N*H, 0, sd_alpha), nrow = N, ncol = H)
  
  # estimation procedure for indirect inference estimator
  if (dim_par == 1) {
    # one-dimensional optimization
    est <- optimize(objective_multi, c(-1, 1), shocks_sim = shocks_sim, 
                    alpha_i = al_sim, coef = coef_cali, reg_form = spec_fe, 
                    X = covariate)
    results <- est$minimum
  } else {
    # use multiple starting points
    # results are sensitive to the starting points and number of H
    est <- multistart(initials, objective_multi, method = "L-BFGS-B", 
                      lower = rep(-1, dim_par), upper = rep(1, dim_par), 
                      shocks_sim = shocks_sim, alpha_i = al_sim, 
                      coef = coef_cali, reg_form = spec_fe, X = covariate)
    # return the estimates that yield the small function value
    ind <- which.min(est$value)
    colnames(est) <- colnames(covariate)
    results <- est[ind, c(1:dim_par)]
  }
  return(results)
}


########### 3. Simulation #############
set.seed(88) # for reproduction
#pm <- as.matrix(rbind(c(0, 0, 0), c(-0.7, 0.2, -0.3), c(0.5, 1, -0.25), 
#                      c(-0.7, 0.3, 0.4), c(-0.5, -1, -0.3)))

pm <- as.matrix(rbind(c(0, 0), c(-0.5, -0.2), c(0.5, -0.25), c(-0.6, 0.2)))

results <- estimate_ii_multi(lfpdata, index, N, T, H, spec_fe, spec2, 
                             lfpdata[, c("kids0_2", "age")], pm)

