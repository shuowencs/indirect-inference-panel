# 11/15/2020 Shuowen Chen
# Calibrated Static profit for the project: Indirect Inference for Nonlinear Panel

library(alpaca)
library(doParallel)
library(doRNG)
library(plm)
library(foreign)
library(speedglm)

#filepath <- "/Users/shuowenchen/desktop/Research/indirect-inference-nonlinear-panel"
#setwd(filepath)

################# Load data
lfpdata <- read.dta("LFPmovers.dta")
lfpdata <- pdata.frame(lfpdata, index = c("id", "year"))
# VALUES OF THE PARAMETERS;
N <- length(unique(lfpdata$id))  # number of individuals
T <- length(unique(lfpdata$year)) # number of time periods
H <- 5
S <- 1

# Regression specification
spec_fe <- lfp ~ age | id
spec <- lfp ~ age + factor(id)

# Check number of missing LFP values
nomissing <- !is.na(lfpdata$lfp)
pop <- length(lfpdata$lfp[nomissing])

# Estimate fitted value from the true dataset to create calibrated dataset
fit <- feglm(spec_fe, data = lfpdata, family = binomial(link = "probit"))
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
  data_s <- data.frame(lfp = lfp_s, kids0_2 = data$kids0_2, kids3_5 = data$kids3_5, 
                       kids6_17 = data$kids6_17, loghusbandincome = data$loghusbandincome, 
                       age = data$age, age2 = data$age2, id = kronecker(sample(N), rep(1, T)), 
                       year = data$year)
  # eliminate obs with all 0's or all 1's
  D <- matrix(data_s$id, ncol = 1)
  td <- kronecker(tapply(lfp_s, D, sum),  matrix(1, T, 1)) 
  insample <- ((td > 0) & (td < T))
  data_s <- data_s[insample, ]
  return(data_s)
}

# A function to simulate data in indirect inference
sim_dat <- function(phi, shocks, unitfe, covariate) {
  T <- dim(shocks)[1] # time series dimension
  N <- dim(shocks)[2] # cross section dimension
  
  lfp_s <- matrix(0, nrow = T, ncol = N)
  # Interpretation: beta*X_{i,t}+\alpha_{i} > rnorm(N)
  temp <- matrix(covariate, nrow = T, ncol = N)
  for (t in 1:T) lfp_s[t, ] <- (phi*temp[t, ] + unitfe > shocks[t, ])
  lfp_s <- matrix(lfp_s, ncol = 1, byrow = FALSE)
  # All other covariates are taken from the original data
  data_s <- data.frame(lfp = lfp_s, age = covariate, 
                       id = kronecker(sample(N), rep(1, T)), 
                       year = kronecker(rep(1, N), c(1:T)))
  # eliminate obs with all 0's or all 1's
  D <- matrix(data_s$id, ncol = 1)
  td <- kronecker(tapply(lfp_s, D, sum),  matrix(1, T, 1)) 
  insample <- ((td > 0) & (td < T))
  data_s <- data_s[insample, ]
  return(data_s)
}

# objective function for indirect inference
# inputs:
#  coef: fixed effect coefficients
#  reg_form: regression formulat
#  X: exogenous covariate
objective <- function(phi, shocks_sim, alpha_i, coef, reg_form, X) {
  # compute simulated beta, returns |beta-beta_sim|
  H <- dim(shocks_sim)[3]
  # loop over sample paths
  coef_sim <- 0
  for (h in 1:H) {
    dats <- sim_dat(phi, shocks_sim[, , h], alpha_i[, h], X)
    fit_sim <- feglm(reg_form, dats, family = binomial(link = "probit"))
    coef_sim <- coef_sim + coef(fit_sim)/H
  }
  # Use identity weighting matrix
  return( sum( (coef - coef_sim)^2 ) )
}


# a function to obtain indirect inference estimator for one simulation
# inputs:
# data: true dataset at hand
# N: number of cross section units
# T: number of time series 
# H: number of simulation paths
# spec_fe: regression specification
# covariate: RHS exogenous covariate

estimate_ii <- function(data, index, N, T, H, spec_fe, covariate) {
  # generate a calibrated panel data
  dat_syn <- cali_dat(data, T, N, index)
  # fixed effect estimate
  fit_cali <- feglm(spec_fe, data = dat_syn, family = binomial(link = "probit"))
  coef_cali <- coef(fit_cali)
  
  # shocks and unit fe for sim data in indirect inference
  shocks_sim <- array(rnorm(N*T*H, 0, 1), c(T, N, H))
  # Try different ways of simulating the individual fixed effect
  # 1. From normal dist.
  alpha_sim <- matrix(rnorm(N*H, 0, 1/4), nrow = N, ncol = H) 
  #alpha_sim2 <- matrix(rexp(N*H), nrow = N, ncol = H) 
  # 2. Use the true DGP's fixed effect
  #alf <- matrix(predict(fit) - coef_fe*lfpdata$age, nrow = T, byrow = FALSE)
  #alpha <- matrix(rep(alf[1, ], H), ncol = H)
  # 3. Use the sd of estimated individual fe in the synthetic data
  # Note: this doesn't work well because individual fe is biased estimate
  #f <- speedglm(lfp ~ age +factor(id), dat_syn, family = binomial(link = "probit"))
  #sd_alpha <- sd(coef(f)[-c(1:2)])
  #al_sim <- matrix(rnorm(N*H, 0, sd_alpha), nrow = N, ncol = H)
  
  # estimation procedure for indirect inference estimator
  # one-dimensional optimization
  est <- optimize(objective, c(-1, 1), shocks_sim = shocks_sim, 
                  alpha_i = alpha_sim, coef = coef_cali, 
                  reg_form = spec_fe, X = covariate)
  
  # return results
  return(est$minimum)
}

########### 3. Simulation #############
set.seed(88)

nCores <- 1   # number of CPUs for parallelization
registerDoParallel(cores = nCores)

# loop over number of simulations
results <- estimate_ii(lfpdata, index, N, T, H, spec_fe, lfpdata$age)
#results_par <- foreach(s = 1:S) %dorng% {
#  results <- estimate_ii(lfpdata, index, N, T, H, spec_fe, lfpdata$kids0_2)
#}

# convert results to a matrix
#sim_results <- t(matrix(unlist(results_par), 1, S)) 

# report the result
#print(c(mean(sim_results), sd(sim_results)), digits = 2)






