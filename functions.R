# 1/17/2020 Shuowen Chen
# This script stores all the function for the indirect inference project

################# Functions
# Create a calibrated panel dataset
# Inputs: 
# data:  true dataset (LFP)
# fit:   regression fit of the true dataset
# model: indicate whether the model is static or dynamic
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
  if (model == "static") {
    insample <- ((td > 0) & (td < T))
  } else {
    td2 <- kronecker(tapply(laglfp_s, D, sum), matrix(1, T, 1))
    insample <- ((td > 0) & (td < T) & (td2 > 0) & (td2 < T))
  }
  # sort the data by id
  data_s <- data_s[order(data_s$id), ]
  data_s <- data_s[insample, ]
  return(data_s)
}

# Simulate data in indirect inference
# Inputs:
# phi:       estimate of parameter of interest
# shocks:    panel shocks
# unitfe:    individual fixed effects
# covariate: regressors
sim_dat <- function(phi, shocks, unitfe, covariate,
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
                         id = kronecker(sample(N), rep(1, T)), 
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
  if (model == "static") {
    insample <- ((td > 0) & (td < T))
  } else {
    td2 <- kronecker(tapply(laglfp_s, D, sum), matrix(1, T, 1))
    insample <- ((td > 0) & (td < T) & (td2 > 0) & (td2 < T))
  }
  data_s <- data_s[order(data_s$id), ] # sort data by id
  data_s <- data_s[insample, ]
  return(data_s)
}

# objective function for indirect inference
# Inputs:
#  phi:      parameters to be optimized
#  coef:     fixed effect coefficients of the auxiliary model
#  reg_form: regression formula
#  X:        exogenous covariate
objective <- function(phi, shocks_sim, alpha_i, coef, reg_form, X,
                      model = c("static", "dynamic")) {
  model <- match.arg(model)
  # compute simulated beta, returns |beta-beta_sim|
  print(phi)
  H <- dim(shocks_sim)[3]
  coef_sim <- rep(0, length(phi))
  # loop over sample paths
  for (h in 1:H) {
    dats <- sim_dat(phi, shocks_sim[, , h], alpha_i[, h], X, model)
    fit_sim <- try(feglm(reg_form, dats, family = binomial(link = "probit")))
    coef_sim <- tryCatch(
      {coef_sim + coef(fit_sim)/H},
      error = function(e) coef_sim
    )
  }
  # Use identity weighting matrix
  return(sum((coef - coef_sim)^2))
}

iise <- function(coef_ii, coef_fe, fit_fe) {
  cov <- vcov(fit_fe)
  jac <- -2*(coef_fe - coef_ii)
  return(sum((coef - phi)^2))
}

# a function to obtain indirect inference estimator for one simulation
# INPUTS:
# data:      true dataset at hand
# fit0:      fitted model for constructing the data for simulation
# H:         number of simulation paths
# spec_fe:   regression specification for feglm
# xnames:    names of RHS exogenous covariate
# initials:  starting points for the optimization routine
# model:     whether the model has a lagged dependent variable
# var:       covariate of interest to compute APE
# var_type:  whether variable is continuous or binary
# at:        if var_type = "continuous", the value at which to compute the APE
# OUTPUTS:
# ii:        indirect inference estimator
# auxiliary: FE estimator without bias correction
estimate_ii <- function(data, fit0, H, spec_fe, xnames, initials,
                        model = c("static", "dynamic"), var, var_type, at) {
  model <- match.arg(model)
  # generate a calibrated panel data
  dat_syn <- cali_dat(data, fit0, model)
  N <- length(unique(dat_syn$id))  # number of individuals
  T <- length(unique(dat_syn$year)) # number of time periods
  
  # coefficients from the auxiliary model
  fit_aux <- feglm(spec_fe, data = dat_syn, family = binomial(link = "probit"))
  coef_aux <- coef(fit_aux)
  #coef_aux <- ape(fit_aux, coef(fit_aux), T, var, var_type, at)
  # Analytical standard errors for the uncorrected FE
  # vcov_fe <- vcov(fit_aux)
  ase <- summary(fit_aux)[1]$cm[, 2]
  # number of pars of interest
  dim_par <- length(coef_aux)
  
  # Shocks and unit fe for sim data in indirect inference
  shocks_sim <- array(rnorm(N*T*H, 0, 1), c(T, N, H))
  # use the fixed effect estimate from fit_cali
  alf <- getFEs(fit_aux)$id 
  al_sim <- matrix(rep(alf, H), nrow = N, ncol = H)
  
  # estimation procedure for indirect inference estimator
  if (dim_par == 1) {
    # one-dimensional optimization
    est <- optimize(objective, c(-1, 1), shocks_sim = shocks_sim, 
                    alpha_i = al_sim, coef = coef_aux, reg_form = spec_fe, 
                    X = dat_syn[, ..xnames], var = var, var_type = var_type,
                    at = at)
    results <- est$minimum
  } else {
    est <- optim(initials, objective, method = "Nelder-Mead",
                 shocks_sim = shocks_sim, alpha_i = al_sim,
                 coef = coef_aux, reg_form = spec_fe, model = model,
                 X = dat_syn[, ..xnames]) 
    # results <- ape(fit_aux, est$par, T, var, var_type, at)
    results <- est$par
  }
  # An iterated estimator: use the indirect inference estimator to 
  # construct the data and reestimate things
  #if (iterated == TRUE) {
  #  est <- optim(initials, objective, method = "Nelder-Mead",
  #               shocks_sim = shocks_sim, alpha_i = al_sim,
  #               coef = coef_aux, reg_form = spec_fe, model = model,
  #               X = dat_syn[, ..xnames]) 
  #}
  return(list(ii = results, auxiliary = coef_aux, ase = ase))
}

# a function to compute average partial effect for panel probit
# NOTE: need to modify the function to compute all APEs
# Inputs:
#  fit:      regression fit of FE estimation
#  T:        time series dimension
#  var:      covariate of interest
#  var_type: whether the support of var is discrete or continuous
#  at:       if var_type is continuous, specify where to evaluate 
#            the partial effect, mean or a quantile
ape <- function(fit, coef_hat, T, var, var_type = c("discrete", "continuous"),
                at = c("mean", 0.5)) {
  alf_hat <- getFEs(fit)$id
  coef_names <- names(coef(fit))
  d0 <- d1 <- fit$data[, ..coef_names]
  if (var_type == "binary") {
    # counterfactual
    d0[, var] <- 0
    d1[, var] <- 1
    ape_est <- mean(pnorm(as.matrix(d1) %*% matrix(coef_hat) + rep(alf_hat, each = T)) - 
                      pnorm(as.matrix(d0) %*% matrix(coef_hat) + rep(alf_hat, each = T)))
  } else if (var_type == "continuous") {
    if (at == "mean") {
      d0[, var] <- apply(d0[, ..var], 2, mean)
    } else {
      d0[, var] <- apply(d0[, ..var], 2, quantile, probs = at)
    }
    ape_est <- coef_hat[var]*mean(dnorm(as.matrix(d0) %*% matrix(coef_hat) + rep(alf_hat, each = T)))
  } 
  names(ape_est) <- "APE"
  return(ape_est)
}

objective_ape <- function(phi, shocks_sim, alpha_i, coef, 
                          reg_form, X, model = c("static", "dynamic"),
                          var, var_type, at) {
  model <- match.arg(model)
  H <- dim(shocks_sim)[3]
  T <- dim(shocks_sim)[1]
  coef_sim <- 0
  # loop over sample paths
  for (h in 1:H) {
    dats <- sim_dat(phi, shocks_sim[, , h], alpha_i[, h], X, model)
    fit_sim <- try(feglm(reg_form, dats, family = binomial(link = "probit")))
    coef_sim <- tryCatch(
      {coef_sim + ape(fit_sim, T, var, var_type, at)/H},
      error = function(e) coef_sim
    )
  }
  print(coef_sim)
  # Use identity weighting matrix
  return(sum((coef - coef_sim)^2))
}


# a function that produces the diagnostic statistics
table_simulation <- function(est, est0, ase) {
  tab <- matrix(0, nrow = 4, ncol = length(est0))
  rownames(tab) <- c('Bias', 'Std Dev', 'RMSE', 'p.95 (ASE)')
  colnames(tab) <- names(est0)
  tab[1, ] <- 100*(apply(est, 1, mean)/est0 - 1)
  tab[2, ] <- 100*(apply(est/est0, 1, sd))
  tab[3, ] <- 100*sqrt((apply((est/est0 - 1)^2, 1, mean)))
  temp <- matrix(rep(est0, 500), ncol = 500)
  tab[4, ] <- apply((est + qnorm(.05/2)*ase <= temp) & 
                      (est + qnorm(1 - .05/2)*ase >= temp), 1, mean)
  return(tab)
}


