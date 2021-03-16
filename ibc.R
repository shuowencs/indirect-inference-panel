# 12/18/2020 Shuowen Chen
# Code for indirect inference bias correction for nonlinear panel
# The script allows for static and dynamic nonlinear panel

library(optimr)
library(alpaca) # for command feglm
library(plm)
library(foreign)
library(doParallel)
library(doRNG)

source("functions.R")
lfpdata <- read.dta("LFPmovers.dta")
lfpdata <- pdata.frame(lfpdata, index = c("id", "year"))

H <- 10 # number of simulation path
S <- 10 # number of simulations

################# DGP Construction
# Regression specification
spec_s <- lfp ~ kids0_2 + kids3_5 + kids6_17 + loghusbandincome + age + age2 | id
#spec_d <- lfp ~ laglfp + kids0_2 + kids3_5 + age | id
spec_d <- lfp ~ laglfp + kids0_2 + kids3_5 + kids6_17 + loghusbandincome + age + age2 | id

# Estimate fitted value from the true dataset to create calibrated dataset
fit_s <- feglm(spec_s, data = lfpdata, family = binomial(link = "probit"))
fit_d <- feglm(spec_d, data = lfpdata, family = binomial(link = "probit"))

# True coefficients for calibration exercise
coef0_s <- coef(fit_s) 
coef0_d <- coef(fit_d)
#ape0_age_s <- ape(fit_s, coef0_s, 9, "age", "con", at = "mean")
#ape0_age_d <- ape(fit_d, coef0_d, 9, "age", "con", at = "mean")


########### 3. Simulation #############
set.seed(88, kind = "L'Ecuyer-CMRG") # for reproduction
nCores <- 1 # number of CPUs for parallelization
registerDoParallel(cores = nCores)


#results_par <- foreach(s = 1:S) %dorng% {
results_dynamic2 <- estimate_ii(lfpdata, fit_d, H, spec_d, 
                                c("laglfp", "kids0_2", "kids3_5", "kids6_17", "loghusbandincome", "age", "age2"), 
                                coef0_d, model = "dynamic", var = "age", var_type = "continuous", at = "mean")
#}
set.seed(88, kind = "L'Ecuyer-CMRG")
results_static <- estimate_ii(lfpdata, fit_s, H, spec_s, 
                              c("kids0_2", "kids3_5", "kids6_17", "loghusbandincome", "age", "age2"), 
                                coef0_s, model = "static", var = "age", var_type = "continuous", at = "mean")
########### 4. Table Results #############
# convert results to a matrix
ii_matrix_d <- sapply(results_pard, function(x) return(x$ii))
fe_matrix_d <- sapply(results_pard, function(x) return(x$auxiliary))
ii_matrix_s <- sapply(results_pars, function(x) return(x$ii))
fe_matrix_s <- sapply(results_pars, function(x) return(x$auxiliary))
se_matrix_d <- sapply(results_pard, function(x) return(x$ase))
se_matrix_s <- sapply(results_pars, function(x) return(x$ase))

table_ii_d <- table_simulation(ii_matrix_d, coef0_d, se_matrix_d)
table_fe_d <- table_simulation(fe_matrix_d, coef0_d, se_matrix_d)

table_ii_s <- table_simulation(ii_matrix_s, coef0_s, se_matrix_s)
table_fe_s <- table_simulation(fe_matrix_s, coef0_s, se_matrix_s)







