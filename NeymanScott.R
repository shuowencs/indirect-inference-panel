# 3/15/2021 Shuowen Chen
# Neyman Scott example

# generate x_it ~ N(mu_i, sigma^2)
set.seed(88)
n <- 250
T <- 5
H <- 1
#mu_i <- 0.1 * seq(n)
mu_i <- rnorm(n, mean = 0, sd = 1)
sigma <- 2

dat <- matrix(0, nrow = n, ncol = T)
for (i in 1:n) dat[i, ] <- rnorm(T, mean = mu_i[i], sd = sigma)

# MLE estimate
hat_mui <- apply(dat, 1, mean)
hat_sigma <- sqrt(sum(sweep(dat, 1, hat_mui)^2)/(n*T)) 

# indirect inference procedure
sim_dat <- function(phi, unitfe){
  # simulate data
  dat <- matrix(0, nrow = n, ncol = T)
  for (i in 1:n) dat[i, ] <- rnorm(T, mean = unitfe[i], sd = phi)
  return(dat)
}

objective <- function(phi, H, alpha_hat, coef) {
  coef_sim <- 0
  for (h in 1:H) {
    dats <- sim_dat(phi, alpha_hat)
    s_mui <- apply(dats, 1, mean)
    # MLE estimate
    s_sigma <- sqrt(sum(sweep(dats, 1, s_mui)^2)/(n*T)) 
    coef_sim <- coef_sim + s_sigma/H
  }
  return( sum( (coef - coef_sim)^2 ) )
}

est <- optimize(objective, c(0, 3), H, hat_mui, hat_sigma)
result <- list(true = sigma, mle = hat_sigma, ii = est$minimum, H = H)
print(result)
