
library(mvtnorm)
library(fields)
library(sp)

source("c:/gitProject/SpaTempoDA/HW2/helper.fns.R")


## Simulate data
set.seed(3)
n <- 200
x <- matrix(runif(n*2, min=0, max=1), ncol = 2) # locations
X <- cbind(1, x) # trend surface model
head(X)

beta <- c(0, 1, 1); sigma2 <- 0.1; rho <- 0.1 # true parameters

d <- rdist(x)
Sigma <- sigma2 * exp(-d/rho) # Exponential covariance model
Y <- drop(rmvnorm(1, mean = X %*% beta, sigma = Sigma))


plot.point.ref(x, Y)



## Profile log-likelihood, ignoring a constant

pll <- function(rho, d, Y, Xmat, verbose = FALSE){
  n <- length(Y)
  K <- exp(-d/rho)
  beta.hat <- solve(t(Xmat) %*% solve(K, Xmat)) %*% t(Xmat) %*% solve(K, Y)
  resid <- drop(Y - Xmat %*% beta.hat)
  sigma2.hat <- t(resid) %*% solve(K, resid) / n
  p <- - 0.5 * n * log(sigma2.hat) - 0.5 * determinant(K, log = TRUE)$modulus
  if (verbose) print(c(rho, pll = p))
  return(p)
}

## below blocks are my code

#1a
rho_profile_log_liklihood_graph_val = rep(0, 1000)
rho_candidate = seq(from=0.005, to=0.5, length.out=1000)
for(i in 1:1000){
  rho_profile_log_liklihood_graph_val[i]= pll(rho_candidate[i], d=d, Y=Y, Xmat=X)
}
plot(rho_candidate, rho_profile_log_liklihood_graph_val, type='l')
rho_critical_pt = optimize(pll, lower=0.005, upper=0.5, maximum=TRUE, d=d, Y=Y, Xmat=X)
lines(c(0.005,0.5), c(rho_critical_pt$objective, rho_critical_pt$objective), col="red")
lines(c(rho_critical_pt$maximum,rho_critical_pt$maximum), c(0, rho_critical_pt$objective), col="red")
critical_point_text = paste("(", round(rho_critical_pt$maximum,5),", ", round(rho_critical_pt$objective,3),")", sep="")
text(rho_critical_pt$maximum, rho_critical_pt$objective, labels = critical_point_text)
dev.copy(png,filename="C:/gitProject/SpaTempoDA/HW2/prob1_rho_profile_log_likelihood.png")
dev.off()

#1b
pll_beta_sigma2 = function(rho, d, Y, Xmat, verbose = FALSE){
  n <- length(Y)
  K <- exp(-d/rho)
  beta.hat <- solve(t(Xmat) %*% solve(K, Xmat)) %*% t(Xmat) %*% solve(K, Y)
  resid <- drop(Y - Xmat %*% beta.hat)
  sigma2.hat <- t(resid) %*% solve(K, resid) / n
  if (verbose) {
    print(c(beta.hat=beta.hat))
    print(c(sigma2.hat=sigma2.hat))
  }
  return (c(beta.hat, sigma2.hat))
}

beta_sigma2_vec = pll_beta_sigma2(rho_critical_pt$maximum, d=d, Y=Y, Xmat=X, verbose=TRUE)

