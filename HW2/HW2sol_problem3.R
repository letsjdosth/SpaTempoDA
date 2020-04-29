# R code created by Choi Seokjun
# UTF-8 CRLF
# windows 10 64-bit.
# R 4.0.0


library(nimble)

library(mvtnorm)
library(fields)
library(sp)
library(MASS)

source("c:/gitProject/SpaTempoDA/HW2/helper.fns.R")


## Simulate data
set.seed(2019311252)
n = 200
sim_x = matrix(runif(n*2, min=0, max=1), ncol = 2) # locations
sim_X = cbind(1, sim_x) # trend surface model

beta = c(0, 1, 1); sigma2 = 0.1; rho = 0.1 # true parameters

dist_mat = rdist(sim_x)
Sigma = sigma2 * exp(-dist_mat/rho) # Exponential covariance model
sim_eta = drop(rmvnorm(1, mean = sim_X %*% beta, sigma = Sigma))
sim_p = exp(sim_eta)/(1+exp(sim_eta))
sim_Y = rep(0, n)
for(i in 1:n){
    sim_Y[i]= rbinom(1,1,sim_p[i])
}

plot.point.ref(sim_x, sim_Y)



expcov <- nimbleFunction(     
  run = function(dists = double(2), phi = double(0)) {
    returnType(double(2))
    n <- dim(dists)[1]
    result <- matrix(nrow = n, ncol = n, init = FALSE)
    for(i in 1:n){
      for(j in 1:n){
        result[i, j] <- exp(-dists[i,j]/phi)
      }
    }
    
    return(result)
  })


# nimble
model_string = nimbleCode({
    
    #Data Model
    for(i in 1:n){
        prob[i] <- exp(W[i]+XB[i])/(1+exp(W[i]+XB[i]))
        Z[i] ~ dbinom(1,prob[i])
    }
    #Constant and Cov Matrix
    XB[1:n] <- beta1*X[,1] + beta2*X[,2]
    # covMat[1:n,1:n] <- expcov(dists[1:n,1:n],phi)
    covMMat[1:n,1:n] <- 1
    fullCovMat[1:n,1:n] <- sigma2 * covMat[1:n,1:n]

    #Process Model
    W[1:n] ~ dmnorm(mean=mn[1:n], cov=fullCovMat[1:n,1:n])

    #Parameter Models
    sigma2 ~ dinvgamma(0.2, 0.2)
    phi ~ dunif(0,1)
    beta1 ~ dnorm(0, sd=sqrt(100))
    beta2 ~ dnorm(0, sd=sqrt(100))
})

niter = 100000
consts = list(n=n, X=sim_X, dists=dist_mat, mn=rep(0,n))
data = list(Z=sim_Y)
inits = list(beta1=rnorm(1), beta2=rnorm(1), phi=0.5, sigma2=2, W=rnorm(n))

samples = nimbleMCMC(model_string, data=data, inits=inits, constants=consts,
    monitors=c("beta1","beta2","phi","sigma2","W"),
    samplesAsCodaMCMC=TRUE, WAIC=FALSE, summary=FALSE,
    niter=niter, nburnin=0, nchains=1)








## Using Ming-Hui Chen's paper in Journal of Computational and Graphical Stats.
hpd <- function(samp,p=0.05){
  ## to find an approximate (1-p)*100% HPD interval from a
  ## given posterior sample vector samp
  
  r <- length(samp)
  samp <- sort(samp)
  rang <- matrix(0,nrow=trunc(p*r),ncol=3)
  dimnames(rang) <- list(NULL,c("low","high","range"))
  for (i in 1:trunc(p*r)) {
    rang[i,1] <- samp[i]
    rang[i,2] <- samp[i+(1-p)*r]
    rang[i,3] <- rang[i,2]-rang[i,1]
  }
  hpd <- rang[order(rang[,3])[1],1:2]
  return(hpd)
}

# Exponential Covariance Function
expCov<-function(distMat,phi){
  exp(-distMat/phi)
}

sqeCov<-function(distMat,phi){
  exp(-0.5*(distMat/phi)^2)
}

matCov<-function(distMat,phi){
  (1+(sqrt(5)*(distMat/phi))+((5*distMat^2)/(3*(phi^2))))*exp(-(sqrt(5)*(distMat/phi)))
}


# Matern Cov Function + Acceptance Rate function
Matern <- function(d, param = c(scale = 1, range = 1, smoothness = 2)) {
  scale <- param[1]
  range <- param[2]
  smoothness <- param[3]
  if (any(d < 0))
    stop("distance argument must be nonnegative")
  d <- d / range
  d[d == 0] <- 1e-10
  rootcon<-sqrt(2*smoothness)
  con <- (2^(smoothness - 1)) * gamma(smoothness)
  con <- 1 / con
  return(scale * con * ((rootcon*d)^smoothness) * besselK(rootcon*d, smoothness))
}
accRateFunc<-function(x){
  accRate<-(length(unique(x))-1)/(length(x)-1)
  return(accRate)
}


# Summary 
summaryFunction<-function(mcmcDat,bmseThresh=0.01,time){
  
  # Parameters
  summaryMat<-rbind(apply(mcmcDat,2,mean),
                          apply(mcmcDat,2,hpd),
                          apply(mcmcDat,2,accRateFunc),
                          bmmat(mcmcDat)[,2],
                          abs(apply(mcmcDat,2,mean))*bmseThresh,
                          apply(mcmcDat,2,ess),
                          apply(mcmcDat,2,ess)/time)
  
  rownames(summaryMat)<-c("Mean","95%CI-Low","95%CI-High",
                                "Accept","BMSE",paste(bmseThresh,"x mean"),
                                "ESS","ESS/sec")
  return(summaryMat)
}
