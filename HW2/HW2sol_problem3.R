# R code created by Choi Seokjun
# UTF-8 CRLF
# windows 10 64-bit.
# R 3.6.3


library(nimble)

library(classInt) # breakpoints for color plotting
library(mvtnorm)
library(fields)
library(sp)
library(MASS)


plot.point.ref <- function(spatialdata, vals) {
  pal <- c("#0000FF","#FF0000")
  ints <- classIntervals(vals, n = 2, style = "pretty") # Determine breakpoints
  # also see style options "quantile" and "fisher"
  intcols <- findColours(ints, pal) # vector of colors
  # if pal doesn't have the same length as # classes, findColours will interpolate

  par(mar = rep(0, 4))
  plot(spatialdata, col = intcols, pch = 19)
  points(spatialdata, pch = 1)
  legend("topleft", fill = attr(intcols, "palette"),
         legend = names(attr(intcols, "table")), bty = "n")
}



## Simulate data
set.seed(2019311252)
n = 400
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
        Z[i] ~ dbern(prob[i])
    }
    #Constant and Cov Matrix
    XB[1:n] <- beta0*X[,1] + beta1*X[,2] + beta2*X[,3]
    covMat[1:n,1:n] <- expcov(dists[1:n,1:n],rho)
    fullCovMat[1:n,1:n] <- sigma2 * covMat[1:n,1:n]

    #Process Model
    W[1:n] ~ dmnorm(mean=mn[1:n], cov=fullCovMat[1:n,1:n])

    #Parameter Models
    sigma2 ~ dinvgamma(0.2, 0.2)
    rho ~ dunif(0,1)
    beta0 ~ dnorm(0, sd=sqrt(100))
    beta1 ~ dnorm(0, sd=sqrt(100))
    beta2 ~ dnorm(0, sd=sqrt(100))
})

niter = 200000
consts = list(n=n, X=sim_X, dists=dist_mat, mn=rep(0,n))
data = list(Z=sim_Y)
inits = list(beta0 = rnorm(1), beta1=rnorm(1), beta2=rnorm(1), rho=0.5, sigma2=2, W=rnorm(n))

pt = proc.time()
samples = nimbleMCMC(model_string, data=data, inits=inits, constants=consts,
    monitors=c("beta1","beta2","rho","sigma2","W"),
    samplesAsCodaMCMC=TRUE, WAIC=FALSE, summary=FALSE,
    niter=niter, nburnin=0, nchains=1)
ptFinal = proc.time()-pt







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

accRateFunc<-function(x){
  accRate<-(length(unique(x))-1)/(length(x)-1)
  return(accRate)
}

# Summary 
summaryFunction<-function(mcmcDat){
  # Parameters
  summaryMat<-rbind(apply(mcmcDat,2,mean),
                          apply(mcmcDat,2,hpd),
                          apply(mcmcDat,2,accRateFunc))
  
  rownames(summaryMat)<-c("Mean","95%CI-Low","95%CI-High",
                                "Accept")
  return(summaryMat)
}

plotFunction = function(mcmcDat, trueParamVal, paramName){
  dev.new(width=15, height=7.5, unit="cm")

  par(mfrow=c(1,3))
  MCsample_mean = mean(mcmcDat)
  MCsample_hpd = hpd(mcmcDat)

  plot.ts(mcmcDat)
  abline(h=trueParamVal, col="red", lwd=2)
  abline(h=MCsample_mean, col="blue", lwd=2)
  title(paste("traceplot of", paramName))
  text(0,trueParamVal,label="TRUE", col="red", pos=3)
  text(0,MCsample_mean,label="MCMC", col="blue", pos=3)


  plot(density(mcmcDat), main=paste("emp. density of", paramName))
  abline(v=trueParamVal, col="red", lwd=2)
  abline(v=MCsample_mean, col="blue", lwd=2)
  text(trueParamVal,0.001,label="TRUE",col="red",pos=3)
  text(MCsample_mean,0.001,label="MCMC",col="blue",pos=3)
  lines(MCsample_hpd,c(0.001,0.001), col="blue",lwd=3)
  text(MCsample_hpd[1],0.001,label="95%HPD",col="blue",pos=3)

  acf(mcmcDat, main="acf")

  plot_filename=paste("C:/gitProject/SpaTempoDA/HW2/prob3_",paramName,".png",sep="")
  dev.copy(png,filename=plot_filename)
  dev.off()
}



SummaryMat = list()
SummaryMat[[1]] = summaryFunction(samples[,c("beta1","beta2","sigma2","rho")])
SummaryMat[[2]] = summaryFunction(samples[,1:n])
# save(SummaryMat,samples,ptFinal,file="C:/gitProject/SpaTempoDA/HW2/prob2_MCMCresult.RData")
print(SummaryMat[[1]])

# beta = c(0, 1, 1); sigma2 = 0.1; rho = 0.1 # true parameters
plotFunction(samples[10000:200000,"sigma2"], 0.1, "sigma2")
plotFunction(samples[10000:200000,"beta1"], 1, "beta1")
plotFunction(samples[10000:200000,"beta2"], 1, "beta2")
plotFunction(samples[10000:200000,"rho"], 0.1, "rho")

