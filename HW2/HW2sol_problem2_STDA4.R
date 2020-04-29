
library(classInt)
library(fields)
library(maps)
library(sp)
library(gstat)
library(geoR)
library(mvtnorm)
library(MCMCpack)
library(coda)

#############################
## California temperatures ##
#############################

load("c:/gitProject/SpaTempoDA/HW2/CAtemps.RData")

ploteqc <- function(spobj, z, breaks, ...){
  pal <- tim.colors(length(breaks)-1)
  fb <- classIntervals(z, n = length(pal), 
                       style = "fixed", fixedBreaks = breaks)
  col <- findColours(fb, pal)
  plot(spobj, col = col, ...)
  image.plot(legend.only = TRUE, zlim = range(breaks), col = pal)
}

## Plotting

range(CAtemp$avgtemp)
breaks <- seq(40, 75, by = 5)
ploteqc(CAtemp, CAtemp$avgtemp, breaks, pch = 19)
map("county", region = "california", add = TRUE)
title(main = "Average Annual Temperatures, 1961-1990, Degrees F")

range(CAgrid$elevation)
breaks <- seq(-100, 3600, by = 100)
ploteqc(CAgrid, CAgrid$elevation, breaks, pch = 19)
map("county", region = "california", add = TRUE)
title(main = "Elevations at prediction locations, m")

## Preliminary model fitting

linmod <- lm(avgtemp~lon+lat+elevation, data = CAtemp)
summary(linmod)
CAtemp$resid <- linmod$resid

range(CAtemp$resid)
breaks <- seq(-7, 7, by = 1)
ploteqc(CAtemp, CAtemp$resid, breaks, pch = 19)
map("county", region = "california", add = TRUE)
title(main = "Residuals")

vg <- variogram(resid ~ 1, data = CAtemp)
plot(vg, xlab = "Distance", ylab = "Semi-variogram estimate")
fitvg <- fit.variogram(vg, 
                       vgm(1, "Exp", range = 200, nugget = 1), 
                       fit.method = 2)
plot(vg, fitvg, xlab = "Distance", 
     ylab = "Semi-variogram estimate")

## Prior parameters

m.beta <- rep(0, 4); V.beta <- 100000 * diag(4)
a.s2 <- 0.001; b.s2 <- 0.001
a.t2 <- 0.001; b.t2 <- 0.001

rhoseq <- seq(0.01, 300, length = 100)
plot(rhoseq, dgamma(rhoseq, shape = 1, scale = 1)) # old prior for rho
m.rho <- 100; v.rho <- 5000
b.rho <- v.rho/m.rho; a.rho <- m.rho/b.rho
plot(rhoseq, dgamma(rhoseq, shape = a.rho, scale = b.rho), type = "l") # new prior for rho

## Setup, storage, and starting values

y <- CAtemp$avgtemp
n <- nrow(CAtemp); m <- nrow(CAgrid)
d <- rdist.earth(coordinates(CAtemp))
X <- cbind(rep(1, n), CAtemp$lon, CAtemp$lat, CAtemp$elevation)
Xpred <- cbind(rep(1, m), CAgrid$lon, CAgrid$lat, CAgrid$elevation)

B <- 1000

beta.samps <- matrix(NA, nrow = 4, ncol = B)
beta.samps[,1] <- coef(linmod)

s2.samps <- t2.samps <- rho.samps <- rep(NA, B)
s2.samps[1] <- fitvg$psill[2]
rho.samps[1] <- fitvg$range[2]
t2.samps[1] <- fitvg$psill[1]

eta.obs.samps <- matrix(NA, nrow = n, ncol = B)

v.prop <- 100^2

## MCMC sampler

Gamma <- exp(-d/rho.samps[1]) # initalize Gamma matrix
Ginv <- solve(Gamma)

for(i in 2:B){
  
  if(i%%100==0) print(i)
  
  ## eta_obs | Rest
  V <- solve(diag(n)/t2.samps[i-1] + Ginv/s2.samps[i-1])
  m <- V %*% (y/t2.samps[i-1] + Ginv %*% X %*% 
                beta.samps[,i-1] / s2.samps[i-1])
  eta.obs.samps[,i] <- rmvnorm(1, mean = m, sigma = V, method = "svd")
  
  ## beta | Rest
  V <- solve(t(X) %*% Ginv %*% X / s2.samps[i-1] + solve(V.beta))
  m <- V %*% (t(X) %*% Ginv %*% eta.obs.samps[,i] / 
                s2.samps[i-1] + solve(V.beta, m.beta))
  beta.samps[,i] <- rmvnorm(1, mean = m, sigma = V, method = "svd")
  
  ## s2 | Rest
  a <- a.s2 + n/2
  resid <- eta.obs.samps[,i] - X %*% beta.samps[,i]
  b <- b.s2 + t(resid) %*% Ginv %*% resid /2
  s2.samps[i] <- rinvgamma(1, a, b)
  
  ## t2 | Rest
  a <- a.t2 + n/2
  resid <- y - eta.obs.samps[,i]
  b <- b.t2 + t(resid) %*% resid / 2
  t2.samps[i] <- rinvgamma(1, a, b)
  
  ## rho | Rest 
  
  # Visualize posterior surface
  # The ratio of this function at rho.cand to rho.samps[i-1] is what determines r
  if(FALSE){
    prho <- sapply(rhoseq, function(rho){
      dmvnorm(eta.obs.samps[,i], mean = X %*% beta.samps[,i], 
              sigma = s2.samps[i] * exp(-d/rho), log = TRUE) +
        dgamma(rho, shape = a.rho, scale = b.rho, log = TRUE)})
    plot(rhoseq, exp(prho), type = "l")
  }
  
  rho.cand <- rnorm(1, mean = rho.samps[i-1], sd = sqrt(v.prop))
  if(rho.cand < 0){ # automatically reject
    rho.samps[i] <- rho.samps[i-1]
  } else {
    lik1 <- dmvnorm(eta.obs.samps[,i], mean = X %*% beta.samps[,i],
                    sigma = s2.samps[i] * exp(-d/rho.cand), log = TRUE)
    lik2 <- dmvnorm(eta.obs.samps[,i], mean = X %*% beta.samps[,i],
                    sigma = s2.samps[i] * exp(-d/rho.samps[i-1]), log = TRUE)
    p1 <- dgamma(rho.cand, shape = a.rho, scale = b.rho, log = TRUE)
    p2 <-   dgamma(rho.samps[i-1], shape = a.rho, scale = b.rho, log = TRUE)
    r <- exp(lik1 + p1 - lik2 - p2)
    if(runif(1) < r){ # accept
      rho.samps[i] <- rho.cand
      Gamma <- exp(-d/rho.cand) 
      Ginv <- solve(Gamma)
    } else { # reject
      rho.samps[i] <- rho.samps[i-1]
    }
  }
  
}

## Diagnostics

plot(beta.samps[1,], type = "l")
plot(s2.samps, type = "l")
plot(rho.samps, type = "l")
plot(t2.samps, type = "l")
plot(eta.obs.samps[1,], type = "l")

length(unique(rho.samps[1:B]))/B # acc rate

burnin <- 100
s2.final <- s2.samps[-(1:burnin)]
t2.final <- t2.samps[-(1:burnin)]
beta.final <- beta.samps[,-(1:burnin)]
eta.obs.final <- eta.obs.samps[,-(1:burnin)]
rho.final <- rho.samps[-(1:burnin)]

acf(s2.final)
acf(t2.final)
acf(beta.final[1,])
acf(eta.obs.final[1,])
acf(rho.final)

effectiveSize(s2.final)
effectiveSize(t2.final)
effectiveSize(beta.final[1,])
effectiveSize(eta.obs.final[1,])
effectiveSize(rho.final)


tie_plot = function(final, param_name){
  dev.new(width=15, height=7.5, unit="cm")
  par(mfrow=c(1,2), pty="s")
  plot(final, type="l")
  title(paste("trace plot of", param_name, sep=" "))
  acf(final, main="acf")
  text(30,0.9,labels=paste("ESS:", effectiveSize(final)),pos=2)
  
  file_name = paste("C:/gitProject/SpaTempoDA/HW2/prob2_",param_name,".png",sep="")
  dev.copy(png,filename=file_name)
  dev.off()
}
tie_plot(s2.final, "sigma_square")
tie_plot(t2.final, "tau_square")
tie_plot(beta.final[1,], "beta0")
tie_plot(beta.final[2,], "beta1")
tie_plot(beta.final[3,], "beta2")
tie_plot(beta.final[4,], "beta3")
tie_plot(eta.obs.final[1,],"eta(first_dim)")
tie_plot(rho.final,"rho")


## Prediction

dcross <- rdist.earth(coordinates(CAtemp), coordinates(CAgrid))
dpred <- rdist.earth(coordinates(CAgrid))

index <- seq(1, B-burnin, by = 20) # which samples to use (thinning)
eta.pred <- matrix(NA, nrow = nrow(CAgrid), ncol = length(index))

for(i in 1:length(index)){
  print(i)
  j <- index[i]
  
  # Construct the covariance matrices
  Gamma <- exp(-d/rho.samps[j]) 
  Ginv <- solve(Gamma)
  g <- exp(-dcross/rho.samps[j])
  Gpred <- exp(-dpred/rho.samps[j])
  
  m <- Xpred %*% beta.final[,j] + t(g) %*% Ginv %*% 
    (y - X %*% beta.final[,j])
  V <- s2.final[j] * (Gpred - t(g)%*%Ginv%*%g)
  eta.pred[,i] <- rmvnorm(1, m, V, method = "svd")
}

## Find pointwise posterior means and sds

eta.pred.m <- apply(eta.pred, 1, mean)
eta.pred.sd <- apply(eta.pred, 1, sd)

range(eta.pred.m)
breaks <- seq(30, 80, by = 5)

ploteqc(CAgrid, eta.pred.m, breaks, pch = 19)
map("county", region = "california", add = TRUE)
title(main = "Posterior Mean")
dev.copy(png,filename="C:/gitProject/SpaTempoDA/HW2/prob2_posterior.png")
dev.off()


range(eta.pred.sd)
breaks <- seq(0, 3.4, by = 0.2)

ploteqc(CAgrid, eta.pred.sd, breaks, pch = 19)
map("county", region = "california", add = TRUE)
points(CAtemp)
title(main = "Posterior Standard Deviation")
dev.copy(png,filename="C:/gitProject/SpaTempoDA/HW2/prob2_posterior_SE.png")
dev.off()

