# 2. The file CAtemps.RData contains two R objects of class SpatialPointsDataFrame, called CAtemp and CAgrid. 
# CAtemp contains average temperatures from 1961-1990 at 200 locations (latitude and longitude) in California in degrees Fahrenheit, 
# along with their elevations in meters. CAgrid contains elevations in meters over a grid of locations. 
# I’ve given you some code to get started with this data in HW1.R.

# Consider the following model for the temperature data.
# Yi = µ(si;β) + e(si;σ2,ρ,τ)
# where µ(si,β) = β0+β1Longitude(s)+β2Latitude(s)+β3Elevation(s) and e(si;σ2,ρ,τ) 
# is a zero mean stationary Gaussian process with exponential covariance function.

# Another way of writing this is as
# Yi = µ(si;β) + e(si;σ2,ρ) + ei
# where now Z is a mean zero Gaussian process like e but without the nugget term, and the ei are iid N(0,τ2), independent of Z. 
# This is important because we want to predict µ(si;beta) + Z(si;σ2,ρ) without the measurement error.

library(sp)
library(gstat)
library(fields)
library(classInt)
library(maps)

ploteqc = function(spobj, z, breaks, ...){
  pal = tim.colors(length(breaks)-1)
  fb = classIntervals(z, n = length(pal), 
                       style = "fixed", fixedBreaks = breaks)
  col = findColours(fb, pal)
  plot(spobj, col = col, ...)
  image.plot(legend.only = TRUE, zlim = range(breaks), col = pal)
}

load("C:/gitProject/SpaTempoDA/HW1/CAtemps.RData")
# head(CAtemp)
# typeof(CAtemp) #S4
# names(CAtemp) #[1] "avgtemp"   "elevation"
## lat, lon, avgtemp, elevation

## Plotting
# range(CAtemp$avgtemp)
breaks = 30:75
x11()
ploteqc(CAtemp, CAtemp$avgtemp, breaks, pch = 19)
map("county", region = "california", add = TRUE)
title(main = "Average Annual Temperatures, 1961-1990, Degrees F")



# (a) Using the CAtemp data, form a preliminary estimate of β using ordinary least squares and make a color plot of the residuals. 
# Include your estimates and plot.

ols = lm(avgtemp~lon+lat+elevation, data=CAtemp)
summary(ols)
# x11()
# plot(ols$residual)
abline(0,0)

CAtemp$ols.residual = ols$residual

## Plotting
# range(CAtemp$ols.residual)
breaks = seq(-7, 7, by = 0.7)
x11()
ploteqc(CAtemp, CAtemp$ols.residual, breaks, pch = 19)
map("county", region = "california", add = TRUE)
title(main = "ols residual")




# (b) Estimate the variogram nonparametrically and then fit the exponential variogram to it using weighted least squares. 
# Make and include a plot of the nonparametric and parametric variogram functions. 
# Also store your parameter estimates and report them.


##  Nonparametric estimation of the variogram

vg = variogram(ols.residual ~ 1, data = CAtemp, width=50) # What happens when we change the bin width?
x11()
plot(vg, xlab = "Distance", ylab = "Semi-variogram estimate", width=5) 


vgangle = variogram(ols.residual ~ 1, data = CAtemp, alpha = c(0, 45, 90, 135))
x11()
plot(vgangle, xlab = "Distance", ylab = "Semi-variogram estimate")

##  Modeling the semivariogram

# show.vgms()

# second argument has starting values
fitvg = fit.variogram(vg, vgm(1, "Exp", range=300, nugget=3))
# print(fitvg)
s2.hat = fitvg$psill[2]
rho.hat = fitvg$range[2]
tau2.hat = fitvg$psill[1]
x11()
plot(vg, fitvg, xlab = "Distance", ylab = "Semi-variogram estimate")


# (c) We will now form the GLS estimate of β by hand, rather than using the gls function. 
# (This function doesn’t handle longitude and latitude well, and I also want to give you some practice with matrix calculations in R.) 
# - Use the rdist function in fields to create a matrix of distances (in miles) between pairs of locations in CAtemp. 
# - Create the covariance matrix, plugging in your estimates from the fitted variogram. 
# (Hint: Sum two matrices, one without a nugget and one using the diag function to create the matrix τ^2I.) 
# - Invert the covariance matrix and store it for later reference. 
# - Create the X matrix. (Hint: Use cbind.) 
# - Put all the pieces together to formb βGLS. 

# lec1.2-40page
# lec2.1-24page

mean_earth_radius_in_miles = 3958.8 #by google
chordaldist_x = mean_earth_radius_in_miles * cos(CAtemp$lat) * cos(CAtemp$lon)
chordaldist_y = mean_earth_radius_in_miles * cos(CAtemp$lat) * sin(CAtemp$lon)
chordaldist_z = mean_earth_radius_in_miles * sin(CAtemp$lat)
chordaldist_mat = cbind(chordaldist_x, chordaldist_y, chordaldist_z)

dist_mat = rdist(chordaldist_mat, chordaldist_mat)
# dim(dist_mat)
# max(dist_mat)
# min(dist_mat)

cov_big_gamma = diag(200)
for(i in 1:200){
  for(j in 1:200){
     if(i!=j){
      cov_big_gamma[i,j] = exp(-dist_mat[i,j] * rho.hat)
    }
  }
}
cov_spatial = s2.hat*cov_big_gamma
cov_nugget = diag(200) * tau2.hat
cov_mat = cov_spatial + cov_nugget

b0 = rep(1,200)
gls.X = as.matrix(data.frame(b0, CAtemp$lon, CAtemp$lat, CAtemp$elevation))
gls.Y = as.matrix(CAtemp$avgtemp)
inv_cov_mat = solve(cov_mat)
gls.beta = solve(t(gls.X) %*% inv_cov_mat %*% gls.X) %*% t(gls.X) %*% inv_cov_mat %*% gls.Y
(gls.beta)

# (d) Calculate and plot the EBLUP of µ + Z at the locations in CAgrid, plugging in your estimates from (b) and (c). 
# Calculate and plot the (estimated) standard error of Z at each prediction location.
#lec2.2 page38

# head(CAgrid)
# dim(CAgrid) #664 1

pred_X = as.matrix(data.frame(rep(1,664), CAgrid$lon, CAgrid$lat, CAgrid$elevation))
pred_mu = pred_X %*% gls.beta
dim(pred_mu)
head(pred_mu)


#cov
#inner
mean_earth_radius_in_miles = 3958.8 #by google
pred_chordaldist_x = mean_earth_radius_in_miles * cos(CAgrid$lat) * cos(CAgrid$lon)
pred_chordaldist_y = mean_earth_radius_in_miles * cos(CAgrid$lat) * sin(CAgrid$lon)
pred_chordaldist_z = mean_earth_radius_in_miles * sin(CAgrid$lat)
pred_chordaldist_mat = cbind(pred_chordaldist_x, pred_chordaldist_y, pred_chordaldist_z)

pred_inner_dist_mat = rdist(pred_chordaldist_mat, pred_chordaldist_mat)

pred_cov_big_gamma = diag(664)
for(i in 1:664){
  for(j in 1:664){
     if(i!=j){
      pred_cov_big_gamma[i,j] = exp(-pred_inner_dist_mat[i,j] * rho.hat)
    }
  }
}
pred_cov_spatial = s2.hat*pred_cov_big_gamma
pred_cov_nugget = diag(664) * tau2.hat
pred_inner_cov_mat = pred_cov_spatial + pred_cov_nugget


#cross
pred_cross_dist_mat = rdist(chordaldist_mat, pred_chordaldist_mat)
dim(pred_cross_dist_mat) # 200, 664
pred_cov_cross_gamma = matrix(0, 200, 664)
for(i in 1:200){
  for(j in 1:664){
      pred_cov_cross_gamma[i,j] = exp(-pred_cross_dist_mat[i,j] * rho.hat)
  }
}
pred_cross_cov_mat = pred_cov_cross_gamma * s2.hat



#krigging
pred_Y = pred_mu + t(pred_cross_cov_mat) %*% inv_cov_mat %*% (gls.Y - gls.X%*%gls.beta)
dim(pred_Y)
CAgrid$pred.temp = pred_Y

# Plotting
range(pred_Y)
breaks = 30:75
x11()
plot.new()

ploteqc(CAgrid, CAgrid$pred.temp, breaks, pch = 19)
map("county", region = "california", add = TRUE)
title(main = "Predicted Average Annual Temperatures, 1961-1990, Degrees F")


