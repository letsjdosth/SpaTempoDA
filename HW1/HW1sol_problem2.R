# R code created by Choi Seokjun
# UTF-8 CRLF
# windows 10 64-bit.
# R 3.5.2

# Note: many parts of this code is reconstructed from Park's code.

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

## function for plotting
ploteqc = function(spobj, z, breaks, ...){
  pal = tim.colors(length(breaks)-1)
  fb = classIntervals(z, n = length(pal), 
                       style = "fixed", fixedBreaks = breaks)
  col = findColours(fb, pal)
  plot(spobj, col = col, ...)
  image.plot(legend.only = TRUE, zlim = range(breaks), col = pal)
}


## Data loading and EDA

load("C:/gitProject/SpaTempoDA/HW1/CAtemps.RData")
# head(CAtemp)
# typeof(CAtemp) #S4
# names(CAtemp) #[1] "avgtemp"   "elevation" ## lat, lon, avgtemp, elevation

# Plotting for EDA
# range(CAtemp$avgtemp)
breaks = 30:75 #<-see range and set by your hand!
x11()
ploteqc(CAtemp, CAtemp$avgtemp, breaks, pch = 19)
map("county", region = "california", add = TRUE)
title(main = "Average Annual Temperatures, 1961-1990, Degrees F")
dev.copy(png,filename="C:/gitProject/SpaTempoDA/HW1/prob2_CAtemp_avgtemp.png")
dev.off()



## (a) Using the CAtemp data, form a preliminary estimate of β using ordinary least squares and make a color plot of the residuals. 
# Include your estimates and plot.

ols = lm(avgtemp~lon+lat+elevation, data=CAtemp)
print("ols coefficients")
print(ols$coeff) #estimates
print('----------------')
# x11()
# plot(ols$residual) # r ordinary plot
# abline(0,0)

CAtemp$ols.residual = ols$residual

# Plotting
# range(CAtemp$ols.residual)
breaks = seq(-7, 7, by = 0.7)
x11()
ploteqc(CAtemp, CAtemp$ols.residual, breaks, pch = 19)
map("county", region = "california", add = TRUE)
title(main = "ols residual")
dev.copy(png,filename="C:/gitProject/SpaTempoDA/HW1/prob2_CAtemp_ols_residual.png")
dev.off()



## (b) Estimate the variogram nonparametrically and then fit the exponential variogram to it using weighted least squares. 
# Make and include a plot of the nonparametric and parametric variogram functions. 
# Also store your parameter estimates and report them.

#  Nonparametric estimation of the variogram
vg = variogram(ols.residual ~ 1, data = CAtemp, width=50) # What happens when we change the bin width?
print(plot(vg, xlab = "Distance", ylab = "Semi-variogram estimate"))
dev.copy(png,filename="C:/gitProject/SpaTempoDA/HW1/prob2_CAtemp_variogram.png")
dev.off()

vgangle = variogram(ols.residual ~ 1, data = CAtemp, alpha = c(0, 45, 90, 135))
x11()
plot.new()
print(plot(vgangle, xlab = "Distance", ylab = "Semi-variogram estimate"))
dev.copy(png,filename="C:/gitProject/SpaTempoDA/HW1/prob2_CAtemp_variogram_angle.png")
dev.off()

#  Modeling the semivariogram
x11()
plot.new()
print(show.vgms())
dev.copy(png,filename="C:/gitProject/SpaTempoDA/HW1/prob2_CAtemp_variogram_various_fit.png")
dev.off()


fitvg = fit.variogram(vg, vgm(1, "Exp", range=300, nugget=3)) # second argument has starting values
# print(fitvg)
# estimated parameters of variance term
s2.hat = fitvg$psill[2]
rho.hat = fitvg$range[2]
tau2.hat = fitvg$psill[1]
x11()
plot.new()
print(plot(vg, fitvg, xlab = "Distance", ylab = "Semi-variogram estimate"))
dev.copy(png,filename="C:/gitProject/SpaTempoDA/HW1/prob2_CAtemp_variogram_exp_fit.png")
dev.off()



## (c) We will now form the GLS estimate of β by hand, rather than using the gls function. 
# (This function doesn’t handle longitude and latitude well, and I also want to give you some practice with matrix calculations in R.) 
# - Use the rdist function in fields to create a matrix of distances (in miles) between pairs of locations in CAtemp. 
# - Create the covariance matrix, plugging in your estimates from the fitted variogram. 
# (Hint: Sum two matrices, one without a nugget and one using the diag function to create the matrix τ^2I.) 
# - Invert the covariance matrix and store it for later reference. 
# - Create the X matrix. (Hint: Use cbind.) 
# - Put all the pieces together to formb βGLS. 

# lec1.2-40page
# lec2.1-24page

make_chordaldist_u1u2_for_rdist = function(lat, lon){
  mean_earth_radius_in_miles = 3958.8 #by google
  chordaldist_x = mean_earth_radius_in_miles * cos(lat) * cos(lon)
  chordaldist_y = mean_earth_radius_in_miles * cos(lat) * sin(lon)
  chordaldist_z = mean_earth_radius_in_miles * sin(lat)
  chordaldist_mat = cbind(chordaldist_x, chordaldist_y, chordaldist_z)
  return(chordaldist_mat)
}

make_EXPcov_mat_without_nugget = function(dist_mat, rho, sigma_sqaure){
  #lecture note 1-2. page 40
  num_row = dim(dist_mat)[1] 
  num_col = dim(dist_mat)[2]
  cov_gamma = matrix(0, num_row, num_col)
  for(i in 1:num_row){
    for(j in 1:num_col){
        cov_gamma[i,j] = exp(-dist_mat[i,j] * rho)
    }
  }
  cov_gamma = cov_gamma * sigma_sqaure
  return(cov_gamma)
}

# make distance matrix
CAtemp_data_u1u2= make_chordaldist_u1u2_for_rdist(CAtemp$lat, CAtemp$lon)
CAtemp_data_dist_mat = rdist(CAtemp_data_u1u2, CAtemp_data_u1u2)
# dim(dist_mat) # 200 200
# max(dist_mat)
# min(dist_mat)


# make covariance matrix
CAtemp_data_cov_spatial = make_EXPcov_mat_without_nugget(CAtemp_data_dist_mat, rho.hat, s2.hat)
CAtemp_data_cov_nugget = diag(200) * tau2.hat
CAtemp_data_cov_mat = CAtemp_data_cov_spatial + CAtemp_data_cov_nugget


# gls fit
b0 = rep(1,200)
gls.X = as.matrix(data.frame(b0, CAtemp$lon, CAtemp$lat, CAtemp$elevation))
gls.Y = as.matrix(CAtemp$avgtemp)
CAtemp_data_inv_cov_mat = solve(CAtemp_data_cov_mat)
gls.beta = solve(t(gls.X) %*% CAtemp_data_inv_cov_mat %*% gls.X) %*% t(gls.X) %*% CAtemp_data_inv_cov_mat %*% gls.Y
print("gls coefficients")
cat('(Intercept)    lon       lat    elevation\n',gls.beta,'\n')
print('----------------')




## (d) Calculate and plot the EBLUP of µ + Z at the locations in CAgrid, plugging in your estimates from (b) and (c). 
# Calculate and plot the (estimated) standard error of Z at each prediction location.
#lec2.2 page38

# head(CAgrid)
# dim(CAgrid) #664 1

# make m0(=pred_mu) term(in lec2-2. page38)
pred_X = as.matrix(data.frame(rep(1,664), CAgrid$lon, CAgrid$lat, CAgrid$elevation))
pred_mu = pred_X %*% gls.beta
dim(pred_mu)
head(pred_mu)


# make covariance matrix
#inner
CAgrid_data_u1u2= make_chordaldist_u1u2_for_rdist(CAgrid$lat, CAgrid$lon)
pred_inner_dist_mat = rdist(CAgrid_data_u1u2, CAgrid_data_u1u2)

pred_cov_spatial = make_EXPcov_mat_without_nugget(pred_inner_dist_mat, rho.hat, s2.hat)
# pred_cov_nugget = diag(664) * tau2.hat #no nugget term!
pred_inner_cov_mat = pred_cov_spatial #+ pred_cov_nugget #no nugget term!


#cross
pred_cross_dist_mat = rdist(CAtemp_data_u1u2, CAgrid_data_u1u2)
dim(pred_cross_dist_mat) # 200, 664
pred_cross_cov_mat = make_EXPcov_mat_without_nugget(pred_cross_dist_mat, rho.hat, s2.hat)


## krigging mean
pred_Y = pred_mu + t(pred_cross_cov_mat) %*% CAtemp_data_inv_cov_mat %*% (gls.Y - gls.X%*%gls.beta)
dim(pred_Y)
CAgrid$pred.temp = pred_Y

# sd (I'll skip this. Instead, I'll find mse.)
pred_Y_var_mat = pred_inner_cov_mat - t(pred_cross_cov_mat) %*% CAtemp_data_inv_cov_mat %*% pred_cross_cov_mat
# dim(pred_Y_var_mat) #664 664
# pred_Y_var = diag(pred_Y_var_mat)
# CAgrid$pred.sd.temp = sqrt(pred_Y_var)

# predicting mse
cal_b = t(pred_X) - t(gls.X) %*% CAtemp_data_inv_cov_mat %*% pred_cross_cov_mat
pred_mse_mat = pred_Y_var_mat + t(cal_b) %*% solve(t(gls.X) %*% CAtemp_data_inv_cov_mat %*% gls.X) %*% cal_b
pred_mse = diag(pred_mse_mat)
CAgrid$pred.mse.temp = pred_mse


## Plotting
# mean
range(pred_Y)
breaks = 30:75
x11()
plot.new()

ploteqc(CAgrid, CAgrid$pred.temp, breaks, pch = 19)
map("county", region = "california", add = TRUE)
title(main = "Predicted Average Annual Temperatures, Degrees F")
dev.copy(png,filename="C:/gitProject/SpaTempoDA/HW1/prob2_CAgrid_predicted_mean.png")
dev.off()

# # #std
# # quantile(sqrt(pred_Y_var))
# # breaks = seq(1.5, 2.15, length.out=40)
# # x11()
# # plot.new()

# # ploteqc(CAgrid, CAgrid$pred.sd.temp, breaks, pch = 19)
# # map("county", region = "california", add = TRUE)
# # title(main = "Std of Predicted Average Annual Temperatures, 1961-1990, Degrees F")


# mse
quantile(pred_mse)
breaks = seq(4.6, 5, length.out=50)
x11()
plot.new()

ploteqc(CAgrid, CAgrid$pred.mse.temp, breaks, pch = 19)
map("county", region = "california", add = TRUE)
title(main = "MSE of Predicted Temperatures")
dev.copy(png,filename="C:/gitProject/SpaTempoDA/HW1/prob2_CAgrid_predicted_MSE.png")
dev.off()
