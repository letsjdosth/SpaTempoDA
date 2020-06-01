library(sp)
library(spdep)
library(classInt)
library(fields)
library(MCMCpack)
library(mvtnorm)

set.seed(2019311252)

## load data
load("C:/gitProject/SpaTempoDA/HW3/munichrents.RData")
head(rents)
names(rents)

## ordinary least square fit (problem 1)
lm.fit = lm(RentPerM2 ~ Year+NoHotWater+NoCentralHeat+NoBathTiles
    +SpecialBathroom+SpecialKitchen+Room2+Room3+Room4+Room5+Room6, data=rents)
summary(lm.fit)
# plot(lm.fit$residual)

# Call:
# lm(formula = RentPerM2 ~ Year + NoHotWater + NoCentralHeat +
#     NoBathTiles + SpecialBathroom + SpecialKitchen + Room2 +
#     Room3 + Room4 + Room5 + Room6, data = rents)

# Residuals:
#     Min      1Q  Median      3Q     Max
# -8.2192 -1.3791 -0.0393  1.3833  9.8943

# Coefficients:
#                   Estimate Std. Error t value Pr(>|t|)
# (Intercept)     -16.012330   4.048974  -3.955 7.93e-05 ***
# Year              0.013331   0.002059   6.473 1.20e-10 ***
# NoHotWater       -1.870505   0.303264  -6.168 8.33e-10 ***
# NoCentralHeat    -1.225761   0.206737  -5.929 3.57e-09 ***
# NoBathTiles      -0.725711   0.123431  -5.879 4.80e-09 ***
# SpecialBathroom   0.660372   0.169562   3.895 0.000102 ***
# SpecialKitchen    1.462887   0.185429   7.889 4.94e-15 ***
# Room2            -1.319372   0.157060  -8.400  < 2e-16 ***
# Room3            -1.886408   0.156998 -12.016  < 2e-16 ***
# Room4            -2.464631   0.193683 -12.725  < 2e-16 ***
# Room5            -2.378158   0.343380  -6.926 5.80e-12 ***
# Room6            -2.482621   0.590678  -4.203 2.75e-05 ***
# ---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

# Residual standard error: 2.132 on 2023 degrees of freedom
# Multiple R-squared:  0.2561,    Adjusted R-squared:  0.2521
# F-statistic: 63.31 on 11 and 2023 DF,  p-value: < 2.2e-16



## plot sp objects (problem 2)
class(districts.sp)
names(districts.sp)
plot(districts.sp, border = "gray", main="grey area : parks")

names(parks.sp)
plot(parks.sp, add=T, col= "gray")
filename = paste("C:/gitProject/SpaTempoDA/HW3/map_parks.png",sep="")
dev.copy(png,filename=filename)
dev.off()


## Creating nb object based on shared boundaries
coords = coordinates(districts.sp)
nb.bound <- poly2nb(districts.sp)

print(nb.bound) # spatial neighbors list object
nb.bound[[1]] # Neighbors; no neighbors indicated by a 0
summary(nb.bound)

plot(districts.sp, border = "gray", main="neighborhood structure")
plot(nb.bound, coords, pch = 19, cex = 0.6, add = TRUE)
filename = paste("C:/gitProject/SpaTempoDA/HW3/map_neighborhood_structure.png",sep="")
dev.copy(png,filename=filename)
dev.off()


## counting observations of each district (problem 3)
dim(H)
dim(rents)
eachdist_obs_num=colSums(H)

pal <- tim.colors(5)
q5 <- classIntervals(eachdist_obs_num, n = 5, style = "quantile")
col <- findColours(q5, pal)
plot(districts.sp, border = "black", col = col, main="number of observations in each district")
plot(parks.sp, add=T, col= "gray")
legend("topright", fill = attr(col, "palette"),
       legend = names(attr(col, "table")),
       bty="n", cex = 0.8, y.intersp = 1.5)
filename = paste("C:/gitProject/SpaTempoDA/HW3/map_obsnum.png",sep="")
dev.copy(png,filename=filename)
dev.off()




## gibbs sampler
# creating model matrix
# help(nb2mat)
model.W = nb2mat(nb.bound, style="B")
model.Dw = diag(x = rowSums(model.W), 380, 380)

model.H = H
model.m = dim(H)[2] #380

# explore data
# dim(X)
# head(X)
# names(rents)

# length(y)
# y==rents$RentPerM2 #all TRUE
model.n = length(y) #2035




# Full conditional distributions
# calculate some expensive matrix
model.inv.tXX = solve(t(X)%*%X)
model.tHH = t(model.H)%*%model.H
model.Dw_minus_W = model.Dw-model.W

#beta
gibbs.cond.beta = function(last_param){
    cond.mean = model.inv.tXX %*% t(X) %*% (y - model.H %*% last_param$eta)
    cond.var = last_param$sigma2 * model.inv.tXX
    new_beta = rmvnorm(1, mean=cond.mean, sigma=cond.var)
    return(new_beta)
}

gibbs.cond.eta = function(last_param){
    cond.precision = (1/last_param$sigma2) * model.tHH + (1/last_param$tau2) * model.Dw_minus_W
    cond.var = solve(cond.precision)
    cond.mean = (1/last_param$sigma2)*cond.var %*% t(model.H) %*% (y - X%*%last_param$beta)
    new_eta = rmvnorm(1, mean=cond.mean, sigma=cond.var)
    new_eta = new_eta - mean(new_eta) # imposing constraint
    return(new_eta)
}

gibbs.cond.sigma2 = function(last_param){
    cond.shape = 0.001 + model.n/2
    diff_fit = y - X%*%last_param$beta - model.H%*%last_param$eta
    cond.scale = 0.001 + 0.5 * t(diff_fit) %*% diff_fit
    new_sigma2 = rinvgamma(1, shape=cond.shape, scale=cond.scale)
    return(new_sigma2)
}

gibbs.cond.tau2 = function(last_param){
    cond.shape = 0.001 + (model.m - 1)/2
    cond.scale = 0.001 + 0.5 * t(last_param$eta) %*% model.Dw_minus_W %*% last_param$eta
    new_tau2 = rinvgamma(1, shape=cond.shape, scale=cond.scale)
    return(new_tau2)
}

# sampler
gibbs.iter.B = 10000
gibbs.param_mat = matrix(0, gibbs.iter.B, 12+380+1+1) # order: beta, eta, sigma2, tau2

colname_set = rep("",12+380+1+1)
for(i in 1:(12+380+1+1)){
    if(i>0 & i<=12) colname_set[i] = paste("beta",i, sep="")
    if(i>12 & i<=12+380) colname_set[i] = paste("eta",i-12, sep="")
    if(i==12+380+1) colname_set[i] = "sigma2"
    if(i==12+380+1+1) colname_set[i] = "tau2"
}
colnames(gibbs.param_mat)=colname_set

gibbs.initial_cond = list(
    beta = c(rnorm(12,0,1)), #12
    eta = c(rnorm(380,0,10)), #380
    sigma2 = 1, #1
    tau2 = 1 #1
)

gibbs.last_param = gibbs.initial_cond
gibbs.param_mat[1,]=unlist(gibbs.last_param)

gibbs.start_time = Sys.time()
print("gibbs: start sampling")
for(i in 2:gibbs.iter.B){
    if(i%%20==0) {
        elapsed = Sys.time()-gibbs.start_time
        est_time_total = elapsed * gibbs.iter.B / i
        est_time_remain = est_time_total - elapsed
        cat("iteration: ", i,"/",gibbs.iter.B,", remain(estimated,sec):", est_time_remain,"\n")
    }
    gibbs.last_param$beta = c(gibbs.cond.beta(gibbs.last_param))
    gibbs.last_param$eta = c(gibbs.cond.eta(gibbs.last_param))
    gibbs.last_param$sigma2 = c(gibbs.cond.sigma2(gibbs.last_param))
    gibbs.last_param$tau2 = c(gibbs.cond.tau2(gibbs.last_param))
    gibbs.param_mat[i,] = unlist(gibbs.last_param)
}
print("gibbs: end sampling")
gibbs.end_time = Sys.time()
cat("gibbs: execution time: ", gibbs.end_time-gibbs.start_time, "\n")


# save(gibbs.param_mat,file="C:/gitProject/SpaTempoDA/HW3/__Gibbsresult.RData")
# load("C:/gitProject/SpaTempoDA/HW3/__Gibbsresult.RData")


# burn-in cut
burn_in_samples_num = 1000
gibbs.param_mat_aftercut = gibbs.param_mat[burn_in_samples_num:gibbs.iter.B,]



# plot
plotting_func = function(param_idx){
    param_name = colname_set[param_idx]
    idx.sample = gibbs.param_mat_aftercut[,param_idx]

    title = paste("posterior samples of",param_name)
    
    #hist
    hist(idx.sample, breaks=100, main=title, col="gray", lty="blank", xlab=param_name)
    
    prob_interval = quantile(idx.sample, c(0.025, 0.975))
    lines(prob_interval, c(0,0), col="blue", lwd=5)
    text(prob_interval[1],0.1, labels=round(prob_interval[1],3),pos=3)
    text(prob_interval[2],0.1, labels=round(prob_interval[2],3),pos=3)
    
    abline(v=mean(idx.sample), col="blue", lwd=3)
    text(mean(idx.sample),0.1, labels=round(mean(idx.sample),3),pos=3)

    filename = paste("C:/gitProject/SpaTempoDA/HW3/",param_name,"_hist.png",sep="")
    dev.copy(png,filename=filename)
    dev.off()

    #traceplot
    idx.ESS = effectiveSize(idx.sample)
    xlab_text = paste("ESS=",round(idx.ESS,2),sep="")
    plot(1:length(idx.sample), idx.sample, main=title, type="l",
        xlab=xlab_text, ylab=param_name)
    abline(h=mean(idx.sample), col="blue", lwd=3)
    
    filename = paste("C:/gitProject/SpaTempoDA/HW3/",param_name,"_traceplot.png",sep="")
    dev.copy(png,filename=filename)
    dev.off()

    #acf
    acf(idx.sample, main=title)
    filename = paste("C:/gitProject/SpaTempoDA/HW3/",param_name,"_acf.png",sep="")
    dev.copy(png,filename=filename)
    dev.off()

}

for(i in c(1:12,12+380+1,12+380+1+1)){
    plotting_func(i)
}


# summary print
summary.mat = matrix(0,14,3)
summary.namevec = rep("",14)
colnames(summary.mat)=c("mean","0.025q","0.975q")
summary.row_iter = 1
for(i in c(1:12,12+380+1,12+380+1+1)){
    summary.idx.sample = gibbs.param_mat_aftercut[,i]
    summary.namevec[summary.row_iter] = colname_set[i]
    summary.mat[summary.row_iter,] = c(mean(summary.idx.sample), quantile(summary.idx.sample, c(0.025, 0.975)))
    summary.row_iter = summary.row_iter + 1
}

rownames(summary.mat)= summary.namevec
print(round(summary.mat,3))
# > print(round(summary.mat,3))
#           mean  0.025q  0.975q
# beta1  -25.533 -34.024 -17.015
# beta2    0.018   0.014   0.022
# beta3   -1.888  -2.464  -1.310
# beta4   -1.320  -1.716  -0.919
# beta5   -0.648  -0.888  -0.407
# beta6    0.623   0.300   0.944
# beta7    1.358   1.006   1.715
# beta8   -1.242  -1.547  -0.943
# beta9   -1.757  -2.056  -1.456
# beta10  -2.276  -2.649  -1.908
# beta11  -2.228  -2.873  -1.580
# beta12  -2.480  -3.612  -1.374
# sigma2   4.077   3.819   4.353
# tau2     1.202   0.672   1.910


#for etas
gibbs.eta.samples = gibbs.param_mat_aftercut[,(12+1):(12+380)]
gibbs.eta.means = colMeans(gibbs.eta.samples)
gibbs.eta.stds = rep(0,380)
for(i in 1:380){
    gibbs.eta.stds[i] = sd(gibbs.eta.samples[,i])
}


pal <- tim.colors(5)
q5 <- classIntervals(gibbs.eta.means, n = 5, style = "quantile")
col <- findColours(q5, pal)
plot(districts.sp, border = "black", col = col, main="eta means")
plot(parks.sp, add=T, col= "gray")
legend("topright", fill = attr(col, "palette"),
       legend = names(attr(col, "table")),
       bty="n", cex = 0.8, y.intersp = 1.5)
filename = paste("C:/gitProject/SpaTempoDA/HW3/map_posterior_eta_means.png",sep="")
dev.copy(png,filename=filename)
dev.off()


pal <- tim.colors(5)
q5 <- classIntervals(gibbs.eta.stds, n = 5, style = "quantile")
col <- findColours(q5, pal)
plot(districts.sp, border = "black", col = col, main="eta stds")
plot(parks.sp, add=T, col= "gray")
legend("topright", fill = attr(col, "palette"),
       legend = names(attr(col, "table")),
       bty="n", cex = 0.8, y.intersp = 1.5)
filename = paste("C:/gitProject/SpaTempoDA/HW3/map_posterior_eta_stds.png",sep="")
dev.copy(png,filename=filename)
dev.off()
