# R code created by Choi Seokjun
# UTF-8 CRLF
# windows 10 64-bit.
# R 3.6.3


library(spatstat)

#
data(bei)
plot(bei)

# bw.scott(bei) #76.80943 41.01130
plot(bei)
filename = paste("C:/gitProject/SpaTempoDA/HW4/prob3_bei.png",sep="")
dev.copy(png,filename=filename)
dev.off()

plot(density(bei, sigma = bw.scott(bei)))
points(bei, pch = 20, cex=0.25)
filename = paste("C:/gitProject/SpaTempoDA/HW4/prob3_intensity_nocovariates.png",sep="")
dev.copy(png,filename=filename)
dev.off()


head(bei.extra)
plot(bei.extra$elev)
filename = paste("C:/gitProject/SpaTempoDA/HW4/prob3_elev.png",sep="")
dev.copy(png,filename=filename)
dev.off()
plot(bei.extra$grad)
filename = paste("C:/gitProject/SpaTempoDA/HW4/prob3_grad.png",sep="")
dev.copy(png,filename=filename)
dev.off()



# model 1
grad = bei.extra$grad
fit.model1 = ppm(bei, ~slope, covariates=list(slope=grad))
# > ppm(bei, ~slope, covariates=list(slope=grad))
# Nonstationary Poisson process

# Log intensity:  ~slope

# Fitted trend coefficients:
# (Intercept)       slope
#   -5.391053    5.026710 

#              Estimate       S.E.   CI95.lo   CI95.hi Ztest      Zval
# (Intercept) -5.391053 0.03001787 -5.449887 -5.332219   *** -179.5948
# slope        5.026710 0.24534296  4.545847  5.507573   ***   20.4885
# => λ(u) = exp(−5.391053 + 5.026710 Z(u)). 

plot(fit.model1, how = "image", se = FALSE)
filename = paste("C:/gitProject/SpaTempoDA/HW4/prob3_model1_fit.png",sep="")
dev.copy(png,filename=filename)
dev.off()

plot(predict(fit.model1, type="cif", ngrid=2^10))
points(bei, pch = 20, cex=0.25)
filename = paste("C:/gitProject/SpaTempoDA/HW4/prob3_model1_predict.png",sep="")
dev.copy(png,filename=filename)
dev.off()



# model 2
fit.model2 = ppm(bei, ~offset(log(slope)), covariates = list(slope = grad))
#              Estimate       S.E.   CI95.lo   CI95.hi Ztest      Zval
# (Intercept) -2.427165 0.01665742 -2.459813 -2.394517   *** -145.7108
# => λ(u) = exp(-2.427165)Z(u) = 0.08828677 Z(u)

plot(fit.model2, how = "image", se = FALSE)
filename = paste("C:/gitProject/SpaTempoDA/HW4/prob3_model2_fit.png",sep="")
dev.copy(png,filename=filename)
dev.off()

plot(predict(fit.model2, type="cif", ngrid=256))
points(bei, pch = 20, cex=0.25)
filename = paste("C:/gitProject/SpaTempoDA/HW4/prob3_model2_predict.png",sep="")
dev.copy(png,filename=filename)
dev.off()
