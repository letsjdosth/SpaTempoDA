# R code created by Choi Seokjun
# UTF-8 CRLF
# windows 10 64-bit.
# R 3.6.3

library(spatstat)
dental = read.table("c:/gitProject/SpaTempoDA/Hw4/dental.reduced.dat", stringsAsFactors=FALSE)
head(dental)
dim(dental)


dental.unaffected = dental[(dental[2]==0),]
dental.affected = dental[(dental[2]==1),]
dim(dental.unaffected)
dim(dental.affected)

head(dental.affected)
range(dental.unaffected[3]) #4492 10324
range(dental.affected[3]) #5035 10373

range(dental.unaffected[4]) #3469 10431
range(dental.affected[4]) #3085 9541


#1. rectangular window
# ?ppp
# ?owin
rec_window = owin(xrange=c(4492, 10373), yrange=c(3085, 10431))
ppp.rec.unaffected = ppp(unlist(dental.unaffected[3]), unlist(dental.unaffected[4]), window=rec_window)
ppp.rec.affected = ppp(unlist(dental.affected[3]), unlist(dental.affected[4]), window=rec_window)
plot(ppp.rec.unaffected, pch = 19, main="rectangular window")
plot(ppp.rec.affected, pch = 19, main = "affected", col = "red", add=TRUE)
filename = paste("C:/gitProject/SpaTempoDA/HW4/prob1_rectangular_window_scatterplot.png",sep="")
dev.copy(png,filename=filename)
dev.off()



#2.
# ?envelope
# The summary statistic fun is applied to each of these simulated patterns. 
# Typically fun is one of the functions Kest, Gest, Fest, Jest, pcf, Kcross, Kdot, Gcross, Gdot, Jcross, Jdot, Kmulti, Gmulti, Jmulti or Kinhom. 
# It may also be a character string containing the name of one of these functions.

evlp.rec.unaffected.F = envelope(ppp.rec.unaffected, Kest, nsim = 999, nrank = 10, global = TRUE)
plot(evlp.rec.unaffected.F)
filename = paste("C:/gitProject/SpaTempoDA/HW4/prob1_rectangular_window_unaffected_evlp_F.png",sep="")
dev.copy(png,filename=filename)
dev.off()

evlp.rec.affected.F = envelope(ppp.rec.affected, Kest, nsim = 999, nrank = 10, global = TRUE)
plot(evlp.rec.affected.F)
filename = paste("C:/gitProject/SpaTempoDA/HW4/prob1_rectangular_window_affected_evlp_F.png",sep="")
dev.copy(png,filename=filename)
dev.off()

evlp.rec.unaffected.G = envelope(ppp.rec.unaffected, Gest, nsim = 999, nrank = 10, global = TRUE)
plot(evlp.rec.unaffected.G)
filename = paste("C:/gitProject/SpaTempoDA/HW4/prob1_rectangular_window_unaffected_evlp_G.png",sep="")
dev.copy(png,filename=filename)
dev.off()
evlp.rec.affected.G = envelope(ppp.rec.affected, Gest, nsim = 999, nrank = 10, global = TRUE)
plot(evlp.rec.affected.G)
filename = paste("C:/gitProject/SpaTempoDA/HW4/prob1_rectangular_window_affected_evlp_G.png",sep="")
dev.copy(png,filename=filename)
dev.off()



#3. window of rough polygon outline to surround the points.
#locator
locator_set_mod = FALSE #for first run, change it TRUE and set polygon. and modify 'else' block
if(locator_set_mod){
    plot(ppp.unaffected, pch = 19, main = "unaffected")
    plot(ppp.affected, pch = 19, main = "affected", col = "red", add=TRUE)
    bdry <- locator()
}else{
    #given : print(poly_window$bdry)
    bdry = list(x=c(10456.320, 9415.281, 4354.227, 4418.291, 5859.730, 6388.258, 8021.889, 9399.265, 10536.400),
        y=c(8333.973, 10560.196, 7677.317, 5595.238, 5579.222, 6828.470, 4089.735, 3064.712, 3016.664))
}

poly_window = owin(poly=bdry)
# print(poly_window$bdry)
ppp.poly.unaffected = ppp(unlist(dental.unaffected[3]), unlist(dental.unaffected[4]), window=poly_window)
ppp.poly.affected = ppp(unlist(dental.affected[3]), unlist(dental.affected[4]), window=poly_window)
plot(ppp.poly.unaffected, pch = 19, main = "polygon-window")
plot(ppp.poly.affected, pch = 19, col = "red", add=TRUE)
filename = paste("C:/gitProject/SpaTempoDA/HW4/prob1_poly_window_scatterplot.png",sep="")
dev.copy(png,filename=filename)
dev.off()


#4.
evlp.poly.unaffected.F = envelope(ppp.poly.unaffected, Kest, nsim = 999, nrank = 10, global = TRUE)
plot(evlp.poly.unaffected.F)
filename = paste("C:/gitProject/SpaTempoDA/HW4/prob1_poly_window_unaffected_evlp_F.png",sep="")
dev.copy(png,filename=filename)
dev.off()

evlp.poly.affected.F = envelope(ppp.poly.affected, Kest, nsim = 999, nrank = 10, global = TRUE)
plot(evlp.poly.affected.F)
filename = paste("C:/gitProject/SpaTempoDA/HW4/prob1_poly_window_affected_evlp_F.png",sep="")
dev.copy(png,filename=filename)
dev.off()

evlp.poly.unaffected.G = envelope(ppp.poly.unaffected, Gest, nsim = 999, nrank = 10, global = TRUE)
plot(evlp.poly.unaffected.G)
filename = paste("C:/gitProject/SpaTempoDA/HW4/prob1_poly_window_unaffected_evlp_G.png",sep="")
dev.copy(png,filename=filename)
dev.off()

evlp.poly.affected.G = envelope(ppp.poly.affected, Gest, nsim = 999, nrank = 10, global = TRUE)
plot(evlp.poly.affected.G)
filename = paste("C:/gitProject/SpaTempoDA/HW4/prob1_poly_window_affected_evlp_G.png",sep="")
dev.copy(png,filename=filename)
dev.off()
