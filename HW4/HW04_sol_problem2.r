# R code created by Choi Seokjun
# UTF-8 CRLF
# windows 10 64-bit.
# R 3.6.3

library(spatstat)

set.seed(2019311252)

hom_pois_process_generator = function(intensity){
    point_num = rpois(1, intensity)
    x_loc = runif(point_num,0,1)
    y_loc = runif(point_num,0,1)
    rec_window = owin(xrange=c(0, 1), yrange=c(0, 1))
    ppp.hom_pois_proc = ppp(x_loc, y_loc, window=rec_window)
    return(ppp.hom_pois_proc)
}

hom_pois_process_generator2 = function(intensity){
    rec_window = owin(xrange=c(0, 1), yrange=c(0, 1))
    return(rpoispp(intensity, lmax = 300, win = rec_window))
}

sim1 = hom_pois_process_generator(100)
sim2 = hom_pois_process_generator(100)
sim3 = hom_pois_process_generator2(100)
sim4 = hom_pois_process_generator2(100)

plot_function = function(ppp.sim, main=""){
    plot(ppp.sim, main=main, pch = 19)
    filename = paste("C:/gitProject/SpaTempoDA/HW4/prob2_",main,"_scatterplot.png",sep="")
    dev.copy(png,filename=filename)
    dev.off()

    plot(density(sim1, sigma = 10), main=main)
    points(sim1, pch = 20, cex=0.25)
    filename = paste("C:/gitProject/SpaTempoDA/HW4/prob2_",main,"_ker.png",sep="")
    dev.copy(png,filename=filename)
    dev.off()
}
plot_function(sim1, "sim1")
plot_function(sim2, "sim2")
plot_function(sim3, "sim3")
plot_function(sim4, "sim4")