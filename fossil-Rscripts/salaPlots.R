rm(list=ls())

library(BAMMtools)

extant <- read.tree("extant.tre")
fossil <- read.tree("salamandridae.tre")

ex_edata <- getEventData(extant, "extant_event_data.txt", burnin=0.1, nsamples=200)
f_edata <- getEventData(fossil, "event_data.txt", burnin=0.1, nsamples=200)

ex_Rates <- getCladeRates(ex_edata)
f_Rates <- getCladeRates(f_edata)

par(mfrow=c(2, 2), bty="n", las=1)
hist(ex_Rates$lambda, col='gray70', border='gray70', main="",  xlab="", ylab="", xlim=c(0, 0.2), ylim=c(0, 80))
hist(f_Rates$lambda, add=T, border=rgb(1, 0, 0, 0.5), col=rgb(1, 0, 0, 0.5))

hist(ex_Rates$mu, col='gray70', border='gray70', main="",  xlab="", ylab="", xlim=c(0, 0.2), ylim=c(0, 80))
hist(f_Rates$mu, add=T, border=rgb(1, 0, 0, 0.5), col=rgb(1, 0, 0, 0.5))

hist(ex_Rates$lambda - ex_Rates$mu, col='gray70', border='gray70', main="",  xlab="", ylab="", xlim=c(0, 0.11), ylim=c(0, 60))
hist(f_Rates$lambda - f_Rates$mu, add=T, border=rgb(1, 0, 0, 0.5), col=rgb(1, 0, 0, 0.5))

hist(ex_Rates$mu/ex_Rates$lambda, col='gray70', border='gray70', main="",  xlab="", ylab="", xlim=c(0, 1), ylim=c(0, 60))
hist(f_Rates$mu/f_Rates$lambda, add=T, border=rgb(1, 0, 0, 0.5), col=rgb(1, 0, 0, 0.5))

