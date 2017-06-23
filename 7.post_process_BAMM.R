



#-------------------------------------------
# BAMM run setup:
# 

library(BAMMtools)
whales <- read.tree("https://macroevolution.github.io/workshop-OSU/data/whales/whaletree.tre")


# Now we use setBAMMpriors to match the scale of priors
#   to the scale of the tree:

setBAMMpriors(whales)


# Now copy and paste the prior block in the relevant location in your control file



##############################################
#---------------------------------------------
# Postprocess 3 : Testing convergence

library(coda)

# If you could not get BAMM to run correctly, 
#  you can simply load some of our example runs
#  distributed with BAMMtools

# data(mcmc.whales)


mcmc.whales <- read.csv("mcmc_out.txt")


# get rid of burnin
dim(mcmc.whales)
head(mcmc.whales)

xx <- mcmc.whales[201:2000, ]

# plot a simple MCMC trace of log-likelihoods:
plot(xx$logLik ~ xx$generation)

effectiveSize(xx$logLik)
effectiveSize(xx$N_shifts)

#---------------------------------------------



##############################################
#---------------------------------------------
# Postprocess 4: how much rate variation?

# Simple summary of posterior
table(xx$N_shifts ) / nrow(xx)

#this is the posterior distribution of shifts!

# How does it compare to the prior?
plotPrior(xx)

# What about a more formal model comparison that takes
#   the prior into account?
# We will compute a matrix of pairwise Bayes factors:

bfmat <- computeBayesFactors(xx, expectedNumberOfShifts = 1)

# comparisons of each model against the null (zero shift) model:
bfmat[, "0"]

 
#---------------------------------------------

##############################################
#---------------------------------------------
# Postprocess 5.3: Reading event data / basic plotting

# If you cannot get BAMM to run:
# data(events.whales)
# ed <- getEventData(whales, events.whales, burnin = 0.1, nsamples = 200)


ed <- getEventData(whales, "event_data.txt", nsamples = 200)

summary(ed)
plot.bammdata(ed)

# lots of options to explore
#
#  tau
#  pal

z <- plot.bammdata(ed, tau = 0.002, lwd=2)
addBAMMlegend(z)


# View an individual shift configuration

x <- subsetEventData(ed, index = 10)
plot.bammdata(x, tau = 0.002, lwd=3)
addBAMMshifts(x, cex=2)

x <- subsetEventData(ed, index = 20)
plot.bammdata(x, tau = 0.002, lwd=3)
addBAMMshifts(x, cex=2)


# Loop over set of samples from the posterior and plot each:

pvec <- 1:9 * 10

quartz.options(height=10, width=10)
par(mar=c(1,1,1,1))

plot.new()
par(mfrow=c(3,3))

# get general color map: 
zcol <- plot.bammdata(ed, tau = 0.002, show=F)

for (ii in pvec ){
	x <- subsetEventData(ed, index = ii)
 	plot.bammdata(x, tau = 0.002, method = "polar", lwd=1.3, colorbreaks = zcol$colorbreaks)
	addBAMMshifts(x, cex=2, par.reset=F, method= "polar")
	mtext(text = paste("sample ", ii, sep=""), side= 1)
	
}




#---------------------------------------------

##############################################
#---------------------------------------------
# Postprocess 6: clade rates and tip rates

# Viewing some node numbers on the tree:
plot.phylo(whales, cex = 0.7)
nodelabels()

tcols <- rep("black", length(whales$tip.label))
dclade <- extract.clade(whales, node = 140)$tip.label
tcols[whales$tip.label %in% dclade] <- "red"

# Highlighting the descendants of the dolphin clade
plot.phylo(whales, tip.color = tcols, cex = 0.7) 
 
# getCladeRates = mean rates across clade for each sample
#                  in posterior

whalerates <- getCladeRates(ed)
names(whalerates)
whalerates$lambda


plot.new()
par(mfrow=c(2,1))
hist(whalerates$lambda, breaks=50, xlim = c(0, 0.2), col = "red")
hist(whalerates$mu, breaks=50, xlim=c(0, 0.2), col="blue")

# Now we are going to get mean clade rates for 3 groups:
#   1. the dolphins only
#   2. the "not-dolphins" only
#   3. everything


# Just focus on node 140 (common ancestor of dolphins)
#      for now:

dolphins1 <- getCladeRates(ed, node = 140)
mean(dolphins1$lambda)

# not dolphins, using nodetype = "exclude"
not_dolphins <- getCladeRates(ed, node = 140, nodetype = "exclude")
mean(not_dolphins$lambda)

# setup 3 panel plot:
plot.new()
par(mfrow=c(3,1))
hist(dolphins1$lambda, breaks=50, col="red", xlim=c(0,0.35))
hist(not_dolphins$lambda, breaks=50, col = "blue", xlim=c(0,0.35))
hist(whalerates$lambda, breaks = 50, col = "lightgreen", xlim=c(0,0.35))

# TIP RATES (think of equal-splits...) 
 
tiprates <- getTipRates(ed)
tip_lambda <- tiprates$lambda.avg

# How do these compare to the equal splits rates?
source("supporting/traitDependent_functions.R")

es <- getEqualSplitsSpeciation(whales)
es <- es[names(tip_lambda)]
cor.test(es, tip_lambda )

# What is the difference between tip rates and clade rates,
#   as estimated by BAMM??? 



#---------------------------------------------

#---------------------------------------------

##############################################
#---------------------------------------------
# Postprocess 7: Visualizing rates through time

# For the whales dataset:
 
plotRateThroughTime(ed, avgCol = "blue")

# this is the overall rate-through-time trajectory

# But we already found that different clades had different dynamics
# so what if we split them into "dolphins" and "not-dolphins"?

ivec <- seq(0.1, 0.9, length.out=100)

plot.new()
par(mfrow=c(3,1))
plotRateThroughTime(ed, node = 140, intervals = ivec, intervalCol = "red", start.time = 36, ylim=c(0, 0.35))
mtext(side=3, text = "Dolphins only", line = -2)

plotRateThroughTime(ed, node = 140, intervals = ivec, avgCol = "blue", start.time = 36, ylim=c(0, 0.35), nodetype = "exclude")
mtext(side=3, text = "Not a dolphin", line = -2)
 
# and what if we used everything together?
plotRateThroughTime(ed, intervals = ivec, intervalCol = "darkgreen", avgCol = "darkgreen", start.time = 36, ylim=c(0, 0.35))
mtext(side=3, text = "all whales, averaged together", line = -2)
 


# Putting together on the same plot:
par(new=T)
plot.new() 
plotRateThroughTime(ed, node = 140, intervals = ivec, intervalCol = "red", start.time = 36, ylim=c(0, 0.35))
 
plotRateThroughTime(ed, node = 140, intervals = ivec, avgCol = "blue", start.time = 36, ylim=c(0, 0.35), nodetype = "exclude", add=T)
 
 
#---------------------------------------------
# Postprocess 8
# 
# Distinct shift configurations
# This plots the 95% credible set 
#      of topologically-distinct shift configurations...

css <- credibleShiftSet(ed, expectedNumberOfShifts = 1)
plot(css, lwd=1.5)

# what is the overall "best" shift configuration?

best <- getBestShiftConfiguration(ed, expectedNumberOfShifts = 1)
z <- plot(best, lwd=2)
addBAMMshifts(best, cex=3)
addBAMMlegend(z)

#---------------------------------------------
# Postprocess 9
#  marginal shift probabilities:

mst <- marginalShiftProbsTree(ed)

plot(mst, cex = 0.7)
add.scale.bar(length = 0.5, lcol="red")

# side-by-side with tree:
par(new=TRUE)
plot.new()
par(mfrow=c(1,2))

plot(whales, cex = 0.7)
plot(mst, cex = 0.7)
axis(1)
mtext(side = 1, text = "probability", line=3)

#---------------------------------------------
# Postprocess 10
# cohorts

cmat <- getCohortMatrix(ed)
cohorts(cmat, ed)

# or including the phylorate plots:
cmat <- getCohortMatrix(ed)
cohorts(cmat, ed, use.plot.bammdata = TRUE, lwd=1.5)

