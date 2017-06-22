



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

# data(mcmc.whales)
# mcmc.whales <- read.csv("mcmc_out.txt")

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
plot.phylo(whales)
nodelabels()

# Just focus on node 140 for now:

dolphins1 <- getCladeRates(ed, node = 140)
mean(dolphins1$lambda)

# not dolphins
not_dolphins <- getCladeRates(ed, node = 140, nodetype = "exclude")
mean(not_dolphins$lambda)

# get tip rates:
tiprates <- getTipRates(ed)
tip_lambda <- tiprates$lambda.avg

# How do these compare to the equal splits rates?

es <- getEqualSplitsSpeciation(whales)
es <- es[names(tip_lambda)]
cor.test(es, tip_lambda )





#---------------------------------------------

#---------------------------------------------

##############################################
#---------------------------------------------
# Postprocess 7: Visualizing rates through time

# For the whales dataset:






 










