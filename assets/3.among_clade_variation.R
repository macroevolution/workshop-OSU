
source("supporting/diversification_functions1.R")
source("supporting/traitDependent_functions.R")
source("supporting/simulate_shift_trees.R")

#----------------------------------------------
#   Exercise 5: Chance!
 
# 	How much variation in species richness can you get from 
#    the same diversification process?

lambda <- 0.2
mu <- 0
max.t <- 25


REPS <- 1000
taxon_count <- numeric(REPS)

pars <- c(lambda = lambda, mu = mu)

simulateTree(pars = pars , max.t = max.t )

for (i in 1:REPS){
	cat(i, '\n')
	tree <- simulateTree(c(0.2, 0), max.t=25)
	taxon_count[i] <- length(tree$tip.label)

}

#hist: plots a histogram
hist(taxon_count, breaks=100)
mean(taxon_count)
 
#  
#----------------------------------------------
#   Exercise 5b: Colless index ** Bonus exercise! **
#      Implement function to compute colless index
#      Apply it to a real dataset
#      Assess significance by simulation

# Here we define the colless imbalance statistic
colless <- function(x){
	
	#N <- length(x$tip.label)
	nn <- balance(x)
	cstat <- sum(abs(nn[,1] - nn[,2]))
	#return(2 * cstat / ((N - 1) * (N - 2)))
	return(cstat)
}

skinks <- read.tree("data/skinks/skinks216.tre")

# this is a tree with 216 tips and is >95% complete at
#   the species level

# Now to assess significance: 
#     Simulate trees under constant-rate model
#     compute colless for each
#     store value
#  This procedure gives us a null distribution that 
#    we can compare to the observed

source("supporting/diversification_functions1.R")

simulateTree(c(1,0), max.taxa = 216)

null_colless <- rep(NA, 1000)

for (ii in 1:1000){
	cat(ii, "\n")
	tree <- simulateTree(c(1,0), max.taxa = 216)
	null_colless[ii] <- colless(tree)
}

# Plot the null distribution
hist(null_colless, breaks=50)

# Visualize: how imbalanced is the skink tree
#    relative to the simulations under constant-rate model?

obs <- colless(skinks)
lines(x=c(obs, obs), y = c(0, 100), lwd=3, col="red")

# and the pvalue, two-tailed:
2* sum(null_colless > obs) / (1 + length(null_colless))



#----------------------------------------------
#   Exercise 6: Tip-specific rates!

rm(list = ls())
source("traitDependent_functions.R")

skinks <- read.tree("data/skinks/skinks216.tre")

t_rates <- getEqualSplitsSpeciation(skinks)

# make a function that interpolates a set of colors
fx <- colorRampPalette(c("blue",  "red"))
colset <- fx(100)

#plot(x = 1:100, y = 1:100, pch = 19, col = colset)

 
colvec <- colorFunction(t_rates, min = quantile(t_rates, 0.25), maxx = quantile(t_rates, 0.75), colset)
 
plot(skinks, type = "fan", show.tip.label=F)
tiplabels(pch=21, bg=colvec, cex=1.5)



#----------------------------------------------
#   Exercise 7: simulate shift trees and compute tip-specific rates

source("simulate_shift_trees.R")
library(BAMMtools)


#  these are exponential distributions with means of 0.15 and 0.05 
lamfx <- function() return(rexp(1, 1/0.15))
mufx <- function() return(rexp(1, 1/0.05))

# rate at which events occur along the phylogeny
trate <- 0.006

tt <- simulateShiftTree(35.5, trate, lamfx, mufx, seed=8)

ed <- getEventData(phy = tt$phy, eventdata = tt$events)
 
z <- plot.bammdata(ed, lwd=2, breaksmethod = "linear")
addBAMMlegend(z)



true_lambda <- getTipRates(ed)$lambda.avg

rates <- getEqualSplitsSpeciation(tt$phy)

plot( rates ~  true_lambda, xlim=c(0,1), ylim=c(0,1))
abline(0,1)
cor.test(rates, true_lambda, method = "spear")



#----------------------------------------------
#   Exercise 8: simulate 
#       a batch of shift trees and compute tip-specific rates
#       and check correlation with true rates!

REPS <- 50

cormat <- matrix(NA, nrow = REPS, ncol = 3)

seedvec <- 1:REPS + 100


for (i in 1:REPS){
	cat(i, "\n")
	tt <- simulateShiftTree(35.5, trate, lamfx, mufx, seed=seedvec[i])
	
	if (!is.na(tt$phy)[1]){
		
		ed <- getEventData(phy = tt$phy, eventdata = tt$events)
		true_lambda <- getTipRates(ed)$lambda.avg
		rates <- getEqualSplitsSpeciation(as.phylo(ed))
	 	cc    <- cor.test(true_lambda, rates )
	 	cormat[i, 1] <- cc$estimate	
		tx <- sort(table(ed$tipStates), decreasing=T)
		cormat[i,2] <- length(ed$tip.label)
		cormat[i,3] <- tx[2]
		
	}
	
	
	
}
















