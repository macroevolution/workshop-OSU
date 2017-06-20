source("diversification_functions1.R")
source("traitDependent_functions.R")

#------------------------------------------------
# Exercise 1: basic maximum likelihood inference
 
library(ape)
skinktree <- read.tree('data/skinks/skinks216.tre')

# dir(...) prints out files in current directory
# getwd(...) tells you what directory you are in

# get number of taxa:
length(skinktree$tip.label)

# get sum of branch lengths:
sum(skinktree$edge.length)
 
# check with lambda = 0.001, 0.01, 0.1, and 1.0

# Implement the likelihood function for the pure-birth model:

# set lambda to 1
lambda <- 0.001
n <- length(skinktree$tip.label)

sum_of_branches <- sum(skinktree$edge.length)

loglik <- (n - 2) * log(lambda) - lambda*sum_of_branches

loglik

#########
ls()  # lists all the things in your R workspace

# Use the loglik_yule function
# from diversification_functions.R
#	This fxn implements the likelihood function we just implemented above


loglike_purebirth(skinktree, lambda = 1)

# look at a lot more values of lambda
lambda_vector <- seq(0.001, 2.0, length.out=1000)
length(lambda_vector)
plot(lambda_vector)

lik_vector <- loglike_purebirth(lambda_vector, skinktree)

# plot lik_vector by the lambda_vector
plot(lik_vector ~ lambda_vector)

# Here we can pull out the value of lambda 
#  for which the value of the log-likelihood was highest
#	e.g., the ML estimate
lambda_vector[which(lik_vector == max(lik_vector)) ]

# OR we can directly obtain the same quantity like this
#	using the analytical solution 
ml_est <- (n - 2) / sum_of_branches
ml_est

##############

# Plot likelihood curve:
plot(x = lambda_vector, y = lik_vector)
lines(x = c(ml_est, ml_est), y= c(-500, -4000), lwd=4, col='red')

 

#------------------------------------------------
# Exercise 2:
# Tree simulation  
# with high and low relative extinction

rm(list = ls())
library(phytools)
source('diversification_functions1.R')
library(diversitree)

ls()

# check arguments names to simulateTreee
args(simulateTree)

# pars: vector of c(lambda, mu)
# max.taxa = number of species we want in tree


# set parameters for simulation
lambda <- 1
mu <- 0

# First, generate tree with no extinction and 100 tips:
tree_pb <- simulateTree(pars = c(lambda, mu), max.taxa = 100)

plot.phylo(tree_pb, show.tip.label=F)

# Now, generate tree with high relative extinction and 100 tips
tree_highE <- simulateTree(pars=c(1, 0.99), max.taxa=100)

## Here we will plot the two trees side by side:
##	one with 0 relative extinction, 
#	one with very high relative extinction


plot.new()
par(mfrow=c(1,2))
par(mar=c(0,0,0,0))
plot.phylo(tree_pb, show.tip.label=F)
plot.phylo(tree_highE, show.tip.label=F)

# Visualize lineage-through-time plots for the 2 trees:

# plot lineage through time plots:
ltt(tree_pb)

ltt(tree_highE)

plot.new()
Le <- ltt(tree_pb, plot = F)
He <- ltt(tree_highE, plot=F)

# Now plotting on same plot with isometric (1:1) line for comparison
Le$times <- Le$times / max(Le$times)
He$times <- He$times / max(He$times)


plot(x=c(0,1), y = c(0, log(100)), lwd=3, col="black", type = "l")
lines(Le$times, log(Le$ltt), col="blue", lwd=0.5) 
points(Le$times, log(Le$ltt), pch=21, bg="blue", cex=1.3)
lines(He$times, log(He$ltt), col="red", lwd=0.5) 
points(He$times, log(He$ltt), pch=21, bg = "red", cex=1.3)

 
#-------------------------------------------------------
# Birth-death models
##### Fit constant-rate birth-death model to simulated tree

# fitCRBD  = fits a constant-rate diversification model

fitCRBD(tree_pb)

fitCRBD(tree_highE)

 

#---------------------------------------------------------------------
# Exercise 3: simulation plus ML parameter estimation
#
# Here we will fit the constant-rate birth-death
#	model to 1000 simulated trees,
#	and where each simulation gets a unique speciation-extinction
#	parameterization

REPS <- 1000
x <- numeric(REPS)
res <- data.frame(true_lambda=x, true_mu=x, est_lambda=x, est_mu=x)


for (i in 1:REPS){
	cat(i, '\n')

#	pick a lambda:	
	lambda <- runif(1, 0, 2)

# pick a relative extinction rate:
	rel_ex <- runif(1  , 0, 0.95)
 
 # calculate mu
	mu <- lambda * rel_ex 
	# simulate tree:	
	tree_sim <- simulateTree(pars=c(lambda, mu), max.taxa=100)
 
	fit <- fitCRBD(tree_sim)
	
	res$true_lambda[i] <- lambda
	res$true_mu[i] <- mu

	res$est_lambda[i] <- fit["lambda"]
	res$est_mu[i] <- fit["mu"]

}


#----------------------------------------------
# do 1000 replicates
# test correlation between true_lambda and est_lambda
# test correlation between true_mu and est_mu

# plot estimated lambda by true lambda
plot(x = res$true_lambda, y = res$est_lambda )

# Plot both:
plot.new()
par(mfrow=c(1,2))
plot(x = res$true_lambda, y = res$est_lambda )
plot(x = res$true_mu, y = res$est_mu , col='red')


cor.test(res$true_lambda, res$est_lambda )
cor.test(res$true_mu, res$est_mu )


 














