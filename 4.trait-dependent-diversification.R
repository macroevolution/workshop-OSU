

#------------------------------------------------
# Testing for trait-dependent diversification using FiSSE (and BiSSE; see below)

# Example script illustrating use of FISSE test

library(ape)
library(phangorn)
library(diversitree)

source("supporting/traitDependent_functions.R")

tree          <- read.tree("example_tree.tre")
xx            <- read.csv("example_trait.csv", header=F)
traits        <- xx[,2]
names(traits) <- xx[,1]
 

traits <- traits[tree$tip.label]

#---- Plot tree w traits before analysis

colvec <- rep("white", length(traits))
colvec[traits == 1] <- "black"

quartz.options(height=12, width=12)
plot.phylo(tree, type = "fan", show.tip.label=F)
tiplabels(pch=21, bg=colvec, cex=0.8)

#------- Arguments to FISSE.binary function	
#
#    phy               = class phylo phylogenetic tree
#    states            = vector of state data (must have names attribute)
#    reps              = Number of simulations to perform when generating null distribution
#    tol               = tolerance value to accept simulation as valid
#                           default value of 0.1 means that null simulations will be accepted 
#                           if parsimony-reconstructed changes are +/- 10% of the observed value
#    qratetype         = How to estimate rates for the mk1 simulations of the null distribution. 
#                            default value of mk1 fits a symmetric 1-rate Mk model to the data
#                            option "parsimony" simply divides the 
#                             number of changes by the summed edge lengths
#                     
 
 
res <- FISSE.binary(tree, traits)


# components of returned object:
#
#    lambda0           = inverse splits rate estimate for state 0
#    lambda1           = inverse splits rate estimate for state 1
#    pval              = proportion of simulations where observed 
#                             test statistic (lambda1 - lambda0)
#                           greater than the simulated value.
#   null_mean_diff     =   average (lambda1 - lambda0) for null distribution           
#   null_sd            =  std deviation of the null distribution
#   nchanges_parsimony = number of parsimony-reconstructed changes in trait values
#   qpars              = transition rate under symmetric (Mk1) model of trait evolution
#	                         (used to simulated null distribution)

# two-tailed pvalue is obtained as
pval_2tailed   <- min(res$pval, 1-res$pval)*2

# this example should give significant evidence for SDD 
#    (this is correct: for this dataset, 
#        state 0 has lambda = 0.05 and state 1 has lambda = 0.1, with mu0 = mu1 = 0.01)
#
#
 
# ------------------------------------------
#     Simulation: generate tree under true trait-dependent diversification

# Here we will set up parameters for simulation under a BiSSE process with
#   a derived state that represents a 3x increase in the rate of speciation.
pars <- c(lambda0 = 0.1, lambda1 = 0.3, mu0=0, mu1=0, q01 = 0.05, q10 = 0.05)

# Simulating the tree, with 100 extant tips:
xtree <- tree.bisse(pars, max.taxa = 100) 
table(xtree$tip.state)

states <- xtree$tip.state
colvec <- rep("white", length(states))
colvec[states == 1] <- "black"

quartz.options(height=12, width=12)
plot.phylo(xtree, type = "fan", show.tip.label=F)
tiplabels(pch=21, bg=colvec, cex=0.8) 
 
# Try FiSSE on this dataset:

res <- FISSE.binary(xtree, states)


 
# -------- EMPIRICAL EXAMPLE ---------------


# ------------------------------------------
# Trait-dependent analysis: speciation in accipiters (hawks) as a
#  function of color polymorphism

rm(list = ls())
library(diversitree)
source("supporting/diversification_functions1.R")
source("supporting/traitDependent_functions.R")


tree <- read.nexus('data/accipitrids/accipiters.nex')
x <- read.csv('data/accipitrids/accip.csv', stringsAsFactors=F)
states <- x$CP
names(states) <- x$species

colvec <- rep('white', length(states))
colvec[states == 1] <- 'black'
plot.phylo(tree, type = 'fan', show.tip.label=F)
tiplabels(pch=21, bg = colvec)

 
#-------------------------------------
#  Formal BiSSE analysis

#     make 5 parameter TRAIT DEPENDENT speciation model:

sampling <- c(1, 1)

likefunc <- make.bisse(tree, states=states, sampling.f = sampling)
likefunc

like_traitd <- constrain(likefunc, mu0 ~ mu1)
argnames(likefunc)
argnames(like_traitd)

# Now make a 4 parameter model where speciation rates are equal
#   by state:

like_con <- constrain(likefunc, mu0~mu1, lambda0~lambda1)
 

res_td <- fitDiversitree(like_traitd, nopt=3)
res_con <- fitDiversitree(like_con, nopt = 3)

lr <- 2 * (res_td$loglik - res_con$loglik)

1- pchisq(lr, df=1)

# Look at parameters from model
res_td$pars
# which suggests much faster speciation rates for species in state 1 (polymorphic) relative 
#   to non-polymorphic



#-------------------------------------
#  Repeat this analysis with FiSSE

res <- FISSE.binary(tree, states)


res

# How do we interpret this?



#-------------------------------------
#  Revised BiSSE analysis, following Beaulieu and O'Meara (2016)
#    Includes "character-independent" CID-2 model
#
#

library(hisse)

# Revisiting the accipiter (hawk) example:

tree <- read.nexus('data/accipitrids/accipiters.nex')
x <- read.csv('data/accipitrids/accip.csv', stringsAsFactors=F)
 
 
# We will fit BiSSE and the CID2 models using a wrapping function
#  to make things a little easier 
 
res_bisse <- fitBiSSE(tree, x)

res_cid2 <- fitCID2(tree, x)


# Look at AICs:

res_bisse

res_cid2
 









 

