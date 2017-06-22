




#----------------------------------------------
# Exercise 4: Gamma statistic
 
warbs <- read.tree("data/warblers/warbs.tre")

gammaStat(warbs) 

skinks <- read.tree("data/skinks/skinks216.tre")
gammaStat(skinks)



#-----------------------------------------------
# Exercise 4b: Simulate trees with time-dependent 
#               speciation and extinction rates
#

library(TESS)

lambda_fx <- function(x) return(1 * exp(x * -0.2))


# look at the rates that this function would generate:
tvec <- seq(0, 10, length.out=100)
rates <- lambda_fx(tvec)
plot(rates ~ tvec)


# We can simulate trees conditioned on a particular age:

tree <- tess.sim.age(n=1, age = 10, lambda = lambda_fx, mu = 0)

plot(tree[[1]])

gammaStat(tree[[1]])

# We can simulate trees dependent on a certain number of tips:
#    argument max is the maximum possible age of the tree we allow:

tree2 <- tess.sim.taxa(n=1, nTaxa = 100, lambda = lambda_fx, mu = 0, max=25)
gammaStat(tree2[[1]])











