library(BAMMtools)
library(paleobioDB)

### Read in tree to get names

Tree <- read.tree("horse.tre")

### Now we're going to use BAMMtools to automatically find which horses are extant

# This function adds a "begin" and "end" element to each edge in the phylogeny
# we can use this to assess when a given lineage went extinct
BTs <- BAMMtools:::NU.branching.times(Tree, return.type="begin.end")

# the longest time since the root for an edge is the "modern day"
endT <- max(BTs$end)

# we can't simply ask which edges end at the modern day, as there will be
# some small amount of numerical error.
# So we instead ask which end *close enough* to the modern that we can count
# them as extant.
tol <- 0.0001
extant_edges <- BTs$edge[which(BTs$end >= endT - tol), 2]

### The names are misspecified relative to what the paleobiology database expects
### so we're going to reformat them

extant_tips <- BTs$tip.label[extant_edges]

spaceNames <- gsub("_", " ", Tree$tip.label)

Names <- paste(spaceNames, collapse=", ")

### Now we're going to download all fossil occurrences of species of horses present
### in this phylogeny

Occs <- pbdb_occurrences(taxon_name=Names, show=c("min_ma", "max_ma"), limit="all")

                   

