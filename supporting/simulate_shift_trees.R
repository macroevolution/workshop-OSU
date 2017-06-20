
# THis file is distributed as part of the "testing BAMM" page on the project website
#   but has the name "R_test_BAMM_functions.R"


require(geiger)
require(diversitree)


####-----------------------------------------
#### The following block of functions:
#
# simulateCBDPTtree
# CPBDPStochasticMap
# pruneCPBDPTree
# computeShiftDescendants
#
# were taken from the Dryad supplement to accompany Moore et al
# PNAS 2016: www.pnas.org/cgi/doi/10.1073/pnas.1518659113
# 
# Dryad doi: http://datadryad.org/resource/doi:10.5061/dryad.mb0sd
#
# We slightly modified SimulateCBDPTree such that it would 
#   generate a larger range of tree sizes by removing the 
#   size selection bias that contributed to ascertainment bias 
#   in the Moore et al article
 
SimulateCBDPTree <- function(time,
                             transition_rate,
                             lambda_function,
                             mu_function,
                             init_lambda = NULL,
                             init_mu = NULL,
                             verbose = TRUE,
                             condition_on_survival = TRUE,
                             condition_on_root = TRUE,
                             NMAX,
                             MAX_FAILS = 1,
                             cur_fails = 0) {

  # If initial lambda and mu are not provided, simulate them from the prior distribution
  if ( is.null(init_lambda) ) {
    initial_speciation <- lambda_function()
  } else {
    initial_speciation <- init_lambda
  }

  if ( is.null(init_mu) ) {
    initial_extinction <- mu_function()
  } else {
    initial_extinction <- init_mu
  }

  # Initialize the number of species
  lineage_index <- 3
  current_num_species <- 2
  total_number_of_transitions <- 0

  # Initialize the edge matrix
  edge <- data.frame(ancestor=1,descendant=2:3,start_time=0,end_time=NA,
                     current_speciation_rate=initial_speciation,current_extinction_rate=initial_extinction,
                     speciation_rates=I(list(initial_speciation,initial_speciation)),extinction_rates=I(list(initial_extinction,initial_extinction)),
                     transition_times=I(list(NA,NA)),status="alive",states=I(list(0,0)),stringsAsFactors=FALSE)

  # Start the simulation
  current_time <- 0

  if(verbose) bar <- txtProgressBar(style=3,width=40)
 

  fail_count = cur_fails

  while( TRUE ) {

    # get the next event time
    these_rates <- (edge$current_speciation_rate + edge$current_extinction_rate + transition_rate) * as.numeric(edge$status == "alive")
    next_times <- suppressWarnings(rexp(length(these_rates),these_rates))
    which_edge_has_event <- which.min(next_times)
    next_event_time <- next_times[which_edge_has_event]

    # increment the current time; if it is greater than max time, break
    current_time <- current_time + next_event_time

    if ( current_time > time ) {
      current_time <- time
      if(verbose) setTxtProgressBar(bar,current_time/time)
      break
    }

    if(verbose) setTxtProgressBar(bar,current_time/time)

    # perform the next event
    this_edge <- edge[which_edge_has_event,]
    these_rates <- c(this_edge$current_speciation_rate,this_edge$current_extinction_rate,transition_rate)
    total_rate <- sum(these_rates)
    probs <- these_rates / total_rate
    probs[these_rates == Inf] <- 1.0
    next_event_type <- sample.int(3,1,prob=probs)

    if ( next_event_type == 1 ) {

      # Speciation event
      edge[which_edge_has_event,]$status <- "speciated"
      edge[which_edge_has_event,]$end_time <- current_time

      this_ancenstor <- edge[which_edge_has_event,]$descendant
      these_states <- edge[which_edge_has_event,]$states[[1]]
      this_state <- these_states[length(these_states)]
      new_descendants <- 1:2 + lineage_index

      these_species <- data.frame(ancestor=this_ancenstor,descendant=new_descendants,start_time=current_time,end_time=NA,
                                  current_speciation_rate=these_rates[1],current_extinction_rate=these_rates[2],
                                  speciation_rates=I(list(these_rates[1],these_rates[1])),extinction_rates=I(list(these_rates[2],these_rates[2])),
                                  transition_times=I(list(NA,NA)),status="alive",states=I(list(this_state,this_state)),stringsAsFactors=FALSE)

      edge <- rbind(edge,these_species)

      current_num_species <- current_num_species + 1
      lineage_index <- lineage_index + 2

      # cat("current time: ",current_time," -- speciation event -- current species: ",current_num_species," -- number of transitions: ",total_number_of_transitions,"\n",sep="")

    } else if ( next_event_type == 2 ) {

      # Extinction event
      edge[which_edge_has_event,]$status <- "extinct"
      edge[which_edge_has_event,]$end_time <- current_time
      current_num_species <- current_num_species - 1

      # cat("current time: ",current_time," -- extinction event -- current species: ",current_num_species," -- number of transitions: ",total_number_of_transitions,"\n",sep="")

    } else if ( next_event_type == 3 ) {

      # Transition event
      total_number_of_transitions <- total_number_of_transitions + 1
      new_speciation_rate <- lambda_function()
      new_extinction_rate <- mu_function()

      these_speciation_rates <- this_edge$speciation_rates
      these_extinction_rates <- this_edge$extinction_rates
      these_transition_times <- this_edge$transition_times
      these_states           <- this_edge$states

      these_speciation_rates[[1]] <- c(these_speciation_rates[[1]],new_speciation_rate)
      these_extinction_rates[[1]] <- c(these_extinction_rates[[1]],new_extinction_rate)
      these_transition_times[[1]] <- c(these_transition_times[[1]],current_time)
      these_states[[1]]           <- c(these_states[[1]],total_number_of_transitions)

      this_edge$current_speciation_rate <- new_speciation_rate
      this_edge$current_extinction_rate <- new_extinction_rate
      this_edge$speciation_rates <- these_speciation_rates
      this_edge$extinction_rates <- these_extinction_rates
      this_edge$transition_times <- these_transition_times
      this_edge$states           <- these_states
      edge[which_edge_has_event,] <- this_edge

      # cat("current time: ",current_time," -- transition event -- current species: ",current_num_species," -- number of transitions: ",total_number_of_transitions,"\n",sep="")

    }

    # if the total number of species is 0, break
    if (current_num_species == 0) {
      if(verbose) {
        flush.console()
        if (verbose) cat("\nTree died!\n")
      }
   
      
      break
    }

    if ( current_num_species > NMAX ) {
      cat("Aborted / N = ", current_num_species, " Fails: ", fail_count, "\t" )
      
      fail_count <- fail_count + 1
      if (fail_count >= MAX_FAILS){
      	return(NA)
      }
      
      break
    }
 
  }

  if ( current_num_species < 1 & condition_on_survival ) {
    return(SimulateCBDPTree(time=time,transition_rate=transition_rate,
                             lambda_function=lambda_function,
                             mu_function=mu_function,
                             init_lambda=init_lambda,
                             init_mu=init_mu,
                             verbose=verbose,
                             condition_on_survival=condition_on_survival, NMAX=NMAX, MAX_FAILS = MAX_FAILS, cur_fails = fail_count))
  }

  if ( current_num_species > NMAX ) {
    return(SimulateCBDPTree(time=time,transition_rate=transition_rate,
                            lambda_function=lambda_function,
                            mu_function=mu_function,
                            init_lambda=init_lambda,
                            init_mu=init_mu,
                            verbose=verbose,
                            condition_on_survival=condition_on_survival, NMAX=NMAX, MAX_FAILS = MAX_FAILS, cur_fails = fail_count))
  }

  edge$end_time[is.na(edge$end_time)] <- current_time

  # Now make the tree ape-friendly
  phy_edge_length <- edge$end_time - edge$start_time
  phy_edge <- matrix(NA,nrow=nrow(edge),ncol=2)

  new_tips <- 1:sum(edge$status != "speciated")
  new_root <- max(new_tips) + 1
  new_nodes <- new_root + 1:(length(new_tips)-2)
  old_nodes <- edge$descendant[edge$status == "speciated"]

  phy_edge[,1] <- new_nodes[match(edge$ancestor,old_nodes)]
  phy_edge[,2] <- new_nodes[match(edge$descendant,old_nodes)]

  phy_edge[edge$status != "speciated",2] <- new_tips
  phy_edge[edge$ancestor == 1,1] <- new_root

  Nnode <- new_root - 2

  phy <- list()
  class(phy) <- "phylo"
  phy$edge <- phy_edge
  phy$edge.length <- phy_edge_length
  phy$tip.label <- new_tips
  phy$Nnode <- Nnode
  phy$full_process <- edge

  if ( condition_on_root ) {

    these_descendants <- phy$edge[phy$edge[,1] == Nnode + 2,2]
    left_descendants <- tips(phy,these_descendants[1])
    right_descendants <- tips(phy,these_descendants[2])
    if ( !any("alive" %in% edge[phy$edge[,2] %in% left_descendants,]$status) | !any("alive" %in% edge[phy$edge[,2] %in% right_descendants,]$status) ){
      if (verbose) cat("\nRoot died!\n")
      return(SimulateCBDPTree(time=time,transition_rate=transition_rate,
                              lambda_function=lambda_function,
                              mu_function=mu_function,
                              init_lambda=init_lambda,
                              init_mu=init_mu,
                              verbose=verbose,
                              condition_on_survival=condition_on_survival, NMAX = NMAX, MAX_FAILS = MAX_FAILS, cur_fails = fail_count))
    }
  }

  return(phy)

}

CPBDPStochasticMap <- function(tree) {

  edge <- tree$full_process

  # Now make the SIMMAP-formated transition events
  total_states <- 1 + max(unlist(edge$states))

  maps <- lapply(1:nrow(edge),function(x){
    row <- edge[x,]
    transition_times <- row$transition_times[[1]]
    durations <- diff(na.omit(c(row$start_time,transition_times,row$end_time)))
    names(durations) <- row$states[[1]]
    return(durations)
  })

  blank_map <- rep(0,total_states)
  names(blank_map) <- 1:total_states-1

  mapped.edge <- do.call(rbind,lapply(maps,function(map){
    for(i in 1:length(map)) blank_map[names(map[i])] <- blank_map[names(map[i])] + map[i]
    return(blank_map)
  }))

  tree$maps <- maps
  tree$mapped.edge <- mapped.edge

  return(tree)

}

pruneCPBDPTree <- function(tree) {

  edge <- tree$full_process

  # these nodes are already sorted from youngest to oldest
  # we reverse them to get the downpass sequence
  # exclude the root
  nodes <- rev(unique(edge$ancestor)[-1])

  for ( node in nodes ) {

    node_ancestor <- edge$ancestor[edge$descendant == node]
    these_edges <- edge[edge$ancestor %in% node,]
    if ( all(these_edges$status == "extinct") ) {
      # if both of the nodes descendants are extinct, drop all of them
      # and set the current node's status to extinct
      edge[edge$descendant %in% node,]$status <- "extinct"
      edge <- edge[!edge$ancestor %in% node,]
    } else if ( any(these_edges$status == "extinct") ) {
      # otherwise, just drop the extinct tip
      extinct_edge <- these_edges[these_edges$status == "extinct",]
      extinct_descendant <- extinct_edge$descendant

      # remove the extinct descendant
      edge <- edge[edge$descendant != extinct_descendant,]

      # relabel the ancestor and descendants
      alive_edge <- these_edges[these_edges$status != "extinct",]
      alive_descendant <- alive_edge$descendant
      which_recipient <- which(edge$descendant == node)
      which_donor <- which(edge$descendant == alive_descendant)
      edge[which_recipient,]$end_time <- edge[which_donor,]$end_time
      edge[which_recipient,]$status <- edge[which_donor,]$status

      if ( length(na.omit(edge[which_donor,]$transition_times[[1]])) > 0 ) {
        # if the donor branch has transitions, glue together the rates
        edge[which_recipient,]$speciation_rates <- list(c(edge[which_recipient,]$speciation_rates[[1]],edge[which_donor,]$speciation_rates[[1]][!is.na(edge[which_donor,]$transition_times[[1]])]))
        edge[which_recipient,]$extinction_rates <- list(c(edge[which_recipient,]$extinction_rates[[1]],edge[which_donor,]$extinction_rates[[1]][!is.na(edge[which_donor,]$transition_times[[1]])]))
        edge[which_recipient,]$current_speciation_rate <- edge[which_donor,]$current_speciation_rate
        edge[which_recipient,]$current_extinction_rate <- edge[which_donor,]$current_extinction_rate
        edge[which_recipient,]$states <- list(c(edge[which_recipient,]$states[[1]],edge[which_donor,]$states[[1]][-1]))
        edge[which_recipient,]$transition_times <- list(as.numeric(c(edge[which_recipient,]$transition_times[[1]],na.omit(edge[which_donor,]$transition_times[[1]]))))
      }

      # now, remove the donor
      edge <- edge[-which_donor,]
      edge$descendant[edge$descendant == node] <- alive_descendant

    }

  }

  # remember to do something with the root
  node <- 1
  these_edges <- edge[edge$ancestor %in% node,]

  if ( all(these_edges$status == "extinct") ) {
    message <- "The whole tree went extint. What happened??"
    cat(message)
    class(message) <- "try-error"
    return(invisible(message))
  } else if ( any(these_edges$status == "extinct" ) ) {
    # We need to remove the root
    warning("Removing root. Was this supposed to happen?")
    edge <- edge[!edge$ancestor %in% node,]
  }

  # Now make the tree ape-friendly
  smallest_node <- min(edge$ancestor)
  phy_edge_length <- edge$end_time - edge$start_time
  phy_edge <- matrix(NA,nrow=nrow(edge),ncol=2)
  new_tips <- 1:sum(edge$status != "speciated")
  new_root <- max(new_tips) + 1
  new_nodes <- new_root + 1:(length(new_tips)-2)
  old_nodes <- edge$descendant[edge$status == "speciated"]
  phy_edge[,1] <- new_nodes[match(edge$ancestor,old_nodes)]
  phy_edge[,2] <- new_nodes[match(edge$descendant,old_nodes)]
  phy_edge[edge$status != "speciated",2] <- new_tips
  phy_edge[edge$ancestor == smallest_node,1] <- new_root
  Nnode <- new_root - 2

  phy <- list()
  class(phy) <- "phylo"
  phy$edge <- phy_edge
  phy$tip.label <- new_tips
  phy$edge.length <- phy_edge_length
  phy$Nnode <- Nnode
  phy$full_process <- edge

  if(!is.null(tree$maps)){
    phy <- CPBDPStochasticMap(phy)
  }

  return(phy)


}

computeShiftDescendants <- function(tree) {

  edge <- tree$full_process
  edge$extant <- apply(edge,1,function(row){
    descendants <- tips(tree,tree$edge[edge$descendant == row$descendant,2])
    any(edge[tree$edge[,2] %in% descendants,]$status == "alive")
  })

  tree$full_process <- edge

  return(tree)

}


#### End MEA block of functions
####-----------------------------------------


#### End block of BAMM-testing code
####-----------------------------------------


# Function to extract basic shift regime data
# including numbers of tips in each rate regime
# and number of branches included in each rate regime
# from a bammdata object.
#  This only extracts from a single indexed shift 
#  configuration. 
shiftRegimeData <- function(ed, index = 1){
	
	ex <- ed$eventData[[index]]
	tipstates <- table(ed$tipStates[[index]])
	tipstates <- tipstates[as.character(ex$index)]
	
	edgestates <- table(ed$eventVectors[[1]])
	edgestates <- edgestates[as.character(ex$index)]
	
	ex$ntips      <- tipstates
	ex$nbranches  <- edgestates 
	return(ex)
}


colorByRegime <- function(ed){
	
	ev <- ed$eventVectors[[1]]
	cols <- sample(rainbow(n = length(unique(ev)) ))
	
	colset <- cols[ev]
	
	return(colset)
}


birthdeath_fit <- function(phy, node = NULL, dropnodes = NULL){
	
	droptips <- NULL
	
	if (!is.null(dropnodes)){
		
		for (i in 1:length(dropnodes)){
			
			if (i <= length(phy$tip.label)){
				droptips <- c(droptips, phy$tip.label[i])
			}else{
				tmp <- extract.clade(phy, node = dropnodes[i])$tip.label
				droptips <- c(droptips, tmp)
			}
			
			
		}
		
		phy <- drop.tip(phy, droptips)
		
	}
	if (!is.null(node)){
	     phy <- extract.clade(phy, node = node)	
	}
	
	lfx <- make.bd(phy)
	lfunc <- function(par){
		 -lfx(exp(par))
	} 
	
	
	ntips <- length(phy$tip.label)
	
	init <- 1.25 * (log(ntips) - log(2)) / max(branching.times(phy))
	init <- log(c(init, 0.5*init))
	names(init) <- c("lambda", "mu")
	res <- optim(par=init, lfunc, method="Nelder" )
	
	rates <- round(exp(res$par), 3)
	names(rates) <- c("lambda", "mu")
	
	return(rates)
	
}



plotSetup <- function(mx = 0.8){
	plot.new()
	par(oma=c(1,1,1,1))
	par(mar=c(5,5,1,1))
	plot.window(xlim=c(0, mx), ylim=c(0, mx), asp=1)	
	abline(0, 1, lwd=3, col="gray50")
	axis(1, at=seq(-0.2, 1.2, by=0.2))
	axis(2, at=seq(-0.2, 1.2, by=0.2), las=1)
	mtext("True rate", side=1, line=3, cex=1.4)
	mtext("BAMM estimated rate", side = 2, line = 3, cex = 1.4)
	
	lines(x=c(0,0.1), y=c(0.6, 0.6), lwd=3, col="gray50")
	text(x=0.1, y=0.6, pos=4, label = "1:1 line", font = 3)
	lines(x=c(0,0.1), y=c(0.55, 0.55), lty="dotted", lwd=3)
	text(x=0.1, y=0.55, pos=4, label = "BAMM fit", font=3)
		
	
}

# There is a huge amount of 
#  overplotting in the following plots, 
#   since most taxa have identical rates.
#  So we will "jitter" the points slightly so you can see
#  roughly how many branches fall out together

j <- function(x){
	return(jitter(x, amount = 0.015))
}
 

 


getSpanningTips <- function(phy, node){
	
	if (node <= length(phy$tip.label)){
		return(phy$tip.label[node])
	}else{
		
		dset <- phy$edge[,2][phy$edge[,1] == node]
		n1 <- dset[1]
		while (n1 > length(phy$tip.label)){
			n1 <- phy$edge[,2][phy$edge[,1] == n1][1]
		}
		n2 <- dset[2]
		while(n2 > length(phy$tip.label)){
			n2 <- phy$edge[,2][phy$edge[,1] == n2][1]
		}
		
		return(phy$tip.label[c(n1,n2)])
		
	}
	
	
}


# This function converts the events from a tree generated using the CPBDP simulator to 
#    BAMM's event data format

cpbdp_to_bammdata <- function(tree)	{
	# May 5th 2017: fixed by adding zero_Ns vector to *force* zero-tip shifts to be properly included
	# May 6th 2017: changed root-branch ID protocol to be based on start time
	# May 8th 2017: force the root regime as first row, even when no edge has root regime as current regime
	# May 9th 2017: redid loop to focus on each transition, instead of looping through states, to avoid
	#                 annoying property of how edges with multiple transitions are recorded
	require(phytools)
	N <- Ntip(tree)
	
	# Find which edge each transition occurs on
	edgeTimes <- sapply(tree$full_process$transition_times, length)
	shiftEdges <- which(edgeTimes > 1)
	
	# identify all regimes
	States <- sapply(tree$full_process$states, function(x) x[length(x)])
	edgeStates <- setNames(States, tree$edge[,2])
	tipStates <- edgeStates[as.character(1:N)]
	allStates <- unique(unlist(tree$full_process$states))
	
	# how many tips in each regime?
	stateN <- tapply(tipStates, tipStates, length)
	toAdd <- setdiff(allStates, names(stateN))
	zero_Ns <- rep(0, length(toAdd))
	names(zero_Ns) <- toAdd
	stateN <- c(stateN, zero_Ns)
	stateN <- stateN[order(as.numeric(names(stateN)))]
	
	# store root & start event data output
	rootN <- stateN["0"]	
	Root <- N + 1
	rootTips <- getSpanningTips(tree, Root)
 	eventData <- data.frame(generation=0, leftchild=rootTips[1], rightchild=rootTips[2], abstime=0, 
 							lambdainit=tree$full_process$speciation_rates[[1]][1], lambdashift=0, 
 							muinit=tree$full_process$extinction_rates[[1]][1], mushift=0, n_taxa=rootN)
			 	
	# go through each edge with at least one rate shift
	for (regime in shiftEdges)	{
		
		# what node does this edge correspond to?
		Node <- tree$edge[regime, 2]
		
		# how many transitions occur on this edge?
		Len <- length(tree$full_process$states[regime][[1]])
		if (Len == 1)	{
			cat("No transition on edge ", regime, "\n")
			stop()
		}
		
		# loop through all transitions on a given edge and save edata
		for (j in 2:Len)		{
			abstime <- tree$full_process$transition_times[[regime]][j]
			Spec <- tree$full_process$speciation_rates[[regime]][j]
			Ext <- tree$full_process$extinction_rates[[regime]][j]
			State <- tree$full_process$states[[regime]][j]
			if (Node >= N)	{
				Descs <- getSpanningTips(tree, Node)
				lchild <- as.character(Descs[1])
				rchild <- as.character(Descs[length(Descs)])
				eventData <- rbind(eventData, 
						data.frame(generation=0, leftchild=lchild, rightchild=rchild, abstime=abstime, 
						lambdainit=Spec, lambdashift=0, muinit=Ext, 
						mushift=0, n_taxa=stateN[as.character(State)]))
			}
			else	 if (Node < N){
				eventData <- rbind(eventData, 
						 data.frame(generation=0, leftchild=as.character(tree$tip.label[Node]), rightchild="NA", 
						 abstime=abstime, lambdainit=Spec, lambdashift=0, muinit=Ext, mushift=0, 
						 n_taxa=stateN[as.character(State)]))
			}
		}	
		#cat(regime, "\n")		
	}

	eventData <- eventData[order(eventData$abstime),]
	if (is.na(sum(eventData$abstime)))	{
		cat("Some tranisition time is NA!\n")
		stop()
	}
	return(eventData)
}




simulateShiftTree <- function(tmax, shiftrate, lamfx, mufx, seed = NULL, NMAX = 1000, MAX_FAILS = 1, NMIN = 25){
	
 
	
	if (!is.null(seed)){
		set.seed(seed)
	}
	
	# Using MEA tree simulator to generate the tree:
	tree <- SimulateCBDPTree(tmax, transition_rate=shiftrate, lamfx, mufx, verbose=F, NMAX = NMAX, MAX_FAILS = MAX_FAILS)

	
	if (length(tree) == 1){
		res <- list(phy = NA, events = NA)
		return(res)
	}
	
	
	
	# Map shifts onto tree
	simTree <- CPBDPStochasticMap(tree)
 
	# Prune extinct taxa 
	prunedTree <- pruneCPBDPTree(simTree)
	
	if (length(prunedTree$tip.label) < NMIN){
		res <- list(phy = NA, events = NA)
		return(res)	
	}
	
 
	# Convert the regime mappings from the simulation to BAMM-style eventData format
	edata <- cpbdp_to_bammdata(prunedTree)[,1:8]	
	
	prunedTree <- read.tree(text = write.tree(prunedTree))

	rownames(edata) <- NULL
	
	res <- list(phy = prunedTree, events = edata)
	
	return(res)
	
}






