

simulateTree <- function(pars, max.taxa = Inf, max.t=Inf, min.taxa = 2, include.extinct=F){
 
	badcount <- 0;
	while (1){
		
		tree <- tree.bd(pars, max.taxa=max.taxa, max.t=max.t, include.extinct=include.extinct);
		if (!is.null(tree)){
			if (length(tree$tip.label) >= min.taxa){
				break;
			}
		}
		
		badcount <- badcount + 1;
		if (badcount > 200){
			stop("Too many trees going extinct\n");
		}
	}
	tree$node.label <- NULL;
	return(tree);
}

# Computes the Colless imbalance statistic 
#	across an entire tree.
colless <- function(phy){
	
	bb <- balance(phy);
	ss <- sum(abs(bb[,1] - bb[,2]));
	n <- length(phy$tip.label);
	return((2 / ((n-1)*(n-2))) * ss);
	
}

logit <- function(x, min=0, max=1){
	p <- (x-min)/(max-min)
	log(p/(1-p))
}

invlogit <- function(x, min=0, max=1)
{
	p <- exp(x)/(1+exp(x))
	p <- ifelse( is.na(p) & !is.na(x), 1, p ) # fix problems with +Inf
	p * (max-min) + min
}

# Gets a vector of initial parameters 
#	for a bisse, birth-death, or mk2 model
#	using the likelihood functions created 
#	in diversitree with make.mk2, make.bisse, or make.bd
# These can be plugged directly into the corresponding
#	likelihood function. However, they are not guaranteed
#	to generate finite log-likelihoods. 
# Arguments:
#	fx: the diversitree likelihood function
#	lmin: the minimum value across all parameters
#	lmax: the maximum value across all parameters

getStartingParamsDiversitree <- function(fx, lmin, lmax){
		
		lamset <- runif(3, lmin, lmax);
		names(lamset) <- c('lambda', paste('lambda', 0:1, sep=''));
		muset <- runif(3, 0, 1) * lamset;
		names(muset) <- c('mu', paste('mu', 0:1, sep=''));
		qset <- runif(4, lmin, lmax * 0.2);
		names(qset) <- c('q01', 'q10', 'q12', 'q21');
		parvec <- c(lamset, muset, qset);
	 
		if (length(setdiff(argnames(fx), names(parvec))) > 0){
			stop("Invalid argnames from function\n");
		}
		
		parset <- intersect(names(parvec), argnames(fx));
				
		return(parvec[parset]);
}


# A general purpose optimization function
# that optimizes parameters of a diversitree likelihood function.
# The likelihood function must correspond to one of the following models:
#	a) BiSSE (or any constrained submodel)
#	b) birth-death
#	c) mk2 (2 state character only model)
fitDiversitree <- function(fx, nopt=1, lmin = 0.001, lmax=20.0, MAXBAD = 100, initscale = 0.1){

	
	for (i in 1:nopt){
		
		badcount <- 0;
		
		iv <- getStartingParamsDiversitree(fx, lmin=lmin, lmax=lmax*initscale);
		
		resx <- try(optim(iv ,fx, method='L-BFGS-B', control=list(maxit=1000, fnscale=-1), lower=lmin, upper=lmax), silent=T);
		while (class(resx) == 'try-error'){
		iv <- getStartingParamsDiversitree(fx, lmin=lmin, lmax=lmax*initscale);
			resx <- try(optim(iv , fx, method='L-BFGS-B', control=list(maxit=1000, fnscale=-1), lower=lmin, upper=lmax), silent=T);
			
			badcount <- badcount + 1;
			if (badcount > MAXBAD){
				stop("Too many fails in fitDiversitree\n");
			}
		}
		
		if (i == 1){
			best <- resx;
		}else{
			if (best$value < resx$value){
				best <- resx;
			}
		}
		
	}
	
	fres <- list(pars=best$par, loglik=best$value);
	fres$AIC <- -2*fres$loglik + 2*length(argnames(fx));
	fres$counts <- best$counts;
	#fres$like_function <- fx;
	fres$convergence <- best$convergence;
	fres$message <- best$message;
	return(fres);
}


# simDiscrete: a simple wrapper function for 
#	geiger::sim.char. Allows user to specify minimum
#	and maximum frequencies of the rare character state.
#	Arguments:
#		phy: 	phylogenetic tree, ape format
#		q:		rate (only allows symmetric rates)
#		minf:	minimum frequency of rarer state
#		maxf:	maximum frequency of rarer state
#		root:	root character state
#	Returns: A vector of character states, coded 
#	as 1 or 0 

simDiscrete <- function(phy, q, minf = 0.05, maxf = 0.5, root=0){
	MAXBADCOUNT <- 200;
	require(geiger);
	mm <- matrix(rep(q, 4), nrow=2);
	diag(mm) <- -1 * diag(mm);
	
	bad <- TRUE;
	badcount <- 0;
	while (bad){
		chars <- sim.char(v, par=mm, nsim=1, model='discrete', root=root+1)[,1,1];
		tx <- table(chars);
		if (length(tx) == 2){
			ff <- min(tx)/length(chars);
			if (ff >= minf & ff <= maxf){
				bad <- FALSE;
			}
		}
		if (badcount > MAXBADCOUNT){
			stop("Exceeded MAXBADCOUNT in simDiscrete\n");
		}
		badcount <- badcount + 1;
	}
	return(chars - 1);
	
}


fitCRBD <- function(phy, nopt=5, lmin=0.001, lmax=5.0, MAXBAD = 200){
	
	if (length(phy$tip.label) < 3){
		pars <- c(0.0001,0)
		names(pars) <- c("lambda", "mu")
		return(pars)
	}
	
	fx <- make.bd(phy)
	
	for (i in 1:nopt){
	
		lam <- runif(1, 0, 0.5)	
	 	mu <- lam * runif(1, 0, 1)
	 
		badcount <- 0
 
		resx <- try(optim(c(lam, mu) ,fx, method='L-BFGS-B', control=list(maxit=1000, fnscale=-1), lower=lmin, upper=lmax), silent=T)
		while (class(resx) == 'try-error'){

			lam <- runif(1, 0, 0.5)	
	 		mu <- lam * runif(1, 0, 1)
			
			resx <- try(optim(c(lam, mu) , fx, method='L-BFGS-B', control=list(maxit=1000, fnscale=-1), lower=lmin, upper=lmax), silent=T);
			
			badcount <- badcount + 1;
			if (badcount > MAXBAD){
				stop("Too many fails in fitDiversitree\n")
			}
		}
		
		if (i == 1){
			best <- resx
		}else{
			if (best$value < resx$value){
				best <- resx
			}
		}
		
	}
	
	fres <- list(pars=best$par, loglik=best$value)
	fres$AIC <- -2*fres$loglik + 2*length(argnames(fx))
	fres$counts <- best$counts
	#fres$like_function <- fx
	fres$convergence <- best$convergence
	fres$message <- best$message
	
	pars <- fres$pars
	names(pars) <- c("lambda", "mu")
	
	return(pars)
}



loglike_purebirth <- function(lambda, phy){
	
	n <- length(phy$tip.label)
	ll <- (n - 2) * log(lambda) - lambda*sum(phy$edge.length) 
	return(ll)
}



colorFunction <- function(x, minn, maxx, colorset){
	
 	
	sq <- seq(from=minn, to = maxx, length.out = length(colorset))
	cols <- rep(colorset[1], length(x))
	
	cols[x <= minn] <- colorset[1]
	cols[x >= maxx] <- colorset[length(colorset)]
	
	for (i in 2:length(sq)){
		isIn <- x >= sq[i-1] & x < sq[i]
		cols[isIn] <- colorset[i]
	}
	cols[is.na(x)] <- NA
	return(cols)

}









 
