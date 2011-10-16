.initial <- function() {
	assign("DOSEEnv", new.env(),.GlobalEnv)
	assign("SemSimCache", new.env(), .GlobalEnv)
	
	tryCatch(utils::data(list="DO2EG", package="DOSE"))
	assign("DO2EG", DO2EG, envir=DOSEEnv)
	
	tryCatch(utils::data(list="EG2DO", package="DOSE"))
	assign("EG2DO", EG2DO, envir=DOSEEnv)
	
}



geneSim <- function(geneID1, geneID2, method="Wang", ont="DO", organism="human", combine="rcmax.avg") {
	DOID1 <- sapply(geneID1, gene2DO)
	DOID2 <- sapply(geneID2, gene2DO)
	m <- length(geneID1)
	n <- length(geneID2)
	scores <- matrix(NA, nrow=m, ncol=n)
	rownames(scores) <- geneID1
	colnames(scores) <- geneID2
	
	for (i in 1:m) {
		for (j in 1:n) {
			if(any(!is.na(DOID1[[i]])) &&  any(!is.na(DOID2[[j]]))) {
				s <- doSim(DOID1[[i]], DOID2[[j]], method, ont, organism)
				scores[i,j] = .combineScores(s, combine)
			}
		}
	}
	return(scores)
}

doSim <- function(DOID1, DOID2, method="Wang", ont="DO", organism="human") {
	m <- length(DOID1)
	n <- length(DOID2)
	scores <- matrix(nrow=m, ncol=n)
	rownames(scores) <- DOID1
	colnames(scores) <- DOID2
	for( i in 1:m) {
		for (j in 1:n) {
			if ( is.na(DOID1[i]) || is.na(DOID2[j]) ) {
					scores[i,j] <- NA
			} else {
				if (method == "Wang") {
					scores[i,j] <- .wangMethod(DOID1[i], DOID2[j], ont=ont)
				} else {
					scores[i,j] <- .infoContentMethod(DOID1[i], DOID2[j], ont=ont, method=method, organism=organism)
				}
			}
		}
	}	
	return(scores)
}


###########################################################
## Method *Wang* for semantic similarity measuring ###
###########################################################
.wangMethod <- function(ID1, ID2, ont="DO", ONTPARENTS=DOPARENTS, weight.isa=0.8, weight.partof=0.6, weight.do=0.7) {
	
	if (ID1 == ID2)
		return (sim=1)		

	sv.a <- 1
	sv.b <- 1
	sw <- 1
	names(sv.a) <- ID1
	names(sv.b) <- ID2 
	
	Parents <- ONTPARENTS
	sv.a <- .SemVal(ID1, ont, Parents, sv.a, sw, weight.isa, weight.partof, weight.do)
	sv.b <- .SemVal(ID2, ont, Parents, sv.b, sw, weight.isa, weight.partof, weight.do)
	
	sv.a <- .uniqsv(sv.a)
	sv.b <- .uniqsv(sv.b)
	
	idx <- intersect(names(sv.a), names(sv.b))
	inter.sva <- unlist(sv.a[idx])
	inter.svb <- unlist(sv.b[idx])
	if (is.null(inter.sva) || is.null(inter.svb) || length(inter.sva) == 0 || length(inter.svb) ==0) {
		sim <- NA
	} else {
		sim <- sum(inter.sva,inter.svb) / sum(sv.a, sv.b)
	}
	return(sim)
}

.uniqsv <- function(sv) {
	sv <- unlist(sv)
	una <- unique(names(sv))
	sv <- unlist(sapply(una, function(x) {max(sv[names(sv)==x])}))
	return (sv)
}

.SemVal_internal <- function(ID, ont, Parents, sv, w, weight.isa, weight.partof, weight.do) {
	if (!exists(ID, Parents)) {
		return(NA)
	}
	p <- get(ID, Parents)
	#p <- unlist(p[[1]])
	if (length(p) == 0 || is.na(p)) {
		#warning(ID, " may not belong to Ontology ", ont)
		return(NA)
	}
	
	old.w <- w
	if (ont == "DO") {
		topNode <- "DOID:4"
	} else {
		relations <- names(p)
		topNode <- "all"
	}

	for (i in 1:length(p)) {
		if (ont == "DO") {
			w <- old.w * weight.do
		} else {
			if (relations[i] == "is_a") {
				w <- old.w * weight.isa
			} else {
				w <- old.w * weight.partof
			}
		}
		names(w) <- p[i]
		sv <- c(sv,w)
		if (p[i] != topNode) {
			sv <- .SemVal_internal(p[i], ont, Parents, sv, w, weight.isa, weight.partof, weight.do)
		}
	}
	return (sv)
}

.SemVal <- function(ID, ont, Parents, sv, w, weight.isa, weight.partof, weight.do) {
#	if(!exists("SemSimCache")) return(.SemVal_internal(ID, ont, Parents, sv, w, weight.isa, weight.partof, weight.do))
	if(!exists("SemSimCache")) {
		.initial()
	}
	ID.ont <- paste(ID, ont, sep=".")
	if (!exists(ID.ont, envir=SemSimCache)) {
	  	value <- .SemVal_internal(ID, ont, Parents, sv, w, weight.isa, weight.partof, weight.do)
	  	assign(ID.ont, value, envir=SemSimCache)
		#cat("recompute ", ID, value, "\n")
	}
	else{
		#cat("cache ", ID, get(ID, envir=SemSimCache), "\n")
	}
	return(get(ID.ont, envir=SemSimCache))
}


computeIC <- function(ont="DO", organism="human"){
		##require(DO.db)
		doids <- toTable(DOTERM)
		doterms <- doids$do_id
		docount <- table(doterms)
		doids <- names(docount) ##unique(doterms)
		cnt <- sapply(doids,function(x){ n=docount[get(x, DOOFFSPRING)]; docount[x]+sum(n[!is.na(n)])})		
		names(cnt) <- doids	
		p <- cnt/sum(docount) 
		# IC of DO terms was quantified as the negative log likelihood. 	
		IC <- -log(p)
		fname <- paste(paste("Info_Contents", organism, ont, sep="_"), ".rda", sep="")
		save(IC, file=fname)
}

loadICdata <- function(organism, ont) {
	if(!exists("DOSEEnv")) .initial()
	fname <- paste("Info_Contents", organism, ont,  sep="_")
	tryCatch(utils::data(list=fname, package="DOSE"))
	IC <- get("IC")
	org.ont.IC <- paste(organism, ont, "IC", sep="")
	assign(eval(org.ont.IC), IC, envir=DOSEEnv)
	rm (IC)
}

###########################################################
## Information Content Based Methods for semantic similarity measuring ###
###########################################################
.infoContentMethod <- function(ID1, ID2, ont="DO", ONTANCESTOR=DOANCESTOR, method, organism="human") {
	if(!exists("DOSEEnv")) {
		.initial()
	} 

	org.ont.IC <- paste(organism, ont, "IC", sep="")
	if(!exists(org.ont.IC, envir=DOSEEnv)) {
		loadICdata(organism, ont)
	}	
	IC <- get(org.ont.IC, envir=DOSEEnv)
	
	# more specific term, larger IC value.
	# Normalized, all divide the most informative IC.
	# all IC values range from 0(root node) to 1(most specific node)	
	mic <- max(IC[IC!=Inf])	
	
	if (ont == "DO") {
		topNode <- "DOID:4"
	} else {
		topNode <- "all"
	}
	
	IC[topNode] = 0	
	
	ic1 <- IC[ID1]/mic
	ic2 <- IC[ID2]/mic
	
	if (ic1 == 0 || ic2 == 0) return (NA)
				
	ancestor1 <- get(ID1, ONTANCESTOR)
	ancestor2 <- get(ID2, ONTANCESTOR)
	if (ID1 == ID2) {
		commonAncestor <- ID1
	} else if (ID1 %in% ancestor2) {
		commonAncestor <- ID1
	} else if (ID2 %in% ancestor1) {
		commonAncestor <- ID2
	} else { 
		commonAncestor <- intersect(ancestor1, ancestor2)
	}
	if (length(commonAncestor) == 0) return (NA)
	
	#Information Content of the most informative common ancestor (MICA)
	mica <- max(IC[commonAncestor])/mic  
	
	## IC is biased
	## because the IC of a term is dependent of its children but not on its parents.
	sim <- switch(method,
   	    Resnik = mica, ## Resnik does not consider how distant the terms are from their common ancestor.
   	    ## Lin and Jiang take that distance into account.
   	    Lin = 2*mica/(ic1+ic2),
   	    Jiang = 1 - min(1, -2*mica + ic1 + ic2), 
   	    Rel = 2*mica/(ic1+ic2)*(1-exp(-mica*mic))  ## mica*mic equals to the original IC value. and exp(-mica*mic) equals to the probability of the term's occurence. 
	)   	
	return (sim)
}

###########################################################
## Function for combine scores ###
###########################################################
.combineScores <- function(SimScores, combine) {

	if (length(combine) == 0) {  #if not define combine
		return(round(SimScores, digits=3))
	} else {
	}
	
	## combine was define...	
    if(!sum(!is.na(SimScores))) return (NA)
    if (is.vector(SimScores) || nrow(SimScores)==1 || ncol(SimScores)==1) {
        if (combine == "avg") {
            return(round(mean(SimScores), digits=3))
        } else {
            return (round(max(SimScores), digits=3)) 
        }
    }
    if (combine == "avg") {
		result <- mean(SimScores, na.rm=TRUE)
	} else if (combine == "max") {
		result <- max(SimScores, na.rm=TRUE)
	} else if (combine == "rcmax") {
		rowScore <- mean(apply(SimScores, 1, max, na.rm=TRUE))
		colScore <- mean(apply(SimScores, 2, max, na.rm=TRUE))
		result <- max(rowScore, colScore)
	} else if (combine == "rcmax.avg") {
		result <- sum( apply(SimScores, 1, max, na.rm=TRUE), apply(SimScores, 2, max, na.rm=TRUE) ) / sum(dim(SimScores))
	}
	
	return (round(result, digits=3))
}

gene2DO <- function(gene) {
	if(!exists("DOSEEnv")) .initial()
	EG2DO <- get("EG2DO", envir=DOSEEnv)
	DO <- EG2DO[[gene]]
	if (is.null(DO)) {
		return(NA)
	} 
	if (sum(!is.na(DO)) == 0) {
		return(NA)
	}
	if (length(DO) == 0) {
		return(NA)
	}
	return(DO)
}

HyperG <- function(numWdrawn, numW, numB, numDrawn) {
	#numWdrawn: number of White balls drawn
	#numW: number of White balls
	#numB: number of Black balls
	#numDrawn: number of balls drawn
	pvalue <- phyper(numWdrawn, numW, numB, numDrawn, lower.tail=FALSE)
	return(pvalue)
}

.yPaste <- function(a, b) {
	x=paste(a, "/", b, sep="", collapse="")
	return(x)
}
