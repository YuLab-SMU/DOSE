setValidity("DOParams",
	function(object) {
		organs <- c("human")  ## only human supported.
		onts <- c("DO")  ## 
		types <- c("DOID", "GeneID")
		mets <- c("Resnik", "Jiang", "Lin", "Rel", "Wang")
		combines <- c("max", "avg", "rcmax", "rcmax.avg")
		if (length(object@IDs) != 2) {
			warning("IDs must set to a list of length 2.")
		}
		if (object@ontology %in% onts) {
		
		} else {
			stop("*ontology* must be one of ", paste(onts, collapse=","))
		}
		if(length(object@organism != 0)) {
			if (object@organism %in% organs) {
			
			} else {
				stop("*organism* must be one of ", paste(organs, collapse=","))
			}
		} else {
		}
		if (object@method %in% mets) {
			
		} else {
			stop("*method* must be one of ", paste(mets, collapse=","))
		}
		if (object@type %in% types) {
			
		} else {
			stop("*type* must be one of ", paste(types, collapse=","))
		}
		if(length(object@combine) != 0) {
			if (object@combine %in% combines) {
				
			} else {
				stop("*combine* must be one of ", paste(combines, collapse=","))
			}
		} else {
		
		}
		return (TRUE)
	}
)

setMethod(
	f= "sim", 
	signature= "DOParams", 
	definition=function(params){
		if (params@type == "DOID") {
			scores <- doSim(DOID1=params@IDs[[1]], DOID2=params@IDs[[2]], method=params@method,ont=params@ontology, organism=params@organism) 
			if (length(params@combine) == 0) {
				result <- round(scores, digits=3)
			} else {
				result <- .combineScores(scores, params@combine)
			}
		} else {  ## type == "GeneID"
			if (length(params@combine) == 0) {
				stop("*combine* must be setting for combining semantic similarity scores of multiple DO terms. ") ##\nUsing setCombineMethod(\"Params\") to specify which method to combine.
			}
			result <- geneSim(geneID1=params@IDs[[1]], geneID2=params@IDs[[2]], method=params@method, ont=params@ontology, organism=params@organism, combine=params@combine)
		}
		return(result)
	}
)


