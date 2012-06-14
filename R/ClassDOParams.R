##' Class "DOParams"
##' This class contains parameters for calculating DO semantic similarity among DO term or Gene list.
##'
##'
##' @name DOParams-class
##' @aliases DOParams-class
##'   sim,DOParams-method
##'
##' @docType class
##' @slot IDs containing a list of DO terms or Gene IDs.
##' @slot type specify the type of IDs, one of "DOID", "GeneID".
##' @slot method Method for calculating DO semantic similarity, one of "Resnik", "Jiang", "Lin", "Rel", "Wang".
##' @slot combine Method for combining DO semantic similarity scores, one of "avg", "max", "rcmax", "rcmax.avg"
##' @slot dropCodes dropCodes for mapping Gene to DO Terms.
##' @slot organism currently, only "human" supported.
##' @exportClass DOParams
##' @author Guangchuang Yu \url{http://ygc.name}
##' @seealso \code{\link{sim}}
##' @keywords classes
setClass(Class="DOParams",
         representation(
                        IDs = "list",
                        type = "character",
                        method="character",
                        combine="character",
                        dropCodes="character",
                        organism="character"
                        ),
         prototype=prototype (
         dropCodes="NULL",
         organism="human")
         )



setValidity("DOParams",
            function(object) {
		organs <- c("human")  ## only human supported.
		types <- c("DOID", "GeneID")
		mets <- c("Resnik", "Jiang", "Lin", "Rel", "Wang")
		combines <- c("max", "avg", "rcmax", "rcmax.avg", "BMA")
		if (length(object@IDs) != 2) {
                    warning("IDs must set to a list of length 2.")
		}
		##if (object@ontology %in% onts) {

		##} else {
                    ##stop("*ontology* must be one of ", paste(onts, collapse=","))
		##}
		if(length(object@organism != 0)) {
                    if (object@organism %in% organs) {

                    } else {
                        stop("*organism* must be one of ",
                             paste(organs, collapse=","))
                    }
		} else {
		}
		if (object@method %in% mets) {

		} else {
                    stop("*method* must be one of ",
                         paste(mets, collapse=","))
		}
		if (object@type %in% types) {

		} else {
                    stop("*type* must be one of ",
                         paste(types, collapse=","))
		}
		if(length(object@combine) != 0) {
                    if (object@combine %in% combines) {

                    } else {
                        stop("*combine* must be one of ",
                             paste(combines, collapse=","))
                    }
		} else {

		}
		return (TRUE)
            }
            )

##' sim method for \code{DOParams} instance
##'
##'
##' @name sim
##' @docType methods
##' @rdname sim-methods
##'
##' @title Methods for calculating semantic similarity
##' @param params A \code{DOParams} instance.
##' @return Semantic similarity value or matrix.
##' @importFrom GOSemSim combineScores
##' @importFrom DO.db DOPARENTS
##' @importFrom DO.db DOANCESTOR
##' @author Guangchuang Yu \url{http://ygc.name}
setMethod(
          f= "sim",
          signature= "DOParams",
          definition=function(params){
              if (params@type == "DOID") {
                  scores <- doSim(DOID1=params@IDs[[1]],
                                  DOID2=params@IDs[[2]],
                                  method=params@method,
                                  organism=params@organism)
                  if (length(params@combine) == 0) {
                      result <- round(scores, digits=3)
                  } else {
                      result <- combineScores(scores, params@combine)
                  }
              } else {  ## type == "GeneID"
                  if (length(params@combine) == 0) {
                      stop("*combine* must be setting for combining semantic similarity scores of multiple DO terms. ")
                      ##\nUsing setCombineMethod(\"Params\") to specify which method to combine.
                  }
                  result <- geneSim(geneID1=params@IDs[[1]],
                                    geneID2=params@IDs[[2]],
                                    method=params@method,
                                    organism=params@organism,
                                    combine=params@combine)
              }
              return(result)
          }
          )
