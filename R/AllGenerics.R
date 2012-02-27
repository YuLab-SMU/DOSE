##' plot method generics
##'
##'
##' @docType methods
##' @name plot
##' @rdname plot-methods
##' @title plot method
##' @param ... Additional argument list
##' @return plot
##' @importFrom BiocGenerics plot
##' @export
##' @author Guangchuang Yu \url{http://ygc.name}
if ( !isGeneric("plot") )
	setGeneric("plot", function(x, ...) standardGeneric("plot"))
	
##' @exportMethod sim
setGeneric(
           name = "sim",
           def=function(params){standardGeneric("sim")}
           )

##' Mapping External ID to Ontology Term ID
##'
##' @S3method EXTID2TERMID DO
##' @export
##' @param gene gene ID vector
##' @param organism organism
EXTID2TERMID <- function(gene, organism) {
    UseMethod("EXTID2TERMID")
}

##' Mapping Ontology Term ID to External ID
##'
##' @S3method TERMID2EXTID DO
##' @export
##' @param term term ID vector
##' @param organism organism
TERMID2EXTID <- function(term, organism) {
    UseMethod("TERMID2EXTID")
}

##' Get all background External ID.
##'
##' @S3method ALLEXTID DO
##' @export
##' @param organism organism
ALLEXTID <- function(organism) {
    UseMethod("ALLEXTID")
}

##' Mapping Ontology Term ID to Name Symbol or Description
##'
##' @S3method TERM2NAME DO
##' @export
##' @param term term ID vector
TERM2NAME <- function(term) {
    UseMethod("TERM2NAME")
}

##' Mapping External ID to Name Symbol.
##'
##' @S3method EXTID2NAME DO
##' @export
##' @param geneID gene ID vector
##' @param organism organism
EXTID2NAME <- function(geneID, organism) {
    UseMethod("EXTID2NAME")
}
