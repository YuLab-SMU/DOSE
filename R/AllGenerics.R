##' @docType methods
##' @name cnetplot
##' @rdname cnetplot-methods
##' @title cnetplot method
##' @param x enrichResult object
##' @param showCategory number of category plotted
##' @param categorySize one of geneNum or pvalue
##' @param foldChange fold change of expression value
##' @param fixed logical
##' @param ... additional parameters
##' @return plot
##' @export
##' @author Guangchuang Yu \url{http://ygc.name}
setGeneric("cnetplot",
           function(x, showCategory=5, categorySize="geneNum", foldChange=NULL, fixed=TRUE, ...)
               standardGeneric("cnetplot"))



##' @docType methods
##' @name dotplot
##' @rdname dotplot-methods
##' @title dotplot method
##' @param ... additional parameter
##' @return plot
##' @export
##' @author Guangchuang Yu
setGeneric("dotplot", function(object, ...) standardGeneric("dotplot"))

## ## @exportMethod "setReadable<-"
## setGeneric(
##            name="setReadable<-",
##            def=function(x, value) {standardGeneric("setReadable<-")}
##            )

##' preparing geneSets for gene set enrichment analysis
##'
## @S3method getGeneSet DO
##' @title getGeneSet
##' @param setType type of gene sets
##' @param organism organism
##' @param ... additional parameter
##' @export
getGeneSet <- function(setType, organism, ...) {
    UseMethod("getGeneSet")
}

##' Mapping External ID to Ontology Term ID
##'
## @S3method EXTID2TERMID DO
##' @title EXTID2TERMID
##' @param gene gene ID vector
##' @param organism organism
##' @param ... additional parameter
##' @export
EXTID2TERMID <- function(gene, organism, ...) {
    UseMethod("EXTID2TERMID")
}

##' Mapping Ontology Term ID to External ID
##'
## @S3method TERMID2EXTID DO
##' @title TERMID2EXTID
##' @param term term ID vector
##' @param organism organism
##' @param ... additional parameter
##' @export
TERMID2EXTID <- function(term, organism, ...) {
    UseMethod("TERMID2EXTID")
}

##' Get all background External ID.
##'
## @S3method ALLEXTID DO
##' @title ALLEXTID
##' @param organism organism
##' @param ... additional parameter
##' @export
ALLEXTID <- function(organism, ...) {
    UseMethod("ALLEXTID")
}

##' Mapping Ontology Term ID to Name Symbol or Description
##'
## @S3method TERM2NAME DO
##' @title TERM2NAME
##' @param term term ID vector
##' @param organism organism
##' @param ... additional parameter
##' @export
TERM2NAME <- function(term, organism, ...) {
    UseMethod("TERM2NAME")
}
