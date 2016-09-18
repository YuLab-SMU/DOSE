##' cnetplot
##'
##' 
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
##' @author Guangchuang Yu \url{http://guangchuangyu.github.io}
setGeneric("cnetplot",
           function(x, showCategory=5, categorySize="geneNum", foldChange=NULL, fixed=TRUE, ...)
               standardGeneric("cnetplot"))


##' dotplot
##'
##' 
##' @docType methods
##' @name dotplot
##' @rdname dotplot-methods
##' @title dotplot method
##' @param ... additional parameter
##' @return plot
##' @export
##' @author Guangchuang Yu
setGeneric("dotplot", function(object, ...) standardGeneric("dotplot"))


##' upsetplot method generics
##'
##'
##' @docType methods
##' @name upsetplot
##' @rdname upsetplot-methods
##' @title upsetplot method
##' @param x object
##' @param ... additional parameters
##' @return plot
##' @export
setGeneric("upsetplot", function(x, ...) standardGeneric("upsetplot"))

#' geneID generic
#'
#' @param x enrichResult
geneID <- function(x) {
   UseMethod("geneID", x)
}

#' geneInCategory generic
#'
#' @param x enrichResult
geneInCategory <- function(x) {
   UseMethod("geneInCategory", x)
}
