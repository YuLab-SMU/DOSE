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


