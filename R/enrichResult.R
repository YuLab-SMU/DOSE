


## ##' plot method generics
## ##'
## ##'
## ##' @docType methods
## ##' @name plot
## ##' @rdname plot-methods
## ##' @aliases plot,enrichResult,ANY-method
## ##' @title plot method
## ##' @param x A \code{enrichResult} instance
## ##' @param type one of bar, cnet or enrichMap
## ##' @param ... Additional argument list
## ##' @return plot
## ##' @importFrom stats4 plot
## ##' @exportMethod plot
## ##' @author Guangchuang Yu \url{http://guangchuangyu.github.io}
## setMethod("plot", signature(x="enrichResult"),
##           function(x, type = "bar", ... ) {
##               if (type == "cnet" || type == "cnetplot") {
##                   cnetplot.enrichResult(x, ...)
##               }
##               if (type == "bar" || type == "barplot") {
##                   barplot(x, ...)
##               }
##               if (type == "enrichMap") {
##                   enrichMap(x, ...)
##               }
##               if (type == "dot" || type == "dotplot") {
##                   dotplot(x, ...)
##               }
##           }
##           )


## ##' dotplot for enrichResult
## ##'
## ##'
## ##' @rdname dotplot-methods
## ##' @aliases dotplot,enrichResult,ANY-method
## ##' @param object an instance of enrichResult
## ##' @param x variable for x axis
## ##' @param colorBy one of 'pvalue', 'p.adjust' and 'qvalue'
## ##' @param showCategory number of category
## ##' @param split separate result by 'category' variable
## ##' @param font.size font size
## ##' @param title plot title
## ##' @exportMethod dotplot
## ##' @author Guangchuang Yu
## setMethod("dotplot", signature(object="enrichResult"),
##           function(object, x="geneRatio", colorBy="p.adjust", showCategory=10, split=NULL, font.size=12, title="") {
##               dotplot_internal(object, x, colorBy, showCategory, split, font.size, title)
##           }
##           )




## ##' @rdname cnetplot-methods
## ##' @exportMethod cnetplot
## setMethod("cnetplot", signature(x="enrichResult"),
##           function(x, showCategory=5, categorySize="pvalue", foldChange=NULL, fixed=TRUE, ...) {
##               cnetplot.enrichResult(x,
##                                     showCategory=showCategory,
##                                     categorySize=categorySize,
##                                     foldChange=foldChange,
##                                     fixed=fixed, ...)
##           }
##           )


