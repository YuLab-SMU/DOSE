## ##' dotplot for gseaResult
## ##'
## ##'
## ##' @rdname dotplot-methods
## ##' @aliases dotplot,gseaResult,ANY-method
## ##' @exportMethod dotplot
## ##' @author Guangchuang Yu
## setMethod("dotplot", signature(object="gseaResult"),
##           function(object, x="geneRatio", colorBy="p.adjust", showCategory=10, split=NULL, font.size=12, title="") {
##               dotplot_internal(object, x, colorBy, showCategory, split, font.size, title)
##           }
##           )


## ##' @rdname cnetplot-methods
## ##' @exportMethod cnetplot
## setMethod("cnetplot", signature(x="gseaResult"),
##           function(x, showCategory=5, categorySize="pvalue", foldChange=NULL, fixed=TRUE, ...) {
##               cnetplot.enrichResult(x,
##                                     showCategory=showCategory,
##                                     categorySize=categorySize,
##                                     foldChange=foldChange,
##                                     fixed=fixed, ...)
##           }
##           )







## ##' plot method for gseaResult
## ##'
## ##'
## ##' @docType methods
## ##' @name plot
## ##' @rdname plot-methods
## ##' @aliases plot,gseaResult,ANY-method
## ##' @title plot method
## ## @param type one of gseaplot or enrichMap
## ## @param ... additional parameter
## ##' @return plot
## ##' @importFrom stats4 plot
## ##' @exportMethod plot
## ##' @author Yu Guangchuang
## setMethod("plot", signature(x="gseaResult"),
##           function(x, type="gseaplot", ...) {
##           ## function(x, geneSetID, by="all", ...) {
##               if (type == "gseaplot") {
##                   gseaplot(x, ...)
##               }
##               if (type == "enrichMap") {
##                   enrichMap(x, ...)
##               }
##           }
##           )

