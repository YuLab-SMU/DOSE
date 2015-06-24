##' Class "enrichResult"
##' This class represents the result of enrichment analysis.
##'
##'
##' @name enrichResult-class
##' @aliases enrichResult-class
##'   show,enrichResult-method summary,enrichResult-method
##'   plot,enrichResult-method
##'
##' @docType class
##' @slot result enrichment analysis
##' @slot pvalueCutoff pvalueCutoff
##' @slot pAdjustMethod pvalue adjust method
##' @slot qvalueCutoff qvalueCutoff
##' @slot organism only "human" supported
##' @slot ontology biological ontology
##' @slot gene Gene IDs
##' @slot universe background gene
##' @slot geneInCategory gene and category association
##' @slot geneSets gene sets
##' @slot readable logical flag of gene ID in symbol or not.
##' @exportClass enrichResult
##' @author Guangchuang Yu \url{http://ygc.name}
##' @seealso \code{\link{enrichDO}}
##' @keywords classes
setClass("enrichResult",
         representation=representation(
             result         ="data.frame",
             pvalueCutoff   ="numeric",
             pAdjustMethod  ="character",
             qvalueCutoff   ="numeric",
             organism       = "character",
             ontology       = "character",
             gene           = "character",
             universe       = "character",
             geneInCategory = "list",
             geneSets       = "list",
             readable       = "logical"
             ),
         prototype=prototype(readable = FALSE)
         )

##' show method for \code{enrichResult} instance
##'
##' @name show
##' @docType methods
##' @rdname show-methods
##' 
##' @title show method
##' @param object A \code{enrichResult} instance.
##' @return message
##' @importFrom methods show
##' @exportMethod show
##' @usage show(object)
##' @author Guangchuang Yu \url{http://ygc.name}
setMethod("show", signature(object="enrichResult"),
          function (object){
              organism = object@organism
              ontology = object@ontology
              geneNum = length(object@gene)
              pvalueCutoff=object@pvalueCutoff
              cat (geneNum, organism, "Genes to ", ontology,
                   " test for over-representation.", "\n",
                   "with p value <", pvalueCutoff, "\n")
          }
          )

##' summary method for \code{enrichResult} instance
##'
##'
##' @name summary
##' @docType methods
##' @rdname summary-methods
##'
##' @title summary method
##' @param object A \code{enrichResult} instance.
##' @param ... additional parameter
##' @return A data frame
##' @importFrom stats4 summary
##' @exportMethod summary
##' @usage summary(object, ...)
##' @author Guangchuang Yu \url{http://ygc.name}
setMethod("summary", signature(object="enrichResult"),
          function(object, ...) {
              return(object@result)
          }
          )

##' plot method generics
##'
##'
##' @docType methods
##' @name plot
##' @rdname plot-methods
##' @aliases plot,enrichResult,ANY-method
##' @title plot method
##' @param x A \code{enrichResult} instance
##' @param type one of bar, cnet or enrichMap
##' @param ... Additional argument list
##' @return plot
##' @importFrom stats4 plot
##' @exportMethod plot
##' @author Guangchuang Yu \url{http://ygc.name}
setMethod("plot", signature(x="enrichResult"),
          function(x, type = "bar", ... ) {
              if (type == "cnet" || type == "cnetplot") {
                  cnetplot.enrichResult(x, ...)
              }
              if (type == "bar" || type == "barplot") {
                  barplot(x, ...)
              }
              if (type == "enrichMap") {
                  enrichMap(x, ...)
              }
              if (type == "dot" || type == "dotplot") {
                  dotplot(x, ...)
              }
          }
          )


##' @rdname dotplot-methods
##' @aliases dotplot,enrichResult,ANY-method
##' @param object an instance of enrichResult
##' @param x variable for x axis
##' @param colorBy one of 'pvalue', 'p.adjust' and 'qvalue'
##' @param showCategory number of category
##' @param font.size font size
##' @param title plot title
##' @exportMethod dotplot
##' @author Guangchuang Yu
setMethod("dotplot", signature(object="enrichResult"),
          function(object, x="geneRatio", colorBy="p.adjust", showCategory=10, font.size=12, title="") {
              dotplot.enrichResult(object, x, colorBy, showCategory, font.size, title)
          }
          )



##' mapping geneID to gene Symbol
##'
##'
##' @title setReadable
##' @param x enrichResult Object
##' @return enrichResult Object
##' @author Yu Guangchuang
##' @export
setReadable <- function(x) {
    if (!(class(x) != "enrichResult" || class(x) != "groupGOResult"))
        stop("input should be an 'enrichResult' object...")
    if (x@readable == FALSE) {
        organism = x@organism
        gc <- x@geneInCategory
        genes <- x@gene
        gn <- EXTID2NAME(genes, organism=organism)
        ##gc <- lapply(gc, EXTID2NAME, organism=organism)
        gc <- lapply(gc, function(i) gn[i])
        x@geneInCategory <- gc

        res <- x@result
        gc <- gc[as.character(res$ID)]
        geneID <- sapply(gc, function(i) paste(i, collapse="/"))
        res$geneID <- unlist(geneID)

        x@result <- res

        x@readable <- TRUE
    }
    return(x)
}

##' @rdname cnetplot-methods
##' @exportMethod cnetplot
setMethod("cnetplot", signature(x="enrichResult"),
          function(x, showCategory=5, categorySize="geneNum", foldChange=NULL, fixed=TRUE, ...) {
              cnetplot.enrichResult(x,
                                    showCategory=showCategory,
                                    categorySize=categorySize,
                                    foldChange=foldChange,
                                    fixed=fixed, ...)
          }
          )
