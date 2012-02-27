##' Class "enrichResult"
##' This class represents the result of DO enrichment analysis.
##'
##'
##' @name enrichResult-class
##' @aliases enrichResult-class
##'   show,enrichResult-method summary,enrichResult-method plot,enrichResult-method
##'
##' @docType class
##' @slot result DO enrichment result
##' @slot pvalueCutoff pvalueCutoff
##' @slot qvalueCutoff qvalueCutoff
##' @slot organism only "human" supported
##' @slot gene Gene IDs
##' @slot geneInCategory gene and category association
##' @exportClass enrichResult
##' @author Guangchuang Yu \url{http://ygc.name}
##' @seealso \code{\link{enrichDO}}
##' @keywords classes
setClass("enrichResult",
         representation=representation(
         result="data.frame",
         pvalueCutoff="numeric",
         qvalueCutoff="numeric",
         organism = "character",
         gene = "character",
         geneInCategory = "list"
         )
         )


##' DO Enrichment Analysis of a gene set.
##'
##' Given a vector of genes, this function will return the enrichment DO
##' categories with FDR control.
##'
##'
##' @param gene a vector of entrez gene id.
##' @param pvalueCutoff Cutoff value of pvalue.
##' @param qvalueCutoff Cutoff value of qvalue.
##' @param readable whether mapping gene ID to gene Name
##' @return A \code{enrichResult} instance.
##' @export
##' @seealso \code{\link{enrichResult-class}}
##' @author Guangchuang Yu \url{http://ygc.name}
##' @keywords manip
##' @examples
##'
##' 	set.seed(123)
##' 	data(EG2DO)
##' 	gene = sample(names(EG2DO), 30)
##' 	yy = enrichDO(gene, pvalueCutoff=0.05)
##' 	summary(yy)
##'
enrichDO <- function(gene, pvalueCutoff=0.05, qvalueCutoff=0.05, readable=F) {
    enrich.internal(gene,
                    organism = "human",
                    pvalueCutoff,
                    qvalueCutoff,
                    ont = "DO",
                    readable)
}

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
##' @author Guangchuang Yu \url{http://ygc.name}
setMethod("show", signature(object="enrichResult"),
	function (object){
		organism = object@organism
		geneNum = length(object@gene)
		pvalueCutoff=object@pvalueCutoff
		cat (geneNum, organism, "Genes to DO test for over-representation.", "\n", "p value <", pvalueCutoff, "\n")
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
##' @return A data frame
##' @importFrom BiocGenerics summary
##' @exportMethod summary
##' @author Guangchuang Yu \url{http://ygc.name}
setMethod("summary", signature(object="enrichResult"),
	function(object) {
		return(object@result)
	}
)

##' @rdname plot-methods
##' @aliases plot,enrichResult,ANY-method
setMethod("plot", signature(x="enrichResult"),
          function(x, type = "cnet", ... ) {
              if (type == "cnet") {
                  cnetplot.enrichResult(x, ...)
              }
          }
          )


