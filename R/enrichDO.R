##' Class "enrichResult"
##' This class represents the result of DO enrichment analysis.
##'
##'
##' @name enrichResult-class
##' @aliases enrichResult-class
##'   show,enrichResult-method summary,enrichResult-method
##'   plot,enrichResult-method setReadable<-,enrichResult-method
##'
##' @docType class
##' @slot result DO enrichment result
##' @slot pvalueCutoff pvalueCutoff
##' @slot qvalueCutoff qvalueCutoff
##' @slot organism only "human" supported
##' @slot ontology biological ontology
##' @slot gene Gene IDs
##' @slot geneInCategory gene and category association
##' @slot readable logical flag of gene ID in symbol or not.
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
         ontology = "character",
         gene = "character",
         geneInCategory = "list",
         readable = "logical"
         ),
         prototype=prototype(readable = FALSE)
         )


##' DO Enrichment Analysis of a gene set.
##'
##' Given a vector of genes, this function will return the enrichment DO
##' categories with FDR control.
##'
##'
##' @param gene a vector of entrez gene id.
##' @param ont one of DO or DOLite.
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
enrichDO <- function(gene, ont="DOLite",
                     pvalueCutoff=0.05, qvalueCutoff=1,
                     readable=F) {

    enrich.internal(gene,
                    organism = "human",
                    pvalueCutoff,
                    qvalueCutoff,
                    ont = ont,
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
##' @exportMethod show
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
##' @return A data frame
##' @importFrom stats4 summary
##' @exportMethod summary
##' @author Guangchuang Yu \url{http://ygc.name}
setMethod("summary", signature(object="enrichResult"),
          function(object) {
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
##' @param ... Additional argument list
##' @return plot
##' @importFrom stats4 plot
##' @exportMethod plot
##' @author Guangchuang Yu \url{http://ygc.name}
setMethod("plot", signature(x="enrichResult"),
          function(x, type = "cnet", ... ) {
              if (type == "cnet") {
                  cnetplot.enrichResult(x, ...)
              }
              if (type == "bar") {
                  barplot.enrichResult(x, ...)
              }
          }
          )

##' setReadable method for \code{enrichResult} instance
##'
##'
##' @name setReadable<-
##' @docType methods
##' @rdname setReadable-methods
##' @aliases setReadable
##' @aliases setReadable-methods
##' @aliases setReadable<-,enrichResult,ANY-method
##'
##' @title Methods mapping gene ID to gene symbol for \code{enrichResult} instance
##' @param x A \code{enrichResult} instance.
##' @param value readable flag.
##' @return A \code{enrichResult} instance.
##' @author Guangchuang Yu \url{http://ygc.name}
setReplaceMethod(
                 f= "setReadable",
                 signature= "enrichResult",
                 definition=function(x, value){
                     if (value == FALSE)
                         return(x)

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
                         res$geneID <- geneID

                         x@result <- res

                         x@readable <- value
                         return(x)
                     }
                 }
                 )

