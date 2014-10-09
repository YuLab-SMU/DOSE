##' Class "enrichResult"
##' This class represents the result of enrichment analysis.
##'
##'
##' @name enrichResult-class
##' @aliases enrichResult-class
##'   show,enrichResult-method summary,enrichResult-method
##'   plot,enrichResult-method cnetplot,enrichResult-method
##'
##' @docType class
##' @slot result enrichment analysis
##' @slot pvalueCutoff pvalueCutoff
##' @slot pAdjustMethod pvalue adjust method
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
         pAdjustMethod="character",
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
##' @param pAdjustMethod one of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"
##' @param universe background genes
##' @param minGSSize minimal size of genes annotated by Ontology term for testing.
##' @param qvalueCutoff qvalue Cutoff
##' @param readable whether mapping gene ID to gene Name
##' @return A \code{enrichResult} instance.
##' @export
##' @seealso \code{\link{enrichResult-class}}
##' @author Guangchuang Yu \url{http://ygc.name}
##' @keywords manip
##' @examples
##'
##'	data(geneList)
##' 	gene = names(geneList)[geneList > 1]
##' 	yy = enrichDO(gene, pvalueCutoff=0.05)
##' 	summary(yy)
##'
enrichDO <- function(gene, ont="DOLite",
                     pvalueCutoff=0.05,
                     pAdjustMethod="BH",
                     universe,
                     minGSSize = 5,
                     qvalueCutoff=0.2,
                     readable=FALSE) {

    enrich.internal(gene,
                    organism = "human",
                    pvalueCutoff=pvalueCutoff,
                    pAdjustMethod=pAdjustMethod,
                    ont = ont,
                    universe = universe,
                    minGSSize = minGSSize,
                    qvalueCutoff = qvalueCutoff,
                    readable = readable)
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
              if (type == "cnet") {
                  cnetplot.enrichResult(x, ...)
              }
              if (type == "bar") {
                  barplot(x, ...)
              }
              if (type == "enrichMap") {
                  enrichMap(x, ...)
              }
          }
          )

##' cnetplot method generics
##'
##'
##' @docType methods
##' @name cnetplot
##' @rdname cnetplot-methods
##' @aliases cnetplot,enrichResult,ANY-method
##' @title cnetplot method
##' @param x enrichResult object
##' @param showCategory number of category plotted
##' @param categorySize one of geneNum or pvalue
##' @param foldChange fold change of expression value
##' @param fixed logical
##' @param ... additional parameter
##' @return plot
##' @exportMethod cnetplot
##' @usage cnetplot(x, showCategory=5, categorySize="geneNum", foldChange=NULL, fixed=TRUE, ...)
##' @author Guangchuang Yu \url{http://ygc.name}
setMethod("cnetplot", signature(x="enrichResult"),
          function(x, showCategory=5, categorySize="geneNum", foldChange=NULL, fixed=TRUE, ...) {
              cnetplot.enrichResult(x,
                                    showCategory=showCategory,
                                    categorySize=categorySize,
                                    foldChange=foldChange,
                                    fixed=fixed, ...)
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

