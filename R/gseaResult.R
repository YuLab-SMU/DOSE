##' Class "gseaResult"
##' This class represents the result of GSEA analysis
##'
##'
##' @name gseaResult-class
##' @aliases gseahResult-class
##'   show,gseaResult-method summary,gseaResult-method
##'   plot,gseaResult-method
##'
##' @docType class
##' @slot result GSEA anaysis
##' @slot organism organism
##' @slot setType setType
##' @slot geneSets geneSets
##' @slot core_enrichment leading genes of enriched sets
##' @slot geneList order rank geneList
##' @slot keytype ID type of gene
##' @slot permScores permutation scores
##' @slot params parameters
##' @slot gene2Symbol gene ID to Symbol
##' @slot readable whether convert gene ID to symbol
##' @exportClass gseaResult
##' @author Guangchuang Yu \url{https://guangchuangyu.github.io}
##' @seealso \code{\link{gseaplot}}
##' @keywords classes
setClass("gseaResult",
         representation   = representation(
             result          = "data.frame",
             organism        = "character",
             setType         = "character", 
             geneSets        = "list",
             core_enrichment = "list",
             geneList        = "numeric",
             keytype         = "character",
             permScores      = "matrix",
             params          = "list",
             gene2Symbol     = "character",
             readable        = "logical"
         )
         )


##' @rdname cnetplot-methods
##' @exportMethod cnetplot
setMethod("cnetplot", signature(x="gseaResult"),
          function(x, showCategory=5, categorySize="pvalue", foldChange=NULL, fixed=TRUE, ...) {
              cnetplot.enrichResult(x,
                                    showCategory=showCategory,
                                    categorySize=categorySize,
                                    foldChange=foldChange,
                                    fixed=fixed, ...)
          }
          )


##' @rdname subset2-methods
##' @exportMethod [[
setMethod("[[", signature(x="gseaResult"),
          function(x, i) {
              if (!i %in% names(x@core_enrichment))
                  stop("input term not found...")
              x@core_enrichment[[i]]
          })


##' @rdname subset-methods
##' @exportMethod [
setMethod("[", signature(x="gseaResult"),
          function(x, i, j) {
              x@result[i,j]
})


##' @rdname subset3-methods
##' @exportMethod $
setMethod("$", signature(x="gseaResult"),
          function(x, name) {
              x@result[, name]
          })


##' @importFrom utils head
##' @method head gseaResult
##' @export
head.gseaResult <- function(x, n=6L, ...) {
    head(x@result, n, ...)
}

##' @importFrom utils tail
##' @method tail gseaResult
##' @export
tail.gseaResult <- function(x, n=6L, ...) {
    tail(x@result, n, ...)
}



##' show method for \code{gseaResult} instance
##'
##' @name show
##' @docType methods
##' @rdname show-methods
##' 
##' @title show method
##' @return message
##' @importFrom methods show
##' @exportMethod show
##' @usage show(object)
##' @author Guangchuang Yu \url{https://guangchuangyu.github.io}
setMethod("show", signature(object="gseaResult"),
          function (object){
              params <- object@params
              cat("#\n# Gene Set Enrichment Analysis\n#\n")
              cat("#...@organism", "\t", object@organism, "\n")
              cat("#...@setType", "\t", object@setType, "\n")
              kt <- object@keytype
              if (kt != "UNKNOWN") {
                  cat("#...@keytype", "\t", kt, "\n")
              }

              cat("#...@geneList", "\t")
              str(object@geneList)
              cat("#...nPerm", "\t", params$nPerm, "\n")
              cat("#...pvalues adjusted by", paste0("'", params$pAdjustMethod, "'"),
                  paste0("with cutoff <", params$pvalueCutoff), "\n")
              cat(paste0("#...", nrow(object@result)), "enriched terms found\n")
              str(object@result)
              cat("#...Citation\n")
              if (object@setType == "DO" || object@setType == "DOLite" || object@setType == "NCG") {
                  citation_msg <- paste("  Guangchuang Yu, Li-Gen Wang, Guang-Rong Yan, Qing-Yu He. DOSE: an",
                                    "  R/Bioconductor package for Disease Ontology Semantic and Enrichment",
                                    "  analysis. Bioinformatics 2015, 31(4):608-609", sep="\n", collapse="\n")
              } else if (object@setType == "Reactome") {
                  citation_msg <- paste("  Guangchuang Yu, Qing-Yu He. ReactomePA: an R/Bioconductor package for",
                                        "  reactome pathway analysis and visualization. Molecular BioSystems",
                                        "  2016, 12(2):477-479", sep="\n", collapse="\n")
              } else {
                  citation_msg <- paste("  Guangchuang Yu, Li-Gen Wang, Yanyan Han and Qing-Yu He.",
                                        "  clusterProfiler: an R package for comparing biological themes among",
                                        "  gene clusters. OMICS: A Journal of Integrative Biology",
                                        "  2012, 16(5):284-287", sep="\n", collapse="\n")
              }
              cat(citation_msg, "\n\n")
          }
          )


##' summary method for \code{gseaResult} instance
##'
##'
##' @name summary
##' @docType methods
##' @rdname summary-methods
##'
##' @title summary method
##' @return A data frame
##' @importFrom stats4 summary
##' @exportMethod summary
##' @usage summary(object, ...)
##' @author Guangchuang Yu \url{https://guangchuangyu.github.io}
setMethod("summary", signature(object="gseaResult"),
          function(object, ...) {
              return(object@result)
          }
          )

##' plot method for gseaResult
##'
##'
##' @docType methods
##' @name plot
##' @rdname plot-methods
##' @aliases plot,gseaResult,ANY-method
##' @title plot method
## @param type one of gseaplot or enrichMap
## @param ... additional parameter
##' @return plot
##' @importFrom stats4 plot
##' @exportMethod plot
##' @author Yu Guangchuang
setMethod("plot", signature(x="gseaResult"),
          function(x, type="gseaplot", ...) {
          ## function(x, geneSetID, by="all", ...) {
              if (type == "gseaplot") {
                  gseaplot(x, ...)
              }
              if (type == "enrichMap") {
                  enrichMap(x, ...)
              }
          }
          )

