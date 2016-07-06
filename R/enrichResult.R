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
##' @slot keytype Gene ID type
##' @slot universe background gene
##' @slot geneInCategory gene and category association
##' @slot gene2Symbol mapping gene to Symbol
##' @slot geneSets gene sets
##' @slot readable logical flag of gene ID in symbol or not.
##' @exportClass enrichResult
##' @author Guangchuang Yu \url{https://guangchuangyu.github.io}
##' @seealso \code{\link{enrichDO}}
##' @keywords classes
setClass("enrichResult",
         representation=representation(
             result         = "data.frame",
             pvalueCutoff   = "numeric",
             pAdjustMethod  = "character",
             qvalueCutoff   = "numeric",
             organism       = "character",
             ontology       = "character",
             gene           = "character",
             keytype        = "character",
             universe       = "character",
             geneInCategory = "list",
             gene2Symbol    = "character",
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
##' @importFrom utils str
##' @importFrom methods show
##' @exportMethod show
##' @usage show(object)
##' @author Guangchuang Yu \url{https://guangchuangyu.github.io}
setMethod("show", signature(object="enrichResult"),
          function (object){
              cat("#\n# over-representation test\n#\n")
              cat("#...@organism", "\t", object@organism, "\n")
              cat("#...@ontology", "\t", object@ontology, "\n")
              kt <- object@keytype
              if (kt != "UNKNOWN") {
                  cat("#...@keytype", "\t", kt, "\n")
              }
              cat("#...@gene", "\t")
              str(object@gene)
              cat("#...pvalues adjusted by", paste0("'", object@pAdjustMethod, "'"),
                  paste0("with cutoff <", object@pvalueCutoff), "\n")
              cat(paste0("#...", nrow(object@result)), "enriched terms found\n")
              str(object@result)
              cat("#...Citation\n")
              if (object@ontology == "DO" || object@ontology == "DOLite" || object@ontology == "NCG") {
                  citation_msg <- paste("  Guangchuang Yu, Li-Gen Wang, Guang-Rong Yan, Qing-Yu He. DOSE: an",
                                    "  R/Bioconductor package for Disease Ontology Semantic and Enrichment",
                                    "  analysis. Bioinformatics 2015, 31(4):608-609", sep="\n", collapse="\n")
              } else if (object@ontology == "Reactome") {
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
          })


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
##' @author Guangchuang Yu \url{http://guangchuangyu.github.io}
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
##' @author Guangchuang Yu \url{http://guangchuangyu.github.io}
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


##' dotplot for enrichResult
##'
##' 
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


##' upsetplot
##'
##' 
##' @rdname upsetplot-methods
##' @aliases upsetplot,enrichResult,ANY-method
##' @param n number of categories to be plotted
##' @author Guangchuang Yu
##' @exportMethod upsetplot
##' @examples
##' \dontrun{
##' require(DOSE)
##' data(geneList)
##' de=names(geneList)[1:100]
##' x <- enrichDO(de)
##' upsetplot(x, 8)
##' }
setMethod("upsetplot", signature(x="enrichResult"),
          function(x, n=10, ...) {
              upsetplot.enrichResult(x, n, ...)
          })


##' mapping geneID to gene Symbol
##'
##'
##' @title setReadable
##' @param x enrichResult Object
##' @param OrgDb OrgDb
##' @param keytype keytype of gene
##' @return enrichResult Object
##' @author Yu Guangchuang
##' @export
setReadable <- function(x, OrgDb, keytype="auto") {
    if (!(is(x, "enrichResult") || is(x, "groupGOResult") || is(x, "gseaResult")))
        stop("input should be an 'enrichResult' or 'gseaResult' object...")
    
    isGSEA <- FALSE
    if (is(x, 'gseaResult')) 
        isGSEA <- TRUE
    
    if (keytype == "auto") {
        keytype <- x@keytype
    }
    
    if (x@readable)
        return(x)
    
    if (isGSEA) {
        gc <- x@core_enrichment
        genes <- names(x@geneList)
    } else {
        gc <- x@geneInCategory
        genes <- x@gene
    }
    
    gn <- EXTID2NAME(OrgDb, genes, keytype)        
    gc <- lapply(gc, function(i) gn[i])
    
    res <- x@result
    gc <- gc[as.character(res$ID)]
    geneID <- sapply(gc, paste0, collapse="/")
    if (isGSEA) {
        x@core_enrichment <- gc
        res$coreEnrichment <- unlist(geneID)
    } else {
        x@geneInCategory <- gc
        res$geneID <- unlist(geneID)
    }

    x@gene2Symbol <- gn
    x@result <- res
    x@keytype <- keytype
    x@readable <- TRUE
    
    return(x)
}

##' @rdname cnetplot-methods
##' @exportMethod cnetplot
setMethod("cnetplot", signature(x="enrichResult"),
          function(x, showCategory=5, categorySize="pvalue", foldChange=NULL, fixed=TRUE, ...) {
              cnetplot.enrichResult(x,
                                    showCategory=showCategory,
                                    categorySize=categorySize,
                                    foldChange=foldChange,
                                    fixed=fixed, ...)
          }
          )

##' accessing gene set
##'
##' 
##' @rdname subset2-methods
##' @title [[ method
##' @param x enrichResult object
##' @param i term to query
##' @exportMethod [[
setMethod("[[", signature(x="enrichResult"),
          function(x, i) {
              if (!i %in% names(x@geneInCategory))
                  stop("input term not found...")
              x@geneInCategory[[i]]
})


##' accessing enriched result
##'
##' 
##' @rdname subset-methods
##' @title [ method
##' @param x enrichResult object
##' @param i index
##' @exportMethod [
setMethod("[", signature(x="enrichResult"),
          function(x, i) {
              x@result[i,]
})


##' @importFrom utils head
##' @method head enrichResult
##' @export
head.enrichResult <- function(x, n=6L, ...) {
    head(x@result, n, ...)
}

##' @importFrom utils tail
##' @method tail enrichResult
##' @export
tail.enrichResult <- function(x, n=6L, ...) {
    tail(x@result, n, ...)
}

