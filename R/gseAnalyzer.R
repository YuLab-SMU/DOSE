##' Class "gseaResult"
##' This class represents the result of GSEA analysis
##'
##'
##' @name gseaResult-class
##' @aliases gseahResult-class
##'   show,gseaResult-method summary,gseaResult-method
##'   plot,gseaResult-method [[,gseaResult-method
##'
##' @docType class
##' @slot result GSEA anaysis
##' @slot organism organism
##' @slot setType setType
##' @slot geneSets geneSets
##' @slot geneList order rank geneList
##' @slot keytype ID type of gene
##' @slot permScores permutation scores
##' @slot params parameters
##' @exportClass gseaResult
##' @author Guangchuang Yu \url{http://guangchuangyu.github.io}
##' @seealso \code{\link{gseaplot}}
##' @keywords classes
setClass("gseaResult",
         representation = representation(
             result     = "data.frame",
             organism   = "character",
             setType    = "character", 
             geneSets   = "list",
             geneList   = "numeric",
             keytype    = "character",
             permScores = "matrix",
             params     = "list"
         )
         )


##' accessing gene set
##'
##' 
##' @rdname subset-methods
##' @title [[ method
##' @exportMethod [[
setMethod("[[", signature(x="gseaResult"),
          function(x, term) {
              if (!term %in% names(x@geneSets))
                  stop("input term not found...")
              x@geneSets[[term]]
          })


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
##' @author Guangchuang Yu \url{http://guangchuangyu.github.io}
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
                                    "  analysis. Bioinformatics 2015 31(4):608-609", sep="\n", collapse="\n")
              } else if (object@setType == "Reactome") {
                  citation_msg <- paste("  Guangchuang Yu, Qing-Yu He. ReactomePA: an R/Bioconductor package for",
                                        "  reactome pathway analysis and visualization. Molecular BioSystems",
                                        "  2015 accepted", sep="\n", collapse="\n")
              } else {
                  citation_msg <- paste("  Guangchuang Yu, Li-Gen Wang, Yanyan Han and Qing-Yu He.",
                                        "  clusterProfiler: an R package for comparing biological themes among",
                                        "  gene clusters. OMICS: A Journal of Integrative Biology 2012,",
                                        "  16(5):284-287", sep="\n", collapse="\n")
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
##' @author Guangchuang Yu \url{http://guangchuangyu.github.io}
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

##' DO Gene Set Enrichment Analysis
##'
##'
##' perform gsea analysis
##' @param geneList order ranked geneList
##' @param exponent weight of each step
##' @param nPerm permutation numbers
##' @param minGSSize minimal size of each geneSet for analyzing
##' @param maxGSSize maximal size of each geneSet for analyzing
##' @param pvalueCutoff pvalue Cutoff
##' @param pAdjustMethod p value adjustment method
##' @param verbose print message or not
##' @param seed logical
##' @return gseaResult object
##' @export
##' @author Yu Guangchuang
##' @keywords manip
gseDO <- function(geneList,
                  exponent=1,
                  nPerm=1000,
                  minGSSize = 10,
                  maxGSSize = 500,
                  pvalueCutoff=0.05,
                  pAdjustMethod="BH",
                  verbose=TRUE,
                  seed=FALSE) {
    
    res <- GSEA_internal(geneList          = geneList,
                         exponent          = exponent,
                         nPerm             = nPerm,
                         minGSSize         = minGSSize,
                         maxGSSize         = maxGSSize,
                         pvalueCutoff      = pvalueCutoff,
                         pAdjustMethod     = pAdjustMethod,
                         verbose           = verbose,
                         seed              = seed,
                         USER_DATA         = get_DO_data())

    if (is.null(res))
        return(res)
    
    res@organism <- "Homo sapiens"
    res@setType <- "DO"
    res@keytype <- "ENTREZID"
    return(res)
}

##' NCG Gene Set Enrichment Analysis
##'
##'
##' perform gsea analysis
##' @inheritParams gseDO
##' @return gseaResult object
##' @export
##' @author Yu Guangchuang
##' @keywords manip
gseNCG <- function(geneList,
                  exponent=1,
                  nPerm=1000,
                  minGSSize = 10,
                  maxGSSize = 500,
                  pvalueCutoff=0.05,
                  pAdjustMethod="BH",
                  verbose=TRUE,
                  seed=FALSE) {
        
    res <- GSEA_internal(geneList          = geneList,
                         exponent          = exponent,
                         nPerm             = nPerm,
                         minGSSize         = minGSSize,
                         maxGSSize         = 500,
                         pvalueCutoff      = pvalueCutoff,
                         pAdjustMethod     = pAdjustMethod,
                         verbose           = verbose,
                         seed              = seed,
                         USER_DATA         = get_NCG_data())

    if (is.null(res))
        return(res)

    res@organism <- "Homo sapiens"
    res@setType <- "NCG"
    res@keytype <- "ENTREZID"
    return(res)
}


##' convert gsea result for ggplot2
##'
##' 
##' @importFrom ggplot2 fortify
## @S3method fortify gseaResult
##' @title fortify.gseaResult
##' @param model gseaResult object
##' @param data not used.
##' @param geneSetID gene set ID
##' @param ... additional parameter
##' @return figure
##' @author G Yu
##' @method fortify gseaResult
##' @export
fortify.gseaResult <- function(model, data, geneSetID, ...) {
    object <- model ## gseaResult object
    geneList <- object@geneList

    if (is.numeric(geneSetID))
        geneSetID <- object@result[geneSetID, "ID"]

    geneSet <- object@geneSets[[geneSetID]]
    exponent <- object@params[["exponent"]]
    df <- gseaScores(geneList, geneSet, exponent, fortify=TRUE)
    df$ymin=0
    df$ymax=0
    pos <- df$position == 1
    h <- diff(range(df$runningScore))/20
    df$ymin[pos] <- -h
    df$ymax[pos] <- h
    df$geneList <- geneList

    return(df)
}


##' visualize analyzing result of GSEA
##'
##' plotting function for gseaResult
##' @importFrom ggplot2 ggplot
##' @importFrom ggplot2 geom_linerange
##' @importFrom ggplot2 geom_line
##' @importFrom ggplot2 geom_vline
##' @importFrom ggplot2 geom_hline
##' @importFrom ggplot2 xlab
##' @importFrom ggplot2 ylab
##' @importFrom ggplot2 aes
##' @importFrom ggplot2 ggplotGrob
##' @importFrom grid grid.newpage
##' @importFrom grid viewport
##' @importFrom grid grid.layout
##' @importFrom grid pushViewport
##' @param gseaResult gseaResult object
##' @param geneSetID geneSet ID
##' @param by one of "runningScore" or "position"
##' @return ggplot2 object
##' @export
##' @author Yu Guangchuang
gseaplot <- function(gseaResult, geneSetID, by="all") {
    by <- match.arg(by, c("runningScore", "position", "all"))

    ## to satisfy codetools
    x <- ymin <- ymax <- runningScore <- es <- pos <- geneList <- NULL
    p <- ggplot(gseaResult,geneSetID=geneSetID,
                aes(x=x, ymin=ymin, ymax=ymax)) +
                    theme_dose() +
                        xlab("Position in the Ranked List of Genes")

    if (by == "runningScore" || by == "all") {
        p.res <- p+geom_linerange(colour="#DAB546")
        p.res <- p.res + geom_line(aes(y=runningScore))

        enrichmentScore <- gseaResult@result[geneSetID, "enrichmentScore"]
        es.df <- data.frame(es = which(p$data$runningScore == enrichmentScore))
        p.res <- p.res + geom_vline(data=es.df, aes(xintercept=es),
                            colour="#FA5860", linetype="dashed")
        p.res <- p.res + ylab("Running Enrichment Score")
        p.res <- p.res + geom_hline(aes(yintercept=0))
    }

    if (by == "position" || by == "all" ) {
        df2 <- data.frame(pos=which(p$data$position==1))
        p.pos <- p + geom_vline(data=df2, aes(xintercept=pos),
                            colour="#DAB546", alpha=.3)
        p.pos <- p.pos + geom_line(aes(y=geneList), colour="red")
        p.pos <- p.pos + ylab("Phenotype")
        p.pos <- p.pos + geom_hline(aes(yintercept=0))
    }

    if (by == "runningScore")
        return (p.res)
    if (by == "position")
        return (p.pos)

    p.pos <- p.pos + xlab("") +
        theme(axis.text.x = element_blank(),
              axis.ticks.x= element_blank())
    p.res <- p.res + theme(axis.title.x=element_text(vjust=-.3))
    ## two plots in one page
    grid.newpage()
    pushViewport(viewport(layout=grid.layout(2,1, heights=c(.3, .7))))
    print(p.pos, vp=viewport(layout.pos.row=1, layout.pos.col=1))
    print(p.res, vp=viewport(layout.pos.row=2, layout.pos.col=1))
    invisible(list(runningScore=p.res, position=p.pos))
}

