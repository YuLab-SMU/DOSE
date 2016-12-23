##' dotplot for gseaResult
##'
##'
##' @rdname dotplot-methods
##' @aliases dotplot,gseaResult,ANY-method
##' @exportMethod dotplot
##' @author Guangchuang Yu
setMethod("dotplot", signature(object="gseaResult"),
          function(object, x="geneRatio", colorBy="p.adjust", showCategory=10, split=NULL, font.size=12, title="") {
              dotplot_internal(object, x, colorBy, showCategory, split, font.size, title)
          }
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
              warning("summary method to convert the object to data.frame is deprecated, please use as.data.frame instead.")
              return(as.data.frame(object, ...))
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

