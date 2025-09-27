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
##' @author Guangchuang Yu \url{https://yulab-smu.top}
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
              print_citation_msg(object@setType)
          }
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
##' @author Guangchuang Yu \url{https://yulab-smu.top}
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
              object <- get_enriched(object)
              n <- nrow(object@result)
              cat(paste0("#...", n), "enriched terms found\n")
              if (n > 0) str(object@result)
              cat("#...Citation\n")
              print_citation_msg(object@ontology)
        }
)


print_citation_msg <- function(ontology) {
    refs <- yulab.utils:::ref_knownledge()

    if (ontology == "HDO" || ontology == "NCG") {
        citation_msg <- refs["DOSE"]
    } else if (ontology == "Reactome") {
        citation_msg <- refs["ReactomePA"]
    } else if (ontology == "MeSH") {
        citation_msg <- refs["meshes"]
    } else {
        citation_msg <- refs["clusterProfiler_NP"]
    }
    cat(citation_msg, "\n\n")
}

