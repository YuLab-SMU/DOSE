#' @importFrom enrichit gsea_gson
gseDisease <- function(geneList,
                       organism = "hsa",
                       exponent=1,
                       nPerm = 1000,
                       minGSSize = 10,
                       maxGSSize = 500,
                       pvalueCutoff=0.05,
                       pAdjustMethod="BH",
                       verbose=TRUE,
                       ontology,
                       method = "multilevel",
                       adaptive = FALSE,
                       minPerm = 1000,
                       maxPerm = 10000,
                       ...) {

    annoData <- get_anno_data(ontology)

    res <- gsea_gson(geneList          = geneList,
                         exponent          = exponent,
                         nPerm             = nPerm,
                         minGSSize         = minGSSize,
                         maxGSSize         = maxGSSize,
                         pvalueCutoff      = pvalueCutoff,
                         pAdjustMethod     = pAdjustMethod,
                         verbose           = verbose,
                         gson              = annoData,
                         method            = method,
                         adaptive          = adaptive,
                         minPerm           = minPerm,
                         maxPerm           = maxPerm,
                         ...)

    if (is.null(res))
        return(res)

    if (organism == "hsa") {
        res@organism <- "Homo sapiens"
    } else {
        res@organism <- "Mus musculus"
    }
    res@setType <- ontology
    res@keytype <- "ENTREZID"
    return(res)
}

#' DO Gene Set Enrichment Analysis
#'
#'
#' perform gsea analysis
#' @inheritParams dose_params
#' @param ... other parameter
#' @return gseaResult object
#' @export
#' @author Guangchuang Yu
#' @keywords manip
gseDO <- function(geneList,
                  ont = "HDO",
                  organism = "hsa",
                  exponent=1,
                  nPerm = 1000,
                  minGSSize = 10,
                  maxGSSize = 500,
                  pvalueCutoff=0.05,
                  pAdjustMethod="BH",
                  verbose=TRUE,
                  method = "multilevel",
                  adaptive = FALSE,
                  minPerm = 1000,
                  maxPerm = 10000,
                  ...) {
     

    gseDisease(geneList          = geneList,
               exponent          = exponent,
               nPerm             = nPerm,
               minGSSize         = minGSSize,
               maxGSSize         = maxGSSize,
               pvalueCutoff      = pvalueCutoff,
               pAdjustMethod     = pAdjustMethod,
               verbose           = verbose,
               ontology          = ont,
               method            = method,
               adaptive          = adaptive,
               minPerm           = minPerm,
               maxPerm           = maxPerm,
               ...)

}

#' NCG Gene Set Enrichment Analysis
#'
#'
#' perform gsea analysis
#' @inheritParams dose_params
#' @param ... other parameter
#' @return gseaResult object
#' @export
#' @author Guangchuang Yu
#' @keywords manip
gseNCG <- function(geneList,
                   exponent=1,
                   nPerm = 1000,
                   minGSSize = 10,
                   maxGSSize = 500,
                   pvalueCutoff=0.05,
                   pAdjustMethod="BH",
                   verbose=TRUE,
                   method = "multilevel",
                   adaptive = FALSE,
                   minPerm = 1000,
                   maxPerm = 10000,
                   ...) {
                  

    gseDisease(geneList          = geneList,
               exponent          = exponent,
               nPerm             = nPerm,
               minGSSize         = minGSSize,
               maxGSSize         = maxGSSize,
               pvalueCutoff      = pvalueCutoff,
               pAdjustMethod     = pAdjustMethod,
               verbose           = verbose,
               ontology          = "NCG",
               method            = method,
               adaptive          = adaptive,
               minPerm           = minPerm,
               maxPerm           = maxPerm,
               ...)
    


}


