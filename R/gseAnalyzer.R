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
#' @param geneList order ranked geneList
#' @param ont one of "HDO", "HPO" or "MPO"
#' @param organism one of "hsa" and "mmu"
#' @param exponent weight of each step
#' @param nPerm permutation numbers
#' @param minGSSize minimal size of each geneSet for analyzing
#' @param maxGSSize maximal size of each geneSet for analyzing
#' @param pvalueCutoff pvalue Cutoff
#' @param pAdjustMethod p value adjustment method
#' @param verbose print message or not
#' @param adaptive logical, use adaptive permutation or not (default: FALSE)
#' @param minPerm minimum number of permutations for adaptive mode (default: 1000)
#' @param maxPerm maximum number of permutations for adaptive mode (default: 10000)
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
               adaptive          = adaptive,
               minPerm           = minPerm,
               maxPerm           = maxPerm,
               ...)

}

#' NCG Gene Set Enrichment Analysis
#'
#'
#' perform gsea analysis
#' @inheritParams gseDO
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
               adaptive          = adaptive,
               minPerm           = minPerm,
               maxPerm           = maxPerm,
               ...)
    


}

#' DisGeNET Gene Set Enrichment Analysis
#'
#'
#' perform gsea analysis
#' @inheritParams gseDO
#' @return gseaResult object
#' @export
#' @author Guangchuang Yu
#' @keywords manip
gseDGN <- function(geneList,
                   exponent=1,
                   nPerm = 1000,
                   minGSSize = 10,
                   maxGSSize = 500,
                   pvalueCutoff=0.05,
                   pAdjustMethod="BH",
                   verbose=TRUE,
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
               ontology          = "DisGeNET",
               adaptive          = adaptive,
               minPerm           = minPerm,
               maxPerm           = maxPerm,
               ...)
}
