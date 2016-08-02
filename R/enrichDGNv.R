##' Enrichment analysis based on the DisGeNET (\url{http://www.disgenet.org/})
##'
##' given a vector of genes, this function will return the enrichment NCG
##' categories with FDR control
##'
##'
##' @title enrichDGN
##' @param snp a vector of SNP
##' @param pvalueCutoff pvalue cutoff
##' @param pAdjustMethod one of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"
##' @param universe background genes
##' @param minGSSize minimal size of genes annotated by NCG category for testing
##' @param maxGSSize maximal size of each geneSet for analyzing
##' @param qvalueCutoff qvalue cutoff
##' @param readable whether mapping gene ID to gene Name
##' @return A \code{enrichResult} instance
##' @export
##' @references Janet et al. (2015) DisGeNET: a discovery platform for the dynamical exploration of human diseases and their genes. \emph{Database} bav028
##' \url{http://database.oxfordjournals.org/content/2015/bav028.long}
##' @author Guangchuang Yu
enrichDGNv <- function(snp,
                      pvalueCutoff = 0.05,
                      pAdjustMethod = "BH",
                      universe,
                      minGSSize = 10,
                      maxGSSize = 500,
                      qvalueCutoff = 0.2,
                      readable = FALSE){
    enrichDisease(gene = snp,
                  pvalueCutoff = pvalueCutoff,
                  pAdjustMethod = pAdjustMethod,
                  universe = universe,
                  minGSSize = minGSSize,
                  maxGSSize = maxGSSize,
                  qvalueCutoff = qvalueCutoff,
                  readable = readable,
                  ontology = "snpDisGeNET")
 
}

get_VDGN_data <- function() {
    if (!exists(".VDGN_DOSE_Env")) {
        tryCatch(utils::data(list="VDGN_DOSE_Env", package="DOSE"))
    }
    get(".VDGN_DOSE_Env", envir = .GlobalEnv)
}


