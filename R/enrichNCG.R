##' Enrichment analysis based on the Network of Cancer Genes database (http://ncg.kcl.ac.uk/)
##'
##' given a vector of genes, this function will return the enrichment NCG
##' categories with FDR control
##'
##' 
##' @title enrichNCG
##' @param gene a vector of entrez gene id
##' @param pvalueCutoff pvalue cutoff
##' @param pAdjustMethod one of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"
##' @param universe background genes
##' @param minGSSize minimal size of genes annotated by NCG category for testing
##' @param qvalueCutoff qvalue cutoff
##' @param readable whether mapping gene ID to gene Name
##' @return A \code{enrichResult} instance
##' @export
##' @author Guangchuang Yu
enrichNCG <- function(gene,
                      pvalueCutoff = 0.05,
                      pAdjustMethod = "BH",
                      universe,
                      minGSSize = 5,
                      qvalueCutoff = 0.2,
                      readable = FALSE){
    
    NCG_DOSE_Env <- get_NCG_data()
    res <- enricher_internal(gene = gene,
                             pvalueCutoff = pvalueCutoff,
                             pAdjustMethod = pAdjustMethod,
                             universe = universe,
                             minGSSize = minGSSize,
                             qvalueCutoff = qvalueCutoff,
                             USER_DATA = NCG_DOSE_Env)
    res@organism <- "Homo sapiens"
    res@keytype <- "ENTREZID"
    res@ontology <- "NCG"

    if(readable) {
        res <- setReadable(res, 'org.Hs.eg.db')
    }
    return(res)
}

get_NCG_data <- function() {
    if (!exists("NCG_DOSE_Env")) {
        tryCatch(utils::data(list="NCG_DOSE_Env", package="DOSE"))
    }
    get("NCG_DOSE_Env", envir = .GlobalEnv)
}


