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
##' @param maxGSSize maximal size of each geneSet for analyzing
##' @param qvalueCutoff qvalue cutoff
##' @param readable whether mapping gene ID to gene Name
##' @return A \code{enrichResult} instance
##' @export
##' @author Guangchuang Yu
enrichNCG <- function(gene,
                      pvalueCutoff = 0.05,
                      pAdjustMethod = "BH",
                      universe,
                      minGSSize = 10,
                      maxGSSize = 500,
                      qvalueCutoff = 0.2,
                      readable = FALSE){
    
    enrichDisease(gene = gene,
                  pvalueCutoff = pvalueCutoff,
                  pAdjustMethod = pAdjustMethod,
                  universe = universe,
                  minGSSize = minGSSize,
                  maxGSSize = maxGSSize,
                  qvalueCutoff = qvalueCutoff,
                  readable = readable,
                  ontology = "NCG")
}

get_NCG_data <- function() {
    if (!exists(".DOSEenv")) .initial()
    .DOSEEnv <- get(".DOSEEnv", envir = .GlobalEnv)
    
    if (!exists(".NCG_DOSE_Env", envir=.DOSEEnv)) {
        tryCatch(utils::data(list="NCG_EXTID2PATHID", package="DOSE"))
        tryCatch(utils::data(list="NCG_PATHID2EXTID", package="DOSE"))
        tryCatch(utils::data(list="NCG_PATHID2NAME", package="DOSE"))
        EXTID2PATHID <- NCG_EXTID2PATHID <- get("NCG_EXTID2PATHID")
        PATHID2EXTID <- NCG_PATHID2EXTID <- get("NCG_PATHID2EXTID")
        PATHID2NAME <- NCG_PATHID2NAME <- get("NCG_PATHID2NAME")

        rm(NCG_EXTID2PATHID, envir = .GlobalEnv)
        rm(NCG_PATHID2EXTID, envir = .GlobalEnv)
        rm(NCG_PATHID2NAME, envir = .GlobalEnv)

        assign(".NCG_DOSE_Env", new.env(), envir = .DOSEEnv)
        .NCG_DOSE_Env <- get(".NCG_DOSE_Env", envir = .DOSEEnv)
        assign("EXTID2PATHID", EXTID2PATHID, envir = .NCG_DOSE_Env)
        assign("PATHID2EXTID", PATHID2EXTID, envir = .NCG_DOSE_Env)
        assign("PATHID2NAME", PATHID2NAME, envir = .NCG_DOSE_Env)
    }
    
    get(".NCG_DOSE_Env", envir = .DOSEEnv)
}


