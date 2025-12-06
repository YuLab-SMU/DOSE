#' Enrichment analysis based on the Network of Cancer Genes database (http://ncg.kcl.ac.uk/)
#'
#' given a vector of genes, this function will return the enrichment NCG
#' categories with FDR control
#'
#' 
#' @title enrichNCG
#' @inheritParams dose_params
#' @return A \code{enrichResult} instance
#' @export
#' @author Guangchuang Yu
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
    
    if (exists(".NCG_DOSE_GSON", envir=.DOSEEnv)) {
        res <- get(".NCG_DOSE_GSON", envir = .DOSEEnv)
        return(res)
    }
    
    tryCatch(utils::data(list="NCG_PATHID2EXTID", package="DOSE"))
    tryCatch(utils::data(list="NCG_PATHID2NAME", package="DOSE"))
    PATHID2EXTID <- get("NCG_PATHID2EXTID")
    PATHID2NAME <- get("NCG_PATHID2NAME")

    rm(NCG_PATHID2EXTID, envir = .GlobalEnv)
    rm(NCG_PATHID2NAME, envir = .GlobalEnv)

    # gsid2gene
    gsid2gene <- stack(PATHID2EXTID)
    colnames(gsid2gene) <- c("gene", "gsid")
    gsid2gene <- gsid2gene[, c("gsid", "gene")]

    # gsid2name
    gsid2name <- data.frame(gsid = names(PATHID2NAME), name = PATHID2NAME)
    rownames(gsid2name) <- NULL

    gson_obj <- gson::gson(gsid2gene = gsid2gene, 
                           gsid2name = gsid2name,
                           species = "Homo sapiens",
                           gsname = "NCG",
                           keytype = "ENTREZID",
                           version = "unknown",
                           accessed_date = as.character(Sys.Date()))

    assign(".NCG_DOSE_GSON", gson_obj, envir = .DOSEEnv)
    return(gson_obj)
}
