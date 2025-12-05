##' Enrichment analysis based on the DisGeNET (\url{http://www.disgenet.org/})
##'
##' given a vector of genes, this function will return the enrichment NCG
##' categories with FDR control
##'
##'
##' @inheritParams enrichNCG
##' @return A \code{enrichResult} instance
##' @export
##' @references Janet et al. (2015) DisGeNET: a discovery platform for the dynamical exploration of human diseases and their genes. \emph{Database} bav028
##' \url{http://database.oxfordjournals.org/content/2015/bav028.long}
##' @author Guangchuang Yu
enrichDGN <- function(gene,
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
                  ontology = "DisGeNET")
    
}

get_DGN_data <- function() {
    if (!exists(".DOSEenv")) .initial()
    .DOSEEnv <- get(".DOSEEnv", envir = .GlobalEnv)
    
    if (exists(".DGN_DOSE_GSON", envir=.DOSEEnv)) {
        res <- get(".DGN_DOSE_GSON", envir = .DOSEEnv)
        return(res)
    }

    tryCatch(utils::data(list="DGN_PATHID2EXTID", package="DOSE"))
    tryCatch(utils::data(list="DGN_PATHID2NAME", package="DOSE"))
    PATHID2EXTID <- get("DGN_PATHID2EXTID")
    PATHID2NAME <- get("DGN_PATHID2NAME")

    rm(DGN_PATHID2EXTID, envir = .GlobalEnv)
    rm(DGN_PATHID2NAME, envir = .GlobalEnv)

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
                           gsname = "DisGeNET",
                           keytype = "ENTREZID",
                           version = "unknown",
                           accessed_date = as.character(Sys.Date()))

    assign(".DGN_DOSE_GSON", gson_obj, envir = .DOSEEnv)
    return(gson_obj)
}


