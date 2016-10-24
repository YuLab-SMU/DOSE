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
    
    if (!exists(".DGN_DOSE_Env", envir=.DOSEEnv)) {
        tryCatch(utils::data(list="DGN_EXTID2PATHID", package="DOSE"))
        tryCatch(utils::data(list="DGN_PATHID2EXTID", package="DOSE"))
        tryCatch(utils::data(list="DGN_PATHID2NAME", package="DOSE"))
        EXTID2PATHID <- DGN_EXTID2PATHID <- get("DGN_EXTID2PATHID")
        PATHID2EXTID <- DGN_PATHID2EXTID <- get("DGN_PATHID2EXTID")
        PATHID2NAME <- DGN_PATHID2NAME <- get("DGN_PATHID2NAME")

        rm(DGN_EXTID2PATHID, envir = .GlobalEnv)
        rm(DGN_PATHID2EXTID, envir = .GlobalEnv)
        rm(DGN_PATHID2NAME, envir = .GlobalEnv)

        assign(".DGN_DOSE_Env", new.env(), envir = .DOSEEnv)
        .DGN_DOSE_Env <- get(".DGN_DOSE_Env", envir = .DOSEEnv)
        assign("EXTID2PATHID", EXTID2PATHID, envir = .DGN_DOSE_Env)
        assign("PATHID2EXTID", PATHID2EXTID, envir = .DGN_DOSE_Env)
        assign("PATHID2NAME", PATHID2NAME, envir = .DGN_DOSE_Env)
    }
    get(".DGN_DOSE_Env", envir = .DOSEEnv)
}


