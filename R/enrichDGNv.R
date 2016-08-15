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
    if (!exists(".DOSEenv")) .initial()
    .DOSEEnv <- get(".DOSEEnv", envir = .GlobalEnv)
    
    if (!exists(".VDGN_DOSE_Env", envir=.DOSEEnv)) {
        tryCatch(utils::data(list="VDGN_EXTID2PATHID", package="DOSE"))
        tryCatch(utils::data(list="VDGN_PATHID2EXTID", package="DOSE"))
        tryCatch(utils::data(list="VDGN_PATHID2NAME", package="DOSE"))
        EXTID2PATHID <- VDGN_EXTID2PATHID <- get("VDGN_EXTID2PATHID")
        PATHID2EXTID <- VDGN_PATHID2EXTID <- get("VDGN_PATHID2EXTID")
        PATHID2NAME <- VDGN_PATHID2NAME <- get("VDGN_PATHID2NAME")

        rm(VDGN_EXTID2PATHID, envir = .GlobalEnv)
        rm(VDGN_PATHID2EXTID, envir = .GlobalEnv)
        rm(VDGN_PATHID2NAME, envir = .GlobalEnv)

        assign(".VDGN_DOSE_Env", new.env(), envir = .DOSEEnv)
        .VDGN_DOSE_Env <- get(".VDGN_DOSE_Env", envir = .DOSEEnv)
        assign("EXTID2PATHID", EXTID2PATHID, envir = .VDGN_DOSE_Env)
        assign("PATHID2EXTID", PATHID2EXTID, envir = .VDGN_DOSE_Env)
        assign("PATHID2NAME", PATHID2NAME, envir = .VDGN_DOSE_Env)
    }
    get(".VDGN_DOSE_Env", envir = .DOSEEnv)
}


