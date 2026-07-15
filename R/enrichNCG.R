#' Enrichment analysis based on the Network of Cancer Genes database (https://www.network-cancer-genes.org/)
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
    .DOSEEnv <- get_dose_env()
    
    if (exists(".NCG_DOSE_GSON", envir=.DOSEEnv)) {
        res <- get(".NCG_DOSE_GSON", envir = .DOSEEnv)
        return(res)
    }

    urls <- c("https://yulab-smu.top/DOSE",
              "https://raw.githubusercontent.com/YuLab-SMU/DOSE/refs/heads/gh-pages")
    
    dbfile <- yulab.utils::download_yulab_file("NCG.gson", urls, 
                                               gzfile = TRUE, appname = "DOSE")
    gson_obj <- gson::read.gson(dbfile)

    assign(".NCG_DOSE_GSON", gson_obj, envir = .DOSEEnv)
    return(gson_obj)
}
