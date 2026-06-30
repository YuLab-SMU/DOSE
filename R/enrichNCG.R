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
#' @importFrom utils read.delim
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

    urls <- c("https://yulab-smu.top/DOSE",
              "https://raw.githubusercontent.com/YuLab-SMU/DOSE/refs/heads/gh-pages")
    
    dbfile <- yulab.utils::download_yulab_file("NCG.tsv.gz", urls, 
                                               gzfile = FALSE, appname = "DOSE")
    
    ncg <- read.delim(gzfile(dbfile), stringsAsFactors = FALSE)
    PATHID2EXTID <- split(as.character(ncg$entrez), as.character(ncg$cancer_type))

    # gsid2gene
    gsid2gene <- stack(PATHID2EXTID)
    colnames(gsid2gene) <- c("gene", "gsid")
    gsid2gene <- gsid2gene[, c("gsid", "gene")]

    # gsid2name: use cancer_type as both ID and name
    gsid2name <- data.frame(gsid = unique(ncg$cancer_type),
                            name = unique(ncg$cancer_type),
                            stringsAsFactors = FALSE)

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
