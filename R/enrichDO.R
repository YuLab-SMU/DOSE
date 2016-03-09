##' DO Enrichment Analysis 
##'
##' Given a vector of genes, this function will return the enrichment DO
##' categories with FDR control.
##'
##'
##' @param ont one of DO or DOLite.
##' @inheritParams enrichNCG
##' @return A \code{enrichResult} instance.
##' @export
##' @seealso \code{\link{enrichResult-class}}
##' @author Guangchuang Yu \url{http://guangchuangyu.github.io}
##' @keywords manip
##' @examples
##'
##'	data(geneList)
##' 	gene = names(geneList)[geneList > 1]
##' 	yy = enrichDO(gene, pvalueCutoff=0.05)
##' 	summary(yy)
##'
enrichDO <- function(gene, ont="DO",
                     pvalueCutoff=0.05,
                     pAdjustMethod="BH",
                     universe,
                     minGSSize = 10,
                     maxGSSize = 500,
                     qvalueCutoff=0.2,
                     readable = FALSE){
    
    res <- enricher_internal(gene,
                             pvalueCutoff=pvalueCutoff,
                             pAdjustMethod=pAdjustMethod,
                             universe = universe,
                             minGSSize = minGSSize,
                             maxGSSize = maxGSSize,
                             qvalueCutoff = qvalueCutoff,
                             USER_DATA = get_DO_data(ont)
                             )

    if (is.null(res))
        return(res)
    
    res@organism <- "Homo sapiens"
    res@keytype <- "ENTREZID"
    res@ontology <- ont
    if(readable) {
        res <- setReadable(res, 'org.Hs.eg.db')
    }
    return(res)
}

get_DO_data <- function(ont="DO") {
    ont <- match.arg(ont, c("DO", "DOLite"))
    if (!exists("DOSEEnv")) {
        tryCatch(utils::data(list="DOSEEnv", package="DOSE"))
    }
    DOSEEnv <- get("DOSEEnv", envir = .GlobalEnv)
    if (ont == "DO") {
        PATHID2EXTID <- get("DO2ALLEG", envir = DOSEEnv)
        EXTID2PATHID <- get("EG2ALLDO", envir = DOSEEnv)
        PATH2NAME.df <- toTable(DOTERM)
        PATH2NAME.df <- PATH2NAME.df[, c("do_id", "Term")]
        PATH2NAME.df <- unique(PATH2NAME.df)
        PATH2NAME <- PATH2NAME.df[,2]
        names(PATH2NAME) <- PATH2NAME.df[,1]
    } else {
        PATHID2EXTID <- get("DOLite2EG", envir=DOSEEnv)
        EXTID2PATHID <- get("EG2DOLite", envir=DOSEEnv)
        PATH2NAME <- get("DOLiteTerm", envir=DOSEEnv)
    }
    
    assign("PATHID2EXTID", PATHID2EXTID, envir = DOSEEnv)
    assign("EXTID2PATHID", EXTID2PATHID, envir = DOSEEnv)
    assign("PATHID2NAME", PATH2NAME, envir = DOSEEnv)
    
    return(DOSEEnv)
}




