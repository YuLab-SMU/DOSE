##' mapping geneID to gene Symbol
##'
##'
##' @title setReadable
##' @param x enrichResult Object
##' @param OrgDb OrgDb
##' @param keyType keyType of gene
##' @return enrichResult Object
##' @author Yu Guangchuang
##' @export
setReadable <- function(x, OrgDb, keyType="auto") {
    OrgDb <- load_OrgDb(OrgDb)
    if (!'SYMBOL' %in% columns(OrgDb)) {
        warning("Fail to convert input geneID to SYMBOL since no SYMBOL information available in the provided OrgDb...")
    }

    if (!(is(x, "enrichResult") || is(x, "groupGOResult") || is(x, "gseaResult") || is(x,"compareClusterResult")))
        stop("input should be an 'enrichResult' , 'gseaResult' or 'compareClusterResult' object...")

    isGSEA <- FALSE
    isCompare <- FALSE
    if (is(x, 'gseaResult'))
        isGSEA <- TRUE
        
    if (is(x, 'compareClusterResult'))
        isCompare <- TRUE
    
    if (keyType == "auto" & isCompare == FALSE) {
        keyType <- x@keytype
        if (keyType == 'UNKNOWN') {
            stop("can't determine keyType automatically; need to set 'keyType' explicitly...")
        }
    }

    if (isCompare == FALSE && x@readable)
        return(x)

    gc <- geneInCategory(x)

    if (isGSEA) {
        genes <- names(x@geneList)
    } else if (isCompare) {
        genes <- NULL
        for(i in seq_len(length(x@geneClusters))) {
            genes <- c(genes,x@geneClusters[[i]])
        }
        genes <- unique(genes)
    } else {
        genes <- x@gene
    }

    gn <- EXTID2NAME(OrgDb, genes, keyType)
    gc <- lapply(gc, function(i) gn[i])

    if(isCompare) {
        res <- x@compareClusterResult
    } else {
        res <- x@result
    }
    
    gc <- gc[as.character(res$ID)]
    geneID <- sapply(gc, paste0, collapse="/")
    if (isGSEA) {
        res$core_enrichment <- unlist(geneID)
    } else {
        res$geneID <- unlist(geneID)
    }

    if(isCompare){
        x@compareClusterResult <- res    
    } else {
        x@gene2Symbol <- gn
        x@result <- res
        x@keytype <- keyType
        x@readable <- TRUE
    }


    return(x)
}
