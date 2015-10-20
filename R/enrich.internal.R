##' interal method for enrichment analysis
##'
##' using the hypergeometric model
##' @title enrich.internal
##' @param gene a vector of entrez gene id.
##' @param organism supported organism.
##' @param pvalueCutoff Cutoff value of pvalue.
##' @param pAdjustMethod one of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"
##' @param ont Ontology
##' @param universe background genes
##' @param minGSSize minimal size of genes annotated by Ontology term for testing.
##' @param qvalueCutoff cutoff of qvalue
##' @param readable whether mapping gene ID to gene Name
##' @param ... additional parameter
##' @return  A \code{enrichResult} instance.
##' @importClassesFrom methods data.frame
##' @importFrom plyr .
##' @importFrom plyr dlply
##' @importFrom qvalue qvalue
##' @importFrom methods new
##' @export
##' @keywords manip
##' @author Guangchuang Yu \url{http://ygc.name}
enrich.internal <- function(gene,
                            organism,
                            pvalueCutoff,
                            pAdjustMethod="BH",
                            ont,
                            universe,
                            minGSSize=5,
                            qvalueCutoff=0.2,
                            readable=FALSE,
                            ...) {

    ## query external ID to Term ID
    gene <- as.character(unique(gene))
    class(gene) <- ont
    qExtID2TermID = EXTID2TERMID(gene, organism, ...)
    qTermID <- unlist(qExtID2TermID)
    if (is.null(qTermID)) {
        return(NA)
    }

    ## Term ID -- query external ID association list.
    qExtID2TermID.df <- data.frame(extID=rep(names(qExtID2TermID),
                                   times=lapply(qExtID2TermID, length)),
                                   termID=qTermID)
    qExtID2TermID.df <- unique(qExtID2TermID.df)

    termID <- NULL ## to satisfy code tools
    qTermID2ExtID <- dlply(qExtID2TermID.df, .(termID),
                           .fun=function(i) as.character(i$extID))


    class(organism) <- ont
    extID <- ALLEXTID(organism, ...)
    if(!missing(universe)) {
        extID <- intersect(extID, universe)
    }

    qTermID2ExtID <- sapply(qTermID2ExtID, intersect, extID)

    ## Term ID annotate query external ID
    qTermID <- unique(names(qTermID2ExtID))

    class(qTermID) <- ont
    termID2ExtID <- TERMID2EXTID(qTermID, organism, ...)
    termID2ExtID <- sapply(termID2ExtID, intersect, extID)
    
    idx <- sapply(termID2ExtID, length) > minGSSize
    if (sum(idx) == 0)
        return (NULL)

    termID2ExtID <- termID2ExtID[idx]
    qTermID2ExtID <- qTermID2ExtID[idx]
    qTermID <- unique(names(qTermID2ExtID))
    
    ## prepare parameter for hypergeometric test
    k <- sapply(qTermID2ExtID, length)
    k <- k[qTermID]
    
    if (length(qTermID)== 1) {
        M <- nrow(termID2ExtID)
    } else {
        M <- sapply(termID2ExtID, length) 
        M <- M[qTermID]
    }

    N <- rep(length(extID), length(M))
    ## n <- rep(length(gene), length(M)) ## those genes that have no annotation should drop.
    n <- rep(length(qExtID2TermID), length(M))
    args.df <- data.frame(numWdrawn=k-1, ## White balls drawn
                          numW=M,        ## White balls
                          numB=N-M,      ## Black balls
                          numDrawn=n)    ## balls drawn


    ## calcute pvalues based on hypergeometric model
    pvalues <- apply(args.df, 1, function(n)
                     phyper(n[1], n[2], n[3], n[4], lower.tail=FALSE)
                     )

    ## gene ratio and background ratio
    GeneRatio <- apply(data.frame(a=k, b=n), 1, function(x)
                       paste(x[1], "/", x[2], sep="", collapse="")
                       )
    BgRatio <- apply(data.frame(a=M, b=N), 1, function(x)
                     paste(x[1], "/", x[2], sep="", collapse="")
                     )


    Over <- data.frame(ID=as.character(qTermID),
                       GeneRatio=GeneRatio,
                       BgRatio=BgRatio,
                       pvalue=pvalues)

    p.adj <- p.adjust(Over$pvalue, method=pAdjustMethod)
    qobj <- tryCatch(qvalue(p=Over$pvalue, lambda=0.05, pi0.method="bootstrap"), error=function(e) NULL)
    
    if (class(qobj) == "qvalue") {
        qvalues <- qobj$qvalues
    } else {
        qvalues <- NA
    }

    geneID <- sapply(qTermID2ExtID, function(i) paste(i, collapse="/"))
    geneID <- geneID[qTermID]
    Over <- data.frame(Over,
                       p.adjust = p.adj,
                       qvalue=qvalues,
                       geneID=geneID,
                       Count=k)


    class(qTermID) <- ont
    Description <- TERM2NAME(qTermID, organism, ...)
    
    if (length(qTermID) != length(Description)) {
        idx <- qTermID %in% names(Description)
        Over <- Over[idx,] 
    }
    Over$Description <- Description
    nc <- ncol(Over)
    Over <- Over[, c(1,nc, 2:(nc-1))]


    Over <- Over[order(pvalues),]

    Over <- Over[ Over$pvalue <= pvalueCutoff, ]
    Over <- Over[ Over$p.adjust <= pvalueCutoff, ]
    if (! any(is.na(Over$qvalue))) {
        Over <- Over[ Over$qvalue <= qvalueCutoff, ]
    }

    Over$ID <- as.character(Over$ID)
    Over$Description <- as.character(Over$Description)

    category <- as.character(Over$ID)

    row.names(Over) <- category

    
    x <- new("enrichResult",
             result         = Over,
             pvalueCutoff   = pvalueCutoff,
             pAdjustMethod  = pAdjustMethod,
             organism       = as.character(organism),
             ontology       = as.character(ont),
             gene           = as.character(gene),
             universe       = extID,
             geneInCategory = as.list(qTermID2ExtID[category]),
             geneSets       = termID2ExtID
             )
    if(readable)
        x <- setReadable(x)

    return (x)
}

