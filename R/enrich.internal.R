##' @importMethodsFrom AnnotationDbi get
##' @importMethodsFrom AnnotationDbi exists
EXTID2TERMID <- function(gene, organism="human", ont="DO") {
    match.arg(ont, "DO")

    if(!exists("DOSEEnv")) .initial()
    EG2DO <- get("EG2DO", envir=DOSEEnv)

    ## query external ID to Ontology ID
    qExtID2Term= EG2DO[gene]
    len <- sapply(qExtID2Term, length)
    notZero.idx <- len != 0
    qExtID2Term <- qExtID2Term[notZero.idx]

    return(qExtID2Term)
}


TERMID2EXTID <- function(TermID, organism) {
    res <- DO2EG[TermID]
    return(res)
}

ALLEXTID <- function(ont="DO") {
    if (ont == "DO") {
        res <- unique(unlist(DO2EG))
    }
    return(res)
}

##' @importMethodsFrom DO.db Term
TERM2NAME <- function(term, ont="DO") {
    if (ont == "DO") {
        desc = sapply(term, Term)
    }
    return(desc)
}

##' convert gene IDs to gene Names.
##'
##'
##' @title convert gene IDs to gene Names
##' @param geneID a vector of gene IDs
##' @param organism only "human" supported
##' @return a vector of gene names.
##' @importFrom org.Hs.eg.db org.Hs.egSYMBOL
##' @author Guangchuang Yu \url{http://ygc.name}
##' @export
EXTID2NAME <- function(geneID, organism) {
    annotation <- switch(organism,
                         human = org.Hs.egSYMBOL
                         )
    gn <- unique(unlist(mget(geneID, annotation, ifnotfound=NA)))
    return(gn)
}

##' interal method for enrichment analysis
##'
##' using the hypergeometric model
##' @title enrich.internal
##' @param gene a vector of entrez gene id.
##' @param organism supported organism.
##' @param pvalueCutoff Cutoff value of pvalue.
##' @param qvalueCutoff Cutoff value of qvalue.
##' @param ont Ontology
##' @param readable whether mapping gene ID to gene Name
##' @return  A \code{enrichDOResult} instance.
##' @importClassesFrom methods data.frame
##' @importFrom plyr .
##' @importFrom plyr dlply
##' @importFrom qvalue qvalue
##' @importFrom methods new
##' @export
##' @keywords manip
##' @author Guangchuang Yu \url{http://ygc.name}
enrich.internal <- function(gene, organism, pvalueCutoff, qvalueCutoff, ont, readable) {
    qExtID2TermID = EXTID2TERMID(gene)
    qTermID <- unlist(qExtID2TermID)

    qExtID2TermID.df <- data.frame(extID=rep(names(qExtID2TermID), times=lapply(qExtID2TermID, length)),
                                   termID=qTermID)
    qExtID2TermID.df <- unique(qExtID2TermID.df)

    termID <- NULL ## to satisfy code tools
    qTermID2ExtID <- dlply(qExtID2TermID.df, .(termID), .fun=function(i) as.character(i$extID))

    qTermID <- unique(qTermID)
    k <- sapply(qTermID2ExtID, length)
    k <- k[qTermID]

    termID2ExtID <- TERMID2EXTID(qTermID)

    M <- sapply(termID2ExtID, length)
    M <- M[qTermID]

    extID <- ALLEXTID(ont)

    N <- rep(length(extID), length(M))
    n <- rep(length(gene), length(M))
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


    Description <- TERM2NAME(qTermID)

    Over <- data.frame(ID=qTermID,
                       Description=Description,
                       GeneRatio=GeneRatio,
                       BgRatio=BgRatio,
                       pvalue=pvalues)


    qobj = qvalue(Over$pvalue, lambda=0.05, pi0.method="bootstrap")
    qvalues <- qobj$qvalues

    if(readable) {
        qTermID2ExtID <- lapply(qTermID2ExtID, EXTID2NAME)
    }

    geneID <- sapply(qTermID2ExtID, function(i) paste(i, collapse="/"))
    geneID <- geneID[qTermID]
    Over <- data.frame(Over,
                       qvalue=qvalues,
                       geneID=geneID,
                       Count=k)

    Over <- Over[order(pvalues),]

    Over <- Over[ Over$pvalue <= pvalueCutoff, ]
    Over <- Over[ Over$qvalue <= qvalueCutoff, ]

    Over$Description <- as.character(Over$Description)



    new("enrichResult",
        result = Over,
        pvalueCutoff=pvalueCutoff,
        qvalueCutoff=qvalueCutoff,
        organism = organism,
        gene = gene,
        geneInCategory = qTermID2ExtID
	)
}
