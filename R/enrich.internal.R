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
                            qvalueCutoff,
                            ont,
                            readable) {

    ## query external ID to Term ID
    gene <- as.character(gene)
    class(gene) <- ont
    qExtID2TermID = EXTID2TERMID(gene, organism)
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

    ## Term ID annotate query external ID
    qTermID <- unique(qTermID)

    ## prepare parameter for hypergeometric test
    k <- sapply(qTermID2ExtID, length)
    k <- k[qTermID]

    class(qTermID) <- ont
    termID2ExtID <- TERMID2EXTID(qTermID, organism)
    M <- sapply(termID2ExtID, length)
    M <- M[qTermID]

    class(organism) <- ont
    extID <- ALLEXTID(organism)
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

    class(qTermID) <- ont
    Description <- TERM2NAME(qTermID)

    Over <- data.frame(ID=as.character(qTermID),
                       Description=Description,
                       GeneRatio=GeneRatio,
                       BgRatio=BgRatio,
                       pvalue=pvalues)

    qobj = qvalue(p=Over$pvalue, lambda=0.05, pi0.method="bootstrap")
    qvalues <- qobj$qvalues



    geneID <- sapply(qTermID2ExtID, function(i) paste(i, collapse="/"))
    geneID <- geneID[qTermID]
    Over <- data.frame(Over,
                       qvalue=qvalues,
                       geneID=geneID,
                       Count=k)


    Over <- Over[order(pvalues),]

    Over <- Over[ Over$pvalue <= pvalueCutoff, ]
    Over <- Over[ Over$qvalue <= qvalueCutoff, ]

    Over$ID <- as.character(Over$ID)
    Over$Description <- as.character(Over$Description)

    category <- as.character(Over$ID)

    rownames(Over) <- category


    x <- new("enrichResult",
             result = Over,
             pvalueCutoff=pvalueCutoff,
             qvalueCutoff=qvalueCutoff,
             organism = as.character(organism),
             ontology = as.character(ont),
             gene = as.character(gene),
             geneInCategory = qTermID2ExtID[category]
             )

    setReadable(x) <- readable

    return (x)
}

##################
##
##     DO
##
##################

##' @importMethodsFrom AnnotationDbi get
##' @importMethodsFrom AnnotationDbi exists
##' @method EXTID2TERMID DO
EXTID2TERMID.DO <- function(gene, organism) {
    if(!exists("DOSEEnv")) .initial()
    EG2ALLDO <- get("EG2ALLDO", envir=DOSEEnv)

    ## query external ID to Ontology ID
    qExtID2Term= EG2ALLDO[gene]
    len <- sapply(qExtID2Term, length)
    notZero.idx <- len != 0
    qExtID2Term <- qExtID2Term[notZero.idx]

    return(qExtID2Term)
}

##' @importMethodsFrom AnnotationDbi get
##' @importMethodsFrom AnnotationDbi exists
##' @method TERMID2EXTID DO
TERMID2EXTID.DO <- function(term, organism) {
    if(!exists("DOSEEnv")) .initial()
    DO2ALLEG <- get("DO2ALLEG", envir=DOSEEnv)
    res <- DO2ALLEG[term]
    return(res)
}

##' @importMethodsFrom DO.db Term
##' @method TERM2NAME DO
TERM2NAME.DO <- function(term) {
    desc = sapply(term, Term)
    return(desc)
}

##' @importMethodsFrom AnnotationDbi get
##' @importMethodsFrom AnnotationDbi exists
##' @method ALLEXTID DO
ALLEXTID.DO <- function(organism) {
    ##match.arg(organism, "human")
    if(!exists("DOSEEnv")) .initial()
    DO2ALLEG <- get("DO2ALLEG", envir=DOSEEnv)
    res <- unique(unlist(DO2ALLEG))
    return(res)
}


##################
##
##     DOLite
##
##################
##' @importMethodsFrom AnnotationDbi get
##' @importMethodsFrom AnnotationDbi exists
##' @method EXTID2TERMID DOLite
EXTID2TERMID.DOLite <- function(gene, organism) {
    if(!exists("DOSEEnv")) .initial()
    EG2DOLite <- get("EG2DOLite", envir=DOSEEnv)

    ## query external ID to Ontology ID
    qExtID2Term= EG2DOLite[gene]
    len <- sapply(qExtID2Term, length)
    notZero.idx <- len != 0
    qExtID2Term <- qExtID2Term[notZero.idx]

    return(qExtID2Term)
}

##' @importMethodsFrom AnnotationDbi get
##' @importMethodsFrom AnnotationDbi exists
##' @method TERMID2EXTID DOLite
TERMID2EXTID.DOLite <- function(term, organism) {
    if(!exists("DOSEEnv")) .initial()
    DOLite2EG <- get("DOLite2EG", envir=DOSEEnv)
    res <- DOLite2EG[term]
    return(res)
}

##' @method TERM2NAME DOLite
TERM2NAME.DOLite <- function(term) {
    if(!exists("DOSEEnv")) .initial()
    DOLiteTerm <- get("DOLiteTerm", envir=DOSEEnv)
    desc = DOLiteTerm[term]
    return(desc)
}

##' @importMethodsFrom AnnotationDbi get
##' @importMethodsFrom AnnotationDbi exists
##' @method ALLEXTID DOLite
ALLEXTID.DOLite <- function(organism) {
    ##match.arg(organism, "human")
    if(!exists("DOSEEnv")) .initial()
    DOLite2EG <- get("DOLite2EG", envir=DOSEEnv)
    res <- unique(unlist(DOLite2EG))
    return(res)
}
