##' DO Enrichment Analysis 
##'
##' Given a vector of genes, this function will return the enrichment DO
##' categories with FDR control.
##'
##'
##' @param gene a vector of entrez gene id.
##' @param ont one of DO or DOLite.
##' @param pvalueCutoff Cutoff value of pvalue.
##' @param pAdjustMethod one of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"
##' @param universe background genes
##' @param minGSSize minimal size of genes annotated by Ontology term for testing.
##' @param qvalueCutoff qvalue Cutoff
##' @param readable whether mapping gene ID to gene Name
##' @return A \code{enrichResult} instance.
##' @export
##' @seealso \code{\link{enrichResult-class}}
##' @author Guangchuang Yu \url{http://ygc.name}
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
                     minGSSize = 5,
                     qvalueCutoff=0.2,
                     readable=FALSE) {

    enrich.internal(gene,
                    organism = "human",
                    pvalueCutoff=pvalueCutoff,
                    pAdjustMethod=pAdjustMethod,
                    ont = ont,
                    universe = universe,
                    minGSSize = minGSSize,
                    qvalueCutoff = qvalueCutoff,
                    readable = readable)
}


##################
##
##     DO
##
##################

##' @importMethodsFrom AnnotationDbi get
##' @importMethodsFrom AnnotationDbi exists
##' @method EXTID2TERMID DO
##' @export
EXTID2TERMID.DO <- function(gene, organism, ...) {
    if(!exists("DOSEEnv")) .initial()
    EG2ALLDO <- get("EG2ALLDO", envir=DOSEEnv)

    ## query external ID to Ontology ID
    qExtID2Term <- EG2ALLDO[gene]
    len <- sapply(qExtID2Term, length)
    notZero.idx <- len != 0
    qExtID2Term <- qExtID2Term[notZero.idx]

    return(qExtID2Term)
}


##' @importMethodsFrom AnnotationDbi get
##' @importMethodsFrom AnnotationDbi exists
##' @method TERMID2EXTID DO
##' @export
TERMID2EXTID.DO <- function(term, organism, ...) {
    if(!exists("DOSEEnv")) .initial()
    DO2ALLEG <- get("DO2ALLEG", envir=DOSEEnv)
    res <- DO2ALLEG[term]
    return(res)
}


##' @importMethodsFrom DO.db Term
##' @method TERM2NAME DO
##' @export
TERM2NAME.DO <- function(term, organism, ...) {
    desc = sapply(term, Term)
    return(desc)
}

##' @importMethodsFrom AnnotationDbi get
##' @importMethodsFrom AnnotationDbi exists
##' @method ALLEXTID DO
##' @export
ALLEXTID.DO <- function(organism, ...) {
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
##' @export
EXTID2TERMID.DOLite <- function(gene, organism, ...) {
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
##' @export
TERMID2EXTID.DOLite <- function(term, organism, ...) {
    if(!exists("DOSEEnv")) .initial()
    DOLite2EG <- get("DOLite2EG", envir=DOSEEnv)
    res <- DOLite2EG[term]
    return(res)
}


##' @method TERM2NAME DOLite
##' @export
TERM2NAME.DOLite <- function(term, organism, ...) {
    if(!exists("DOSEEnv")) .initial()
    DOLiteTerm <- get("DOLiteTerm", envir=DOSEEnv)
    desc = DOLiteTerm[term]
    return(desc)
}


##' @importMethodsFrom AnnotationDbi get
##' @importMethodsFrom AnnotationDbi exists
##' @method ALLEXTID DOLite
##' @export
ALLEXTID.DOLite <- function(organism, ...) {
    ##match.arg(organism, "human")
    if(!exists("DOSEEnv")) .initial()
    DOLite2EG <- get("DOLite2EG", envir=DOSEEnv)
    res <- unique(unlist(DOLite2EG))
    return(res)
}
