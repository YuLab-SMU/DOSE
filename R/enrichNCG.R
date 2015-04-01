##' NCG enrichment analysis
##'
##' given a vector of genes, this function will return the enrichment NCG
##' categories with FDR control
##'
##' 
##' @title enrichNCG
##' @param gene a vector of entrez gene id
##' @param pvalueCutoff pvalue cutoff
##' @param pAdjustMethod one of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"
##' @param universe background genes
##' @param minGSSize minimal size of genes annotated by NCG category for testing
##' @param qvalueCutoff qvalue cutoff
##' @param readable whether mapping gene ID to gene Name
##' @return A \code{enrichResult} instance
##' @export
##' @author Guangchuang Yu
enrichNCG <- function(gene,
                      pvalueCutoff=0.05,
                      pAdjustMethod="BH",
                      universe,
                      minGSSize=5,
                      qvalueCutoff=0.2,
                      readable=FALSE) {
    
    NCG_DOSE_Env <- get_NCG_data()
    enrich.internal(gene = gene,
                    organism = "human",
                    pvalueCutoff = pvalueCutoff,
                    pAdjustMethod = pAdjustMethod,
                    ont = "NCG",
                    universe = universe,
                    minGSSize = minGSSize,
                    qvalueCutoff = qvalueCutoff,
                    readable = readable,
                    USER_DATA = NCG_DOSE_Env)
}

get_NCG_data <- function() {
    if (!exists("NCG_DOSE_Env")) {
        tryCatch(utils::data(list="NCG_DOSE_Env", package="DOSE"))
    }
    get("NCG_DOSE_Env", envir = .GlobalEnv)
}


##' @method EXTID2TERMID NCG
##' @export
EXTID2TERMID.NCG <- function(gene, organism, ...) {
    EXTID2TERMID.USER_DEFINED(gene, organism, ...)
}

##' @method TERMID2EXTID NCG
##' @export
TERMID2EXTID.NCG <- function(term, organism, ...) {
    TERMID2EXTID.USER_DEFINED(term, organism, ...)
}

##' @method TERM2NAME NCG
##' @export
TERM2NAME.NCG <- function(term, organism, ...) {
    params <- as.list(match.call()[-1])
    if ("USER_DATA" %in% names(params)) {
        TERM2NAME.USER_DEFINED(term, organism, ...)
    } else {
        USER_DATA <- get_NCG_data()
        TERM2NAME.USER_DEFINED(term, organism, USER_DATA)        
    }
}

##' @method ALLEXTID NCG
##' @export
ALLEXTID.NCG <- function(organism, ...) {
    ALLEXTID.USER_DEFINED(organism, ...)
}


##' @method EXTID2TERMID USER_DEFINED
##' @export
EXTID2TERMID.USER_DEFINED <- function(gene, organism, ...) {
    EXTID2TERMID.USER_DEFINED.internal(gene, organism, ...)
}

##' @method TERMID2EXTID USER_DEFINED
##' @export
TERMID2EXTID.USER_DEFINED <- function(term, organism, ...) {
    TERMID2EXTID.USER_DEFINED.internal(term, organism, ...)
}

##' @method TERM2NAME USER_DEFINED
##' @export
TERM2NAME.USER_DEFINED <- function(term, organism, ...) {
    TERM2NAME.USER_DEFINED.internal(term, organism, ...)
}

##' @method ALLEXTID USER_DEFINED
##' @export
ALLEXTID.USER_DEFINED <- function(organism, ...) {
    ALLEXTID.USER_DEFINED.internal(organism, ...)
}


EXTID2TERMID.USER_DEFINED.internal <- function(gene, organism, USER_DATA, ...) {
    EXTID2PATHID <- get("EXTID2PATHID", envir = USER_DATA)

    qExtID2Path <- EXTID2PATHID[gene]
    len <- sapply(qExtID2Path, length)
    notZero.idx <- len != 0
    qExtID2Path <- qExtID2Path[notZero.idx]

    return(qExtID2Path)
}

TERMID2EXTID.USER_DEFINED.internal <- function(term, organism, USER_DATA, ...) {
    PATHID2EXTID <- get("PATHID2EXTID", envir = USER_DATA)
    res <- PATHID2EXTID[term]
    return(res)
}

TERM2NAME.USER_DEFINED.internal <- function(term, organism, USER_DATA, ...) {
    PATHID2NAME <- get("PATHID2NAME", envir = USER_DATA)
    if (is.null(PATHID2NAME) || is.na(PATHID2NAME)) {
        return(as.character(term))
    }
    return(PATHID2NAME[term])
}

ALLEXTID.USER_DEFINED.internal <- function(organism, USER_DATA, ...) {
    PATHID2EXTID <- get("PATHID2EXTID", envir = USER_DATA)
    res <- unique(unlist(PATHID2EXTID))
    return(res)
}
