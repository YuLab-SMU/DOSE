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

    enrichDisease(gene = gene,
                  pvalueCutoff = pvalueCutoff,
                  pAdjustMethod = pAdjustMethod,
                  universe = universe,
                  minGSSize = minGSSize,
                  maxGSSize = maxGSSize,
                  qvalueCutoff = qvalueCutoff,
                  readable = readable,
                  ontology = ont)
}

get_DO_data <- function(ont="DO") {
    ont <- match.arg(ont, c("DO", "DOLite"))
    if (!exists(".DOSEEnv")) {
        .initial()
    }
    DOSEEnv <- get(".DOSEEnv", envir = .GlobalEnv)
    if (ont == "DO") {
        if (!exists("DO2ALLEG", envir=DOSEEnv)) {
            tryCatch(utils::data(list="DO2ALLEG", package="DOSE"))
            assign("DO2ALLEG", DO2ALLEG, envir = DOSEEnv)
            DO2ALLEG <- get("DO2ALLEG")
            rm(DO2ALLEG, envir = .GlobalEnv)
        }

        if (!exists("EG2ALLDO", envir = DOSEEnv)) {
            tryCatch(utils::data(list="EG2ALLDO", package="DOSE"))
            assign("EG2ALLDO", EG2ALLDO, envir = DOSEEnv)
            EG2ALLDO <- get("EG2ALLDO")
            rm(EG2ALLDO, envir = .GlobalEnv)
        }
            
        PATHID2EXTID <- get("DO2ALLEG", envir = DOSEEnv)
        EXTID2PATHID <- get("EG2ALLDO", envir = DOSEEnv)
        
        PATH2NAME.df <- toTable(DOTERM)
        PATH2NAME.df <- PATH2NAME.df[, c("do_id", "Term")]
        PATH2NAME.df <- unique(PATH2NAME.df)
        PATH2NAME <- PATH2NAME.df[,2]
        names(PATH2NAME) <- PATH2NAME.df[,1]
    } else {
        if (!exists("DOLite2EG", envir = DOSEEnv)) {
            tryCatch(utils::data(list="DOLite2EG", package="DOSE"))
            assign("DOLite2EG", DOLite2EG, envir = DOSEEnv)
            DOLite2EG <- get("DOLite2EG")
            rm(DOLite2EG, envir = .GlobalEnv)
        }

        if (!exists("EG2DOLite", envir = DOSEEnv)) {
            tryCatch(utils::data(list="EG2DOLite", package="DOSE"))
            assign("EG2DOLite", EG2DOLite, envir = DOSEEnv)
            EG2DOLite <- get("EG2DOLite")
            rm(EG2DOLite, envir = .GlobalEnv)
        }

        if (!exists("DOLiteTerm", envir = DOSEEnv)) {
            tryCatch(utils::data(list="DOLiteTerm", package="DOSE"))
            assign("DOLiteTerm", DOLiteTerm, envir = DOSEEnv)
            DOLiteTerm <- get("DOLiteTerm")
            rm(DOLiteTerm, envir = .GlobalEnv)
        }
        
        PATHID2EXTID <- get("DOLite2EG")
        EXTID2PATHID <- get("EG2DOLite")
        PATH2NAME <- get("DOLiteTerm")
    }
    
    assign("PATHID2EXTID", PATHID2EXTID, envir = DOSEEnv)
    assign("EXTID2PATHID", EXTID2PATHID, envir = DOSEEnv)
    assign("PATHID2NAME", PATH2NAME, envir = DOSEEnv)

    return(DOSEEnv)
}




