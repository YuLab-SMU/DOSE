.initial <- function() {
    assign("DOSEEnv", new.env(),.GlobalEnv)
    assign("SemSimCache", new.env(), .GlobalEnv)
    assign("ICEnv", new.env(), .GlobalEnv)

    tryCatch(utils::data(list="DO2ALLEG", package="DOSE"))
    assign("DO2ALLEG", DO2ALLEG, envir=DOSEEnv)

    tryCatch(utils::data(list="EG2ALLDO", package="DOSE"))
    assign("EG2ALLDO", EG2ALLDO, envir=DOSEEnv)

}

##' provide gene ID, this function will convert to the corresponding DO Terms
##'
##'
##' @title convert Gene ID to DO Terms
##' @param gene entrez gene ID
##' @return DO Terms
##' @importMethodsFrom AnnotationDbi get
##' @importMethodsFrom AnnotationDbi exists
##' @author Guangchuang Yu \url{http://ygc.name}
gene2DO <- function(gene) {
    if(!exists("DOSEEnv")) .initial()
    EG2ALLDO <- get("EG2ALLDO", envir=DOSEEnv)
    DO <- EG2ALLDO[[gene]]
    if (is.null(DO)) {
        return(NA)
    }
    if (sum(!is.na(DO)) == 0) {
        return(NA)
    }
    if (length(DO) == 0) {
        return(NA)
    }
    return(DO)
}

##' rebuilding entrez and DO mapping datasets
##'
##'
##' @title rebuiding annotation data
##' @param file do_rif.human.txt
##' @return NULL
##' @importFrom plyr dlply
##' @importFrom plyr .
##' @importFrom DO.db DOANCESTOR
##' @importFrom DO.db DOTERM
##' @importFrom AnnotationDbi mget
##' @author Guangchuang Yu \url{http://ygc.name}
rebuildAnnoData <- function(file) {
                                        #
                                        # do_rif.human.txt was downloaded from
                                        # http://projects.bioinformatics.northwestern.edu/do_rif/
                                        #
    do.rif <- read.delim2(file, sep="\t", stringsAsFactors=F, header=F)
    eg.do <- do.rif[,c(1,5)]
    colnames(eg.do) <- c("eg", "doid")

    eg <- doid <- NULL # to satisfy codetools

    DO2EG <- dlply(eg.do, .(doid), .fun=function(i) i$eg)
    doids <- toTable(DOTERM)
    doterms <- doids$do_id
    idx <- names(DO2EG) %in% doterms
    DO2EG <- DO2EG[idx]
    DO2EG <- lapply(DO2EG, function(i) unique(i))
    save(DO2EG, file="DO2EG.rda")


    EG2DO <- dlply(eg.do, .(eg), .fun=function(i) i$doid)
    EG2DO <- lapply(EG2DO, function(i) unique(i[ i %in% doterms ]))

    i <- unlist(lapply(EG2DO, function(i) length(i) != 0))
    EG2DO <- EG2DO[i]
    save(EG2DO, file="EG2DO.rda")

    EG2ALLDO <- lapply(EG2DO,
                       function(i) {
                           ans <- unlist(mget(i, DOANCESTOR))
                           ans <- ans[ !is.na(ans) ]
                           ans <- c(i, ans)
                           ans <- unique(ans)
                           return(ans)
                       })
    save(EG2ALLDO, file="EG2ALLDO.rda")

    len <- lapply(EG2ALLDO,length)
    EG2ALLDO.df <- data.frame(EG=rep(names(EG2ALLDO), times=len),
                              DO=unlist(EG2ALLDO))
    DO <- NULL ## satisfy code tools
    DO2ALLEG <- dlply(EG2ALLDO.df, .(DO), function(i) as.character(i$EG))
    DO2ALLEG <- lapply(DO2ALLEG, unique)
    save(DO2ALLEG, file="DO2ALLEG.rda")
}
