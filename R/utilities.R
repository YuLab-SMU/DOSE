.initial <- function() {
    assign("DOSEEnv", new.env(),.GlobalEnv)
    assign("SemSimCache", new.env(), .GlobalEnv)
    assign("ICEnv", new.env(), .GlobalEnv)

    tryCatch(utils::data(list="DO2EG", package="DOSE"))
    assign("DO2EG", DO2EG, envir=DOSEEnv)

    tryCatch(utils::data(list="EG2DO", package="DOSE"))
    assign("EG2DO", EG2DO, envir=DOSEEnv)

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
    EG2DO <- get("EG2DO", envir=DOSEEnv)
    DO <- EG2DO[[gene]]
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
    idx <- grep("DOID", names(DO2EG))
    DO2EG <- DO2EG[idx]
    save(DO2EG, file="DO2EG.rda")

    EG2DO <- dlply(eg.do, .(eg), .fun=function(i) i$doid)
    i=unlist(lapply(EG2DO, function(i) length(grep("DOID", i)))) != 0
    EG2DO <- EG2DO[i]
    save(EG2DO, file="EG2DO.rda")
}
