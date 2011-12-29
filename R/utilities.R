.initial <- function() {
    assign("DOSEEnv", new.env(),.GlobalEnv)
    assign("SemSimCache", new.env(), .GlobalEnv)
    assign("ICEnv", new.env(), .GlobalEnv)

    tryCatch(utils::data(list="DO2EG", package="DOSE"))
    assign("DO2EG", DO2EG, envir=DOSEEnv)

    tryCatch(utils::data(list="EG2DO", package="DOSE"))
    assign("EG2DO", EG2DO, envir=DOSEEnv)

}


##' measuring similarities bewteen two gene vectors.
##'
##' provide two entrez gene vectors, this function will calculate their similarity.
##' @title geneSim
##' @param geneID1 entrez gene vector
##' @param geneID2 entrez gene vector
##' @param method one of "Wang", "Resnik", "Rel", "Jiang", and "Lin".
##' @param organism only "human" supported
##' @param combine One of "max", "average", "rcmax", "rcmax.avg" methods, for combining semantic similarity scores of multiple DO terms associated with gene/protein.
##' @return score matrix
##' @author Guangchuang Yu \url{http://ygc.name}
geneSim <- function(geneID1,
                    geneID2,
                    method="Wang",
                    organism="human",
                    combine="rcmax.avg") {

    DOID1 <- sapply(geneID1, gene2DO)
    DOID2 <- sapply(geneID2, gene2DO)
    m <- length(geneID1)
    n <- length(geneID2)
    scores <- matrix(NA, nrow=m, ncol=n)
    rownames(scores) <- geneID1
    colnames(scores) <- geneID2

    for (i in 1:m) {
        for (j in 1:n) {
            if(any(!is.na(DOID1[[i]])) &&  any(!is.na(DOID2[[j]]))) {
                s <- doSim(DOID1[[i]],
                           DOID2[[j]],
                           method,
                           organism
                           )
                scores[i,j] = combineScores(s, combine)
            }
        }
    }
    return(scores)
}


##' measuring similarities between two DO term vectors.
##'
##' provide two DO term vectors, this function will calculate their similarities.
##' @title doSim
##' @param DOID1 DO term vector
##' @param DOID2 DO term vector
##' @param method one of "Wang", "Resnik", "Rel", "Jiang", and "Lin".
##' @param organism only "human" supported
##' @return score matrix
##' @author Guangchuang Yu \url{http://ygc.name}
doSim <- function(DOID1,
                  DOID2,
                  method="Wang",
                  organism="human") {

    ont <- "DO"
    m <- length(DOID1)
    n <- length(DOID2)
    scores <- matrix(nrow=m, ncol=n)
    rownames(scores) <- DOID1
    colnames(scores) <- DOID2
    for( i in 1:m) {
        for (j in 1:n) {
            if ( is.na(DOID1[i]) || is.na(DOID2[j]) ) {
                scores[i,j] <- NA
            } else {
                if (method == "Wang") {
                    scores[i,j] <- wangMethod(DOID1[i],
                                              DOID2[j],
                                              ont=ont)
                } else {
                    scores[i,j] <- infoContentMethod(DOID1[i],
                                                     DOID2[j],
                                                     ont=ont,
                                                     method=method,
                                                     organism=organism)
                }
            }
        }
    }
    return(scores)
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


##' hypergeometric test for over-representation.
##'
##'
##' @title Hypergeometric test
##' @param numWdrawn number of White balls drawn
##' @param numW number of White balls
##' @param numB number of Black balls
##' @param numDrawn number of balls drawn
##' @return p value
##' @author Guangchuang Yu \url{http://ygc.name}
HyperG <- function(numWdrawn, numW, numB, numDrawn) {
    pvalue <- phyper(numWdrawn,
                     numW,
                     numB,
                     numDrawn,
                     lower.tail=FALSE)
    return(pvalue)
}

##' provide numerator and denominator, return numerator/denominator
##'
##'
##' @title getRatio
##' @param a numerator
##' @param b denominator
##' @return numerator/denominator
##' @author Guangchuang Yu \url{http://ygc.name}
getRatio <- function(a, b) {
    x=paste(a, "/", b, sep="", collapse="")
    return(x)
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

