##' measuring similarities bewteen two gene vectors.
##'
##' provide two entrez gene vectors, this function will calculate their similarity.
##' @title geneSim
##' @param geneID1 entrez gene vector
##' @param geneID2 entrez gene vector
##' @param measure one of "Wang", "Resnik", "Rel", "Jiang", and "Lin".
##' @param combine One of "max", "avg", "rcmax", "BMA" methods, for combining semantic similarity scores of multiple DO terms associated with gene/protein.
##' @return score matrix
##' @importFrom DO.db DOPARENTS
##' @importFrom DO.db DOANCESTOR
##' @importFrom GOSemSim combineScores
##' @export
##' @author Guangchuang Yu \url{http://ygc.name}
geneSim <- function(geneID1,
                    geneID2=NULL,
                    measure="Wang",
                    combine="BMA") {

    DOID1 <- lapply(geneID1, gene2DO)
    if (is.null(geneID2)) {
        geneID2 <- geneID1
        DOID2 <- DOID1
    } else {
        DOID2 <- lapply(geneID2, gene2DO)
    }

    m <- length(geneID1)
    n <- length(geneID2)
    scores <- matrix(NA, nrow=m, ncol=n)
    rownames(scores) <- geneID1
    colnames(scores) <- geneID2

    for (i in 1:m) {
        if (length(geneID1) == length(geneID2) && all(geneID1 == geneID2)) {
           nn <- i
           flag <- TRUE
        } else {
            flag <- FALSE
            nn <- n
        }

        for (j in 1:nn) {
            if(any(!is.na(DOID1[[i]])) &&  any(!is.na(DOID2[[j]]))) {
                s <- doSim(DOID1[[i]],
                           DOID2[[j]],
                           measure
                           )
                scores[i,j] = combineScores(s, combine)
                if (flag == TRUE && j != i) {
                    scores[j, i] <- scores[i,j]
                }
            }
        }
    }
    if (nrow(scores) == 1 & ncol(scores) == 1)
        scores = as.numeric(scores)
    return(scores)
}
