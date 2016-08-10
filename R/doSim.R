##' measuring similarities between two DO term vectors.
##'
##' provide two DO term vectors, this function will calculate their similarities.
##' @title doSim
##' @param DOID1 DO term vector
##' @param DOID2 DO term vector
##' @param measure one of "Wang", "Resnik", "Rel", "Jiang", and "Lin".
##' @return score matrix
##' @importFrom GOSemSim termSim
##' @export
##' @author Guangchuang Yu \url{https://guangchuangyu.github.io}
doSim <- function(DOID1,
                  DOID2,
                  measure="Wang") {
    scores <- termSim(DOID1,DOID2, dodata(), measure)
    if(nrow(scores) == 1 & ncol(scores) == 1)
        scores <- as.numeric(scores)
    return(scores)
}


