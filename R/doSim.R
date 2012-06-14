##' measuring similarities between two DO term vectors.
##'
##' provide two DO term vectors, this function will calculate their similarities.
##' @title doSim
##' @param DOID1 DO term vector
##' @param DOID2 DO term vector
##' @param method one of "Wang", "Resnik", "Rel", "Jiang", and "Lin".
##' @param organism only "human" supported
##' @return score matrix
##' @importFrom GOSemSim termSim
##' @export
##' @author Guangchuang Yu \url{http://ygc.name}
doSim <- function(DOID1,
                  DOID2,
                  method="Wang",
                  organism="human") {

    ont <- "DO"
    scores <- termSim(DOID1,DOID2, method, organism, ont)
    return(scores)
}
